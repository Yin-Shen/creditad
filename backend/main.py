import os
import sys
import io

SERVER_PORT = int(os.environ.get("CREDITAD_PORT", "5001"))

# Ensure backend dir is on sys.path so embedded Python finds local modules
_BACKEND_DIR = os.path.dirname(os.path.abspath(__file__))
if _BACKEND_DIR not in sys.path:
    sys.path.insert(0, _BACKEND_DIR)
os.chdir(_BACKEND_DIR)

# 核心修复：强制标准输出和错误流使用 UTF-8 编码 (Win console encoding)
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

import startup_log
startup_log.banner()

import asyncio
import uuid
import json
from concurrent.futures import ThreadPoolExecutor

from tad_vci.runtime_contract import MethodSetContractError, resolve_method_set

with startup_log.step("deps", "Loading numerical stack (numpy / pandas / h5py)"):
    import pandas as pd
    import numpy as np
    import h5py

with startup_log.step("deps", "Loading aiohttp web framework"):
    from aiohttp import web
    import aiohttp_cors

with startup_log.step("deps", "Loading TAD detection (mcool_bed)"):
    from mcool_bed import detect_TAD_boundaries_v2, detect_TAD_boundaries_tadvci

with startup_log.step("deps", "Loading SQLite layer"):
    from sql import init_db, add_token_record, update_tad_outputs, mark_running, mark_failed, get_record

with startup_log.step("deps", "Loading optional Hi-C reader (.hic via hicstraw, optional)"):
    try:
        from hic_reader import HicFileReader, HIC_AVAILABLE
    except ImportError:
        HIC_AVAILABLE = False
        startup_log.emit("deps", "info", ".hic format support not available (hicstraw missing - .mcool still works)")

# Legacy Keras/PyTorch/GB boundary classifiers belong to retired prediction and
# active-learning directions. They are not imported and have no server routes.
# The active workflow is five-method detection -> named evidence fields ->
# deterministic tier -> optional local review record.

# BigWig support is provided via bigwig_adapter (winbbi on Win, pyBigWig elsewhere).
BIGWIG_AVAILABLE = True

# =========================================================
# 1. BED 渲染引擎 (智能解析版)
# =========================================================
class BedTileRenderer:
    def __init__(self, bed_path, tile_size=1024,prediction_path=None):
        self.path = bed_path
        self.tile_size = tile_size
        self.prediction_map = {}
        try:
            # 1. 智能读取：sep=r'\s+' 兼容空格和 Tab，comment='#' 跳过注释
            # 先不指定 dtype，防止表头导致 int 转换崩溃
            self.df = pd.read_csv(
                bed_path,
                sep=r'\s+',
                header=None,
                comment='#',
                engine='python'
            )

            # 2. 自动识别列数并处理
            col_count = len(self.df.columns)
            if col_count >= 4:
                self.df = self.df.iloc[:, :4]
                self.df.columns = ['chrom', 'start', 'end', 'name']
            elif col_count == 3:
                self.df = self.df.iloc[:, :3]
                self.df.columns = ['chrom', 'start', 'end']
                self.df['name'] = ""  # 补齐列
            else:
                raise ValueError(f"BED 文件至少需要 3 列 (chrom, start, end)，检测到 {col_count} 列")

            # 3. 强制类型清洗 (关键：处理掉可能存在的 'start', 'end' 字符串表头)
            self.df['start'] = pd.to_numeric(self.df['start'], errors='coerce')
            self.df['end'] = pd.to_numeric(self.df['end'], errors='coerce')

            # 删除无法转换为数字的行（即表头行或脏数据）
            self.df = self.df.dropna(subset=['start', 'end'])

            self.df['start'] = self.df['start'].astype(int)
            self.df['end'] = self.df['end'].astype(int)
            self.df['chrom'] = self.df['chrom'].astype(str).str.strip().str.replace('^chr', '', case=False, regex=True)
            self.df = self.df.fillna("")

            # 4. 建立索引映射
            self.chrom_groups = {name: group for name, group in self.df.groupby('chrom')}

            # 这里的分辨率需要与你的 Hi-C mcool 保持一致或自定义层级
            self.resolutions = [1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000]

            print(f" BED 加载成功: {bed_path}, 解析列数: {col_count}, 总记录: {len(self.df)}")
            if prediction_path and os.path.exists(prediction_path):
                pred_df = pd.read_csv(prediction_path, sep=r'\s+', header=0)
                pred_df['chrom'] = pred_df['chrom'].astype(str).str.replace('^chr', '', regex=True)
                for _, row in pred_df.iterrows():
                    key = (str(row['chrom']), int(row['start']), int(row['end']))
                    self.prediction_map[key] = float(row['prediction'])
                print(f" 预测概率加载成功: {len(self.prediction_map)} 条")
        except Exception as e:
            print(f" BED 加载失败: {str(e)}")
            raise e

    def _normalize_chrom(self, chrom):
        return str(chrom).strip().replace('chr', '').replace('Chr', '')

    def get_all_chroms(self):
        chrom_info = []
        offset = 0
        # 按照 DataFrame 中的染色体顺序获取
        for name in self.df['chrom'].unique():
            group = self.chrom_groups[name]
            length = int(group['end'].max())
            chrom_info.append({"name": name, "length": length, "offset": offset})
            offset += length
        return chrom_info

    def get_info(self, chrom):
        chrom_key = self._normalize_chrom(chrom)
        if chrom_key not in self.chrom_groups:
            return {"error": f"chromosome not found: {chrom_key}",
                    "available": list(self.chrom_groups.keys())}
        group = self.chrom_groups[chrom_key]
        return {
            "chromosome": chrom_key,
            "max_zoom": len(self.resolutions) - 1,
            "available_resolutions": self.resolutions,
            "genome_length": int(group['end'].max()),
            "tile_size": self.tile_size
        }

    async def fetch_tile(self, chrom, zoom, x):
        chrom_key = self._normalize_chrom(chrom)
        if chrom_key not in self.chrom_groups: return []

        # 确保 zoom 不越界
        z = min(max(0, zoom), len(self.resolutions) - 1)
        res = self.resolutions[z]

        start_bp, end_bp = x * self.tile_size * res, (x + 1) * self.tile_size * res
        df = self.chrom_groups[chrom_key]

        # 范围查询优化
        subset = df[(df['end'] >= start_bp) & (df['start'] <= end_bp)].copy()

        # 抽样逻辑：防止低倍率下数据量过大导致前端卡死
        if zoom < 3 and len(subset) > 1000:
            subset = subset.sample(n=1000)

        return subset.to_dict(orient='records')


# =========================================================
# 2. HiC 渲染引擎
# =========================================================
# Two pre-read bounds on one range request. An unclamped genome-wide span at 25 kb reads
# ~450 M pixels and materialises tens of GB of int64 row/column arrays before any
# filtering — enough for the kernel to OOM-kill the backend and take every open session
# with it. _MAX_RANGE_BINS bounds each axis; _MAX_RANGE_CELLS bounds the RESULT, which is
# a Python dict per surviving pixel and is what actually dominates (a 4096 x 4096 window
# on the chr7 demo answered 335 MB of JSON). Both internal callers sit far below: tiles
# are 256 x 256 = 65 k cells and the colour-scale probe is at most 768 x 768 = 590 k.
_MAX_RANGE_BINS = 4096
_MAX_RANGE_CELLS = 2_000_000

# h5py slicing and the numpy reduction below are synchronous and CPU/IO bound. Running
# them on the event loop serialises the whole server behind one matrix read (measured:
# /api/health went from 0.3 ms to 1.1 s while a first /api/hic/info ran, and six tiles
# took 0.075/0.151/0.222/0.289/0.355/0.420 s — perfectly additive). A small dedicated
# pool keeps the loop responsive while bounding how many big reads run at once; h5py
# serialises HDF5 calls internally, so concurrent readers of one file are safe.
_MATRIX_POOL = ThreadPoolExecutor(max_workers=4, thread_name_prefix="hic-io")


class HiCRangeRenderer:
    def __init__(self, mcool_path):
        self.path = mcool_path
        self.f = h5py.File(mcool_path, "r")
        raw_res = [int(k) for k in self.f["resolutions"].keys() if k.isdigit()]
        self.resolutions = sorted(raw_res, reverse=True)
        self.data_map = self._load_metadata()

    def _load_metadata(self):
        data = {}
        for res in self.resolutions:
            g = self.f[f"resolutions/{res}"]
            chrom_names = g["chroms"]["name"][:].astype(str)
            chrom_lengths = g["chroms"]["length"][:]
            chrom_offsets = g["indexes"]["chrom_offset"][:]
            chrom_info = {
                name.replace("chr", ""): {
                    "length": int(length),
                    "start_bin": int(chrom_offsets[i]),
                    "n_bins": int(np.ceil(length / res))
                } for i, (name, length) in enumerate(zip(chrom_names, chrom_lengths))
            }
            data[res] = {
                "chrom_info": chrom_info,
                "pixels": g["pixels"],
                "bin_offsets": g["indexes"]["bin1_offset"][:],
                "weights": g["bins"]["weight"][:] if "weight" in g["bins"] else None
            }
        return data

    def get_all_chroms(self):
        res = self.resolutions[-1]
        info_map = self.data_map[res]["chrom_info"]
        chrom_info = []
        offset = 0
        for name, info in info_map.items():
            chrom_info.append({"name": name, "length": info["length"], "offset": offset})
            offset += info["length"]
        return chrom_info

    async def fetch_range_data(self, res, b1_range_abs, b2_range_abs):
        d = self.data_map.get(res)
        if not d: return [], res
        for lo, hi in (b1_range_abs, b2_range_abs):
            if hi - lo > _MAX_RANGE_BINS:
                raise ValueError(
                    f"range too large: {hi - lo} bins at {res} bp requested, "
                    f"maximum is {_MAX_RANGE_BINS} per axis")
        cells = max(0, b1_range_abs[1] - b1_range_abs[0]) * max(0, b2_range_abs[1] - b2_range_abs[0])
        if cells > _MAX_RANGE_CELLS:
            raise ValueError(
                f"range too large: {cells} bin pairs requested at {res} bp, maximum is "
                f"{_MAX_RANGE_CELLS}; request a smaller window or a coarser resolution")
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(
            _MATRIX_POOL, self._fetch_range_sync, res, b1_range_abs, b2_range_abs)

    def _fetch_range_sync(self, res, b1_range_abs, b2_range_abs):
        d = self.data_map.get(res)
        if not d: return [], res

        p = d["pixels"]
        offsets = d["bin_offsets"]
        w = d["weights"]

        total_bins = len(offsets) - 1
        abs_s1, abs_e1 = max(0, min(total_bins, b1_range_abs[0])), max(0, min(total_bins, b1_range_abs[1]))
        abs_s2, abs_e2 = max(0, min(total_bins, b2_range_abs[0])), max(0, min(total_bins, b2_range_abs[1]))

        if abs_s1 >= abs_e1: return [], res

        idx_start, idx_end = int(offsets[abs_s1]), int(offsets[abs_e1])
        if idx_start >= idx_end: return [], res

        bin2_ids_abs = p["bin2_id"][idx_start:idx_end]
        counts = p["count"][idx_start:idx_end]

        diffs = np.diff(offsets[abs_s1:abs_e1 + 1])
        rows_abs = np.repeat(np.arange(abs_s1, abs_e1), diffs)

        mask = (bin2_ids_abs >= abs_s2) & (bin2_ids_abs < abs_e2)
        if not np.any(mask): return [], res

        f_rows, f_cols, f_vals = rows_abs[mask], bin2_ids_abs[mask], counts[mask].astype(float)

        if w is not None:
            valid = (f_rows < len(w)) & (f_cols < len(w))
            f_rows, f_cols, f_vals = f_rows[valid], f_cols[valid], f_vals[valid]
            norm_vals = w[f_rows.astype(int)] * w[f_cols.astype(int)]
            f_vals = np.nan_to_num(f_vals * norm_vals)

        pos_mask = f_vals > 0
        f_rows, f_cols, f_vals = f_rows[pos_mask], f_cols[pos_mask], f_vals[pos_mask]

        results = [{"x": int(c * res), "y": int(r * res), "v": float(v)}
                   for r, c, v in zip(f_rows, f_cols, f_vals)]
        return results, res

    def close(self):
        self.f.close()

    def get_genome_index(self):
        # 1. 选一个分辨率（通常用最小分辨率获取 chrom 信息）
        res = self.resolutions[-1]
        d = self.data_map[res]

        chrom_info = []
        current_offset = 0

        # 2. 遍历 .mcool 内部存储的染色体顺序
        for name, info in d["chrom_info"].items():
            chrom_info.append({
                "name": name,
                "length": info["length"],
                "offset": current_offset
            })
            current_offset += info["length"]

        return {
            "total_length": current_offset,
            "chromosomes": chrom_info
        }


# =========================================================
# 2b. BigWig 渲染引擎
# Cross-platform: uses bigwig_adapter (winbbi on Win / pyBigWig elsewhere)
# =========================================================
from bigwig_adapter import BigWigReader as _BigWigReader

# Upper bound for /api/bigwig/signal?bins=. 4096 bins over the longest human
# chromosome is ~61 kb per bin — finer than any zoom the viewer draws — while the
# unbounded parameter let one request occupy the process for half a minute.
_MAX_SIGNAL_BINS = 4096

_GAPS_CACHE: dict[str, dict] = {}  # assembly → {chrom: [(s,e), ...]}


def _load_gaps(assembly: str = "hg19") -> dict:
    """Return {raw_chrom: [(start, end), ...]} of UCSC gap intervals (centromere,
    heterochromatin, contig, clone, telomere) for the requested assembly.

    Cached per-assembly after first call.  Both `chr1` and `1` keys are populated
    so callers using either chromosome-naming convention work.
    """
    if assembly in _GAPS_CACHE:
        return _GAPS_CACHE[assembly]
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "extra_mode", f"{assembly}_gaps.bed")
    result: dict = {}
    if not os.path.isfile(path):
        _GAPS_CACHE[assembly] = result
        return result
    try:
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom, s, e = parts[0], int(parts[1]), int(parts[2])
                result.setdefault(chrom, []).append((s, e))
                if chrom.startswith("chr"):
                    result.setdefault(chrom[3:], []).append((s, e))
        for k in result:
            result[k].sort()
    except Exception as exc:
        print(f"[gap-mask] failed to load {path}: {exc}")
    _GAPS_CACHE[assembly] = result
    return result


def _detect_assembly_from_chroms(chrom_to_size: dict) -> str:
    """hg38 chr1 = 248,956,422 bp; hg19 chr1 = 249,250,621 bp.  Difference is
    ~300 kb but reliable across all human references.  Falls back to hg19 if
    we can't tell."""
    chr1_len = chrom_to_size.get("chr1") or chrom_to_size.get("1")
    if chr1_len is None:
        return "hg19"
    return "hg38" if chr1_len < 249_000_000 else "hg19"


# Backward-compat shim for any call site still using the old name
def _load_hg19_gaps() -> dict:
    return _load_gaps("hg19")


class BigWigRenderer:
    def __init__(self, bw_path):
        self.path = bw_path
        self.chroms = {}
        self.name_map = {}  # clean_name -> raw_name from file
        # Read header eagerly with a short-lived reader (avoids long-held DLL on Win)
        try:
            with _BigWigReader() as reader:
                reader.open(bw_path)
                for raw in reader.get_chromosomes():
                    clean = raw.replace("chr", "").replace("Chr", "")
                    self.name_map[clean] = raw
                    self.chroms[raw] = int(reader.get_chrom_size(raw))
        except Exception as exc:
            print(f"[BigWigRenderer] header read failed: {exc}")

    def get_all_chroms(self):
        chrom_info = []
        offset = 0
        for raw, length in self.chroms.items():
            clean = raw.replace("chr", "").replace("Chr", "")
            chrom_info.append({"name": clean, "length": length, "offset": offset})
            offset += length
        return chrom_info

    def get_genome_index(self):
        cl = self.get_all_chroms()
        return {"total_length": sum(c["length"] for c in cl), "chromosomes": cl}

    def _resolve(self, chrom):
        if chrom in self.chroms:
            return chrom
        clean = chrom.replace("chr", "").replace("Chr", "")
        if clean in self.name_map:
            return self.name_map[clean]
        with_chr = f"chr{clean}"
        if with_chr in self.chroms:
            return with_chr
        return None

    # -- per-chromosome baseline cache (lazily populated) --
    #   Some ENCODE ChIP-seq bigWigs fill every bin with a constant baseline
    #   (fold-change-over-control or signal-p-value files) so "gap" regions
    #   show up as a flat line instead of blanks.  We detect that baseline
    #   once per chromosome and mask it out of the fetched signal.
    _baseline_cache: dict = None

    def _chromosome_baseline(self, raw_chrom: str) -> float | None:
        """Return an estimated baseline value for a chromosome (for masking).

        Some ENCODE ChIP-seq bigWigs (fold-change-over-control, signal p-value)
        fill gap regions with a constant non-zero baseline; classic raw-signal
        bigWigs leave gaps as ``None``.  We detect the first type by:

          1. Probing 2000 evenly-spaced bins across the chromosome.
          2. Computing the 25th percentile (p25).
          3. If p25 > 0.05 (meaningfully above zero), we treat p25 × 1.1 as
             the mask threshold — bins below it are gap / low-signal and get
             hidden.  This correctly keeps raw-signal bigWigs (p25 ≈ 0) intact
             while masking the 0.56 baseline in fold-change/p-value bigWigs.

        Cached per chromosome.
        """
        if self._baseline_cache is None:
            self._baseline_cache = {}
        if raw_chrom in self._baseline_cache:
            return self._baseline_cache[raw_chrom]

        clen = int(self.chroms.get(raw_chrom, 0))
        if clen <= 0:
            self._baseline_cache[raw_chrom] = None
            return None

        try:
            with _BigWigReader() as r:
                r.open(self.path)
                probes = r.read_zoom_signal(
                    chrom=raw_chrom, start=0, end=clen, num_bins=2000, use_closest=True,
                )
        except Exception:
            self._baseline_cache[raw_chrom] = None
            return None

        vals = [float(v) for v in (probes or []) if v is not None]
        if len(vals) < 200:
            self._baseline_cache[raw_chrom] = None
            return None

        vals_sorted = sorted(vals)
        n = len(vals_sorted)
        p25 = vals_sorted[n // 4]
        # Only mask if p25 is clearly above zero — indicates baseline-filled file.
        baseline = p25 if p25 > 0.05 else None
        self._baseline_cache[raw_chrom] = baseline
        return baseline

    async def fetch_signal(self, chrom, start, end, bins=500):
        raw = self._resolve(chrom)
        if raw is None:
            return []
        chrom_len = self.chroms[raw]
        start = max(0, min(int(start), chrom_len - 1))
        end = max(start + 1, min(int(end), chrom_len))
        bins = max(1, int(bins))

        bw_chrom_sizes: dict = {}
        def _read():
            try:
                with _BigWigReader() as reader:
                    reader.open(self.path)
                    try:
                        bw_chrom_sizes.update(reader.chroms())
                    except Exception:
                        pass
                    vals = reader.read_zoom_signal(
                        chrom=raw, start=start, end=end, num_bins=bins, use_closest=True,
                    )
                return [None if v is None else float(v) for v in vals]
            except Exception:
                return [None] * bins

        loop = asyncio.get_event_loop()
        vals = await loop.run_in_executor(_MATRIX_POOL, _read)

        # Authoritative gap/centromere mask from UCSC gap.txt + centromeres.txt.
        # Auto-detected assembly via chr1 length so hg19 and hg38 BigWigs both
        # mask their respective gap intervals.  Any bin whose CENTER falls
        # inside a centromere / heterochromatin / contig / clone / telomere
        # interval is dropped — no "fake" baseline values leak through.
        assembly = _detect_assembly_from_chroms(bw_chrom_sizes)
        gap_intervals = _load_gaps(assembly).get(raw, [])
        def _in_gap(pos: int) -> bool:
            # binary search since gap_intervals is sorted
            lo, hi = 0, len(gap_intervals)
            while lo < hi:
                mid = (lo + hi) // 2
                gs, ge = gap_intervals[mid]
                if pos < gs:
                    hi = mid
                elif pos >= ge:
                    lo = mid + 1
                else:
                    return True
            return False

        step = (end - start) / max(1, bins)

        def _mask():
            out = []
            for i, v in enumerate(vals):
                if v is None or v <= 0:
                    continue
                pos = int(start + i * step)
                if _in_gap(pos):
                    continue
                out.append({"pos": pos, "val": round(v, 4)})
            return out

        # Also off the loop: the read was already offloaded, but this per-bin binary
        # search over the gap track ran inline and is the part that scales with bins.
        return await loop.run_in_executor(_MATRIX_POOL, _mask)

    async def fetch_tile_data(self, chrom, start, end, bins=256):
        """Used by /api/v1/tiles/1d/* — returns flat list of values."""
        raw = self._resolve(chrom)
        if raw is None:
            return [0.0] * bins

        def _read():
            try:
                with _BigWigReader() as reader:
                    reader.open(self.path)
                    vals = reader.read_zoom_signal(
                        chrom=raw, start=int(start), end=int(end),
                        num_bins=int(bins), use_closest=True,
                    )
                return [0.0 if v is None else float(v) for v in vals]
            except Exception:
                return [0.0] * bins

        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(_MATRIX_POOL, _read)

    def close(self):
        # No persistent reader to close (opened per-request)
        pass


# =========================================================
# 3. 统一文件管理器
# =========================================================
# How many contact matrices may stay open at once. Each open .mcool costs one file
# descriptor plus ~100 MB of resident metadata (measured: registering a third full-genome
# matrix moved RSS 831 MB -> 939 MB), and nothing ever released them before process exit,
# so a long desktop session that imported ten matrices carried >1 GB it no longer used.
# Eviction is invisible: the token keeps working and its renderer is rebuilt on next use,
# so a view that is still on screen can never be broken by another tab's import.
_MAX_OPEN_MATRICES = 8


class GlobalManager:
    def __init__(self):
        self.renderers = {}
        self._spec = {}          # token -> (path, file_type, prediction_path), kept forever
        self._matrix_lru = []    # tokens of open matrices, least-recently-used first

    def _build(self, path, file_type, prediction_path=None):
        if file_type == "bed":
            return BedTileRenderer(path, prediction_path=prediction_path)
        if file_type == "bigwig" and BIGWIG_AVAILABLE:
            return BigWigRenderer(path)
        if file_type == "hic_juicer" and HIC_AVAILABLE:
            return HicFileReader(path)
        return HiCRangeRenderer(path)

    def _is_matrix(self, file_type):
        return file_type not in ("bed", "bigwig")

    def _touch(self, token, file_type):
        if not self._is_matrix(file_type):
            return
        if token in self._matrix_lru:
            self._matrix_lru.remove(token)
        self._matrix_lru.append(token)
        while len(self._matrix_lru) > _MAX_OPEN_MATRICES:
            victim = self._matrix_lru.pop(0)
            # Drop the reference rather than calling close(): a read already running in
            # the matrix pool still holds the renderer, so the HDF5 handle closes when
            # that read finishes instead of being pulled out from under it.
            if self.renderers.pop(victim, None) is not None:
                print(f"[manager] closing idle matrix session {victim} "
                      f"(over {_MAX_OPEN_MATRICES} open); it reopens on next use")

    def register(self, path, file_type, prediction_path=None):
        for token, spec in self._spec.items():
            if spec[0] == path:
                # For BED files, reload to pick up changes
                if file_type == "bed":
                    self.renderers[token] = BedTileRenderer(path, prediction_path=prediction_path)
                elif token not in self.renderers:
                    self.renderers[token] = self._build(*spec)
                self._touch(token, spec[1])
                return token
        token = str(uuid.uuid4())
        self._spec[token] = (path, file_type, prediction_path)
        self.renderers[token] = self._build(path, file_type, prediction_path)
        self._touch(token, file_type)
        return token

    def get(self, token):
        r = self.renderers.get(token)
        if r is not None:
            self._touch(token, self._spec.get(token, (None, "bed"))[1])
            return r
        spec = self._spec.get(token)
        if spec is None:
            return None
        try:
            r = self._build(*spec)
        except Exception as exc:
            print(f"[manager] reopening {spec[0]} for token {token} failed: {exc}")
            return None
        self.renderers[token] = r
        self._touch(token, spec[1])
        return r

    def close_all(self):
        for r in self.renderers.values():
            if hasattr(r, 'close'): r.close()

    def get_abs_info(self, token):
        renderer = self.get(token)
        if not renderer:
            return {}, 0, []

        # 情况 A: 如果是 Hi-C 渲染器 (基于 h5py，有 .f 属性)
        if hasattr(renderer, 'f'):
            res = renderer.resolutions[-1]
            # 访问 HDF5 内部结构
            chroms = renderer.f[f"resolutions/{res}/chroms"]
            names = chroms["name"][:].astype(str)
            lengths = chroms["length"][:]

            offsets = {}
            current_total = 0
            chrom_list = []
            for name, length in zip(names, lengths):
                clean_name = name.replace("chr", "").replace("Chr", "")
                offsets[clean_name] = current_total
                chrom_list.append({"name": clean_name, "length": int(length), "offset": current_total})
                current_total += int(length)
            return offsets, current_total, chrom_list

        # 情况 B: 如果是 BED 渲染器 (基于 Pandas，使用 get_all_chroms)
        else:
            chrom_list = renderer.get_all_chroms()
            # 这里的 chrom_list 已经是归一化后的 [{"name": "1", "offset": 0, ...}]
            offsets = {item['name']: item['offset'] for item in chrom_list}
            total_len = sum(item['length'] for item in chrom_list)
            return offsets, total_len, chrom_list

    def abs_to_rel(self, abs_pos, offsets):
        """将全局绝对位置转回 (染色体, 相对位置)"""
        last_chrom = None
        for name, offset in offsets.items():
            if abs_pos < offset:
                break
            last_chrom = name
        if last_chrom:
            return last_chrom, abs_pos - offsets[last_chrom]
        return None, 0


# =========================================================
# 4. API 路由处理函数
# =========================================================

async def handle_register(request):
    data = {}
    try:
        # 1. 强制使用 UTF-8 解码，防止 Windows 路径编码问题
        body = await request.read()
        data = json.loads(body.decode('utf-8'))
        if not isinstance(data, dict):
            return web.json_response({"error": "body must be a JSON object"}, status=400)

        path = data.get("path")
        if not path:
            return web.json_response({"error": "path is required"}, status=400)

        # 2. 强力清洗路径（后端双保险）
        path = path.replace('\u202a', '').replace('\u202b', '').strip()

        # friendly existence check (avoids surfacing a raw h5py/library traceback to the user)
        if not os.path.exists(path):
            return web.json_response(
                {"error": f"File not found: {path} — check the path is correct and readable by the backend."},
                status=404)

        f_type = data.get("type", "hic")
        lower_path = path.lower()
        if f_type == "hic" and lower_path.endswith('.hic'):
            f_type = "hic_juicer"
        elif lower_path.endswith(('.bigwig', '.bw')) or f_type == "bigwig":
            f_type = "bigwig"
        prediction_path = data.get("prediction_path", None)
        token = request.app['manager'].register(path, f_type, prediction_path=prediction_path)

        if f_type == "hic":
            # 3. 使用 try-except 包裹数据库操作，防止因为参数类型不匹配导致 500
            try:
                add_token_record(
                    token=token,
                    mcool_path=path,
                    resolution=data.get("resolution"),
                    is_threshold=data.get("is_threshold"),
                    overlap=data.get("overlap"),
                    hash_tag=data.get("hash_tag"),
                    fast_mode=data.get("fast_mode")
                )
            except Exception as db_err:
                print(f"Database error (ignored): {db_err}")

        return web.json_response({"token": token, "type": f_type})

    except json.JSONDecodeError as e:
        return web.json_response({"error": f"Invalid JSON format: {str(e)}"}, status=400)
    except Exception as e:
        # Path-not-found (404) and bad-JSON (400) are handled above; reaching here means the file
        # exists but could not be opened/parsed. Return a generic friendly message and log the raw
        # error SERVER-SIDE only — never leak the library traceback / server path to the client.
        # `data` is initialised to {} above: this handler used to call .get() on whatever
        # json.loads produced, so a JSON array made the except branch raise AttributeError
        # of its own and the friendly 400 below could never be returned.
        print(f"[register] open/parse failed for path={data.get('path', '?')!r}: {e}")
        return web.json_response(
            {"error": "Could not open this file as a valid mcool/.hic/bigWig — "
                      "check the file type and that it is not corrupted."},
            status=400)


async def get_genome_index(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)
    data = r.get_all_chroms()
    return web.json_response({"total_length": sum(i['length'] for i in data), "chromosomes": data})


# --- BED 路由 ---
async def bed_info(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)
    info = r.get_info(req.match_info['chrom'])
    # A missing chromosome used to come back as HTTP 200 carrying an error field, so a
    # client that branches on the status code read the failure as a success. Every other
    # not-found on this server is a 404.
    if isinstance(info, dict) and info.get("error"):
        return web.json_response(info, status=404)
    return web.json_response(info)


async def handle_bed_tile(req):
    token = req.match_info['token']

    # 1. 安全获取参数，防止索引越界
    try:
        z = int(req.match_info['z'])
        x = int(req.match_info['x'])
    except ValueError:
        return web.json_response({"error": "Invalid coordinates"}, status=400)

    r = req.app['manager'].get(token)
    if not r:
        return web.json_response({"error": "Invalid token"}, status=404)

    # 2. 获取基因组偏移量
    offsets, total_len, _ = req.app['manager'].get_abs_info(token)

    # 3. 边界检查：防止 z 超过 resolutions 数组长度导致 500
    safe_z = min(max(0, z), len(r.resolutions) - 1)
    res = r.resolutions[safe_z]
    tile_size = 256

    # 3b. CHROM-RELATIVE tiling. When ?chrom is supplied the tile is indexed within that
    # ONE chromosome and features are returned in chromosome-relative coordinates — matching
    # the Hi-C heatmap (handle_hic_tile is chrom-relative) and the TAD-VCI tier track. Without
    # this, a genome-wide BED is served in ABSOLUTE genome coordinates, so only chr1 (offset 0)
    # aligns and chr2..22 method tracks are misplaced/empty. (Single-chrom BEDs like the chr7
    # demo are unaffected: offset 0 → absolute == chrom-relative.)
    chrom_q = req.query.get('chrom')
    if chrom_q is not None and hasattr(r, 'chrom_groups'):
        ck = str(chrom_q).replace('chr', '').replace('Chr', '')
        rel_start = x * tile_size * res
        rel_end = (x + 1) * tile_size * res
        df = r.chrom_groups.get(ck)
        out = []
        if df is not None and not df.empty:
            subset = df[(df['end'] >= rel_start) & (df['start'] <= rel_end)].copy()
            subset['raw_chrom'] = ck
            subset['raw_start'] = subset['start']
            subset['raw_end'] = subset['end']
            out = subset.to_dict(orient='records')
        return web.json_response({"tile_id": f"{z}.{x}", "abs_range": [rel_start, rel_end],
                                  "resolution": res, "data": out})

    # 4. 计算该瓦片的全局范围
    abs_start = x * tile_size * res
    abs_end = (x + 1) * tile_size * res

    results = []

    # 5. 遍历染色体偏移表
    for mcool_chrom_name, offset in offsets.items():
        # --- 核心修复：归一化染色体名称，确保 "chr7" 和 "7" 能匹配 ---
        clean_name = str(mcool_chrom_name).replace('chr', '').replace('Chr', '')

        # 检查 BED 中是否有这条染色体的数据
        df = r.chrom_groups.get(clean_name)
        if df is None or df.empty:
            continue

        # 获取该染色体的最大长度（基于 BED 数据）
        chrom_len = int(df['end'].max())
        chrom_end_abs = offset + chrom_len

        # 6. 判断染色体是否在当前瓦片的绝对范围内
        if chrom_end_abs >= abs_start and offset <= abs_end:

            # 计算在该染色体内的相对查询范围
            rel_start = max(0, abs_start - offset)
            rel_end = abs_end - offset

            # 筛选重叠区间
            subset = df[(df['end'] >= rel_start) & (df['start'] <= rel_end)].copy()

            if not subset.empty:
                # 转换为绝对坐标供前端对齐 Hi-C 矩阵
                # handle_bed_tile 里 subset 转 dict 前加两列
                subset['raw_chrom'] = clean_name
                subset['raw_start'] = subset['start']  # 加偏移前
                subset['raw_end'] = subset['end']
                subset['start'] = subset['start'] + offset
                subset['end'] = subset['end'] + offset

                # 将该段数据加入结果集
                results.extend(subset.to_dict(orient='records'))

    # 7. 返回结果
    return web.json_response({
        "tile_id": f"{z}.{x}",
        "abs_range": [abs_start, abs_end],
        "resolution": res,
        "data": results
    })

async def bed_tiles(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)
    chrom, zoom, x = req.match_info['chrom'], int(req.match_info['zoom']), int(req.match_info['x'])
    data = await r.fetch_tile(chrom, zoom, x)
    return web.json_response({"data": data, "zoom": zoom, "x": x})


async def handle_mcool_chroms(request):
    """获取 mcool 文件中包含的染色体列表"""
    token = request.match_info['token']
    renderer = request.app['manager'].get(token)

    if not renderer or not hasattr(renderer, 'data_map'):
        return web.json_response({"error": "Invalid token or not a HiC file"}, status=404)

    # 从最低分辨率中获取染色体信息（通常所有分辨率的染色体是一致的）
    res = renderer.resolutions[-1]
    chrom_info = renderer.data_map[res]["chrom_info"]

    # 返回染色体名称列表，例如 ["1", "2", "X", ...]
    return web.json_response({
        "chromosomes": list(chrom_info.keys())
    })

# --- HiC 路由 ---
# Heatmap intensity = log1p(v * scale) / log1p(100). The curve denominator and the
# anchor value are fixed so EVERY matrix renders on one perceptual scale: a raw-count
# matrix (scale=1) puts its median contact (~3 counts) at ~0.30 intensity, and we map
# each balanced level onto that same curve by choosing scale so its median lands at 3.0.
_COLOR_ANCHOR_V = 3.0          # raw-count value GM12878's median sits at → ~0.30 intensity


async def _compute_color_scale(r):
    """Per-RESOLUTION intensity pre-scale for an mcool, cached on the renderer.

    Heatmap intensity is log1p(v * scale) / log1p(100). scale=1 reproduces the original
    raw-count rendering exactly (GM12878's raw levels stay pixel-identical to the
    known-good demo). The problem: ICE-BALANCED levels have contact values ~1e-3, where
    log1p(v) is in its linear regime, so a plain ceiling renders them almost white and
    only the top ~2% of contacts show — the matrix LOOKS sparse even though it is dense
    (4DN K562/IMR90/HepG2 have 92-99% non-zero pixels vs GM12878's 65%). We instead map
    each balanced level's MEDIAN contact onto the same intensity GM12878's raw median
    renders at (~0.30) by setting scale = 3.0 / median. That lifts the bulk of contacts
    into the visible gradient, so balanced matrices show full TAD structure like the raw
    demo instead of a sparse scatter.

    The same mcool can be raw at one resolution and balanced at another, so scale is
    computed PER resolution by sampling that level's near-diagonal block. Within a zoom
    every tile shares one scale → no per-tile discontinuity; changing zoom re-normalises,
    as HiGlass/Juicebox do. Cached.

    Two caveats on the probe, measured 2026-07-31 on the shipped GM12878 chr7 matrix:
      * It samples a fixed BIN COUNT (3*T = 768), so it covers a different genomic span at
        each level — 19.2 Mb at 25 kb but 7.68 Mb at 10 kb. Contact frequency decays with
        distance, so the wider block carries proportionally more long-range low pixels and
        its median is pulled down.
      * Over the SAME span the coarser level's median is genuinely higher (25 kb 3.624e-4
        vs 10 kb 1.711e-4 over 7.68 Mb) because a coarser bin aggregates more contacts.
    The two effects act in opposite directions and here largely cancel: as shipped the
    probe reports 1.627e-4 at 25 kb against 1.711e-4 at 10 kb, i.e. near-equal scales.
    That cancellation is a coincidence of this data, not a property of the design. Fixing
    the span would change every heatmap's appearance, so it is left as-is and recorded
    here rather than silently 'improved'.
    """
    cached = getattr(r, '_color_scale_arr', None)
    if cached is not None:
        return cached
    import numpy as _np
    out = []
    failed = []
    T = 256
    for res in r.resolutions:
        scale = 1.0  # raw-count default → identical to the original rendering
        try:
            d = r.data_map.get(res)
            first_chrom = next(iter(d['chrom_info']))
            ci = d['chrom_info'][first_chrom]
            off, nb = ci['start_bin'], ci['n_bins']
            n = min(nb, 3 * T)                       # near-diagonal block of the first chrom
            if n > 1:
                data, _ = await r.fetch_range_data(res, (off, off + n), (off, off + n))
                arr = _np.asarray([p['v'] for p in data if p['v'] and p['v'] > 0], dtype=float)
                if arr.size >= 50 and float(_np.percentile(arr, 99)) < 1.0:
                    med = float(_np.median(arr))      # balanced level → anchor median to ~0.30
                    if med > 0:
                        scale = _COLOR_ANCHOR_V / med
        except Exception as exc:
            # Falling back to scale=1.0 renders a BALANCED matrix almost blank. That
            # looked identical to "this region really is empty", so the failure has to
            # be visible: logged here and reported to the client as color_scale_ok.
            failed.append(res)
            print(f"[color-scale] {r.path}: resolution {res} probe failed "
                  f"({type(exc).__name__}: {exc}); falling back to raw-count scale 1.0")
        out.append(scale)
    r._color_scale_arr = out
    r._color_scale_failed = failed
    return out


def _chrom_safe_span(r, chrom, length):
    """Largest gap-free interval [lo, hi] of the chromosome, so the initial random view
    lands in a data-containing region instead of a telomere/centromere gap (where the
    heatmap is empty). Gaps come from the UCSC gap track; assembly is inferred from chr1's
    length so it works for any token. Falls back to [0, length] if gaps are unavailable."""
    cache = getattr(r, '_safe_span_cache', None)
    if cache is None:
        cache = {}; r._safe_span_cache = cache
    if chrom in cache:
        return cache[chrom]
    span = [0, int(length)]
    try:
        ci = r.data_map[r.resolutions[-1]]["chrom_info"]
        c1 = (ci.get("1") or ci.get("chr1") or {}).get("length")
        assembly = "hg19" if c1 == 249250621 else ("hg38" if c1 == 248956422 else "hg19")
        gaps = sorted(_load_gaps(assembly).get(str(chrom), []))
        # largest interval between consecutive gaps within [0, length]
        best_lo, best_hi, prev = 0, int(length), 0
        best = 0
        for gs, ge in list(gaps) + [(int(length), int(length))]:
            if gs - prev > best:
                best = gs - prev; best_lo, best_hi = prev, gs
            prev = max(prev, ge)
        if best > 0:
            span = [int(best_lo), int(best_hi)]
    except Exception as exc:
        print(f"[safe-span] {getattr(r, 'path', '?')} chr{chrom}: gap lookup failed "
              f"({type(exc).__name__}: {exc}); using the whole chromosome")
    cache[chrom] = span
    return span


async def hic_info(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)
    chrom = req.match_info['chrom'].replace("chr", "")
    res_key = r.resolutions[-1]
    info = r.data_map[res_key]["chrom_info"].get(chrom, {})
    scale_by_zoom = await _compute_color_scale(r)
    length = info.get("length", 0)
    return web.json_response({"length": length,
                              "available_resolutions": r.resolutions,
                              # intensity pre-scale per resolution level (index aligns with available_resolutions)
                              "color_scale_by_zoom": scale_by_zoom,
                              # False when at least one resolution fell back to the raw-count
                              # scale because its probe failed; a balanced matrix then renders
                              # near-blank, which is indistinguishable from "no contacts here".
                              "color_scale_ok": not getattr(r, "_color_scale_failed", []),
                              # largest gap-free interval -> initial random view avoids empty telomere/centromere
                              "safe_span": _chrom_safe_span(r, chrom, length)})


def _bp_window(start, end, limit):
    """Validate a base-pair window against ``limit`` (the length it must lie inside).

    Returns (start, end) or raises ValueError -> 400 via json_error_mw. Negative
    coordinates used to be clamped silently: on a genome-wide matrix a negative start
    still landed at a positive ABSOLUTE bin, so a request for chr7 was answered with
    contacts from the tail of chr6 and nothing said so.
    """
    if start < 0 or end < 0:
        raise ValueError("start and end must be >= 0")
    if end <= start:
        raise ValueError("end must be greater than start")
    if start >= limit:
        raise ValueError(f"start {start} is past the end of this sequence ({limit} bp)")
    return start, min(end, limit)


async def hic_range(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r: return web.json_response({"error": "Invalid token"}, status=404)

    chrom = req.match_info['chrom'].replace("chr", "")
    res = int(req.query.get("res", r.resolutions[0]))
    s1, e1 = int(req.query.get("start", 0)), int(req.query.get("end", 0))
    s2, e2 = int(req.query.get("start2", s1)), int(req.query.get("end2", e1))

    d = r.data_map.get(res)
    if not d or chrom not in d["chrom_info"]:
        return web.json_response({"error": "Chrom not found"}, status=404)

    # Clamp to THIS chromosome, not to the genome. fetch_range_data only clamped to the
    # total bin count, so an out-of-range end read straight past the chromosome into the
    # rest of the genome — the whole-genome case reads ~450 M pixels and OOM-kills the
    # process, taking every open session with it.
    chrom_len = int(d["chrom_info"][chrom]["length"])
    s1, e1 = _bp_window(s1, e1, chrom_len)
    s2, e2 = _bp_window(s2, e2, chrom_len)

    offset = d["chrom_info"][chrom]["start_bin"]
    data, actual_res = await r.fetch_range_data(
        res,
        (offset + s1 // res, offset + e1 // res),
        (offset + s2 // res, offset + e2 // res)
    )
    return web.json_response({"resolution": actual_res, "data": data})


async def hic_global_range(req):
    r = req.app['manager'].get(req.match_info.get('token'))
    if not r: return web.json_response({"error": "Invalid token"}, status=404)

    res = int(req.query.get("res", r.resolutions[0]))
    s1, e1 = int(req.query.get("start", 0)), int(req.query.get("end", 0))
    s2, e2 = int(req.query.get("start2", s1)), int(req.query.get("end2", e1))

    d = r.data_map.get(res)
    if not d:
        return web.json_response({"error": f"resolution {res} not in this matrix"}, status=404)
    genome_len = sum(int(ci["length"]) for ci in d["chrom_info"].values())
    s1, e1 = _bp_window(s1, e1, genome_len)
    s2, e2 = _bp_window(s2, e2, genome_len)

    data, actual_res = await r.fetch_range_data(res, (s1 // res, e1 // res), (s2 // res, e2 // res))
    return web.json_response({"resolution": actual_res, "data": data})


async def handle_hic_tile(req):
    r = req.app['manager'].get(req.match_info['token'])
    if not r or not hasattr(r, 'resolutions'):
        return web.json_response({"error": "Invalid token"}, status=404)
    try:
        z, x, y = int(req.match_info['z']), int(req.match_info['x']), int(req.match_info['y'])
    except (ValueError, TypeError):
        return web.json_response({"error": "Invalid tile coordinates"}, status=400)
    chrom = req.query.get('chrom', '').replace('chr', '').replace('Chr', '')

    tile_size = 256
    if z < 0 or z >= len(r.resolutions):
        return web.json_response({"error": "zoom out of range"}, status=400)
    res = r.resolutions[z]

    # Chromosome bin offset (for chrom-relative coords) + n_bins (to CLAMP the tile to this
    # chromosome — without the clamp a coarse-zoom tile spans many chromosomes and bleeds
    # neighbouring-chromosome contacts into the displayed chromosome).
    d = r.data_map.get(res)
    chrom_offset = 0
    chrom_nbins = None
    if chrom and d and chrom in d['chrom_info']:
        chrom_offset = d['chrom_info'][chrom]['start_bin']
        chrom_nbins = d['chrom_info'][chrom]['n_bins']

    x0, x1 = x * tile_size, (x + 1) * tile_size
    y0, y1 = y * tile_size, (y + 1) * tile_size
    if chrom_nbins is not None:
        x1 = min(x1, chrom_nbins)
        y1 = min(y1, chrom_nbins)
    if x1 <= x0 or y1 <= y0:
        return web.json_response({"tile_id": f"{z}.{x}.{y}", "data": []})

    data, actual_res = await r.fetch_range_data(
        res,
        (chrom_offset + x0, chrom_offset + x1),
        (chrom_offset + y0, chrom_offset + y1),
    )

    # Convert to chromosome-relative coords + drop anything beyond the chromosome (defensive).
    #
    # `v` is serialised to 6 significant digits. It is a float64 whose full 17-digit decimal
    # expansion was ~60% of the wire size of a tile (the initial view of a 10 kb package
    # fetches 60 tiles / 38 MB), and the client's only use of it is
    #   ci = int( min(1, log1p(v*scale)/log1p(100)) * 255 )
    # i.e. one of 256 colour indices. The colour index moves by at most
    #   255 * (v*scale/(1+v*scale)) * 5e-6 / log1p(100)  <  3e-4
    # for any v, so 6 digits cannot change a rendered pixel. Verified by pixel-diffing the
    # heatmap before and after (0 differing pixels), not only by that bound.
    offset_bp = chrom_offset * res
    chrom_len_bp = chrom_nbins * res if chrom_nbins is not None else None
    out = []
    for p in data:
        p['x'] -= offset_bp
        p['y'] -= offset_bp
        if chrom_len_bp is None or (0 <= p['x'] < chrom_len_bp and 0 <= p['y'] < chrom_len_bp):
            v = p.get('v')
            if isinstance(v, float):
                p['v'] = float(f"{v:.6g}")
            out.append(p)

    return web.json_response({"tile_id": f"{z}.{x}.{y}", "data": out})


async def bed_get_prediction(req):
    """
    GET /api/bed/prediction/{token}?chrom=7&start=65000&end=75000
    返回该区间的预测概率
    """
    token = req.match_info['token']
    r = req.app['manager'].get(token)
    if not r:
        return web.json_response({"error": "Invalid token"}, status=404)

    chrom = req.query.get('chrom', '').replace('chr', '')
    try:
        start = int(req.query.get('start', 0))
        end   = int(req.query.get('end', 0))
    except ValueError:
        return web.json_response({"error": "Invalid coordinates"}, status=400)

    prob = r.prediction_map.get((chrom, start, end), None)
    return web.json_response({
        "chrom": chrom,
        "start": start,
        "end":   end,
        "prediction": prob  # None 表示没有预测值
    })

async def handle_tad_run(request):
    """
    核心修改：只对最终重合的 TAD_final 文件进行预测
    """
    try:
        data = await request.json()
    except Exception:
        return web.json_response({"error": "Invalid JSON body"}, status=400)
    if not isinstance(data, dict):
        return web.json_response({"error": "body must be a JSON object"}, status=400)
    token = data.get('token')
    if not token:
        return web.json_response({"error": "Token is required"}, status=400)

    renderer = request.app['manager'].get(token)
    if not renderer or not hasattr(renderer, 'path'):
        return web.json_response({"error": "Invalid token or not a HiC file"}, status=400)

    try:
        method_set = resolve_method_set(data)
    except MethodSetContractError as exc:
        return web.json_response({"error": str(exc)}, status=exc.status_code)

    # 标记任务开始
    mark_running(token)

    # 获取参数
    record = get_record(token) or {}
    resolution = data.get('resolution', record.get('resolution', 10000))
    is_threshold = data.get('is_threshold', record.get('is_threshold', 0.3))
    overlap = data.get('overlap', record.get('overlap', 0.7))
    hash_tag = data.get('hash_tag', record.get('hash_tag', token[:8]))
    fast_mode = data.get('fast_mode', record.get('fast_mode', True))
    output_dir = data.get('output_dir', './tad_results')
    # method_set is validated before task state is mutated. The active default is
    # tadvci; v2 remains an explicit compatibility path and legacy is retired.

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    selected_chroms = data.get('selected_chroms', [])  # 例如 ["1", "7"]
    try:
        # 1. 运行 TAD 检测。method_set:
        #    'tadvci' (5-method built-in panel + per-import TAD-VCI evidence tier) — the
        #             unified path: 5 callers on the user's mcool -> votes (BD1) + optional
        #             CTCF/RAD21 (BD2/BD3) -> D5..D1 tier annotation, no "Disputed" track.
        #    'v2' is the explicit clean-room 4-method compatibility path.
        if method_set == 'tadvci':
            ctcf_bw = data.get('ctcf_bw') or (data.get('bigwig_paths') or [None])[0]
            rad21_bw = data.get('rad21_bw') or (data.get('bigwig_paths') or [None, None])[1] \
                if isinstance(data.get('bigwig_paths'), list) and len(data.get('bigwig_paths', [])) > 1 else data.get('rad21_bw')
            user_bed_paths = data.get('user_bed_paths') or []
            # BD4 only when the client supplies loops_path. Never auto-mount demo
            # HiCCUPS / ChIP files on user imports (demos come from use_multicell_example).
            loops_path = data.get('loops_path') or None
            output_paths = await asyncio.get_event_loop().run_in_executor(
                None,
                lambda: detect_TAD_boundaries_tadvci(
                    mcool_file=renderer.path, selected_chroms=selected_chroms,
                    resolution=resolution, hash_tag=hash_tag, output_dir=output_dir,
                    ctcf_bw=ctcf_bw, rad21_bw=rad21_bw, user_bed_paths=user_bed_paths,
                    loops_path=loops_path, cell_line=data.get('cell_line'),
                    assembly=data.get('assembly'),
                    # 'contribute' = uploaded BED joins the consensus as an extra voter;
                    # 'evaluate' (default) = uploaded BED is scored + flagged but does not vote.
                    upload_mode=(data.get('upload_mode') or 'evaluate'))
            )
            # register the per-import tier annotation so the visualization + sign-out read THIS data
            try:
                from tad_vci.api_routes import set_imported_annotation
                set_imported_annotation(token, output_paths.get('TADVCI_annotation'))
            except Exception as reg_err:
                print(f"tadvci annotation registration failed (non-fatal): {reg_err}")
        else:
            detector = detect_TAD_boundaries_v2
            output_paths = await asyncio.get_event_loop().run_in_executor(
                None,
                lambda: detector(
                    mcool_file=renderer.path,
                    selected_chroms=selected_chroms,
                    resolution=resolution,
                    is_threshold=is_threshold,
                    overlap=overlap,
                    hash_tag=hash_tag,
                    fast_mode=fast_mode,
                    output_dir=output_dir
                )
            )

        # The retired Keras boundary-prediction branch used to live here. It called a
        # module named `test` that is neither imported nor defined anywhere in the tree,
        # so it was a NameError guarded only by the hardcoded DACTOR_AVAILABLE = False —
        # flipping that flag would have crashed every detection run.
        prediction_output = None

        update_tad_outputs(
            token=token,
            output_paths=output_paths,
            prediction_output=prediction_output
        )

        return web.json_response({
            "status": "done",
            "output_paths": output_paths,
            "prediction_output": prediction_output,
            "message": "TAD detection completed (only final overlapping TADs predicted)"
        })

    except Exception as e:
        error_msg = str(e)
        mark_failed(token, error_msg)
        return web.json_response({
            "status": "failed",
            "error": error_msg
        }, status=500)

async def handle_tad_status(request):
    """
    GET /api/tad/status/{token}
    查询 TAD 任务状态及所有入库字段。
    """
    token  = request.match_info['token']
    record = get_record(token)
    if not record:
        return web.json_response({"error": "Not found"}, status=404)
    return web.json_response(record)


# =========================================================
# 5. 服务启动
# =========================================================
async def on_shutdown(app):
    app['manager'].close_all()
    _MATRIX_POOL.shutdown(wait=False, cancel_futures=True)


# Subsystems whose absence leaves the process listening but the product unusable.
# Registration of the TAD-VCI route table can genuinely fail (a packaging layout where
# api_routes' parents[3] walks off the drive root is enough) and the failure was only a
# startup warn line, after which the banner still said "All systems operational" and
# /api/health still answered a flat 200 — so the launcher reported a healthy start for a
# backend on which every /api/tadvci/* call 404s.
DEGRADED: dict[str, str] = {}


# --- 在 handle_functions 区域添加 ---
async def handle_health(request):
    """用于 Electron 探测后端是否启动完成"""
    ok = not DEGRADED
    return web.json_response(
        {"status": "ok" if ok else "degraded", "service": "creditad-api",
         "tadvci": "tadvci" not in DEGRADED,
         "errors": [f"{k}: {v}" for k, v in sorted(DEGRADED.items())]},
        status=200 if ok else 503,
    )

# --- 查询所有 token ---
async def handle_get_all_records(request):
    """
    GET /api/records
    返回数据库中所有 token 记录。
    """
    from sql import get_all_records
    records = get_all_records()
    return web.json_response({"records": records})


async def handle_genome_index(request):
    token = request.match_info['token']
    manager = request.app['manager']
    renderer = manager.get(token)

    if not renderer:
        return web.json_response({"error": "Invalid token"}, status=404)

    # 调用上面的解析逻辑
    index_data = renderer.get_genome_index()
    return web.json_response(index_data)




# --- 删除单条记录 ---
async def handle_delete_record(request):
    """
    DELETE /api/records/{token}
    删除指定 token 的记录。
    """
    from sql import delete_record

    token = request.match_info['token']
    success = delete_record(token)
    if success:
        return web.json_response({"status": "deleted", "token": token})
    else:
        return web.json_response({"error": "Token not found"}, status=404)


async def bed_global_info(req):
    """
    GET /api/bed/global-info/{token}
    返回全局 BED 信息，不需要染色体参数
    供 TAD BED（绝对坐标）使用
    """
    token = req.match_info['token']
    r = req.app['manager'].get(token)
    if not r:
        return web.json_response({"error": "Invalid token"}, status=404)

    # 获取所有染色体及总长度
    _, total_len, chrom_list = req.app['manager'].get_abs_info(token)

    return web.json_response({
        "genome_length": total_len,
        "available_resolutions": r.resolutions,
        "chromosomes": chrom_list
    })


# --- 删除所有记录 ---
async def handle_delete_all_records(request):
    """
    DELETE /api/records
    删除数据库中所有 token 记录。
    """
    from sql import delete_all_records

    count = delete_all_records()
    return web.json_response({"status": "all_deleted", "deleted_count": count})

_labels_db_path = os.path.join(os.path.dirname(__file__), 'extra_mode', 'labels.json')

def _load_labels():
    if os.path.exists(_labels_db_path):
        try:
            with open(_labels_db_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception:
            pass
    return []

def _save_labels(labels):
    os.makedirs(os.path.dirname(_labels_db_path), exist_ok=True)
    with open(_labels_db_path, 'w', encoding='utf-8') as f:
        json.dump(labels, f, ensure_ascii=False, indent=2)

_labels_db = _load_labels()

async def handle_labels_save(request):
    global _labels_db
    try:
        data = await request.json()
    except Exception:
        return web.json_response({"error": "Invalid JSON body"}, status=400)
    if not isinstance(data, dict) or not data.get('id'):
        return web.json_response({"error": "label must be a JSON object with an 'id'"}, status=400)
    idx = next((i for i, l in enumerate(_labels_db) if l.get('id') == data.get('id')), None)
    if idx is not None:
        _labels_db[idx] = data
    else:
        _labels_db.append(data)
    _save_labels(_labels_db)
    return web.json_response({"status": "saved"})

async def handle_labels_list(request):
    return web.json_response({"labels": _labels_db})

async def handle_labels_delete(request):
    global _labels_db
    label_id = request.match_info['id']
    _labels_db = [l for l in _labels_db if l.get('id') != label_id]
    _save_labels(_labels_db)
    return web.json_response({"status": "deleted"})


async def handle_bigwig_signal(request):
    token = request.match_info['token']
    chrom = request.match_info['chrom']
    r = request.app['manager'].get(token)
    if not r or not isinstance(r, BigWigRenderer):
        return web.json_response({"error": "Invalid bigwig token"}, status=404)
    start = int(request.query.get('start', 0))
    end = int(request.query.get('end', 10000000))
    bins = int(request.query.get('bins', 500))
    raw = r._resolve(chrom)
    if raw is None:
        return web.json_response({"error": f"chromosome not found: {chrom}",
                                  "available": sorted(r.name_map)}, status=404)
    # Reject instead of silently clamping: start=-100&end=-1 used to answer HTTP 200 and
    # echo the impossible window back, so a caller could not tell a bad query from an
    # empty region. bins was unbounded — bins=200000 took 31.6 s and returned 6.5 MB
    # (bins=500 takes 0.07 s), and the client cancelling does not stop the work.
    if start < 0 or end < 0:
        return web.json_response({"error": "start and end must be >= 0"}, status=400)
    if end <= start:
        return web.json_response({"error": "end must be greater than start"}, status=400)
    if bins < 1 or bins > _MAX_SIGNAL_BINS:
        return web.json_response(
            {"error": f"bins must be between 1 and {_MAX_SIGNAL_BINS}"}, status=400)
    data = await r.fetch_signal(chrom, start, end, bins)
    return web.json_response({"chrom": chrom, "start": start, "end": end, "data": data})


def main():
    with startup_log.step("db", "Initializing SQLite tables (idempotent)"):
        init_db()

    with startup_log.step("bigwig", "Validating canonical BigWigs (file integrity)"):
        try:
            from validate_bigwig import validate_bigwig as _vb
            # Validate the tracks the product actually grades against: the per-cell hg38
            # fold-change-over-control pair used by the six demo packages. The old hg19
            # copies under example_data/ and extra_mode/ are legacy and unread; scanning
            # them made startup look healthy while saying nothing about the live path.
            # Follow the same root the resolver and the Data Manager guide use, so a
            # machine with CREDITAD_SHARED_RAW set validates the files it will actually
            # read — and a packaged build does not print a developer path six times.
            from data_manager import shared_raw_root
            _shared = str(shared_raw_root())
            checks = {
                "GM12878 CTCF (hg38 fc)":  f"{_shared}/gm12878/ctcf/ENCFF734CUT.bigWig",
                "GM12878 RAD21 (hg38 fc)": f"{_shared}/gm12878/rad21/ENCFF571ZJJ.bigWig",
                "IMR90 CTCF (hg38 fc)":    f"{_shared}/imr90/ctcf/ENCFF105FHL.bigWig",
                "IMR90 RAD21 (hg38 fc)":   f"{_shared}/imr90/rad21/ENCFF048PZI.bigWig",
                "HepG2 CTCF (hg38 fc)":    f"{_shared}/hepg2/ctcf/ENCFF357NFO.bigWig",
                "HepG2 RAD21 (hg38 fc)":   f"{_shared}/hepg2/rad21/ENCFF972ODZ.bigWig",
            }
            for label, path in checks.items():
                if not os.path.exists(path):
                    startup_log.emit("bigwig", "info", f"{label}: missing — skipped")
                    continue
                r = _vb(path)
                if r.verdict == "ok":
                    startup_log.emit("bigwig", "info",
                                     f"{label}: max={r.max_val:.0f} sumData={r.sum_data:.1e} OK")
                elif r.verdict == "empty":
                    startup_log.emit("bigwig", "warn",
                                     f"{label}: file looks empty — {r.detail}")
                else:  # error
                    startup_log.emit("bigwig", "warn",
                                     f"{label}: {r.detail}")
        except Exception as e:
            startup_log.emit("bigwig", "warn",
                             f"BigWig validation skipped ({type(e).__name__}: {e})")

    @web.middleware
    async def json_error_mw(request, handler):
        # Never let an endpoint emit a plain-text/HTML 500 — the frontend does response.json()
        # on tile/range responses, so an HTML 500 throws on the client and blanks the view.
        # Bad input (e.g. non-int query params) → 400 JSON; anything else → 500 JSON.
        try:
            return await handler(request)
        except web.HTTPException:
            raise
        except (ValueError, TypeError) as e:
            return web.json_response({"error": f"Bad request: {e}"}, status=400)
        except Exception as e:
            print(f"[api-error] {request.method} {request.path}: {e!r}")
            return web.json_response({"error": "Internal error processing request"}, status=500)

    app = web.Application(client_max_size=1024 ** 2 * 100, middlewares=[json_error_mw])
    app['manager'] = GlobalManager()
    app.on_shutdown.append(on_shutdown)

    # Register data download manager routes
    try:
        from data_manager import register_data_routes
        register_data_routes(app)
        startup_log.emit("deps", "info", "Data Manager API registered (/api/data/*)")
    except ImportError:
        startup_log.emit("deps", "warn", "data_manager.py not found - download feature unavailable")

    app.add_routes([
        web.get("/api/health", handle_health),
        web.post("/api/register",              handle_register),
        web.get('/api/genome/index/{token}',   get_genome_index),
        web.get("/api/records", handle_get_all_records),
        web.get('/api/bed/info/{token}/{chrom}',          bed_info),
        web.get('/api/bed/tiles/{token}/{chrom}/{zoom}/{x}', bed_tiles),
        web.delete("/api/records/{token}", handle_delete_record),
        web.delete("/api/records", handle_delete_all_records),
        web.get('/api/v1/tiles/1d/{token}/{z}/{x}', handle_bed_tile),
        web.get('/api/bed/global-info/{token}', bed_global_info),
        web.get('/api/bed/prediction/{token}', bed_get_prediction),
        web.get("/api/mcool/chroms/{token}", handle_mcool_chroms),
        web.get("/api/hic/info/{token}/{chrom}",  hic_info),
        web.get("/api/hic/range/{token}/{chrom}", hic_range),
        web.get("/api/hic/global-range/{token}",  hic_global_range),
        web.get('/api/v1/tiles/2d/{token}/{z}/{x}/{y}', handle_hic_tile),
        web.post("/api/tad/run",               handle_tad_run),
        web.get("/api/tad/status/{token}",     handle_tad_status),
        web.post("/api/labels/save", handle_labels_save),
        web.get("/api/labels/list", handle_labels_list),
        web.delete("/api/labels/{id}", handle_labels_delete),
        web.get("/api/bigwig/signal/{token}/{chrom}", handle_bigwig_signal),
    ])

    # TAD-VCI graded-evidence curation routes (/api/tadvci/*)
    try:
        from tad_vci.api_routes import register_tadvci_routes
        register_tadvci_routes(app)
        startup_log.emit("deps", "info", "TAD-VCI curation API registered (/api/tadvci/*)")
    except Exception as _e:  # noqa: BLE001
        DEGRADED["tadvci"] = f"{type(_e).__name__}: {_e}"
        startup_log.emit("deps", "error", f"TAD-VCI API not registered: {_e}")

    cors = aiohttp_cors.setup(app, defaults={
        "*": aiohttp_cors.ResourceOptions(allow_headers="*", allow_methods="*")
    })
    for route in list(app.router.routes()):
        cors.add(route)

    # Bind ONCE, here, and hand the bound socket to run_app. The old code probed with a
    # throwaway socket, closed it, and let run_app bind again — two bugs in one:
    #
    #  * the probe set SO_REUSEADDR, which on WINDOWS means "bind even if another socket
    #    already holds this port" (the opposite of the POSIX meaning). The guard was
    #    therefore a no-op on the only platform that ships this code to end users;
    #  * "listening"/"All systems operational" were emitted from an on_startup hook, which
    #    aiohttp runs BEFORE the site binds. A bind failure after that produced a log that
    #    reported a fully healthy backend and a process that then exited — exactly what a
    #    user reported on 2026-07-31 (full startup log, then "exited with code 2").
    #
    # Binding once removes the race and makes the banner true: if this raises, nothing has
    # claimed to be ready.
    sock = _reserve_port(SERVER_PORT)
    startup_log.emit("port", "done", f"TCP port {SERVER_PORT} bound on 127.0.0.1")

    async def _announce_listen(app_):
        startup_log.emit("listen", "done", f"HTTP server listening on http://127.0.0.1:{SERVER_PORT}")
        startup_log.ready(SERVER_PORT, DEGRADED)

    app.on_startup.append(_announce_listen)
    print(f" CrediTAD API ready [Port: {SERVER_PORT}]")
    web.run_app(app, sock=sock, print=lambda *_: None)


def _reserve_port(port: int):
    """Return a socket bound to ``port``, or raise with an actionable message.

    The returned socket is handed straight to ``web.run_app(sock=…)``, so the port is
    held continuously from this check to serving — there is no window in which another
    process can take it, and no second bind that can fail after the app has announced
    itself as ready.

    The previous implementation ran ``lsof -ti:<port> | xargs -r kill -9`` (and the
    taskkill equivalent on Windows), which terminates ANY process on that port —
    including an unrelated service belonging to another user on a shared host. A tool
    that grades other people's data must not do that. We probe, and on conflict we exit
    with an actionable message naming the port and the override.
    """
    import socket

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        # NEVER SO_REUSEADDR on Windows: there it permits binding a port another socket
        # already holds, so the conflict this function exists to detect goes unnoticed.
        # SO_EXCLUSIVEADDRUSE is the Windows way to ask for a genuinely exclusive bind.
        if hasattr(socket, "SO_EXCLUSIVEADDRUSE"):
            s.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        elif os.name == "posix":
            # POSIX: this only relaxes TIME_WAIT reuse, which is what a restarting server
            # wants; it does not let two live listeners share the port.
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        s.bind(("127.0.0.1", port))
        s.listen(128)
        return s
    except OSError:
        s.close()
        raise OSError(
            f"Port {port} on 127.0.0.1 is already in use. CrediTAD will not terminate the "
            f"process holding it. Stop that process yourself (another CrediTAD window, or "
            f"an installed copy running alongside this one), or start on another port with "
            f"CREDITAD_PORT=<port> (currently {port})."
        ) from None


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        startup_log.fatal("listen", f"Backend crashed during startup: {type(exc).__name__}: {exc}")
        raise
