#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TAD边界检测综合脚本 - 集成 DI (Dixon 2012) + TopDom (Shin 2016) + IS (Crane 2015) + cooltools (Venev 2024) 四种方法
纯 h5py + numpy + pandas 实现，无需 cooler / pybedtools / bedtools，Windows 可用。
用法: python tad_detection.py -i <mcool文件> -r <分辨率> [选项]
"""

import argparse
import os
import sys
import time

import numpy as np
import h5py
import pandas as pd
from tqdm import tqdm
from hmmlearn import hmm

try:
    from hic_reader import HicFileReader, HIC_AVAILABLE
except ImportError:
    HIC_AVAILABLE = False


# ═══════════════════════════════════════════════
# .hic 读取适配层
# ═══════════════════════════════════════════════

class HicReaderAdapter:
    def __init__(self, hic_path: str, resolution: int, chroms: list = None):
        if not HIC_AVAILABLE:
            raise ImportError("hicstraw is required for .hic files. Install with: pip install hic-straw")
        self._hic = HicFileReader(hic_path)
        available_res = self._hic.resolutions
        self.resolution = min(available_res, key=lambda x: abs(x - resolution))
        all_chroms = list(self._hic._chrom_info.keys())
        if chroms:
            self._target_chroms = [c for c in chroms if c in all_chroms]
        else:
            self._target_chroms = [c for c in all_chroms if c not in ('M', 'MT')]
        print(f" HicReaderAdapter | resolution: {self.resolution} | chromosomes: {len(self._target_chroms)}")

    @property
    def chromnames(self) -> list:
        return self._target_chroms

    def bins_for(self, chrom: str) -> pd.DataFrame:
        info = self._hic._chrom_info.get(chrom)
        if not info:
            return pd.DataFrame(columns=["chrom", "start", "end"])
        length = info["length"]
        starts = np.arange(0, length, self.resolution)
        ends = np.minimum(starts + self.resolution, length)
        return pd.DataFrame({"chrom": [chrom] * len(starts), "start": starts, "end": ends})

    def matrix(self, chrom: str) -> np.ndarray:
        print(f"  Extracting matrix: {chrom} ...")
        t0 = time.time()
        info = self._hic._chrom_info.get(chrom)
        if not info:
            raise ValueError(f"Chromosome {chrom} not found")
        length = info["length"]
        mat = self._hic.get_contact_matrix(chrom, 0, length, self.resolution)
        n = mat.shape[0]
        print(f"  Done | size: {n}x{n} | time: {time.time() - t0:.2f}s")
        return mat


# ═══════════════════════════════════════════════
# mcool 读取层（替代 cooler）
# ═══════════════════════════════════════════════

import h5py
import pandas as pd
import numpy as np
import time


class McoolReader:
    """高性能 .mcool 读取器：支持染色体筛选与切片读取。"""

    _DIVISIVE_BALANCE_NAMES = {"KR", "VC", "VC_SQRT"}

    def __init__(self, cool_file: str, resolution: int, chroms: list = None,
                 balance_name: str = "weight", strict_balance: bool = False,
                 balance_weights_path: str | None = None):
        self.cool_file = cool_file
        self.resolution = self._resolve_resolution(resolution)
        self._root_path = f"resolutions/{self.resolution}"
        self.balance_name = str(balance_name) if balance_name is not None else None
        self.strict_balance = bool(strict_balance)
        self.balance_weights_path = str(balance_weights_path) if balance_weights_path else None

        # 1. 加载基础元数据
        self._bins_df = self._load_bins()
        self._external_balance_weights = self._load_external_balance_weights()

        # 2. 确定目标染色体列表（并排序）
        all_available = self._bins_df["chrom"].unique().tolist()
        if chroms:
            # 兼容两种命名风格：传入 '7' 或 'chr7' 都能匹配 mcool 中
            # 的 'chr7' / '7'。先按 "去掉 chr 前缀" 归一化到 ``key``，
            # 然后从 all_available 里挑出 key 相同的真实名字。
            def _key(c):
                s = str(c)
                return s[3:] if s.lower().startswith('chr') else s
            wanted_keys = {_key(c) for c in chroms}
            self._target_chroms = sorted(
                [c for c in all_available if _key(c) in wanted_keys],
                key=lambda c: self._chrom_sort_key(c)
            )
        else:
            self._target_chroms = sorted(
                all_available,
                key=lambda c: self._chrom_sort_key(c)
            )

        print(f" McoolReader 初始化完成 | 分辨率: {self.resolution} | 目标染色体: {len(self._target_chroms)} 条")

    # ── 内部辅助 ──────────────────────────────

    @staticmethod
    def _chrom_sort_key(name: str):
        n = str(name).replace("chr", "")
        if n.isdigit(): return (0, int(n))
        return (1, {"X": 100, "Y": 101, "M": 102, "MT": 102}.get(n.upper(), 200))

    def _resolve_resolution(self, resolution: int) -> int:
        with h5py.File(self.cool_file, "r") as f:
            available = [int(r) for r in f["resolutions"].keys() if r.isdigit()]
        if resolution in available: return resolution
        closest = min(available, key=lambda x: abs(x - resolution))
        return closest

    def _load_bins(self) -> pd.DataFrame:
        """加载 bins 信息并预计算每个染色体的 bin 范围。"""
        with h5py.File(self.cool_file, "r") as f:
            grp = f[self._root_path]
            chrom_ids = grp["bins"]["chrom"][:]
            starts = grp["bins"]["start"][:]
            ends = grp["bins"]["end"][:]
            chrom_names = [c.decode() if isinstance(c, bytes) else c for c in grp["chroms"]["name"][:]]

        chrom_col = [chrom_names[i] for i in chrom_ids]
        return pd.DataFrame({"chrom": chrom_col, "start": starts, "end": ends})

    def _load_external_balance_weights(self) -> np.ndarray | None:
        """Load and validate a manifest-bound multiplicative ICE sidecar."""
        if not self.balance_weights_path:
            return None
        path = self.balance_weights_path
        if not os.path.isfile(path):
            raise ValueError(f"balance sidecar does not exist: {path}")
        try:
            with np.load(path, allow_pickle=False) as sidecar:
                required = {
                    "weights",
                    "resolution_bp",
                    "balance_name",
                    "balance_semantics",
                }
                missing = sorted(required.difference(sidecar.files))
                if missing:
                    raise ValueError(f"missing arrays {missing}")
                weights = np.asarray(sidecar["weights"], dtype=float)
                resolution = int(np.asarray(sidecar["resolution_bp"]).item())
                balance_name = str(np.asarray(sidecar["balance_name"]).item())
                semantics = str(np.asarray(sidecar["balance_semantics"]).item())
        except (OSError, ValueError, KeyError) as exc:
            raise ValueError(f"invalid balance sidecar {path}: {exc}") from exc
        if resolution != self.resolution:
            raise ValueError(
                f"balance sidecar resolution {resolution} != mcool resolution {self.resolution}"
            )
        if semantics != "multiplicative":
            raise ValueError(
                f"balance sidecar semantics must be multiplicative, got {semantics!r}"
            )
        if balance_name != "ICE_weight" or self.balance_name != balance_name:
            raise ValueError(
                "balance sidecar name must be ICE_weight and match the requested "
                f"balance name, got sidecar={balance_name!r}, requested={self.balance_name!r}"
            )
        if weights.ndim != 1 or len(weights) != len(self._bins_df):
            raise ValueError(
                f"balance sidecar weights shape {weights.shape} != ({len(self._bins_df)},)"
            )
        if np.isinf(weights).any():
            raise ValueError("balance sidecar contains infinite weights")
        finite = np.isfinite(weights)
        if np.any(weights[finite] <= 0):
            raise ValueError("balance sidecar contains non-positive finite weights")
        return weights

    # ── 公开接口 ──────────────────────────────

    @property
    def chromnames(self) -> list:
        """返回本次任务需要处理的染色体列表"""
        return self._target_chroms

    def bins_for(self, chrom: str) -> pd.DataFrame:
        return self._bins_df[self._bins_df["chrom"] == chrom].reset_index(drop=True)

    def matrix(self, chrom: str, balanced: bool = True,
               balance_name: str | None = None) -> np.ndarray:
        """
        利用索引（bin1_offset）实现 O(1) 内存开销的切片读取。

        当 ``balanced=True``（默认）时，应用构造器声明的 ``balance_name``
        （默认 ``weight``）。Cooler/ICE 的 ``weight`` 是乘法型；4DN/Juicer 的
        ``KR``、``VC``、``VC_SQRT`` 是除法型，与
        ``cooler.matrix(balance=<name>)`` 的语义一致：
            multiplicative: balanced[i, j] = raw[i, j] * w[i] * w[j]
            divisive:       balanced[i, j] = raw[i, j] / w[i] / w[j]
        权重为 NaN 的 bin（cooler balance 标记的低复杂度 / unmappable
        bin）在 balanced 矩阵里也是 NaN —— TAD caller 必须能 propagate
        NaN，不能用 0 替代（否则会污染局部统计量）。

        若声明的列不存在，普通交互模式会警告并回退 raw counts；
        ``strict_balance=True`` 的科研重建模式直接报错，禁止静默混用 normalization。
        """
        print(f"  提取矩阵: {chrom} ...")
        start_time = time.time()

        # 1. 获取该染色体对应的全局 bin ID 范围
        chrom_bins = self._bins_df[self._bins_df["chrom"] == chrom]
        if chrom_bins.empty:
            raise ValueError(f"染色体 {chrom} 不存在")

        id_min = int(chrom_bins.index[0])
        id_max = int(chrom_bins.index[-1])
        n = id_max - id_min + 1

        weight = None
        requested_balance = self.balance_name if balance_name is None else str(balance_name)
        with h5py.File(self.cool_file, "r") as f:
            root = f[self._root_path]
            pix_grp = root["pixels"]
            # 关键：获取 bin1 的偏移索引，实现精准定位
            bin1_offsets = root["indexes"]["bin1_offset"][:]

            # 2. 确定该染色体在 pixels 数组中的物理起始位置
            idx_start = bin1_offsets[id_min]
            idx_end = bin1_offsets[id_max + 1]

            # 3. 仅读取该染色体相关的像素块（极大节省内存！）
            b1 = pix_grp["bin1_id"][idx_start:idx_end]
            b2 = pix_grp["bin2_id"][idx_start:idx_end]
            cnt = pix_grp["count"][idx_start:idx_end]

            # 4. 平衡权重（如果存在）
            bins_grp = root["bins"]
            if (balanced and requested_balance and self._external_balance_weights is not None
                    and requested_balance == self.balance_name):
                weight = self._external_balance_weights[id_min:id_max + 1]
            elif balanced and requested_balance and requested_balance in bins_grp:
                weight = bins_grp[requested_balance][id_min:id_max + 1].astype(float)
            elif balanced and requested_balance and self.strict_balance:
                available = sorted(k for k in bins_grp.keys()
                                   if k not in {"chrom", "start", "end"})
                raise ValueError(
                    f"requested balance column bins/{requested_balance!s} is missing at "
                    f"{self._root_path}; available={available}"
                )

        # 5. 构建密集矩阵
        # 注意：由于 .cool 是对称矩阵存储的上三角，我们只需要过滤 b2 越界的情况
        mask = (b2 >= id_min) & (b2 <= id_max)
        row = b1[mask] - id_min
        col = b2[mask] - id_min
        data = cnt[mask].astype(float)

        mat = np.zeros((n, n), dtype=float)
        mat[row, col] = data
        mat[col, row] = data  # 填充对称部分

        if weight is not None:
            # Match cooler's 4DN convention exactly: KR/VC/VC_SQRT are divisive;
            # ordinary cooler `weight` columns are multiplicative. NaN weights
            # propagate in either case so callers can mask low-complexity bins.
            wcol = weight[:, None]
            wrow = weight[None, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                if requested_balance in self._DIVISIVE_BALANCE_NAMES:
                    mat = mat / wcol / wrow
                else:
                    mat = mat * wcol * wrow
            n_bad = int((~np.isfinite(weight)).sum())
            source = "sidecar" if self._external_balance_weights is not None else "mcool"
            print(f"  应用 {requested_balance} 平衡权重 ({source}) | bad bins: {n_bad}/{n}")
        elif balanced:
            print(
                f"  [警告] mcool 该分辨率没有 {requested_balance} 列，使用原始 counts；"
                "科研重建请启用 strict_balance 并声明统一 normalization"
            )

        print(f"  完成 | 大小: {n}x{n} | 耗时: {time.time() - start_time:.2f}s")
        return mat


# ═══════════════════════════════════════════════
# 公共工具
# ═══════════════════════════════════════════════

def fmt_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def boundary_bp_to_bin_interval(chrom: str, boundary_bp: int, bins: pd.DataFrame,
                                resolution: int) -> tuple[str, int, int] | None:
    """Convert a boundary coordinate to the containing 0-based half-open bin."""
    if resolution <= 0:
        raise ValueError("resolution must be positive")
    idx = int(round(boundary_bp / resolution))
    if idx < 0 or idx >= len(bins):
        return None
    row = bins.iloc[idx]
    return fmt_chrom(str(row["chrom"])), int(row["start"]), int(row["end"])


# ═══════════════════════════════════════════════
# 方法一：DI (Directionality Index)
# ═══════════════════════════════════════════════

def calculate_DI(hic_matrix: np.ndarray, window_size: int = 200) -> np.ndarray:
    print("  计算 Directionality Index (DI)...")
    n  = hic_matrix.shape[0]
    DI = np.zeros(n)
    for i in tqdm(range(window_size, n - window_size), desc="  DI"):
        A = np.sum(hic_matrix[i, i - window_size:i])
        B = np.sum(hic_matrix[i, i + 1:i + window_size + 1])
        E = (A + B) / 2
        if E == 0:
            DI[i] = 0
        else:
            sign  = (B - A) / abs(B - A) if B != A else 0
            DI[i] = sign * (((A - E) ** 2 / E) + ((B - E) ** 2 / E))
    return DI


def detect_boundaries_DI(DI_values: np.ndarray, resolution: int,
                         min_tad_size: int = 5) -> list:
    """Canonical Dixon et al. 2012 DI boundary caller.

    Per the original Nature paper: feed DI values into a 3-state Gaussian HMM
    (biased upstream / neutral / biased downstream) and declare a boundary at
    every state transition.  The only filter applied is a minimum TAD length
    between consecutive state changes, which is the standard post-hoc
    requirement used by every downstream analysis.
    """
    print("  HMM boundary detection (Dixon 2012)...")
    vals = DI_values.reshape(-1, 1).copy()
    vals = np.nan_to_num(vals)
    vals[np.isinf(vals)] = 0
    model  = hmm.GaussianHMM(
        n_components=3,
        covariance_type="full",
        n_iter=100,
        random_state=0,
    )
    model.fit(vals)
    states = model.predict(vals)

    raw_bds = [i for i in range(1, len(states)) if states[i] != states[i - 1]]

    filtered = []
    for i, b in enumerate(raw_bds):
        next_b = raw_bds[i + 1] if i + 1 < len(raw_bds) else len(states)
        if next_b - b >= min_tad_size:
            filtered.append(b * resolution)

    print(f"  Raw state changes: {len(raw_bds)}, after min-size filter ({min_tad_size} bins): {len(filtered)}")
    return filtered


def run_DI(reader: McoolReader, output_file: str):
    print("\n[1/4] 运行 DI 方法...")
    res    = reader.resolution
    chroms = reader.chromnames
    open(output_file, "w").close()
    ok = 0
    for i, chrom in enumerate(chroms):
        print(f"\n  进度: {i+1}/{len(chroms)} — {chrom}")
        try:
            mat = reader.matrix(chrom)
            bins = reader.bins_for(chrom)
            DI  = calculate_DI(mat)
            bds = detect_boundaries_DI(DI, res)
            with open(output_file, "a") as f:
                for b in bds:
                    interval = boundary_bp_to_bin_interval(chrom, b, bins, res)
                    if interval is not None:
                        f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")
            ok += 1
        except Exception as exc:
            print(f"  !! 染色体 {chrom} 出错: {exc}")
    print(f"  DI 完成，成功处理 {ok}/{len(chroms)} 条染色体 → {output_file}")


# ═══════════════════════════════════════════════
# 方法二：TopDom (Shin et al. 2016, Nucleic Acids Research)
# ═══════════════════════════════════════════════

def calculate_TopDom(hic_matrix: np.ndarray, window_size: int = 20) -> np.ndarray:
    """TopDom binSignal (Shin et al. 2016, NAR 44:e70).

    For each bin i, compute the mean contact frequency inside the W × W
    off-diagonal square spanning [i-W, i-1] × [i+1, i+W].  This is the "diamond
    above the bin" region whose signal drops at a true TAD boundary.  W=20 bins
    at 10 kb resolution recovers the original paper's 200 kb window scale.
    """
    print(f"  TopDom binSignal 计算（窗口={window_size} bins，Shin 2016）...")
    n = hic_matrix.shape[0]
    binSignal = np.full(n, np.nan)
    for i in tqdm(range(window_size, n - window_size), desc="  TopDom"):
        diamond = hic_matrix[i - window_size:i, i + 1:i + window_size + 1]
        if diamond.size == 0:
            continue
        binSignal[i] = np.nanmean(diamond)
    return binSignal


def detect_boundaries_TopDom(binSignal: np.ndarray, hic_matrix: np.ndarray,
                             resolution: int,
                             window_size: int = 20,
                             p_thresh: float = 0.05,
                             min_tad_size: int = 10) -> list:
    """Canonical Shin 2016 TopDom boundary caller.

    1. Detect local minima of binSignal (candidate boundaries).
    2. For each candidate, run a one-sided Wilcoxon rank-sum test comparing
       within-TAD contacts (upper-triangular of the two flanking W×W on-diagonal
       blocks) against between-TAD contacts (the cross rectangle [i-W:i]×[i+1:i+W]).
       Accept the boundary iff the within-TAD distribution is significantly
       greater than the between-TAD distribution (p ≤ p_thresh).
    3. Enforce a minimum TAD size between consecutive accepted boundaries.
    """
    from scipy.stats import ranksums
    n = len(binSignal)

    # --- Step 1: local minima ---
    candidates = []
    for i in range(1, n - 1):
        bs = binSignal[i]
        if np.isnan(bs):
            continue
        left, right = binSignal[i - 1], binSignal[i + 1]
        if np.isnan(left) or np.isnan(right):
            continue
        if bs <= left and bs <= right and (bs < left or bs < right):
            candidates.append(i)

    # --- Step 2: Wilcoxon rank-sum test (Shin 2016 Eq. 4) ---
    accepted = []
    for i in candidates:
        lo = max(0, i - window_size)
        hi = min(n, i + window_size + 1)
        up = hic_matrix[lo:i, lo:i]
        dn = hic_matrix[i + 1:hi, i + 1:hi]
        between = hic_matrix[lo:i, i + 1:hi]
        if up.shape[0] < 2 or dn.shape[0] < 2 or between.size < 3:
            continue
        within = np.concatenate([
            up[np.triu_indices_from(up, k=1)],
            dn[np.triu_indices_from(dn, k=1)],
        ])
        within = within[~np.isnan(within)]
        bet    = between[~np.isnan(between)].ravel()
        if within.size < 3 or bet.size < 3:
            continue
        try:
            _, p = ranksums(within, bet, alternative="greater")
        except TypeError:
            _, p = ranksums(within, bet)
        if p <= p_thresh:
            accepted.append(i)

    # --- Step 3: min TAD size filter ---
    filtered = []
    for idx, b in enumerate(accepted):
        next_b = accepted[idx + 1] if idx + 1 < len(accepted) else n
        if next_b - b >= min_tad_size:
            filtered.append(b * resolution)

    print(f"  Local minima: {len(candidates)}, "
          f"after Wilcoxon (p≤{p_thresh}): {len(accepted)}, "
          f"after min-size filter ({min_tad_size} bins): {len(filtered)}")
    return filtered


def run_TopDom(reader: McoolReader, output_file: str):
    print("\n[2/4] 运行 TopDom 方法 (Shin 2016)...")
    res    = reader.resolution
    chroms = reader.chromnames
    open(output_file, "w").close()
    ok = 0
    for i, chrom in enumerate(chroms):
        print(f"\n  进度: {i+1}/{len(chroms)} — {chrom}")
        try:
            mat = reader.matrix(chrom)
            bins = reader.bins_for(chrom)
            bs  = calculate_TopDom(mat)
            bds = detect_boundaries_TopDom(bs, mat, res)
            with open(output_file, "a") as f:
                for b in bds:
                    interval = boundary_bp_to_bin_interval(chrom, b, bins, res)
                    if interval is not None:
                        f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")
            ok += 1
        except Exception as exc:
            print(f"  !! 染色体 {chrom} 出错: {exc}")
    print(f"  TopDom 完成，成功处理 {ok}/{len(chroms)} 条染色体 → {output_file}")


# ═══════════════════════════════════════════════
# 方法三：IS (Insulation Score, Crane et al. 2015 Nature 523:240)
# ═══════════════════════════════════════════════

def calculate_IS(hic_matrix: np.ndarray, window_size: int) -> np.ndarray:
    """Canonical Crane 2015 insulation score.

    For each bin i, score = mean contact inside the W × W off-diagonal square
    [i-W, i-1] × [i+1, i+W] (the "diamond above the bin" region) — i.e. the
    average contact between the upstream window and the downstream window.
    Wider, symmetric (and with no self-diagonal) than a raw 2W×2W sum.
    """
    print(f"  计算 Insulation Score（窗口={window_size} bins，Crane 2015）...")
    n  = hic_matrix.shape[0]
    IS = np.full(n, np.nan)
    for i in tqdm(range(window_size, n - window_size), desc="  IS"):
        block = hic_matrix[i - window_size:i, i + 1:i + window_size + 1]
        if block.size == 0:
            continue
        IS[i] = np.nanmean(block)
    return IS


def detect_boundaries_IS(is_values: np.ndarray, threshold: float,
                         min_tad_size: int = 5) -> list:
    """Canonical Crane 2015 boundary detection.

    1. Normalize: score = log2(IS / mean(IS)) — mean over non-NaN bins.
    2. Detect local minima of the normalized score below ``-threshold``
       (threshold is given as a positive "delta" number in log2 units,
       with the published default 0.1).
    3. Enforce min TAD size between consecutive boundaries.
    """
    print(f"  IS boundary detection (log2 ratio, delta≤-{threshold})...")
    mean_is = np.nanmean(is_values)
    if mean_is <= 0 or np.isnan(mean_is):
        return []
    with np.errstate(divide="ignore", invalid="ignore"):
        log2_ratio = np.log2(np.where(is_values > 0, is_values, np.nan) / mean_is)

    raw_bds = []
    for i in range(1, len(log2_ratio) - 1):
        v = log2_ratio[i]
        if np.isnan(v):
            continue
        left, right = log2_ratio[i - 1], log2_ratio[i + 1]
        if np.isnan(left) or np.isnan(right):
            continue
        if v < left and v < right and v < -threshold:
            raw_bds.append(i)

    filtered = []
    for i, b in enumerate(raw_bds):
        next_b = raw_bds[i + 1] if i + 1 < len(raw_bds) else len(log2_ratio)
        if next_b - b >= min_tad_size:
            filtered.append(b)

    print(f"  Raw: {len(raw_bds)}, after min-size filter ({min_tad_size} bins): {len(filtered)}")
    return filtered


def run_IS(reader: McoolReader, threshold: float, output_file: str):
    print("\n[3/4] 运行 IS 方法...")
    res         = reader.resolution
    window_bins = 500_000 // res
    print(f"  窗口大小: {window_bins} bins (500kb)")
    chroms = reader.chromnames
    open(output_file, "w").close()
    ok = 0
    for i, chrom in enumerate(chroms):
        print(f"\n  进度: {i+1}/{len(chroms)} — {chrom}")
        try:
            mat     = reader.matrix(chrom)
            bins    = reader.bins_for(chrom)
            is_vals = calculate_IS(mat, window_bins)
            bds     = detect_boundaries_IS(is_vals, threshold)
            with open(output_file, "a") as f:
                for idx in bds:
                    if idx < len(bins):
                        row = bins.iloc[idx]
                        f.write(f"{fmt_chrom(row['chrom'])}\t{row['start']}\t{row['end']}\n")
            ok += 1
        except Exception as exc:
            print(f"  !! 染色体 {chrom} 出错: {exc}")
    print(f"  IS 完成，成功处理 {ok}/{len(chroms)} 条染色体 → {output_file}")


# ═══════════════════════════════════════════════
# 快速模式：每条染色体只加载一次矩阵，三方法共用
# ═══════════════════════════════════════════════

def run_all_methods(reader: McoolReader, is_threshold: float,
                    di_file: str, topdom_file: str, is_file: str,
                    min_tad_bins: int = 5):
    """
    每条染色体只读取一次矩阵，依次运行 DI / TopDom / IS，
    结果分别追加到各自的输出文件。三种都是严格的文献原版方法：
      DI      — Dixon et al. 2012 Nature (window=2Mb + HMM)
      TopDom  — Shin et al. 2016 NAR    (binSignal + Wilcoxon)
      IS      — Crane et al. 2015 Nature (log2 ratio + local minima)
    """
    print("\n[快速模式] 每条染色体单次加载，依次运行 DI / TopDom / IS ...")
    res         = reader.resolution
    window_bins = 500_000 // res
    chroms      = reader.chromnames

    for path in (di_file, topdom_file, is_file):
        open(path, "w").close()

    ok = 0
    t_total = time.time()

    for idx, chrom in enumerate(chroms):
        print(f"\n{'─'*50}")
        print(f"  [{idx+1}/{len(chroms)}] 染色体 {chrom}")
        print(f"{'─'*50}")

        try:
            mat  = reader.matrix(chrom)
            bins = reader.bins_for(chrom)

            # ── DI (Dixon 2012) ──
            print("  → DI")
            DI  = calculate_DI(mat)
            bds = detect_boundaries_DI(DI, res, min_tad_size=min_tad_bins)
            with open(di_file, "a") as f:
                for b in bds:
                    interval = boundary_bp_to_bin_interval(chrom, b, bins, res)
                    if interval is not None:
                        f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")

            # ── TopDom (Shin 2016) ──
            print("  → TopDom")
            td_signal = calculate_TopDom(mat)
            bds = detect_boundaries_TopDom(td_signal, mat, res, min_tad_size=min_tad_bins)
            with open(topdom_file, "a") as f:
                for b in bds:
                    interval = boundary_bp_to_bin_interval(chrom, b, bins, res)
                    if interval is not None:
                        f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")

            # ── IS (Crane 2015) ──
            print("  → IS")
            is_vals = calculate_IS(mat, window_bins)
            bds     = detect_boundaries_IS(is_vals, is_threshold, min_tad_size=min_tad_bins)
            with open(is_file, "a") as f:
                for i in bds:
                    if i < len(bins):
                        row = bins.iloc[i]
                        f.write(f"{fmt_chrom(row['chrom'])}\t{row['start']}\t{row['end']}\n")

            ok += 1
        except Exception as exc:
            print(f"  !! 染色体 {chrom} 出错: {exc}")

    print(f"\n[快速模式] 完成，成功处理 {ok}/{len(chroms)} 条染色体，"
          f"总用时: {time.time()-t_total:.1f}s")
    print(f"  DI     → {di_file}")
    print(f"  TopDom → {topdom_file}")
    print(f"  IS     → {is_file}")


# ═══════════════════════════════════════════════
# 步骤四：并集 - 交集（纯 pandas，替代 pybedtools/bedtools）
# ═══════════════════════════════════════════════

def load_bed(path: str) -> pd.DataFrame:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end"], usecols=[0, 1, 2])
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    return df


def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """等价于 bedtools merge。"""
    if df.empty:
        return df.copy()
    out = []
    for chrom, grp in df.groupby("chrom", sort=False):
        grp = grp.sort_values("start")
        cs, ce = int(grp.iloc[0]["start"]), int(grp.iloc[0]["end"])
        for _, row in grp.iloc[1:].iterrows():
            if int(row["start"]) <= ce:
                ce = max(ce, int(row["end"]))
            else:
                out.append((chrom, cs, ce))
                cs, ce = int(row["start"]), int(row["end"])
        out.append((chrom, cs, ce))
    return pd.DataFrame(out, columns=["chrom", "start", "end"])


def _has_overlap(as_, ae, bs, be, f: float) -> bool:
    inter = max(0, min(ae, be) - max(as_, bs))
    if inter == 0:
        return False
    return (inter / (ae - as_) >= f) and (inter / (be - bs) >= f)


def intersect_rows(df_a: pd.DataFrame, df_b: pd.DataFrame,
                   f: float) -> pd.DataFrame:
    """返回 df_a 中与 df_b 至少有 f 比例重叠的行（等价 bedtools intersect -u）。"""
    if df_a.empty or df_b.empty:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    b_by_chrom = {c: g.values for c, g in df_b.groupby("chrom")}
    hits = []
    for row in df_a.itertuples(index=False):
        grp = b_by_chrom.get(row.chrom)
        if grp is None:
            continue
        for rb in grp:                      # rb: (chrom, start, end)
            if _has_overlap(row.start, row.end, int(rb[1]), int(rb[2]), f):
                hits.append((row.chrom, row.start, row.end))
                break
    return pd.DataFrame(hits, columns=["chrom", "start", "end"])


def subtract_rows(df_a: pd.DataFrame, df_b: pd.DataFrame,
                  f: float) -> pd.DataFrame:
    """返回 df_a 中 *不* 与 df_b 重叠的行（等价 bedtools intersect -v）。"""
    if df_b.empty:
        return df_a.copy()
    overlap_set = set(
        map(tuple, intersect_rows(df_a, df_b, f)[["chrom", "start", "end"]].values)
    )
    mask = [tuple(r) not in overlap_set
            for r in df_a[["chrom", "start", "end"]].values]
    return df_a[mask].reset_index(drop=True)


def _pad_boundaries(df: pd.DataFrame, pad: int) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    out = df.copy()
    out["start"] = (out["start"] - pad).clip(lower=0)
    out["end"]   = out["end"] + pad
    return merge_intervals(out)


def _nearby_match(row, other_df, tolerance):
    for _, orow in other_df.iterrows():
        if row["chrom"] != orow["chrom"]:
            continue
        if abs(row["start"] - orow["start"]) <= tolerance:
            return True
    return False


def run_intersection(di_file: str, topdom_file: str, is_file: str,
                     output_file: str, overlap: float,
                     tolerance: int = 50000):
    """3-method legacy consensus/disputed caller (DI + TopDom + IS).

    Per-original-bin vote with a `tolerance * 2` bp matching window:
      • votes == 1  → Disputed (only this method calls a boundary here)
      • votes >= 2  → Consensus (at least one other method agrees)

    Outputs both files as **4-column BED**: ``chrom\\tstart\\tend\\tmethod``.
    For Consensus rows the 4th column is the comma-joined set of methods
    that hit the same site; for Disputed it is the single supporting
    method name.  The frontend BedTrack reads this 4th column to colour
    each Disputed box by its contributing method so it visually matches
    the upstream method row.
    """
    print(f"\n[4/4] Computing consensus & disputed boundaries (tolerance={tolerance//1000}kb)...")

    di_df     = load_bed(di_file)
    topdom_df = load_bed(topdom_file)
    is_df     = load_bed(is_file)

    all_boundaries = pd.concat([
        di_df.assign(method="DI"),
        topdom_df.assign(method="TopDom"),
        is_df.assign(method="IS"),
    ], ignore_index=True)

    method_dfs = {"DI": di_df, "TopDom": topdom_df, "IS": is_df}
    consensus = []
    disputed = []
    for _, row in all_boundaries.iterrows():
        own_method = row["method"]
        supporters = [own_method]
        for m, odf in method_dfs.items():
            if m == own_method or odf.empty:
                continue
            chrom_match = odf[odf["chrom"] == row["chrom"]]
            if len(chrom_match) == 0:
                continue
            start_close = (chrom_match["start"] - row["start"]).abs()
            end_close = (chrom_match["end"] - row["end"]).abs()
            if (start_close + end_close).min() <= tolerance * 2:
                supporters.append(m)
        entry = {
            "chrom":  row["chrom"],
            "start":  int(row["start"]),
            "end":    int(row["end"]),
            "method": own_method,                # method that contributed THIS bin
            "support": ",".join(sorted(set(supporters))),
        }
        if len(supporters) >= 2:
            consensus.append(entry)
        else:
            disputed.append(entry)

    cons_cols = ["chrom", "start", "end", "support"]
    disp_cols = ["chrom", "start", "end", "method"]
    if consensus:
        consensus_df = (
            pd.DataFrame(consensus)
              .drop_duplicates(subset=["chrom", "start", "end"])
              .sort_values(["chrom", "start"])
              [cons_cols]
              .rename(columns={"support": "method"})  # for 4-col BED uniformity
        )
    else:
        consensus_df = pd.DataFrame(columns=disp_cols)
    if disputed:
        disputed_df = (
            pd.DataFrame(disputed)
              .drop_duplicates(subset=["chrom", "start", "end"])
              .sort_values(["chrom", "start"])
              [disp_cols]
        )
    else:
        disputed_df = pd.DataFrame(columns=disp_cols)

    consensus_df.to_csv(output_file, sep="\t", index=False, header=False)
    disputed_out = output_file.replace("TAD_final", "Disputed")
    disputed_df.to_csv(disputed_out, sep="\t", index=False, header=False)

    print(f"\n===== Boundary Statistics =====")
    print(f"  DI     boundaries: {len(di_df)}")
    print(f"  TopDom boundaries: {len(topdom_df)}")
    print(f"  IS     boundaries: {len(is_df)}")
    print(f"  Consensus (>=2 agree): {len(consensus_df)}  -> {output_file}")
    print(f"  Disputed  (1 only)   : {len(disputed_df)}  -> {disputed_out}")
    print(f"===============================")


# ═══════════════════════════════════════════════
# 主入口
# ═══════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description="TAD边界检测（DI + TopDom + IS + cooltools），纯h5py/pandas实现，Windows可用",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例: python tad_detection.py -i sample.mcool -r 10000 -t 0.1 -o 0.7 -b myhash"
    )
    p.add_argument("-i", "--input",        required=True,  help="输入 .mcool 文件路径")
    p.add_argument("-r", "--resolution",   required=True,  type=int,   help="分析分辨率 (bp)")
    p.add_argument("-t", "--is-threshold", default=0.5,    type=float, help="IS方法阈值 (默认: 0.5)")
    p.add_argument("-o", "--overlap",      default=0.5,    type=float, help="交集重叠度 0-1 (默认: 0.5)")
    p.add_argument("-b", "--hash",         default="",                 help="输出文件名标识符")
    p.add_argument("--fast",  action="store_true",
                   help="快速模式：每条染色体只加载一次矩阵，三方法共用（推荐）")
    return p.parse_args()


def detect_TAD_boundaries(
        mcool_file: str,
        resolution: int,
        selected_chroms: list = None,
        is_threshold: float = 0.5,
        overlap: float = 0.5,
        hash_tag: str = "",
        fast_mode: bool = True,
        output_dir: str = "./tad_results",
        min_tad_bins: int = 10,
        enable_cooltools: bool = False,
):
    """3-method legacy TAD boundary pipeline (DI + TopDom + IS).

    Restored as the default visualisation pipeline on 2026-05-14: matches the
    demo BED files shipped in ``backend/tad_results`` (chr7canon / chr17canon)
    and the historical 4-track UI layout (DI / TopDom / IS / Disputed).

    ``Disputed`` here is the *legacy* definition produced by
    :func:`run_intersection`: a per-original-bin vote with a 100 kb tolerance
    window — a method bin is Disputed when no other method has any boundary
    within ~100 kb of it.  Output rows are 10 kb single bins, aligned with
    every method track row pixel-for-pixel.

    Set ``enable_cooltools=True`` to additionally run the cooltools
    insulation caller as a 4th method and emit ``Disputed_4methods`` /
    ``Consensus_4methods`` files; the returned ``Disputed`` key, however,
    keeps pointing to the legacy 3-method Disputed for UI compatibility.
    """
    os.makedirs(output_dir, exist_ok=True)

    # 如果没有传入，默认跑全基因组（保持向下兼容）
    if not selected_chroms:
        selected_chroms = []

        # 文件命名加入标识，方便区分不同染色体的运行结果
    chrom_suffix = "_all" if not selected_chroms else f"_{len(selected_chroms)}chroms"

    di_out     = f"{output_dir}/DI_res{resolution}_{hash_tag}{chrom_suffix}.bed"
    topdom_out = f"{output_dir}/TopDom_res{resolution}_{hash_tag}{chrom_suffix}.bed"
    is_out     = f"{output_dir}/IS_res{resolution}_{hash_tag}{chrom_suffix}.bed"
    union_out  = f"{output_dir}/TAD_final_res{resolution}_{hash_tag}{chrom_suffix}.bed"
    # run_intersection writes the legacy 3-method Disputed alongside TAD_final
    disp3_out  = union_out.replace("TAD_final", "Disputed")

    if mcool_file.lower().endswith('.hic'):
        reader = HicReaderAdapter(mcool_file, resolution, chroms=selected_chroms)
    else:
        reader = McoolReader(mcool_file, resolution, chroms=selected_chroms)

    if fast_mode:
        run_all_methods(reader, is_threshold, di_out, topdom_out, is_out, min_tad_bins=min_tad_bins)
    else:
        run_DI(reader, di_out)
        run_TopDom(reader, topdom_out)
        run_IS(reader, is_threshold, is_out)

    # Legacy 3-method consensus/disputed via per-bin vote (tolerance=50 kb*2)
    run_intersection(di_out, topdom_out, is_out, union_out, overlap)

    cooltools_out = None
    disp4_out = None
    cons4_out = None
    if enable_cooltools:
        cooltools_out = f"{output_dir}/cooltools_res{resolution}_{hash_tag}{chrom_suffix}.bed"
        cooltools_ok = _run_cooltools_if_available(
            mcool_file, resolution, cooltools_out,
            selected_chroms=selected_chroms,
        )
        method_beds = {"DI": di_out, "TopDom": topdom_out, "IS": is_out}
        if cooltools_ok:
            method_beds["cooltools"] = cooltools_out
        disp4_out = f"{output_dir}/Disputed_{len(method_beds)}methods_{hash_tag}{chrom_suffix}.bed"
        cons4_out = f"{output_dir}/Consensus_{len(method_beds)}methods_{hash_tag}{chrom_suffix}.bed"
        _build_consensus_disputed(method_beds, disp4_out, cons4_out, tol=50000)

    print(f"\n===== TAD检测完成 ({'Partial' if selected_chroms else 'Full'}) =====")
    return {
        "DI": di_out,
        "TopDom": topdom_out,
        "IS": is_out,
        "cooltools": cooltools_out,
        "TAD_final": union_out,
        # Default Disputed = 3-method legacy definition (matches demo BEDs).
        "Disputed": disp3_out,
        "Disputed_4methods": disp4_out,
        "Consensus_4methods": cons4_out,
    }


def _run_cooltools_if_available(mcool_path: str, resolution: int, out_bed: str,
                                window_bp: int = 250000,
                                selected_chroms: list | None = None) -> bool:
    """Strict Open2C cooltools.insulation caller.

    Uses the library's own ``threshold='Li'`` boundary selection (Li &
    Lee 1993 minimum cross-entropy threshold, the default published in
    Venev et al. 2024 eLife).  We report every bin where
    ``is_boundary_{window_bp}`` is True and apply **no** additional custom
    strength filter — staying faithful to the canonical library output.

    If ``selected_chroms`` is provided, only boundaries on those chromosomes
    are written to the BED (matches the chromosome filter applied to the
    classical callers via ``McoolReader``).
    """
    try:
        import cooler
        import cooltools
    except ImportError:
        print("[cooltools] package not installed — skipping 4th method")
        return False
    try:
        clr = cooler.Cooler(f"{mcool_path}::/resolutions/{resolution}")
        if "weight" not in clr.bins().columns:
            try:
                import cooler.balance  # noqa
                print("[cooltools] mcool not balanced, running cooler balance...")
                import subprocess
                subprocess.run(
                    ["cooler", "balance", f"{mcool_path}::/resolutions/{resolution}"],
                    check=True, capture_output=True, timeout=900,
                )
                clr = cooler.Cooler(f"{mcool_path}::/resolutions/{resolution}")
            except Exception as exc:
                print(f"[cooltools] balance failed: {exc}")
                return False
        ins_tab = cooltools.insulation(
            clr, [window_bp],
            ignore_diags=2, min_dist_bad_bin=3,
            threshold="Li", append_raw_scores=False, verbose=False,
        )
        col_bnd = f"is_boundary_{window_bp}"
        if col_bnd not in ins_tab.columns:
            return False
        b = ins_tab[ins_tab[col_bnd]].copy()

        # Match the chromosome filter applied to DI/TopDom/IS via McoolReader.
        if selected_chroms:
            # Accept both "chr7" and "7" forms: normalize both sides to the
            # set that the library's `ins_tab["chrom"]` column uses.
            chrom_values = set(b["chrom"].astype(str).unique())
            wanted = set()
            for c in selected_chroms:
                cs = str(c)
                if cs in chrom_values:
                    wanted.add(cs)
                else:
                    alt = cs[3:] if cs.startswith("chr") else f"chr{cs}"
                    if alt in chrom_values:
                        wanted.add(alt)
            if wanted:
                b = b[b["chrom"].astype(str).isin(wanted)].copy()

        b[["chrom", "start", "end"]].to_csv(out_bed, sep="\t", header=False, index=False)
        chrom_desc = f"chroms={sorted(b['chrom'].unique().tolist())}" if selected_chroms else "all chroms"
        print(f"[cooltools] wrote {len(b)} boundaries → {out_bed} "
              f"(threshold=Li, window={window_bp} bp, {chrom_desc})")
        return True
    except Exception as exc:
        print(f"[cooltools] failed: {exc}")
        return False


def _build_consensus_disputed(method_beds: dict, out_disputed: str, out_consensus: str,
                              tol: int = 50000) -> None:
    """Aggregate boundary BEDs from N methods and split into Consensus vs Disputed.

    Semantics (user-requested):
      • Consensus = sites where ALL N methods agree (within ``tol`` bp).
      • Disputed  = sites where EXACTLY ONE method calls a boundary AND no
        other method has a boundary within ``tol`` bp of it.  These are the
        truly contested sites that benefit most from human annotation.
      • Sites with 2 to N-1 methods agreeing are "partial agreement" — they
        are real-ish but not annotation targets; we skip them.

    ``tol`` is the merge tolerance used when grouping nearby boundaries across
    methods into a single site.  50 kb = 5 bins at 10 kb resolution, a standard
    window used in TAD-boundary consensus papers (e.g. Zufferey et al 2018).
    """
    import pandas as pd
    all_rows = []
    for name, path in method_beds.items():
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            continue
        df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                         names=["chrom", "start", "end"])
        df["method"] = name
        df["chrom"] = df["chrom"].astype(str).str.replace(r"^chr", "", regex=True)
        all_rows.append(df)
    if not all_rows:
        pd.DataFrame(columns=["chrom", "start", "end"]).to_csv(out_disputed, sep="\t", header=False, index=False)
        pd.DataFrame(columns=["chrom", "start", "end"]).to_csv(out_consensus, sep="\t", header=False, index=False)
        return
    merged = pd.concat(all_rows, ignore_index=True).sort_values(["chrom", "start"]).reset_index(drop=True)

    # ─── Two-stage algorithm (locked, matches user expectation) ──────────
    # Stage 1 (vote judgement):
    #   Group neighbouring boundary intervals (across methods) that are
    #   within ``tol`` bp into a single "site".  vote(site) = number of
    #   distinct methods that contributed any bin to that site.
    # Stage 2 (per-bin output):
    #   Walk every original method bin and tag it with its site's vote.
    #   Disputed  = bins whose site has vote == 1 (only one method).
    #   Consensus = bins whose site has vote == N (all methods).
    # Each output row is a real 10 kb method bin (same width as the
    # method tracks above), de-duplicated across methods.
    # ─────────────────────────────────────────────────────────────────────
    site_assignments: list[tuple[str, int, int, int]] = []  # (chrom, start, end, vote)
    for chrom, g in merged.groupby("chrom", sort=False):
        cur_end = None
        cur_methods: set = set()
        cur_bins: list[tuple[str, int, int]] = []
        def _flush():
            vote = len(cur_methods)
            for cbin in cur_bins:
                site_assignments.append((cbin[0], cbin[1], cbin[2], vote))
        for _, row in g.iterrows():
            if cur_end is None:
                cur_end = int(row.end)
                cur_methods = {row.method}
                cur_bins = [(chrom, int(row.start), int(row.end))]
                continue
            if int(row.start) <= cur_end + tol:
                cur_end = max(cur_end, int(row.end))
                cur_methods.add(row.method)
                cur_bins.append((chrom, int(row.start), int(row.end)))
            else:
                _flush()
                cur_end = int(row.end)
                cur_methods = {row.method}
                cur_bins = [(chrom, int(row.start), int(row.end))]
        if cur_end is not None:
            _flush()

    assigned = pd.DataFrame(site_assignments, columns=["chrom", "start", "end", "vote"])
    n_total = len(method_beds)
    # De-duplicate per (chrom, start, end): if multiple methods contribute
    # the same exact bin, that bin appears once.  Vote is the max across
    # those rows (equivalent because all contributors share the same site).
    grouped = (
        assigned.groupby(["chrom", "start", "end"], sort=False)["vote"].max().reset_index()
    )
    cons = grouped[grouped["vote"] == n_total][["chrom", "start", "end"]].copy()
    disp = grouped[grouped["vote"] == 1       ][["chrom", "start", "end"]].copy()
    cons = cons.assign(chrom="chr" + cons["chrom"].astype(str)).sort_values(["chrom", "start"])
    disp = disp.assign(chrom="chr" + disp["chrom"].astype(str)).sort_values(["chrom", "start"])
    cons.to_csv(out_consensus, sep="\t", header=False, index=False)
    disp.to_csv(out_disputed,  sep="\t", header=False, index=False)

    # Breakdown by agreement level for logging (over unique bins)
    vote_dist = grouped["vote"].value_counts().sort_index()
    print(f"[consensus/disputed] tol={tol//1000}kb  total_unique_bins={len(grouped)}")
    for k, v in vote_dist.items():
        tag = "Disputed" if k == 1 else ("Consensus" if k == n_total else "partial")
        print(f"    bins at vote {k}/{n_total}: {v:4d}  ({tag})")
    print(f"  → Consensus (vote=={n_total}): {len(cons)} bins  → {out_consensus}")
    print(f"  → Disputed  (vote==1 ): {len(disp)} bins  → {out_disputed}")

# ═══════════════════════════════════════════════
# v2 4-method pipeline (insulation / topdom_like / contact_contrast / spectral_profile)
# ═══════════════════════════════════════════════

V2_METHOD_IDS = ("insulation", "topdom_like", "contact_contrast", "spectral_profile")


def detect_TAD_boundaries_v2(
        mcool_file: str,
        resolution: int,
        selected_chroms: list | None = None,
        is_threshold: float = 0.5,
        overlap: float = 0.5,
        hash_tag: str = "",
        fast_mode: bool = True,
        output_dir: str = "./tad_results",
        min_tad_bins: int = 10,
) -> dict:
    """builtin_callers_v2 4-method TAD boundary pipeline.

    Replaces the legacy DI / TopDom / IS / cooltools chain with the
    clean-room v2 callers (insulation, topdom_like, contact_contrast,
    spectral_profile).  DI was dropped 2026-05-13 after a sensitivity
    audit produced all-zero calls under cooler-balanced inputs.

    Signature mirrors :func:`detect_TAD_boundaries` so :mod:`main.py` can
    swap entry points without touching the request handler.  Parameters
    ``is_threshold``, ``overlap``, ``min_tad_bins`` are accepted but
    ignored (each v2 caller uses its own dataclass defaults).

    Returns a dict with the four per-method BED paths plus
    ``Disputed`` / ``Consensus`` / ``TAD_final`` (TAD_final == Consensus
    for downstream-prediction backward compatibility).
    """
    from builtin_callers_v2 import run_all_callers
    from builtin_callers_v2.insulation import InsulationParams
    from builtin_callers_v2.thresholds import QuantileThreshold

    os.makedirs(output_dir, exist_ok=True)
    if not selected_chroms:
        selected_chroms = []
    chrom_suffix = "_all" if not selected_chroms else f"_{len(selected_chroms)}chroms"

    if mcool_file.lower().endswith(".hic"):
        reader = HicReaderAdapter(mcool_file, resolution, chroms=selected_chroms)
    else:
        reader = McoolReader(mcool_file, resolution, chroms=selected_chroms)

    # Web pipeline override: stricter insulation threshold so the visualised
    # boundaries are the top-quartile strongest valleys (cleaner heatmap
    # overlay).  The cleanroom evaluation pipeline still uses the unmodified
    # builtin_callers_v2 default (QuantileThreshold(0.5)).
    per_method_params = {
        "insulation": InsulationParams(threshold_strategy=QuantileThreshold(0.75)),
    }

    method_beds = {
        m: f"{output_dir}/{m}_res{resolution}_{hash_tag}{chrom_suffix}.bed"
        for m in V2_METHOD_IDS
    }
    for path in method_beds.values():
        open(path, "w").close()

    chroms = reader.chromnames
    print(f"\n[v2 4-method] resolution={resolution}  chroms={len(chroms)}  "
          f"methods={list(V2_METHOD_IDS)}")
    t_total = time.time()
    ok = 0
    per_method_counts = {m: 0 for m in V2_METHOD_IDS}

    for idx, chrom in enumerate(chroms):
        print(f"\n  [{idx+1}/{len(chroms)}] chromosome {chrom}")
        try:
            mat = reader.matrix(chrom)
        except Exception as exc:
            print(f"  !! load matrix failed: {exc}")
            continue
        chrom_label = fmt_chrom(str(chrom))
        try:
            results = run_all_callers(
                mat, chrom_label, int(resolution),
                params_by_method=per_method_params,
            )
        except Exception as exc:
            print(f"  !! run_all_callers failed on {chrom}: {exc}")
            continue
        any_success = False
        for m in V2_METHOD_IDS:
            res = results.get(m)
            if res is None:
                continue
            bdf = res.boundaries
            if bdf is None or bdf.empty:
                continue
            with open(method_beds[m], "a") as fh:
                for row in bdf.itertuples(index=False):
                    fh.write(f"{row.chrom}\t{int(row.start)}\t{int(row.end)}\n")
            per_method_counts[m] += len(bdf)
            any_success = True
        if any_success:
            ok += 1

    print(
        f"\n[v2 4-method] {ok}/{len(chroms)} chromosomes produced "
        f"boundaries  total_time={time.time()-t_total:.1f}s"
    )
    for m in V2_METHOD_IDS:
        print(f"    {m}: {per_method_counts[m]} boundaries  → {method_beds[m]}")

    disp_out = f"{output_dir}/Disputed_4methods_v2_{hash_tag}{chrom_suffix}.bed"
    cons_out = f"{output_dir}/Consensus_4methods_v2_{hash_tag}{chrom_suffix}.bed"
    _build_consensus_disputed(method_beds, disp_out, cons_out, tol=50000)

    tad_final = f"{output_dir}/TAD_final_v2_{hash_tag}{chrom_suffix}.bed"
    if os.path.isfile(cons_out):
        import shutil
        shutil.copyfile(cons_out, tad_final)
    else:
        open(tad_final, "w").close()

    return {
        "insulation":       method_beds["insulation"],
        "topdom_like":      method_beds["topdom_like"],
        "contact_contrast": method_beds["contact_contrast"],
        "spectral_profile": method_beds["spectral_profile"],
        "TAD_final":        tad_final,
        "Disputed":         disp_out,
        "Consensus":        cons_out,
    }


# ═══════════════════════════════════════════════
# Active 5-method built-in panel + a per-import evidence tier computed on the
# uploaded data (not a static precomputed annotation).
#   BD1 = votes among {insulation, topdom_like, contact_contrast, spectral_profile,
#         network_modularity} on the user's own mcool
#   BD2/BD3 = CTCF / RAD21 ChIP fold IF the user supplied bigwig paths, else N/A
#   -> combine_tier -> D5..D1, written as a token-scoped annotation TSV the
#      visualization + sign-out read back. No "Disputed" track.
# ═══════════════════════════════════════════════
TADVCI_METHOD_IDS = ("insulation", "topdom_like", "contact_contrast",
                     "spectral_profile", "network_modularity")


def detect_TAD_boundaries_tadvci(
        mcool_file: str,
        resolution: int,
        selected_chroms: list | None = None,
        hash_tag: str = "",
        output_dir: str = "./tad_results",
        ctcf_bw: str | None = None,
        rad21_bw: str | None = None,
        user_bed_paths: list | None = None,
        loops_path: str | None = None,
        cell_line: str | None = None,
        assembly: str | None = None,
        upload_mode: str = "evaluate",
        tol: int = 50000,
        balance_name: str = "weight",
        strict_balance: bool = False,
        balance_weights_path: str | None = None,
        balance_id: str | None = None,
        enable_bcp: bool = False,
        **_ignored,
) -> dict:
    """Five-method built-in panel plus a deterministic TAD-VCI evidence tier on the
    uploaded mcool. Emits 5 per-method BEDs, a tier annotation TSV, and a tier BED.
    BD2/BD3 are measured only if ctcf_bw/rad21_bw are given (per-cell background).

    BYO: boundaries from user-supplied BED files (other algorithms) are folded into
    the candidate set so they get the same deterministic evidence tier; each candidate
    carries a `user_bed` flag (1 = a user boundary lies within tol)."""
    import numpy as np, pandas as pd
    from builtin_callers_v2 import run_all_callers
    from builtin_callers_v2.insulation import InsulationParams
    from builtin_callers_v2.thresholds import QuantileThreshold
    from tad_vci import annotate

    os.makedirs(output_dir, exist_ok=True)
    selected_chroms = selected_chroms or []
    suffix = "_all" if not selected_chroms else f"_{len(selected_chroms)}chroms"
    if mcool_file.lower().endswith(".hic"):
        reader = HicReaderAdapter(mcool_file, resolution, chroms=selected_chroms)
    else:
        reader = McoolReader(mcool_file, resolution, chroms=selected_chroms,
                             balance_name=balance_name, strict_balance=strict_balance,
                             balance_weights_path=balance_weights_path)
    # Production caller parameters (scientific review 2026-06-02). Set EXPLICITLY (not defaults)
    # so the values are reproducible + paper-documentable (see analyses/.../CALLER_PARAM_SPEC.md):
    #   insulation: prominence_window_bins=20 -> LOCAL (500 kb) valley-drop instead of full-arm.
    #     Fixes the edge bias that systematically dropped boundaries near chromosome starts
    #     (the former None = whole-arm prominence gave start-proximal minima artificially low
    #     prominence, so e.g. chr13's first boundary appeared only at 35 Mb instead of ~20 Mb).
    #     Verified: chr13 first boundary 35.75->20.0 Mb, chr7 6.85->0.93 Mb, total count
    #     essentially unchanged (chr13 154->152, chr7 244->243), so votes/tier shift minimally.
    #   spectral_profile: aggregate="mean" (not max) + min_distance_bins=4 -> removes the
    #     single-bin-spike false positives the max aggregation produced (~50% on synthetic).
    from builtin_callers_v2.spectral_profile import SpectralProfileParams
    per_method_params = {
        "insulation": InsulationParams(threshold_strategy=QuantileThreshold(0.75),
                                       prominence_window_bins=20),
        "spectral_profile": SpectralProfileParams(aggregate="mean", min_distance_bins=4),
    }

    method_beds = {m: f"{output_dir}/{m}_res{resolution}_{hash_tag}{suffix}.bed" for m in TADVCI_METHOD_IDS}
    for p in method_beds.values():
        open(p, "w").close()

    # per (chrom) -> per method -> sorted boundary start positions (for vote counting)
    by_chrom: dict = {}
    counts = {m: 0 for m in TADVCI_METHOD_IDS}
    chroms = reader.chromnames
    print(f"\n[TAD-VCI 5-method] resolution={resolution} chroms={len(chroms)} methods={list(TADVCI_METHOD_IDS)}")
    t0 = time.time(); ok = 0
    for idx, chrom in enumerate(chroms):
        try:
            mat = reader.matrix(chrom)
        except Exception as exc:
            print(f"  !! load matrix failed on {chrom}: {exc}"); continue
        clabel = fmt_chrom(str(chrom))
        try:
            results = run_all_callers(mat, clabel, int(resolution), params_by_method=per_method_params)
        except Exception as exc:
            print(f"  !! run_all_callers failed on {chrom}: {exc}"); continue
        per = {}
        for m in TADVCI_METHOD_IDS:
            res = results.get(m)
            bdf = None if res is None else res.boundaries
            if bdf is None or bdf.empty:
                per[m] = np.array([], int); continue
            starts = bdf["start"].astype(int).values
            with open(method_beds[m], "a") as fh:
                for row in bdf.itertuples(index=False):
                    fh.write(f"{row.chrom}\t{int(row.start)}\t{int(row.end)}\n")
            counts[m] += len(bdf); per[m] = np.sort(starts)
        by_chrom[clabel] = per
        if any(len(v) for v in per.values()):
            ok += 1
    print(f"[TAD-VCI 5-method] {ok}/{len(chroms)} chroms; counts={counts}; {time.time()-t0:.1f}s")

    # ---- fold in user-uploaded boundaries (BYO from other algorithms) ----
    user_by_cnum: dict = {}
    for ubed in (user_bed_paths or []):
        if not ubed or not os.path.isfile(ubed):
            continue
        try:
            udf = pd.read_csv(ubed, sep=r"\s+", header=None, comment="#",
                              usecols=[0, 1], names=["chrom", "start"], engine="python")
            udf = udf[pd.to_numeric(udf["start"], errors="coerce").notna()].copy()
            udf["chrom"] = udf["chrom"].astype(str).str.replace("^chr", "", regex=True)
            udf["start"] = udf["start"].astype(int)
            for cnum, g in udf.groupby("chrom"):
                user_by_cnum.setdefault(cnum, []).extend(g["start"].tolist())
        except Exception as exc:
            print(f"  !! user BED parse failed {ubed}: {exc}")
    user_by_cnum = {c: np.sort(np.array(v, int)) for c, v in user_by_cnum.items()}
    if user_by_cnum:
        print(f"[TAD-VCI BYO] folded user boundaries: {{{', '.join(f'{c}:{len(a)}' for c, a in user_by_cnum.items())}}}")

    # ---- consensus clustering + votes (BD1): cluster the 5-method ∪ user boundaries within
    #      `tol` so each method is counted at most once per consensus boundary (no vote
    #      double-counting across boundaries within tol of each other). ----
    # 'contribute' mode: the uploaded method joins the panel as an extra voter (BD1 graded
    # against n_methods = 5 + 1). 'evaluate' (default): upload is flagged + scored but never
    # votes, so the built-in consensus is unchanged. Mode only matters if uploads exist.
    contribute = (str(upload_mode).lower() == "contribute") and bool(user_by_cnum)
    n_methods = len(TADVCI_METHOD_IDS) + (1 if contribute else 0)
    rows = []
    by_cnum = {clabel.replace("chr", ""): per for clabel, per in by_chrom.items()}
    cnums = set(by_cnum) | set(user_by_cnum)
    for cnum in cnums:
        per = by_cnum.get(cnum, {m: np.array([], int) for m in TADVCI_METHOD_IDS})
        upos = user_by_cnum.get(cnum)
        # Per-chrom BD1 denominator: a method with zero boundaries on this chrom can never
        # vote here (e.g. Arrowhead when KR fails). Counting it in n_methods deflates grades.
        present = [m for m in TADVCI_METHOD_IDS if len(per.get(m, np.array([], int)))]
        n_eff = max(1, len(present) + (1 if contribute and upos is not None and len(upos) else 0))
        absent = [m for m in TADVCI_METHOD_IDS if m not in present]
        for cl in cluster_consensus({m: per.get(m, np.array([], int)) for m in TADVCI_METHOD_IDS},
                                    tol, user_pos=upos, count_user_as_vote=contribute):
            rows.append({
                "chrom": cnum,
                "pos": cl["pos"],
                "votes": cl["votes"],
                "user_bed": cl["user_bed"],
                "n_methods": int(n_methods),
                "n_methods_effective": int(n_eff),
                "callers_absent_on_chrom": ";".join(absent) if absent else "",
            })
    df = pd.DataFrame(rows)

    # ---- optional BD4 loop-anchor count from a HiCCUPS-style loops file ----
    if len(df) and loops_path and os.path.isfile(loops_path):
        anchors = _parse_loop_anchors(loops_path)
        def _lc(c, p):
            k = str(c).replace("chr", "")
            return int(np.sum(np.abs(anchors[k] - p) <= tol)) if k in anchors else 0
        df["loop_anchor_count"] = [_lc(c, p) for c, p in zip(df["chrom"], df["pos"])]

    # ---- optional BD2/BD3 from user-supplied ChIP (per-cell background + PER-CELL bands) ----
    bg_ctcf = bg_rad21 = None
    from tad_vci.graded_evidence import load_bands, load_null_quantiles
    ctcf_bands, rad21_bands, bands_calibrated, bands_cell = load_bands(cell_line)
    bands_source = "packaged" if bands_calibrated else "in_situ"
    # Where each observed fold sits in this cell's random-window null (packaged, or
    # sampled in situ below). Without it "weak" cannot be told from "just cleared it".
    null_quantiles: dict = dict(load_null_quantiles(cell_line) or {}) if bands_calibrated else {}
    if len(df):
        if ctcf_bw and os.path.isfile(ctcf_bw):
            df["ctcf"], bg_ctcf = _bw_window_signal(ctcf_bw, df, resolution)
        if rad21_bw and os.path.isfile(rad21_bw):
            df["rad21"], bg_rad21 = _bw_window_signal(rad21_bw, df, resolution)
        # Unify the BD2/BD3 background to the CANONICAL random-window mean (the value the calibration
        # bands are percentiles of). The served bigWig header covered-mean drifts ~4% from it and biased
        # RAD21 grades one level high for ~12% of boundaries (audit #3 -> bd3_unify_report.json). For
        # calibrated cells, divide by the stored per-cell random-window mean instead of the header mean.
        if bands_calibrated:
            try:
                import json as _json
                from tad_vci.graded_evidence import _BANDS_JSON
                _cal = _json.loads(_BANDS_JSON.read_text())
                if bands_cell in _cal:
                    if bg_ctcf is not None and "ctcf" in df:
                        bg_ctcf = float(_cal[bands_cell]["CTCF"]["mean"])
                    if bg_rad21 is not None and "rad21" in df:
                        bg_rad21 = float(_cal[bands_cell]["RAD21"]["mean"])
            except Exception as _bg_err:
                print(f"[TAD-VCI] BD2/BD3 canonical-background override skipped (non-fatal): {_bg_err}")
        else:
            # UNCALIBRATED: never grade against GM12878 cutpoints with the user's own
            # track as numerator (one-directional inflation). Estimate in situ, or drop.
            from tad_vci.insitu_bands import estimate_bands_and_null
            got_ctcf = got_rad21 = False
            if ctcf_bw and os.path.isfile(ctcf_bw) and "ctcf" in df:
                est, est_bg, est_q = estimate_bands_and_null(ctcf_bw)
                if est:
                    ctcf_bands, bg_ctcf, got_ctcf = est, est_bg, True
                    if est_q:
                        null_quantiles["CTCF"] = est_q
            if rad21_bw and os.path.isfile(rad21_bw) and "rad21" in df:
                est, est_bg, est_q = estimate_bands_and_null(rad21_bw)
                if est:
                    rad21_bands, bg_rad21, got_rad21 = est, est_bg, True
                    if est_q:
                        null_quantiles["RAD21"] = est_q
            if "ctcf" in df and not got_ctcf:
                df = df.drop(columns=["ctcf"])
            if "rad21" in df and not got_rad21:
                df = df.drop(columns=["rad21"])
            bands_source = "in_situ"
            bands_cell = f"in_situ:{cell_line or 'unspecified'}"
        # The retired BD5 motif experiment is outside the evidence-card contract.
        # BD1 continuous companion: insulation boundary-strength from the Hi-C map (informational,
        # refines the discrete vote count; not a separate criterion / not in the tier).
        try:
            from builtin_callers_v2.insulation import (
                compute_insulation_production, InsulationParams, _resolve_params)
            _rd = McoolReader(mcool_file, resolution, balance_name=balance_name,
                              strict_balance=strict_balance,
                              balance_weights_path=balance_weights_path)
            _p = _resolve_params(InsulationParams(), resolution)
            _ins = {}
            for _ch in df["chrom"].astype(str).unique():
                for _name in (_ch, "chr" + _ch, _ch.replace("chr", "")):   # mcools differ: '1' vs 'chr1'
                    try:
                        _mx = np.nan_to_num(np.asarray(_rd.matrix(_name), float))
                        _ins[_ch] = compute_insulation_production(_mx, _p)["boundary_strength"]
                        break
                    except Exception:
                        continue
            df["insulation_strength"] = [
                (round(float(_ins[str(c)][int(p) // resolution]), 4)
                 if (str(c) in _ins and 0 <= int(p) // resolution < len(_ins[str(c)])
                     and np.isfinite(_ins[str(c)][int(p) // resolution])) else np.nan)
                for c, p in zip(df["chrom"], df["pos"])]
        except Exception as _i_err:
            print(f"[TAD-VCI] BD1 insulation companion skipped (non-fatal): {_i_err}")
        full = annotate(df, bg_ctcf=bg_ctcf or 1.0, bg_rad21=bg_rad21 or 1.0,
                        ctcf_bands=ctcf_bands, rad21_bands=rad21_bands, n_methods=n_methods,
                        null_quantiles=null_quantiles)
        full["bands_cell"] = bands_cell
        full["bands_source"] = bands_source
        full["bands_calibrated"] = bool(bands_calibrated)
        full["data_cell_line"] = str(cell_line or "")
        full["hic_balance"] = ("hic_reader_default" if mcool_file.lower().endswith(".hic")
                               else str(balance_name))
        full["hic_balance_source"] = (
            "sidecar" if balance_weights_path else
            ("hic_reader_default" if mcool_file.lower().endswith(".hic") else "mcool_bins")
        )
        full["hic_normalization_id"] = str(balance_id or balance_name)
        full["chip_window_bp"] = 50000
        full["tier_rule_version"] = "max_support_v2"
        full["vote_tol_bp"] = int(tol)
        # BD1 vote model, stated per row because the two product paths differ and the
        # difference is consequential. Here positions are CLUSTERED within tol first
        # (cluster_consensus), so each method contributes at most one vote per cluster and
        # each cluster is one candidate. The BED path (annotate_from_beds) grades the raw
        # union of caller calls instead, where one call can support several neighbouring
        # candidates -- see supporting_offsets_bp there.
        full["bd1_vote_model"] = "cluster_consensus"
        if not bands_calibrated and (ctcf_bw or rad21_bw):
            print(
                f"[TAD-VCI] BD2/BD3 bands NOT calibrated for '{cell_line}' "
                f"-> estimated in situ from the supplied track (bands_source=in_situ)"
            )
        # Fail closed if a historical caller tries to reactivate BCP.  Silent
        # omission would let a run appear to contain an audited score when the
        # entire label/split/model chain has been retired.
        if enable_bcp:
            raise RuntimeError(
                "enable_bcp=True is retired after the July 2026 provenance audit; "
                "no active pipeline may emit BCP or conformal fields"
            )
    else:
        full = df

    ann_tsv = f"{output_dir}/tadvci_annotation_{hash_tag}{suffix}.tsv"
    full.to_csv(ann_tsv, sep="\t", index=False)
    tier_bed = f"{output_dir}/tadvci_tier_{hash_tag}{suffix}.bed"
    with open(tier_bed, "w") as fh:
        if len(full):
            for r in full.itertuples(index=False):
                fh.write(f"{r.chrom}\t{int(r.pos)}\t{int(r.pos)+resolution}\t{r.D_tier}\n")
    if len(full):
        dist = full["D_tier"].value_counts().reindex(["D5","D4","D3","D2","D1"]).fillna(0).astype(int).to_dict()
        print(f"[TAD-VCI tier] {len(full)} boundaries; dist={dist}; -> {ann_tsv}")

    out = {m: method_beds[m] for m in TADVCI_METHOD_IDS}
    out["TADVCI_annotation"] = ann_tsv
    out["TADVCI_tier"] = tier_bed
    return out


def cluster_consensus(per_method: dict, tol: int, user_pos=None, count_user_as_vote: bool = False) -> list:
    """Cluster boundary positions within `tol` into consensus boundaries, counting each
    method AT MOST ONCE per cluster.

    Fixes a vote double-count: when the union kept boundaries that sit within `tol` of each
    other as separate entries (finer 25 kb dedup vs a 50 kb vote tolerance), each entry
    credited the SAME method's nearby call, inflating votes and creating phantom 'strong'
    boundaries supported by only one method (e.g. chr7:50,075,000 borrowed 3 votes from a
    real boundary 50 kb away at 50,125,000). Here every method-boundary point is assigned to
    exactly ONE cluster, so no call is counted toward more than one consensus boundary.

    Greedy clustering anchored to each cluster's start bounds cluster width to `tol` (no
    runaway single-linkage chaining). Deterministic (sorted input, median-member position).
    Returns [{pos, votes, methods(set), user_bed}], one per cluster."""
    pts = []
    for m, arr in per_method.items():
        for p in np.asarray(arr).tolist():
            pts.append((int(p), m))
    if user_pos is not None:
        for p in np.asarray(user_pos).tolist():
            pts.append((int(p), "__user__"))
    if not pts:
        return []
    pts.sort(key=lambda x: x[0])
    out, cur, anchor = [], [pts[0]], pts[0][0]

    def flush(cl):
        methods = set(m for _, m in cl if m != "__user__")
        user = any(m == "__user__" for _, m in cl)
        # In 'contribute' mode the uploaded method counts as ONE additional voter (BD1 grading
        # then uses n_methods = 5 + 1 so the bands stay comparable). In 'evaluate' mode the
        # upload is only a flag and never alters the consensus of the built-in panel.
        votes = len(methods) + (1 if (user and count_user_as_vote) else 0)
        positions = sorted(p for p, _ in cl)
        rep = int(positions[len(positions) // 2])    # median member -> on-grid, deterministic
        out.append({"pos": rep, "votes": votes, "methods": methods, "user_bed": int(user)})

    for p, m in pts[1:]:
        if p - anchor <= tol:
            cur.append((p, m))
        else:
            flush(cur); cur = [(p, m)]; anchor = p
    flush(cur)
    return out


def _parse_loop_anchors(loops_path: str) -> dict:
    """Parse a HiCCUPS-style loops file -> {chrom: sorted np.array of anchor midpoints}.

    Handles BOTH formats:
      * header form (HiCCUPS .txt): named columns chr1,x1,x2,chr2,y1,y2 (+ extras);
      * headerless 6-column .bedpe: chr1 x1 x2 chr2 y1 y2.
    Both loop ends are anchors. Chrom keys are normalised to bare numbers ('7', 'X')."""
    import numpy as np
    by: dict[str, list] = {}
    try:
        with open(loops_path) as _fh:
            first = _fh.readline().rstrip("\n").split("\t")
    except Exception:
        return {}
    has_header = any(t.strip().lower() in ("chr1", "x1", "chrom1", "#chr1") for t in first)
    try:
        if has_header:
            ldf = pd.read_csv(loops_path, sep="\t")
            cols = {c.lower(): c for c in ldf.columns}
            need = ("chr1", "x1", "x2", "chr2", "y1", "y2")
            if not all(k in cols for k in need):
                return {}
            recs = ldf[[cols["chr1"], cols["x1"], cols["x2"],
                        cols["chr2"], cols["y1"], cols["y2"]]].itertuples(index=False, name=None)
        else:
            # headerless .bedpe: first 6 columns are chr1 x1 x2 chr2 y1 y2
            ldf = pd.read_csv(loops_path, sep="\t", header=None, comment="#")
            if ldf.shape[1] < 6:
                return {}
            recs = ldf.iloc[:, :6].itertuples(index=False, name=None)
        for c1, x1, x2, c2, y1, y2 in recs:
            for cc, a, b in ((c1, x1, x2), (c2, y1, y2)):
                key = str(cc).replace("chr", "")
                by.setdefault(key, []).append((int(a) + int(b)) // 2)
    except Exception:
        return {}
    return {k: np.sort(np.asarray(v)) for k, v in by.items()}


def _bw_window_signal(bw_path: str, df, resolution: int, flank: int = 25000):
    """Mean ChIP signal in a SYMMETRIC +/-flank window per boundary + a random-window
    background mean.

    CANONICAL WINDOW (train==serve, fixed 2026-06-04): exactly [pos-flank, pos+flank]
    -- a single symmetric +/-25 kb window with NO +resolution term. This is the same
    window definition used to build the calibration null (calibrate_bands.random_windows:
    [p-FLANK, p+FLANK]). The previous asymmetric [pos-flank, pos+resolution+flank]
    window was 75 kb wide (vs the 50 kb calibration window), so a fold graded
    against the random-window bands mixed a 75 kb numerator with a 50 kb null -- a
    calibration/serving mismatch on top of the denominator mismatch. `resolution` is kept
    in the signature for backward compatibility but no longer affects the window.
    """
    import numpy as np, pyBigWig
    bw = pyBigWig.open(bw_path)
    def win(chrom, pos):
        # Accept either naming convention. ENCODE bigWigs ship both "chr7" and bare "7"
        # keys depending on the processing pipeline; hard-coding the "chr" prefix made
        # every lookup on a bare-named track miss, and each miss returned NaN, so BD2/BD3
        # silently graded not_assessable for the whole run instead of failing loudly.
        # Measured on the shipped chr7 tracks (keys are {"7": 159345973}): 1,528 of 1,528
        # candidates returned NaN before this fix.
        raw = str(chrom)
        stripped = raw[3:] if raw.startswith("chr") else raw
        avail = bw.chroms()
        c = next((k for k in (raw, stripped, "chr" + stripped) if k in avail), None)
        if c is None:
            return np.nan
        L = avail[c]
        if not L:
            return np.nan
        a, b = max(0, pos - flank), min(L, pos + flank)
        try:
            v = bw.stats(c, a, b, type="mean")[0]
        except (RuntimeError, KeyError):
            return np.nan
        return float(v) if v is not None else 0.0
    sig = np.array([win(c, p) for c, p in zip(df["chrom"], df["pos"])], float)
    # background = genome-wide covered mean (header) as a stable per-cell reference.
    # NOTE: for CALIBRATED cells the served path overrides this header mean with the
    # per-cell random-window mean (mcool_bed detect_TAD_boundaries_tadvci), so the served
    # fold = raw +/-25 kb mean / random-window mean -- the exact canonical fold the bands
    # are percentiles of.
    h = bw.header(); bg = h["sumData"] / h["nBasesCovered"] if h["nBasesCovered"] else 1.0
    bw.close()
    return sig, float(bg)


def main():
    results = detect_TAD_boundaries(
        mcool_file="example_data/GSE63525_GM12878_insitu_DpnII_combined.mcool",
        resolution=10000,
        is_threshold=0.3,
        overlap=0.7,
        hash_tag="sample1",
        fast_mode=True,
        output_dir="./"
    )

    print(results)


if __name__ == "__main__":
    main()
