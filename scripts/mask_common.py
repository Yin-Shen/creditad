"""Shared masking predicates for referee-2 Major 5 (A1).

The gap / blacklist / off-end predicates are lifted VERBATIM in behaviour from
review/referee2/r2_04_unmappable_and_windows.py so the v1.1 masking targets exactly the
rows the referee counted (4,366 gap / 416 blacklist / 12 off-end). Nothing is
re-derived here: same three reference files, same point-in-interval test on `pos`,
same chrom-name normalisation (leading "chr" stripped).

The tier engine is the SHIPPED one (backend/tad_vci/graded_evidence.py); combine_tier
and grade_votes are imported, never reimplemented.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260825

ROOT = Path("/home/coder/TAD/creditad")
V10 = ROOT / "example_data" / "multicell"
WORK = ROOT / "rebuild_2026-08-25_fixes" / "masking"
OUT = WORK / "out"
V11 = WORK / "multicell_v1.1"

sys.path.insert(0, str(ROOT / "backend"))
from tad_vci.graded_evidence import grade_votes, combine_tier  # noqa: E402

PACKAGES = ["GM12878_25kb", "GM12878_10kb",
            "IMR90_25kb", "IMR90_10kb",
            "HepG2_25kb", "HepG2_10kb"]

CHROM_SIZES_F = Path("/home/coder/_shared_data/raw/reference/hg38/hg38.chrom.sizes")
GAP_F = Path("/home/coder/_shared_data/raw/reference/hg38/hg38_gap_regions.bed")
BLACKLIST_F = Path("/home/coder/_shared_data/raw/reference/blacklist/ENCFF356LFX.bed.gz")

TIERS = ["D1", "D2", "D3", "D4", "D5"]
GRADES_NA = ["none", "weak", "moderate", "strong", "not_assessable"]


def load_chrom_sizes() -> dict[str, int]:
    sizes = {}
    for line in CHROM_SIZES_F.open():
        c, L = line.split()[:2]
        sizes[c.replace("chr", "")] = int(L)
    return sizes


def _intervals(path: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"], compression="infer",
                     dtype={0: str})
    df["chrom"] = df["chrom"].str.replace("^chr", "", regex=True)
    return {c: g[["start", "end"]].to_numpy() for c, g in df.groupby("chrom")}


def load_gap() -> dict[str, np.ndarray]:
    return _intervals(GAP_F)


def load_blacklist() -> dict[str, np.ndarray]:
    return _intervals(BLACKLIST_F)


def in_any(iv, pos: int) -> bool:
    """Point-in-interval, half-open [start, end) -- identical to r2_04.in_any."""
    if iv is None or not len(iv):
        return False
    return bool(np.any((iv[:, 0] <= pos) & (pos < iv[:, 1])))


def flag_rows(df: pd.DataFrame, sizes, gap_by_c, bl_by_c) -> pd.DataFrame:
    """Return a frame of boolean flags aligned to df.index."""
    chrom = df["chrom"].astype(str).to_numpy()
    pos = df["pos"].to_numpy()
    return pd.DataFrame({
        "off_end": [int(p) >= sizes.get(c, 0) for c, p in zip(chrom, pos)],
        "in_gap": [in_any(gap_by_c.get(c), int(p)) for c, p in zip(chrom, pos)],
        "in_blacklist": [in_any(bl_by_c.get(c), int(p)) for c, p in zip(chrom, pos)],
    }, index=df.index)


def mask_reason(in_gap: bool, in_bl: bool) -> str:
    if in_gap and in_bl:
        return "assembly_gap+encode_blacklist"
    if in_gap:
        return "assembly_gap"
    if in_bl:
        return "encode_blacklist"
    return ""


def load_package(pkg: str) -> pd.DataFrame:
    """Load v1.0 annotation.tsv with EVERY column kept as raw text.

    dtype=str is deliberate: the masking pass must be able to re-emit unmodified rows
    byte-for-byte, so no float parsing / re-formatting may touch them.
    """
    p = V10 / pkg / "annotation.tsv"
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    return df


def load_manifest(pkg: str, base: Path = V10) -> dict:
    with (base / pkg / "manifest.json").open() as fh:
        return json.load(fh)


def tier_vec(bd1, bd2, bd3, bd4):
    """Memoised shipped combine_tier over 4 grade sequences -> (tiers, na_strings)."""
    cache: dict[tuple, tuple[str, str]] = {}
    tiers, nas = [], []
    for k in zip(bd1, bd2, bd3, bd4):
        v = cache.get(k)
        if v is None:
            t, meta = combine_tier({"BD1": k[0], "BD2": k[1], "BD3": k[2], "BD4": k[3]})
            v = (t, ",".join(meta["not_assessable"]))
            cache[k] = v
        tiers.append(v[0])
        nas.append(v[1])
    return tiers, nas
