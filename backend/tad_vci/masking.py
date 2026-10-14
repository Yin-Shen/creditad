"""Position masks for the ChIP axes (BD2, BD3).

Inside an assembly gap, a centromeric or heterochromatic band, or an ENCODE exclusion-list
region, a ChIP track carries no interpretable signal: a low fold there is an absence of
measurement, not a measured absence. A mask is a BED file of such regions. A boundary
whose position (BED column 2, the bin start) lies inside any mask interval, half-open
[start, end), has BD2 and BD3 recorded as `not_assessable`, and its tier and
`not_assessable_criteria` are re-derived with the released rule (`combine_tier`). BD1
(caller support) and BD4 (loop anchors) are not changed: neither is a ChIP measurement at
the boundary coordinate. The raw `ctcf` / `rad21` / fold / null-percentile columns keep
their values, so the withheld measurement stays visible; the grade, not the number, is
withheld.

This is the rule the v1.1 evidence tables (data/multicell_v1.1/) were built with, applied
after grading, so the cutpoints, the background and every unmasked row are exactly what
an unmasked run gives. Chromosome names are compared without a leading "chr".

`hg38` names the mask bundled with the package (masks/hg38_mask.bed.gz: UCSC hg38 gap
table, UCSC cytoBand acen/gvar/stalk bands and ENCODE ENCFF356LFX, chr1-22, X, Y).
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, Sequence

import numpy as np
import pandas as pd

from tad_vci.graded_evidence import combine_tier

MASK_DIR = Path(__file__).resolve().parent / "masks"
BUNDLED_MASKS = {"hg38": "hg38_mask.bed.gz"}
MASKED_COLUMNS = {"BD2": "BD2_ctcf", "BD3": "BD3_rad21"}
RECORD_COLUMNS = [
    "chrom", "pos", "mask_source", "BD1_caller_support",
    "BD2_ctcf_unmasked", "BD3_rad21_unmasked", "BD4_loop_anchor",
    "D_tier_unmasked", "D_tier", "not_assessable_criteria",
]


@dataclass(frozen=True)
class Mask:
    path: str
    n_intervals: int
    merged: dict          # chrom -> (starts, ends), sorted, overlapping/adjacent merged
    raw: dict             # chrom -> (starts, ends, labels) as read


def resolve_mask(spec: str, assembly: Optional[str] = None) -> str:
    """Return the file a `--mask` value names: a bundled mask name or a path.

    A bundled name takes precedence over a file of the same name in the working directory
    (pass ./hg38 for such a file). A bundled mask is tied to its assembly, so asking for
    it on a run declared as another assembly is refused."""
    key = str(spec).strip()
    if not key:
        raise ValueError("mask is empty")
    bundled = BUNDLED_MASKS.get(key.lower())
    if bundled is not None:
        if assembly and str(assembly).strip().lower() != key.lower():
            raise ValueError(
                f"mask '{key}' is the bundled {key.lower()} mask, but the run is declared "
                f"as assembly {assembly}")
        path = MASK_DIR / bundled
        if not path.is_file():
            raise FileNotFoundError(f"bundled mask missing from this install: {path}")
        return str(path)
    path = Path(key).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"mask BED not found: {key}")
    return str(path)


def _lines(path: str) -> Iterator[str]:
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    opener = gzip.open if gz else open
    with opener(path, "rt") as fh:
        yield from fh


def _merge(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.lexsort((ends, starts))
    s, e = starts[order], ends[order]
    out_s, out_e = [int(s[0])], [int(e[0])]
    for a, b in zip(s[1:], e[1:]):
        if a <= out_e[-1]:
            out_e[-1] = max(out_e[-1], int(b))
        else:
            out_s.append(int(a))
            out_e.append(int(b))
    return np.array(out_s, dtype=np.int64), np.array(out_e, dtype=np.int64)


def load_mask(path: str) -> Mask:
    """Read a BED mask (plain or gzip). Columns 1-3 are required; column 4, when present,
    names the source of the interval. Comment, `track` and `browser` lines are skipped,
    as is a header line whose coordinates are not integers."""
    default_label = Path(path).name
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    n = 0
    for lineno, line in enumerate(_lines(path), start=1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.rstrip("\r\n").split("\t") if "\t" in line else line.split()
        if len(f) < 3:
            raise ValueError(f"mask {path} line {lineno}: fewer than 3 columns")
        try:
            start, end = int(f[1]), int(f[2])
        except ValueError:
            continue
        if start < 0 or end <= start:
            raise ValueError(f"mask {path} line {lineno}: invalid interval {start}-{end}")
        chrom = f[0][3:] if f[0].startswith("chr") else f[0]
        label = f[3].strip() if len(f) > 3 and f[3].strip() else default_label
        by_chrom.setdefault(chrom, []).append((start, end, label))
        n += 1
    if n == 0:
        raise ValueError(f"mask {path} has no parseable intervals")
    merged, raw = {}, {}
    for chrom, ivs in by_chrom.items():
        s = np.array([i[0] for i in ivs], dtype=np.int64)
        e = np.array([i[1] for i in ivs], dtype=np.int64)
        raw[chrom] = (s, e, np.array([i[2] for i in ivs], dtype=object))
        merged[chrom] = _merge(s, e)
    return Mask(path=str(path), n_intervals=n, merged=merged, raw=raw)


def _norm_chrom(chrom: Sequence) -> np.ndarray:
    return np.array([c[3:] if c.startswith("chr") else c for c in map(str, chrom)], dtype=object)


def flag_positions(mask: Mask, chrom: Sequence, pos: Sequence) -> np.ndarray:
    """Boolean array: is each (chrom, pos) inside any mask interval, half-open."""
    ch = _norm_chrom(chrom)
    p = np.asarray(pos, dtype=np.int64)
    out = np.zeros(len(p), dtype=bool)
    for c in set(ch):
        iv = mask.merged.get(c)
        if iv is None:
            continue
        sel = ch == c
        starts, ends = iv
        i = np.searchsorted(starts, p[sel], side="right") - 1
        out[sel] = (i >= 0) & (p[sel] < ends[np.maximum(i, 0)])
    return out


def mask_sources(mask: Mask, chrom: str, pos: int) -> str:
    """'+'-joined, sorted labels of the mask intervals containing (chrom, pos)."""
    c = chrom[3:] if str(chrom).startswith("chr") else str(chrom)
    iv = mask.raw.get(c)
    if iv is None:
        return ""
    s, e, lab = iv
    hit = (s <= int(pos)) & (int(pos) < e)
    return "+".join(sorted(set(lab[hit])))


def apply_mask(table: pd.DataFrame, mask: Mask) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (masked copy of an evidence table, one record per masked row).

    Only BD2_ctcf, BD3_rad21, D_tier and not_assessable_criteria of masked rows change;
    every other cell, and every unmasked row, is returned as given."""
    hit = flag_positions(mask, table["chrom"], pd.to_numeric(table["pos"]))
    out = table.copy()
    if not hit.any():
        return out, pd.DataFrame(columns=RECORD_COLUMNS)
    for col in ("BD2_ctcf", "BD3_rad21", "D_tier", "not_assessable_criteria"):
        # A table read back with NA parsing can hold an all-empty grade column as float;
        # pandas 3 refuses to put strings into it.
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].astype(object)
    before = out.loc[hit, ["chrom", "pos", "BD1_caller_support", "BD2_ctcf", "BD3_rad21",
                           "BD4_loop_anchor", "D_tier"]].copy()
    for col in MASKED_COLUMNS.values():
        out.loc[hit, col] = "not_assessable"
    tiers, na = [], []
    for b1, b2, b3, b4 in zip(out.loc[hit, "BD1_caller_support"], out.loc[hit, "BD2_ctcf"],
                              out.loc[hit, "BD3_rad21"], out.loc[hit, "BD4_loop_anchor"]):
        tier, meta = combine_tier({"BD1": b1, "BD2": b2, "BD3": b3, "BD4": b4})
        tiers.append(tier)
        na.append(",".join(meta["not_assessable"]))
    out.loc[hit, "D_tier"] = tiers
    out.loc[hit, "not_assessable_criteria"] = na
    records = pd.DataFrame({
        "chrom": before["chrom"].to_numpy(),
        "pos": before["pos"].to_numpy(),
        "mask_source": [mask_sources(mask, c, p)
                        for c, p in zip(before["chrom"], pd.to_numeric(before["pos"]))],
        "BD1_caller_support": before["BD1_caller_support"].to_numpy(),
        "BD2_ctcf_unmasked": before["BD2_ctcf"].to_numpy(),
        "BD3_rad21_unmasked": before["BD3_rad21"].to_numpy(),
        "BD4_loop_anchor": before["BD4_loop_anchor"].to_numpy(),
        "D_tier_unmasked": before["D_tier"].to_numpy(),
        "D_tier": tiers,
        "not_assessable_criteria": na,
    }, columns=RECORD_COLUMNS)
    return out, records
