"""Shared helpers for deriving the v1.1 evidence tables (gap and blacklist masking).

The masking itself is the package's own `tad_vci.masking` (the code behind
`creditad annotate --mask`): a candidate whose `pos` lies inside a mask interval,
half-open [start, end), with chromosome names compared without a leading "chr", gets
BD2/BD3 = not_assessable and its tier re-derived with the released `combine_tier`. Nothing
is reimplemented here; these helpers only read the v1.0 packages, find the candidates
beyond their chromosome end, and translate mask labels into the v1.1 reason names.

Every location is an argument (see mask_and_retier.py); the defaults point into this
repository.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))
from tad_vci.masking import apply_mask, load_mask, resolve_mask  # noqa: E402,F401

SEED = 20260825

DEFAULT_PACKAGES_DIR = REPO / "data" / "multicell"
DEFAULT_CHROM_SIZES = REPO / "data" / "masks" / "hg38.chrom.sizes"
DEFAULT_MASK = "hg38"

PACKAGES = ["GM12878_25kb", "GM12878_10kb",
            "IMR90_25kb", "IMR90_10kb",
            "HepG2_25kb", "HepG2_10kb"]


def load_chrom_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    for line in Path(path).open():
        c, length = line.split()[:2]
        sizes[c.replace("chr", "")] = int(length)
    return sizes


def off_end(df: pd.DataFrame, sizes: dict[str, int]) -> pd.Series:
    """Candidates at or beyond the end of their chromosome (a BED-generation artefact).
    A chromosome missing from the sizes file counts as length 0, as in the v1.1 build."""
    return pd.Series([int(p) >= sizes.get(str(c), 0) for c, p in zip(df["chrom"], df["pos"])],
                     index=df.index)


def v11_reason(mask_source: str) -> str:
    """Translate the bundled mask's interval labels into the v1.1 reason names.

    The v1.1 tables call the union of UCSC gap rows and the centromeric/heterochromatic
    cytobands (acen, gvar, stalk) the assembly-gap set. Labels of any other mask are
    passed through unchanged."""
    labels = set(mask_source.split("+"))
    gap = any(lab == "assembly_gap" or lab.startswith("cytoband_") for lab in labels)
    bl = "encode_blacklist" in labels
    if gap and bl:
        return "assembly_gap+encode_blacklist"
    if gap:
        return "assembly_gap"
    if bl:
        return "encode_blacklist"
    return mask_source


def load_package(packages_dir: Path, pkg: str) -> pd.DataFrame:
    """Load a v1.0 annotation table with EVERY column kept as raw text, so unmasked rows
    are re-emitted byte for byte (no float parsing or re-formatting)."""
    base = Path(packages_dir) / pkg
    path = base / "annotation.tsv"
    if not path.is_file():
        path = base / "annotation.tsv.gz"
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)


def load_manifest(packages_dir: Path, pkg: str) -> dict:
    with (Path(packages_dir) / pkg / "manifest.json").open() as fh:
        return json.load(fh)
