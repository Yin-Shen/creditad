#!/usr/bin/env python3
"""Quick start: annotate the chr7 candidates of the GM12878 25 kb package from files that
ship in this repository, and check the result against the expected output byte for byte.

Inputs (all in the repository): the package's candidate BED and six caller BEDs, the chr7
preview ChIP tracks and the GM12878 HiCCUPS loop calls. The preview tracks are 250 bp
re-binned display copies, so the records here are NOT the released GM12878 25 kb records
(see data/chr7_preview_tracks/README.md); the released tables are regenerated exactly by
scripts/reproduce_packages.py with the genome-wide ENCODE tracks. This example checks that
an installation runs the documented command and produces the documented output.

    pip install .
    python examples/quickstart/run_quickstart.py          # exit 0 and PASS
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PKG = REPO / "data" / "multicell" / "GM12878_25kb"
TRACKS = REPO / "data" / "chr7_preview_tracks"
LOOPS = REPO / "data" / "loops" / "GM12878_HiCCUPS_loops_hg38.bedpe.gz"
EXPECTED = Path(__file__).resolve().parent / "expected" / "tadvci_annotation_beds_quickstart_chr7.tsv"
CALLERS = ["TopDom", "SpectralTAD", "OnTAD", "MSTD", "arrowhead", "DI"]


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def chr7_candidates(dst: Path) -> int:
    with gzip.open(PKG / "union_boundaries.bed.gz", "rt") as fh:
        rows = [ln for ln in fh if ln.split("\t", 1)[0] in ("chr7", "7")]
    dst.write_text("".join(rows))
    return len(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=None, help="output directory (default: a temporary one)")
    ap.add_argument("--write-expected", action="store_true", help="maintainers only: refresh expected/")
    args = ap.parse_args()
    out = args.out or Path(tempfile.mkdtemp(prefix="creditad_quickstart_"))
    out.mkdir(parents=True, exist_ok=True)
    cand = out / "GM12878_25kb_chr7_candidates.bed"
    n = chr7_candidates(cand)
    cmd = [sys.executable, "-m", "tad_vci.cli", "annotate", "--boundaries", str(cand)]
    for c in CALLERS:
        cmd += ["--panel", f"{c}={PKG / f'{c}_boundaries.bed.gz'}"]
    cmd += ["--ctcf", str(TRACKS / "GM12878_CTCF_ENCFF734CUT_hg38_chr7.bigWig"),
            "--rad21", str(TRACKS / "GM12878_RAD21_ENCFF571ZJJ_hg38_chr7.bigWig"),
            "--loops", str(LOOPS), "--cell-line", "GM12878", "--assembly", "hg38",
            "--resolution", "25000", "--tag", "quickstart_chr7", "--out", str(out)]
    print("$", " ".join(cmd[1:]))
    subprocess.run(cmd, check=True)
    table = out / "tadvci_annotation_beds_quickstart_chr7.tsv"
    if args.write_expected:
        EXPECTED.parent.mkdir(parents=True, exist_ok=True)
        EXPECTED.write_bytes(table.read_bytes())
        print(f"expected output refreshed: {EXPECTED} ({n} candidates)")
        return 0
    ok = sha256(table) == sha256(EXPECTED)
    print(f"{n} chr7 candidates annotated -> {table}")
    print(f"sha256 {sha256(table)}")
    print("PASS: identical to examples/quickstart/expected/" if ok else
          "FAIL: differs from examples/quickstart/expected/ (diff the two tables)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
