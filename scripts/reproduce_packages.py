#!/usr/bin/env python3
"""reproduce_packages.py -- regenerate the six delivered CrediTAD evidence tables and
verify them against the shipped ones, column by column.

WHAT THIS REPRODUCES
    The six delivered annotation tables (GM12878 / IMR90 / HepG2 x 25 kb / 10 kb;
    284,744 candidate boundaries in total), hg38, tier rule max_support_v2 -- i.e. exactly
    the tables the manuscript reports on. It re-runs the annotation path from each
    package's own candidate BED and six caller BEDs plus the molecular inputs, and
    compares EVERY column of the result against the shipped annotation.tsv.

WHAT IT DOES NOT DO
    It does not re-call TADs. The per-caller boundary BEDs are inputs, produced upstream
    by six external callers; regenerating those requires the raw Hi-C contact matrices and
    the callers themselves, and is out of scope here (see REPRODUCE.md).

INPUTS YOU MUST FETCH YOURSELF (not redistributable: ENCODE/4DN terms + size)
    $CREDITAD_DATA_ROOT/<cell>/ctcf/<ACC>.bigWig      ENCODE CTCF fold-change-over-control
    $CREDITAD_DATA_ROOT/<cell>/rad21/<ACC>.bigWig     ENCODE RAD21 fold-change-over-control
        GM12878  CTCF ENCFF734CUT   RAD21 ENCFF571ZJJ
        IMR90    CTCF ENCFF105FHL   RAD21 ENCFF048PZI
        HepG2    CTCF ENCFF357NFO   RAD21 ENCFF972ODZ
    $CREDITAD_LOOPS_DIR/<CELL>_HiCCUPS_loops_hg38.bedpe   HiCCUPS loop calls (Zenodo)
    $CREDITAD_PACKAGES_DIR/<CELL>_<RES>kb/   the delivered packages -- these SHIP in this
        repository under data/multicell/ (gzipped), so this default already works.
    $CREDITAD_LOOPS_DIR                      HiCCUPS loop calls -- also shipped, under
        data/loops/ (gzipped).
    Only the two ChIP bigWigs per cell line must be fetched from ENCODE.
    Both plain and .gz inputs are accepted.

IMPORTANT -- USE THE GENOME-WIDE ChIP TRACKS.
    The chr7-only bigWigs distributed for browser display are 250 bp re-binned previews.
    They are NOT value-identical to the genome-wide tracks (measured: median window ratio
    1.12 for CTCF, 1.11 for RAD21; 0 of 1,528 chr7 windows agree to within 1e-6) and they
    do NOT reproduce the delivered grades. Point CREDITAD_DATA_ROOT at the full tracks.

Exit code 0 iff every column of every available package matches.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

CALLERS = ["TopDom", "SpectralTAD", "OnTAD", "MSTD", "arrowhead", "DI"]
PACKAGES = [("GM12878", 25000), ("GM12878", 10000), ("IMR90", 25000),
            ("IMR90", 10000), ("HepG2", 25000), ("HepG2", 10000)]
ACC = {
    "GM12878": {"ctcf": "ENCFF734CUT", "rad21": "ENCFF571ZJJ"},
    "IMR90": {"ctcf": "ENCFF105FHL", "rad21": "ENCFF048PZI"},
    "HepG2": {"ctcf": "ENCFF357NFO", "rad21": "ENCFF972ODZ"},
}
EXPECTED_ROWS = {"GM12878_25kb": 28259, "GM12878_10kb": 69479,
                 "IMR90_25kb": 28476, "IMR90_10kb": 65157,
                 "HepG2_25kb": 29513, "HepG2_10kb": 63860}


def paths():
    """Defaults point at the tier-1 data bundled in this repository."""
    return (Path(os.environ.get("CREDITAD_PACKAGES_DIR", REPO / "data" / "multicell")),
            Path(os.environ.get("CREDITAD_DATA_ROOT", REPO / "data" / "raw")),
            Path(os.environ.get("CREDITAD_LOOPS_DIR", REPO / "data" / "loops")))


def resolve(base: Path, name: str) -> Path | None:
    """Return `name` or `name.gz` under `base`, whichever exists."""
    for cand in (base / name, base / (name + ".gz")):
        if cand.exists():
            return cand
    return None


def materialise(path: Path, tmp: Path) -> str:
    """Give the engine a plain-text path: decompress .gz into `tmp` if needed."""
    if path.suffix != ".gz":
        return str(path)
    out = tmp / path.name[:-3]
    if not out.exists():
        tmp.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "rb") as fin, out.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
    return str(out)


def compare(fresh: pd.DataFrame, shipped: pd.DataFrame) -> tuple[int, int, list]:
    shared = [c for c in fresh.columns if c in shipped.columns]
    bad = []
    for c in shared:
        a, b = fresh[c], shipped[c]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            eq = np.isclose(a.astype(float), b.astype(float),
                            rtol=1e-9, atol=1e-12, equal_nan=True)
        else:
            eq = (a.fillna("").astype(str).to_numpy()
                  == b.fillna("").astype(str).to_numpy())
        if not eq.all():
            bad.append({"column": c, "n_mismatch": int((~eq).sum()),
                        "agreement": float(eq.mean())})
    return len(shared), len(shared) - len(bad), bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(REPO / "reproduce_out"))
    ap.add_argument("--packages", nargs="*", default=None,
                    help="subset, e.g. GM12878_25kb (default: all six)")
    args = ap.parse_args()

    PKGS, DATA, LOOPS = paths()
    from tad_vci.annotate_from_beds import annotate_from_beds

    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)
    report, skipped, failures = [], [], []
    for cell, res in PACKAGES:
        name = f"{cell}_{res // 1000}kb"
        if args.packages and name not in args.packages:
            continue
        pkg = PKGS / name
        ctcf = DATA / cell.lower() / "ctcf" / f"{ACC[cell]['ctcf']}.bigWig"
        rad21 = DATA / cell.lower() / "rad21" / f"{ACC[cell]['rad21']}.bigWig"
        shipped_tsv = resolve(pkg, "annotation.tsv")
        union_bed = resolve(pkg, "union_boundaries.bed")
        loops = resolve(LOOPS, f"{cell}_HiCCUPS_loops_hg38.bedpe")
        caller_beds = {c: resolve(pkg, f"{c}_boundaries.bed") for c in CALLERS}
        missing = [str(x) for x in (
            shipped_tsv or pkg / "annotation.tsv(.gz)",
            union_bed or pkg / "union_boundaries.bed(.gz)",
            loops or LOOPS / f"{cell}_HiCCUPS_loops_hg38.bedpe(.gz)",
            ctcf, rad21) if x is None or not Path(x).exists()]
        missing += [f"{pkg}/{c}_boundaries.bed(.gz)"
                    for c, v in caller_beds.items() if v is None]
        if missing:
            skipped.append({"package": name, "missing": missing})
            print(f"SKIP {name}: missing {len(missing)} input(s); first = {missing[0]}")
            continue

        tmp = Path(tempfile.mkdtemp(prefix=f"creditad_{name}_"))
        panel = [(materialise(caller_beds[c], tmp), c) for c in CALLERS]
        r = annotate_from_beds(
            materialise(union_bed, tmp), method_beds=panel, resolution=res,
            tol=res, loop_tol=50000, ctcf_bw=str(ctcf), rad21_bw=str(rad21),
            loops_path=materialise(loops, tmp), cell_line=cell, assembly="hg38",
            hash_tag=name, output_dir=str(outdir / name), origin="builtin_demo",
        )
        fresh = pd.read_csv(r["TADVCI_annotation"], sep="\t", low_memory=False)
        shipped = pd.read_csv(shipped_tsv, sep="\t", low_memory=False)  # pandas reads .gz
        if len(fresh) != len(shipped):
            failures.append({"package": name, "reason":
                             f"row count {len(fresh)} != shipped {len(shipped)}"})
            print(f"FAIL {name}: row count {len(fresh)} != {len(shipped)}")
            continue
        n_shared, n_ok, bad = compare(fresh, shipped)
        exp = EXPECTED_ROWS.get(name)
        row_ok = exp is None or len(fresh) == exp
        ok = not bad and row_ok
        report.append({"package": name, "n_rows": len(fresh),
                       "n_rows_expected": exp, "columns_compared": n_shared,
                       "columns_identical": n_ok, "mismatched_columns": bad, "pass": ok})
        print(f"{'PASS' if ok else 'FAIL'} {name}: {n_ok}/{n_shared} columns identical, "
              f"{len(fresh)} rows")
        if not ok:
            failures.append({"package": name, "bad_columns": bad})

    total_rows = sum(x["n_rows"] for x in report)
    summary = {"packages_checked": len(report), "packages_skipped": skipped,
               "total_rows_compared": total_rows, "failures": failures,
               "all_pass": bool(report) and not failures}
    (outdir / "reproduction_report.json").write_text(json.dumps(
        {"summary": summary, "per_package": report}, indent=2))
    print("=" * 70)
    if not report:
        print("NOTHING CHECKED - no package had all its inputs present.")
        print("See the docstring above and REPRODUCE.md for what to fetch.")
        return 2
    print(f"{len(report)} package(s), {total_rows} rows compared")
    if len(report) == 6 and total_rows != 284744:
        print(f"WARNING: all six packages present but total rows {total_rows} != 284744")
    print("REPRODUCTION: " + ("PASS" if summary["all_pass"] else "FAIL"))
    print(f"report -> {outdir / 'reproduction_report.json'}")
    return 0 if summary["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
