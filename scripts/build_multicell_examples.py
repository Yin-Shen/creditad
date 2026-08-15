#!/usr/bin/env python3
"""Build the CrediTAD demo packages: 3 cell lines x 2 resolutions, no consensus track.

Candidate set (the boundaries CrediTAD grades) = the deduplicated **union** of the six
callers' boundary positions. That is the honest candidate definition for a downstream
evidence-card tool: every position any caller proposed is put on trial, and BD1 records
how many of the six panel methods support it. No ConsensusTAD output is used, shipped,
or displayed.

Inputs  : $CREDITAD_CALLS_ROOT/<CELL>_<RES>/percaller.tsv  (upstream per-caller calls;
          see REPRODUCE.md -- not redistributed with this repository)
Outputs : creditad/example_data/multicell/<CELL>_<RES>/
            <caller>_boundaries.bed   (6 files, one per panel method)
            union_boundaries.bed      (candidate positions)
            annotation.tsv            (CrediTAD BD1-BD3 + D-tier)
            manifest.json
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

# Raw ENCODE/4DN inputs are NOT redistributed (licence + size). Point this at your own
# download directory; the layout below mirrors the accessions listed in REPRODUCE.md.
DATA = Path(os.environ.get("CREDITAD_DATA_ROOT", ROOT / "data" / "raw"))

SRC = Path(os.environ.get("CREDITAD_CALLS_ROOT", ROOT / "upstream_calls"))
OUT = ROOT / "example_data" / "multicell"

CALLERS = ["TopDom", "SpectralTAD", "OnTAD", "MSTD", "arrowhead", "DI"]

# Raw inputs are NOT redistributed with this repository (ENCODE/4DN licence + size:
# the six ChIP bigWigs alone are ~2.3 GB). Fetch them by the accessions below and point
# $CREDITAD_DATA_ROOT at the download directory; see REPRODUCE.md for the exact URLs.
# Loop calls live in $CREDITAD_LOOPS_DIR; they ship in this repository under data/loops/
# (~5 MB, gzipped), so the default already works.
LOOPS = Path(os.environ.get("CREDITAD_LOOPS_DIR", ROOT / "data" / "loops"))

CELLS = {
    "GM12878": {
        # 2026-07-25 unification: every demo cell line is hg38 with
        # fold-change-over-control ChIP and hg38 HiCCUPS loops.
        "assembly": "hg38",
        "mcool": str(DATA / "gm12878/hic/4DNFIXP4QG5B.mcool"),
        "ctcf": str(DATA / "gm12878/ctcf/ENCFF734CUT.bigWig"),
        "rad21": str(DATA / "gm12878/rad21/ENCFF571ZJJ.bigWig"),
        "ctcf_acc": "ENCFF734CUT",
        "rad21_acc": "ENCFF571ZJJ",
        "loops": str(LOOPS / "GM12878_HiCCUPS_loops_hg38.bedpe"),
        "loops_src": "HiCCUPS (GM12878, hg38)",
    },
    "IMR90": {
        "assembly": "hg38",
        "mcool": str(DATA / "imr90/hic/4DNFIJTOIGOI.mcool"),
        "ctcf": str(DATA / "imr90/ctcf/ENCFF105FHL.bigWig"),
        "rad21": str(DATA / "imr90/rad21/ENCFF048PZI.bigWig"),
        "ctcf_acc": "ENCFF105FHL",
        "rad21_acc": "ENCFF048PZI",
        "loops": str(LOOPS / "IMR90_HiCCUPS_loops_hg38.bedpe"),
        "loops_src": "HiCCUPS (IMR90, hg38)",
    },
    "HepG2": {
        "assembly": "hg38",
        "mcool": str(DATA / "hepg2/hic/4DNFIS6HAUPP.mcool"),
        "ctcf": str(DATA / "hepg2/ctcf/ENCFF357NFO.bigWig"),
        "rad21": str(DATA / "hepg2/rad21/ENCFF972ODZ.bigWig"),
        "ctcf_acc": "ENCFF357NFO",
        "rad21_acc": "ENCFF972ODZ",
        "loops": str(LOOPS / "HepG2_HiCCUPS_loops_hg38.bedpe"),
        "loops_src": "HiCCUPS (HepG2, hg38)",
    },
}

RESOLUTIONS = [25000, 10000]

# Package-specific residual DATA defects only. Method-level notes (tier≠predictor,
# published-parameter gaps, ChIP not re-BAM'd) live in docs/MEASURED_LIMITATIONS.md
# and must NOT be listed here — the Home page shows this list as "known caveats for
# this dataset", so empty means the delivered package has no open data defect.
#
# History (closed):
# - GM12878 Arrowhead chrom gaps: closed after SCALE/96g rebuild (all 23 chroms present).
# - IMR90_10kb Arrowhead chr7: closed 2026-07-26 by replacing KR-only run with juicer
#   SCALE output (out10000_scale; chr7=421 domains). KR still fails on that chrom.
# - HepG2_10kb DI +1 bin (10 kb) vs TopDom/SpectralTAD: closed 2026-07-26 by shifting
#   DI domains +10 kb in calling_repro_final/HepG2_10kb/percaller.tsv (was documented as
#   +2 bin / 20 kb on an earlier build; remeasured mode offset was +1 bin).
KNOWN_CAVEATS: dict[str, list[str]] = {}
GLOBAL_CAVEATS: list[str] = []


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def boundaries_of(df: pd.DataFrame, res: int) -> set[tuple[str, int]]:
    """Domain endpoints snapped to the bin grid, 0-based BED starts."""
    out: set[tuple[str, int]] = set()
    for c, s, e in zip(df["chr"], df["start"].astype(int), df["end"].astype(int)):
        for p in (s - 1, e):                 # one-based start -> 0-based; end is a boundary too
            out.add((str(c), (p // res) * res))
    return out


def write_bed(path: Path, positions: set[tuple[str, int]], res: int) -> int:
    def key(x):
        c = x[0][3:] if x[0].startswith("chr") else x[0]
        return ((0, int(c)) if c.isdigit() else (1, c), x[1])

    rows = sorted(positions, key=key)
    with path.open("w") as f:
        for c, p in rows:
            f.write(f"{c}\t{p}\t{p + res}\n")
    return len(rows)


def build(cell: str, res: int) -> dict:
    combo = f"{cell}_{res // 1000}kb"
    src = SRC / combo / "percaller.tsv"
    if not src.is_file():
        raise SystemExit(f"missing input: {src}")
    pc = pd.read_csv(src, sep="\t")
    d = OUT / combo
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)

    per_caller_counts, method_beds = {}, []
    union: set[tuple[str, int]] = set()
    for caller in CALLERS:
        sub = pc[pc["meta.tool"] == caller]
        pos = boundaries_of(sub, res) if len(sub) else set()
        bed = d / f"{caller}_boundaries.bed"
        per_caller_counts[caller] = write_bed(bed, pos, res)
        union |= pos
        if per_caller_counts[caller]:
            method_beds.append((str(bed), caller))

    union_bed = d / "union_boundaries.bed"
    n_union = write_bed(union_bed, union, res)

    meta = CELLS[cell]
    from tad_vci.annotate_from_beds import annotate_from_beds

    result = annotate_from_beds(
        str(union_bed),
        method_beds=method_beds,
        resolution=res,
        # BD1 vote matching: +/-1 bin (2026-07-29). At 2 bins one caller call could be
        # counted toward up to 5 neighbouring candidates. BD4 keeps its own published
        # 50 kb window and no longer scales with the resolution.
        tol=res * 1,
        loop_tol=50000,
        ctcf_bw=meta["ctcf"],
        rad21_bw=meta["rad21"],
        loops_path=meta.get("loops"),
        cell_line=cell,
        assembly=meta["assembly"],
        hash_tag=f"{combo}_union",
        origin="builtin_demo",
        output_dir=str(d / "_work"),
    )
    ann_src = Path(result["TADVCI_annotation"])
    if not ann_src.is_file():
        raise SystemExit(f"{combo}: annotate_from_beds produced no annotation; keys={list(result)}")
    shutil.copy(ann_src, d / "annotation.tsv")
    tier_src = Path(result.get("TADVCI_tier", ""))
    if tier_src.is_file():
        shutil.copy(tier_src, d / "tier_boundaries.bed")
    shutil.rmtree(d / "_work", ignore_errors=True)

    ann = pd.read_csv(d / "annotation.tsv", sep="\t")
    tiers = ann["D_tier"].value_counts().to_dict() if "D_tier" in ann else {}
    manifest = {
        "name": f"{cell}_{res // 1000}kb_multicaller_union",
        "cell_line": cell,
        "assembly": meta["assembly"],
        "resolution_bp": res,
        "role": "CrediTAD downstream evidence-card demo — boundaries NOT called by CrediTAD",
        "candidate_definition": ("deduplicated union of the six panel callers' domain endpoints, "
                                 "snapped to the bin grid; no ConsensusTAD involved"),
        "consensus_used": False,
        "upstream_pipeline": ("local multi-caller run (OnTAD level-2, SpectralTAD level-1, TopDom, "
                              "MSTD v1, Arrowhead Juicer 1.22.01 with SCALE balancing, DI domaincaller); "
                              "cooler ICE for Hi-C. Parameter CHOICES follow Li et al. Nat Methods 2026 "
                              "where obtainable; this is a local panel, not a byte-reproduction of TADShop. "
                              "Method-level limitations: docs/MEASURED_LIMITATIONS.md."),
        # Host paths are not provenance a reader can act on; record the pipeline stage
        # and the file identity instead so a clone can be checked without this machine.
        "source_percaller": f"calling_repro_final/{combo}/percaller.tsv",
        "source_percaller_sha256": _sha256(src),
        "n_candidates": int(n_union),
        "per_caller_boundary_counts": per_caller_counts,
        "bd1_panel": [n for _, n in method_beds],
        "bd1_n_methods": len(method_beds),
        "chip_tracks": {"CTCF": meta["ctcf_acc"], "RAD21": meta["rad21_acc"]},
        "chip_output_type": "fold change over control",
        "cross_cell_comparable": True,
        "loop_calls": meta.get("loops_src"),
        "tier_counts": tiers,
        "caveats": KNOWN_CAVEATS.get(combo, []) + GLOBAL_CAVEATS,
        "reference_note": ("TAD boundaries have no ground truth; this package is caller output "
                           "graded by evidence, not a validated boundary set."),
    }
    (d / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{combo}: candidates={n_union} panel={len(method_beds)} tiers={tiers}")
    return manifest


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    index = []
    only = sys.argv[1:] or None
    for cell in CELLS:
        for res in RESOLUTIONS:
            combo = f"{cell}_{res // 1000}kb"
            if only and combo not in only:
                continue
            index.append(build(cell, res))
    if index:
        (OUT / "index.json").write_text(json.dumps(
            {"datasets": [{k: m[k] for k in ("name", "cell_line", "assembly", "resolution_bp",
                                             "n_candidates", "bd1_n_methods", "consensus_used")}
                          for m in index]}, indent=2) + "\n")
    print(f"\nwrote {len(index)} demo packages under {OUT}")


if __name__ == "__main__":
    main()
