#!/usr/bin/env python3
"""Command-line interface for CrediTAD evidence grading.

Grades already-called TAD boundaries and writes the evidence-card table plus the
tier BED. This is the same code path the HTTP API uses (`tad_vci.annotate_from_beds`),
so a scripted run and a GUI run on identical inputs produce identical output.

Examples
--------
Grade one caller's boundaries against a five-method panel, with CTCF/RAD21 and loops:

    creditad annotate \\
        --boundaries candidates.bed \\
        --panel TopDom=TopDom.bed --panel OnTAD=OnTAD.bed --panel DI=DI.bed \\
        --ctcf CTCF.bigWig --rad21 RAD21.bigWig --loops loops.bedpe \\
        --cell-line GM12878 --assembly hg19 --resolution 25000 \\
        --out results/

Grade a single caller with no panel (BD1 is floored to `weak`, which the card states):

    creditad annotate --boundaries mycalls.bed --resolution 10000 --out results/

Report the grading rule that produced a run:

    creditad criteria
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# The package lives beside the rest of the backend; make sibling modules importable
# whether this is run from an installed wheel or straight from the source tree.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def _parse_panel(values: list[str]) -> list[tuple[str, str]]:
    """`--panel NAME=path.bed` (repeatable) -> [(path, name), ...] for annotate_from_beds."""
    out: list[tuple[str, str]] = []
    for raw in values or []:
        if "=" in raw:
            name, path = raw.split("=", 1)
        else:
            path = raw
            name = Path(raw).stem.replace("_boundaries", "")
        p = Path(path).expanduser()
        if not p.is_file():
            raise SystemExit(f"--panel file not found: {p}")
        out.append((str(p), name))
    return out


def cmd_annotate(args: argparse.Namespace) -> int:
    from tad_vci.annotate_from_beds import annotate_from_beds

    primary = Path(args.boundaries).expanduser()
    if not primary.is_file():
        raise SystemExit(f"--boundaries file not found: {primary}")
    for label, val in (("--ctcf", args.ctcf), ("--rad21", args.rad21), ("--loops", args.loops)):
        if val and not Path(val).expanduser().is_file():
            raise SystemExit(f"{label} file not found: {val}")

    panel = _parse_panel(args.panel)
    result = annotate_from_beds(
        str(primary),
        method_beds=panel,
        resolution=args.resolution,
        tol=args.vote_tolerance,   # None -> annotate_from_beds derives resolution*1
        ctcf_bw=str(Path(args.ctcf).expanduser()) if args.ctcf else None,
        rad21_bw=str(Path(args.rad21).expanduser()) if args.rad21 else None,
        loops_path=str(Path(args.loops).expanduser()) if args.loops else None,
        cell_line=args.cell_line,
        assembly=args.assembly,
        hash_tag=args.tag,
        output_dir=str(Path(args.out).expanduser()),
        origin="user",
        uncalibrated_bands=args.uncalibrated_bands,
    )
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"boundaries graded : {result['n_boundaries']}")
        print(f"BD1 panel (n={result['n_methods']}) : {', '.join(result['bd1_panel']) or '(single method)'}")
        print(f"tier distribution : {result['tier_dist']}")
        print(f"evidence table    : {result['TADVCI_annotation']}")
        print(f"tier BED          : {result['TADVCI_tier']}")
        if not args.ctcf:
            print("note: no --ctcf given, so BD2 is not_assessable (distinct from 'depleted')")
        if not args.rad21:
            print("note: no --rad21 given, so BD3 is not_assessable")
        if not args.loops:
            print("note: no --loops given, so BD4 is not_assessable")
        if not args.cell_line:
            print("note: no --cell-line given, so BD2/BD3 are not treated as calibrated; "
                  f"bands policy = {args.uncalibrated_bands}")
    return 0


def cmd_criteria(args: argparse.Namespace) -> int:
    from tad_vci.graded_evidence import CRITERIA

    if args.json:
        print(json.dumps({k: {
            "name": c.name, "assay": c.assay, "availability": c.availability,
            "role": c.role, "grading": c.grading,
        } for k, c in CRITERIA.items()}, indent=2))
        return 0
    for key, c in CRITERIA.items():
        print(f"{key}  {c.name}  [{c.availability} / {c.role}]")
        print(f"    assay   : {c.assay}")
        print(f"    grading : {c.grading}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="creditad",
        description="Evidence cards and tiers for already-called TAD boundaries. "
                    "CrediTAD does not call TADs and does not rank callers.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = p.add_subparsers(dest="command", required=True)

    a = sub.add_parser("annotate", help="grade a boundary BED and write the evidence table + tier BED")
    a.add_argument("--boundaries", required=True,
                   help="BED of candidate boundary positions to grade (column 2 is the position)")
    a.add_argument("--panel", action="append", default=[], metavar="NAME=FILE.bed",
                   help="a caller in the BD1 voting panel; repeat once per caller. "
                        "With no panel, BD1 is floored to 'weak'.")
    a.add_argument("--ctcf", help="CTCF bigWig for BD2 (omit -> BD2 not_assessable)")
    a.add_argument("--rad21", help="RAD21/cohesin bigWig for BD3 (omit -> BD3 not_assessable)")
    a.add_argument("--loops", help="loop calls (.bedpe or HiCCUPS .txt) for BD4 (omit -> BD4 not_assessable)")
    a.add_argument("--resolution", type=int, default=25000, help="bin size in bp (default 25000)")
    a.add_argument("--vote-tolerance", type=int, default=None, metavar="BP",
                   help="BD1 vote matching tolerance in bp (default: 1 x resolution)")
    a.add_argument("--cell-line", default=None,
                   help="cell line, used to pick the packaged BD2/BD3 bands. Unknown or "
                        "omitted is NOT treated as calibrated; see --uncalibrated-bands "
                        "for what happens instead (default: estimate from your own track)")
    a.add_argument("--assembly", default=None, help="hg19 or hg38, recorded in the output provenance")
    a.add_argument("--out", default="./creditad_results", help="output directory")
    a.add_argument("--tag", default="cli", help="suffix for the output filenames")
    a.add_argument("--uncalibrated-bands", default="in_situ",
                   choices=["in_situ", "not_assessable", "gm12878_fallback"],
                   help="what to do when the cell line has no packaged BD2/BD3 calibration. "
                        "in_situ (default): estimate this track's own p75/p90 and use its own "
                        "background. not_assessable: refuse to grade BD2/BD3. gm12878_fallback: "
                        "the legacy behaviour, which inflated grades one-directionally "
                        "(measured BD3 +31%% on HepG2) — opt in only if you know why.")
    a.add_argument("--json", action="store_true", help="print the result summary as JSON")
    a.set_defaults(func=cmd_annotate)

    c = sub.add_parser("criteria", help="print the BD1-BD4 grading rules actually compiled in")
    c.add_argument("--json", action="store_true")
    c.set_defaults(func=cmd_criteria)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
