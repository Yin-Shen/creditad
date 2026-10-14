"""Build the v1.1 multicell packages (gap and blacklist masking) from the v1.0 packages.

WHAT CHANGES relative to v1.0
  (a) Every candidate whose `pos` falls inside the mask -- by default the bundled hg38 mask
      (UCSC hg38 gap table, UCSC cytoBand acen/gvar/stalk bands, ENCODE hg38 blacklist
      ENCFF356LFX) -- has BD2_ctcf and BD3_rad21 re-graded to `not_assessable`. At such a
      position the ChIP readout is not a measured absence, it is an absence of measurement.
      BD1 and BD4 are NOT masked:
        BD1 is caller agreement over the Hi-C contact map, a property of the callers'
        output at that coordinate, not of the ChIP mappability at it;
        BD4 counts loop-call anchors within +/-50 kb, and a loop call whose anchor lands
        in a gap is still a call that exists in the supplied loop file.
      The masking is tad_vci.masking.apply_mask, the code behind
      `creditad annotate --mask hg38`; it is not reimplemented here.
  (b) The candidates whose `pos` is at or beyond the end of their chromosome are DELETED
      (BED-generation artefact; 12 in the six packages); they are listed in
      <out>/impact/removed_off_chromosome_end.csv.
  (c) D_tier and not_assessable_criteria are recomputed for masked rows ONLY, with the
      released backend/tad_vci/graded_evidence.combine_tier. Every other row is passed
      through byte-for-byte.
  (d) union_boundaries.bed / tier_boundaries.bed are regenerated (derived from
      annotation.tsv, one line per candidate). The per-caller BEDs are not copied: they are
      upstream caller output, unchanged from v1.0, and are read from data/multicell/.

RETAINED ON PURPOSE: the raw ctcf / rad21 / *_fold / *_null_pct columns keep their v1.0
values on masked rows. The 38-column schema is fixed, and a reader is better served by
seeing the number that was withheld than by a blank.

Usage:
    python scripts/mask_and_retier.py --out DIR [--packages-dir data/multicell]
        [--mask hg38|MASK.bed] [--chrom-sizes data/masks/hg38.chrom.sizes]
        [--check data/multicell_v1.1] [PKG ...]          (default: all six packages)

With --check, every file written is compared with the released v1.1 package of the same
name (annotation.tsv, union_boundaries.bed, tier_boundaries.bed, masking_records.tsv byte
for byte after decompression; the three summary tables in impact/; manifest.json field by
field). The exit code is 0 iff every table is identical and every manifest field this
script computes matches; manifest fields that name file locations or the generating
script, and fields copied unchanged from the v1.0 manifests, are listed when they differ.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mask_common as M  # noqa: E402

MASK_CRITERIA = ["BD2", "BD3"]
PROVENANCE_FIELDS = ("reference_files", "engine_caveat", "generator")
COMPUTED_FIELDS = ("data_version", "n_candidates", "tier_counts", "supersedes")

MASK_NOTE = (
    "v1.1 mappability masking. A candidate whose bin start falls "
    "inside an hg38 assembly gap or inside the ENCODE hg38 blacklist has its CTCF (BD2) "
    "and cohesin (BD3) grades set to not_assessable, because at such a position the ChIP "
    "readout is unmeasurable rather than measured-negative. BD1 (cross-caller support) "
    "and BD4 (loop anchors) are left as graded: neither depends on ChIP mappability at "
    "the candidate coordinate. D_tier and not_assessable_criteria are recomputed for "
    "masked rows with the shipped combine_tier; all other rows are byte-identical to "
    "v1.0. The raw ctcf/rad21/fold/null-percentile columns keep their v1.0 values so the "
    "withheld measurement stays inspectable; the grade, not the number, is what is "
    "withheld."
)


def _display(path: Path) -> str:
    p = Path(path).resolve()
    try:
        return p.relative_to(M.REPO).as_posix()
    except ValueError:
        return str(path)


def union_bytes(chrom, pos, res: int) -> bytes:
    buf = io.StringIO()
    for c, p in zip(chrom, pos):
        buf.write(f"chr{c}\t{p}\t{p + res}\n")
    return buf.getvalue().encode()


def tier_bytes(chrom, pos, tier, res: int) -> bytes:
    buf = io.StringIO()
    for c, p, t in zip(chrom, pos, tier):
        buf.write(f"{c}\t{p}\t{p + res}\t{t}\n")
    return buf.getvalue().encode()


def build_manifest(pkg: str, man: dict, stats: dict, tier_counts: dict, n_rows: int,
                   refs: dict) -> dict:
    out: dict = {}
    for k, v in man.items():
        out[k] = v
        if k == "name":
            out["data_version"] = "1.1"
    out["n_candidates"] = n_rows
    out["tier_counts"] = tier_counts
    out["masking"] = {
        "applied": True,
        "description": MASK_NOTE,
        "predicate": "point-in-interval on the candidate `pos`, half-open [start, end); "
                     "identical to the gap and blacklist predicate of r2_04_unmappable_and_windows.py",
        "reference_files": refs,
        "criteria_masked": MASK_CRITERIA,
        "criteria_left_graded": ["BD1", "BD4"],
        "criteria_left_graded_reason":
            "BD1 is caller agreement on the Hi-C map and BD4 is loop-anchor overlap in the "
            "supplied loop file; neither is a ChIP measurement at the candidate coordinate, "
            "so neither is invalidated by low mappability there.",
        "n_rows_masked": stats["n_masked"],
        "n_rows_masked_assembly_gap_only": stats["n_gap_only"],
        "n_rows_masked_blacklist_only": stats["n_bl_only"],
        "n_rows_masked_both": stats["n_both"],
        "n_rows_removed_off_chromosome_end": stats["n_removed"],
        "raw_signal_columns_retained": True,
        "engine_caveat":
            "Masking is applied after graded_evidence.annotate(), which has no mappability "
            "model of its own. `creditad annotate --mask hg38` applies the same masking to "
            "a new annotation run; scripts/mask_and_retier.py derives these tables from the "
            "v1.0 tables.",
        "records_file": "masking_records.tsv",
        "generator": "scripts/mask_and_retier.py",
        "seed": M.SEED,
    }
    out["supersedes"] = {
        "data_version": "1.0",
        "package_dir": f"example_data/multicell/{pkg}",
        "n_candidates": man["n_candidates"],
        "tier_counts": man["tier_counts"],
    }
    return out


def process(pkg: str, args, mask, sizes: dict, refs: dict) -> dict:
    df = M.load_package(args.packages_dir, pkg)
    cols = list(df.columns)
    res = int(df["resolution_bp"].iloc[0])

    # (a)+(c) mask BD2/BD3 and re-derive the tier of masked rows (package code)
    masked, rec = M.apply_mask(df, mask)
    assert list(masked.columns) == cols, "38-column schema changed"

    # (b) delete off-chromosome-end rows (from the table and from the mask records)
    gone = M.off_end(df, sizes)
    removed = df[gone].copy()
    removed.insert(0, "package", pkg)
    removed_out = removed[["package", "chrom", "pos", "votes", "BD1_caller_support",
                           "BD2_ctcf", "BD3_rad21", "BD4_loop_anchor", "D_tier"]].copy()
    removed_out["chrom_size_hg38"] = [sizes.get(str(c), -1) for c in removed["chrom"]]
    removed_out["reason"] = "pos_at_or_beyond_chromosome_end"
    out = masked[~gone.to_numpy()].reset_index(drop=True)
    gone_keys = set(zip(removed["chrom"], removed["pos"]))
    rec = rec[[k not in gone_keys for k in zip(rec["chrom"], rec["pos"])]].reset_index(drop=True)
    ipos = out["pos"].astype("int64")

    d = args.out / pkg
    d.mkdir(parents=True, exist_ok=True)
    (d / "annotation.tsv").write_bytes(
        out.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    (d / "union_boundaries.bed").write_bytes(union_bytes(out["chrom"], ipos, res))
    (d / "tier_boundaries.bed").write_bytes(tier_bytes(out["chrom"], ipos, out["D_tier"], res))

    reasons = [M.v11_reason(s) for s in rec["mask_source"]]
    records = pd.DataFrame({
        "chrom": rec["chrom"].to_numpy(),
        "pos": rec["pos"].astype("int64").to_numpy(),
        "reason": reasons,
        "BD1_caller_support": rec["BD1_caller_support"].to_numpy(),
        "BD2_ctcf_v10": rec["BD2_ctcf_unmasked"].to_numpy(),
        "BD3_rad21_v10": rec["BD3_rad21_unmasked"].to_numpy(),
        "BD4_loop_anchor": rec["BD4_loop_anchor"].to_numpy(),
        "D_tier_v10": rec["D_tier_unmasked"].to_numpy(),
        "D_tier_v11": rec["D_tier"].to_numpy(),
        "not_assessable_criteria_v11": rec["not_assessable_criteria"].to_numpy(),
    })
    records.to_csv(d / "masking_records.tsv", sep="\t", index=False)

    tier_counts = dict(sorted(Counter(out["D_tier"]).items(), key=lambda kv: -kv[1]))
    stats = dict(
        package=pkg,
        n_v10=int(len(df)),
        n_removed=int(len(removed)),
        n_v11=int(len(out)),
        n_masked=int(len(records)),
        n_gap_only=reasons.count("assembly_gap"),
        n_bl_only=reasons.count("encode_blacklist"),
        n_both=reasons.count("assembly_gap+encode_blacklist"),
    )
    man = build_manifest(pkg, M.load_manifest(args.packages_dir, pkg), stats, tier_counts,
                         len(out), refs)
    with (d / "manifest.json").open("w") as fh:
        json.dump(man, fh, indent=2)
        fh.write("\n")

    stats["tier_counts_v11"] = json.dumps(tier_counts)
    return {"stats": stats, "records": records.assign(package=pkg), "removed": removed_out}


def _ref_bytes(ref_dir: Path, name: str) -> bytes | None:
    for cand, gz in ((ref_dir / name, False), (ref_dir / (name + ".gz"), True)):
        if cand.is_file():
            raw = cand.read_bytes()
            return gzip.decompress(raw) if gz else raw
    return None


def _first_diff(a: bytes, b: bytes) -> str:
    la, lb = a.splitlines(), b.splitlines()
    for i, (x, y) in enumerate(zip(la, lb), start=1):
        if x != y:
            return f"line {i}"
    return f"line count {len(la)} vs {len(lb)}"


def check(args, pkgs: list[str]) -> bool:
    ok = True
    files = ["annotation.tsv", "union_boundaries.bed", "tier_boundaries.bed",
             "masking_records.tsv"]
    for pkg in pkgs:
        for name in files:
            new = (args.out / pkg / name).read_bytes()
            ref = _ref_bytes(args.check / pkg, name)
            if ref is None:
                print(f"[check] {pkg}/{name}: no released file to compare")
                ok = False
            elif new == ref:
                print(f"[check] {pkg}/{name}: identical ({len(new.splitlines())} lines)")
            else:
                print(f"[check] {pkg}/{name}: DIFFERS at {_first_diff(new, ref)}")
                ok = False
        new_m = json.loads((args.out / pkg / "manifest.json").read_text())
        ref_m = json.loads((args.check / pkg / "manifest.json").read_text())
        diff = sorted(k for k in set(new_m) | set(ref_m)
                      if k != "masking" and new_m.get(k) != ref_m.get(k))
        mdiff = sorted(k for k in set(new_m["masking"]) | set(ref_m["masking"])
                       if new_m["masking"].get(k) != ref_m["masking"].get(k))
        # Fields this script writes must match; the rest are copied from the v1.0 manifest
        # in --packages-dir and are reported, not judged.
        unexpected = ([k for k in diff if k in COMPUTED_FIELDS]
                      + [f"masking.{k}" for k in mdiff if k not in PROVENANCE_FIELDS])
        named = [f"masking.{k}" for k in mdiff if k in PROVENANCE_FIELDS]
        copied = [k for k in diff if k not in COMPUTED_FIELDS]
        print(f"[check] {pkg}/manifest.json: "
              + ("identical apart from " + ", ".join(named) if named else "identical")
              + (f"; copied from the v1.0 manifest and different: {', '.join(copied)}"
                 if copied else "")
              + (f"; DIFFERS in {', '.join(unexpected)}" if unexpected else ""))
        ok = ok and not unexpected
    if len(pkgs) == len(M.PACKAGES):
        for name in ("mask_package_stats.csv", "masking_records_all.csv",
                     "removed_off_chromosome_end.csv"):
            new = (args.out / "impact" / name).read_bytes()
            ref = _ref_bytes(args.check / "impact", name)
            same = ref is not None and new == ref
            print(f"[check] impact/{name}: " + ("identical" if same else
                  "no released file" if ref is None else f"DIFFERS at {_first_diff(new, ref)}"))
            ok = ok and same
    print("[check] " + ("ALL IDENTICAL" if ok else "DIFFERENCE FOUND"))
    return ok


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("packages", nargs="*", help=f"subset of {', '.join(M.PACKAGES)}")
    ap.add_argument("--out", type=Path, required=True, help="output directory (new tables)")
    ap.add_argument("--packages-dir", type=Path, default=M.DEFAULT_PACKAGES_DIR,
                    help="the v1.0 packages (default: data/multicell)")
    ap.add_argument("--mask", default=M.DEFAULT_MASK,
                    help="mask BED, or 'hg38' for the bundled hg38 mask (default)")
    ap.add_argument("--chrom-sizes", type=Path, default=M.DEFAULT_CHROM_SIZES,
                    help="chromosome sizes for the off-end deletion "
                         "(default: data/masks/hg38.chrom.sizes)")
    ap.add_argument("--check", type=Path, default=None, metavar="V11_DIR",
                    help="compare every file written with the released v1.1 packages")
    args = ap.parse_args(argv)
    pkgs = args.packages or M.PACKAGES
    unknown = [p for p in pkgs if p not in M.PACKAGES]
    if unknown:
        ap.error(f"unknown package(s): {', '.join(unknown)}")

    mask_path = Path(M.resolve_mask(args.mask, "hg38"))
    mask = M.load_mask(str(mask_path))
    sizes = M.load_chrom_sizes(args.chrom_sizes)
    refs = {
        "mask": _display(mask_path),
        "mask_sha256": hashlib.sha256(mask_path.read_bytes()).hexdigest(),
        "chrom_sizes": _display(args.chrom_sizes),
        "blacklist_accession": "ENCFF356LFX",
    }

    all_stats, all_rec, all_rem = [], [], []
    for pkg in pkgs:
        r = process(pkg, args, mask, sizes, refs)
        all_stats.append(r["stats"])
        all_rec.append(r["records"])
        all_rem.append(r["removed"])
        s = r["stats"]
        print(f"[mask] {pkg:14s} v1.0={s['n_v10']:6d} removed={s['n_removed']:2d} "
              f"v1.1={s['n_v11']:6d} masked={s['n_masked']:5d} "
              f"(gap_only={s['n_gap_only']}, bl_only={s['n_bl_only']}, both={s['n_both']})",
              flush=True)

    impact = args.out / "impact"
    impact.mkdir(parents=True, exist_ok=True)
    st = pd.DataFrame(all_stats)
    st.to_csv(impact / "mask_package_stats.csv", index=False)
    pd.concat(all_rec, ignore_index=True).to_csv(impact / "masking_records_all.csv",
                                                 index=False)
    pd.concat(all_rem, ignore_index=True).to_csv(impact / "removed_off_chromosome_end.csv",
                                                 index=False)
    print(st.drop(columns=["tier_counts_v11"]).to_string(index=False))
    print(f"TOTAL v1.0={st.n_v10.sum()} removed={st.n_removed.sum()} "
          f"v1.1={st.n_v11.sum()} masked={st.n_masked.sum()}")
    print("->", args.out)
    if args.check is not None:
        return 0 if check(args, pkgs) else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
