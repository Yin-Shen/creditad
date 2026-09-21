"""Referee 2 / Major 5 fix (A1): build the v1.1 multicell packages.

WHAT CHANGES relative to v1.0
  (a) Every candidate whose `pos` falls inside an hg38 assembly gap or inside the ENCODE
      hg38 blacklist (ENCFF356LFX) has BD2_ctcf and BD3_rad21 re-graded to
      `not_assessable`. At such a position the ChIP readout is not a measured absence, it
      is an absence of measurement, and the paper's whole distinction is that one.
      BD1 and BD4 are NOT masked:
        BD1 is caller agreement over the Hi-C contact map, a property of the callers'
        output at that coordinate, not of the ChIP mappability at it;
        BD4 counts loop-call anchors within +/-50 kb, and a loop call whose anchor lands
        in a gap is still a call that exists in the supplied loop file.
      Only mappability-limited *ChIP* axes are withheld.
  (b) The 12 candidates whose `pos` is at or beyond the end of their chromosome are
      DELETED (BED-generation artefact); they are listed in
      out/removed_off_chromosome_end.csv.
  (c) D_tier and not_assessable_criteria are recomputed for affected rows ONLY, with the
      SHIPPED backend/tad_vci/graded_evidence.combine_tier. Every other row is passed
      through byte-for-byte.
  (d) union_boundaries.bed / tier_boundaries.bed are regenerated (they are derived from
      annotation.tsv, one line per candidate, and step2 proves the derivation is exact).
      The six per-caller BEDs are copied verbatim -- they are upstream caller output and
      the masking pass has no business editing them.

RETAINED ON PURPOSE: the raw ctcf / rad21 / *_fold / *_null_pct columns keep their v1.0
values on masked rows. The 38-column schema is fixed, so dropping them was not an option,
and a reader is better served by seeing the number that was withheld than by a blank. The
consequence is stated in the report: re-running graded_evidence.annotate() on a v1.1 table
regrades those axes back to `none`; masking is a documented post-engine layer, applied by
this script, not something the engine knows about.

Usage:  python mask_and_retier.py [PKG ...]      (default: all six)
"""
from __future__ import annotations

import io
import json
import shutil
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mask_common as M

CALLER_BEDS = ["TopDom", "SpectralTAD", "OnTAD", "MSTD", "arrowhead", "DI"]
MASK_CRITERIA = ["BD2", "BD3"]
MASK_COLS = {"BD2": "BD2_ctcf", "BD3": "BD3_rad21"}

MASK_NOTE = (
    "v1.1 mappability masking (referee 2, Major 5). A candidate whose bin start falls "
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


def build_manifest(pkg: str, man: dict, stats: dict, tier_counts: dict, n_rows: int) -> dict:
    out: dict = {}
    for k, v in man.items():
        out[k] = v
        if k == "name":
            out["data_version"] = "1.1"
    out["n_candidates"] = n_rows
    out["tier_counts"] = tier_counts
    out["masking"] = {
        "applied": True,
        "referee_item": "referee 2, Major 5",
        "description": MASK_NOTE,
        "predicate": "point-in-interval on the candidate `pos`, half-open [start, end); "
                     "identical to review/referee2/r2_04_unmappable_and_windows.py",
        "reference_files": {
            "assembly_gap": str(M.GAP_F),
            "encode_blacklist": str(M.BLACKLIST_F),
            "chrom_sizes": str(M.CHROM_SIZES_F),
            "blacklist_accession": "ENCFF356LFX",
        },
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
            "Masking is applied AFTER graded_evidence.annotate(); the engine itself has no "
            "mappability model. Re-running annotate() on this table would regrade the "
            "masked axes back to `none`. Reproduce v1.1 with "
            "rebuild_2026-08-25_fixes/masking/mask_and_retier.py.",
        "records_file": "masking_records.tsv",
        "generator": "rebuild_2026-08-25_fixes/masking/mask_and_retier.py",
        "seed": M.SEED,
    }
    out["supersedes"] = {
        "data_version": "1.0",
        "package_dir": f"example_data/multicell/{pkg}",
        "n_candidates": man["n_candidates"],
        "tier_counts": man["tier_counts"],
    }
    return out


def process(pkg: str, outroot: Path) -> dict:
    sizes = M.load_chrom_sizes()
    gap, bl = M.load_gap(), M.load_blacklist()

    df = M.load_package(pkg)
    cols = list(df.columns)
    res = int(df["resolution_bp"].iloc[0])
    ipos = df["pos"].astype("int64")
    flags = M.flag_rows(df.assign(pos=ipos), sizes, gap, bl)

    # (b) delete off-chromosome-end rows
    removed = df[flags.off_end].copy()
    removed.insert(0, "package", pkg)
    removed_out = removed[["package", "chrom", "pos", "votes", "BD1_caller_support",
                           "BD2_ctcf", "BD3_rad21", "BD4_loop_anchor", "D_tier"]].copy()
    removed_out["chrom_size_hg38"] = [sizes.get(str(c), -1) for c in removed["chrom"]]
    removed_out["reason"] = "pos_at_or_beyond_chromosome_end"

    keep = ~flags.off_end
    df = df[keep].reset_index(drop=True)
    fl = flags[keep].reset_index(drop=True)
    ipos = ipos[keep].reset_index(drop=True)

    # (a) mask BD2/BD3
    mask = (fl.in_gap | fl.in_blacklist).to_numpy()
    n_masked = int(mask.sum())
    before = df.loc[mask, ["BD2_ctcf", "BD3_rad21", "BD4_loop_anchor", "BD1_caller_support",
                           "D_tier", "not_assessable_criteria"]].copy()
    for c in MASK_COLS.values():
        df.loc[mask, c] = "not_assessable"

    # (c) recompute tier for affected rows only, with the shipped engine
    sub = df.loc[mask]
    new_tier, new_na = M.tier_vec(sub["BD1_caller_support"], sub["BD2_ctcf"],
                                  sub["BD3_rad21"], sub["BD4_loop_anchor"])
    df.loc[mask, "D_tier"] = new_tier
    df.loc[mask, "not_assessable_criteria"] = new_na
    assert list(df.columns) == cols, "38-column schema changed"

    # ---- outputs
    d = outroot / pkg
    d.mkdir(parents=True, exist_ok=True)
    (d / "annotation.tsv").write_bytes(
        df.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    (d / "union_boundaries.bed").write_bytes(union_bytes(df["chrom"], ipos, res))
    (d / "tier_boundaries.bed").write_bytes(tier_bytes(df["chrom"], ipos, df["D_tier"], res))
    for c in CALLER_BEDS:
        shutil.copy2(M.V10 / pkg / f"{c}_boundaries.bed", d / f"{c}_boundaries.bed")

    rec = pd.DataFrame({
        "chrom": sub["chrom"].to_numpy(),
        "pos": ipos[mask].to_numpy(),
        "reason": [M.mask_reason(g, b) for g, b in zip(fl.in_gap[mask], fl.in_blacklist[mask])],
        "BD1_caller_support": before["BD1_caller_support"].to_numpy(),
        "BD2_ctcf_v10": before["BD2_ctcf"].to_numpy(),
        "BD3_rad21_v10": before["BD3_rad21"].to_numpy(),
        "BD4_loop_anchor": before["BD4_loop_anchor"].to_numpy(),
        "D_tier_v10": before["D_tier"].to_numpy(),
        "D_tier_v11": new_tier,
        "not_assessable_criteria_v11": new_na,
    })
    rec.to_csv(d / "masking_records.tsv", sep="\t", index=False)

    tier_counts = dict(sorted(Counter(df["D_tier"]).items(), key=lambda kv: -kv[1]))
    stats = dict(
        package=pkg,
        n_v10=int(len(df) + len(removed)),
        n_removed=int(len(removed)),
        n_v11=int(len(df)),
        n_masked=n_masked,
        n_gap_only=int((fl.in_gap & ~fl.in_blacklist).sum()),
        n_bl_only=int((~fl.in_gap & fl.in_blacklist).sum()),
        n_both=int((fl.in_gap & fl.in_blacklist).sum()),
    )
    man = build_manifest(pkg, M.load_manifest(pkg), stats, tier_counts, len(df))
    with (d / "manifest.json").open("w") as fh:
        json.dump(man, fh, indent=2)
        fh.write("\n")

    stats["tier_counts_v11"] = json.dumps(tier_counts)
    return {"stats": stats, "records": rec.assign(package=pkg),
            "removed": removed_out, "before": before, "mask": mask, "df": df}


def main(argv):
    pkgs = argv[1:] or M.PACKAGES
    outroot = M.V11
    outroot.mkdir(parents=True, exist_ok=True)
    M.OUT.mkdir(parents=True, exist_ok=True)

    all_stats, all_rec, all_rem = [], [], []
    for pkg in pkgs:
        r = process(pkg, outroot)
        all_stats.append(r["stats"])
        all_rec.append(r["records"])
        all_rem.append(r["removed"])
        s = r["stats"]
        print(f"[mask] {pkg:14s} v1.0={s['n_v10']:6d} removed={s['n_removed']:2d} "
              f"v1.1={s['n_v11']:6d} masked={s['n_masked']:5d} "
              f"(gap_only={s['n_gap_only']}, bl_only={s['n_bl_only']}, both={s['n_both']})",
              flush=True)

    st = pd.DataFrame(all_stats)
    st.to_csv(M.OUT / "mask_package_stats.csv", index=False)
    rec = pd.concat(all_rec, ignore_index=True)
    rec.to_csv(M.OUT / "masking_records_all.csv", index=False)
    rem = pd.concat(all_rem, ignore_index=True)
    rem.to_csv(M.OUT / "removed_off_chromosome_end.csv", index=False)

    print(st.drop(columns=["tier_counts_v11"]).to_string(index=False))
    print(f"TOTAL v1.0={st.n_v10.sum()} removed={st.n_removed.sum()} "
          f"v1.1={st.n_v11.sum()} masked={st.n_masked.sum()}")
    print("->", outroot)


if __name__ == "__main__":
    main(sys.argv)
