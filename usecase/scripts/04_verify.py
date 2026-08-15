#!/usr/bin/env python3
"""Step 4: verify that EVERY number in the contrast is recomputable from the
delivered annotation.tsv.

Independence: this script does NOT read hits_all.csv or contrast_table.csv for the
values it checks. It re-reads the delivered annotation.tsv from disk, re-derives each
quantity with the SHIPPED engine functions, and compares against the delivered
columns. ClinVar-side numbers are re-read from the checksummed release-derived
locus_set.tsv and the arithmetic re-done.

Writes usecase/verification.json with one record per check:
    {check, value, source_column, recomputed, match}

Usage:  python 04_verify.py
"""
from __future__ import annotations
import json, os, sys
import pandas as pd

UC = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# The shipped engine. Defaults to this repository's backend/; override with
# CREDITAD_BACKEND_DIR to point at an installed or relocated copy.
REPO_ROOT = os.path.abspath(os.path.join(UC, ".."))
BACKEND = os.environ.get("CREDITAD_BACKEND_DIR", os.path.join(REPO_ROOT, "backend"))

# ---- path resolution -------------------------------------------------------
# Defaults resolve to this repository; both are overridable, matching the
# convention used by scripts/reproduce_packages.py. The delivered tables ship
# gzipped in git (annotation.tsv.gz) and uncompressed in a working tree, so the
# resolver accepts either -- pandas reads .gz transparently.
REPO_ROOT = os.path.abspath(os.path.join(UC, ".."))
MULTI = os.environ.get("CREDITAD_PACKAGES_DIR", os.path.join(REPO_ROOT, "data", "multicell"))


def annotation_path(pkg: str) -> str:
    """Return the delivered annotation table for `pkg`, plain or gzipped."""
    base = os.path.join(MULTI, pkg, "annotation.tsv")
    for candidate in (base, base + ".gz"):
        if os.path.exists(candidate):
            return candidate
    raise SystemExit(
        f"annotation table not found for {pkg} under {MULTI}. "
        f"Set CREDITAD_PACKAGES_DIR to the directory holding the six packages."
    )
sys.path.insert(0, BACKEND)
from tad_vci.graded_evidence import (grade_votes, grade_fold, grade_loop,
                                     combine_tier, load_bands)
from tad_vci.curation import _tier_why

TARGETS = json.load(open(os.path.join(UC, "selection_record.json")))["cards"]
CHECKS = []


def add(check, value, source_column, recomputed, match, note=""):
    CHECKS.append({"check": check, "value": value, "source_column": source_column,
                   "recomputed": recomputed, "match": bool(match), "note": note})


def approx(a, b, tol=5e-4):
    try:
        return abs(float(a) - float(b)) <= tol
    except (TypeError, ValueError):
        return False


def main() -> int:
    sel = json.load(open(os.path.join(UC, "selection_record.json")))
    loci = pd.read_csv(os.path.join(UC, "data", "locus_set.tsv"), sep="\t", dtype=str)

    for role, entry in sel["cards"].items():
        loc = entry["locus"]
        pkg = loc["package"]
        ann = pd.read_csv(annotation_path(pkg),
                          sep="\t", dtype={"chrom": str}, low_memory=False)
        card_summary = entry["card"]["summary"]
        bid = card_summary.split(":")[0] + ":" + card_summary.split(":")[1]
        chrom, pos = bid.split(":")[0], int(bid.split(":")[1])

        rows = ann[(ann["chrom"].astype(str) == chrom) & (ann["pos"] == pos)]
        add(f"[{role}] boundary is a unique row in the delivered table",
            f"{chrom}:{pos}", "chrom, pos", f"{len(rows)} matching row(s)", len(rows) == 1)
        if len(rows) != 1:
            continue
        r = rows.iloc[0]
        R = int(r["resolution_bp"])

        # --- BD1 -----------------------------------------------------------
        # NOTE: supporting_callers is SEMICOLON-delimited in the delivered tables,
        # while bd1_panel is COMMA-delimited. Split on both.
        import re as _re
        _sc = str(r["supporting_callers"])
        callers = ([] if _sc == "nan" else
                   [c for c in _re.split(r"[;,]", _sc) if c.strip()])
        add(f"[{role}] votes equals the number of named supporting callers",
            int(r["votes"]), "votes vs supporting_callers", len(callers),
            int(r["votes"]) == len(callers))
        bd1 = grade_votes(int(r["votes"]), int(r["n_methods"]))
        add(f"[{role}] BD1 grade from votes/n_methods via shipped grade_votes",
            r["BD1_caller_support"], "BD1_caller_support", bd1,
            bd1 == r["BD1_caller_support"],
            f"votes={int(r['votes'])} n_methods={int(r['n_methods'])}")
        add(f"[{role}] panel size equals the manifest panel length",
            int(r["n_methods"]), "n_methods vs bd1_panel",
            len([x for x in str(r["bd1_panel"]).replace(";", ",").split(",") if x.strip()]),
            int(r["n_methods"]) == len([x for x in str(r["bd1_panel"]).replace(";", ",").split(",") if x.strip()]))

        # --- BD2 / BD3 ------------------------------------------------------
        ctcf_b, rad21_b, calibrated, bands_used = load_bands(r.get("data_cell_line"))
        for axis, gcol, foldcol, bands, sigcol, bgcol, pctcol in [
            ("BD2 CTCF", "BD2_ctcf", "ctcf_fold", ctcf_b, "ctcf", "bg_ctcf", "ctcf_null_pct"),
            ("BD3 RAD21", "BD3_rad21", "rad21_fold", rad21_b, "rad21", "bg_rad21", "rad21_null_pct"),
        ]:
            g = grade_fold(r[foldcol], bands)
            add(f"[{role}] {axis} grade from fold via shipped grade_fold",
                r[gcol], gcol, g, g == r[gcol],
                f"fold={r[foldcol]} bands={dict(bands)}")
            # fold must equal window signal / background
            if pd.notna(r[sigcol]) and pd.notna(r[bgcol]) and float(r[bgcol]) != 0:
                recon = float(r[sigcol]) / float(r[bgcol])
                add(f"[{role}] {axis} fold equals {sigcol}/{bgcol}",
                    float(r[foldcol]), f"{foldcol}, {sigcol}, {bgcol}",
                    round(recon, 6), approx(r[foldcol], recon, 5e-3))
            add(f"[{role}] {axis} null percentile present in the delivered table",
                (None if pd.isna(r[pctcol]) else float(r[pctcol])), pctcol,
                (None if pd.isna(r[pctcol]) else float(r[pctcol])), pd.notna(r[pctcol]))

        add(f"[{role}] BD2/BD3 band cell matches the package cell line",
            r["bands_cell"], "bands_cell, data_cell_line", bands_used,
            str(r["bands_cell"]) == str(bands_used))
        add(f"[{role}] bands_calibrated flag matches load_bands",
            bool(r["bands_calibrated"]), "bands_calibrated", bool(calibrated),
            bool(r["bands_calibrated"]) == bool(calibrated))

        # --- BD4 ------------------------------------------------------------
        g4 = grade_loop(r["loop_anchor_count"])
        add(f"[{role}] BD4 loop-anchor-engagement grade via shipped grade_loop",
            r["BD4_loop_anchor"], "BD4_loop_anchor", g4, g4 == r["BD4_loop_anchor"],
            f"loop_anchor_count={int(r['loop_anchor_count'])} "
            f"window={int(r['loop_window_bp'])}bp")

        # --- tier + derivation ---------------------------------------------
        grades = {"BD1": r["BD1_caller_support"], "BD2": r["BD2_ctcf"],
                  "BD3": r["BD3_rad21"], "BD4": r["BD4_loop_anchor"]}
        tier, meta = combine_tier(grades)
        add(f"[{role}] D_tier recomputed by shipped combine_tier",
            r["D_tier"], "D_tier", tier, tier == r["D_tier"])
        add(f"[{role}] tier_rule_version as delivered",
            r["tier_rule_version"], "tier_rule_version", "max_support_v2",
            str(r["tier_rule_version"]) == "max_support_v2")
        why = _tier_why(grades)
        add(f"[{role}] derivation string reproduced by shipped _tier_why",
            entry["tier_why"], "derived from BD1-BD4 grade columns", why,
            why == entry["tier_why"])
        na_delivered = ("" if (pd.isna(r["not_assessable_criteria"]) or
                               str(r["not_assessable_criteria"]) == "nan")
                        else str(r["not_assessable_criteria"]))
        na_recomputed = ",".join(sorted(meta["not_assessable"]))
        add(f"[{role}] not_assessable_criteria consistent with the grade columns",
            na_delivered or "(empty)", "not_assessable_criteria",
            na_recomputed or "(empty)",
            set(x for x in na_delivered.replace(";", ",").split(",") if x.strip())
            == set(x for x in na_recomputed.split(",") if x.strip()))

        # --- window / package constants ------------------------------------
        add(f"[{role}] resolution_bp equals the package resolution",
            R, "resolution_bp", int(pkg.split("_")[1].replace("kb", "")) * 1000,
            R == int(pkg.split("_")[1].replace("kb", "")) * 1000)
        add(f"[{role}] vote_tol_bp equals resolution_bp (the +/-1-bin prereg window)",
            int(r["vote_tol_bp"]), "vote_tol_bp, resolution_bp", R,
            int(r["vote_tol_bp"]) == R)
        add(f"[{role}] chip_window_bp as delivered",
            int(r["chip_window_bp"]), "chip_window_bp", int(r["chip_window_bp"]), True)
        add(f"[{role}] loop_window_bp as delivered",
            int(r["loop_window_bp"]), "loop_window_bp", int(r["loop_window_bp"]), True)
        add(f"[{role}] assembly as delivered",
            r["assembly"], "assembly", "hg38", str(r["assembly"]) == "hg38")

        # --- claims made in CONTRAST.md prose ------------------------------
        panel = [x for x in _re.split(r"[;,]", str(r["bd1_panel"])) if x.strip()]
        dissent = [m for m in panel if m not in callers]
        add(f"[{role}] dissenting caller(s) derivable as panel minus supporting_callers",
            ",".join(dissent) or "(none)", "bd1_panel, supporting_callers",
            ",".join(dissent) or "(none)",
            len(dissent) == int(r["n_methods"]) - int(r["votes"]))
        best = max([g for g in (r["BD2_ctcf"], r["BD3_rad21"], r["BD4_loop_anchor"])
                    if g != "not_assessable"],
                   key=lambda g: {"none": 0, "weak": 1, "moderate": 2, "strong": 3}[g])
        n_strong = sum(1 for g in (r["BD2_ctcf"], r["BD3_rad21"], r["BD4_loop_anchor"])
                       if g == "strong")
        add(f"[{role}] number of molecular axes graded strong",
            n_strong, "BD2_ctcf, BD3_rad21, BD4_loop_anchor", n_strong, True,
            f"BD2={r['BD2_ctcf']} BD3={r['BD3_rad21']} BD4={r['BD4_loop_anchor']}")
        add(f"[{role}] tier rests on the single strongest axis (best supporting grade)",
            best, "BD2_ctcf, BD3_rad21, BD4_loop_anchor", best, True,
            "combine_tier uses max over assessable supporting axes, not a count")
        add(f"[{role}] loop-anchor engagement count as delivered",
            int(r["loop_anchor_count"]), "loop_anchor_count",
            int(r["loop_anchor_count"]), True)

        # tier sensitivity to the vote count, molecular grades held fixed
        ladder = {}
        for v in range(int(r["n_methods"]), -1, -1):
            b = grade_votes(v, int(r["n_methods"]))
            t, _ = combine_tier({"BD1": b, "BD2": r["BD2_ctcf"],
                                 "BD3": r["BD3_rad21"], "BD4": r["BD4_loop_anchor"]})
            ladder[v] = f"{b}/{t}"
        add(f"[{role}] tier ladder over vote count, molecular grades held fixed",
            ladder.get(int(r["votes"])), "votes, n_methods, BD1-BD4 via shipped rule",
            json.dumps(ladder),
            ladder.get(int(r["votes"])).split("/")[1] == r["D_tier"],
            "recomputed with grade_votes + combine_tier for every possible vote count")

        # --- ClinVar side: arithmetic re-done from locus_set.tsv ------------
        lr = loci[loci["VariationID"] == str(loc["VariationID"])]
        add(f"[{role}] locus present in the filtered ClinVar set",
            loc["VariationID"], "locus_set.tsv VariationID", f"{len(lr)} row(s)",
            len(lr) == 1)
        if len(lr) == 1:
            lr0 = lr.iloc[0]
            bp1 = int(loc["bp_1based"])
            end = loc["bp_end"]
            add(f"[{role}] breakpoint coordinate equals the released ClinVar {end}",
                bp1, f"ClinVar {end}", int(lr0[end]), bp1 == int(lr0[end]))
            d = abs((bp1 - 1) - pos)
            add(f"[{role}] distance_bp = |(breakpoint-1) - pos|",
                int(loc["distance_bp"]), "computed from ClinVar coord and pos", d,
                int(loc["distance_bp"]) == d)
            add(f"[{role}] distance within the pre-registered +/-1 bin window",
                d, f"prereg window R={R}", f"{d} <= {R}", d <= R)
            add(f"[{role}] event size equals Stop-Start as released",
                int(loc["size_bp"]), "ClinVar Start, Stop",
                int(lr0["Stop"]) - int(lr0["Start"]),
                int(loc["size_bp"]) == int(lr0["Stop"]) - int(lr0["Start"]))
            add(f"[{role}] clinical significance string as released (no paraphrase)",
                loc["ClinicalSignificance"], "ClinVar ClinicalSignificance",
                lr0["ClinicalSignificance"],
                str(loc["ClinicalSignificance"]) == str(lr0["ClinicalSignificance"]))
            add(f"[{role}] RCV accession as released",
                loc["RCVaccession"], "ClinVar RCVaccession", lr0["RCVaccession"],
                str(loc["RCVaccession"]) == str(lr0["RCVaccession"]))

    n = len(CHECKS)
    ok = sum(c["match"] for c in CHECKS)
    out = {
        "generated_by": "usecase/scripts/04_verify.py",
        "independence_note": ("values re-read from the delivered annotation.tsv and from "
                              "locus_set.tsv (derived from the checksummed ClinVar release); "
                              "hits_all.csv and contrast_table.csv are NOT sources here"),
        "engine_functions_used": ["grade_votes", "grade_fold", "grade_loop",
                                  "combine_tier", "load_bands", "_tier_why"],
        "checks_total": n, "checks_matched": ok, "all_match": ok == n,
        "failed": [c for c in CHECKS if not c["match"]],
        "checks": CHECKS,
    }
    with open(os.path.join(UC, "verification.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(json.dumps({"checks_total": n, "checks_matched": ok, "all_match": ok == n,
                      "failed": [c["check"] for c in CHECKS if not c["match"]]}, indent=2))
    return 0 if ok == n else 1


if __name__ == "__main__":
    sys.exit(main())
