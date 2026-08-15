#!/usr/bin/env python3
"""Step 3: apply the pre-registered deterministic choice rule and build the
two-column contrast.

Implements PREREGISTERED_RULE.md section 5 EXACTLY:
  reference package GM12878_10kb; primary = highest tier (D5>D4>D3>D2>D1);
  tie-break: chrom (1..22,X) -> pos -> VariationID int -> bp_end (Start<Stop);
  companion = lowest tier under the identical chain.
The record side is rendered by the SHIPPED engine functions
(tad_vci.curation.evidence_card / _tier_why), not by a local paraphrase.

Usage:  python 03_select_and_contrast.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import pandas as pd

UC = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# The shipped engine. Defaults to this repository's backend/; override with
# CREDITAD_BACKEND_DIR to point at an installed or relocated copy.
REPO_ROOT = os.path.abspath(os.path.join(UC, ".."))
BACKEND = os.environ.get("CREDITAD_BACKEND_DIR", os.path.join(REPO_ROOT, "backend"))
sys.path.insert(0, BACKEND)
from tad_vci.curation import evidence_card, _tier_why          # shipped engine
from tad_vci.graded_evidence import load_bands, combine_tier   # shipped engine

REF_PKG = "GM12878_10kb"
FALLBACK = ["GM12878_10kb", "IMR90_10kb", "HepG2_10kb",
            "GM12878_25kb", "IMR90_25kb", "HepG2_25kb"]
TIER_RANK = {"D1": 1, "D2": 2, "D3": 3, "D4": 4, "D5": 5}


def ckey(c: str) -> int:
    return 23 if str(c) == "X" else int(c)


def pick(df: pd.DataFrame, direction: str) -> pd.Series:
    """direction 'max' -> highest tier; 'min' -> lowest tier. Identical tie chain."""
    d = df.copy()
    d["_tier"] = d["ev_D_tier"].map(TIER_RANK)
    target = d["_tier"].max() if direction == "max" else d["_tier"].min()
    d = d[d["_tier"] == target]
    d["_ck"] = d["Chromosome"].map(ckey)
    d["_be"] = d["bp_end"].map({"Start": 0, "Stop": 1})
    d = d.sort_values(["_ck", "ev_pos", "VariationID_int", "_be"], kind="mergesort")
    return d.iloc[0]


def main() -> int:
    hits = pd.read_csv(os.path.join(UC, "hits_all.csv"), dtype={"Chromosome": str,
                                                               "ev_chrom": str})
    if not len(hits):
        raise SystemExit("hits_all.csv is empty — NULL_RESULT branch applies")
    hits["VariationID_int"] = hits["VariationID"].astype(np.int64)

    chosen_pkg = None
    for p in FALLBACK:
        if len(hits[hits["package"] == p]):
            chosen_pkg = p
            break
    ref = hits[hits["package"] == chosen_pkg].copy()

    tiers_present = sorted(ref["ev_D_tier"].unique(), key=lambda t: TIER_RANK[t])
    primary = pick(ref, "max")
    companion = pick(ref, "min") if len(tiers_present) > 1 else None

    sel = {
        "reference_package_declared": REF_PKG,
        "reference_package_used": chosen_pkg,
        "package_fallback_invoked": chosen_pkg != REF_PKG,
        "hits_in_reference_package": int(len(ref)),
        "tiers_present_in_reference_package": tiers_present,
        "primary_rule": "highest tier, ties: chrom -> pos -> VariationID -> bp_end",
        "companion_rule": "lowest tier, identical tie chain",
        "companion_omitted_reason": (None if companion is not None else
                                     "all hits in the reference package share one tier"),
        "n_hits_at_primary_tier": int((ref["ev_D_tier"] == primary["ev_D_tier"]).sum()),
        "n_hits_at_companion_tier": (int((ref["ev_D_tier"] == companion["ev_D_tier"]).sum())
                                     if companion is not None else None),
    }

    rows = {"primary": primary}
    if companion is not None:
        rows["companion"] = companion

    # ---------- build the two-column contrast --------------------------------
    contrast_records = []
    cards = {}
    for role, r in rows.items():
        ev = {k[3:]: r[k] for k in r.index if k.startswith("ev_")}
        card = evidence_card(ev)
        grades = {"BD1": ev["BD1_caller_support"], "BD2": ev["BD2_ctcf"],
                  "BD3": ev["BD3_rad21"], "BD4": ev["BD4_loop_anchor"]}
        why = _tier_why(grades)
        tier_recomputed, meta = combine_tier(grades)
        ctcf_b, rad21_b, calibrated, bands_used = load_bands(ev.get("data_cell_line"))
        cards[role] = {"card": card, "tier_why": why,
                       "tier_recomputed": tier_recomputed,
                       "tier_delivered": ev["D_tier"],
                       "bands_cell_used": bands_used, "bands_calibrated": bool(calibrated),
                       "ctcf_bands": dict(ctcf_b), "rad21_bands": dict(rad21_b),
                       "meta": meta, "ev": ev,
                       "locus": {k: r[k] for k in
                                 ["VariationID", "Type", "ClinicalSignificance",
                                  "ReviewStatus", "GeneSymbol", "PhenotypeList",
                                  "RCVaccession", "dbVar", "Name", "Cytogenetic",
                                  "Chromosome", "Start", "Stop", "size_bp",
                                  "bp_end", "bp_1based", "distance_bp", "package"]}}

        pos = int(ev["pos"]); R = int(ev["resolution_bp"])
        na = ev.get("not_assessable_criteria")
        na = "" if (na is None or str(na) == "nan") else str(na)

        def row(field, coord_col, rec_col, source_column):
            contrast_records.append({
                "role": role,
                "boundary_id": f"chr{ev['chrom']}:{pos}",
                "field": field,
                "coordinate_list_shows": coord_col,
                "evidence_record_shows": rec_col,
                "source_column": source_column,
            })

        row("genomic position",
            f"chr{ev['chrom']}:{pos}-{pos + R} (bin start {pos}, {R} bp bin)",
            f"chr{ev['chrom']}:{pos}-{pos + R} (identical)",
            "chrom, pos, resolution_bp")
        row("caller support (BD1)", "not shown",
            f"{int(ev['votes'])}/{int(ev['n_methods'])} callers: {ev['supporting_callers']}"
            f"; grade {ev['BD1_caller_support']}",
            "votes, n_methods, supporting_callers, BD1_caller_support")
        row("caller offsets", "not shown",
            f"{ev['supporting_offsets_bp']} bp from the candidate position",
            "supporting_offsets_bp")
        row("CTCF (BD2)", "not shown",
            (f"grade {ev['BD2_ctcf']}; fold {float(ev['ctcf_fold']):.4f} vs background "
             f"{float(ev['bg_ctcf']):.5f}; null percentile {float(ev['ctcf_null_pct']):.2f}"
             if str(ev['BD2_ctcf']) != 'not_assessable' else "grade not_assessable"),
            "BD2_ctcf, ctcf_fold, bg_ctcf, ctcf_null_pct")
        row("cohesin/RAD21 (BD3)", "not shown",
            (f"grade {ev['BD3_rad21']}; fold {float(ev['rad21_fold']):.4f} vs background "
             f"{float(ev['bg_rad21']):.5f}; null percentile {float(ev['rad21_null_pct']):.2f}"
             if str(ev['BD3_rad21']) != 'not_assessable' else "grade not_assessable"),
            "BD3_rad21, rad21_fold, bg_rad21, rad21_null_pct")
        row("loop-anchor engagement (BD4)", "not shown",
            f"grade {ev['BD4_loop_anchor']}; loop_anchor_count = "
            f"{int(ev['loop_anchor_count'])} loop ends in a {int(ev['loop_window_bp'])} bp window",
            "BD4_loop_anchor, loop_anchor_count, loop_window_bp")
        row("evidence tier", "not shown",
            f"{ev['D_tier']} (rule {ev['tier_rule_version']})",
            "D_tier, tier_rule_version")
        row("derivation of the tier", "not shown", why,
            "reproduced by tad_vci.curation._tier_why from BD1-BD4 grade columns")
        row("what was NOT measured", "not shown",
            (f"criteria not assessable: {na}" if na else
             "all four criteria were assessable at this position"),
            "not_assessable_criteria")
        row("threshold provenance", "not shown",
            f"BD2/BD3 bands from {ev['bands_cell']} ({ev['bands_source']}), "
            f"calibrated={ev['bands_calibrated']}; ChIP window {int(ev['chip_window_bp'])} bp",
            "bands_cell, bands_source, bands_calibrated, chip_window_bp")
        row("panel definition", "not shown",
            f"{ev['bd1_panel']}; rule: {ev['bd1_rule']}",
            "bd1_panel, bd1_rule")
        row("why this boundary is being looked at", "nothing in the list says",
            f"ClinVar {r['VariationID']} ({r['RCVaccession']}), {r['Type']}, "
            f"{r['ClinicalSignificance']}, chr{r['Chromosome']}:{int(r['Start'])}-{int(r['Stop'])} "
            f"({int(r['size_bp'])} bp); breakpoint {r['bp_end']} at {int(r['bp_1based'])} "
            f"lies {int(r['distance_bp'])} bp from the candidate position",
            "external: ClinVar variant_summary (see SOURCE_PROVENANCE.json)")

    cdf = pd.DataFrame(contrast_records)
    cdf.to_csv(os.path.join(UC, "contrast_table.csv"), index=False)

    with open(os.path.join(UC, "selection_record.json"), "w") as fh:
        json.dump({"selection": sel,
                   "cards": {k: {kk: vv for kk, vv in v.items() if kk != "ev"}
                             for k, v in cards.items()}},
                  fh, indent=2, default=str)

    print(json.dumps({"selection": sel,
                      "primary": {"boundary": f"chr{cards['primary']['ev']['chrom']}:"
                                              f"{int(cards['primary']['ev']['pos'])}",
                                  "tier": cards['primary']['ev']['D_tier'],
                                  "VariationID": cards['primary']['locus']['VariationID'],
                                  "why": cards['primary']['tier_why']},
                      "companion": ({"boundary": f"chr{cards['companion']['ev']['chrom']}:"
                                                 f"{int(cards['companion']['ev']['pos'])}",
                                     "tier": cards['companion']['ev']['D_tier'],
                                     "VariationID": cards['companion']['locus']['VariationID'],
                                     "why": cards['companion']['tier_why']}
                                    if 'companion' in cards else None)},
                     indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
