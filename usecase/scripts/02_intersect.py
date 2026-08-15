#!/usr/bin/env python3
"""Step 2 of the CrediTAD use-case demonstration: intersect the pre-registered
locus set with all six delivered evidence packages.

Implements PREREGISTERED_RULE.md section 3 EXACTLY:
    breakpoint b (1-based) hits candidate (chrom, pos) at resolution R iff
        | (b - 1) - pos | <= R
Both breakpoints (Start, Stop) of every locus are tested independently.
Every hit is written. No curation, no post-hoc exclusion.

Usage:  python 02_intersect.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import pandas as pd

UC = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

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
PACKAGES = ["GM12878_10kb", "GM12878_25kb", "HepG2_10kb", "HepG2_25kb",
            "IMR90_10kb", "IMR90_25kb"]
EXPECTED_TOTAL_ROWS = 284_744          # GROUP_BRIEF check value
TIERS = ["D1", "D2", "D3", "D4", "D5"]


def ckey(c: str) -> int:
    return 23 if c == "X" else int(c)


def main() -> int:
    loci = pd.read_csv(os.path.join(UC, "data", "locus_set.tsv"), sep="\t", dtype=str)
    loci["Start"] = loci["Start"].astype(np.int64)
    loci["Stop"] = loci["Stop"].astype(np.int64)
    loci["VariationID_int"] = loci["VariationID"].astype(np.int64)
    n_loci = len(loci)

    # long form: one row per breakpoint (2 per locus), both tested independently
    bps = []
    for end_name in ("Start", "Stop"):
        t = loci.copy()
        t["bp_end"] = end_name
        t["bp_1based"] = t[end_name]
        bps.append(t)
    bp = pd.concat(bps, ignore_index=True)
    bp["bp_0based"] = bp["bp_1based"] - 1

    loaded_rows = 0
    hits_frames = []
    per_package = {}

    for pkg in PACKAGES:
        ann_path = annotation_path(pkg)
        ann = pd.read_csv(ann_path, sep="\t", dtype={"chrom": str}, low_memory=False)
        loaded_rows += len(ann)
        assert len(ann.columns) == 38, f"{pkg}: expected 38 columns, got {len(ann.columns)}"

        res_vals = ann["resolution_bp"].unique()
        assert len(res_vals) == 1, f"{pkg}: non-unique resolution_bp {res_vals}"
        R = int(res_vals[0])
        assert R == (10_000 if pkg.endswith("10kb") else 25_000), f"{pkg}: R={R}"

        ann["chrom"] = ann["chrom"].astype(str).str.replace(r"^chr", "", regex=True)

        pkg_hits = []
        for chrom, g in ann.groupby("chrom", sort=False):
            sub = bp[bp["Chromosome"] == chrom]
            if sub.empty:
                continue
            g = g.sort_values("pos", kind="mergesort").reset_index(drop=True)
            pos = g["pos"].to_numpy(dtype=np.int64)
            b0 = sub["bp_0based"].to_numpy(dtype=np.int64)
            lo = np.searchsorted(pos, b0 - R, side="left")
            hi = np.searchsorted(pos, b0 + R, side="right")
            counts = hi - lo
            if counts.sum() == 0:
                continue
            bp_idx = np.repeat(np.arange(len(sub)), counts)
            cand_idx = np.concatenate([np.arange(l, h) for l, h in zip(lo, hi) if h > l])
            left = sub.iloc[bp_idx].reset_index(drop=True)
            right = g.iloc[cand_idx].reset_index(drop=True)
            right.columns = [f"ev_{c}" for c in right.columns]
            m = pd.concat([left.reset_index(drop=True), right], axis=1)
            m["package"] = pkg
            m["resolution_bp_pkg"] = R
            m["distance_bp"] = (m["bp_0based"] - m["ev_pos"]).abs()
            pkg_hits.append(m)

        if pkg_hits:
            ph = pd.concat(pkg_hits, ignore_index=True)
        else:
            ph = pd.DataFrame()
        per_package[pkg] = ph
        hits_frames.append(ph)

    assert loaded_rows == EXPECTED_TOTAL_ROWS, \
        f"corpus row total {loaded_rows} != expected {EXPECTED_TOTAL_ROWS}"

    hits = pd.concat([h for h in hits_frames if len(h)], ignore_index=True) \
        if any(len(h) for h in hits_frames) else pd.DataFrame()

    # ---- deterministic output order (PREREG §4) ---------------------------
    if len(hits):
        hits["_ck"] = hits["Chromosome"].map(ckey)
        hits["_pk"] = hits["package"].map({p: i for i, p in enumerate(PACKAGES)})
        hits["_be"] = hits["bp_end"].map({"Start": 0, "Stop": 1})
        hits = hits.sort_values(["_pk", "_ck", "ev_pos", "VariationID_int", "_be"],
                                kind="mergesort").reset_index(drop=True)
        hits = hits.drop(columns=["_ck", "_pk", "_be"])
        front = ["package", "resolution_bp_pkg", "VariationID", "Type",
                 "ClinicalSignificance", "ReviewStatus", "GeneSymbol", "PhenotypeList",
                 "RCVaccession", "dbVar", "Name", "Cytogenetic", "Chromosome",
                 "Start", "Stop", "size_bp", "bp_end", "bp_1based", "bp_0based",
                 "distance_bp"]
        ev = [c for c in hits.columns if c.startswith("ev_")]
        rest = [c for c in hits.columns if c not in front + ev]
        hits = hits[front + ev + rest]
        hits.to_csv(os.path.join(UC, "hits_all.csv"), index=False)
    else:
        pd.DataFrame().to_csv(os.path.join(UC, "hits_all.csv"), index=False)

    funnel = json.load(open(os.path.join(UC, "data", "locus_funnel.json")))
    summary = {
        "prereg_file": "usecase/PREREGISTERED_RULE.md",
        "prereg_sha256": json.load(open(os.path.join(UC, "PREREG_SEAL.json")))["sha256"],
        "data_source": json.load(open(os.path.join(UC, "SOURCE_PROVENANCE.json")))["source_name"],
        "corpus_rows_loaded": int(loaded_rows),
        "corpus_rows_expected": EXPECTED_TOTAL_ROWS,
        "corpus_row_check_pass": bool(loaded_rows == EXPECTED_TOTAL_ROWS),
        "locus_funnel": funnel["funnel"],
        "duplicate_tuples_collapsed": funnel["duplicate_tuples_collapsed"],
        "loci_tested": int(n_loci),
        "breakpoints_tested": int(len(bp)),
        "total_hits": int(len(hits)),
        "hits_per_package": {p: int(len(per_package[p])) for p in PACKAGES},
        "distinct_loci_with_hit": int(hits["VariationID"].nunique()) if len(hits) else 0,
        "distinct_candidates_with_hit": int(
            hits.groupby(["package", "Chromosome", "ev_pos"]).ngroups) if len(hits) else 0,
        "loci_with_hit_per_package": {
            p: int(per_package[p]["VariationID"].nunique()) if len(per_package[p]) else 0
            for p in PACKAGES},
        "breakpoint_end_split": (hits["bp_end"].value_counts().to_dict() if len(hits) else {}),
        "tier_distribution_of_hits_overall": (
            hits["ev_D_tier"].value_counts().reindex(TIERS).fillna(0).astype(int).to_dict()
            if len(hits) else {}),
        "tier_distribution_of_hits_per_package": {
            p: (per_package[p]["ev_D_tier"].value_counts()
                .reindex(TIERS).fillna(0).astype(int).to_dict()
                if len(per_package[p]) else {}) for p in PACKAGES},
        "distance_bp_distribution": (
            hits["distance_bp"].value_counts().sort_index().to_dict() if len(hits) else {}),
        "script": "usecase/scripts/02_intersect.py",
        "output": "usecase/hits_all.csv",
    }
    with open(os.path.join(UC, "summary_counts.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=int)

    print(json.dumps({k: summary[k] for k in
                      ["corpus_rows_loaded", "corpus_row_check_pass", "loci_tested",
                       "breakpoints_tested", "total_hits", "hits_per_package",
                       "distinct_loci_with_hit", "distinct_candidates_with_hit",
                       "tier_distribution_of_hits_overall", "breakpoint_end_split"]},
                     indent=2, default=int))
    return 0


if __name__ == "__main__":
    sys.exit(main())
