"""Structured evidence-card engine for already-called TAD boundaries.

The engine records four named criteria and applies the inspectable ``max_support_v2``
rule to produce a five-level *evidence tier* (D1--D5). The tier is not a calibrated
probability, a truth label or an expert decision.

BD1 is caller support from the supplied multi-caller panel and is the required primary field.
BD2 (CTCF), BD3 (RAD21) and BD4 (loop anchors) are optional supporting observations.
Measured negative evidence (``none``) and an unavailable assay (``not_assessable``)
are distinct on the card AND, since ``max_support_v2``, in the tier: a boundary whose
every measured supporting axis came back depleted cannot reach D4/D5 on caller
agreement alone, while an axis nobody measured stays neutral (see the depletion clause
in ``combine_tier``). Only the strongest assessable supporting grade enters the tier,
so correlated supporting assays are not summed.

BD5, the experimental BCP/conformal model and all external validation analyses are
outside this engine. Because BD2--BD4 are tier inputs, enrichment of those same
signals cannot validate the tier. Scientific performance claims require separate,
manifest-bound analyses; this module enforces only the data and rule contract.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import pandas as pd

Grade = Literal["strong", "moderate", "weak", "none", "not_assessable"]
Tier = Literal["D5", "D4", "D3", "D2", "D1"]

# Each criterion is classified on TWO independent, non-overlapping axes (a reviewer
# must be able to read both unambiguously; see BD_SPEC.md for the canonical table):
#   AVAILABILITY — how often the required data exists for an arbitrary cell line:
#     intrinsic        always (computed from the Hi-C contact map itself)
#     standard_assay   usually (a common public ChIP-seq / loop call)
#     reference_panel  needs a panel of OTHER cell lines
#   ROLE — how the criterion enters the evidence-tier rule:
#     primary          the obligatory anchor of the tier (must be assessable)
#     supporting       refines the tier when present; an UNMEASURED axis never
#                      downgrades, but all-measured-and-depleted does (v2 clause)
# The former held-out-validator role and BD5 experiment are outside this contract.
# External analyses live in manifest-gated analysis modules, never in the tier engine.
Availability = Literal["intrinsic", "standard_assay", "reference_panel", "genome_sequence"]
Role = Literal["primary", "supporting"]

# Default bands = the canonical GM12878 calibration, DERIVED reproducibly by
# analyses/tad_vci_20260529/calibrate_bands.py (random +/-25 kb windows, fold =
# window/background-mean, averaged over 5 seeds). These constants MUST equal
# calibrated_bands.json["GM12878"] (a unit test asserts it, so engine and calibration
# cannot drift).
#
# BAND LADDER v2 (2026-07-29) — every rung is a quantile of the cell's OWN random-window
# fold null: weak = p75, moderate = p90, strong = p97.5; anything below the weak cut
# grades "none". v1 pinned the weak floor to fold = 1.0 ("at least background"), which
# is ~the MEDIAN of that null: 50.01% of random IMR90 windows cleared the RAD21 weak
# band and 44.78% cleared CTCF (measured; stored per track as
# p_random_at_or_above_background). A rung half of random positions reach is not
# evidence, and it was load-bearing — it was the only reason some packages had any D2 at
# all. Raising only the weak floor to p75 was impossible without colliding with v1's
# moderate cut (also p75), so the whole ladder moved up and strong gained a new p97.5
# cut. Migration measured per package in docs/MEASURED_LIMITATIONS.md.
#
# "none" therefore now means "not in the top quartile of this cell's random windows",
# NOT "depleted below background". Read the depletion clause in combine_tier with that
# meaning; the served card carries the raw fold, the background and the cutpoints so the
# distinction is visible per boundary.
#
# Bands are PER-CELL: other cell lines pass their own bands to annotate() (the
# bg-invariant, correct recipe — a boundary is graded against the percentiles of its OWN
# cell's random windows), since fold distributions differ across cells.
# WINDOW STATISTIC: the fold numerator is the MEAN over the +/-25 kb window, and the bands
# below are percentiles of that same mean over random windows — calibrate and serve agree by
# construction. `max` or `median` would be defensible alternatives but would need their own
# recalibration; no sensitivity analysis over the choice ships with this version.
CTCF_BANDS = (("strong", 2.141), ("moderate", 1.564), ("weak", 1.198))
RAD21_BANDS = (("strong", 1.641), ("moderate", 1.394), ("weak", 1.212))

import json as _json
from pathlib import Path as _Path

_BANDS_JSON = _Path(__file__).resolve().parent / "calibrated_bands.json"


def load_bands(cell_line: Optional[str]):
    """Return (ctcf_bands, rad21_bands, calibrated: bool, used_cell: str) for a cell line.

    BD2/BD3 fold thresholds are PER-CELL (p75/p90/p97.5 of that cell's own random-window fold).
    calibrated_bands.json holds regenerated measurement bands per cell. If the requested
    cell line is not calibrated, return calibrated=False. Callers MUST NOT grade BD2/BD3
    against the GM12878 cutpoints on a foreign track: that mixes two scales and only
    inflates grades (see tad_vci.insitu_bands). Prefer estimate_bands() in situ, or drop
    the assay so the card reports not_assessable.

    calibrated=True means "these cutpoints were derived from THIS cell line's own random
    windows". An absent or unrecognised cell line must therefore return calibrated=False.
    The tuple still carries GM12878 cutpoints only as a last-resort placeholder that
    annotate_from_beds / mcool_bed refuse to apply unless uncalibrated_bands='gm12878_fallback'."""
    uncal = (CTCF_BANDS, RAD21_BANDS, False, "GM12878")
    if not cell_line:
        return uncal
    try:
        cal = _json.loads(_BANDS_JSON.read_text())
    except Exception:
        return uncal
    # case-insensitive match against the json keys
    key = next((k for k in cal if not k.startswith("_") and k.lower() == str(cell_line).lower()), None)
    if key is None:
        return uncal
    def _b(d):
        # No default for weak_fold. A pre-v2 artifact carries weak_fold = 1.0 (the retired
        # background floor ~50% of random windows cleared); silently defaulting a missing
        # key would resurrect exactly that floor under a v2 label.
        return (("strong", float(d["strong_fold"])), ("moderate", float(d["moderate_fold"])),
                ("weak", float(d["weak_fold"])))
    try:
        return (_b(cal[key]["CTCF"]), _b(cal[key]["RAD21"]), True, key)
    except Exception:
        return (CTCF_BANDS, RAD21_BANDS, False, "GM12878")


def load_null_quantiles(cell_line: Optional[str]) -> Optional[dict]:
    """Return {"CTCF": {"pct": [...], "fold": [...]}, "RAD21": {...}} for a calibrated cell.

    The stored null lets a served card say WHERE an observed fold sits in that cell's
    random-window distribution ("fold 1.065 = 57.7th percentile of random windows"), which
    is the only form in which a band cutpoint is interpretable. Returns None when the cell
    is not calibrated or the artifact predates the stored quantiles.
    """
    if not cell_line:
        return None
    try:
        cal = _json.loads(_BANDS_JSON.read_text())
    except Exception:
        return None
    key = next((k for k in cal if not k.startswith("_") and k.lower() == str(cell_line).lower()), None)
    if key is None:
        return None
    out = {}
    for axis in ("CTCF", "RAD21"):
        rec = cal[key].get(axis) or {}
        pct, fold = rec.get("null_quantiles_pct"), rec.get("null_quantiles")
        if isinstance(pct, list) and isinstance(fold, list) and len(pct) == len(fold) >= 2:
            out[axis] = {"pct": [float(p) for p in pct], "fold": [float(f) for f in fold]}
    return out or None


def null_percentile(fold: Optional[float], quant: Optional[dict]) -> Optional[float]:
    """Percentile of `fold` in a stored random-window null (linear interpolation).

    `quant` is one axis of load_null_quantiles(). Returns None when either is missing.
    Values outside the tabulated range clamp to the end percentiles rather than
    extrapolating a number the null cannot support."""
    if fold is None or quant is None:
        return None
    try:
        f = float(fold)
    except (TypeError, ValueError):
        return None
    if f != f:
        return None
    import numpy as _np
    xs, ys = quant["fold"], quant["pct"]
    return round(float(_np.interp(f, xs, ys)), 2)


@dataclass(frozen=True)
class Criterion:
    code: str
    name: str
    availability: Availability   # how often the required data exists
    role: Role                   # how it enters the evidence-tier rule
    assay: str                   # data required to assess it
    grading: str                 # exact rule mapping the measurement to a grade


# Canonical criteria registry. ONE row per criterion; both axes explicit; grading
# rule stated. This is the single source of truth — BD_SPEC.md and the manuscript
# are generated to match it, and the combining rule below is DERIVED from `role`
# (so the taxonomy and the rule can never silently drift).
CRITERIA: dict[str, Criterion] = {
    "BD1": Criterion("BD1", "cross-caller support", "intrinsic", "primary",
                     "caller boundary sets (votes)",
                     "votes over the provided caller panel (n_methods): "
                     "strong if votes >= ceil(0.8*n), moderate if votes >= ceil(0.6*n), "
                     "else weak (never none for a called boundary). "
                     "A single-method panel (n=1) cannot express cross-caller "
                     "agreement, so BD1 is floored to weak. "
                     "n is the size of the supplied panel and is reported per boundary. "
                     "A method votes if ANY of its calls lies within the position-matching "
                     "tolerance (vote_tol_bp, default 1 bin), so one caller call can "
                     "support several neighbouring candidates: votes are correlated across "
                     "nearby candidates and are not independent observations"),
    "BD2": Criterion("BD2", "CTCF ChIP enrichment", "standard_assay", "supporting",
                     "CTCF ChIP-seq",
                     "+/-25 kb mean fold vs the cell's own random-window background: "
                     ">= p97.5 strong, >= p90 moderate, >= p75 weak, below p75 none "
                     "(percentiles of that cell's random-window fold null)"),
    "BD3": Criterion("BD3", "RAD21 / cohesin enrichment", "standard_assay", "supporting",
                     "RAD21 ChIP-seq",
                     "+/-25 kb mean fold vs the cell's own random-window background: "
                     ">= p97.5 strong, >= p90 moderate, >= p75 weak, below p75 none "
                     "(percentiles of that cell's random-window fold null)"),
    "BD4": Criterion("BD4", "loop-anchor overlap", "standard_assay", "supporting",
                     "loop calls (e.g. HiCCUPS)",
                     "loop anchors within 50 kb: >=2 strong, 1 moderate, 0 none"),
}
# The in-tier supporting set is DERIVED from role, not hard-coded, so it always
# matches the registry. Primary (BD1) is handled explicitly.
IN_TIER_SUPPORTING = [c for c, cr in CRITERIA.items() if cr.role == "supporting"]
SUPPORTING_PLUS_SECONDARY_CORE = IN_TIER_SUPPORTING   # backward-compatible alias


PANEL_METHODS = ["insulation", "topdom_like", "contact_contrast", "spectral_profile", "network_modularity"]


import math as _math


def grade_votes(votes: int, n_methods: int = 5) -> Grade:
    """BD1 over the caller panel, FRACTION-based so it stays comparable when the panel
    size changes (e.g. a user uploads an extra method in 'contribute' mode):
        strong  if votes >= ceil(0.8 * n_methods)
        moderate if votes >= ceil(0.6 * n_methods)
        weak    otherwise (BD1 is never 'none' — it stays assessable so the tier
                engine's "BD1 must be assessable" contract holds even at votes=0).
    With the default 5-method panel this reproduces the original rule exactly
    (>=4 strong, 3 moderate, <=2 weak).

    SINGLE-METHOD FLOOR: with a one-method panel (n_methods < 2) cross-caller
    agreement cannot be expressed — every called boundary is trivially 1/1, which
    is NOT evidence of cross-caller support. BD1 is floored to 'weak' rather than
    reporting a vacuous 'strong'."""
    n = max(1, int(n_methods))
    if n < 2:
        return "weak"
    if votes >= _math.ceil(0.8 * n):
        return "strong"
    if votes >= _math.ceil(0.6 * n):
        return "moderate"
    return "weak"


def grade_fold(fold: Optional[float], bands: tuple) -> Grade:
    if fold is None or pd.isna(fold):
        return "not_assessable"
    for name, cut in bands:
        if fold >= cut:
            return name
    return "none"


def grade_loop(n_anchors: Optional[float]) -> Grade:
    """BD4: boundary coincidence with HiCCUPS loop anchors (CTCF/cohesin corner
    domains). >=2 anchors -> strong, 1 -> moderate, 0 -> none."""
    if n_anchors is None or pd.isna(n_anchors):
        return "not_assessable"
    n = int(n_anchors)
    if n >= 2:
        return "strong"
    if n == 1:
        return "moderate"
    return "none"


_RANK = {"none": 0, "weak": 1, "moderate": 2, "strong": 3}


def _rank(g: Grade) -> int:
    return _RANK.get(g, 0)          # not_assessable -> 0, but excluded upstream


def combine_tier(grades: dict[str, Grade]) -> tuple[Tier, dict]:
    """Transparent evidence-tier rule using only assessable evidence.

    Returns (tier, meta) where meta carries the provenance of which criteria
    were assessed vs not_assessable.

    MONOTONICITY CONTRACT (enforced by test_tier_monotone_in_every_grade): the
    tier is NON-DECREASING in every per-criterion grade — raising any single grade
    (e.g. one supporting criterion weak->moderate->strong, or BD1 weak->moderate->
    strong) can never lower the tier.

    NA-NEUTRALITY: deliberately NOT part of the contract since max_support_v2. A
    supporting axis moving from not_assessable to an assessed `none` CAN lower the tier,
    but only in the single case covered by the depletion clause below (BD1 strong with
    every measured supporting axis depleted). Grade monotonicity still holds throughout:
    within a fixed assessability pattern, raising any grade never lowers the tier. This
    is what makes the missing-vs-negative distinction consequential rather than cosmetic;
    before v2 the two were provably indistinguishable in the tier (measured: 0.00% of
    boundaries changed tier between the two, all 6 released packages).

    To guarantee grade monotonicity, ALL branches aggregate the
    supporting evidence through ONE consistent magnitude statistic `best` (the
    strongest assessable supporting grade). The previous rule mixed `best` with a
    COUNT statistic (`n_mod_plus >= 2`), which made the rule non-monotone: one
    STRONG support (best=strong, count=1) ranked BELOW two MODERATE supports
    (best=moderate, count=2), so upgrading a moderate support to strong could DROP
    the tier. We resolve this by anchoring every branch on `best`, with the
    monotone ordering "one strong support >= two moderate supports" (magnitude
    dominates count). See BD_SPEC.md "Monotonicity property".
    """
    bd1 = grades.get("BD1", "not_assessable")
    if bd1 == "not_assessable":
        raise ValueError("BD1 (caller support) is core and must be assessable for a called boundary")

    assessable = {c: g for c in IN_TIER_SUPPORTING
                  if (g := grades.get(c, "not_assessable")) != "not_assessable"}
    best = max((_rank(g) for g in assessable.values()), default=0)

    if bd1 == "strong" and best >= _RANK["strong"]:
        tier: Tier = "D5"
    elif bd1 == "strong" and assessable and best == 0:
        # NO-SUPPORT CLAUSE (tier_rule_version max_support_v2). Caller agreement alone is
        # not enough for the top two tiers when NO supporting axis that was measured
        # reached even its weak band. This is the one place where "measured and unsupported"
        # is treated differently from "never measured" — the `assessable` guard keeps an
        # unmeasured axis neutral, so a boundary is never penalised for evidence nobody
        # collected.
        # WORDING, corrected 2026-07-29 with the v2 band ladder: `none` on a ChIP axis now
        # means "below the p75 of this cell's random-window null", i.e. NOT ENRICHED, which
        # is weaker than the "depleted below background" this clause used to assume (BD4's
        # `none` was always plain absence: zero loop anchors). The clause is therefore about
        # absence of positive support on every axis we measured, and is stated that way on
        # the card. ASSUMPTION, for the record: a boundary that is unenriched on CTCF,
        # cohesin AND loops is read as unsupported rather than unmeasured. It is a weaker
        # objection than for any single axis (CTCF-independent boundaries exist).
        tier = "D3"
    elif bd1 == "strong" or (bd1 == "moderate" and best >= _RANK["moderate"]):
        tier = "D4"
    elif bd1 == "moderate" or (bd1 == "weak" and best >= _RANK["moderate"]):
        tier = "D3"
    elif bd1 == "weak" and best >= _RANK["weak"]:
        tier = "D2"
    else:
        tier = "D1"

    assessed = ["BD1"] + list(assessable.keys())
    not_assessable = [c for c in CRITERIA if c not in assessed]
    meta = {
        "assessed": assessed,
        "not_assessable": not_assessable,
    }
    return tier, meta


def annotate(df: pd.DataFrame, bg_ctcf: float, bg_rad21: float,
             ctcf_bands: tuple = CTCF_BANDS, rad21_bands: tuple = RAD21_BANDS,
             n_methods: int = 5, null_quantiles: Optional[dict] = None) -> pd.DataFrame:
    """Grade whatever is present in df and combine. Required: votes. Optional
    columns (graded only if present): ctcf, rad21, loop_anchor_count(BD4).
    Missing -> not_assessable. External analysis columns are not graded here and simply
    ride through untouched; the served evidence-card contract carries only BD1-BD4.

    bg_ctcf/bg_rad21 and ctcf_bands/rad21_bands are PER-CELL: pass the calling
    cell's own background means and its own (p75,p90,p97.5) bands from calibrate_bands.py.
    Defaults are the GM12878 calibration.

    null_quantiles (optional), as returned by load_null_quantiles() or built in situ from
    the user's own track, adds ctcf_null_pct / rad21_null_pct: where each observed fold
    sits in that cell's random-window null. Without it a reader cannot tell a fold that
    beats 99% of random windows from one that beats 55% — both used to read simply
    "weak"."""
    out = df.copy()
    # Grade against the panel that could actually vote on each row's chromosome when the
    # caller supplied it (annotate_from_beds writes n_methods_effective); fall back to the
    # scalar for older tables so existing artifacts keep grading identically.
    if "n_methods_effective" in out:
        out["BD1_caller_support"] = [
            grade_votes(int(v), int(n) if n == n and n else n_methods)
            for v, n in zip(out["votes"], out["n_methods_effective"])
        ]
    else:
        out["BD1_caller_support"] = out["votes"].apply(lambda v: grade_votes(v, n_methods))

    nq = null_quantiles or {}
    if "ctcf" in out:
        out["ctcf_fold"] = out["ctcf"] / bg_ctcf
        out["BD2_ctcf"] = out["ctcf_fold"].apply(lambda f: grade_fold(f, ctcf_bands))
        if nq.get("CTCF"):
            out["ctcf_null_pct"] = out["ctcf_fold"].apply(
                lambda f: null_percentile(f, nq["CTCF"]))
    else:
        out["BD2_ctcf"] = "not_assessable"
    if "rad21" in out:
        out["rad21_fold"] = out["rad21"] / bg_rad21
        out["BD3_rad21"] = out["rad21_fold"].apply(lambda f: grade_fold(f, rad21_bands))
        if nq.get("RAD21"):
            out["rad21_null_pct"] = out["rad21_fold"].apply(
                lambda f: null_percentile(f, nq["RAD21"]))
    else:
        out["BD3_rad21"] = "not_assessable"
    out["BD4_loop_anchor"] = (out["loop_anchor_count"].apply(grade_loop)
                              if "loop_anchor_count" in out else "not_assessable")

    tiers, na = [], []
    for _, r in out.iterrows():
        grades = {"BD1": r["BD1_caller_support"], "BD2": r["BD2_ctcf"], "BD3": r["BD3_rad21"],
                  "BD4": r["BD4_loop_anchor"]}
        tier, meta = combine_tier(grades)
        tiers.append(tier)
        na.append(",".join(meta["not_assessable"]))
    out["D_tier"] = tiers
    out["not_assessable_criteria"] = na
    return out
