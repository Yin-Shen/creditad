"""Unit tests for the TAD-VCI graded-evidence engine (3-tier ACMG-style)."""
from __future__ import annotations

import sys

import pandas as pd
import pytest

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2]))
from tad_vci import (  # noqa: E402
    CRITERIA, annotate, combine_tier, grade_fold, grade_votes,
    CTCF_BANDS, RAD21_BANDS,
)
NA = "not_assessable"
TIER_ORDER = {"D1": 1, "D2": 2, "D3": 3, "D4": 4, "D5": 5}


def g(b1, b2=NA, b3=NA, b4=NA):
    return {"BD1": b1, "BD2": b2, "BD3": b3, "BD4": b4}


def tier(b1, b2=NA, b3=NA, b4=NA):
    return combine_tier(g(b1, b2, b3, b4))[0]


# ---- per-criterion grading ----
def test_grade_votes():
    # 5-caller panel: >=4 strong, 3 moderate, <=2 weak
    assert grade_votes(5) == "strong"
    assert grade_votes(4) == "strong"
    assert grade_votes(3) == "moderate"
    assert grade_votes(2) == "weak"
    assert grade_votes(1) == "weak"


def test_grade_fold_bands():
    # Band ladder v2: every cut is a percentile of the cell's own random-window fold null
    # (weak p75, moderate p90, strong p97.5). Assert against the live constants so the test
    # tracks the single source of truth rather than a frozen literal.
    strong_cut = CTCF_BANDS[0][1]
    mod_cut = CTCF_BANDS[1][1]
    weak_cut = CTCF_BANDS[2][1]
    assert grade_fold(strong_cut + 1.0, CTCF_BANDS) == "strong"
    assert grade_fold(strong_cut, CTCF_BANDS) == "strong"           # exactly the strong cut
    assert grade_fold(strong_cut - 1e-3, CTCF_BANDS) == "moderate"  # just below strong
    assert grade_fold(mod_cut, CTCF_BANDS) == "moderate"
    assert grade_fold(mod_cut - 1e-3, CTCF_BANDS) == "weak"
    assert grade_fold(weak_cut, CTCF_BANDS) == "weak"
    # The v1 floor: "at or above background" is roughly the MEDIAN of the null, so it must
    # NOT reach the weak band any more. This is the regression that retired v1.
    assert weak_cut > 1.0
    assert grade_fold(1.00, CTCF_BANDS) == "none"
    assert grade_fold(weak_cut - 1e-3, CTCF_BANDS) == "none"
    assert grade_fold(0.80, CTCF_BANDS) == "none"


def test_grade_fold_missing_is_not_assessable():
    assert grade_fold(None, RAD21_BANDS) == NA
    assert grade_fold(float("nan"), RAD21_BANDS) == NA


def test_grade_conservation_removed_from_engine():
    """BD5 cross-cell conservation was removed from the core engine 2026-06-04; the grader
    must NOT be importable from the tad_vci package (it lives only in the analysis module)."""
    import tad_vci
    assert not hasattr(tad_vci, "grade_conservation")
    assert "grade_conservation" not in tad_vci.__all__


# ---- combining rule: ACMG branches ----
def test_tier_branches():
    assert tier("strong", b2="strong") == "D5"
    # max_support_v2 depletion clause: strong caller agreement with every MEASURED
    # supporting axis depleted no longer reaches D4 (it was D4 up to v1).
    assert tier("strong", b2="none", b3="none") == "D3"
    assert tier("strong", b2="none", b3="weak") == "D4"   # any real support restores D4
    assert tier("moderate", b2="strong") == "D4"
    assert tier("moderate", b2="weak") == "D3"
    assert tier("weak", b2="moderate", b3="moderate") == "D3"   # two moderate supports
    assert tier("weak", b2="strong") == "D3"                    # one STRONG support (>= 2 moderate)
    assert tier("weak", b2="moderate") == "D3"                  # one moderate support lifts weak-BD1
    assert tier("weak", b2="weak") == "D2"
    assert tier("weak", b2="none", b3="none") == "D1"


def test_tier_monotone_in_every_grade():
    """The tier must be NON-DECREASING in every per-criterion grade: raising any single
    grade one step up the ordered scale (none<weak<moderate<strong, and BD1 weak<moderate<
    strong) can never lower the tier. Also: a not_assessable supporting criterion becoming
    an assessable depleted ('none') one must never lower the tier. This is the monotonicity
    contract documented in BD_SPEC.md; it directly rules out the prior non-monotone bug where
    one STRONG support ranked below two MODERATE supports."""
    GRADES = ["none", "weak", "moderate", "strong", NA]
    SUP_SCALE = ["none", "weak", "moderate", "strong"]
    BD1_SCALE = ["weak", "moderate", "strong"]   # BD1 is never none/NA for a called boundary

    def t(b1, b2, b3, b4):
        return TIER_ORDER[combine_tier({"BD1": b1, "BD2": b2, "BD3": b3, "BD4": b4})[0]]

    checked = 0
    for b1 in BD1_SCALE:
        for b2 in GRADES:
            for b3 in GRADES:
                for b4 in GRADES:
                    base = t(b1, b2, b3, b4)
                    # raise BD1 one step
                    i = BD1_SCALE.index(b1)
                    if i + 1 < len(BD1_SCALE):
                        assert t(BD1_SCALE[i + 1], b2, b3, b4) >= base
                        checked += 1
                    # raise each supporting criterion one step (or NA -> none)
                    for idx, bx in ((0, b2), (1, b3), (2, b4)):
                        args = [b2, b3, b4]
                        if bx in SUP_SCALE and SUP_SCALE.index(bx) + 1 < len(SUP_SCALE):
                            args[idx] = SUP_SCALE[SUP_SCALE.index(bx) + 1]
                            assert t(b1, *args) >= base
                            checked += 1
                        elif bx == NA:
                            # Since max_support_v2 an axis becoming assessable-but-depleted
                            # may LOWER the tier, but ONLY via the depletion clause: BD1
                            # strong with every measured supporting axis depleted. Anywhere
                            # else NA -> none must still be neutral or better.
                            args[idx] = "none"
                            after = t(b1, *args)
                            others = [g for j, g in enumerate(args) if j != idx]
                            depletion_case = (
                                b1 == "strong"
                                and all(g in ("none", NA) for g in others)
                            )
                            if depletion_case:
                                assert after <= base, "depletion clause may only lower or hold"
                            else:
                                assert after >= base, (
                                    f"NA->none lowered the tier outside the depletion clause: "
                                    f"BD1={b1} args={args} {base}->{after}")
                            checked += 1
    assert checked > 1000   # full single-step coverage of the grade lattice


def test_depletion_clause_distinguishes_measured_depletion_from_missing():
    """max_support_v2: caller agreement alone cannot reach D4/D5 once every supporting
    axis that was measured came back depleted — while an axis nobody measured stays
    neutral. Before v2 these two states were indistinguishable in the tier."""
    strong_only = {"BD1": "strong", "BD2": NA, "BD3": NA, "BD4": NA}
    all_depleted = {"BD1": "strong", "BD2": "none", "BD3": "none", "BD4": "none"}
    one_measured_depleted = {"BD1": "strong", "BD2": "none", "BD3": NA, "BD4": NA}
    partial_support = {"BD1": "strong", "BD2": "none", "BD3": "weak", "BD4": "none"}

    assert combine_tier(strong_only)[0] == "D4", "never-measured must not be penalised"
    assert combine_tier(all_depleted)[0] == "D3", "measured-and-empty must not reach D4"
    assert combine_tier(one_measured_depleted)[0] == "D3", "the only measured axis is empty"
    assert combine_tier(partial_support)[0] == "D4", "any real support keeps D4"
    # the distinction must be visible, i.e. the two states differ
    assert combine_tier(strong_only)[0] != combine_tier(all_depleted)[0]


def test_one_strong_support_at_least_two_moderate():
    """The chosen monotonicity ordering: one STRONG support is at least as good as two
    MODERATE supports (magnitude dominates count). Was VIOLATED before the fix."""
    assert TIER_ORDER[tier("weak", b2="strong")] >= TIER_ORDER[tier("weak", b2="moderate", b3="moderate")]
    assert TIER_ORDER[tier("moderate", b2="strong")] >= TIER_ORDER[tier("moderate", b2="moderate", b3="moderate")]


def test_absence_is_not_negative():
    # weak BD1 with all support not_assessable -> D1 (not penalized below floor)
    assert tier("weak") == "D1"
    # strong BD1 with no support still reaches D4 (absence != downgrade)
    assert tier("strong") == "D4"


def test_bd1_must_be_assessable():
    with pytest.raises(ValueError):
        combine_tier(g(NA, b2="strong"))


def test_provenance_lists_not_assessable():
    _, meta = combine_tier(g("moderate", b2="strong"))
    assert "BD1" in meta["assessed"] and "BD2" in meta["assessed"]
    assert set(["BD3", "BD4"]).issubset(set(meta["not_assessable"]))


# ---- criterion registry (availability classes) ----
def test_criteria_registry_two_axes():
    # axis 1: availability
    assert CRITERIA["BD1"].availability == "intrinsic"
    assert CRITERIA["BD2"].availability == "standard_assay"
    assert CRITERIA["BD3"].availability == "standard_assay"
    assert CRITERIA["BD4"].availability == "standard_assay"
    # axis 2: role (drives the rule)
    assert CRITERIA["BD1"].role == "primary"
    assert CRITERIA["BD2"].role == CRITERIA["BD3"].role == CRITERIA["BD4"].role == "supporting"
    # the in-tier supporting set is DERIVED from role and must match exactly
    from tad_vci.graded_evidence import IN_TIER_SUPPORTING
    assert set(IN_TIER_SUPPORTING) == {"BD2", "BD3", "BD4"}
    assert "BD5" not in CRITERIA          # retired motif experiment is outside the evidence-card contract
    assert "BD6" not in CRITERIA          # degron criterion removed (no data on any cell line)
    assert len(CRITERIA) == 4


def test_bands_match_calibration():
    """Engine bands MUST equal the reproducible GM12878 calibration (single source
    of truth: calibrate_bands.py -> calibrated_bands.json). Prevents engine/calibration drift."""
    import json
    from importlib.resources import files
    j = json.loads(files("tad_vci").joinpath("calibrated_bands.json").read_text())["GM12878"]
    assert CTCF_BANDS[0][1] == j["CTCF"]["strong_fold"] and CTCF_BANDS[1][1] == j["CTCF"]["moderate_fold"]
    assert RAD21_BANDS[0][1] == j["RAD21"]["strong_fold"] and RAD21_BANDS[1][1] == j["RAD21"]["moderate_fold"]
    assert CTCF_BANDS[2][1] == j["CTCF"]["weak_fold"] and RAD21_BANDS[2][1] == j["RAD21"]["weak_fold"]
    # v2: the weak rung is the null's p75, not the background floor. A stored 1.0 would
    # mean a pre-v2 artifact had been shipped under a v2 engine.
    for axis in ("CTCF", "RAD21"):
        assert j[axis]["band_percentiles"] == {"weak": 75.0, "moderate": 90.0, "strong": 97.5}
        assert j[axis]["weak_fold"] > 1.0
        # the measurement that retired the v1 floor must stay in the artifact
        assert 0.3 < j[axis]["p_random_at_or_above_background"] < 0.7
        assert len(j[axis]["null_quantiles"]) == len(j[axis]["null_quantiles_pct"]) == 101


# ---- annotate end-to-end ----
def test_annotate_minimal_and_optional_columns():
    df = pd.DataFrame({
        "chrom": ["1", "1", "1"], "pos": [100000, 200000, 300000],
        "votes": [4, 3, 1], "ctcf": [1.2, 0.8, 0.4], "rad21": [1.5, 0.9, 0.7],
    })
    out = annotate(df, bg_ctcf=0.5, bg_rad21=0.9)
    # 5-method panel, v2 GM12878 bands (CTCF weak 1.198 / moderate 1.564 / strong 2.141):
    #   row1 votes4=strong + ctcf fold 2.4=strong                    -> D5
    #   row2 votes3=moderate + ctcf fold 1.6=moderate                -> D4
    #   row3 votes1=weak + ctcf 0.8=none + rad21 0.78=none           -> D1
    assert list(out["BD2_ctcf"]) == ["strong", "moderate", "none"]
    assert list(out["D_tier"]) == ["D5", "D4", "D1"]
    # BD4 not wired here -> listed not_assessable
    assert "BD4" in out["not_assessable_criteria"].iloc[0]
    assert "BD5" not in out.columns and "BD5_ctcf_motif" not in out.columns   # BD5 removed 2026-06-04
    assert "BD6" not in out.columns and "perturbation_confirmation" not in out.columns
    # cross-cell conservation is NO LONGER graded by the engine (moved to the analysis module).
    # A conservation_count input column just rides through untouched; NO cross_cell_conservation
    # grade column is produced, and the tier is unaffected by it.
    df3 = df.assign(conservation_count=[3, 0, 2])
    out3 = annotate(df3, 0.5, 0.9)
    assert "cross_cell_conservation" not in out3.columns
    assert out3["conservation_count"].tolist() == [3, 0, 2]   # input column preserved
    assert list(out3["D_tier"]) == list(out["D_tier"])        # conservation does not move the tier


# ---- served per-boundary TSV exposes exactly BD1-BD4 (no BD5/conservation re-introduction) ----
def test_served_tsv_header_is_bd1_to_bd4_only():
    """The per-boundary annotation the engine emits (what gets written to the served TSV)
    must carry exactly the BD1-BD4 grade columns and never a BD5/conservation grade column,
    so an accidental BD5 re-introduction is caught here."""
    df = pd.DataFrame({
        "chrom": ["1", "1"], "pos": [100000, 200000], "votes": [4, 2],
        "ctcf": [1.0, 0.6], "rad21": [1.5, 0.9], "loop_anchor_count": [2, 0],
        "conservation_count": [3, 1],   # passed in but must NOT yield a grade column
    })
    out = annotate(df, bg_ctcf=0.5, bg_rad21=0.9)
    grade_cols = {"BD1_caller_support", "BD2_ctcf", "BD3_rad21", "BD4_loop_anchor"}
    present_bd = {c for c in out.columns if c.startswith("BD")}
    assert present_bd == grade_cols
    # no BD5/conservation/perturbation/degron grade columns anywhere
    forbidden = ["BD5", "BD6", "cross_cell_conservation", "BD5_cross_cell",
                 "BD5_ctcf_motif", "perturbation_confirmation", "BD6_degron"]
    assert all(f not in out.columns for f in forbidden)
    # the evidence-card column map (curation.py) must agree with this BD1-BD4 set
    from tad_vci.curation import _CARD_COLS
    assert set(_CARD_COLS.values()) == grade_cols
    assert set(_CARD_COLS.keys()) == {"BD1", "BD2", "BD3", "BD4"}


def test_grade_loop():
    from tad_vci import grade_loop
    assert grade_loop(2) == "strong"
    assert grade_loop(5) == "strong"
    assert grade_loop(1) == "moderate"
    assert grade_loop(0) == "none"
    assert grade_loop(None) == NA
    assert grade_loop(float("nan")) == NA


def test_bd1_denominator_is_the_panel_that_could_vote_on_that_chromosome():
    """A caller that produced no boundaries on a chromosome cannot vote there, so counting
    it in the BD1 denominator deflates every grade on that chromosome. annotate() must use
    the per-row n_methods_effective when the caller supplies it (7,321 candidates across 4
    chromosomes were affected in the shipped packages before this was fixed), and fall back
    to the scalar when the column is absent."""
    import pandas as pd
    from tad_vci.graded_evidence import annotate

    # 4 votes out of a 6-caller panel = moderate; out of the 5 that could vote = strong.
    df = pd.DataFrame([
        {"chrom": "7", "pos": 1_000_000, "votes": 4, "n_methods": 6, "n_methods_effective": 5},
        {"chrom": "8", "pos": 2_000_000, "votes": 4, "n_methods": 6, "n_methods_effective": 6},
    ])
    out = annotate(df, bg_ctcf=1.0, bg_rad21=1.0, n_methods=6)
    assert out.loc[0, "BD1_caller_support"] == "strong", "chr7 lost a caller: denominator is 5"
    assert out.loc[1, "BD1_caller_support"] == "moderate", "chr8 had the full panel"

    # Backward compatibility: no effective column -> scalar path, unchanged behaviour.
    legacy = pd.DataFrame([{"chrom": "7", "pos": 1_000_000, "votes": 4, "n_methods": 6}])
    out2 = annotate(legacy, bg_ctcf=1.0, bg_rad21=1.0, n_methods=6)
    assert out2.loc[0, "BD1_caller_support"] == "moderate"


def test_tier_why_string_always_states_the_tier_the_rule_produced():
    """_tier_why is written into the signed decision record, so a mismatch between the
    stated derivation and the actual tier is a provenance defect, not a cosmetic one.
    Before max_support_v2 was wired into it, a depleted boundary was described as '-> D4'
    while the rule returned D3."""
    import re
    from tad_vci.curation import _tier_why

    cases = [
        {"BD1": "strong", "BD2": "none", "BD3": "none", "BD4": "none"},      # depletion -> D3
        {"BD1": "strong", "BD2": NA, "BD3": NA, "BD4": NA},                  # unmeasured -> D4
        {"BD1": "strong", "BD2": "none", "BD3": "weak", "BD4": "none"},      # partial -> D4
        {"BD1": "strong", "BD2": "strong", "BD3": NA, "BD4": NA},            # -> D5
        {"BD1": "moderate", "BD2": "moderate", "BD3": NA, "BD4": NA},        # -> D4
        {"BD1": "weak", "BD2": "none", "BD3": "none", "BD4": NA},            # -> D1
    ]
    for grades in cases:
        tier = combine_tier(grades)[0]
        why = _tier_why(grades)
        m = re.search(r"→\s*(D[1-5])", why)
        assert m, f"derivation string states no tier: {why!r}"
        assert m.group(1) == tier, f"{why!r} states {m.group(1)} but the rule returned {tier}"
