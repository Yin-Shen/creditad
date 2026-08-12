"""The evidence card must carry what was measured, not only the grade word.

Written after the 2026-07-29 audit, in which a card read `BD1: 2/6 methods` while only one
caller had drawn a block at the position (the other vote came from a call two bins away,
and the tolerance appeared nowhere in the product), and `BD3: weak` on a RAD21 fold of
1.065 that 42% of random windows in the same cell would beat.

Each test below pins one of the disclosures that were missing:
  * the BD1 position-matching tolerance and each supporter's offset;
  * that a vote is not a call at the position, and that neighbouring candidates share votes;
  * the observed fold, the background it was divided by, the cutpoints, and where the fold
    sits in that cell's random-window null.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tad_vci.curation import evidence_card  # noqa: E402
from tad_vci.graded_evidence import load_bands, load_null_quantiles, null_percentile  # noqa: E402

# The boundary from the audit: IMR90 25 kb chr1:45,950,000. votes=2, but only arrowhead
# called this position; OnTAD's vote comes from a call one tolerance-width away.
ROW = {
    "chrom": "1", "pos": 45950000, "D_tier": "D1",
    "BD1_caller_support": "weak", "BD2_ctcf": "none", "BD3_rad21": "none",
    "BD4_loop_anchor": "none",
    "votes": 2, "n_methods": 6, "n_methods_effective": 6,
    "supporting_callers": "OnTAD;arrowhead",
    "supporting_offsets_bp": "50000;0",
    "vote_tol_bp": 50000, "resolution_bp": 25000,
    "bd1_vote_model": "nearest_within_tol",
    "ctcf_fold": 0.887133, "rad21_fold": 1.065147,
    "bg_ctcf": 0.6141, "bg_rad21": 0.6142,
    "ctcf_null_pct": 33.41, "rad21_null_pct": 57.74,
    "loop_anchor_count": 0, "chip_window_bp": 50000,
    "data_cell_line": "IMR90", "bands_calibrated": True,
}


def test_bd1_reports_the_match_window_and_every_supporter_offset():
    card = evidence_card(ROW)
    m = card["criteria"]["BD1"]["measured"]
    assert m["vote_tol_bp"] == 50000 and m["vote_tol_bins"] == 2.0
    assert card["vote_tol_bp"] == 50000
    assert m["supporters"] == [
        {"caller": "OnTAD", "offset_bp": 50000, "offset_bins": 2.0},
        {"caller": "arrowhead", "offset_bp": 0, "offset_bins": 0.0},
    ]
    # the headline the audit found missing: 2 votes, 1 call here
    assert m["n_supporters_at_this_position"] == 1


def test_bd1_states_that_votes_are_position_matched_and_reused():
    card = evidence_card(ROW)
    rule = card["bd1_match_rule"]
    assert "50,000 bp" in rule and "2 bins" in rule
    assert "not a call at this position" in rule
    # the union candidate set shares one call between neighbouring candidates
    assert "correlated" in rule and "independent" in rule


def test_bd1_match_rule_distinguishes_the_clustered_vote_model():
    clustered = dict(ROW, bd1_vote_model="cluster_consensus")
    rule = evidence_card(clustered)["bd1_match_rule"]
    assert "at most one vote per candidate" in rule
    assert "correlated" not in rule


def test_bd1_measured_survives_an_annotation_without_offsets():
    """Pre-2026-07-29 tables have supporting_callers but no offsets: name the callers,
    say the offset is unknown, and do not invent a zero (which would read as
    'called this exact position')."""
    old = {k: v for k, v in ROW.items() if k != "supporting_offsets_bp"}
    m = evidence_card(old)["criteria"]["BD1"]["measured"]
    assert m["supporters"] == [{"caller": "OnTAD"}, {"caller": "arrowhead"}]
    assert m["n_supporters_at_this_position"] == 0   # none are known to be at the position


def test_chip_axes_report_fold_background_cutpoints_and_null_rank():
    card = evidence_card(ROW)
    bd3 = card["criteria"]["BD3"]["measured"]
    assert bd3["fold"] == round(ROW["rad21_fold"], 4)
    assert bd3["background_mean"] == ROW["bg_rad21"]
    assert bd3["window_bp"] == 50000
    _, rad21_bands, calibrated, cell = load_bands("IMR90")
    assert calibrated and cell == "IMR90"
    assert bd3["cutpoints"] == {name: round(cut, 4) for name, cut in dict(rad21_bands).items()}
    assert bd3["cutpoint_percentiles"] == {"weak": 75, "moderate": 90, "strong": 97.5}
    # the number that makes the grade legible: 42% of random windows reach this fold
    assert bd3["null_percentile"] == 57.74
    assert bd3["random_windows_at_or_above"] == round(100 - 57.74, 2)


def test_the_audited_fold_grades_none_and_sits_below_the_weak_cut():
    """Regression on the v1 band ladder: RAD21 fold 1.065 in IMR90 used to grade `weak`
    (the floor was fold >= 1.0, ~the median of the null). It must now be `none`, and the
    card must show why."""
    card = evidence_card(ROW)
    assert card["criteria"]["BD3"]["grade"] == "none"
    m = card["criteria"]["BD3"]["measured"]
    assert m["fold"] < m["cutpoints"]["weak"]
    assert m["fold"] > 1.0          # still above background — `none` is not `depleted`
    q = load_null_quantiles("IMR90")
    assert abs(null_percentile(ROW["rad21_fold"], q["RAD21"]) - m["null_percentile"]) < 0.01


def test_grading_text_quotes_the_percentile_each_cutpoint_is():
    card = evidence_card(ROW)
    for code in ("BD2", "BD3"):
        text = card["criteria"][code]["grading"]
        assert "p97.5" in text and "p90" in text and "p75" in text
        assert ">=1.00 weak" not in text        # the retired background floor
    assert card["criteria"]["BD4"]["measured"] == {
        "loop_anchor_count": 0, "window_bp": 50000,
        "cutpoints": {"strong": 2, "moderate": 1, "none": 0},
    }


def test_no_measurement_block_when_the_annotation_has_no_numbers():
    """An older or partial table must not grow an empty measurement block that reads as
    'measured and zero'."""
    bare = {"chrom": "1", "pos": 1000, "D_tier": "D1", "BD1_caller_support": "weak",
            "BD2_ctcf": "not_assessable", "BD3_rad21": "not_assessable",
            "BD4_loop_anchor": "not_assessable"}
    card = evidence_card(bare)
    assert "measured" not in card["criteria"]["BD2"]
    assert "measured" not in card["criteria"]["BD4"]
    assert card["criteria"]["BD1"]["measured"]["supporters"] is None
    assert card["bd1_match_rule"] is None
