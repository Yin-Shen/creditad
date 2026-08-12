"""Tests for evidence cards and application-level versioned records."""
from __future__ import annotations

import sys

import pytest

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2]))
from tad_vci.curation import (  # noqa: E402
    append_snapshot, content_hash, evidence_card, load_history, make_snapshot,
    next_version, verify_snapshot,
)

ROW = {
    "chrom": "1", "pos": 1000000, "D_tier": "D4",
    "BD1_caller_support": "strong", "BD2_ctcf": "moderate", "BD3_rad21": "weak",
    "BD4_loop_anchor": "moderate",
}
TS = "2026-05-30T00:00:00Z"


def test_evidence_card_structure():
    card = evidence_card(ROW)
    assert card["boundary_id"] == "1:1000000"
    assert card["auto_tier"] == "D4"
    assert card["criteria"]["BD1"]["role"] == "primary"
    assert card["criteria"]["BD1"]["availability"] == "intrinsic"
    assert "BD5" not in card["criteria"] and len(card["criteria"]) == 4   # BD5 removed 2026-06-04
    assert "BD1" in card["assessed"]


def test_snapshot_without_override():
    card = evidence_card(ROW)
    s = make_snapshot(card, "hg19", "reviewerA", TS)
    assert s.version == 1 and s.final_tier == "D4" and s.review_decision is None
    assert len(s.content_hash) == 64
    assert s.hash_schema == "creditad-decision-checksum-v3"


def test_override_requires_rationale_and_tier():
    card = evidence_card(ROW)
    with pytest.raises(ValueError):
        make_snapshot(card, "hg19", "c", TS, decision={"final_tier": "D5", "rationale": "   "})
    with pytest.raises(ValueError):
        make_snapshot(card, "hg19", "c", TS, decision={"final_tier": "", "rationale": "review note"})


def test_application_append_preserves_prior_line_and_history(tmp_path):
    log = tmp_path / "curation.jsonl"
    card = evidence_card(ROW)
    s1 = make_snapshot(card, "hg19", "reviewerA", TS, prior_version=0)
    append_snapshot(s1, log)
    first_line = log.read_text().splitlines()[0]            # capture v1 line verbatim

    # local reviewer changes the tier with a rationale -> new version, old line untouched
    v = next_version(card["boundary_id"], log)
    s2 = make_snapshot(card, "hg19", "reviewerB", "2026-05-30T01:00:00Z",
                       decision={"final_tier": "D5", "verdict": "supported_for_review",
                                 "rationale": "flanking CTCF peak + known enhancer"},
                       prior_version=v)
    append_snapshot(s2, log)

    assert log.read_text().splitlines()[0] == first_line     # this API did not rewrite v1
    hist = load_history(card["boundary_id"], log)
    assert [h["version"] for h in hist] == [1, 2]
    assert hist[0]["final_tier"] == "D4" and hist[1]["final_tier"] == "D5"
    assert hist[1]["review_decision"]["rationale"].startswith("flanking")
    assert hist[0]["content_hash"] != hist[1]["content_hash"]


def test_verify_snapshot_detects_selected_record_mismatch(tmp_path):
    log = tmp_path / "c.jsonl"
    card = evidence_card(ROW)
    append_snapshot(make_snapshot(card, "hg19", "c", TS), log)
    rec = load_history(card["boundary_id"], log)[0]
    assert verify_snapshot(rec) is True
    rec["final_tier"] = "D5"
    assert verify_snapshot(rec) is False


def test_v3_checksum_covers_actor_time_and_context():
    card = evidence_card(ROW)
    rec = make_snapshot(card, "hg19", "reviewer_label", TS,
                        resolution=25000, cell_line="GM12878").__dict__.copy()
    assert verify_snapshot(rec) is True
    for field, changed in (
        ("actor_label", "different_label"),
        ("created_at", "2026-05-31T00:00:00Z"),
        ("resolution", 50000),
        ("cell_line", "K562"),
    ):
        altered = dict(rec)
        altered[field] = changed
        assert verify_snapshot(altered) is False


def test_v2_record_remains_verifiable_for_backward_reading_only():
    card = evidence_card(ROW)
    payload = {
        "boundary_id": card["boundary_id"], "assembly": "hg19",
        "auto_tier": card["auto_tier"], "final_tier": card["auto_tier"],
        "evidence": card["criteria"], "expert_override": None,
        "version": 1, "curator": "legacy_demo_label", "created_at": TS,
        "resolution": 25000, "cell_line": "GM12878",
        "hash_schema": "creditad-decision-checksum-v2",
    }
    rec = {**payload, "content_hash": content_hash(payload)}
    assert verify_snapshot(rec) is True


def test_two_boundaries_isolated(tmp_path):
    log = tmp_path / "c.jsonl"
    for chrom in ("1", "2"):
        row = dict(ROW, chrom=chrom)
        append_snapshot(make_snapshot(evidence_card(row), "hg19", "c", TS), log)
    assert len(load_history("1:1000000", log)) == 1
    assert len(load_history("2:1000000", log)) == 1
