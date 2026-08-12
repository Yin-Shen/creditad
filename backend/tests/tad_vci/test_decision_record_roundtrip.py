"""Decision-record round trip and tamper-evidence — closes inventory item H5.

Full path exercised here: ``evidence_card -> make_snapshot -> append_snapshot ->
load_history -> verify_snapshot``, on a JSONL log inside ``tmp_path``. The repository's
own ``backend/tad_vci/decision_records_v3.jsonl`` is never written to; one test asserts
that explicitly, because a test that appended to it would both dirty the shipped clean
state and pass for the wrong reason.

What the checksum does and does not claim is part of the contract: it detects a change
to any checksummed field of a record you already hold. It does NOT prove the log is
complete, and ``actor_label`` is a caller-supplied string, not an authenticated identity.
The tests are written to pin exactly that scope.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tad_vci.curation import (  # noqa: E402
    CurationSnapshot, append_snapshot, evidence_card, load_history, make_snapshot,
    next_version, verify_snapshot,
)

REPO_LOG = Path(__file__).resolve().parents[2] / "tad_vci" / "decision_records_v3.jsonl"

ROW = {
    "chrom": "7", "pos": 1_250_000, "D_tier": "D4",
    "BD1_caller_support": "strong", "BD2_ctcf": "moderate",
    "BD3_rad21": "weak", "BD4_loop_anchor": "not_assessable",
    "votes": 4, "n_methods": 5, "data_cell_line": "GM12878", "bands_calibrated": True,
}
TS = "2026-07-26T08:00:00Z"

# Every field that make_snapshot folds into the checksum payload, except hash_schema —
# that one selects the verification branch and is covered by its own test below.
CHECKSUMMED_FIELDS = (
    "boundary_id", "assembly", "auto_tier", "final_tier", "evidence",
    "review_decision", "version", "actor_label", "created_at",
    "resolution", "cell_line",
)


@pytest.fixture
def log_path(tmp_path: Path) -> Path:
    return tmp_path / "decision_records_v3.jsonl"


def _appended(log: Path, **kwargs) -> CurationSnapshot:
    card = evidence_card(ROW)
    snap = make_snapshot(card, "hg38", "reviewer-A", TS, resolution=10_000,
                         cell_line="GM12878", **kwargs)
    append_snapshot(snap, log)
    return snap


def test_round_trip_reloads_an_identical_verifying_record(log_path: Path):
    snap = _appended(log_path)

    hist = load_history("7:1250000", log_path)
    assert len(hist) == 1
    rec = hist[0]

    assert verify_snapshot(rec) is True
    assert rec["content_hash"] == snap.content_hash
    assert len(rec["content_hash"]) == 64
    assert rec["hash_schema"] == "creditad-decision-checksum-v3"
    assert rec["version"] == 1
    assert rec["auto_tier"] == "D4" and rec["final_tier"] == "D4"
    assert rec["review_decision"] is None
    assert rec["assembly"] == "hg38"
    assert rec["resolution"] == 10_000 and rec["cell_line"] == "GM12878"
    assert set(rec["evidence"]) == {"BD1", "BD2", "BD3", "BD4"}
    assert rec["evidence"]["BD1"]["grade"] == "strong"

    # one JSON object per line, byte-identical to what the snapshot serialises
    lines = log_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1 and lines[0] == snap.to_json()
    assert json.loads(lines[0]) == rec


def test_reviewer_override_appends_a_second_verifying_version(log_path: Path):
    _appended(log_path)
    card = evidence_card(ROW)
    prior = next_version(card["boundary_id"], log_path)
    assert prior == 1

    snap2 = make_snapshot(
        card, "hg38", "reviewer-B", "2026-07-26T09:00:00Z",
        decision={"final_tier": "D5", "verdict": "supported_for_review",
                  "rationale": "convergent CTCF motif pair at the anchor"},
        prior_version=prior, resolution=10_000, cell_line="GM12878",
    )
    append_snapshot(snap2, log_path)

    hist = load_history(card["boundary_id"], log_path)
    assert [h["version"] for h in hist] == [1, 2]
    assert all(verify_snapshot(h) for h in hist)
    assert hist[0]["final_tier"] == "D4" and hist[1]["final_tier"] == "D5"
    assert hist[1]["auto_tier"] == "D4", "the automatic tier must survive the override"
    assert hist[0]["content_hash"] != hist[1]["content_hash"]


@pytest.mark.parametrize("field", CHECKSUMMED_FIELDS)
def test_tampering_with_any_checksummed_field_fails_verification(field: str, log_path: Path):
    _appended(log_path)
    rec = load_history("7:1250000", log_path)[0]
    assert verify_snapshot(rec) is True

    original = rec[field]
    if field == "evidence":
        rec = json.loads(json.dumps(rec))
        rec["evidence"]["BD2"]["grade"] = "strong"          # quiet upgrade of one axis
    elif field == "review_decision":
        rec["review_decision"] = {"final_tier": "D5", "verdict": "supported_for_review",
                                  "rationale": "inserted after the fact"}
    elif field == "version":
        rec["version"] = 99
    elif field == "resolution":
        rec["resolution"] = 25_000
    else:
        rec[field] = str(original) + "-tampered"

    assert rec[field] != original
    assert verify_snapshot(rec) is False, f"tampering with {field!r} went undetected"


@pytest.mark.parametrize(
    "schema", ["creditad-decision-checksum-v2", "creditad-decision-checksum-v1", "", None]
)
def test_downgrading_the_hash_schema_never_verifies(schema, log_path: Path):
    """`hash_schema` selects the verification branch, so it is the field an attacker would
    rewrite to reach a weaker checksum. The property that matters is that a v3 record with
    a rewritten schema NEVER verifies.

    Current implementation detail, asserted as-is rather than papered over: the legacy
    branches index `rec["expert_override"]` / `rec["curator"]`, which a v3 record does not
    carry, so the call raises KeyError instead of returning False. Either outcome is a
    refusal to verify; it is never a silent True.
    """
    _appended(log_path)
    rec = load_history("7:1250000", log_path)[0]
    rec["hash_schema"] = schema

    try:
        assert verify_snapshot(rec) is False
    except KeyError as exc:
        assert str(exc).strip("'") in {"expert_override", "curator"}


def test_tampering_with_the_checksum_itself_fails_verification(log_path: Path):
    _appended(log_path)
    rec = load_history("7:1250000", log_path)[0]
    rec["content_hash"] = "0" * 64
    assert verify_snapshot(rec) is False


def test_history_is_scoped_to_one_boundary_and_ordered_by_version(log_path: Path):
    _appended(log_path)
    other = evidence_card({**ROW, "pos": 2_000_000})
    append_snapshot(make_snapshot(other, "hg38", "reviewer-A", TS), log_path)

    assert [h["boundary_id"] for h in load_history("7:1250000", log_path)] == ["7:1250000"]
    assert [h["boundary_id"] for h in load_history("7:2000000", log_path)] == ["7:2000000"]
    assert load_history("7:9999999", log_path) == []
    assert next_version("7:9999999", log_path) == 0


def test_round_trip_never_writes_to_the_repository_decision_log(log_path: Path):
    """The shipped log must stay at the clean state it is released in."""
    before = (REPO_LOG.stat().st_size, REPO_LOG.stat().st_mtime_ns) if REPO_LOG.is_file() else None

    _appended(log_path)
    assert log_path.is_file() and log_path.stat().st_size > 0

    after = (REPO_LOG.stat().st_size, REPO_LOG.stat().st_mtime_ns) if REPO_LOG.is_file() else None
    assert after == before
    if REPO_LOG.is_file():
        assert REPO_LOG.stat().st_size == 0, (
            "backend/tad_vci/decision_records_v3.jsonl is shipped empty; a non-empty file "
            "here means a run leaked local review records into the release tree"
        )


def test_make_snapshot_refuses_an_unusable_review_record():
    """An override without a rationale, an override without a tier, or an unlabelled
    actor would all produce a record nobody can audit later."""
    card = evidence_card(ROW)
    with pytest.raises(ValueError):
        make_snapshot(card, "hg38", "reviewer-A", TS,
                      decision={"final_tier": "D5", "rationale": "   "})
    with pytest.raises(ValueError):
        make_snapshot(card, "hg38", "reviewer-A", TS,
                      decision={"final_tier": "", "rationale": "looks right"})
    with pytest.raises(ValueError):
        make_snapshot(card, "hg38", "   ", TS)
