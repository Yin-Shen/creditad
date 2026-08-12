"""Scope-aware decision-record history — closes inventory item B14.

A boundary id is only ``chrom:pos``. The same string is produced by every dataset:
GM12878 hg38 10 kb and K562 hg38 10 kb both contain ``3:1234567``. Before this change
``GET /api/tadvci/history/3:1234567`` returned every cell line's record for that
coordinate, so the evidence panel showed a colleague's GM12878 verdict while the user
was looking at K562 — the cross-cell collision the inventory flags.

The lookup is now qualified with an ``@cell/assembly/resolution`` suffix that
``curation.parse_boundary_selector`` understands, and the front end sends it
(``TADVCIEvidence.vue``: ``historyScope``). This file pins both halves:

  * the pure filter (``parse_boundary_selector`` / ``scope_matches`` / ``load_history``)
  * the real HTTP path — ``handle_tadvci_signout`` writes the JSONL, then
    ``handle_tadvci_history`` reads it back through the same route the UI calls, and the
    record's own checksum is re-verified from what actually landed on disk.

Nothing here writes to the shipped ``backend/tad_vci/decision_records_v3.jsonl``; one
test asserts that explicitly.
"""
from __future__ import annotations

import asyncio
import json
import sys
from pathlib import Path

import pandas as pd
import pytest
from aiohttp.test_utils import make_mocked_request

BACKEND = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BACKEND))

from tad_vci import api_routes  # noqa: E402
from tad_vci.curation import (  # noqa: E402
    append_snapshot, evidence_card, load_history, make_snapshot, next_version,
    parse_boundary_selector, record_scope, scope_matches, verify_snapshot,
)

REPO_LOG = BACKEND / "tad_vci" / "decision_records_v3.jsonl"

BID = "3:1234567"
ROW = {
    "chrom": "3", "pos": 1_234_567, "D_tier": "D4",
    "BD1_caller_support": "strong", "BD2_ctcf": "moderate",
    "BD3_rad21": "weak", "BD4_loop_anchor": "not_assessable",
    "votes": 4, "n_methods": 5,
}
TS = "2026-07-26T09:00:00Z"


def _write(log: Path, cell: str, assembly: str, resolution: int, actor: str) -> dict:
    card = evidence_card(dict(ROW))
    snap = make_snapshot(
        card, assembly, actor, TS,
        decision={"final_tier": "D4", "verdict": "supported_for_review",
                  "rationale": f"{cell} {resolution}"},
        prior_version=next_version(card["boundary_id"], log),
        resolution=resolution, cell_line=cell,
    )
    append_snapshot(snap, log)
    return json.loads(snap.to_json())


# ── pure filter ───────────────────────────────────────────────────────────────

def test_selector_parsing_round_trip():
    assert parse_boundary_selector("3:1234567") == ("3:1234567", None)
    bid, scope = parse_boundary_selector("3:1234567@GM12878/GRCh38/10000")
    assert bid == "3:1234567"
    # cell upper-cased, GRCh38 aliased to hg38, resolution an int — so the UI may send
    # whatever casing/alias the manifest happens to carry.
    assert scope == {"cell_line": "GM12878", "assembly": "hg38", "resolution": 10000}


def test_malformed_suffix_is_not_silently_widened():
    """A half-written selector must not degrade into 'show me every cell line'."""
    assert parse_boundary_selector("3:1234567@GM12878") == ("3:1234567@GM12878", None)
    assert parse_boundary_selector("3:1234567@GM12878/hg38") == ("3:1234567@GM12878/hg38", None)


def test_scope_matches_is_strict_on_every_field():
    rec = {"cell_line": "K562", "assembly": "hg38", "resolution": 10000}
    assert record_scope(rec) == {"cell_line": "K562", "assembly": "hg38", "resolution": 10000}
    assert scope_matches(rec, None) is True
    assert scope_matches(rec, record_scope(rec)) is True
    assert scope_matches(rec, {**record_scope(rec), "cell_line": "GM12878"}) is False
    assert scope_matches(rec, {**record_scope(rec), "assembly": "hg19"}) is False
    assert scope_matches(rec, {**record_scope(rec), "resolution": 25000}) is False


def test_legacy_record_without_a_cell_line_is_not_attributed_to_one():
    """A pre-v3 line carries no cell_line. Showing it under a cell-scoped query would be
    a guess about which map it was reviewed on, so it is withheld there and stays
    reachable only through an unscoped query / the all-records view."""
    legacy = {"boundary_id": BID, "assembly": "hg19", "version": 1}
    assert scope_matches(legacy, None) is True
    assert scope_matches(legacy, {"cell_line": "GM12878", "assembly": "hg19",
                                  "resolution": 25000}) is False


# ── log-level behaviour ───────────────────────────────────────────────────────

@pytest.fixture
def log_path(tmp_path: Path) -> Path:
    return tmp_path / "decision_records_v3.jsonl"


def test_history_scope_separates_cell_lines_at_the_same_coordinate(log_path):
    gm = _write(log_path, "GM12878", "hg38", 10_000, "reviewer_gm")
    k5 = _write(log_path, "K562", "hg38", 10_000, "reviewer_k5")
    assert gm["boundary_id"] == k5["boundary_id"] == BID

    gm_hist = load_history(f"{BID}@GM12878/hg38/10000", log_path)
    k5_hist = load_history(f"{BID}@K562/hg38/10000", log_path)
    assert [r["actor_label"] for r in gm_hist] == ["reviewer_gm"]
    assert [r["actor_label"] for r in k5_hist] == ["reviewer_k5"]
    # and the un-scoped query still sees the whole log for this boundary
    assert len(load_history(BID, log_path)) == 2


def test_history_scope_separates_resolution_and_assembly(log_path):
    _write(log_path, "GM12878", "hg38", 10_000, "at_10kb")
    _write(log_path, "GM12878", "hg38", 25_000, "at_25kb")
    _write(log_path, "GM12878", "hg19", 10_000, "on_hg19")
    assert [r["actor_label"] for r in load_history(f"{BID}@GM12878/hg38/10000", log_path)] \
        == ["at_10kb"]
    assert [r["actor_label"] for r in load_history(f"{BID}@GM12878/hg38/25000", log_path)] \
        == ["at_25kb"]
    assert [r["actor_label"] for r in load_history(f"{BID}@GM12878/hg19/10000", log_path)] \
        == ["on_hg19"]
    assert len(load_history(BID, log_path)) == 3


def test_versions_are_unique_across_scopes(log_path):
    """next_version stays deliberately unscoped: two records may not both claim v1 for
    one boundary_id, or an export cannot tell them apart. Per-scope history therefore
    starts wherever the log had got to."""
    _write(log_path, "GM12878", "hg38", 10_000, "a")
    _write(log_path, "K562", "hg38", 10_000, "b")
    versions = [r["version"] for r in load_history(BID, log_path)]
    assert versions == [1, 2]
    assert load_history(f"{BID}@K562/hg38/10000", log_path)[0]["version"] == 2
    # a scoped selector must not reset the counter either
    assert next_version(f"{BID}@K562/hg38/10000", log_path) == 2


def test_every_written_record_verifies_from_disk(log_path):
    _write(log_path, "GM12878", "hg38", 10_000, "a")
    _write(log_path, "K562", "hg38", 25_000, "b")
    lines = [json.loads(x) for x in log_path.read_text(encoding="utf-8").splitlines() if x.strip()]
    assert len(lines) == 2
    for rec in lines:
        assert verify_snapshot(rec) is True
        # the scope fields are inside the checksum payload, so they cannot be edited
        # afterwards to move a record between cell lines without detection
        assert verify_snapshot({**rec, "cell_line": "HeLa"}) is False
        assert verify_snapshot({**rec, "resolution": 5000}) is False


# ── the real HTTP path the UI uses ────────────────────────────────────────────

def _signout(body: dict, log_path: Path) -> tuple[int, dict]:
    req = make_mocked_request("POST", "/api/tadvci/signout")
    async def _json():
        return body
    req.json = _json
    orig = api_routes._LOG
    api_routes._LOG = log_path
    try:
        resp = asyncio.run(api_routes.handle_tadvci_signout(req))
    finally:
        api_routes._LOG = orig
    return resp.status, json.loads(resp.body.decode())


def _history(selector: str, log_path: Path) -> dict:
    req = make_mocked_request("GET", f"/api/tadvci/history/{selector}",
                              match_info={"boundary_id": selector})
    orig = api_routes._LOG
    api_routes._LOG = log_path
    try:
        resp = asyncio.run(api_routes.handle_tadvci_history(req))
    finally:
        api_routes._LOG = orig
    return json.loads(resp.body.decode())


@pytest.fixture
def imported_token(tmp_path: Path) -> str:
    """An imported annotation that owns BID, so the card lookup is deterministic and does
    not depend on the built-in GM12878 bundle being present on this host."""
    df = pd.DataFrame({
        "chrom": ["3"], "pos": [1_234_567], "votes": [4], "n_methods": [5],
        "BD1_caller_support": ["strong"], "BD2_ctcf": ["moderate"],
        "BD3_rad21": ["weak"], "BD4_loop_anchor": ["not_assessable"],
        "D_tier": ["D4"], "assembly": ["hg38"],
    })
    tsv = tmp_path / "scope_upload.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    api_routes.set_imported_annotation("tk_scope", str(tsv))
    return "tk_scope"


def _ui_body(cell: str, resolution: int, token: str, actor: str) -> dict:
    """Exactly the payload TADVCIEvidence.vue posts from the sign-out form."""
    return {
        "boundary_id": BID, "final_tier": "D4", "verdict": "supported_for_review",
        "rationale": f"recorded from the UI for {cell}", "actor_label": actor,
        "resolution": resolution, "cell_line": cell, "assembly": "hg38", "token": token,
    }


def test_ui_signout_lands_in_the_jsonl_with_a_verifying_checksum(log_path, imported_token):
    status, payload = _signout(_ui_body("K562", 10_000, imported_token, "ui_reviewer"), log_path)
    assert status == 200, payload
    assert payload["checksum_matches"] is True

    on_disk = [json.loads(x) for x in log_path.read_text(encoding="utf-8").splitlines() if x.strip()]
    assert len(on_disk) == 1
    rec = on_disk[0]
    assert rec["boundary_id"] == BID
    assert rec["cell_line"] == "K562"
    assert rec["resolution"] == 10_000
    assert rec["assembly"] == "hg38"
    assert rec["actor_label"] == "ui_reviewer"
    assert rec["review_decision"]["rationale"] == "recorded from the UI for K562"
    # the checksum is recomputed from the bytes that actually landed, not from the
    # in-memory snapshot the handler returned
    assert verify_snapshot(rec) is True
    assert rec["content_hash"] == payload["snapshot"]["content_hash"]


def test_history_route_honours_the_scope_suffix(log_path, imported_token):
    assert _signout(_ui_body("GM12878", 10_000, imported_token, "gm_reviewer"), log_path)[0] == 200
    assert _signout(_ui_body("K562", 10_000, imported_token, "k5_reviewer"), log_path)[0] == 200

    gm = _history(f"{BID}@GM12878/hg38/10000", log_path)
    k5 = _history(f"{BID}@K562/hg38/10000", log_path)
    assert [r["actor_label"] for r in gm["history"]] == ["gm_reviewer"]
    assert [r["actor_label"] for r in k5["history"]] == ["k5_reviewer"]
    # the bare id is still accepted and still returns the full log for that coordinate
    assert len(_history(BID, log_path)["history"]) == 2
    # a scope with no records is empty, not a fallback to "everything"
    assert _history(f"{BID}@HeLa/hg38/10000", log_path)["history"] == []


def test_tests_never_touch_the_shipped_decision_log(log_path, imported_token):
    before = REPO_LOG.read_text(encoding="utf-8") if REPO_LOG.exists() else None
    _signout(_ui_body("GM12878", 10_000, imported_token, "gm_reviewer"), log_path)
    after = REPO_LOG.read_text(encoding="utf-8") if REPO_LOG.exists() else None
    assert after == before
