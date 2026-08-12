"""The tier track is served slim; the full card is served one boundary at a time.

Serving a whole evidence card per boundary made one chromosome of a 10 kb package a
21.2 MB / 22.5 s download and held ~250 MB of card dicts per token, to paint ticks that
need five fields. These tests pin the split so it cannot regress:

* the track payload carries exactly the fields the track draws with — and nothing that
  only the panel reads (that is what made it big);
* the per-boundary endpoint still returns the complete card, including the measurement
  blocks added on 2026-07-29;
* `view=full` still reproduces the old shape for anything that needs it;
* the two agree on tier and votes for every boundary, so slimming cannot desynchronise
  the tick from the card behind it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "backend"))
from tad_vci import api_routes  # noqa: E402

PKG = REPO / "example_data" / "multicell" / "IMR90_25kb" / "annotation.tsv"
TOKEN = "test-token"


@pytest.fixture()
def token(tmp_path):
    if not PKG.exists():
        pytest.skip(f"demo package missing: {PKG}")
    # a small slice keeps the test fast but exercises the real table's columns
    df = pd.read_csv(PKG, sep="\t", low_memory=False)
    df = df[df["chrom"].astype(str) == "1"].head(200)
    out = tmp_path / "annotation.tsv"
    df.to_csv(out, sep="\t", index=False)
    api_routes.set_imported_annotation(TOKEN, str(out))
    yield TOKEN
    api_routes.IMPORTED_ANNOT.pop(TOKEN, None)
    api_routes._imported_cards_cache.pop(TOKEN, None)
    api_routes._imported_track_cache.pop(TOKEN, None)
    api_routes._imported_df_cache.pop(TOKEN, None)


def test_track_payload_carries_only_what_the_track_draws(token):
    rows = api_routes._track_for_token(token, "1")
    assert rows, "the slice has boundaries on chr1"
    assert set(rows[0]) == set(api_routes._TRACK_FIELDS)
    # the fields the panel needs must NOT be in the list payload
    for heavy in ("criteria", "summary", "tier_why", "bd1_threshold", "supporters"):
        assert heavy not in rows[0]
    # a track row is a fraction of a card
    import json
    card = api_routes._card_for_position(token, "1", rows[0]["pos"])
    assert len(json.dumps(rows[0])) * 8 < len(json.dumps(card))


def test_card_endpoint_returns_the_complete_card(token):
    rows = api_routes._track_for_token(token, "1")
    row = rows[0]
    card = api_routes._card_for_position(token, "1", row["pos"])
    assert card is not None
    assert card["boundary_id"] == row["boundary_id"]
    assert card["auto_tier"] == row["auto_tier"]
    assert card["votes"] == row["votes"]
    assert set(card["criteria"]) == {"BD1", "BD2", "BD3", "BD4"}
    # the 2026-07-29 disclosures must survive the split
    assert card["criteria"]["BD1"]["measured"]["vote_tol_bp"]
    assert card["criteria"]["BD3"]["measured"]["cutpoints"]
    assert card["bd1_match_rule"]


def test_missing_position_is_a_miss_not_a_wrong_card(token):
    assert api_routes._card_for_position(token, "1", 12_345) is None
    assert api_routes._card_for_position(token, "22", 1_000_000) is None


def test_track_and_full_view_agree_boundary_for_boundary(token):
    slim = api_routes._track_for_token(token, "1")
    full = api_routes._cards_for_token(token, "1")
    assert len(slim) == len(full)
    for s, f in zip(slim, full):
        assert s["boundary_id"] == f["boundary_id"]
        assert s["auto_tier"] == f["auto_tier"]
        assert s["pos"] == f["pos"]
        assert s["votes"] == f["votes"]


def test_repointing_a_token_invalidates_every_cache(token, tmp_path):
    first = api_routes._track_for_token(token, "1")
    assert first
    empty = tmp_path / "other.tsv"
    pd.DataFrame({"chrom": ["1"], "pos": [999_000_000], "votes": [1],
                  "n_methods": [6], "n_methods_effective": [6],
                  "BD1_caller_support": ["weak"], "BD2_ctcf": ["not_assessable"],
                  "BD3_rad21": ["not_assessable"], "BD4_loop_anchor": ["not_assessable"],
                  "D_tier": ["D1"], "user_bed": [0]}).to_csv(empty, sep="\t", index=False)
    api_routes.set_imported_annotation(token, str(empty))
    again = api_routes._track_for_token(token, "1")
    assert [r["pos"] for r in again] == [999_000_000]
