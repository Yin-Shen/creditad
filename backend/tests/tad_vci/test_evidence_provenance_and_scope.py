"""Provenance and scope guarantees on the TAD-VCI routes.

The tool's deliverable is a signed local review record, so the two things that may never
go wrong are (a) which dataset a record's evidence came from and (b) which genome build a
contact map is drawn on. Each test below pins a case that used to fail silently — the
request succeeded and the picture looked right.
"""
import asyncio
import json
import os
import sys
from pathlib import Path

import pandas as pd
import pytest
from aiohttp.test_utils import make_mocked_request

BACKEND = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BACKEND))

from tad_vci import api_routes  # noqa: E402
from tad_vci.curation import evidence_card  # noqa: E402


def _json(resp) -> dict:
    return json.loads(resp.body.decode())


def _annotation(path: Path, tier: str, support: str) -> Path:
    """One boundary at 7:875000, graded differently in each file."""
    pd.DataFrame({
        "chrom": ["7"], "pos": [875_000], "votes": [5 if support == "strong" else 2],
        "n_methods": [6], "BD1_caller_support": [support],
        "BD2_ctcf": ["not_assessable"], "BD3_rad21": ["not_assessable"],
        "BD4_loop_anchor": ["none"], "D_tier": [tier], "assembly": ["hg38"],
    }).to_csv(path, sep="\t", index=False)
    return path


def _signout(body: dict, log_path: Path):
    req = make_mocked_request("POST", "/api/tadvci/signout")

    async def _j():
        return body
    req.json = _j
    original = api_routes._LOG
    api_routes._LOG = log_path
    try:
        resp = asyncio.run(api_routes.handle_tadvci_signout(req))
    finally:
        api_routes._LOG = original
    return resp.status, _json(resp)


@pytest.fixture
def two_datasets(tmp_path):
    """Two datasets that disagree about the SAME boundary_id.

    boundary_id is only "chrom:pos" and the shipped demo packages share thousands of
    positions at different tiers, so a lookup that is not pinned to one dataset can
    produce a record labelled with one cell line carrying another's evidence.
    """
    a = _annotation(tmp_path / "a.tsv", "D5", "strong")
    b = _annotation(tmp_path / "b.tsv", "D3", "weak")
    api_routes.set_imported_annotation("tok-a", str(a))
    api_routes.set_imported_annotation("tok-b", str(b))
    yield
    for tok in ("tok-a", "tok-b"):
        api_routes.IMPORTED_ANNOT.pop(tok, None)
        api_routes._imported_cards_cache.pop(tok, None)
        api_routes._imported_df_cache.pop(tok, None)


_BASE = {"boundary_id": "7:875000", "final_tier": "D3", "rationale": "review note",
         "verdict": "uncertain", "actor_label": "test_actor", "resolution": 25_000}


@pytest.mark.parametrize("token,expected_tier,expected_support", [
    ("tok-a", "D5", "strong"),
    ("tok-b", "D3", "weak"),
])
def test_signout_records_the_evidence_of_the_named_dataset(
        two_datasets, tmp_path, token, expected_tier, expected_support):
    status, payload = _signout({**_BASE, "token": token, "cell_line": "IMR90"},
                               tmp_path / "log.jsonl")
    assert status == 200, payload
    snap = payload["snapshot"]
    assert snap["auto_tier"] == expected_tier
    assert snap["evidence"]["BD1"]["grade"] == expected_support
    assert payload["checksum_matches"] is True


def test_signout_refuses_to_guess_which_dataset_a_decision_is_about(two_datasets, tmp_path):
    """Without a token the handler used to scan every registered import and take the
    first hit, producing a record whose label and evidence came from different datasets
    while its own checksum still verified."""
    status, payload = _signout({**_BASE, "cell_line": "IMR90"}, tmp_path / "log.jsonl")
    assert status == 422
    assert "token required" in payload["error"]


def test_signout_does_not_fall_back_to_the_default_cell_line(two_datasets, tmp_path):
    """An unknown cell_line was silently downgraded to GM12878 and its built-in cards
    were searched BEFORE the caller's own token."""
    status, payload = _signout({**_BASE, "cell_line": "K562"}, tmp_path / "log.jsonl")
    assert status == 422


def test_signout_rejects_an_unknown_token(two_datasets, tmp_path):
    status, payload = _signout({**_BASE, "token": "tok-gone"}, tmp_path / "log.jsonl")
    assert status == 404
    assert "unknown token" in payload["error"]


def test_signout_reports_which_source_a_missing_boundary_was_looked_for_in(
        two_datasets, tmp_path):
    status, payload = _signout({**_BASE, "boundary_id": "7:1", "token": "tok-a"},
                               tmp_path / "log.jsonl")
    assert status == 404
    assert "tok-a" in payload["error"]


# ── genome-build gate (was: "could not verify" passed as "matches") ──


def _profile_gate(profile: dict, pkg_assembly: str) -> bool:
    """True when the map is served. Mirrors handle_tadvci_use_multicell_example."""
    if profile["assembly"] is None:
        return False
    return not (pkg_assembly and profile["assembly"] != pkg_assembly)


def test_an_unverifiable_matrix_is_not_served(tmp_path):
    """The signature table only knows chr1/chr7/chrX, so a subset matrix profiles as
    assembly=None. That must block the map, not pass the build check."""
    cooler = pytest.importorskip("cooler")
    import numpy as np

    path = tmp_path / "hg19_chr2_only.cool"
    res, length = 1_000_000, 243_199_373          # hg19 chr2
    n = int(np.ceil(length / res))
    bins = pd.DataFrame({"chrom": ["chr2"] * n, "start": np.arange(n) * res,
                         "end": np.minimum((np.arange(n) + 1) * res, length)})
    pixels = pd.DataFrame({"bin1_id": np.arange(n - 1), "bin2_id": np.arange(1, n),
                           "count": np.ones(n - 1, dtype=int)})
    cooler.create_cooler(str(path), bins, pixels, ordered=True)

    profile = api_routes._matrix_profile(path)
    assert profile["assembly"] is None
    assert profile["assembly_source"] is None
    assert profile["unverified_reason"], "an unverified matrix must say why"
    assert _profile_gate(profile, "hg38") is False


def test_a_matrix_that_declares_its_build_is_served(tmp_path):
    """Fail-closed must not brick a legitimate subset matrix: cooler's own metadata is
    accepted as a second, weaker source once no length identifies the build."""
    cooler = pytest.importorskip("cooler")
    import numpy as np

    path = tmp_path / "hg38_chr2_declared.cool"
    res, length = 1_000_000, 242_193_529          # hg38 chr2
    n = int(np.ceil(length / res))
    bins = pd.DataFrame({"chrom": ["chr2"] * n, "start": np.arange(n) * res,
                         "end": np.minimum((np.arange(n) + 1) * res, length)})
    pixels = pd.DataFrame({"bin1_id": np.arange(n - 1), "bin2_id": np.arange(1, n),
                           "count": np.ones(n - 1, dtype=int)})
    cooler.create_cooler(str(path), bins, pixels, ordered=True, assembly="hg38")

    profile = api_routes._matrix_profile(path)
    assert profile["assembly"] == "hg38"
    assert profile["assembly_source"] == "file_metadata"
    assert _profile_gate(profile, "hg38") is True
    assert _profile_gate(profile, "hg19") is False


def test_a_missing_matrix_profiles_as_unverified():
    profile = api_routes._matrix_profile(None)
    assert profile["assembly"] is None and profile["unverified_reason"]


# ── request validation ──


@pytest.mark.parametrize("raw,expected", [(None, 25_000), ("", 25_000),
                                          (10_000, 10_000), (10_000.0, 10_000),
                                          ("10000", 10_000)])
def test_resolution_accepts_an_absent_or_whole_value(raw, expected):
    assert api_routes._parse_resolution(raw, 25_000) == expected


@pytest.mark.parametrize("raw", [0, -5, 10_000.7, "abc", True])
def test_resolution_rejects_every_value_it_cannot_honour(raw):
    """0 used to be swallowed by `int(raw or default)` and silently answered as 25 kb,
    while -5 was rejected — the same field failed loudly or quietly by accident.
    Resolution is hashed into every sign-out record, so it may not be guessed."""
    with pytest.raises(ValueError):
        api_routes._parse_resolution(raw, 25_000)


@pytest.mark.parametrize("tag", ["annotate", "run-1", "a.b_c", "E2E"])
def test_hash_tag_accepts_ordinary_labels(tag):
    assert api_routes._SAFE_TAG_RE.fullmatch(tag)


@pytest.mark.parametrize("tag", ["../../evil", "a/b", "a\\b", "", "x" * 65, "a b"])
def test_hash_tag_rejects_anything_that_can_escape_the_output_directory(tag):
    """hash_tag is concatenated into output FILE NAMES, so a separator in it walks out
    of the chosen directory; unvalidated it produced a 500 that echoed the resolved
    server path back to the caller."""
    assert api_routes._SAFE_TAG_RE.fullmatch(tag) is None


def test_output_dir_is_normalised_and_created(tmp_path):
    target = tmp_path / "nested" / "out"
    resolved = api_routes._resolved_output_dir(str(target))
    assert Path(resolved).is_dir()
    assert resolved == os.path.abspath(str(target))


def test_output_dir_of_a_callers_choosing_is_honoured(tmp_path):
    """Not confined to a fixed root: the GUI, the CLI and the parity check all name
    their own results directory."""
    assert api_routes._resolved_output_dir(str(tmp_path)) == os.path.abspath(str(tmp_path))


def test_output_dir_that_cannot_be_created_is_a_client_error(tmp_path):
    blocker = tmp_path / "afile"
    blocker.write_text("not a directory")
    with pytest.raises(ValueError):
        api_routes._resolved_output_dir(str(blocker / "under"))


def test_empty_output_dir_is_rejected():
    with pytest.raises(ValueError):
        api_routes._resolved_output_dir("   ")


# ── review-record location (was: inside the install tree, wiped by every upgrade) ──


def test_records_log_follows_the_configured_directory(tmp_path, monkeypatch):
    monkeypatch.setenv("CREDITAD_RECORDS_DIR", str(tmp_path / "userdata"))
    log = api_routes._records_log()
    assert log.parent == tmp_path / "userdata"
    assert log.parent.is_dir()
    assert log.name == "decision_records_v3.jsonl", "an existing log must keep working"


def test_records_log_defaults_beside_the_module(monkeypatch):
    monkeypatch.delenv("CREDITAD_RECORDS_DIR", raising=False)
    log = api_routes._records_log()
    assert log.parent == Path(api_routes.__file__).resolve().parent
    assert log.name == "decision_records_v3.jsonl"


def test_records_log_falls_back_when_the_configured_directory_is_unusable(monkeypatch):
    monkeypatch.setenv("CREDITAD_RECORDS_DIR", "/proc/definitely/not/writable")
    assert api_routes._records_log().name == "decision_records_v3.jsonl"


# ── optional full-genome search order (was: a developer path baked in unconditionally) ──


def test_the_development_data_tree_is_not_searched_off_posix(monkeypatch):
    """A packaged Windows build resolved host files through wine's Z: mapping during the
    smoke test, so the branch a real user takes was never exercised."""
    import data_manager

    monkeypatch.setattr(data_manager.os, "name", "nt")
    assert data_manager.dev_shared_raw_default() is None
    for root in api_routes._shared_raw_roots():
        assert "_shared_data" not in str(root)


def test_the_development_data_tree_is_skipped_when_it_does_not_exist(monkeypatch):
    import data_manager

    monkeypatch.setattr(data_manager.os, "name", "posix")
    monkeypatch.setenv("CREDITAD_DEV_SHARED_RAW", "/nonexistent/dev/tree")
    assert data_manager.dev_shared_raw_default() is None


def test_an_explicit_shared_raw_root_is_searched_first(monkeypatch, tmp_path):
    monkeypatch.setenv("CREDITAD_SHARED_RAW", str(tmp_path))
    assert api_routes._shared_raw_roots()[0] == tmp_path


# ── per-cell threshold disclosure (was: an uncalibrated card claimed per-cell numbers) ──


_CARD_ROW = {"chrom": "7", "pos": 1_000_000, "votes": 3, "n_methods": 6,
             "BD1_caller_support": "moderate", "BD2_ctcf": "moderate",
             "BD3_rad21": "weak", "BD4_loop_anchor": "none", "D_tier": "D3"}


def test_a_calibrated_cell_line_reports_its_own_thresholds():
    card = evidence_card({**_CARD_ROW, "data_cell_line": "K562"})
    assert "per-cell K562" in card["criteria"]["BD2"]["grading"]
    assert card["bands_disclosure"] is None


@pytest.mark.parametrize("cell", ["HeLa", None])
def test_an_uncalibrated_cell_line_says_the_thresholds_are_not_its_own(cell):
    """load_bands returns the GM12878 placeholder with calibrated=False for an unknown
    cell line, but the card still read "per-cell GM12878" and, with no bands_calibrated
    column to trigger the disclosure, nothing on the card said the numbers were borrowed."""
    card = evidence_card({**_CARD_ROW, "data_cell_line": cell})
    for code in ("BD2", "BD3"):
        assert "NOT calibrated" in card["criteria"][code]["grading"]
        assert "per-cell" not in card["criteria"][code]["grading"]
    assert card["bands_disclosure"] and "NOT calibrated" in card["bands_disclosure"]


def test_a_band_table_failure_is_disclosed_rather_than_swallowed(monkeypatch):
    """The rewrite was wrapped in `except Exception: pass`, so a failure left the static
    GM12878 registry numbers on the card with nothing saying they were a fallback."""
    import tad_vci.graded_evidence as ge

    def _boom(_cell):
        raise RuntimeError("band table unreadable")

    monkeypatch.setattr(ge, "load_bands", _boom)
    card = evidence_card({**_CARD_ROW, "data_cell_line": "GM12878"})
    assert card["bands_disclosure"] and "could not be resolved" in card["bands_disclosure"]


# ── one session per run (was: two runs of the same BED shared a session token) ──


class _Manager:
    """GlobalManager's register/get contract, de-duplicating by path as it does."""

    def __init__(self):
        self.by_path = {}

    def register(self, path, file_type, prediction_path=None):
        return self.by_path.setdefault(path, f"tok-{len(self.by_path)}")

    def get(self, token):
        return object() if token in self.by_path.values() else None


def _annotate(body: dict, manager):
    req = make_mocked_request("POST", "/api/tadvci/annotate", app={"manager": manager})

    async def _j():
        return body
    req.json = _j
    resp = asyncio.run(api_routes.handle_tadvci_annotate(req))
    return resp.status, _json(resp)


def _boundaries(path: Path, positions, resolution=25_000) -> Path:
    path.write_text("".join(f"7\t{p}\t{p + resolution}\n" for p in positions))
    return path


def test_two_runs_of_one_candidate_bed_get_separate_sessions(tmp_path):
    """register() de-duplicates by path, so binding the session to the primary BED handed
    the second run the first run's token and silently re-pointed that session's tier
    table — an already-open tab then showed another run's tiers with no indication."""
    primary = _boundaries(tmp_path / "primary.bed", [1_000_000, 2_000_000, 3_000_000])
    panel = _boundaries(tmp_path / "panel.bed", [1_000_000])
    manager = _Manager()

    status_a, run_a = _annotate({"primary_bed": str(primary), "hash_tag": "runA",
                                 "output_dir": str(tmp_path / "outA")}, manager)
    assert status_a == 200, run_a
    status_b, run_b = _annotate({"primary_bed": str(primary), "hash_tag": "runB",
                                 "method_beds": [{"path": str(panel), "name": "panel"}],
                                 "output_dir": str(tmp_path / "outB")}, manager)
    assert status_b == 200, run_b

    assert run_a["session_token"] != run_b["session_token"]
    assert api_routes.IMPORTED_ANNOT[run_a["session_token"]] == run_a["TADVCI_annotation"]
    assert api_routes.IMPORTED_ANNOT[run_b["session_token"]] == run_b["TADVCI_annotation"]
    for token in (run_a["session_token"], run_b["session_token"]):
        api_routes.IMPORTED_ANNOT.pop(token, None)
        api_routes._imported_cards_cache.pop(token, None)


def test_annotate_rejects_a_hash_tag_that_escapes_the_output_directory(tmp_path):
    primary = _boundaries(tmp_path / "primary.bed", [1_000_000])
    outside = tmp_path / "should_not_exist"
    status, payload = _annotate({"primary_bed": str(primary), "hash_tag": "../../evil",
                                 "output_dir": str(outside)}, _Manager())
    assert status == 400
    assert "hash_tag" in payload["error"]
    assert not outside.exists(), "a rejected request must not create anything"


@pytest.mark.parametrize("handler,path", [
    ("handle_tadvci_annotate", "/api/tadvci/annotate"),
    ("handle_tadvci_signout", "/api/tadvci/signout"),
    ("handle_tadvci_use_multicell_example", "/api/tadvci/use_multicell_example"),
])
@pytest.mark.parametrize("body", [[1, 2, 3], "a string"])
def test_a_non_object_body_is_a_client_error(handler, path, body):
    """These reached body.get() on a list and raised AttributeError, which the middleware
    turned into an opaque 500 for what is plainly a bad request."""
    req = make_mocked_request("POST", path)

    async def _j():
        return body
    req.json = _j
    resp = asyncio.run(getattr(api_routes, handler)(req))
    assert resp.status == 400
    assert "JSON object" in _json(resp)["error"]
