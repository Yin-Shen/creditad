"""Input guards and resource bounds on the raw-data routes in backend/main.py.

Every case here corresponds to a request that used to be answered with an opaque 500,
a silently wrong window, or unbounded work on the event loop. The handlers are aiohttp
coroutines and the repo ships no pytest-asyncio plugin, so each is driven directly with
asyncio.run() and a mocked request, as backend/tests/tad_vci/test_api_triage_signout.py
already does.
"""
import asyncio
import json
import os
import sys
from pathlib import Path

import pytest
from aiohttp.test_utils import make_mocked_request

# `main` is imported from the installed package (pip install '.[app]'); the source tree
# no longer needs to be prepended to sys.path.

# main.py chdir()s into the backend directory at import time so the embedded Python
# resolves its relative data paths. That is a process-wide side effect, and sibling test
# modules address fixtures relative to the repository root, so put the cwd back.
_CWD = os.getcwd()
import main  # noqa: E402
os.chdir(_CWD)


def _json(resp) -> dict:
    return json.loads(resp.body.decode())


# ── window validation (was: silent clamp, negative start read another chromosome) ──


def test_bp_window_accepts_an_in_range_request():
    assert main._bp_window(0, 1_000, 10_000) == (0, 1_000)


def test_bp_window_truncates_an_overlong_end_to_the_sequence():
    assert main._bp_window(500, 99_999, 10_000) == (500, 10_000)


@pytest.mark.parametrize("start,end", [(-1_000_000, -1), (-1, 100), (0, -5)])
def test_bp_window_rejects_negative_coordinates(start, end):
    """Negative bp used to be clamped to bin 0 and then re-offset by the chromosome's
    start bin, so a request for chr7 was answered with contacts from the end of chr6."""
    with pytest.raises(ValueError, match=">= 0"):
        main._bp_window(start, end, 10_000)


@pytest.mark.parametrize("start,end", [(100, 100), (200, 100)])
def test_bp_window_rejects_an_empty_or_reversed_window(start, end):
    with pytest.raises(ValueError, match="greater than start"):
        main._bp_window(start, end, 10_000)


def test_bp_window_rejects_a_start_past_the_end_of_the_sequence():
    with pytest.raises(ValueError, match="past the end"):
        main._bp_window(20_000, 30_000, 10_000)


# ── range size cap (was: an unclamped span read the whole genome and OOM-killed us) ──


class _StubMatrix:
    """Minimal stand-in for HiCRangeRenderer's data_map contract."""

    path = "<stub>"
    resolutions = [25_000]

    def __init__(self):
        self.data_map = {25_000: {"chrom_info": {"7": {"length": 159_345_973,
                                                       "start_bin": 0, "n_bins": 6374}}}}
        self.sync_calls = 0

    fetch_range_data = main.HiCRangeRenderer.fetch_range_data

    def _fetch_range_sync(self, res, b1, b2):
        self.sync_calls += 1
        return [], res


def test_fetch_range_data_refuses_a_span_wider_than_the_cap():
    stub = _StubMatrix()
    over = main._MAX_RANGE_BINS + 1
    with pytest.raises(ValueError, match="range too large"):
        asyncio.run(stub.fetch_range_data(25_000, (0, over), (0, 10)))
    assert stub.sync_calls == 0, "the oversized read must be refused before any IO"


def test_fetch_range_data_allows_a_span_at_the_cap():
    stub = _StubMatrix()
    asyncio.run(stub.fetch_range_data(25_000, (0, main._MAX_RANGE_BINS), (0, 10)))
    assert stub.sync_calls == 1


def test_tile_and_colour_probe_spans_stay_under_the_cap():
    """The two internal callers must never be able to trip the guard: tiles are 256 bins
    and the colour-scale probe samples 3 x 256."""
    assert 256 <= main._MAX_RANGE_BINS
    assert 3 * 256 <= main._MAX_RANGE_BINS


# ── bigwig signal guards (was: HTTP 200 for a negative window, unbounded bins) ──


class _StubBigWig(main.BigWigRenderer):
    def __init__(self):
        self.path = "<stub>"
        self.chroms = {"chr7": 159_345_973}
        self.name_map = {"7": "chr7"}
        self.fetched = []

    async def fetch_signal(self, chrom, start, end, bins=500):
        self.fetched.append((chrom, start, end, bins))
        return []


def _signal(bw, query: str, chrom: str = "7"):
    req = make_mocked_request(
        "GET", f"/api/bigwig/signal/tok/{chrom}?{query}",
        match_info={"token": "tok", "chrom": chrom},
        app={"manager": type("M", (), {"get": staticmethod(lambda _t: bw)})()},
    )
    return asyncio.run(main.handle_bigwig_signal(req))


def test_bigwig_signal_serves_a_reasonable_request():
    bw = _StubBigWig()
    resp = _signal(bw, "start=0&end=10000000&bins=500")
    assert resp.status == 200
    assert bw.fetched == [("7", 0, 10_000_000, 500)]


@pytest.mark.parametrize("query,message", [
    ("start=-100&end=-1&bins=10", ">= 0"),
    ("start=500&end=100&bins=10", "greater than start"),
    ("start=0&end=1000&bins=0", "between 1 and"),
    ("start=0&end=1000&bins=2000000", "between 1 and"),
])
def test_bigwig_signal_rejects_impossible_windows_and_bin_counts(query, message):
    bw = _StubBigWig()
    resp = _signal(bw, query)
    assert resp.status == 400
    assert message in _json(resp)["error"]
    assert bw.fetched == [], "a rejected request must not start any work"


def test_bigwig_signal_allows_the_maximum_bin_count():
    bw = _StubBigWig()
    assert _signal(bw, f"start=0&end=159000000&bins={main._MAX_SIGNAL_BINS}").status == 200


def test_bigwig_signal_reports_an_unknown_chromosome_as_404():
    bw = _StubBigWig()
    resp = _signal(bw, "start=0&end=1000", chrom="ZZZ")
    assert resp.status == 404
    assert "chromosome not found" in _json(resp)["error"]


# ── malformed bodies (was: AttributeError -> opaque 500) ──


def _post(handler, body, path):
    req = make_mocked_request("POST", path)

    async def _j():
        return body
    req.json = _j
    req.read = _j
    return asyncio.run(handler(req))


@pytest.mark.parametrize("body", [[1, 2, 3], "a string", 7])
def test_tad_run_rejects_a_non_object_body_as_400(body):
    resp = _post(main.handle_tad_run, body, "/api/tad/run")
    assert resp.status == 400
    assert "JSON object" in _json(resp)["error"]


def test_register_rejects_a_non_object_body_as_400():
    req = make_mocked_request("POST", "/api/register")

    async def _read():
        return b"[1,2,3]"
    req.read = _read
    resp = asyncio.run(main.handle_register(req))
    assert resp.status == 400
    assert "JSON object" in _json(resp)["error"]


# ── BED chromosome lookup (was: HTTP 200 carrying a Chinese error string) ──


def _bed_renderer(tmp_path):
    bed = tmp_path / "b.bed"
    bed.write_text("7\t1000000\t1025000\n7\t2000000\t2025000\n")
    return main.BedTileRenderer(str(bed))


def _bed_info(renderer, chrom):
    req = make_mocked_request(
        "GET", f"/api/bed/info/tok/{chrom}", match_info={"token": "tok", "chrom": chrom},
        app={"manager": type("M", (), {"get": staticmethod(lambda _t: renderer)})()},
    )
    return asyncio.run(main.bed_info(req))


def test_bed_info_reports_an_unknown_chromosome_as_404_in_english(tmp_path):
    resp = _bed_info(_bed_renderer(tmp_path), "ZZZ")
    assert resp.status == 404
    body = _json(resp)
    assert body["error"] == "chromosome not found: ZZZ"
    assert body["error"].isascii(), "server errors are English; a client may log them anywhere"


def test_bed_info_still_serves_a_known_chromosome(tmp_path):
    resp = _bed_info(_bed_renderer(tmp_path), "7")
    assert resp.status == 200
    assert _json(resp)["chromosome"] == "7"


# ── open-matrix budget (was: every registered mcool stayed open until process exit) ──


class _FakeMatrix:
    def __init__(self, path):
        self.path = path
        self.closed = False

    def close(self):
        self.closed = True


def _manager_with_fake_matrices(monkeypatch, n):
    mgr = main.GlobalManager()
    monkeypatch.setattr(mgr, "_build",
                        lambda path, file_type, prediction_path=None: _FakeMatrix(path))
    return mgr, [mgr.register(f"/tmp/m{i}.mcool", "hic") for i in range(n)]


def test_only_a_bounded_number_of_matrices_stays_open(monkeypatch):
    mgr, tokens = _manager_with_fake_matrices(monkeypatch, main._MAX_OPEN_MATRICES + 4)
    assert len(mgr.renderers) == main._MAX_OPEN_MATRICES
    assert len(mgr._spec) == len(tokens), "every session must stay addressable"


def test_an_evicted_session_reopens_transparently(monkeypatch):
    """Eviction must never invalidate a token: another tab importing a matrix cannot be
    allowed to break a view that is still on screen."""
    mgr, tokens = _manager_with_fake_matrices(monkeypatch, main._MAX_OPEN_MATRICES + 1)
    oldest = tokens[0]
    assert oldest not in mgr.renderers
    reopened = mgr.get(oldest)
    assert reopened is not None and reopened.path == "/tmp/m0.mcool"
    assert mgr.get(oldest) is reopened


def test_registering_the_same_path_twice_returns_one_session(monkeypatch):
    mgr = main.GlobalManager()
    monkeypatch.setattr(mgr, "_build",
                        lambda path, file_type, prediction_path=None: _FakeMatrix(path))
    assert mgr.register("/tmp/same.mcool", "hic") == mgr.register("/tmp/same.mcool", "hic")
    assert len(mgr._spec) == 1


def test_unknown_tokens_stay_unknown():
    assert main.GlobalManager().get("no-such-token") is None


# ── readiness (was: a flat 200 even with the whole /api/tadvci route table missing) ──


def _health():
    return asyncio.run(main.handle_health(make_mocked_request("GET", "/api/health")))


def test_health_reports_ok_when_every_subsystem_registered():
    main.DEGRADED.clear()
    resp = _health()
    assert resp.status == 200
    assert _json(resp)["status"] == "ok"
    assert _json(resp)["tadvci"] is True


def test_health_fails_closed_when_the_tadvci_routes_are_missing():
    main.DEGRADED.clear()
    main.DEGRADED["tadvci"] = "ImportError: simulated packaging fault"
    try:
        resp = _health()
        body = _json(resp)
        assert resp.status == 503, "a launcher that polls response.ok must not see a ready app"
        assert body["status"] == "degraded"
        assert body["tadvci"] is False
        assert any("tadvci" in e for e in body["errors"])
    finally:
        main.DEGRADED.clear()


def test_fetch_range_data_refuses_a_window_with_too_many_bin_pairs():
    """The response is one dict per surviving pixel, so the RESULT size — not the bin
    span alone — is what has to be bounded: a 4096 x 4096 window inside the per-axis cap
    still answered 335 MB of JSON."""
    stub = _StubMatrix()
    side = int(main._MAX_RANGE_CELLS ** 0.5) + 64
    assert side <= main._MAX_RANGE_BINS, "this case must be inside the per-axis cap"
    with pytest.raises(ValueError, match="bin pairs"):
        asyncio.run(stub.fetch_range_data(25_000, (0, side), (0, side)))
    assert stub.sync_calls == 0


def test_internal_callers_stay_inside_the_bin_pair_budget():
    """A 256-bin tile and the 3 x 256 colour-scale probe must never trip the guard."""
    assert 256 * 256 <= main._MAX_RANGE_CELLS
    assert (3 * 256) ** 2 <= main._MAX_RANGE_CELLS
