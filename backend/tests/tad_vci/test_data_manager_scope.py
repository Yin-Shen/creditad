from __future__ import annotations

import asyncio
import json
import os
import re
from pathlib import Path

import pytest

from data_manager import (
    DATA_FILES,
    VISIBLE_CATEGORIES,
    VISIBLE_FILE_IDS,
    _md5_file,
    _visible_entries,
    _visible_entry,
    handle_data_download,
    handle_data_status,
    handle_data_verify,
)


class _JsonRequest:
    def __init__(self, body: object) -> None:
        self._body = body

    async def json(self) -> object:
        return self._body


DEMO_CELLS = ("GM12878", "IMR90", "HepG2")

# One concrete entry for the negative verify cases; any visible id would do.
_PROBE_ID = "gm12878_chr7_mcool"


def test_data_manager_lists_every_shipped_track_for_all_three_cell_lines() -> None:
    assert VISIBLE_CATEGORIES == {"shipped_chr7_demo"}
    # One contact matrix + CTCF + RAD21 per cell line. The inventory used to hold
    # GM12878 only, so Data Manager showed one cell line while three shipped.
    assert len(VISIBLE_FILE_IDS) == 3 * len(DEMO_CELLS)
    for cell in DEMO_CELLS:
        for kind in ("chr7_mcool", "ctcf_chr7_bigwig", "rad21_chr7_bigwig"):
            entry = _visible_entry(f"{cell.lower()}_{kind}")
            assert entry is not None, f"{cell} {kind} missing from the inventory"
            assert entry["category"] == "shipped_chr7_demo"
            assert entry["cell_line"] == cell
            assert cell in entry["filename"], (
                f"{entry['filename']} is not {cell}'s own file"
            )

    assert {entry["id"] for entry in _visible_entries()} == VISIBLE_FILE_IDS
    # Every shipped name says chr7, and none of them is an ENCODE full download.
    for entry in _visible_entries():
        assert "chr7" in entry["filename"]
        assert not entry.get("url"), "shipped chr7 tracks are not ENCODE full downloads"

    # The retired GRCh37 matrix was shared by all three cell lines.
    assert _visible_entry("chr7_mcool") is None
    assert all(entry["filename"] != "chr7.mcool" for entry in _visible_entries())
    assert _visible_entry("smc3_bigwig") is None
    assert _visible_entry("ctcf_bigwig") is None  # old id retired
    assert _visible_entry("hg19_fa") is None


def test_status_returns_exact_visible_inventory() -> None:
    response = asyncio.run(handle_data_status(None))
    payload = json.loads(response.text)
    assert {entry["id"] for entry in payload["files"]} == VISIBLE_FILE_IDS
    assert len(payload["files"]) == len(VISIBLE_FILE_IDS)


def test_hidden_entry_is_not_exposed_by_any_data_endpoint(monkeypatch) -> None:
    hidden = {
        "id": "legacy_model",
        "category": "legacy_model",
        "filename": "retired.joblib",
        "dir": "/tmp",
        "url": "https://invalid.example/retired.joblib",
        "md5": "0" * 32,
        "size_mb": 1,
        "required": False,
        "description": "must remain hidden",
    }
    monkeypatch.setattr("data_manager.DATA_FILES", [*DATA_FILES, hidden])

    status_response = asyncio.run(handle_data_status(None))
    status_payload = json.loads(status_response.text)
    assert "legacy_model" not in {entry["id"] for entry in status_payload["files"]}

    download_response = asyncio.run(handle_data_download(_JsonRequest({"id": "legacy_model"})))
    verify_response = asyncio.run(handle_data_verify(_JsonRequest({"id": "legacy_model"})))
    assert download_response.status == 404
    assert verify_response.status == 404


# ======================================================================================
# H9 — the declared md5/size of every visible entry is checked against the real bytes.
#
# The Data Manager refuses a download whose md5 does not match, which is worth nothing
# if the declared md5 was itself copied wrong. These tests recompute it from the local
# copy. When a file is not on this host the check skips instead of passing vacuously.
#
# Two locations count as "the local copy":
#   1. the download destination the entry itself declares (dir/filename);
#   2. the canonical raw copy of the same ENCODE accession already on the host, under
#      $CREDITAD_SHARED_RAW (no default; the test skips when unset) — the demo packages
#      are graded against that copy, so it is the same bytes the entry is a pointer to.
# ======================================================================================

_SHARED_RAW = Path(os.environ.get("CREDITAD_SHARED_RAW", "")) if os.environ.get(
    "CREDITAD_SHARED_RAW") else Path("__unset__")
_ACCESSION = re.compile(r"(ENCFF[0-9A-Z]{6})")


def _local_copies(entry: dict) -> list[Path]:
    candidates = [Path(entry["dir"]) / entry["filename"]]
    match = _ACCESSION.search(entry.get("url") or "")
    if match and _SHARED_RAW.is_dir():
        candidates += sorted(_SHARED_RAW.glob(f"*/*/{match.group(1)}.bigWig"))
    return [p for p in candidates if p.is_file() and p.stat().st_size > 0]


def _visible_ids() -> list[str]:
    return sorted(VISIBLE_FILE_IDS)


@pytest.mark.parametrize("file_id", _visible_ids())
def test_declared_checksum_is_well_formed(file_id: str) -> None:
    """Runs everywhere, file or no file: a malformed md5 can never match anything, so it
    would turn the download integrity check into a permanent failure."""
    entry = _visible_entry(file_id)
    assert re.fullmatch(r"[0-9a-f]{32}", entry["md5"]), entry["md5"]
    assert isinstance(entry["size_mb"], int) and entry["size_mb"] > 0
    assert entry["filename"] and entry["dir"]


@pytest.mark.parametrize("file_id", _visible_ids())
def test_declared_md5_matches_the_local_copy(file_id: str) -> None:
    entry = _visible_entry(file_id)
    copies = _local_copies(entry)
    if not copies:
        pytest.skip(f"{file_id}: no local copy on this host ({entry['filename']})")
    path = copies[0]
    assert _md5_file(path) == entry["md5"], (
        f"{file_id}: declared md5 {entry['md5']} does not match {path}"
    )


@pytest.mark.parametrize("file_id", _visible_ids())
def test_declared_size_matches_the_local_copy(file_id: str) -> None:
    """size_mb is decimal MB (10^6 bytes) and displayed as a whole number, so the
    declaration is accepted within 1 MB of the real length."""
    entry = _visible_entry(file_id)
    copies = _local_copies(entry)
    if not copies:
        pytest.skip(f"{file_id}: no local copy on this host ({entry['filename']})")
    actual_mb = copies[0].stat().st_size / 1e6
    assert abs(actual_mb - entry["size_mb"]) <= 1.0, (
        f"{file_id}: declares {entry['size_mb']} MB, {copies[0]} is {actual_mb:.1f} MB"
    )


@pytest.mark.parametrize("file_id", _visible_ids())
def test_verify_endpoint_really_hashes_a_present_file(file_id: str) -> None:
    """End-to-end through the API the UI calls, not just the helper."""
    entry = _visible_entry(file_id)
    # Resolve through the module's own resolver, not a hand-built path: since the shared
    # raw tree is consulted first, hardcoding `entry["dir"]` here skipped the two ChIP
    # entries even though their bytes were on disk under the canonical accession path.
    from data_manager import _file_path
    declared = _file_path(entry)
    if not (declared.is_file() and declared.stat().st_size > 0):
        pytest.skip(f"{file_id}: no local copy resolved (looked at {declared})")
    payload = json.loads(asyncio.run(handle_data_verify(_JsonRequest({"id": file_id}))).text)
    assert payload["verified"] is True, payload
    assert payload["actual"] == payload["expected"] == entry["md5"]


def test_verify_endpoint_reports_a_corrupted_copy(tmp_path, monkeypatch) -> None:
    """The negative case, with no large file involved: same id, wrong bytes on disk."""
    entry = dict(_visible_entry(_PROBE_ID))
    entry["dir"] = str(tmp_path)
    (tmp_path / entry["filename"]).write_bytes(b"not the real matrix")
    monkeypatch.setattr(
        "data_manager.DATA_FILES",
        [entry if e["id"] == _PROBE_ID else e for e in DATA_FILES],
    )

    payload = json.loads(asyncio.run(handle_data_verify(_JsonRequest({"id": _PROBE_ID}))).text)
    assert payload["verified"] is False
    assert payload["expected"] == entry["md5"]
    assert payload["actual"] != payload["expected"]
    assert "mismatch" in payload["reason"].lower()


def test_verify_endpoint_reports_a_missing_file(tmp_path, monkeypatch) -> None:
    entry = dict(_visible_entry(_PROBE_ID))
    entry["dir"] = str(tmp_path / "absent")
    monkeypatch.setattr(
        "data_manager.DATA_FILES",
        [entry if e["id"] == _PROBE_ID else e for e in DATA_FILES],
    )

    payload = json.loads(asyncio.run(handle_data_verify(_JsonRequest({"id": _PROBE_ID}))).text)
    assert payload["verified"] is False
    assert "not found" in payload["reason"].lower()
