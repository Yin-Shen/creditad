"""The demo packages ship their six caller BEDs gzipped, and the demo must show them.

This repository ships `data/multicell/<CELL>_<RES>/<caller>_boundaries.bed.gz`. Release
v1.0.4's use_multicell_example registered only uncompressed `<caller>_boundaries.bed`, so
a demo loaded from the released tree came up with no caller tracks at all; the desktop
staging worked around it by decompressing the BEDs first.
"""
import asyncio
import gzip
import json
import os
import shutil
from pathlib import Path

from aiohttp.test_utils import make_mocked_request

from tad_vci import api_routes

REPO = Path(__file__).resolve().parents[3]
PACKAGES = REPO / "data" / "multicell"
CALLERS = ["TopDom", "SpectralTAD", "OnTAD", "MSTD", "arrowhead", "DI"]


class _Manager:
    """GlobalManager's register/get contract; records what was registered."""

    def __init__(self):
        self.registered = []

    def register(self, path, file_type, prediction_path=None):
        self.registered.append((path, file_type))
        return f"tok-{len(self.registered)}"

    def get(self, token):
        return None


def _load(monkeypatch, packages_dir: Path, cell="GM12878", res=25000):
    monkeypatch.setattr(api_routes, "_MULTICELL_DIR", packages_dir)
    monkeypatch.setattr(api_routes, "_multicell_track_paths",
                        lambda _cell: {"mcool": None, "ctcf": None, "rad21": None})
    mgr = _Manager()
    req = make_mocked_request("POST", "/api/tadvci/use_multicell_example",
                              app={"manager": mgr})

    async def _j():
        return {"cell_line": cell, "resolution": res}
    req.json = _j
    try:
        resp = asyncio.run(api_routes.handle_tadvci_use_multicell_example(req))
    finally:
        api_routes.IMPORTED_ANNOT.pop(f"multicell:{cell}_{res // 1000}kb", None)
    return resp.status, json.loads(resp.body.decode()), mgr


def _bed_names(mgr):
    return [Path(p).name for p, kind in mgr.registered if kind == "bed"]


def test_shipped_gzipped_caller_beds_become_tracks(monkeypatch):
    for cell, res in (("GM12878", 25000), ("HepG2", 10000)):
        status, payload, mgr = _load(monkeypatch, PACKAGES, cell, res)
        assert status == 200, payload
        assert [t["fileName"] for t in payload["bed_tokens"]] == CALLERS
        assert _bed_names(mgr) == [f"{c}_boundaries.bed.gz" for c in CALLERS]


def test_an_uncompressed_bed_is_still_preferred(monkeypatch, tmp_path):
    """A staged package with decompressed BEDs keeps loading exactly what v1.0.4 loaded."""
    pkg = tmp_path / "GM12878_25kb"
    shutil.copytree(PACKAGES / "GM12878_25kb", pkg)
    plain = pkg / "DI_boundaries.bed"
    with gzip.open(pkg / "DI_boundaries.bed.gz", "rb") as fin:
        plain.write_bytes(fin.read())
    status, payload, mgr = _load(monkeypatch, tmp_path)
    assert status == 200, payload
    names = _bed_names(mgr)
    assert names[-1] == "DI_boundaries.bed"
    assert names[:-1] == [f"{c}_boundaries.bed.gz" for c in CALLERS[:-1]]


def test_bed_track_renderer_reads_gz_like_the_plain_file(tmp_path):
    cwd = os.getcwd()
    import main  # main.py changes into the backend directory at import
    os.chdir(cwd)
    gz = PACKAGES / "GM12878_25kb" / "TopDom_boundaries.bed.gz"
    plain = tmp_path / "TopDom_boundaries.bed"
    with gzip.open(gz, "rb") as fin:
        plain.write_bytes(fin.read())
    a, b = main.BedTileRenderer(str(plain)), main.BedTileRenderer(str(gz))
    assert len(a.df) > 0
    assert a.df.equals(b.df)
    assert a.get_all_chroms() == b.get_all_chroms()
