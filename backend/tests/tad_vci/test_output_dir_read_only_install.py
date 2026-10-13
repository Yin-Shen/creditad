"""A relative results directory must be writable in a read-only install.

main.py changes into the backend directory at start-up, so a relative output_dir resolves
inside the install tree. The Linux AppImage mounts that tree read-only, and the Import
page's default "./tad_results" was answered with HTTP 400 "output_dir cannot be created:
Read-only file system". When the
working directory is not writable, a relative output_dir now resolves under the per-user
data directory that sql.py already uses for its token database. With a writable working
directory nothing changes.
"""
import asyncio
import json
import os
from pathlib import Path

import pytest
from aiohttp.test_utils import make_mocked_request

from tad_vci import api_routes


@pytest.fixture
def user_data(tmp_path, monkeypatch):
    base = tmp_path / "xdg_data"
    monkeypatch.setenv("XDG_DATA_HOME", str(base))
    monkeypatch.setenv("LOCALAPPDATA", str(base))
    return base / "CrediTAD"


@pytest.fixture
def read_only_cwd(tmp_path, monkeypatch):
    """The working directory reports not writable (as in the mounted AppImage)."""
    install = tmp_path / "install" / "backend"
    install.mkdir(parents=True)
    monkeypatch.chdir(install)
    real = api_routes._dir_writable
    monkeypatch.setattr(api_routes, "_dir_writable",
                        lambda p: False if Path(p) == install else real(p))
    return install


def test_relative_output_dir_goes_to_the_user_data_dir(read_only_cwd, user_data):
    resolved = api_routes._resolved_output_dir("./tad_results")
    assert resolved == str(user_data / "tad_results")
    assert Path(resolved).is_dir()
    assert not (read_only_cwd / "tad_results").exists()


def test_absolute_output_dir_is_never_moved(read_only_cwd, user_data, tmp_path):
    target = tmp_path / "chosen" / "out"
    assert api_routes._resolved_output_dir(str(target)) == str(target)
    assert not user_data.exists()


def test_writable_working_directory_keeps_the_v1_0_4_behaviour(tmp_path, monkeypatch,
                                                                user_data):
    monkeypatch.chdir(tmp_path)
    assert api_routes.relocate_relative_output_dir("./tad_results") == "./tad_results"
    resolved = api_routes._resolved_output_dir("./tad_results")
    assert resolved == os.path.abspath(str(tmp_path / "tad_results"))
    assert not user_data.exists()


@pytest.mark.skipif(not hasattr(os, "geteuid") or os.geteuid() == 0,
                    reason="permission bits do not restrict root")
def test_a_really_read_only_working_directory(tmp_path, monkeypatch, user_data):
    install = tmp_path / "ro"
    install.mkdir()
    monkeypatch.chdir(install)
    install.chmod(0o555)
    try:
        resolved = api_routes._resolved_output_dir("tad_results")
    finally:
        install.chmod(0o755)
    assert resolved == str(user_data / "tad_results")
    assert Path(resolved).is_dir()


class _Manager:
    def __init__(self):
        self.by_path = {}

    def register(self, path, file_type, prediction_path=None):
        return self.by_path.setdefault(path, f"tok-{len(self.by_path)}")

    def get(self, token):
        return None


def test_import_route_default_output_dir_in_a_read_only_install(read_only_cwd, user_data,
                                                                 tmp_path):
    """The request the Import page sends when the user leaves output_dir at its default
    returned HTTP 400 from the mounted AppImage."""
    primary = tmp_path / "primary.bed"
    primary.write_text("7\t1000000\t1025000\n7\t2000000\t2025000\n")
    req = make_mocked_request("POST", "/api/tadvci/annotate", app={"manager": _Manager()})

    async def _j():
        return {"primary_bed": str(primary), "hash_tag": "ro", "output_dir": "./tad_results"}
    req.json = _j
    resp = asyncio.run(api_routes.handle_tadvci_annotate(req))
    payload = json.loads(resp.body.decode())
    try:
        assert resp.status == 200, payload
        assert Path(payload["TADVCI_annotation"]).parent == user_data / "tad_results"
        assert Path(payload["TADVCI_annotation"]).is_file()
    finally:
        api_routes.IMPORTED_ANNOT.pop(payload.get("session_token"), None)


def test_tad_run_route_uses_the_same_location(read_only_cwd, user_data, monkeypatch):
    """/api/tad/run (the .mcool import) wrote to the same relative default."""
    cwd = os.getcwd()
    import main  # main.py changes into the backend directory at import
    os.chdir(cwd)
    seen = {}

    def _detector(**kw):
        seen["output_dir"] = kw["output_dir"]
        return {}

    monkeypatch.setattr(main, "detect_TAD_boundaries_tadvci", _detector)
    for name in ("mark_running", "mark_failed", "update_tad_outputs"):
        monkeypatch.setattr(main, name, lambda *a, **k: None)
    monkeypatch.setattr(main, "get_record", lambda token: {})

    class _Renderer:
        path = "matrix.mcool"

    class _Mgr:
        def get(self, token):
            return _Renderer()

    req = make_mocked_request("POST", "/api/tad/run", app={"manager": _Mgr()})

    async def _j():
        return {"token": "t1", "output_dir": "./tad_results"}
    req.json = _j
    resp = asyncio.run(main.handle_tad_run(req))
    assert resp.status == 200, resp.body
    assert seen["output_dir"] == str(user_data / "tad_results")
    assert (user_data / "tad_results").is_dir()
    assert not (read_only_cwd / "tad_results").exists()
