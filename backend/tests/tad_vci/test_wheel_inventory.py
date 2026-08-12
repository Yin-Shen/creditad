from __future__ import annotations

import subprocess
import sys
import tomllib
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BACKEND = ROOT / "backend"
EXPECTED_INCLUDE = [
    "tad_vci",
    "tad_vci.*",
    "builtin_callers_v2",
    "builtin_callers_v2.*",
]
EXPECTED_EXCLUDE = [
    "builtin_callers_v2_pre_4method_backup_20260513",
    "builtin_callers_v2_pre_4method_backup_20260513.*",
]
FORBIDDEN_WHEEL_PARTS = (
    "analyses/",
    "builtin_callers/",
    "builtin_callers_v2_pre_4method_backup_20260513/",
    "pre_critfix",
    ".disabled",
    "/tests/",
)


def test_package_discovery_uses_exact_active_names() -> None:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        config = tomllib.load(handle)
    package_find = config["tool"]["setuptools"]["packages"]["find"]
    assert package_find["namespaces"] is False
    assert package_find["include"] == EXPECTED_INCLUDE
    assert package_find["exclude"] == EXPECTED_EXCLUDE

    discovered = {
        str(path.parent.relative_to(BACKEND)).replace("/", ".")
        for path in BACKEND.rglob("__init__.py")
        if path.parent.relative_to(BACKEND).parts[0]
        in {"tad_vci", "builtin_callers_v2"}
    }
    assert discovered == {"tad_vci", "builtin_callers_v2"}


def test_built_wheel_contains_only_active_runtime_packages(tmp_path: Path) -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            str(ROOT),
            "--no-deps",
            "--no-cache-dir",
            "--wheel-dir",
            str(tmp_path),
        ],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    wheels = list(tmp_path.glob("creditad-*.whl"))
    assert len(wheels) == 1

    with zipfile.ZipFile(wheels[0]) as archive:
        names = set(archive.namelist())

    assert "mcool_bed.py" in names
    assert "tad_vci/__init__.py" in names
    assert "tad_vci/calibrated_bands.json" in names
    assert "tad_vci/normalization_spec.json" in names
    assert "builtin_callers_v2/__init__.py" in names
    for name in names:
        normalized = f"/{name}"
        assert not any(part in normalized for part in FORBIDDEN_WHEEL_PARTS), name
