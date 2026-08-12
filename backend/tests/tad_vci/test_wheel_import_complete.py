"""The built wheel must be IMPORT-COMPLETE.

This is the regression test for BUG-1 (2026-08-12): `tad_vci/annotate_from_beds.py`
imported the top-level module `validate_bigwig` unguarded whenever `--assembly` was
passed, but `[tool.setuptools] py-modules` shipped only `mcool_bed`. A core
`pip install .` therefore produced a CLI that raised
`ModuleNotFoundError: No module named 'validate_bigwig'` on the documented command.

The sibling test module `test_wheel_inventory.py` asserts what must NOT be in the wheel.
Nothing asserted that everything the wheel IMPORTS is actually IN it, which is why the
defect reached a release candidate.

Two checks, both static (no build, no network, fast):

1. `test_declared_py_modules_cover_top_level_imports` — parse every packaged source file,
   collect its absolute top-level imports, and require that any name which resolves to a
   sibling module in `backend/` is declared in `py-modules`. This is the check that fails
   on the BUG-1 configuration.
2. `test_third_party_imports_are_declared_or_guarded` — an import of a third-party
   distribution that is not in `[project] dependencies` must be inside a try/except or a
   function body, so it degrades or names its extra rather than crashing at import time.
   This is the check that fails on the BUG-2 configuration (`pyBigWig` used on the
   documented path while declared only in an extra).
"""
from __future__ import annotations

import ast
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
BACKEND = ROOT / "backend"

# Third-party names that are legitimately optional: they are declared in an extra AND
# every import site is guarded (try/except) or lazy (inside a function).
# Deliberately optional: declared in an extra AND the code copes when absent.
# pyBigWig is NOT in this register on purpose — BD2/BD3 are on the documented CLI path,
# so it must stay a CORE dependency. Moving it here would silence the BUG-2 check.
OPTIONAL_THIRD_PARTY = {
    "aiohttp", "aiohttp_cors", "sqlalchemy", "cooltools", "hicstraw", "straw",
    "winbbi", "pytest", "sklearn", "joblib",
}
STDLIB = set(sys.stdlib_module_names)


def _config() -> dict:
    with (ROOT / "pyproject.toml").open("rb") as fh:
        return tomllib.load(fh)


def _packaged_sources() -> list[Path]:
    """Every .py file that ends up inside the wheel."""
    out: list[Path] = []
    for pkg in ("tad_vci", "builtin_callers_v2"):
        out.extend(sorted((BACKEND / pkg).rglob("*.py")))
    cfg = _config()
    for mod in cfg["tool"]["setuptools"].get("py-modules", []):
        p = BACKEND / f"{mod}.py"
        if p.is_file():
            out.append(p)
    return out


def _toplevel_imports(tree: ast.AST) -> set[tuple[str, bool]]:
    """(module_name, is_module_scope_and_unguarded) for each absolute import."""
    found: set[tuple[str, bool]] = set()

    class V(ast.NodeVisitor):
        def __init__(self) -> None:
            self.depth = 0          # inside try / function / class

        def _add(self, name: str) -> None:
            found.add((name.split(".")[0], self.depth == 0))

        def visit_Try(self, node: ast.Try) -> None:
            self.depth += 1
            for child in ast.iter_child_nodes(node):
                self.visit(child)
            self.depth -= 1

        def visit_FunctionDef(self, node) -> None:
            self.depth += 1
            for child in ast.iter_child_nodes(node):
                self.visit(child)
            self.depth -= 1

        visit_AsyncFunctionDef = visit_FunctionDef

        def visit_Import(self, node: ast.Import) -> None:
            for a in node.names:
                self._add(a.name)

        def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
            if node.level == 0 and node.module:
                self._add(node.module)

    V().visit(tree)
    return found


def test_declared_py_modules_cover_top_level_imports() -> None:
    """Any sibling module the wheel imports must itself be shipped. (BUG-1)"""
    cfg = _config()
    declared = set(cfg["tool"]["setuptools"].get("py-modules", []))
    packages = {"tad_vci", "builtin_callers_v2"}
    siblings = {p.stem for p in BACKEND.glob("*.py")}

    missing: dict[str, list[str]] = {}
    for src in _packaged_sources():
        try:
            tree = ast.parse(src.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):  # pragma: no cover
            continue
        for name, _unguarded in _toplevel_imports(tree):
            if name in declared or name in packages or name in STDLIB:
                continue
            if name in siblings:
                missing.setdefault(name, []).append(
                    str(src.relative_to(ROOT)))

    assert not missing, (
        "the wheel imports sibling modules it does not ship — add them to "
        "[tool.setuptools] py-modules in pyproject.toml:\n"
        + "\n".join(f"  {name}  <- imported by {', '.join(sorted(set(where)))}"
                    for name, where in sorted(missing.items()))
    )


def test_third_party_imports_are_declared_or_deliberately_optional() -> None:
    """Every third-party import in the wheel is either a core dependency or a KNOWN
    optional. (BUG-2)

    Guarding is not sufficient on its own. BUG-2 was a *lazy* import — `pyBigWig` is
    imported inside `mcool_bed._bw_window_signal`, so it is never a module-scope crash —
    yet it still killed the documented `creditad annotate --ctcf` run the moment that
    function was called, because the package was declared only in an extra. A test that
    looked only at module-scope imports would have PASSED on the BUG-2 configuration.

    So: a third-party name imported anywhere in packaged source must be in
    `[project] dependencies`, or be listed in OPTIONAL_THIRD_PARTY — the explicit register
    of "this may legitimately be absent and the code copes". Adding a name there is a
    deliberate act; forgetting to declare a dependency is not.
    """
    cfg = _config()
    deps = cfg["project"]["dependencies"]
    declared = {d.split(";")[0].split("[")[0]
                 .split(">")[0].split("<")[0].split("=")[0].split("!")[0].strip().lower()
                for d in deps}
    # distribution name -> import name, where they differ
    # distribution name -> import name, where PyPI name != module name
    alias = {"pybigwig": "pyBigWig", "scikit-learn": "sklearn",
             "scikit-image": "skimage", "pyyaml": "yaml", "pillow": "PIL",
             "opencv-python": "cv2", "attrs": "attr", "protobuf": "google"}
    declared |= {alias[d] for d in list(declared) if d in alias}
    declared_import_names = {d.replace("-", "_") for d in declared} | {
        v for k, v in alias.items() if k in declared}

    packages = {"tad_vci", "builtin_callers_v2"}
    siblings = {p.stem for p in BACKEND.glob("*.py")}

    offenders: dict[str, list[str]] = {}
    for src in _packaged_sources():
        try:
            tree = ast.parse(src.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):  # pragma: no cover
            continue
        for name, _unguarded in _toplevel_imports(tree):
            low = name.lower()
            if (name in STDLIB or name in packages or name in siblings
                    or low in declared or name in declared_import_names
                    or low in OPTIONAL_THIRD_PARTY):
                continue
            offenders.setdefault(name, []).append(str(src.relative_to(ROOT)))

    assert not offenders, (
        "packaged modules import third-party packages that are neither in "
        "[project] dependencies nor registered as deliberately optional — declare them in "
        "core dependencies, or add them to OPTIONAL_THIRD_PARTY and make the code cope "
        "with their absence:\n"
        + "\n".join(f"  {name}  <- {', '.join(sorted(set(where)))}"
                    for name, where in sorted(offenders.items()))
    )


def test_documented_cli_path_imports_are_shippable() -> None:
    """The modules the documented `creditad annotate` path needs are all shipped."""
    cfg = _config()
    declared = set(cfg["tool"]["setuptools"].get("py-modules", []))
    for required in ("mcool_bed", "validate_bigwig"):
        assert required in declared, (
            f"{required} is on the documented `creditad annotate --assembly` path "
            f"but is not in py-modules; this is exactly BUG-1")
