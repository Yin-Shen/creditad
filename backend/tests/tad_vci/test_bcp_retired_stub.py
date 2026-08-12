"""Active-package contract for the retired BCP/conformal component."""
from __future__ import annotations

from pathlib import Path

import pytest

from tad_vci import bcp_infer


def test_bcp_active_module_is_fail_closed() -> None:
    assert bcp_infer.available() is False
    status = bcp_infer.status_for_score(None)
    assert status["status"] == "retired_provenance_failure"
    assert status["probability_or_coverage_claim_allowed"] is False
    for function in (
        bcp_infer.score_features,
        bcp_infer.annotate_bcp,
        bcp_infer.bcp_level,
        bcp_infer.uncertain,
        bcp_infer.concordance,
    ):
        with pytest.raises(RuntimeError, match="retired"):
            function(None)


def test_active_server_has_no_legacy_prediction_or_active_learning_routes() -> None:
    backend = Path(__file__).resolve().parents[2]
    source = (backend / "main.py").read_text()
    assert "handle_deeptad" not in source
    assert "/api/deeptad/" not in source
    assert "ENABLE_LEGACY_AI" not in source
    assert "BoundaryNetPredictor" not in source


def test_backend_root_contains_only_active_runtime_python_modules() -> None:
    backend = Path(__file__).resolve().parents[2]
    observed = {path.name for path in backend.glob("*.py")}
    assert observed == {
        "bigwig_adapter.py",
        # Windows-only child process for winbbi: its Go DLL aborts the interpreter on
        # some calls, so it must not share the server's process (see the module docstring).
        "bigwig_helper.py",
        "data_manager.py",
        "hic_reader.py",
        "main.py",
        "mcool_bed.py",
        "sql.py",
        "startup_log.py",
        "validate_bigwig.py",
    }
