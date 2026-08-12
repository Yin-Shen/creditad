from __future__ import annotations

from pathlib import Path

import pytest

from tad_vci.runtime_contract import MethodSetContractError, resolve_method_set


BACKEND = Path(__file__).resolve().parents[2]


def test_method_set_defaults_to_active_tadvci_not_record_state() -> None:
    assert resolve_method_set({}) == "tadvci"
    assert resolve_method_set({"method_set": None}) == "tadvci"


def test_only_active_and_explicit_v2_method_sets_are_accepted() -> None:
    assert resolve_method_set({"method_set": "tadvci"}) == "tadvci"
    assert resolve_method_set({"method_set": " V2 "}) == "v2"
    with pytest.raises(MethodSetContractError, match="allowed values: tadvci, v2") as exc:
        resolve_method_set({"method_set": "unknown"})
    assert exc.value.status_code == 400
    with pytest.raises(MethodSetContractError, match="non-empty string") as exc:
        resolve_method_set({"method_set": ""})
    assert exc.value.status_code == 400


def test_explicit_legacy_request_is_fail_closed() -> None:
    with pytest.raises(MethodSetContractError, match="retired") as exc:
        resolve_method_set({"method_set": " legacy "})
    assert exc.value.status_code == 410


def test_server_validates_method_before_mutating_task_state() -> None:
    source = (BACKEND / "main.py").read_text(encoding="utf-8")
    handler = source[source.index("async def handle_tad_run"):]
    assert handler.index("method_set = resolve_method_set(data)") < handler.index(
        "mark_running(token)"
    )
    assert "from mcool_bed import detect_TAD_boundaries," not in source
    assert "detector = detect_TAD_boundaries_v2" in handler
