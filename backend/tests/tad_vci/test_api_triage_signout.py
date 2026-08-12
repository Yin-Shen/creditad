"""API-side tests for retired triage behavior and versioned decision records.

Handlers are aiohttp coroutines; the repo ships no pytest-asyncio / pytest-aiohttp
plugin, so each handler is driven directly with asyncio.run() + a mocked request
(query params via the path, JSON body via a stubbed .json()). The handlers only touch
request.match_info / request.query / request.json() — never request.app — so a mocked
request is sufficient and no live server is needed.

Covered regressions:
  1. The scientifically invalid legacy BCP/conformal API is fail-closed even when a
     historical annotation still carries score columns.
  2. No coverage or auto-confirmation claim is returned for any sample.
  3. Sign-out stamps the assembly from the imported annotation TSV (not a fixed hg19
     default), so an hg38 upload cannot be hashed as hg19.
"""
import asyncio
import json
import sys
from pathlib import Path

import pandas as pd
from aiohttp.test_utils import make_mocked_request

BACKEND = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BACKEND))

from tad_vci import api_routes    # noqa: E402
from mcool_bed import detect_TAD_boundaries_tadvci  # noqa: E402

def _point_token(token: str, tsv: Path) -> None:
    api_routes.set_imported_annotation(token, str(tsv))


def _triage(token: str, chrom: str | None = None) -> dict:
    path = f"/api/tadvci/triage/{token}" + (f"?chrom={chrom}" if chrom else "")
    req = make_mocked_request("GET", path, match_info={"token": token})
    resp = asyncio.run(api_routes.handle_tadvci_triage(req))
    return json.loads(resp.body.decode())


def _signout(body: dict, log_path: Path | None = None) -> tuple[int, dict]:
    req = make_mocked_request("POST", "/api/tadvci/signout")
    async def _json():
        return body
    req.json = _json
    # Never append to the application decision log from a test; redirect to a temp log.
    orig_log = api_routes._LOG
    if log_path is not None:
        api_routes._LOG = log_path
    try:
        resp = asyncio.run(api_routes.handle_tadvci_signout(req))
    finally:
        api_routes._LOG = orig_log
    return resp.status, json.loads(resp.body.decode())


def test_triage_three_buckets_partition_n(tmp_path):
    """Historical score columns cannot reactivate retired triage claims."""
    df = _hic_only_df().assign(bcp_score=[0.1, 0.5, 0.9], bcp_uncertain=[False, True, False])
    path = tmp_path / "legacy_score_columns.tsv"
    df.to_csv(path, sep="\t", index=False)
    n_real = len(df)
    _point_token("tk_chr7", path)
    j = _triage("tk_chr7")
    assert j["n"] == n_real
    assert j["scored"] == 0
    assert j["not_scored"] == n_real
    assert j["auto_confident"] == j["routed_to_expert"] == 0
    assert j["bcp_available"] is False
    assert j["coverage_guarantee"] is None
    assert j["bcp_status"].startswith("retired_")


def test_triage_routed_frac_uses_scored_denominator(tmp_path):
    """Chromosome scoping remains available, but no retired ML quantity is returned."""
    path = tmp_path / "legacy_score_columns.tsv"
    _hic_only_df().assign(bcp_score=[0.1, 0.5, 0.9]).to_csv(path, sep="\t", index=False)
    _point_token("tk_chr7b", path)
    j = _triage("tk_chr7b", chrom="7")
    assert j["scope"] == "chr7"
    assert j["n"] > 0
    assert j["scored"] == 0
    assert j["routed_frac"] is None
    assert j["coverage_guarantee"] is None


def _hic_only_df() -> pd.DataFrame:
    """A Hi-C-only evidence annotation without retired model columns."""
    return pd.DataFrame({
        "chrom": ["7", "7", "7"],
        "pos": [1_000_000, 2_000_000, 3_000_000],
        "votes": [4, 2, 5],
        "loop_anchor_count": [1, 0, 2],
        "D_tier": ["D4", "D2", "D5"],
    })


def test_tadvci_detection_disables_legacy_bcp_by_default():
    """The active scientific import path cannot silently reactivate the retired model."""
    import inspect
    assert inspect.signature(detect_TAD_boundaries_tadvci).parameters["enable_bcp"].default is False


def test_builtin_registry_and_decision_log_do_not_serve_legacy_artifacts(monkeypatch):
    """Static demos fail closed until the manifest-gated ICE annotations exist."""
    assert set(api_routes.CELL_REGISTRY) == {"GM12878"}
    for record in api_routes.CELL_REGISTRY.values():
        assert "tad_vci_rebuild_20260710" in record["annotation"]
        assert "cells_ice_v2" in record["annotation"]
        assert "tad_vci_20260529" not in record["annotation"]
    assert api_routes._LOG.name == "decision_records_v3.jsonl"

    def _reject(_root):
        raise api_routes.BundleEligibilityError("canonical manifest is missing")

    monkeypatch.setattr(api_routes, "validate_gm12878_bundle", _reject)
    for cell, record in api_routes.CELL_REGISTRY.items():
        eligible, reason = api_routes._cell_manifest_status(cell, record)
        assert eligible is False
        assert "canonical bundle verification failed" in reason


def _verified_api_bundle(root: Path, annotation: Path) -> dict:
    manifest_path = root / "GM12878_run_manifest.json"
    manifest_path.write_text("{}")
    return {
        "root": root,
        "manifest_path": manifest_path,
        "manifest": {
            "schema": "creditad-gm12878-cell-run-manifest-v2",
            "status": "passed_all_eligibility_gates",
            "cell_line": "GM12878",
            "assembly": "hg19",
            "resolution_bp": 25_000,
            "full_autosome_validation_required": True,
            "exact_reproducibility_required": True,
        },
        "annotation_path": annotation,
        "method_paths": {},
        "inputs": {
            "ice_sidecar_validation": root / "validation-v2.json",
            "ice_reproducibility": root / "reproducibility.json",
            "calibration": root / "calibration.json",
            "calibration_manifest": root / "calibration.manifest.json",
        },
        "validation": {
            "schema": "creditad-ice-sidecar-validation-v2",
            "pass": True,
        },
        "reproducibility": {
            "schema": "creditad-ice-reproducibility-check-v1",
            "pass": True,
        },
    }


def test_cell_registry_uses_complete_shared_v2_verifier(tmp_path, monkeypatch):
    root = tmp_path / "cells_ice_v2"
    root.mkdir()
    monkeypatch.setattr(api_routes, "_CELLS_DIR", root)
    annotation = root / "GM12878_annotation.tsv"
    annotation.write_text("chrom\tpos\tvotes\n1\t25000\t1\n")
    verified = _verified_api_bundle(root, annotation)
    calls = []

    def _verify(requested_root):
        calls.append(Path(requested_root))
        return verified

    monkeypatch.setattr(api_routes, "validate_gm12878_bundle", _verify)
    reg = {"annotation": str(annotation), "assembly": "hg19", "resolution": 25000}
    assert api_routes._cell_manifest_status("GM12878", reg) == (True, "eligible")
    assert calls == [root]
    assert verified["manifest"]["schema"] == "creditad-gm12878-cell-run-manifest-v2"
    assert verified["validation"]["schema"] == "creditad-ice-sidecar-validation-v2"
    assert verified["reproducibility"]["pass"] is True
    assert {"calibration", "calibration_manifest"}.issubset(verified["inputs"])


def test_cell_registry_fails_closed_on_v2_provenance_gate(tmp_path, monkeypatch):
    root = tmp_path / "cells_ice_v2"
    root.mkdir()
    monkeypatch.setattr(api_routes, "_CELLS_DIR", root)
    annotation = root / "GM12878_annotation.tsv"
    annotation.write_text("chrom\tpos\tvotes\n1\t25000\t1\n")

    def _reject(_root):
        raise api_routes.BundleEligibilityError(
            "full-autosome ICE validation is stale; reproducibility/calibration not trusted"
        )

    monkeypatch.setattr(api_routes, "validate_gm12878_bundle", _reject)
    reg = {"annotation": str(annotation), "assembly": "hg19", "resolution": 25000}
    eligible, reason = api_routes._cell_manifest_status("GM12878", reg)
    assert eligible is False
    assert "validation is stale" in reason


def test_cell_registry_rejects_non_gm_and_registry_drift(tmp_path, monkeypatch):
    root = tmp_path / "cells_ice_v2"
    root.mkdir()
    annotation = root / "GM12878_annotation.tsv"
    annotation.write_text("chrom\tpos\tvotes\n1\t25000\t1\n")
    verified = _verified_api_bundle(root, annotation)
    monkeypatch.setattr(api_routes, "_CELLS_DIR", root)
    monkeypatch.setattr(api_routes, "validate_gm12878_bundle", lambda _root: verified)
    reg = {"annotation": str(annotation), "assembly": "hg19", "resolution": 25000}
    assert api_routes._cell_manifest_status("K562", reg)[0] is False
    drifted = dict(reg, assembly="hg38")
    eligible, reason = api_routes._cell_manifest_status("GM12878", drifted)
    assert eligible is False and "assembly differs" in reason


def test_public_record_normalizes_legacy_expert_field_names():
    public = api_routes._public_record({
        "boundary_id": "1:25000",
        "expert_override": {"final_tier": "D3", "verdict": "legacy", "rationale": "old"},
        "curator": "legacy_demo_label",
    })
    assert "expert_override" not in public and "curator" not in public
    assert public["review_decision"]["rationale"] == "old"
    assert public["actor_label"] == "legacy_demo_label"
    assert public["legacy_field_names_normalized"] is True


def test_decision_requires_neutral_verdict_and_explicit_actor_label():
    status, payload = _signout({
        "boundary_id": "1:25000",
        "final_tier": "D3",
        "rationale": "review note",
        "verdict": "supported_for_review",
    })
    assert status == 422
    assert "actor_label required" in payload["error"]

    status, payload = _signout({
        "boundary_id": "1:25000",
        "final_tier": "D3",
        "rationale": "review note",
        "verdict": "confirmed_real",
        "actor_label": "local_user",
    })
    assert status == 422
    assert "verdict must be one of" in payload["error"]


def test_triage_hic_only_reports_bcp_unavailable(tmp_path):
    """Issue 2 (API side): a Hi-C-only annotation makes the triage endpoint report
    bcp_available == False + an explicit 'BCP not available' status, not a 90% badge."""
    out = _hic_only_df()
    tsv = tmp_path / "hic_only.tsv"
    out.to_csv(tsv, sep="\t", index=False)
    _point_token("tk_hic_only", tsv)
    j = _triage("tk_hic_only")
    assert j["scored"] == 0
    assert j["not_scored"] == j["n"] == 3
    assert j["bcp_available"] is False
    assert j["bcp_status"].startswith("retired_")
    assert j["routed_frac"] is None
    assert j["coverage_guarantee"] is None
    assert j["auto_confident"] + j["routed_to_expert"] + j["not_scored"] == j["n"]


def test_signout_assembly_from_annotation_not_hardcoded(tmp_path):
    """Issue 3: sign-out reads assembly from the imported annotation TSV. A token whose
    annotation carries assembly='hg38' must be hashed as hg38 even though the served cell
    resolves to GM12878 (hg19) by default."""
    sub = pd.DataFrame({
        "chrom": ["1"], "pos": [1_000_000], "votes": [3],
        "BD1_caller_support": ["moderate"], "BD2_ctcf": ["not_assessable"],
        "BD3_rad21": ["not_assessable"], "BD4_loop_anchor": ["none"],
        "D_tier": ["D3"], "assembly": ["hg38"],
    })
    tsv = tmp_path / "hg38_upload.tsv"
    sub.to_csv(tsv, sep="\t", index=False)
    _point_token("tk_hg38", tsv)
    bid = f"1:{int(sub.iloc[0]['pos'])}"
    status, j = _signout({
        "boundary_id": bid, "final_tier": "D3",
        "rationale": "test hg38 provenance", "token": "tk_hg38",
        "verdict": "supported_for_review", "actor_label": "test_actor",
    }, log_path=tmp_path / "log.jsonl")
    assert status == 200, j
    assert j["ok"] is True
    assert j["snapshot"]["assembly"] == "hg38"     # NOT the hg19 default
    assert j["checksum_matches"] is True
    assert j["snapshot"]["actor_label"] == "test_actor"
    assert "expert_override" not in j["snapshot"]


def test_signout_assembly_payload_fallback(tmp_path):
    """Issue 3 fallback: when the annotation has no assembly column, an explicit assembly in
    the sign-out payload is honoured (frontend-supplied fallback), still never silent hg19."""
    df = pd.DataFrame({
        "chrom": ["3"], "pos": [5_000_000], "votes": [3],
        "BD1_caller_support": ["moderate"], "D_tier": ["D4"],
    })
    tsv = tmp_path / "no_assembly.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    _point_token("tk_pl", tsv)
    status, j = _signout({
        "boundary_id": "3:5000000", "final_tier": "D4",
        "rationale": "payload assembly fallback", "token": "tk_pl",
        "assembly": "GRCh38", "verdict": "uncertain",
        "actor_label": "test_actor",
    }, log_path=tmp_path / "log.jsonl")
    assert status == 200, j
    assert j["snapshot"]["assembly"] == "hg38"     # GRCh38 alias normalised, not hg19


def test_norm_assembly_helper():
    """_norm_assembly: GRCh37/GRCh38 alias + case-insensitive; unknown -> None (fall through)."""
    assert api_routes._norm_assembly("hg38") == "hg38"
    assert api_routes._norm_assembly("GRCh38") == "hg38"
    assert api_routes._norm_assembly("GRCh37") == "hg19"
    assert api_routes._norm_assembly("HG19") == "hg19"
    assert api_routes._norm_assembly("mm10") is None
    assert api_routes._norm_assembly(None) is None
