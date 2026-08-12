"""Shared fail-closed verifier for the canonical GM12878 rebuild bundle.

This module is deliberately Python-3.8-compatible because the pinned ConsTADs
runtime imports it.  It performs no analysis; it only verifies current files
against the staged-build manifest and its normalization/calibration attestations.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from tad_vci.normalization_spec import (
    ACCEPTANCE_GATES,
    ACTIVE_SPEC_IDENTITY,
    NORMALIZATION_ID,
    PARAMETERS,
    REQUIRED_COOLER_VERSION,
    RESOLUTION_BP,
)


BACKEND = Path(__file__).resolve().parents[1]
REBUILD_ROOT = BACKEND / "analyses" / "tad_vci_rebuild_20260710"
CANONICAL_CELL_ROOT = REBUILD_ROOT / "cells_ice_v2"
BUILDER = REBUILD_ROOT / "build_gm12878_annotations.py"
COMPUTE_SCRIPT = REBUILD_ROOT / "compute_ice_sidecar.py"
COMPARE_SCRIPT = REBUILD_ROOT / "compare_ice_replicates.py"
VALIDATION_SCRIPT = REBUILD_ROOT / "validate_ice_sidecar.py"
MCOOL_READER_SCRIPT = BACKEND / "mcool_bed.py"
PACKAGED_BANDS = BACKEND / "tad_vci" / "calibrated_bands.json"
AUTOSOMES = [str(value) for value in range(1, 23)]
METHODS = {
    "insulation",
    "topdom_like",
    "contact_contrast",
    "spectral_profile",
    "network_modularity",
}


class BundleEligibilityError(ValueError):
    """Raised when the canonical rebuilt bundle is missing, stale, or altered."""


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256": sha256_file(path),
    }


def _reject_constant(value):
    raise BundleEligibilityError("non-finite JSON constant is forbidden: %s" % value)


def _require_finite_json(value, location="$"):
    if isinstance(value, float) and not math.isfinite(value):
        raise BundleEligibilityError("non-finite JSON number at %s" % location)
    if isinstance(value, dict):
        for key, child in value.items():
            _require_finite_json(child, "%s.%s" % (location, key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _require_finite_json(child, "%s[%d]" % (location, index))


def load_json_strict(path):
    try:
        payload = json.loads(
            Path(path).read_text(encoding="utf-8"),
            parse_constant=_reject_constant,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise BundleEligibilityError("invalid strict JSON %s: %s" % (path, exc))
    if not isinstance(payload, dict):
        raise BundleEligibilityError("expected JSON object: %s" % path)
    _require_finite_json(payload)
    return payload


def _identity_from_metadata(metadata, label, require_hash=True):
    if not isinstance(metadata, dict) or not metadata.get("path"):
        raise BundleEligibilityError("%s identity is missing" % label)
    path = Path(str(metadata["path"]))
    if not path.is_file():
        raise BundleEligibilityError("%s file is missing: %s" % (label, path))
    actual = file_identity(path)
    for key in ("path", "size_bytes", "mtime_ns"):
        if key in metadata and metadata.get(key) != actual[key]:
            raise BundleEligibilityError("%s %s mismatch" % (label, key))
    if require_hash and metadata.get("sha256") != actual["sha256"]:
        raise BundleEligibilityError("%s SHA-256 mismatch" % label)
    return actual


def _selected_identity(recorded, actual, label):
    for key in ("path", "size_bytes", "mtime_ns", "sha256"):
        if recorded.get(key) != actual[key]:
            raise BundleEligibilityError("%s %s mismatch" % (label, key))


def _verify_reproducibility(record, identities):
    if (
        record.get("schema") != "creditad-ice-reproducibility-check-v1"
        or record.get("cell_line") != "GM12878"
        or record.get("pass") is not True
    ):
        raise BundleEligibilityError("ICE reproducibility record is failing")
    required_checks = {
        "all_invariant_manifest_fields_equal",
        "weight_shape_equal",
        "primary_has_no_infinite_weights",
        "repeat_has_no_infinite_weights",
        "weight_nan_mask_equal",
        "finite_weight_values_bitwise_equal",
        "weight_arrays_equal_including_nan",
        "finite_weight_max_abs_difference_zero",
    }
    checks = record.get("checks", {})
    if set(checks) != required_checks or not all(value is True for value in checks.values()):
        raise BundleEligibilityError("ICE reproducibility checks are incomplete")
    _selected_identity(
        record.get("primary", {}).get("sidecar", {}),
        identities["ice_sidecar"],
        "reproducibility primary sidecar",
    )
    _selected_identity(
        record.get("primary", {}).get("manifest", {}),
        identities["ice_sidecar_manifest"],
        "reproducibility primary manifest",
    )
    repeat_sidecar = _identity_from_metadata(
        record.get("repeat", {}).get("sidecar", {}), "repeat sidecar"
    )
    _identity_from_metadata(
        record.get("repeat", {}).get("manifest", {}), "repeat manifest"
    )
    _selected_identity(
        record.get("comparison_implementation", {}),
        file_identity(COMPARE_SCRIPT),
        "reproducibility comparator",
    )
    weights = record.get("weights", {})
    if (
        weights.get("n_infinite_primary") != 0
        or weights.get("n_infinite_repeat") != 0
        or weights.get("max_abs_difference") != 0.0
        or repeat_sidecar["sha256"] != identities["ice_sidecar"]["sha256"]
    ):
        raise BundleEligibilityError("reproducibility weights are not exact and finite")


def _verify_validation(record, identities, sidecar_manifest):
    per_autosome = record.get("per_autosome", {})
    if (
        record.get("schema") != "creditad-ice-sidecar-validation-v2"
        or record.get("cell_line") != "GM12878"
        or record.get("pass") is not True
        or record.get("validated_autosomes") != AUTOSOMES
        or set(per_autosome) != set(AUTOSOMES)
        or not all(row.get("pass") is True for row in per_autosome.values())
    ):
        raise BundleEligibilityError("full-autosome ICE validation is incomplete")
    _selected_identity(
        record.get("source", {}), identities["mcool"], "validation source mcool"
    )
    compute_sha = sidecar_manifest.get("implementation", {}).get("compute_script", {}).get("sha256")
    if (
        record.get("source", {}).get("pre_run_sha256_match") is not True
        or record.get("source", {}).get("unchanged_during_validation") is not True
        or record.get("sidecar_sha256") != identities["ice_sidecar"]["sha256"]
        or record.get("sidecar_manifest_sha256")
        != identities["ice_sidecar_manifest"]["sha256"]
        or record.get("reproducibility_record_sha256")
        != identities["ice_reproducibility"]["sha256"]
        or record.get("normalization_id") != NORMALIZATION_ID
        or record.get("normalization_spec_sha256") != ACTIVE_SPEC_IDENTITY["sha256"]
        or record.get("compute_script_sha256") != compute_sha
        or compute_sha != sha256_file(COMPUTE_SCRIPT)
        or record.get("comparison_script_sha256") != sha256_file(COMPARE_SCRIPT)
        or record.get("validation_script_sha256") != sha256_file(VALIDATION_SCRIPT)
        or record.get("mcool_reader_sha256") != sha256_file(MCOOL_READER_SCRIPT)
        or record.get("compute_implementation_matches_current") is not True
        or record.get("algorithm_parameters_match_active_spec") is not True
        or record.get("full_manifest_acceptance_gates_rechecked") is not True
        or record.get("actual_sidecar_coverage_rechecked") is not True
        or record.get("reproducibility_gate_rechecked") is not True
        or record.get("acceptance", {}).get("all_22_autosomes_pass") is not True
        or record.get("dense_reference_check", {}).get("matrix_allclose") is not True
    ):
        raise BundleEligibilityError("full-autosome ICE validation is stale")


def _verify_calibration(calibration, manifest, packaged, identities):
    expected = {
        "GM12878": packaged["GM12878"],
        "_engine_default": packaged["_engine_default"],
        "_calibration_spec": packaged["_calibration_spec"],
    }
    comparison = manifest.get("reference_comparison", {})
    if calibration != expected or (
        manifest.get("schema") != "creditad-gm12878-calibration-provenance-v1"
        or manifest.get("status") != "passed_exact_reference_match"
        or manifest.get("cell_line") != "GM12878"
        or comparison.get("gm_record_exact_equal") is not True
        or comparison.get("calibration_spec_exact_equal") is not True
        or comparison.get("engine_default_exact_equal") is not True
        or comparison.get("no_parameter_tuning_performed") is not True
    ):
        raise BundleEligibilityError("GM calibration provenance is incomplete")
    output = manifest.get("output", {})
    for key in ("path", "size_bytes", "sha256"):
        if output.get(key) != identities["calibration"][key]:
            raise BundleEligibilityError("calibration output %s mismatch" % key)
    _selected_identity(
        manifest.get("implementation", {}).get("generator", {}),
        identities["calibration_generator"],
        "calibration generator",
    )
    _selected_identity(
        manifest.get("implementation", {}).get("packaged_reference", {}),
        identities["packaged_calibrated_bands"],
        "packaged calibration",
    )
    for assay, key in (("CTCF", "ctcf_bigwig"), ("RAD21", "rad21_bigwig")):
        recorded = manifest.get("inputs", {}).get(assay, {})
        _selected_identity(recorded, identities[key], "%s calibration input" % assay)
        if recorded.get("stable_before_and_after") is not True:
            raise BundleEligibilityError("%s calibration input was unstable" % assay)


def validate_gm12878_bundle(cell_root=CANONICAL_CELL_ROOT):
    root = Path(cell_root).resolve()
    if root != CANONICAL_CELL_ROOT.resolve() or root.name != "cells_ice_v2":
        raise BundleEligibilityError("only the canonical cells_ice_v2 root is eligible")
    manifest_path = root / "GM12878_run_manifest.json"
    if not manifest_path.is_file():
        raise BundleEligibilityError("canonical GM12878 run manifest is missing")
    manifest = load_json_strict(manifest_path)
    required = {
        "schema": "creditad-gm12878-cell-run-manifest-v2",
        "status": "passed_all_eligibility_gates",
        "cell_line": "GM12878",
        "assembly": "hg19",
        "resolution_bp": RESOLUTION_BP,
        "hic_balance": "ICE_weight",
        "hic_balance_source": "read_only_sidecar",
        "hic_normalization_id": NORMALIZATION_ID,
        "hic_normalization_spec_sha256": ACTIVE_SPEC_IDENTITY["sha256"],
        "full_autosome_validation_required": True,
        "exact_reproducibility_required": True,
        "chip_window_bp": 50_000,
        "tier_rule_version": "max_support_v1",
        "bcp_enabled": False,
    }
    for key, expected in required.items():
        if manifest.get(key) != expected:
            raise BundleEligibilityError("cell manifest %s mismatch" % key)
    if manifest.get("chromosomes") != AUTOSOMES:
        raise BundleEligibilityError("cell manifest must contain ordered autosomes 1..22")

    inputs = manifest.get("inputs", {})
    required_inputs = {
        "mcool",
        "ctcf_bigwig",
        "rad21_bigwig",
        "loops",
        "normalization_spec",
        "ice_sidecar",
        "ice_sidecar_manifest",
        "ice_sidecar_validation",
        "ice_reproducibility",
        "calibration",
        "calibration_manifest",
        "calibration_generator",
        "packaged_calibrated_bands",
    }
    if set(inputs) != required_inputs:
        raise BundleEligibilityError("cell manifest scientific input set mismatch")
    identities = {
        key: _identity_from_metadata(inputs[key], key, require_hash=True)
        for key in sorted(required_inputs)
    }
    if (
        identities["normalization_spec"]["path"] != ACTIVE_SPEC_IDENTITY["path"]
        or identities["normalization_spec"]["sha256"] != ACTIVE_SPEC_IDENTITY["sha256"]
    ):
        raise BundleEligibilityError("normalization spec is not the active artifact")

    code = manifest.get("code", {})
    _selected_identity(code.get("builder", {}), file_identity(BUILDER), "canonical builder")
    _selected_identity(
        code.get("mcool_bed", {}), file_identity(MCOOL_READER_SCRIPT), "mcool reader"
    )
    callers = code.get("builtin_callers_v2", {})
    actual_caller_paths = sorted(
        path for path in (BACKEND / "builtin_callers_v2").rglob("*.py")
        if "__pycache__" not in path.parts
    )
    expected_callers = {str(path.relative_to(BACKEND)) for path in actual_caller_paths}
    if set(callers) != expected_callers:
        raise BundleEligibilityError("builtin caller implementation set mismatch")
    for path in actual_caller_paths:
        relative = str(path.relative_to(BACKEND))
        _selected_identity(callers[relative], file_identity(path), relative)

    outputs = manifest.get("outputs", {})
    annotation_meta = outputs.get("annotation", {})
    annotation_identity = _identity_from_metadata(
        annotation_meta, "canonical annotation", require_hash=True
    )
    method_meta = outputs.get("method_beds", {})
    if set(method_meta) != METHODS:
        raise BundleEligibilityError("five-method BED set mismatch")
    method_identities = {
        method: _identity_from_metadata(metadata, "%s BED" % method)
        for method, metadata in method_meta.items()
    }

    sidecar_manifest = load_json_strict(identities["ice_sidecar_manifest"]["path"])
    convergence = sidecar_manifest.get("autosome_convergence", {})
    if (
        sidecar_manifest.get("schema") != "creditad-ice-sidecar-v1"
        or sidecar_manifest.get("cell_line") != "GM12878"
        or sidecar_manifest.get("normalization_id") != NORMALIZATION_ID
        or sidecar_manifest.get("algorithm") != "cooler.balance_cooler"
        or sidecar_manifest.get("parameters") != PARAMETERS
        or sidecar_manifest.get("acceptance_gates") != ACCEPTANCE_GATES
        or sidecar_manifest.get("software", {}).get("cooler") != REQUIRED_COOLER_VERSION
        or sidecar_manifest.get("balance_semantics") != "multiplicative"
        or sidecar_manifest.get("normalization_spec") != ACTIVE_SPEC_IDENTITY
        or sidecar_manifest.get("output", {}).get("sha256")
        != identities["ice_sidecar"]["sha256"]
        or set(convergence) != set(AUTOSOMES)
        or not all(value is True for value in convergence.values())
    ):
        raise BundleEligibilityError("ICE sidecar manifest is stale or incomplete")
    reproducibility = load_json_strict(identities["ice_reproducibility"]["path"])
    _verify_reproducibility(reproducibility, identities)
    validation = load_json_strict(identities["ice_sidecar_validation"]["path"])
    _verify_validation(validation, identities, sidecar_manifest)
    calibration = load_json_strict(identities["calibration"]["path"])
    calibration_manifest = load_json_strict(identities["calibration_manifest"]["path"])
    packaged = load_json_strict(identities["packaged_calibrated_bands"]["path"])
    _verify_calibration(calibration, calibration_manifest, packaged, identities)

    annotation = pd.read_csv(annotation_identity["path"], sep="\t")
    required_columns = {
        "chrom",
        "pos",
        "votes",
        "hic_balance",
        "hic_balance_source",
        "hic_normalization_id",
        "chip_window_bp",
        "tier_rule_version",
        "data_cell_line",
        "D_tier",
        "assembly",
    }
    if required_columns - set(annotation.columns):
        raise BundleEligibilityError("annotation trust columns are incomplete")
    if annotation[list(required_columns)].isna().any().any():
        raise BundleEligibilityError("annotation contains missing trust values")
    if any(str(column).lower().startswith("bcp") for column in annotation.columns):
        raise BundleEligibilityError("retired BCP columns are forbidden")
    if len(annotation) != int(annotation_meta.get("n_rows", -1)):
        raise BundleEligibilityError("annotation row count differs from manifest")
    normalized_chromosomes = annotation["chrom"].astype(str).str.replace(
        r"^chr", "", regex=True
    )
    if set(normalized_chromosomes) != set(AUTOSOMES):
        raise BundleEligibilityError("annotation chromosome set is incomplete")
    expected_columns = {
        "hic_balance": "ICE_weight",
        "hic_balance_source": "sidecar",
        "hic_normalization_id": NORMALIZATION_ID,
        "chip_window_bp": 50_000,
        "tier_rule_version": "max_support_v1",
        "data_cell_line": "GM12878",
        "assembly": "hg19",
    }
    for column, expected in expected_columns.items():
        if not annotation[column].eq(expected).all():
            raise BundleEligibilityError("annotation %s contract mismatch" % column)

    return {
        "root": root,
        "manifest_path": manifest_path,
        "manifest": manifest,
        "annotation_path": Path(annotation_identity["path"]),
        "method_paths": {
            method: Path(identity["path"])
            for method, identity in method_identities.items()
        },
        "inputs": {key: Path(identity["path"]) for key, identity in identities.items()},
        "sidecar_manifest": sidecar_manifest,
        "validation": validation,
        "reproducibility": reproducibility,
    }
