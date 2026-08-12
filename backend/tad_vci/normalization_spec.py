"""Fail-closed loader for the canonical ICE-v2 normalization contract.

The JSON file is package data so analysis scripts, the API eligibility gate,
and an installed wheel all resolve the same inspectable scientific contract.
Generated artifacts bind its exact bytes by SHA-256.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SPEC_PATH = Path(__file__).resolve().with_name("normalization_spec.json")
_TOP_LEVEL_KEYS = {
    "schema",
    "spec_version",
    "normalization_id",
    "required_software",
    "resolution_bp",
    "algorithm",
    "balance_name",
    "balance_semantics",
    "parameters",
    "acceptance_gates",
}
_PARAMETER_KEYS = {
    "cis_only",
    "trans_only",
    "ignore_diags",
    "mad_max",
    "min_nnz",
    "min_count",
    "rescale_marginals",
    "tol",
    "max_iters",
    "chunksize",
}
_GATE_KEYS = {
    "required_autosomes",
    "max_normalized_marginal_cv",
    "min_finite_weight_fraction",
    "max_nonpositive_finite_weights",
    "independent_validation_chromosome",
    "matrix_rtol",
    "matrix_atol",
}


class NormalizationSpecError(ValueError):
    """The normalization contract is missing, malformed, or self-inconsistent."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_exact_keys(mapping: Mapping[str, Any], expected: set[str], label: str) -> None:
    observed = set(mapping)
    if observed != expected:
        raise NormalizationSpecError(
            f"{label} keys differ: missing={sorted(expected - observed)}; "
            f"unexpected={sorted(observed - expected)}"
        )


def validate_spec(spec: Mapping[str, Any]) -> None:
    if not isinstance(spec, Mapping):
        raise NormalizationSpecError("normalization spec must be a JSON object")
    _require_exact_keys(spec, _TOP_LEVEL_KEYS, "top-level normalization spec")
    if spec["schema"] != "creditad-ice-normalization-spec-v1":
        raise NormalizationSpecError("unsupported normalization spec schema")
    if spec["spec_version"] != "v2":
        raise NormalizationSpecError("only the frozen ICE-v2 contract is eligible")
    software = spec["required_software"]
    if not isinstance(software, Mapping) or set(software) != {"cooler"}:
        raise NormalizationSpecError("required_software must name only the Cooler version")
    cooler_version = software["cooler"]
    expected_id = f"ICE_cis_cooler_{cooler_version}_{spec['spec_version']}"
    if spec["normalization_id"] != expected_id:
        raise NormalizationSpecError("normalization_id disagrees with software/version")
    if spec["resolution_bp"] != 25_000:
        raise NormalizationSpecError("the frozen ICE-v2 resolution must be 25 kb")
    if spec["algorithm"] != "cooler.balance_cooler":
        raise NormalizationSpecError("unexpected normalization algorithm")
    if spec["balance_name"] != "ICE_weight":
        raise NormalizationSpecError("unexpected balance name")
    if spec["balance_semantics"] != "multiplicative":
        raise NormalizationSpecError("ICE sidecar weights must be multiplicative")

    parameters = spec["parameters"]
    if not isinstance(parameters, Mapping):
        raise NormalizationSpecError("parameters must be a JSON object")
    _require_exact_keys(parameters, _PARAMETER_KEYS, "normalization parameters")
    if parameters["cis_only"] is not True or parameters["trans_only"] is not False:
        raise NormalizationSpecError("ICE-v2 must be cis-only")
    if parameters["rescale_marginals"] is not True:
        raise NormalizationSpecError("ICE-v2 must rescale marginals")
    for key in ("ignore_diags", "mad_max", "min_nnz", "min_count", "max_iters", "chunksize"):
        if type(parameters[key]) is not int or parameters[key] < 0:
            raise NormalizationSpecError(f"{key} must be a non-negative integer")
    if parameters["ignore_diags"] < 1 or parameters["max_iters"] < 1:
        raise NormalizationSpecError("ignore_diags and max_iters must be positive")
    if not isinstance(parameters["tol"], (int, float)) or not 0 < parameters["tol"] < 1:
        raise NormalizationSpecError("tol must be between zero and one")

    gates = spec["acceptance_gates"]
    if not isinstance(gates, Mapping):
        raise NormalizationSpecError("acceptance_gates must be a JSON object")
    _require_exact_keys(gates, _GATE_KEYS, "acceptance gates")
    autosomes = [str(value) for value in gates["required_autosomes"]]
    if autosomes != [str(value) for value in range(1, 23)]:
        raise NormalizationSpecError("required_autosomes must be ordered 1..22")
    if gates["independent_validation_chromosome"] not in autosomes:
        raise NormalizationSpecError("independent validation chromosome must be an autosome")
    if not 0 < gates["max_normalized_marginal_cv"] <= 1:
        raise NormalizationSpecError("marginal-CV gate must be in (0, 1]")
    if not 0 < gates["min_finite_weight_fraction"] <= 1:
        raise NormalizationSpecError("finite-weight gate must be in (0, 1]")
    if type(gates["max_nonpositive_finite_weights"]) is not int:
        raise NormalizationSpecError("nonpositive-weight gate must be an integer")
    for key in ("matrix_rtol", "matrix_atol"):
        if not isinstance(gates[key], (int, float)) or gates[key] <= 0:
            raise NormalizationSpecError(f"{key} must be positive")


def load_spec(path: Path = SPEC_PATH) -> dict[str, Any]:
    try:
        spec = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise NormalizationSpecError(f"cannot read normalization spec: {path}") from exc
    validate_spec(spec)
    return spec


def spec_identity(path: Path = SPEC_PATH) -> dict[str, Any]:
    path = Path(path).resolve()
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "sha256": sha256_file(path),
    }


ACTIVE_SPEC = load_spec()
ACTIVE_SPEC_IDENTITY = spec_identity()
NORMALIZATION_ID = str(ACTIVE_SPEC["normalization_id"])
SPEC_VERSION = str(ACTIVE_SPEC["spec_version"])
REQUIRED_COOLER_VERSION = str(ACTIVE_SPEC["required_software"]["cooler"])
RESOLUTION_BP = int(ACTIVE_SPEC["resolution_bp"])
PARAMETERS = dict(ACTIVE_SPEC["parameters"])
ACCEPTANCE_GATES = dict(ACTIVE_SPEC["acceptance_gates"])
