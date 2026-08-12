from __future__ import annotations

import copy
import json

import pytest

from tad_vci.normalization_spec import (
    ACCEPTANCE_GATES,
    ACTIVE_SPEC,
    ACTIVE_SPEC_IDENTITY,
    NORMALIZATION_ID,
    PARAMETERS,
    RESOLUTION_BP,
    NormalizationSpecError,
    load_spec,
    sha256_file,
)


def test_canonical_ice_v2_contract_is_hash_bound():
    assert NORMALIZATION_ID == "ICE_cis_cooler_0.10.4_v2"
    assert RESOLUTION_BP == 25_000
    assert PARAMETERS["tol"] == 1e-5
    assert PARAMETERS["max_iters"] == 1_000
    assert ACCEPTANCE_GATES["required_autosomes"] == [str(i) for i in range(1, 23)]
    assert ACTIVE_SPEC_IDENTITY["sha256"] == sha256_file(
        ACTIVE_SPEC_IDENTITY["path"]
    )


@pytest.mark.parametrize(
    ("path", "value", "message"),
    [
        (("parameters", "max_iters"), -1, "max_iters"),
        (("normalization_id",), "ICE_cis_cooler_0.10.4_v3", "normalization_id"),
        (("acceptance_gates", "required_autosomes"), ["1", "2"], "autosomes"),
    ],
)
def test_malformed_or_semantically_changed_spec_fails_closed(
    tmp_path, path, value, message
):
    spec = copy.deepcopy(ACTIVE_SPEC)
    target = spec
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    candidate = tmp_path / "normalization_spec.json"
    candidate.write_text(json.dumps(spec))
    with pytest.raises(NormalizationSpecError, match=message):
        load_spec(candidate)


def test_valid_parameter_change_has_a_different_artifact_hash(tmp_path):
    spec = copy.deepcopy(ACTIVE_SPEC)
    spec["parameters"]["max_iters"] = 2_000
    candidate = tmp_path / "normalization_spec.json"
    candidate.write_text(json.dumps(spec))
    assert load_spec(candidate)["parameters"]["max_iters"] == 2_000
    assert sha256_file(candidate) != ACTIVE_SPEC_IDENTITY["sha256"]
