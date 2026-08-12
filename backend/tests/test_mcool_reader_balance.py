"""Normalization contract for the scientific mcool reader."""
from __future__ import annotations

import shutil

import cooler
import h5py
import numpy as np
import pytest

from mcool_bed import McoolReader


def _copy_with_balance(tmp_path, source, name: str, values: np.ndarray):
    target = tmp_path / f"toy_{name}.mcool"
    shutil.copyfile(source, target)
    with h5py.File(target, "r+") as h5:
        bins = h5["resolutions/25000/bins"]
        if name in bins:
            del bins[name]
        bins.create_dataset(name, data=np.asarray(values, dtype=float))
    return target


@pytest.mark.parametrize("name", ["KR", "VC", "VC_SQRT"])
def test_4dn_divisive_balance_matches_cooler(tmp_path, name):
    source = "backend/tests/builtin_callers/toy_8mb.mcool"
    n = cooler.Cooler(f"{source}::resolutions/25000").shape[0]
    values = np.linspace(0.75, 1.25, n)
    target = _copy_with_balance(tmp_path, source, name, values)

    expected = cooler.Cooler(f"{target}::resolutions/25000").matrix(
        balance=name
    ).fetch("chr_toy")
    observed = McoolReader(str(target), 25000, chroms=["chr_toy"],
                           balance_name=name, strict_balance=True).matrix("chr_toy")
    np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12,
                               equal_nan=True)


def test_multiplicative_weight_matches_cooler(tmp_path):
    source = "backend/tests/builtin_callers/toy_8mb.mcool"
    n = cooler.Cooler(f"{source}::resolutions/25000").shape[0]
    values = np.linspace(0.8, 1.2, n)
    target = _copy_with_balance(tmp_path, source, "weight", values)

    expected = cooler.Cooler(f"{target}::resolutions/25000").matrix(
        balance="weight"
    ).fetch("chr_toy")
    observed = McoolReader(str(target), 25000, chroms=["chr_toy"],
                           balance_name="weight", strict_balance=True).matrix("chr_toy")
    np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12,
                               equal_nan=True)


def test_strict_balance_rejects_missing_column():
    source = "backend/tests/builtin_callers/toy_8mb.mcool"
    reader = McoolReader(source, 25000, chroms=["chr_toy"],
                         balance_name="KR", strict_balance=True)
    with pytest.raises(ValueError, match="requested balance column bins/KR is missing"):
        reader.matrix("chr_toy")


def _write_sidecar(
    path,
    weights,
    resolution=25000,
    semantics="multiplicative",
    balance_name="ICE_weight",
):
    np.savez_compressed(
        path,
        weights=np.asarray(weights, dtype=float),
        resolution_bp=np.asarray(resolution, dtype=np.int64),
        balance_name=np.asarray(balance_name),
        balance_semantics=np.asarray(semantics),
    )


def test_external_ice_sidecar_matches_independent_raw_multiplication(tmp_path):
    source = "backend/tests/builtin_callers/toy_8mb.mcool"
    clr = cooler.Cooler(f"{source}::resolutions/25000")
    values = np.linspace(0.8, 1.2, clr.shape[0])
    values[3] = np.nan
    sidecar = tmp_path / "weights.npz"
    _write_sidecar(sidecar, values)

    raw = clr.matrix(balance=False).fetch("chr_toy").astype(float)
    expected = raw * values[:, None] * values[None, :]
    observed = McoolReader(
        source,
        25000,
        chroms=["chr_toy"],
        balance_name="ICE_weight",
        strict_balance=True,
        balance_weights_path=str(sidecar),
    ).matrix("chr_toy")
    np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12,
                               equal_nan=True)


@pytest.mark.parametrize(
    ("values_delta", "resolution", "semantics", "match"),
    [
        (-1, 25000, "multiplicative", "weights shape"),
        (0, 50000, "multiplicative", "resolution"),
        (0, 25000, "divisive", "semantics must be multiplicative"),
    ],
)
def test_external_ice_sidecar_rejects_incompatible_metadata(
    tmp_path, values_delta, resolution, semantics, match
):
    source = "backend/tests/builtin_callers/toy_8mb.mcool"
    n = cooler.Cooler(f"{source}::resolutions/25000").shape[0] + values_delta
    sidecar = tmp_path / "bad_weights.npz"
    _write_sidecar(sidecar, np.ones(n), resolution=resolution, semantics=semantics)
    with pytest.raises(ValueError, match=match):
        McoolReader(
            source,
            25000,
            balance_name="ICE_weight",
            strict_balance=True,
            balance_weights_path=str(sidecar),
        )


@pytest.mark.parametrize(
    ("bad_value", "match"),
    [(np.inf, "infinite weights"), (-np.inf, "infinite weights")],
)
def test_external_ice_sidecar_rejects_infinite_weights(tmp_path, bad_value, match):
    source = "backend/tests/builtin_callers/toy_8mb.mcool"
    n = cooler.Cooler(f"{source}::resolutions/25000").shape[0]
    values = np.ones(n)
    values[3] = bad_value
    sidecar = tmp_path / "infinite_weights.npz"
    _write_sidecar(sidecar, values)
    with pytest.raises(ValueError, match=match):
        McoolReader(
            source,
            25000,
            balance_name="ICE_weight",
            strict_balance=True,
            balance_weights_path=str(sidecar),
        )


def test_external_ice_sidecar_rejects_mismatched_balance_name(tmp_path):
    source = "backend/tests/builtin_callers/toy_8mb.mcool"
    n = cooler.Cooler(f"{source}::resolutions/25000").shape[0]
    sidecar = tmp_path / "wrong_name.npz"
    _write_sidecar(sidecar, np.ones(n), balance_name="weight")
    with pytest.raises(ValueError, match="name must be ICE_weight"):
        McoolReader(
            source,
            25000,
            balance_name="ICE_weight",
            strict_balance=True,
            balance_weights_path=str(sidecar),
        )
