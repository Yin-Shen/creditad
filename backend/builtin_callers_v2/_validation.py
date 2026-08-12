"""Validation helpers for clean-room v2 boundary callers.

This module is referenced by every caller via `from ._validation import ...`.
It implements NaN-preserving input validation:
    * Inf is converted to NaN; the conversion is recorded via the
      `inf_to_nan:<count>` warning string.
    * NaN entries are preserved (NOT zero-imputed) so that downstream
      `np.isfinite` masks and per-window `valid_fraction` gating remain
      correct. Zero-imputation of missing pixels shifts means and would
      break the mean-based callers (insulation, topdom_like,
      contact_contrast) and the sum-based caller (directionality_index).
    * Negative values are rejected when `allow_negative=False` (the default
      for production cooler-balanced/raw matrices).
    * `require_symmetric=True` triggers an average-symmetrization with a
      warning record when the input is not perfectly symmetric within
      tolerance.
    * `window_bins` enforces a minimum matrix size.

The returned `preprocess_warnings` list is appended into
`BoundaryCallResult.warnings` by every caller.
"""

from __future__ import annotations

from typing import Iterable, List, Tuple

import numpy as np

try:
    from scipy import sparse
except Exception:  # pragma: no cover
    sparse = None


def _to_dense_float(matrix) -> np.ndarray:
    if sparse is not None and sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=float)


def _is_approximately_symmetric(arr: np.ndarray, rtol: float = 1e-6, atol: float = 1e-6) -> bool:
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        return False
    finite = np.isfinite(arr) & np.isfinite(arr.T)
    diff = np.abs(arr - arr.T)
    diff[~finite] = 0.0
    if not np.array_equal(np.isnan(arr), np.isnan(arr.T)):
        return False
    scale = np.maximum(np.abs(arr), np.abs(arr.T))
    scale[~finite] = 1.0
    return bool(np.all(diff <= atol + rtol * scale))


def _validate_common(
    matrix,
    *,
    allow_negative: bool,
    require_symmetric: bool,
    window_bins: int | None,
) -> Tuple[np.ndarray, List[str]]:
    arr = _to_dense_float(matrix)
    warnings: List[str] = []
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        raise ValueError("matrix must be a two-dimensional square matrix")
    inf_count = int(np.isinf(arr).sum())
    if inf_count:
        arr = arr.copy()
        arr[np.isinf(arr)] = np.nan
        warnings.append(f"inf_to_nan:{inf_count}")
    nan_count = int(np.isnan(arr).sum())
    if nan_count:
        warnings.append(f"nan_count:{nan_count}")
    finite_mask = np.isfinite(arr)
    if np.any(finite_mask):
        neg_count = int((arr[finite_mask] < 0).sum())
        if neg_count:
            if not allow_negative:
                raise ValueError(
                    f"matrix contains {neg_count} negative finite values but allow_negative=False"
                )
            warnings.append(f"negative_count:{neg_count}")
    if require_symmetric and arr.shape[0] == arr.shape[1] and arr.shape[0] > 0:
        if not _is_approximately_symmetric(arr):
            arr = arr.copy()
            arr_T = arr.T
            both_nan = np.isnan(arr) & np.isnan(arr_T)
            mean = (np.where(np.isnan(arr), 0.0, arr) + np.where(np.isnan(arr_T), 0.0, arr_T)) / 2.0
            mean[both_nan] = np.nan
            arr = mean
            warnings.append("symmetrized_input_was_not_symmetric")
    if window_bins is not None and arr.shape[0] < int(window_bins):
        raise ValueError(
            f"matrix size {arr.shape[0]} smaller than required window_bins {int(window_bins)}"
        )
    finite_zero = np.where(np.isfinite(arr), arr, 0.0)
    if arr.size and np.all(finite_zero == 0):
        warnings.append("all_zero_after_nan_strip")
    zero_rows = int((np.sum(np.abs(finite_zero), axis=1) == 0).sum())
    if zero_rows:
        warnings.append(f"zero_row_count:{zero_rows}")
    return arr, warnings


def validate_for_mean_method(
    matrix,
    *,
    allow_negative: bool = False,
    require_symmetric: bool = True,
    window_bins: int | None = None,
) -> Tuple[np.ndarray, List[str]]:
    """Validate a matrix for mean-based callers (insulation, topdom_like, contact_contrast, spectral_profile)."""

    arr, warnings = _validate_common(
        matrix,
        allow_negative=allow_negative,
        require_symmetric=require_symmetric,
        window_bins=window_bins,
    )
    return arr, warnings


def validate_for_sum_method(
    matrix,
    *,
    allow_negative: bool = False,
    require_symmetric: bool = True,
    window_bins: int | None = None,
) -> Tuple[np.ndarray, List[str]]:
    """Validate a matrix for sum-based callers (directionality_index)."""

    arr, warnings = _validate_common(
        matrix,
        allow_negative=allow_negative,
        require_symmetric=require_symmetric,
        window_bins=window_bins,
    )
    return arr, warnings


def collect_qc_warnings(matrix) -> List[str]:
    """Return a JSON-serializable list of QC warning strings about a contact matrix."""

    arr = _to_dense_float(matrix)
    warnings: List[str] = []
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        warnings.append("not_square_or_not_2d")
        return warnings
    if arr.size == 0:
        warnings.append("empty_matrix")
        return warnings
    nan_count = int(np.isnan(arr).sum())
    inf_count = int(np.isinf(arr).sum())
    if nan_count:
        warnings.append(f"nan_count:{nan_count}")
    if inf_count:
        warnings.append(f"inf_count:{inf_count}")
    finite_zero = np.where(np.isfinite(arr), arr, 0.0)
    if np.all(finite_zero == 0):
        warnings.append("all_zero_matrix")
    zero_rows = int((np.sum(np.abs(finite_zero), axis=1) == 0).sum())
    if zero_rows:
        warnings.append(f"zero_row_count:{zero_rows}")
    finite = np.isfinite(arr)
    if np.any(finite):
        finite_vals = arr[finite]
        neg = int((finite_vals < 0).sum())
        if neg:
            warnings.append(f"negative_count:{neg}")
    sparsity = float(np.sum(arr == 0) / arr.size)
    if sparsity > 0.90:
        warnings.append(f"high_sparsity:{sparsity:.3f}")
    return warnings


__all__ = [
    "validate_for_mean_method",
    "validate_for_sum_method",
    "collect_qc_warnings",
]
