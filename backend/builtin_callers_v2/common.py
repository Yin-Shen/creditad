"""Common utilities for clean-room v2 TAD boundary callers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

try:  # optional; tests and scripts can pass dense arrays directly
    from scipy import sparse
except Exception:  # pragma: no cover
    sparse = None


def _to_dense_float(matrix) -> np.ndarray:
    if sparse is not None and sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=float)


def compute_matrix_qc(
    matrix,
    *,
    chrom: str | None = None,
    region: str | None = None,
    resolution: int | None = None,
    balanced_or_raw: str | None = None,
    balance_available: bool | None = None,
    selected_reason: str | None = None,
    file_path: str | None = None,
    selected_uri: str | None = None,
    dataset_id: str | None = None,
    cell_line: str | None = None,
    organism: str | None = None,
    level: str | None = None,
) -> dict:
    """Return JSON-serializable QC for a contact matrix."""

    arr = _to_dense_float(matrix)
    finite = np.isfinite(arr)
    n_total = int(arr.size)
    finite_vals = arr[finite]
    is_square = arr.ndim == 2 and arr.shape[0] == arr.shape[1]
    diag = np.diag(arr) if is_square else np.array([], dtype=float)
    offdiag = arr.copy() if arr.ndim == 2 else np.array([], dtype=float)
    if is_square:
        np.fill_diagonal(offdiag, 0.0)
    zero_rows = []
    zero_cols = []
    if arr.ndim == 2:
        finite_zero = np.where(np.isfinite(arr), arr, 0.0)
        zero_rows = np.where(np.sum(np.abs(finite_zero), axis=1) == 0)[0].tolist()
        zero_cols = np.where(np.sum(np.abs(finite_zero), axis=0) == 0)[0].tolist()
    qc = {
        "dataset_id": dataset_id,
        "cell_line": cell_line,
        "organism": organism,
        "level": level,
        "file_path": file_path,
        "selected_uri": selected_uri,
        "matrix_shape": list(arr.shape),
        "chrom": chrom,
        "region": region,
        "resolution": resolution,
        "balanced_or_raw": balanced_or_raw,
        "balance_available": balance_available,
        "nan_count": int(np.isnan(arr).sum()),
        "inf_count": int(np.isinf(arr).sum()),
        "zero_row_count": int(len(zero_rows)),
        "zero_col_count": int(len(zero_cols)),
        "negative_count": int((arr < 0).sum()) if arr.size else 0,
        "total_sum": float(np.nansum(arr)) if arr.size else 0.0,
        "sparsity": float(np.sum(arr == 0) / n_total) if n_total else 0.0,
        "min": float(np.nanmin(arr)) if finite_vals.size else None,
        "max": float(np.nanmax(arr)) if finite_vals.size else None,
        "mean": float(np.nanmean(arr)) if finite_vals.size else None,
        "median": float(np.nanmedian(arr)) if finite_vals.size else None,
        "diag_sum": float(np.nansum(diag)) if diag.size else 0.0,
        "offdiag_sum": float(np.nansum(offdiag)) if offdiag.size else 0.0,
        "all_zero": bool(arr.size > 0 and np.all(np.where(np.isfinite(arr), arr, 0.0) == 0)),
        "is_square": bool(is_square),
        "is_symmetric": bool(is_square and np.allclose(arr, arr.T, equal_nan=True)),
        "selected_reason": selected_reason,
        "warnings": [],
    }
    if qc["nan_count"]:
        qc["warnings"].append("NaN values present")
    if qc["inf_count"]:
        qc["warnings"].append("Inf values present")
    if qc["negative_count"]:
        qc["warnings"].append("negative values present")
    if qc["zero_row_count"]:
        qc["warnings"].append("zero rows present")
    if qc["zero_col_count"]:
        qc["warnings"].append("zero columns present")
    if qc["all_zero"]:
        qc["warnings"].append("all-zero matrix")
    if qc["sparsity"] > 0.90:
        qc["warnings"].append("high sparsity")
    if is_square and np.nansum(diag) == 0:
        qc["warnings"].append("zero diagonal sum")
    return qc


def validate_matrix(
    matrix,
    *,
    allow_negative: bool = False,
    require_symmetric: bool = True,
    window_bins: int | None = None,
    allow_nan: bool = False,
    allow_inf: bool = False,
) -> dict:
    """Validate a matrix and return QC metadata.

    The function rejects structural errors and unsafe numeric values by default.
    All-zero, sparse, and zero-row matrices are recorded as warnings, not fatal
    errors, so negative controls can be processed and shown to produce no calls.
    """

    arr = _to_dense_float(matrix)
    qc = compute_matrix_qc(arr)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        raise ValueError("matrix must be a two-dimensional square matrix")
    if qc["nan_count"] and not allow_nan:
        raise ValueError("matrix contains NaN values")
    if qc["inf_count"] and not allow_inf:
        raise ValueError("matrix contains Inf values")
    if require_symmetric and not qc["is_symmetric"]:
        raise ValueError("matrix must be symmetric")
    if qc["negative_count"] and not allow_negative:
        raise ValueError("matrix contains negative values but allow_negative=False")
    if window_bins is not None and arr.shape[0] < int(window_bins):
        raise ValueError("matrix is smaller than required window")
    return qc


def symmetrize_matrix(matrix, mode: str = "average") -> np.ndarray:
    arr = _to_dense_float(matrix)
    if mode == "average":
        return (arr + arr.T) / 2.0
    if mode == "upper":
        return np.triu(arr) + np.triu(arr, 1).T
    if mode == "lower":
        return np.tril(arr) + np.tril(arr, -1).T
    raise ValueError(f"unknown symmetrization mode: {mode}")


def safe_nanmean(values: Iterable[float], min_valid_fraction: float = 0.0) -> float:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0:
        return float("nan")
    valid = np.isfinite(arr)
    if valid.sum() / arr.size < min_valid_fraction:
        return float("nan")
    if valid.sum() == 0:
        return float("nan")
    return float(np.nanmean(arr))


def bin_to_interval(
    local_bin_index: int,
    resolution: int,
    region_start: int,
    chrom: str | None = None,
) -> dict:
    start = int(region_start) + int(local_bin_index) * int(resolution)
    end = start + int(resolution)
    out = {"local_bin_index": int(local_bin_index), "start": start, "end": end}
    if chrom is not None:
        out = {"chrom": chrom, **out}
    return out


def enforce_min_spacing(boundaries: Sequence[int], scores: Sequence[float], min_distance: int) -> np.ndarray:
    b = np.asarray(boundaries, dtype=int)
    s = np.asarray(scores, dtype=float)
    if b.size == 0:
        return b
    if min_distance <= 0:
        return np.sort(b)
    order = np.lexsort((b, -np.nan_to_num(s, nan=-np.inf)))
    kept: list[int] = []
    for idx in order:
        cand = int(b[idx])
        if all(abs(cand - old) >= min_distance for old in kept):
            kept.append(cand)
    return np.asarray(sorted(kept), dtype=int)


def one_to_one_match(pred: Sequence[int], ref: Sequence[int], tolerance_bp: int) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int, int]] = []
    for p in pred:
        for r in ref:
            d = abs(int(p) - int(r))
            if d <= tolerance_bp:
                pairs.append((d, int(p), int(r)))
    pairs.sort()
    used_p: set[int] = set()
    used_r: set[int] = set()
    matches: list[tuple[int, int]] = []
    for _, p, r in pairs:
        if p in used_p or r in used_r:
            continue
        used_p.add(p)
        used_r.add(r)
        matches.append((p, r))
    return sorted(matches)


def precision_recall_f1_jaccard(matches, pred_count: int, ref_count: int) -> dict:
    matched = len(matches)
    precision = matched / pred_count if pred_count else 0.0
    recall = matched / ref_count if ref_count else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    union = pred_count + ref_count - matched
    jaccard = matched / union if union else 0.0
    return {
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "jaccard": float(jaccard),
        "matched_count": int(matched),
        "false_positive_count": int(pred_count - matched),
        "false_negative_count": int(ref_count - matched),
    }


def make_block_matrix(
    block_sizes: Sequence[int],
    *,
    within: float = 10.0,
    between: float = 1.0,
    diag: float = 0.0,
) -> tuple[np.ndarray, list[int]]:
    n = int(sum(block_sizes))
    mat = np.full((n, n), float(between), dtype=float)
    truth: list[int] = []
    start = 0
    for size in block_sizes:
        end = start + int(size)
        mat[start:end, start:end] = float(within)
        start = end
        if start < n:
            truth.append(start)
    np.fill_diagonal(mat, float(diag))
    return mat, truth


def add_distance_decay(matrix, strength: float = 1.0, power: float = 1.0) -> np.ndarray:
    arr = _to_dense_float(matrix).copy()
    n = arr.shape[0]
    idx = np.arange(n)
    dist = np.abs(idx[:, None] - idx[None, :])
    decay = strength / np.power(dist + 1.0, power)
    return arr + decay


def add_noise(matrix, scale: float = 0.1, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    arr = _to_dense_float(matrix).copy()
    noise = rng.normal(0.0, scale, size=arr.shape)
    noise = (noise + noise.T) / 2.0
    out = arr + noise
    out[out < 0] = 0.0
    return out


def add_loops(matrix, loops: Sequence[tuple[int, int, float]]) -> np.ndarray:
    arr = _to_dense_float(matrix).copy()
    for i, j, value in loops:
        arr[int(i), int(j)] += float(value)
        arr[int(j), int(i)] += float(value)
    return arr


def add_compartment_pattern(matrix, amplitude: float = 1.0, period: int = 10) -> np.ndarray:
    arr = _to_dense_float(matrix).copy()
    n = arr.shape[0]
    x = np.sin(np.arange(n) / max(period, 1) * 2 * np.pi)
    plaid = np.outer(x, x) * amplitude
    out = arr + plaid
    out[out < 0] = 0.0
    return (out + out.T) / 2.0


def save_matrix_qc_json(qc: dict, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(qc, indent=2, sort_keys=True))


def local_extrema_mask(values: Sequence[float], mode: str, window: int = 1) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    mask = np.zeros(arr.shape, dtype=bool)
    for i, val in enumerate(arr):
        if not np.isfinite(val):
            continue
        lo = max(0, i - window)
        hi = min(arr.size, i + window + 1)
        neigh = arr[lo:hi]
        finite = neigh[np.isfinite(neigh)]
        if finite.size <= 1:
            continue
        if mode == "min":
            mask[i] = bool(val <= np.min(finite) and np.any(val < finite))
        elif mode == "max":
            mask[i] = bool(val >= np.max(finite) and np.any(val > finite))
        else:
            raise ValueError("mode must be min or max")
    return mask
