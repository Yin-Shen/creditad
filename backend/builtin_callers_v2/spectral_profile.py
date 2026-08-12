"""Clean-room v2 spectral boundary profile caller.

Identity: custom/simplified spectral profile. This is NOT exact SpectralTAD.

Implementation notes:
    Production currently calls slow implementation. No optimized path
    exists. Manuscript and runtime tables MUST NOT claim an optimized
    production for this method.

Per-window quantities:
    For each sliding window starting at index s of size window_bins, we
    compute a normalized graph Laplacian, take eigenvectors at positions
    1..eigen_k (skipping the trivial 0-th), row-normalize them to unit
    sphere, then compute consecutive-row jumps in this embedding. These
    jumps are written to per-bin score lists keyed by global bin index
    s + j (for j = 1, ..., window_bins - 1).

eigen_gap:
    For each window with at least eigen_k + 2 finite eigenvalues, we
    compute gap_s = vals[eigen_k + 1] - vals[eigen_k] (the gap between
    the last selected eigenvalue and the first non-selected one). Per-bin
    eigen_gap is the mean across windows that cover that bin.

fiedler_jump:
    Computed using ONLY the first non-trivial eigenvector (the Fiedler
    vector, vecs[:, 1]). Distinct from spectral_score, which uses all
    eigen_k selected eigenvectors. Aggregated per bin using the same
    rule as spectral_score (max or mean).

quality_flag (per bin):
    Aggregated as the worst severity flag among all windows covering the
    bin. Severity order: ok < constant < bad_degree < too_small.

Tail coverage:
    The sliding-window loop appends a final window at s = n - window_bins
    if not already in the start list, so the chunk tail is covered.

Non-finite handling:
    Inputs containing non-finite values are flagged in qflags as
    constant or bad_degree depending on which condition is hit during
    Laplacian construction.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .common import bin_to_interval, enforce_min_spacing, local_extrema_mask, symmetrize_matrix
from .results import BoundaryCallResult
from .thresholds import QuantileThreshold, ThresholdStrategy, fit_threshold, threshold_name
from ._validation import collect_qc_warnings, validate_for_mean_method


_FLAG_SEVERITY = {"ok": 0, "constant": 1, "bad_degree": 2, "too_small": 3}


@dataclass
class SpectralProfileParams:
    window_bins: int = 20
    eigen_k: int = 2
    min_size_bins: int = 5
    step_bins: Optional[int] = None
    regularization_eps: float = 1e-9
    min_degree: float = 0.0
    threshold_strategy: ThresholdStrategy = QuantileThreshold(0.75)
    min_distance_bins: int = 2
    zero_diagonal: bool = True
    aggregate: str = "max"
    allow_negative: bool = False


def _normalized_laplacian(adjacency: np.ndarray, eps: float = 1e-9) -> tuple[np.ndarray, np.ndarray]:
    A = np.asarray(adjacency, dtype=float)
    degree = A.sum(axis=1)
    inv = 1.0 / np.sqrt(degree + eps)
    L = np.eye(A.shape[0]) - (inv[:, None] * A * inv[None, :])
    return degree, L


def _align_eigenvectors(vecs: np.ndarray) -> np.ndarray:
    out = vecs.copy()
    for k in range(out.shape[1]):
        col = out[:, k]
        idx = int(np.argmax(np.abs(col)))
        if col[idx] < 0:
            out[:, k] *= -1
    return out


def _flag_suffix(flag: str) -> str:
    if flag.endswith("_too_small") or flag == "too_small":
        return "too_small"
    if flag.endswith("_bad_degree") or flag == "bad_degree":
        return "bad_degree"
    if flag.endswith("_constant") or flag == "constant":
        return "constant"
    if flag.endswith("_ok") or flag == "ok":
        return "ok"
    return "ok"


def _aggregate_bin_quality(
    n: int,
    window_starts: list[int],
    window_size: int,
    qflags: list[str],
) -> list[str]:
    per_bin = [""] * n
    per_bin_severity = np.full(n, -1, dtype=int)
    for w_idx, s in enumerate(window_starts):
        if w_idx >= len(qflags):
            break
        flag = qflags[w_idx]
        suffix = _flag_suffix(flag)
        severity = _FLAG_SEVERITY.get(suffix, 0)
        end = min(n, s + window_size)
        for bin_idx in range(s, end):
            if severity > per_bin_severity[bin_idx]:
                per_bin_severity[bin_idx] = severity
                per_bin[bin_idx] = flag
    return per_bin


def _finish_spectral(
    score_lists: list[list[float]],
    fiedler_lists: list[list[float]],
    eigen_gap_lists: list[list[float]],
    n: int,
    p: SpectralProfileParams,
    eigenvalues_list: list[np.ndarray],
    eigenvectors_list: list[np.ndarray],
    jumps_list: list[np.ndarray],
    degree_list: list[np.ndarray],
    qflags: list[str],
    window_starts: list[int],
) -> dict:
    spectral = np.full(n, np.nan)
    fiedler = np.full(n, np.nan)
    eigen_gap = np.full(n, np.nan)
    for i in range(n):
        s_vals = [v for v in score_lists[i] if np.isfinite(v)]
        if s_vals:
            spectral[i] = max(s_vals) if p.aggregate == "max" else float(np.mean(s_vals))
        f_vals = [v for v in fiedler_lists[i] if np.isfinite(v)]
        if f_vals:
            fiedler[i] = max(f_vals) if p.aggregate == "max" else float(np.mean(f_vals))
        g_vals = [v for v in eigen_gap_lists[i] if np.isfinite(v)]
        if g_vals:
            eigen_gap[i] = float(np.mean(g_vals))
    candidate = local_extrema_mask(spectral, "max", 1) & np.isfinite(spectral) & (spectral > 0)
    vals = spectral[candidate]
    cutoff = fit_threshold(p.threshold_strategy, vals)
    per_bin_quality = _aggregate_bin_quality(n, window_starts, p.window_bins, qflags)
    return {
        "degree": np.asarray(degree_list, dtype=object),
        "eigenvalues": np.asarray(eigenvalues_list, dtype=object),
        "selected_eigenvectors": np.asarray(eigenvectors_list, dtype=object),
        "local_jump_per_window": np.asarray(jumps_list, dtype=object),
        "spectral_score": spectral,
        "fiedler_jump": fiedler,
        "eigen_gap": eigen_gap,
        "quality_flags": np.asarray(qflags, dtype=object),
        "per_bin_quality_flag": np.asarray(per_bin_quality, dtype=object),
        "window_starts": np.asarray(window_starts, dtype=int),
        "candidate_mask": candidate,
        "threshold_cutoff": np.asarray(cutoff),
    }


def _empty_result(n: int, p: SpectralProfileParams, reason: str) -> dict:
    result = _finish_spectral(
        score_lists=[[] for _ in range(n)],
        fiedler_lists=[[] for _ in range(n)],
        eigen_gap_lists=[[] for _ in range(n)],
        n=n,
        p=p,
        eigenvalues_list=[],
        eigenvectors_list=[],
        jumps_list=[],
        degree_list=[],
        qflags=[reason],
        window_starts=[],
    )
    result["per_bin_quality_flag"] = np.asarray([reason] * n, dtype=object)
    return result


def compute_spectral_profile_slow(matrix: np.ndarray, params: SpectralProfileParams | None = None) -> dict:
    p = params or SpectralProfileParams()
    step = p.step_bins or max(1, p.window_bins // 2)
    canonical, preprocess_warnings = validate_for_mean_method(
        matrix,
        allow_negative=p.allow_negative,
        require_symmetric=True,
        window_bins=p.min_size_bins,
    )
    arr = symmetrize_matrix(canonical)
    n = arr.shape[0]
    if p.window_bins < p.min_size_bins or n < p.min_size_bins or n < p.window_bins:
        result = _empty_result(n, p, "too_small")
        result["validation_warnings"] = list(preprocess_warnings)
        return result
    starts = list(range(0, n - p.window_bins + 1, step))
    if not starts:
        starts = [0]
    last_required = n - p.window_bins
    if last_required >= 0 and starts[-1] != last_required:
        starts.append(last_required)
    score_lists: list[list[float]] = [[] for _ in range(n)]
    fiedler_lists: list[list[float]] = [[] for _ in range(n)]
    eigen_gap_lists: list[list[float]] = [[] for _ in range(n)]
    eigenvalues_list: list[np.ndarray] = []
    eigenvectors_list: list[np.ndarray] = []
    jumps_list: list[np.ndarray] = []
    degree_list: list[np.ndarray] = []
    qflags: list[str] = []
    window_starts_used: list[int] = []
    for s in starts:
        window_starts_used.append(s)
        sub = symmetrize_matrix(arr[s : s + p.window_bins, s : s + p.window_bins])
        if sub.shape[0] < p.min_size_bins:
            qflags.append(f"window_{s}_too_small")
            degree_list.append(np.asarray([], dtype=float))
            eigenvalues_list.append(np.asarray([], dtype=float))
            eigenvectors_list.append(np.asarray([], dtype=float))
            jumps_list.append(np.asarray([], dtype=float))
            continue
        if p.zero_diagonal:
            sub = sub.copy()
            np.fill_diagonal(sub, 0.0)
        offdiag = sub[~np.eye(sub.shape[0], dtype=bool)]
        if offdiag.size == 0 or np.nanstd(offdiag) <= p.regularization_eps:
            qflags.append(f"window_{s}_constant")
            degree_list.append(sub.sum(axis=1))
            eigenvalues_list.append(np.asarray([], dtype=float))
            eigenvectors_list.append(np.asarray([], dtype=float))
            jumps_list.append(np.asarray([], dtype=float))
            continue
        degree, lap = _normalized_laplacian(sub, p.regularization_eps)
        degree_list.append(degree)
        if np.any(degree <= p.min_degree) or not np.isfinite(lap).all():
            qflags.append(f"window_{s}_bad_degree")
            eigenvalues_list.append(np.asarray([], dtype=float))
            eigenvectors_list.append(np.asarray([], dtype=float))
            jumps_list.append(np.asarray([], dtype=float))
            continue
        vals, vecs = np.linalg.eigh(lap)
        order = np.argsort(vals)
        vals = vals[order]
        vecs = vecs[:, order]
        selected = vecs[:, 1 : 1 + p.eigen_k]
        selected = _align_eigenvectors(selected)
        norms = np.linalg.norm(selected, axis=1)
        Z = np.divide(selected, norms[:, None], out=np.zeros_like(selected), where=norms[:, None] > 0)
        fiedler_col = Z[:, 0] if Z.shape[1] >= 1 else np.zeros(sub.shape[0], dtype=float)
        if 1 + p.eigen_k < vals.size and np.isfinite(vals[p.eigen_k]) and np.isfinite(vals[1 + p.eigen_k]):
            window_gap = float(vals[1 + p.eigen_k] - vals[p.eigen_k])
        else:
            window_gap = np.nan
        jumps = np.full(sub.shape[0], np.nan)
        for j in range(1, sub.shape[0]):
            jumps[j] = np.linalg.norm(Z[j] - Z[j - 1])
            fiedler_jump = float(abs(fiedler_col[j] - fiedler_col[j - 1]))
            global_bin = s + j
            if global_bin < n:
                score_lists[global_bin].append(float(jumps[j]))
                fiedler_lists[global_bin].append(fiedler_jump)
        if np.isfinite(window_gap):
            end = min(n, s + p.window_bins)
            for bin_idx in range(s, end):
                eigen_gap_lists[bin_idx].append(window_gap)
        eigenvalues_list.append(vals)
        eigenvectors_list.append(selected)
        jumps_list.append(jumps)
        qflags.append(f"window_{s}_ok")
    result = _finish_spectral(
        score_lists=score_lists,
        fiedler_lists=fiedler_lists,
        eigen_gap_lists=eigen_gap_lists,
        n=n,
        p=p,
        eigenvalues_list=eigenvalues_list,
        eigenvectors_list=eigenvectors_list,
        jumps_list=jumps_list,
        degree_list=degree_list,
        qflags=qflags,
        window_starts=window_starts_used,
    )
    result["validation_warnings"] = list(preprocess_warnings)
    return result


def compute_spectral_profile_production(matrix: np.ndarray, params: SpectralProfileParams | None = None) -> dict:
    return compute_spectral_profile_slow(matrix, params)


def call_boundaries(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    region_start: int = 0,
    params: Optional[SpectralProfileParams] = None,
) -> BoundaryCallResult:
    p = params or SpectralProfileParams()
    qc_warnings = collect_qc_warnings(matrix)
    out = compute_spectral_profile_production(matrix, p)
    n = np.asarray(matrix).shape[0]
    bins = np.arange(n)
    region_end = region_start + n * resolution
    region = f"{chrom}:{region_start}-{region_end}"
    base = (
        out["candidate_mask"]
        & np.isfinite(out["spectral_score"])
        & (out["spectral_score"] >= float(out["threshold_cutoff"]))
    )
    called = bins[base]
    kept = enforce_min_spacing(called, out["spectral_score"][called], p.min_distance_bins)
    threshold_pass = np.zeros(n, dtype=bool)
    threshold_pass[kept] = True
    score_rows = []
    per_bin_flag = out["per_bin_quality_flag"]
    for i in bins:
        iv = bin_to_interval(int(i), resolution, region_start, chrom)
        score_rows.append(
            {
                **iv,
                "resolution": resolution,
                "region_start": region_start,
                "region": region,
                "method": "spectral_profile",
                "score": out["spectral_score"][i],
                "boundary_candidate": bool(out["candidate_mask"][i]),
                "threshold_pass": bool(threshold_pass[i]),
                "warnings": "",
                "spectral_score": out["spectral_score"][i],
                "eigen_gap": out["eigen_gap"][i],
                "fiedler_jump": out["fiedler_jump"][i],
                "quality_flag": str(per_bin_flag[i]) if i < len(per_bin_flag) else "",
            }
        )
    b_rows = []
    for rank, i in enumerate(kept, start=1):
        iv = bin_to_interval(int(i), resolution, region_start, chrom)
        b_rows.append(
            {
                **iv,
                "resolution": resolution,
                "region_start": region_start,
                "region": region,
                "boundary_score": out["spectral_score"][i],
                "method": "spectral_profile",
                "rank": rank,
                "threshold": float(out["threshold_cutoff"]),
                "warnings": "",
            }
        )
    warnings = [
        "SpectralProfile is not exact SpectralTAD",
        "production calls slow implementation: no optimized path",
        *qc_warnings,
        *list(out.get("validation_warnings", [])),
    ]
    return BoundaryCallResult(
        method_name="spectral_profile",
        chrom=chrom,
        resolution=resolution,
        region_start=region_start,
        region_end=region_end,
        local_bin_indices=bins,
        boundaries=pd.DataFrame(b_rows),
        scores=pd.DataFrame(score_rows),
        parameters={
            **asdict(p),
            "threshold_strategy": threshold_name(p.threshold_strategy),
            "threshold_kind": "per_chunk_quantile",
            "identity": "custom/simplified spectral profile, not SpectralTAD",
            "original_boundary_semantics": "boundary before bin i",
        },
        warnings=warnings,
        intermediates=out,
    )
