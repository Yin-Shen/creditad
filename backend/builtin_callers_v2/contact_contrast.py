"""Clean-room v2 local contact contrast boundary caller.

Identity: custom in-house transparent caller, not an external TAD caller.

Score selection (NOTE):
    score_mode controls the final score column:
        - "normalized" (default): score = raw_contrast / (mean_within + cross_mean + eps)
        - "raw": score = raw_contrast
    allow_negative defaults to False. When allow_negative=True, the score
    is FORCED to raw_contrast regardless of score_mode. This is because
    raw_contrast can take negative values whereas the normalized form is
    bounded to [-1, 1] under non-negative inputs. Users wanting raw
    signed contrast should set allow_negative=True; users wanting the
    bounded form should keep allow_negative=False AND score_mode in
    {"normalized"}. The override is reported in the result warnings when
    triggered.

Non-finite handling:
    Slow oracle drops non-finite pixels per off-diagonal collection.
    Production uses np.isfinite masks. Equivalent at any chunk position
    (no edge-denominator bug). NaN-to-zero imputation upstream MUST NOT
    be applied: mean-based statistics are not invariant to zero
    substitution of missing values.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .common import bin_to_interval, enforce_min_spacing, local_extrema_mask
from .results import BoundaryCallResult
from .thresholds import QuantileThreshold, ThresholdStrategy, fit_threshold, threshold_name
from ._validation import collect_qc_warnings, validate_for_mean_method


@dataclass
class ContactContrastParams:
    window_bins: int = 5
    min_valid_fraction: float = 0.5
    score_mode: str = "normalized"
    threshold_strategy: ThresholdStrategy = QuantileThreshold(0.75)
    min_distance_bins: int = 2
    allow_negative: bool = False
    eps: float = 1e-12


def _offdiag_values(block: np.ndarray) -> np.ndarray:
    if block.size == 0:
        return np.asarray([], dtype=float)
    m = np.asarray(block, dtype=float).copy()
    if m.shape[0] == m.shape[1]:
        mask = ~np.eye(m.shape[0], dtype=bool)
        vals = m[mask]
    else:
        vals = m.ravel()
    return vals[np.isfinite(vals)]


def _finish_contrast(
    left_mean: np.ndarray,
    right_mean: np.ndarray,
    cross_mean: np.ndarray,
    vf_left: np.ndarray,
    vf_right: np.ndarray,
    vf_cross: np.ndarray,
    p: ContactContrastParams,
) -> dict:
    raw = 0.5 * (left_mean + right_mean) - cross_mean
    denom = 0.5 * (left_mean + right_mean) + cross_mean + p.eps
    with np.errstate(divide="ignore", invalid="ignore"):
        norm = raw / denom
    forced_raw = bool(p.allow_negative)
    use_raw = forced_raw or p.score_mode == "raw"
    score = raw if use_raw else norm
    candidate = local_extrema_mask(score, "max", 1) & np.isfinite(score) & (score > 0)
    vals = score[candidate & np.isfinite(score) & (score > 0)]
    cutoff = fit_threshold(p.threshold_strategy, vals)
    return {
        "left_mean": left_mean,
        "right_mean": right_mean,
        "cross_mean": cross_mean,
        "raw_contrast": raw,
        "normalized_contrast": norm,
        "valid_fraction_left": vf_left,
        "valid_fraction_right": vf_right,
        "valid_fraction_cross": vf_cross,
        "candidate_mask": candidate,
        "threshold_cutoff": np.asarray(cutoff),
        "boundary_score": score,
        "score_used": np.asarray("raw" if use_raw else "normalized"),
        "allow_negative_override": np.asarray(forced_raw and p.score_mode != "raw"),
    }


def compute_contact_contrast_slow(
    matrix: np.ndarray, params: ContactContrastParams | None = None
) -> dict:
    p = params or ContactContrastParams()
    arr, preprocess_warnings = validate_for_mean_method(
        matrix,
        allow_negative=p.allow_negative,
        require_symmetric=True,
        window_bins=p.window_bins,
    )
    n = arr.shape[0]
    L = np.full(n, np.nan)
    R = np.full(n, np.nan)
    C = np.full(n, np.nan)
    vfl = np.zeros(n)
    vfr = np.zeros(n)
    vfc = np.zeros(n)
    w = p.window_bins
    for i in range(n):
        if i - w < 0 or i + w > n:
            continue
        left_vals = []
        right_vals = []
        cross_vals = []
        for r in range(i - w, i):
            for c in range(i - w, i):
                if r != c and np.isfinite(arr[r, c]):
                    left_vals.append(arr[r, c])
        for r in range(i, i + w):
            for c in range(i, i + w):
                if r != c and np.isfinite(arr[r, c]):
                    right_vals.append(arr[r, c])
        for r in range(i - w, i):
            for c in range(i, i + w):
                if np.isfinite(arr[r, c]):
                    cross_vals.append(arr[r, c])
        vfl[i] = len(left_vals) / max(w * (w - 1), 1)
        vfr[i] = len(right_vals) / max(w * (w - 1), 1)
        vfc[i] = len(cross_vals) / max(w * w, 1)
        if vfl[i] >= p.min_valid_fraction and left_vals:
            L[i] = np.mean(left_vals)
        if vfr[i] >= p.min_valid_fraction and right_vals:
            R[i] = np.mean(right_vals)
        if vfc[i] >= p.min_valid_fraction and cross_vals:
            C[i] = np.mean(cross_vals)
    result = _finish_contrast(L, R, C, vfl, vfr, vfc, p)
    result["validation_warnings"] = list(preprocess_warnings)
    return result


def compute_contact_contrast_production(
    matrix: np.ndarray, params: ContactContrastParams | None = None
) -> dict:
    p = params or ContactContrastParams()
    arr, preprocess_warnings = validate_for_mean_method(
        matrix,
        allow_negative=p.allow_negative,
        require_symmetric=True,
        window_bins=p.window_bins,
    )
    n = arr.shape[0]
    L = np.full(n, np.nan)
    R = np.full(n, np.nan)
    C = np.full(n, np.nan)
    vfl = np.zeros(n)
    vfr = np.zeros(n)
    vfc = np.zeros(n)
    w = p.window_bins
    for i in range(w, n - w + 1):
        lv = _offdiag_values(arr[i - w : i, i - w : i])
        rv = _offdiag_values(arr[i : i + w, i : i + w])
        cv = arr[i - w : i, i : i + w].ravel()
        cv = cv[np.isfinite(cv)]
        vfl[i] = lv.size / max(w * (w - 1), 1)
        vfr[i] = rv.size / max(w * (w - 1), 1)
        vfc[i] = cv.size / max(w * w, 1)
        if vfl[i] >= p.min_valid_fraction and lv.size:
            L[i] = np.mean(lv)
        if vfr[i] >= p.min_valid_fraction and rv.size:
            R[i] = np.mean(rv)
        if vfc[i] >= p.min_valid_fraction and cv.size:
            C[i] = np.mean(cv)
    result = _finish_contrast(L, R, C, vfl, vfr, vfc, p)
    result["validation_warnings"] = list(preprocess_warnings)
    return result


def call_boundaries(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    region_start: int = 0,
    params: Optional[ContactContrastParams] = None,
) -> BoundaryCallResult:
    p = params or ContactContrastParams()
    qc_warnings = collect_qc_warnings(matrix)
    out = compute_contact_contrast_production(matrix, p)
    n = np.asarray(matrix).shape[0]
    bins = np.arange(n)
    region_end = region_start + n * resolution
    region = f"{chrom}:{region_start}-{region_end}"
    base = (
        out["candidate_mask"]
        & np.isfinite(out["boundary_score"])
        & (out["boundary_score"] >= float(out["threshold_cutoff"]))
    )
    called = bins[base]
    kept = enforce_min_spacing(called, out["boundary_score"][called], p.min_distance_bins)
    threshold_pass = np.zeros(n, dtype=bool)
    threshold_pass[kept] = True
    score_rows = []
    for i in bins:
        iv = bin_to_interval(int(i), resolution, region_start, chrom)
        valid_frac = np.nanmin(
            [
                out["valid_fraction_left"][i],
                out["valid_fraction_right"][i],
                out["valid_fraction_cross"][i],
            ]
        )
        score_rows.append(
            {
                **iv,
                "resolution": resolution,
                "region_start": region_start,
                "region": region,
                "method": "contact_contrast",
                "score": out["boundary_score"][i],
                "boundary_candidate": bool(out["candidate_mask"][i]),
                "threshold_pass": bool(threshold_pass[i]),
                "warnings": "",
                "left_mean": out["left_mean"][i],
                "right_mean": out["right_mean"][i],
                "cross_mean": out["cross_mean"][i],
                "contrast_score": out["boundary_score"][i],
                "valid_fraction": valid_frac,
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
                "boundary_score": out["boundary_score"][i],
                "method": "contact_contrast",
                "rank": rank,
                "threshold": float(out["threshold_cutoff"]),
                "warnings": "",
            }
        )
    warnings = list(qc_warnings) + list(out.get("validation_warnings", []))
    if bool(out["allow_negative_override"]):
        warnings.append("allow_negative=True_forced_score=raw_overriding_score_mode")
    return BoundaryCallResult(
        method_name="contact_contrast",
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
            "identity": "custom in-house",
            "score_used": str(out["score_used"]),
            "original_boundary_semantics": "boundary before bin i",
        },
        warnings=warnings,
        intermediates=out,
    )
