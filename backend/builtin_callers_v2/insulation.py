"""Clean-room v2 insulation score boundary caller.

Identity: cooltools-compatible target. This module implements explicit diamond
pixel enumeration, local-minimum candidate generation, and valley-drop
prominence. It does not call cooltools or any external TAD caller.

Default ignore_diags is 0, which differs from cooltools' typical default of 2.
Override via InsulationParams.ignore_diags when comparing against cooltools.

valley_drop semantics (intermediate key `boundary_strength`):
    For a local minimum at bin i, valley_drop = min(L, R) where L and R are
    the maximum scores on the left and right side of i. If
    prominence_window_bins is None, L/R use the FULL upstream/downstream
    arrays of the chunk (global definition; chunk-size sensitive). If
    prominence_window_bins is a positive integer, L/R use only that many
    bins on each side (local definition; recommended for cross-chunk
    comparability). This metric is intentionally NOT called peak
    prominence in the topographic sense.

bad_bin handling:
    A bin is flagged bad if its diamond has valid_pixel_fraction below
    min_frac_valid_pixels or its raw_insulation is non-finite. If
    min_dist_bad_bin > 0, every bin within that distance of a bad bin is
    excluded from boundary calling.

Non-finite handling:
    Slow oracle drops non-finite values from each diamond pixel set.
    Production masks non-finite via np.isfinite. Both produce identical
    raw_insulation values when valid pixel counts match. NaN-to-zero
    imputation upstream (e.g. in genome-wide runner) MUST NOT be applied:
    zero imputation systematically lowers raw_insulation and shifts the
    log2 normalization, which is not equivalent to NaN-preserving.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Optional

import numpy as np
import pandas as pd

from .common import bin_to_interval, enforce_min_spacing, local_extrema_mask
from .results import BoundaryCallResult
from .thresholds import FixedThreshold, QuantileThreshold, ThresholdStrategy, fit_threshold, threshold_name
from ._validation import collect_qc_warnings, validate_for_mean_method


@dataclass
class InsulationParams:
    window_bp: Optional[int] = None
    window_bins: Optional[int] = None
    ignore_diags: int = 0
    min_frac_valid_pixels: float = 0.66
    min_dist_bad_bin: int = 0
    smooth_window_bins: int = 1
    local_min_window_bins: int = 1
    prominence_window_bins: Optional[int] = None
    threshold_strategy: ThresholdStrategy = QuantileThreshold(0.50)
    min_distance_bins: int = 2
    log_transform: bool = True
    min_score: float | None = None
    allow_negative: bool = False


def _resolve_params(params: InsulationParams | None, resolution: int | None = None) -> InsulationParams:
    p = params or InsulationParams()
    if p.window_bins is None:
        if p.window_bp is not None and resolution is not None:
            return replace(p, window_bins=max(1, int(p.window_bp) // int(resolution)))
        return replace(p, window_bins=10)
    return p


def _moving_nanmean(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return values.copy()
    out = np.full(values.shape, np.nan, dtype=float)
    half = window // 2
    for i in range(values.size):
        lo = max(0, i - half)
        hi = min(values.size, i + half + 1)
        if np.isfinite(values[lo:hi]).any():
            out[i] = np.nanmean(values[lo:hi])
    return out


def _insulation_pixel_values_slow(
    matrix: np.ndarray, i: int, window_bins: int, ignore_diags: int
) -> tuple[np.ndarray, int]:
    vals = []
    possible = 0
    n = matrix.shape[0]
    for r in range(i - window_bins, i):
        for c in range(i, i + window_bins):
            if 0 <= r < i and i <= c < n and c > r and c - r > ignore_diags and c - r <= window_bins:
                possible += 1
                v = matrix[r, c]
                if np.isfinite(v):
                    vals.append(float(v))
    return np.asarray(vals, dtype=float), int(possible)


def _valley_drop(score: np.ndarray, i: int, window: Optional[int] = None) -> float:
    if not np.isfinite(score[i]):
        return 0.0
    n = score.size
    if window is None:
        left = score[:i]
        right = score[i + 1 :]
    else:
        w = int(window)
        lo = max(0, i - w)
        hi = min(n, i + w + 1)
        left = score[lo:i]
        right = score[i + 1 : hi]
    left = left[np.isfinite(left)]
    right = right[np.isfinite(right)]
    if left.size == 0 or right.size == 0:
        return 0.0
    return float(max(0.0, min(np.max(left) - score[i], np.max(right) - score[i])))


def _expand_bad_bins(bad_bin_mask: np.ndarray, dist: int) -> np.ndarray:
    if dist <= 0:
        return bad_bin_mask.copy()
    n = bad_bin_mask.size
    excluded = np.zeros(n, dtype=bool)
    bad_positions = np.where(bad_bin_mask)[0]
    for pos in bad_positions:
        lo = max(0, int(pos) - int(dist))
        hi = min(n, int(pos) + int(dist) + 1)
        excluded[lo:hi] = True
    return excluded


def _finish_insulation(
    raw: np.ndarray,
    valid_fraction: np.ndarray,
    possible: np.ndarray,
    valid: np.ndarray,
    p: InsulationParams,
) -> dict:
    med = np.nanmedian(raw)
    if not np.isfinite(med) or med <= 0:
        insulation_score = np.full(raw.shape, np.nan)
    else:
        norm = raw / med
        with np.errstate(divide="ignore", invalid="ignore"):
            insulation_score = np.log2(norm) if p.log_transform else norm
        insulation_score[~np.isfinite(insulation_score)] = np.nan
    smoothed = _moving_nanmean(insulation_score, p.smooth_window_bins)
    candidate = local_extrema_mask(smoothed, "min", p.local_min_window_bins)
    strength = np.zeros(raw.shape, dtype=float)
    for i in np.where(candidate)[0]:
        strength[i] = _valley_drop(smoothed, int(i), p.prominence_window_bins)
    bad_bin_mask = (valid_fraction < p.min_frac_valid_pixels) | ~np.isfinite(raw)
    excluded_by_bad = _expand_bad_bins(bad_bin_mask, p.min_dist_bad_bin)
    valid_strengths = strength[np.isfinite(strength) & (strength > 0)]
    cutoff = fit_threshold(p.threshold_strategy, valid_strengths)
    if p.min_score is not None:
        cutoff = max(cutoff, float(p.min_score))
    threshold_pass = (
        candidate
        & np.isfinite(strength)
        & (strength >= cutoff)
        & (valid_fraction >= p.min_frac_valid_pixels)
        & ~excluded_by_bad
    )
    return {
        "raw_insulation": raw,
        "valid_pixel_fraction": valid_fraction,
        "possible_pixel_count": possible,
        "valid_pixel_count": valid,
        "insulation_score": smoothed,
        "candidate_mask": candidate,
        "boundary_strength": strength,
        "valley_drop": strength,
        "bad_bin_mask": bad_bin_mask,
        "excluded_by_bad_bin": excluded_by_bad,
        "threshold_values": valid_strengths,
        "threshold_cutoff": np.asarray(cutoff),
        "threshold_strategy_name": np.asarray(threshold_name(p.threshold_strategy)),
        "prominence_window_bins": np.asarray(
            -1 if p.prominence_window_bins is None else int(p.prominence_window_bins)
        ),
        "ignored_diags": np.asarray(p.ignore_diags),
        "threshold_pass": threshold_pass,
    }


def compute_insulation_slow(matrix: np.ndarray, params: InsulationParams | None = None) -> dict:
    p = _resolve_params(params)
    arr, preprocess_warnings = validate_for_mean_method(
        matrix,
        allow_negative=p.allow_negative,
        require_symmetric=True,
        window_bins=p.window_bins,
    )
    n = arr.shape[0]
    raw = np.full(n, np.nan)
    valid_fraction = np.zeros(n, dtype=float)
    possible_counts = np.zeros(n, dtype=int)
    valid_counts = np.zeros(n, dtype=int)
    for i in range(n):
        vals, possible = _insulation_pixel_values_slow(arr, i, p.window_bins, p.ignore_diags)
        possible_counts[i] = possible
        valid_counts[i] = vals.size
        valid_fraction[i] = vals.size / possible if possible else 0.0
        if possible and valid_fraction[i] >= p.min_frac_valid_pixels and vals.size:
            raw[i] = float(np.mean(vals))
    result = _finish_insulation(raw, valid_fraction, possible_counts, valid_counts, p)
    result["validation_warnings"] = list(preprocess_warnings)
    return result


def compute_insulation_production(matrix: np.ndarray, params: InsulationParams | None = None) -> dict:
    p = _resolve_params(params)
    arr, preprocess_warnings = validate_for_mean_method(
        matrix,
        allow_negative=p.allow_negative,
        require_symmetric=True,
        window_bins=p.window_bins,
    )
    n = arr.shape[0]
    raw = np.full(n, np.nan)
    valid_fraction = np.zeros(n, dtype=float)
    possible_counts = np.zeros(n, dtype=int)
    valid_counts = np.zeros(n, dtype=int)
    for i in range(n):
        rows = np.arange(max(0, i - p.window_bins), i)
        cols = np.arange(i, min(n, i + p.window_bins))
        if rows.size == 0 or cols.size == 0:
            continue
        rr, cc = np.meshgrid(rows, cols, indexing="ij")
        mask = (cc > rr) & ((cc - rr) > p.ignore_diags) & ((cc - rr) <= p.window_bins)
        possible = int(mask.sum())
        possible_counts[i] = possible
        if possible == 0:
            continue
        vals = arr[rr[mask], cc[mask]]
        vals = vals[np.isfinite(vals)]
        valid_counts[i] = vals.size
        valid_fraction[i] = vals.size / possible
        if valid_fraction[i] >= p.min_frac_valid_pixels and vals.size:
            raw[i] = float(np.mean(vals))
    result = _finish_insulation(raw, valid_fraction, possible_counts, valid_counts, p)
    result["validation_warnings"] = list(preprocess_warnings)
    return result


def call_boundaries(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    region_start: int = 0,
    params: Optional[InsulationParams] = None,
) -> BoundaryCallResult:
    p = _resolve_params(params, resolution)
    qc_warnings = collect_qc_warnings(matrix)
    out = compute_insulation_production(matrix, p)
    n = np.asarray(matrix).shape[0]
    region_end = int(region_start) + n * int(resolution)
    bins = np.arange(n)
    region = f"{chrom}:{region_start}-{region_end}"
    called = bins[out["threshold_pass"]]
    kept = enforce_min_spacing(called, out["boundary_strength"][called], p.min_distance_bins)
    threshold_pass = np.zeros(n, dtype=bool)
    threshold_pass[kept] = True
    score_rows = []
    for i in bins:
        iv = bin_to_interval(int(i), resolution, region_start, chrom)
        score_rows.append(
            {
                **iv,
                "resolution": resolution,
                "region_start": region_start,
                "region": region,
                "method": "insulation",
                "score": out["boundary_strength"][i],
                "boundary_candidate": bool(out["candidate_mask"][i]),
                "threshold_pass": bool(threshold_pass[i]),
                "warnings": "",
                "insulation_score": out["insulation_score"][i],
                "boundary_strength": out["boundary_strength"][i],
                "valley_drop": out["valley_drop"][i],
                "valid_pixel_fraction": out["valid_pixel_fraction"][i],
                "bad_bin": bool(out["bad_bin_mask"][i]),
                "excluded_by_bad_bin": bool(out["excluded_by_bad_bin"][i]),
            }
        )
    scores = pd.DataFrame(score_rows)
    b_rows = []
    for rank, i in enumerate(kept, start=1):
        iv = bin_to_interval(int(i), resolution, region_start, chrom)
        b_rows.append(
            {
                **iv,
                "resolution": resolution,
                "region_start": region_start,
                "region": region,
                "boundary_score": out["boundary_strength"][i],
                "method": "insulation",
                "rank": rank,
                "threshold": float(out["threshold_cutoff"]),
                "warnings": "",
            }
        )
    warnings = list(qc_warnings) + list(out.get("validation_warnings", []))
    if p.prominence_window_bins is None:
        warnings.append("valley_drop_uses_global_extrema_chunk_size_sensitive")
    if p.ignore_diags == 0:
        warnings.append("ignore_diags=0_differs_from_cooltools_default")
    return BoundaryCallResult(
        method_name="insulation",
        chrom=chrom,
        resolution=resolution,
        region_start=region_start,
        region_end=region_end,
        local_bin_indices=bins,
        boundaries=pd.DataFrame(b_rows),
        scores=scores,
        parameters={
            **asdict(p),
            "threshold_strategy": threshold_name(p.threshold_strategy),
            "threshold_kind": "per_chunk_quantile",
            "original_boundary_semantics": "boundary before bin i",
        },
        warnings=warnings,
        intermediates=out,
    )
