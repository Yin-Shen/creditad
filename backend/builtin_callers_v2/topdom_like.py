"""Clean-room v2 TopDom-like boundary caller.

Identity: TopDom-like transparent implementation. This is NOT exact TopDom R.

Key fix vs prior version:
    The production path previously used `possible = w * w` as the
    denominator for valid_fraction at every bin. At edges (i < w or
    i > n - w) the actual block shape is (r1 - r0) * (c1 - c0) which is
    strictly less than w * w. The current implementation computes the
    correct shape-based denominator so slow and production now agree on
    mean_cf and valid_fraction for every bin position.

scipy dependency:
    The ranksum p-value computation requires scipy.stats.ranksums. If
    scipy is unavailable AND stat_filter is True, call_boundaries raises
    RuntimeError at entry. If scipy is unavailable AND stat_filter is
    False, p-values are filled with NaN and a warning is added to the
    result.

Boundary callable range:
    The statistical test requires 2*w bins of flanking context on each
    side, so boundary_score is finite only for i in [2w, n - 2w]. Outside
    that band, mean_cf and local_ext are still reported but boundary_score
    is NaN and the bin cannot be called. This is by design and predates
    the denominator fix.

Non-finite handling:
    NaN-to-zero imputation upstream MUST NOT be applied: it dilutes
    mean_cf and shifts both the cross/flank Wilcoxon test and the median
    difference boundary_score. The validator and runner should preserve
    NaN at missing pixels and rely on valid_fraction gating instead.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .common import bin_to_interval, enforce_min_spacing, local_extrema_mask
from .results import BoundaryCallResult
from ._validation import collect_qc_warnings, validate_for_mean_method

try:
    from scipy.stats import ranksums

    _SCIPY_AVAILABLE = True
except Exception:
    ranksums = None
    _SCIPY_AVAILABLE = False


@dataclass
class TopDomLikeParams:
    window_bins: int = 5
    stat_filter: bool = True
    pvalue_threshold: float = 0.05
    min_valid_fraction: float = 0.5
    local_ext_window_bins: int = 1
    min_distance_bins: int = 2
    allow_negative: bool = False


def _ranksum_pvalue(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return np.nan
    if ranksums is None:
        return np.nan
    return float(ranksums(a, b).pvalue)


def _finish_topdom(
    mean_cf: np.ndarray,
    valid_fraction: np.ndarray,
    matrix: np.ndarray,
    p: TopDomLikeParams,
) -> dict:
    n = matrix.shape[0]
    local_ext = np.zeros(n, dtype=float)
    local_ext[~np.isfinite(mean_cf)] = -0.5
    local_ext[local_extrema_mask(mean_cf, "min", p.local_ext_window_bins)] = -1
    local_ext[local_extrema_mask(mean_cf, "max", p.local_ext_window_bins)] = 1
    pvals = np.full(n, np.nan)
    stat_mask = np.zeros(n, dtype=bool)
    boundary_score = np.full(n, np.nan)
    cross_mean = np.full(n, np.nan)
    flank_mean = np.full(n, np.nan)
    w = p.window_bins
    cross_n = np.zeros(n, dtype=int)
    flank_n = np.zeros(n, dtype=int)
    for i in range(n):
        if i - 2 * w < 0 or i + 2 * w > n:
            continue
        cross = matrix[i - w : i, i : i + w].ravel()
        left = matrix[i - 2 * w : i - w, i - w : i].ravel()
        right = matrix[i : i + w, i + w : i + 2 * w].ravel()
        flank = np.concatenate([left, right])
        cross = cross[np.isfinite(cross)]
        flank = flank[np.isfinite(flank)]
        cross_n[i] = cross.size
        flank_n[i] = flank.size
        if cross.size:
            cross_mean[i] = np.mean(cross)
        if flank.size:
            flank_mean[i] = np.mean(flank)
        if cross.size and flank.size:
            pvals[i] = _ranksum_pvalue(flank, cross)
            boundary_score[i] = np.median(flank) - np.median(cross)
            stat_mask[i] = (not p.stat_filter) or (
                np.isfinite(pvals[i]) and pvals[i] <= p.pvalue_threshold
            )
    return {
        "mean_cf": mean_cf,
        "local_ext": local_ext,
        "pvalues": pvals,
        "stat_filter_mask": stat_mask,
        "boundary_score": boundary_score,
        "valid_fraction": valid_fraction,
        "cross_mean": cross_mean,
        "flank_mean": flank_mean,
        "cross_values_count": cross_n,
        "flank_values_count": flank_n,
    }


def compute_topdom_like_slow(matrix: np.ndarray, params: TopDomLikeParams | None = None) -> dict:
    p = params or TopDomLikeParams()
    if p.stat_filter and not _SCIPY_AVAILABLE:
        raise RuntimeError(
            "topdom_like.compute_topdom_like_slow with stat_filter=True requires "
            "scipy.stats.ranksums; install scipy or set stat_filter=False."
        )
    arr, preprocess_warnings = validate_for_mean_method(
        matrix,
        allow_negative=p.allow_negative,
        require_symmetric=True,
        window_bins=p.window_bins,
    )
    n = arr.shape[0]
    mean_cf = np.full(n, np.nan)
    valid_fraction = np.zeros(n, dtype=float)
    for i in range(n):
        vals = []
        possible = 0
        for r in range(i - p.window_bins, i):
            for c in range(i, i + p.window_bins):
                if 0 <= r < n and 0 <= c < n:
                    possible += 1
                    if np.isfinite(arr[r, c]):
                        vals.append(arr[r, c])
        valid_fraction[i] = len(vals) / possible if possible else 0.0
        if possible and valid_fraction[i] >= p.min_valid_fraction and vals:
            mean_cf[i] = float(np.mean(vals))
    result = _finish_topdom(mean_cf, valid_fraction, arr, p)
    result["validation_warnings"] = list(preprocess_warnings)
    return result


def compute_topdom_like_production(matrix: np.ndarray, params: TopDomLikeParams | None = None) -> dict:
    p = params or TopDomLikeParams()
    if p.stat_filter and not _SCIPY_AVAILABLE:
        raise RuntimeError(
            "topdom_like.compute_topdom_like_production with stat_filter=True requires "
            "scipy.stats.ranksums; install scipy or set stat_filter=False."
        )
    arr, preprocess_warnings = validate_for_mean_method(
        matrix,
        allow_negative=p.allow_negative,
        require_symmetric=True,
        window_bins=p.window_bins,
    )
    n = arr.shape[0]
    mean_cf = np.full(n, np.nan)
    valid_fraction = np.zeros(n, dtype=float)
    w = p.window_bins
    for i in range(n):
        r0 = max(0, i - w)
        r1 = i
        c0 = i
        c1 = min(n, i + w)
        if r0 >= r1 or c0 >= c1:
            continue
        block = arr[r0:r1, c0:c1]
        possible = (r1 - r0) * (c1 - c0)
        valid = block[np.isfinite(block)]
        valid_fraction[i] = valid.size / possible if possible else 0.0
        if valid_fraction[i] >= p.min_valid_fraction and valid.size:
            mean_cf[i] = float(np.mean(valid))
    result = _finish_topdom(mean_cf, valid_fraction, arr, p)
    result["validation_warnings"] = list(preprocess_warnings)
    return result


def call_boundaries(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    region_start: int = 0,
    params: Optional[TopDomLikeParams] = None,
) -> BoundaryCallResult:
    p = params or TopDomLikeParams()
    qc_warnings = collect_qc_warnings(matrix)
    if p.stat_filter and not _SCIPY_AVAILABLE:
        raise RuntimeError(
            "topdom_like.call_boundaries with stat_filter=True requires scipy.stats.ranksums; "
            "install scipy or set TopDomLikeParams.stat_filter=False to skip the p-value filter."
        )
    out = compute_topdom_like_production(matrix, p)
    n = np.asarray(matrix).shape[0]
    bins = np.arange(n)
    region_end = region_start + n * resolution
    region = f"{chrom}:{region_start}-{region_end}"
    base = (
        (out["local_ext"] == -1)
        & np.isfinite(out["boundary_score"])
        & (out["boundary_score"] > 0)
    )
    if p.stat_filter:
        base &= out["stat_filter_mask"]
    called = bins[base]
    kept = enforce_min_spacing(called, out["boundary_score"][called], p.min_distance_bins)
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
                "method": "topdom_like",
                "score": out["boundary_score"][i],
                "boundary_candidate": bool(base[i]),
                "threshold_pass": bool(threshold_pass[i]),
                "warnings": "",
                "bin_signal": out["mean_cf"][i],
                "pvalue": out["pvalues"][i],
                "local_extreme_flag": out["local_ext"][i],
                "stat_filter_pass": bool(out["stat_filter_mask"][i]),
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
                "method": "topdom_like",
                "rank": rank,
                "threshold": p.pvalue_threshold if p.stat_filter else np.nan,
                "warnings": "",
            }
        )
    warnings = [
        "TopDom-like only; not exact TopDom R",
        *qc_warnings,
        *list(out.get("validation_warnings", [])),
    ]
    if not _SCIPY_AVAILABLE:
        warnings.append("scipy_unavailable_pvalues_set_to_nan")
    return BoundaryCallResult(
        method_name="topdom_like",
        chrom=chrom,
        resolution=resolution,
        region_start=region_start,
        region_end=region_end,
        local_bin_indices=bins,
        boundaries=pd.DataFrame(b_rows),
        scores=pd.DataFrame(score_rows),
        parameters={
            **asdict(p),
            "threshold_kind": "fixed_pvalue_threshold" if p.stat_filter else "none",
            "identity": "TopDom-like, not exact TopDom",
            "original_boundary_semantics": "boundary before bin i",
        },
        warnings=warnings,
        intermediates=out,
    )
