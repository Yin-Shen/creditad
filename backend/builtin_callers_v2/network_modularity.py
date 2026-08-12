"""network_modularity — Network-features family TAD boundary caller.

MrTADFinder-style, SIMPLIFIED (not the exact MrTADFinder agglomerative algorithm,
analogous to how spectral_profile is "not exact SpectralTAD"). Treats the
balanced contact map as a weighted graph and finds contiguous communities that
are denser than expected GIVEN GENOMIC DISTANCE:

    B_ij = A_ij - gamma * E(|i-j|)      E(d) = mean balanced contact at separation d

The distance-corrected null E(d) (NOT a degree-product / Newman null) is the
defining choice: a degree-product null on a Hi-C map just re-discovers
near-diagonal contact density, collapsing the method back onto the insulation /
Linear-score family. Removing the distance-decay expectation makes the signal
about long-range within-domain enrichment instead of the local diamond drop.

Domains are contiguous genomic intervals, so the modularity objective is
maximized EXACTLY by an O(n^2) dynamic program over the chain (deterministic,
no Louvain seed). Boundaries are the domain starts.

Exposes the standard caller API:
    call_boundaries(matrix, chrom, resolution, region_start=0, params=None)
        -> BoundaryCallResult
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np
import pandas as pd

from ._validation import collect_qc_warnings
from .common import bin_to_interval
from .results import BoundaryCallResult


@dataclass
class NetworkModularityParams:
    gamma: float = 1.0            # resolution parameter of the distance-corrected modularity
    min_size_bins: int = 3        # >= 75 kb @ 25 kb
    max_size_bins: int = 240      # <= 6 Mb; bounds DP cost and rules out implausible TADs
    score_flank_bins: int = 10    # window for per-boundary cross-cut strength
    # Prominence filter (consistency with the other callers, which all threshold their
    # score). The chain-DP, run against a CHROMOSOME-GLOBAL distance expectation E, tiles
    # minimum-size domains across any sub-region whose contacts sit BELOW the global E
    # (B = A - gamma*E uniformly negative there), producing runs of spurious 75-kb
    # boundaries (e.g. GM12878 chr7 ~109-113/79/119 Mb). A real boundary is a LOCAL
    # PROMINENCE of the cross-cut strength, not a member of such a flat plateau, so we
    # keep a DP start only if its strength exceeds its local neighbourhood by an adaptive
    # floor. Set prominence_frac=0.0 to recover the old (unfiltered) behaviour.
    prominence_window_bins: int = 5
    prominence_frac: float = 0.3


def _resolve_params(params: Optional[NetworkModularityParams]) -> NetworkModularityParams:
    return params if params is not None else NetworkModularityParams()


def _expected_by_distance(matrix: np.ndarray) -> np.ndarray:
    """E(d) = nanmean over diagonal d, broadcast to full matrix."""
    n = matrix.shape[0]
    ed = np.zeros(n, dtype=np.float64)
    for d in range(n):
        diag = np.diagonal(matrix, d)
        m = np.nanmean(diag) if np.isfinite(diag).any() else 0.0
        ed[d] = m if np.isfinite(m) else 0.0
    idx = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    return ed[idx]


def _modularity_dp(B: np.ndarray, min_size: int, max_size: int) -> np.ndarray:
    """Exact chain-constrained modularity maximization. Returns domain-start bins."""
    n = B.shape[0]
    P = np.zeros((n + 1, n + 1), dtype=np.float64)
    P[1:, 1:] = np.cumsum(np.cumsum(B, axis=0), axis=1)

    def seg(a: int, b: int) -> float:                       # within [a,b) upper-tri sum
        rect = P[b, b] - P[a, b] - P[b, a] + P[a, a]
        return rect / 2.0                                   # symmetric, diagonal zeroed

    NEG = -1e18
    best = np.full(n + 1, NEG)
    best[0] = 0.0
    prev = np.zeros(n + 1, dtype=int)
    for b in range(1, n + 1):
        lo = max(0, b - max_size)
        hi = b - min_size + 1
        for a in range(lo, max(hi, lo + 1)):
            if a > b - 1:
                break
            if best[a] == NEG:
                continue
            v = best[a] + seg(a, b)
            if v > best[b]:
                best[b] = v
                prev[b] = a
        if best[b] == NEG:                                  # tail shorter than min_size
            a = max(0, b - 1)
            best[b] = best[a] + seg(a, b)
            prev[b] = a
    starts = []
    b = n
    while b > 0:
        a = prev[b]
        starts.append(a)
        b = a
    return np.array(sorted(s for s in starts if s > 0), dtype=int)


def _boundary_strength(matrix: np.ndarray, E: np.ndarray, gamma: float, flank: int) -> np.ndarray:
    """Per-bin cross-cut modularity deficit: mean(gamma*E - A) over the off-diagonal
    block straddling a cut before bin i. High => contacts across the cut fall below
    distance expectation => strong boundary."""
    n = matrix.shape[0]
    deficit = gamma * E - matrix                            # >0 where below expectation
    strength = np.zeros(n, dtype=np.float64)
    for i in range(n):
        a0, a1 = max(0, i - flank), i
        b0, b1 = i, min(n, i + flank)
        if a1 <= a0 or b1 <= b0:
            continue
        block = deficit[a0:a1, b0:b1]
        strength[i] = np.nanmean(block) if np.isfinite(block).any() else 0.0
    return strength


def _prominence_filter(starts: np.ndarray, strength: np.ndarray,
                       window: int, frac: float) -> np.ndarray:
    """Keep a DP domain-start only if its cross-cut strength is a LOCAL PROMINENCE:
    strength[i] minus the median strength of its +/-window neighbours must exceed an
    adaptive floor (frac * median strength over all DP starts). This removes the
    min-size 'plateau tiling' the chain DP produces in below-global-expectation regions,
    while keeping genuinely prominent boundaries. frac<=0 disables the filter."""
    if frac <= 0 or len(starts) == 0:
        return starts
    floor = frac * float(np.median(strength[starts]))
    n = len(strength)
    kept = []
    for s in starts:
        lo, hi = max(0, s - window), min(n, s + window + 1)
        nb = np.concatenate([strength[lo:s], strength[s + 1:hi]])
        base = float(np.median(nb)) if len(nb) else 0.0
        if strength[s] - base >= floor:
            kept.append(int(s))
    return np.array(sorted(kept), dtype=int)


def call_boundaries(
    matrix: np.ndarray,
    chrom: str,
    resolution: int,
    region_start: int = 0,
    params: Optional[NetworkModularityParams] = None,
) -> BoundaryCallResult:
    p = _resolve_params(params)
    qc_warnings = collect_qc_warnings(matrix)
    A = np.asarray(matrix, dtype=np.float64)
    n = A.shape[0]
    region_end = int(region_start) + n * int(resolution)
    region = f"{chrom}:{region_start}-{region_end}"

    E = _expected_by_distance(A)
    B = A - p.gamma * E
    B[~np.isfinite(B)] = 0.0                                # unknown pairs contribute nothing
    np.fill_diagonal(B, 0.0)
    dp_starts = _modularity_dp(B, p.min_size_bins, p.max_size_bins)
    strength = _boundary_strength(A, E, p.gamma, p.score_flank_bins)
    # prominence filter removes the below-global-expectation min-size tiling artifact
    starts = _prominence_filter(dp_starts, strength, p.prominence_window_bins, p.prominence_frac)

    bins = np.arange(n)
    candidates = set(int(s) for s in dp_starts)   # raw DP domain starts
    called = set(int(s) for s in starts)          # prominence-passed boundaries
    score_rows = []
    for i in bins:
        iv = bin_to_interval(int(i), resolution, region_start, chrom)
        score_rows.append({
            **iv,
            "resolution": resolution,
            "region_start": region_start,
            "region": region,
            "method": "network_modularity",
            "score": float(strength[i]),
            "boundary_candidate": bool(i in candidates),
            "threshold_pass": bool(i in called),
            "warnings": "",
            "modularity_boundary_strength": float(strength[i]),
        })
    scores = pd.DataFrame(score_rows)

    ranked = sorted(starts, key=lambda i: strength[i], reverse=True)
    rank_of = {int(i): r for r, i in enumerate(ranked, start=1)}
    b_rows = []
    for i in starts:
        iv = bin_to_interval(int(i), resolution, region_start, chrom)
        b_rows.append({
            **iv,
            "resolution": resolution,
            "region_start": region_start,
            "region": region,
            "boundary_score": float(strength[i]),
            "method": "network_modularity",
            "rank": rank_of[int(i)],
            "threshold": float(p.gamma),
            "warnings": "",
        })
    b_df = pd.DataFrame(b_rows, columns=[
        "chrom", "local_bin_index", "start", "end", "resolution", "region_start",
        "region", "boundary_score", "method", "rank", "threshold", "warnings",
    ])

    warnings = list(qc_warnings)
    warnings.append("MrTADFinder-style, simplified; not exact MrTADFinder")
    warnings.append("distance-corrected null E(d); deterministic chain modularity DP")

    return BoundaryCallResult(
        method_name="network_modularity",
        chrom=chrom,
        resolution=resolution,
        region_start=region_start,
        region_end=region_end,
        local_bin_indices=bins,
        boundaries=b_df,
        scores=scores,
        parameters={
            **asdict(p),
            "null_model": "distance_decay_expectation",
            "optimizer": "exact_chain_modularity_dp",
            "boundary_semantics": "domain start bin",
        },
        warnings=warnings,
        intermediates={"expected_by_distance_d0": float(E[0, 0]) if n else 0.0,
                       "n_domains": len(starts) + 1},
    )
