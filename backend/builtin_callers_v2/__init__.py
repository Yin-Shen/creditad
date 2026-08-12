"""CrediTAD built-in clean-room TAD boundary callers (5-method active set).

Active method set:
    insulation        — Crane-style insulation score (diamond pixel, valley drop).
    topdom_like       — TopDom-like Wilcoxon ranksum / median-difference statistic.
    contact_contrast  — custom intra-domain vs cross-boundary contact contrast.
    spectral_profile  — spectral-eigenvector-based domain boundary profile.
    network_modularity — graph modularity-transition boundary profile.

DI (directionality_index) was dropped on 2026-05-13 after a sensitivity
audit produced all-zero boundary calls under cooler-balanced inputs.

Source is a clean-room implementation. No external TAD caller is wrapped.  The
separate explicit ``v2`` compatibility route deliberately selects the historical
four-method subset; the active default ``tadvci`` route uses all five.
"""

from __future__ import annotations

from .common import (
    bin_to_interval,
    enforce_min_spacing,
    local_extrema_mask,
    safe_nanmean,
    symmetrize_matrix,
    validate_matrix,
)
from .results import BoundaryCallResult
from .thresholds import (
    FixedThreshold,
    LiThreshold,
    NoThreshold,
    OtsuThreshold,
    QuantileThreshold,
    ThresholdStrategy,
    fit_threshold,
    threshold_name,
)

from .contact_contrast import ContactContrastParams
from .insulation import InsulationParams
from .network_modularity import NetworkModularityParams
from .spectral_profile import SpectralProfileParams
from .topdom_like import TopDomLikeParams

from .registry import CALLERS, get_caller, list_callers, run_all_callers

__all__ = [
    "BoundaryCallResult",
    "ContactContrastParams",
    "InsulationParams",
    "NetworkModularityParams",
    "SpectralProfileParams",
    "TopDomLikeParams",
    "FixedThreshold",
    "LiThreshold",
    "NoThreshold",
    "OtsuThreshold",
    "QuantileThreshold",
    "ThresholdStrategy",
    "fit_threshold",
    "threshold_name",
    "bin_to_interval",
    "enforce_min_spacing",
    "local_extrema_mask",
    "safe_nanmean",
    "symmetrize_matrix",
    "validate_matrix",
    "CALLERS",
    "get_caller",
    "list_callers",
    "run_all_callers",
]
