"""Fail-closed compatibility stub for the retired BCP/conformal component.

The full pre-audit implementation is preserved under
``backend/analyses/tad_vci_20260529/_quarantined_active_components/``.  Its
labels, split, feature alignment, calibration and checkpoint provenance failed
the July 2026 audit, so an installed active package must not load or score it.
"""
from __future__ import annotations


RETIRED_REASON = (
    "BCP/conformal inference is retired after the July 2026 label, split, "
    "feature-alignment and checkpoint-provenance audit"
)


def available() -> bool:
    """The retired model is never available through the active package."""
    return False


def _retired(*_args, **_kwargs):
    raise RuntimeError(RETIRED_REASON)


score_features = _retired
annotate_bcp = _retired
bcp_level = _retired
uncertain = _retired
concordance = _retired


def status_for_score(_value=None) -> dict[str, object]:
    """Machine-readable retired status for legacy callers that only need scope."""
    return {
        "available": False,
        "status": "retired_provenance_failure",
        "reason": RETIRED_REASON,
        "probability_or_coverage_claim_allowed": False,
    }


__all__ = [
    "RETIRED_REASON",
    "available",
    "score_features",
    "annotate_bcp",
    "bcp_level",
    "uncertain",
    "concordance",
    "status_for_score",
]
