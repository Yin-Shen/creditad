"""Registry for the 5-method clean-room production caller set.

The 5 methods (insulation, corner_ranksum, contact_contrast, laplacian_profile,
network_modularity) expose a function-based API: ``call_boundaries(matrix, chrom,
resolution, region_start=0, params=None) -> BoundaryCallResult``.

NAMING (2026-06-06): two methods were renamed so their public IDs no longer
collide with the official upstream tools used in the TADShop ConsensusTAD
reproduction (the CRAN ``TopDom`` package and the Bioconductor ``SpectralTAD``
package). These five are clean-room in-house re-implementations / custom
statistics, NOT those official tools:
    topdom_like     -> corner_ranksum     (corner-block ranksum / Wilcoxon; TopDom-style statistic, our code)
    spectral_profile -> laplacian_profile  (normalised-Laplacian eigenvector profile; spectral, our code)
The implementations still live in topdom_like.py / spectral_profile.py (internal
filenames); only the public method ID changed. Old IDs remain accepted as
ALIASES so existing data / scripts keep resolving.

DI was removed 2026-05-13 after a sensitivity audit (all-zero calls under
cooler-balanced inputs).
"""

from __future__ import annotations

from typing import Callable

from . import contact_contrast, insulation, network_modularity, spectral_profile, topdom_like


# INTERNAL implementation keys are kept stable (topdom_like / spectral_profile) so
# the whole pipeline, saved data and legacy consumers (sql.py, the v2 demo) keep
# working unchanged. The PUBLIC / paper-facing names are the non-colliding ones
# (corner_ranksum / laplacian_profile) — use PUBLIC_NAMES / public_id() for any
# figure, manuscript, API label, or doc. get_caller() accepts EITHER name.
CALLERS: dict[str, Callable] = {
    "insulation": insulation.call_boundaries,
    "topdom_like": topdom_like.call_boundaries,
    "contact_contrast": contact_contrast.call_boundaries,
    "spectral_profile": spectral_profile.call_boundaries,
    "network_modularity": network_modularity.call_boundaries,
}

# internal id -> public / paper id (renamed 2026-06-06 to avoid colliding with the
# official CRAN TopDom and Bioconductor SpectralTAD used in the TADShop reproduction;
# our two are clean-room re-implementations, NOT those tools)
PUBLIC_NAMES: dict[str, str] = {
    "insulation": "insulation",
    "topdom_like": "corner_ranksum",
    "contact_contrast": "contact_contrast",
    "spectral_profile": "laplacian_profile",
    "network_modularity": "network_modularity",
}
# public / paper id -> internal id  (so get_caller accepts the new names too)
_PUBLIC_TO_INTERNAL: dict[str, str] = {v: k for k, v in PUBLIC_NAMES.items()}


def public_id(method_id: str) -> str:
    """internal id -> paper/display id (e.g. topdom_like -> corner_ranksum)."""
    return PUBLIC_NAMES.get(method_id, method_id)


def internal_id(method_id: str) -> str:
    """accept either name -> internal implementation key."""
    return _PUBLIC_TO_INTERNAL.get(method_id, method_id)


def list_callers() -> list[str]:
    return sorted(CALLERS)


def list_public_callers() -> list[str]:
    return sorted(PUBLIC_NAMES.values())


def get_caller(method_id: str) -> Callable:
    mid = internal_id(method_id)
    if mid not in CALLERS:
        raise KeyError(
            f"unknown built-in caller: {method_id!r}; available: {list_callers()} "
            f"(public: {list_public_callers()})"
        )
    return CALLERS[mid]


def run_all_callers(
    matrix,
    chrom: str,
    resolution: int,
    region_start: int = 0,
    params_by_method: dict | None = None,
):
    params_by_method = params_by_method or {}
    out = {}
    for method_id, caller in CALLERS.items():
        out[method_id] = caller(
            matrix,
            chrom,
            resolution,
            region_start=region_start,
            params=params_by_method.get(method_id),
        )
    return out
