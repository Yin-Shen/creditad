"""Small, dependency-free contracts shared by the active HTTP runtime."""
from __future__ import annotations


DEFAULT_METHOD_SET = "tadvci"
SUPPORTED_METHOD_SETS = frozenset({"tadvci", "v2"})
RETIRED_METHOD_SETS = frozenset({"legacy"})


class MethodSetContractError(ValueError):
    """A public ``method_set`` request that must fail before task mutation."""

    def __init__(self, message: str, status_code: int = 400) -> None:
        super().__init__(message)
        self.status_code = status_code


def resolve_method_set(data: dict) -> str:
    """Resolve the public TAD pipeline without falling back to legacy state."""
    requested = data.get("method_set")
    if requested is None:
        return DEFAULT_METHOD_SET
    if not isinstance(requested, str) or not requested.strip():
        raise MethodSetContractError("method_set must be a non-empty string")
    method_set = requested.strip().lower()
    if method_set in RETIRED_METHOD_SETS:
        raise MethodSetContractError(
            "method_set 'legacy' is retired; use 'tadvci' for the active five-caller workflow",
            status_code=410,
        )
    if method_set not in SUPPORTED_METHOD_SETS:
        allowed = ", ".join(sorted(SUPPORTED_METHOD_SETS))
        raise MethodSetContractError(
            f"Unsupported method_set {requested!r}; allowed values: {allowed}"
        )
    return method_set
