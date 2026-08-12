"""Threshold strategies for clean-room v2 boundary scores."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


class ThresholdStrategy(Protocol):
    name: str

    def fit(self, values) -> float:
        ...


def _finite_values(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


@dataclass(frozen=True)
class FixedThreshold:
    value: float
    name: str = "fixed"

    def fit(self, values) -> float:
        return float(self.value)


@dataclass(frozen=True)
class QuantileThreshold:
    quantile: float = 0.75
    name: str = "quantile"

    def fit(self, values) -> float:
        arr = _finite_values(values)
        if arr.size == 0:
            return float("inf")
        return float(np.quantile(arr, self.quantile))


@dataclass(frozen=True)
class NoThreshold:
    name: str = "none"

    def fit(self, values) -> float:
        return float("-inf")


@dataclass(frozen=True)
class LiThreshold:
    name: str = "li"

    def fit(self, values) -> float:
        arr = _finite_values(values)
        if arr.size == 0:
            return float("inf")
        try:
            from skimage.filters import threshold_li

            return float(threshold_li(arr))
        except Exception:
            return float(np.quantile(arr, 0.75))


@dataclass(frozen=True)
class OtsuThreshold:
    name: str = "otsu"

    def fit(self, values) -> float:
        arr = _finite_values(values)
        if arr.size == 0:
            return float("inf")
        try:
            from skimage.filters import threshold_otsu

            return float(threshold_otsu(arr))
        except Exception:
            return float(np.quantile(arr, 0.75))


def threshold_name(strategy: ThresholdStrategy | None) -> str:
    if strategy is None:
        return "none"
    return getattr(strategy, "name", strategy.__class__.__name__)


def fit_threshold(strategy: ThresholdStrategy | None, values) -> float:
    if strategy is None:
        return NoThreshold().fit(values)
    return float(strategy.fit(values))
