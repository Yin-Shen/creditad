"""Result container and file I/O for clean-room v2 boundary callers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _jsonify(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonify(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonify(v) for v in value]
    if isinstance(value, np.ndarray):
        return _jsonify(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if np.isnan(value):
            return None
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


@dataclass
class BoundaryCallResult:
    method_name: str
    chrom: str
    resolution: int
    region_start: int
    region_end: int
    local_bin_indices: np.ndarray
    boundaries: pd.DataFrame
    scores: pd.DataFrame
    parameters: dict
    warnings: list[str]
    intermediates: dict[str, Any]

    def save_boundaries_tsv(self, outdir: str | Path) -> Path:
        path = Path(outdir) / "boundaries.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        self.boundaries.to_csv(path, sep="\t", index=False)
        return path

    def save_scores_tsv(self, outdir: str | Path) -> Path:
        path = Path(outdir) / "scores.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        self.scores.to_csv(path, sep="\t", index=False)
        return path

    def save_parameters_json(self, outdir: str | Path) -> Path:
        path = Path(outdir) / "parameters.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "method_name": self.method_name,
            "chrom": self.chrom,
            "resolution": self.resolution,
            "region_start": self.region_start,
            "region_end": self.region_end,
            "parameters": self.parameters,
            "warnings": self.warnings,
        }
        path.write_text(json.dumps(_jsonify(payload), indent=2, sort_keys=True))
        return path

    def save_intermediates_npz(self, outdir: str | Path) -> Path:
        path = Path(outdir) / "intermediates.npz"
        path.parent.mkdir(parents=True, exist_ok=True)
        arrays = {}
        for key, value in self.intermediates.items():
            if isinstance(value, str):
                arrays[key] = np.asarray(value)
            elif value is None:
                arrays[key] = np.asarray(np.nan)
            else:
                arrays[key] = np.asarray(value)
        np.savez_compressed(path, **arrays)
        return path

    def save_all(self, outdir: str | Path) -> None:
        self.save_boundaries_tsv(outdir)
        self.save_scores_tsv(outdir)
        self.save_parameters_json(outdir)
        self.save_intermediates_npz(outdir)
