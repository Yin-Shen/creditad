"""One-shot script to create a small synthetic mcool for tests.

Generates a 320-bin x 320-bin block-diagonal contact matrix at 25 kb
resolution (i.e. 8 Mb of synthetic chromosome). Eight TAD-like blocks
of 40 bins each. Saved as a real cooler so cooltools / our Insulation
implementation can both read it.

Run once: `python _make_toy_mcool.py`. Output `toy_8mb.mcool` checked
in to the repo for deterministic CI.
"""
from __future__ import annotations

from pathlib import Path

import cooler
import numpy as np
import pandas as pd

OUT_PATH = Path(__file__).parent / "toy_8mb.mcool"
RESOLUTION = 25_000
N_BINS = 320  # 8 Mb
N_TADS = 8
TAD_SIZE = N_BINS // N_TADS  # 40 bins = 1 Mb each

CHROM = "chr_toy"
CHROM_LEN = N_BINS * RESOLUTION  # 8 Mb


def build_synthetic_matrix() -> np.ndarray:
    rng = np.random.default_rng(seed=42)
    m = np.zeros((N_BINS, N_BINS), dtype=np.float64)
    for t in range(N_TADS):
        lo, hi = t * TAD_SIZE, (t + 1) * TAD_SIZE
        intra = rng.gamma(shape=2.0, scale=10.0, size=(TAD_SIZE, TAD_SIZE))
        m[lo:hi, lo:hi] = (intra + intra.T) / 2.0
    bg = rng.gamma(shape=1.0, scale=1.0, size=(N_BINS, N_BINS))
    bg = (bg + bg.T) / 2.0
    distance_decay = np.exp(
        -np.abs(np.arange(N_BINS)[:, None] - np.arange(N_BINS)[None, :]) / 30.0
    )
    m = m + bg * distance_decay
    np.fill_diagonal(m, 0.0)
    m = (m + m.T) / 2.0
    return m


def main() -> None:
    matrix = build_synthetic_matrix()

    bins = pd.DataFrame(
        {
            "chrom": [CHROM] * N_BINS,
            "start": np.arange(0, CHROM_LEN, RESOLUTION, dtype=np.int64),
            "end": np.arange(RESOLUTION, CHROM_LEN + 1, RESOLUTION, dtype=np.int64),
        }
    )

    pixels = []
    for i in range(N_BINS):
        for j in range(i, N_BINS):
            v = matrix[i, j]
            if v > 0:
                pixels.append({"bin1_id": i, "bin2_id": j, "count": v})
    pixels_df = pd.DataFrame(pixels).astype(
        {"bin1_id": "int64", "bin2_id": "int64", "count": "float64"}
    )

    if OUT_PATH.exists():
        OUT_PATH.unlink()
    cooler.create_cooler(
        cool_uri=f"{OUT_PATH}::resolutions/{RESOLUTION}",
        bins=bins,
        pixels=pixels_df,
        ordered=True,
        symmetric_upper=True,
        dtypes={"count": "float64"},
    )
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
