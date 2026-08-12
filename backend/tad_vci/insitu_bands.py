"""Derive BD2/BD3 cutpoints from the user's OWN ChIP track.

Why this exists
---------------
`graded_evidence.load_bands` can only return cutpoints for cell lines that were
pre-calibrated. For anything else the old behaviour was to hand back the GM12878
cutpoints while the fold numerator and the background denominator were both computed
from the user's track — a mixed scale nobody calibrated. Because every GM12878 cutpoint
is the lowest of the calibrated set, the error is strictly one-directional: grades
inflate and never deflate. Measured on HepG2 chr1 (2,351 boundaries, same input, only
the declared cell line changed): BD2 up 17.4%, BD3 up 31.2%, tier up 3.5%, zero
downgrades. With a single-caller panel BD3 moved 37.4%.

The honest alternative is to measure the same statistic on the user's own track: the
p75/p90/p97.5 of the fold over random symmetric +/-25 kb windows, which is exactly how the
packaged bands were produced (analyses/tad_vci_20260529/calibrate_bands.py). Costs about
2.5 s per track and reproduces the stored HepG2 calibration to under 1%.

BAND LADDER v2 (2026-07-29): weak = p75, moderate = p90, strong = p97.5; below p75 grades
none. The retired v1 floor pinned weak to fold = 1.0, which is roughly the MEDIAN of this
null — around half of RANDOM windows cleared it. This module must track the packaged
ladder exactly, or an in-situ track would be graded on a different scale from a packaged
one. It also returns the sampled null so the served card can report where a fold sits in
it, the same as for a calibrated cell.

The window MUST stay identical to the one the served numerator uses
(`mcool_bed._bw_window_signal`, symmetric +/-25 kb); a calibrate/serve window mismatch is
the exact bug this file's ancestor was written to remove.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

FLANK = 25_000          # matches mcool_bed._bw_window_signal
N_WINDOWS = 20_000      # matches the packaged calibration
SEED = 0                # single seed: the packaged bands use 5, the spread is <0.4%
AUTOSOMES = [str(i) for i in range(1, 23)]


def random_window_stats(path: str, seed: int = SEED, n: int = N_WINDOWS):
    """(folds, background_mean) over random +/-FLANK autosomal windows, or (None, None).

    The background is returned as well because the SERVED fold must be divided by the same
    denominator these percentiles were computed against. Using the bigWig header's covered
    mean instead leaves the numerator and the cutpoints on different scales, which is what
    made the first version of this fix still inflate every grade in one direction.
    """
    try:
        import pyBigWig
    except ImportError:
        return None, None
    try:
        bw = pyBigWig.open(path)
    except (RuntimeError, OSError):
        return None, None
    try:
        avail = bw.chroms()
        chroms = {}
        for c in AUTOSOMES:                       # accept both 'chr1' and '1' naming
            for name in (f"chr{c}", c):
                L = avail.get(name)
                if L and L > 2 * FLANK:
                    chroms[name] = L
                    break
        if not chroms:
            return None, None
        rng = np.random.default_rng(seed)
        names = list(chroms)
        vals: list[float] = []
        for _ in range(n):
            name = names[rng.integers(len(names))]
            p = int(rng.integers(FLANK, chroms[name] - FLANK))
            try:
                v = bw.stats(name, p - FLANK, p + FLANK, type="mean")[0]
            except (RuntimeError, KeyError):
                v = None
            if v is not None:
                vals.append(float(v))
    finally:
        bw.close()
    if len(vals) < n // 10:                       # too sparse to characterise the null
        return None, None
    arr = np.asarray(vals, dtype=float)
    bg = float(arr.mean())
    if not np.isfinite(bg) or bg <= 0:
        return None, None
    return arr / bg, bg


BAND_PERCENTILES = {"weak": 75.0, "moderate": 90.0, "strong": 97.5}
QUANTILE_GRID = list(range(0, 101))


def estimate_bands(path: str, seed: int = SEED, n: int = N_WINDOWS):
    """((bands), background_mean) from this track, or (None, None).

    bands = (('strong', p97.5), ('moderate', p90), ('weak', p75)) of the track's own
    random-window fold null — the v2 ladder, identical to the packaged calibration.
    """
    bands, bg, _ = estimate_bands_and_null(path, seed=seed, n=n)
    return bands, bg


def estimate_bands_and_null(path: str, seed: int = SEED, n: int = N_WINDOWS):
    """(bands, background_mean, null_quantiles) — as estimate_bands, plus the sampled null.

    null_quantiles is {"pct": [0..100], "fold": [...]}, the same shape stored per track in
    calibrated_bands.json, so a card built from an in-situ track can report the observed
    fold's percentile exactly as a packaged one does."""
    folds, bg = random_window_stats(path, seed=seed, n=n)
    if folds is None or folds.size == 0:
        return None, None, None
    cuts = {k: float(np.percentile(folds, q)) for k, q in BAND_PERCENTILES.items()}
    if not all(np.isfinite(v) for v in cuts.values()):
        return None, None, None
    # Two floors, both retained from v1: a band may never sit below the track's own
    # background (that would grade depleted signal as positive evidence), and the ladder
    # must stay ordered. All four packaged cells have p75 >= 1.19, so this floor is
    # inactive for them; it only guards a pathologically skewed user track.
    cuts["weak"] = max(cuts["weak"], 1.0)
    if not (cuts["weak"] <= cuts["moderate"] <= cuts["strong"]) or cuts["strong"] <= 1.0:
        return None, None, None
    quant = {"pct": [float(p) for p in QUANTILE_GRID],
             "fold": [float(v) for v in np.percentile(folds, QUANTILE_GRID)]}
    return (("strong", cuts["strong"]), ("moderate", cuts["moderate"]),
            ("weak", cuts["weak"])), bg, quant
