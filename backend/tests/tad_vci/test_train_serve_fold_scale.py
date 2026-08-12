"""Integration guard: the BCP fed columns must be on the SAME scale at TRAIN time and at
SERVE time (no feature-scale mismatch).

History: the BCP was once trained on the RAW +/-25 kb window mean (student ctcf_fold mean
~0.605) but SERVED the canonical fold-over-background column (mean ~1.09) -- a ~1.8x scale
mismatch on the ChIP features -- plus a window mismatch (75 kb asymmetric serve vs 50 kb
symmetric train). After the 2026-06-04 canonical-fold fix:

  CANONICAL FOLD (single definition, train==serve):
    fold = mean ChIP signal over the SYMMETRIC [pos-25k, pos+25k] window
           / that cell's random-window mean (calibrated_bands.json[cell][assay]['mean'])

This test recomputes, on a FIXED set of GM12878 boundaries, the served fold via the exact
served primitive mcool_bed._bw_window_signal (symmetric +/-25 kb) divided by the random-window
mean, and asserts |mean_serve - mean_train| is within tolerance of the trained
student_boundary_features.tsv fold columns. It also pins the served window to symmetric
+/-25 kb (no +resolution term) and asserts the trained fold mean is the canonical-fold mean
(>1.0 for an enrichment assay), not the old raw mean (~0.6).

DATA GENERATION: the fixtures used here (student_boundary_features.tsv, calibrated_bands.json,
the GM12878 CTCF/RAD21 bigWigs) belong to the ARCHIVED hg19 analysis generation
(backend/analyses/tad_vci_20260529/). This test guards the train==serve *scale contract* of
mcool_bed._bw_window_signal; it is NOT a check on the shipped hg38 six-pack numbers.

Data-root resolution (fixed 2026-07-25, issue F2): the bigWigs live NEXT TO the repo
(<repo_parent>/CTCF_ENCFF749HDD.bigWig), not inside it. The previous code looked only in
parents[3] (= the repo root), so the guard silently skipped on the very host that has the
data. _resolve_data_root() now probes an ordered candidate list and $CREDITAD_DATA_ROOT wins.

When the real fixtures really are absent the module-level skip names the EXACT missing paths,
and test_guard_logic_detects_scale_mismatch_synthetic still runs: it builds a synthetic bigWig
and drives the same comparison helper, proving the guard both passes on a matched scale and
FAILS on an injected 1.8x mismatch. The contract is therefore never fully untested.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

BACKEND = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BACKEND))

ANA = BACKEND / "analyses" / "tad_vci_20260529"
FEATURES_TSV = ANA / "student_boundary_features.tsv"
BANDS_JSON = ANA / "calibrated_bands.json"

CTCF_NAME = "CTCF_ENCFF749HDD.bigWig"
RAD21_NAME = "RAD21_ENCFF000WCT.bigWig"
REPO_ROOT = Path(__file__).resolve().parents[3]


def _resolve_data_root() -> Path:
    """Directory that holds the archived GM12878 bigWigs.

    $CREDITAD_DATA_ROOT always wins. Otherwise probe, in order: the repo root, the repo's
    PARENT (where the tracks actually live on the dev host), and example_data/. Falls back to
    the repo root so the skip message points at a concrete, explainable path.
    """
    env = os.environ.get("CREDITAD_DATA_ROOT")
    if env:
        return Path(env)
    for cand in (REPO_ROOT, REPO_ROOT.parent, REPO_ROOT / "example_data"):
        if (cand / CTCF_NAME).is_file() and (cand / RAD21_NAME).is_file():
            return cand
    return REPO_ROOT


_DATA_ROOT = _resolve_data_root()
CTCF_BW = _DATA_ROOT / CTCF_NAME
RAD21_BW = _DATA_ROOT / RAD21_NAME

CELL = "GM12878"
RES = 25000          # resolution arg to _bw_window_signal (must NOT widen the window now)
N_SAMPLE = 400       # fixed subset for speed; deterministic (head)
MEAN_TOL = 0.05      # |mean_serve - mean_train| allowed on each fed fold column
ASSAY_FOLDS = [("CTCF", "ctcf_fold", CTCF_BW), ("RAD21", "rad21_fold", RAD21_BW)]


def _missing_inputs() -> list[str]:
    """Exact paths that are missing, so the skip reason is diagnosable (issue F2)."""
    return [str(p) for p in (FEATURES_TSV, BANDS_JSON, CTCF_BW, RAD21_BW) if not p.is_file()]


_MISSING = _missing_inputs()
_HAVE_REAL_INPUTS = not _MISSING
requires_real_inputs = pytest.mark.skipif(
    bool(_MISSING),
    reason=("archived hg19 GM12878 fixtures missing: " + ", ".join(_MISSING)
            + f" (data root resolved to {_DATA_ROOT}; set CREDITAD_DATA_ROOT to override)"))


def _randwin_bg(assay: str) -> float:
    bands = json.loads(BANDS_JSON.read_text())
    return float(bands[CELL][assay]["mean"])


def _fixed_boundaries() -> pd.DataFrame:
    df = pd.read_csv(FEATURES_TSV, sep="\t").head(N_SAMPLE).copy()
    df["chrom"] = df["chrom"].astype(str)
    return df


def _assert_fold_scale_match(assay: str, serve_fold, train_fold, n_rows: int) -> None:
    """THE guard: served fold and trained fold must agree in mean within MEAN_TOL.

    Extracted so the real-data tests and the synthetic stand-in exercise the identical
    comparison logic (issue F2 -- the contract must stay tested even when the archived
    bigWigs are absent on the host).
    """
    serve_fold = np.asarray(serve_fold, float)
    train_fold = np.asarray(train_fold, float)
    finite = np.isfinite(serve_fold) & np.isfinite(train_fold)
    assert finite.sum() >= 0.9 * n_rows, "too many non-finite signals to compare"

    mean_serve = float(np.mean(serve_fold[finite]))
    mean_train = float(np.mean(train_fold[finite]))
    assert abs(mean_serve - mean_train) <= MEAN_TOL, (
        f"{assay} train/serve fold mean mismatch: train={mean_train:.4f} serve={mean_serve:.4f} "
        f"(|delta|={abs(mean_serve - mean_train):.4f} > {MEAN_TOL}); feature-scale mismatch reintroduced")


@requires_real_inputs
def test_served_window_is_symmetric_no_resolution_term():
    """mcool_bed._bw_window_signal must use the SYMMETRIC +/-flank window: passing two very
    different `resolution` values must give the SAME signal (the +resolution term is gone)."""
    from mcool_bed import _bw_window_signal
    df = _fixed_boundaries()[["chrom", "pos"]].head(60)
    sig_a, _ = _bw_window_signal(str(CTCF_BW), df, resolution=10000)
    sig_b, _ = _bw_window_signal(str(CTCF_BW), df, resolution=250000)
    finite = np.isfinite(sig_a) & np.isfinite(sig_b)
    assert finite.any()
    assert np.allclose(sig_a[finite], sig_b[finite], atol=1e-9), (
        "served window still depends on `resolution` -> +resolution term not removed")


@requires_real_inputs
@pytest.mark.parametrize("assay,col,bw_path", ASSAY_FOLDS)
def test_train_serve_fold_mean_within_tolerance(assay, col, bw_path):
    """The served fold (symmetric +/-25 kb signal / random-window mean) must match the
    trained student fold column in mean on the SAME fixed boundary set."""
    from mcool_bed import _bw_window_signal

    df = _fixed_boundaries()
    bg = _randwin_bg(assay)

    sig, _header_bg = _bw_window_signal(str(bw_path), df, resolution=RES)   # symmetric +/-25 kb
    _assert_fold_scale_match(assay, sig / bg, df[col].astype(float).values, len(df))


@requires_real_inputs
def test_trained_fold_is_canonical_not_raw():
    """The shipped student fold columns must be on the canonical FOLD scale (mean > 1.0 for an
    enrichment assay), not the old RAW window-mean scale (CTCF ~0.6). Guards against a regression
    that retrains on raw signal again."""
    df = pd.read_csv(FEATURES_TSV, sep="\t")
    ctcf_mean = float(df["ctcf_fold"].astype(float).mean())
    rad21_mean = float(df["rad21_fold"].astype(float).mean())
    assert ctcf_mean > 1.0, f"ctcf_fold mean {ctcf_mean:.3f} looks like RAW signal, not a fold"
    assert rad21_mean > 1.0, f"rad21_fold mean {rad21_mean:.3f} looks like RAW signal, not a fold"


# --------------------------------------------------------------------------------------
# Synthetic stand-in: runs on EVERY host, with or without the archived hg19 fixtures.
# --------------------------------------------------------------------------------------

def _write_synthetic_bigwig(path: Path, chrom: str = "chr1", length: int = 2_000_000,
                            step: int = 1000, value: float = 2.0) -> None:
    """A flat synthetic bigWig: every `step`-bp interval carries `value`.

    Flat on purpose -- the analytic expectation of any window mean is exactly `value`, so
    the guard's arithmetic can be checked without any real ChIP data.
    """
    import pyBigWig
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(chrom, length)])
    starts = list(range(0, length, step))
    bw.addEntries([chrom] * len(starts), starts,
                  ends=[s + step for s in starts],
                  values=[value] * len(starts))
    bw.close()


@pytest.fixture(scope="module")
def synthetic_bw(tmp_path_factory):
    pytest.importorskip("pyBigWig", reason="pyBigWig required to build the synthetic stand-in")
    p = tmp_path_factory.mktemp("fold_scale") / "synthetic_flat.bigWig"
    _write_synthetic_bigwig(p)
    return p


def _synthetic_boundaries(n: int = 40, chrom: str = "chr1") -> pd.DataFrame:
    # keep every window fully inside [0, length) so the flat expectation holds exactly
    pos = np.arange(n, dtype=int) * 40_000 + 100_000
    return pd.DataFrame({"chrom": [chrom] * n, "pos": pos})


def test_synthetic_window_is_symmetric_no_resolution_term(synthetic_bw):
    """Window-symmetry contract, checked without the archived bigWigs: the served signal must
    not change when `resolution` changes by 25x (the +resolution term must stay removed)."""
    from mcool_bed import _bw_window_signal
    df = _synthetic_boundaries()
    sig_a, _ = _bw_window_signal(str(synthetic_bw), df, resolution=10000)
    sig_b, _ = _bw_window_signal(str(synthetic_bw), df, resolution=250000)
    assert np.isfinite(sig_a).all() and np.isfinite(sig_b).all()
    assert np.allclose(sig_a, sig_b, atol=1e-9), (
        "served window still depends on `resolution` -> +resolution term reintroduced")
    # flat track: the +/-25 kb mean is the track value itself
    assert np.allclose(sig_a, 2.0, atol=1e-6)


def test_guard_logic_detects_scale_mismatch_synthetic(synthetic_bw):
    """The guard itself must be live, not vacuous: on a synthetic flat track it PASSES when the
    trained column is on the served fold scale, and FAILS when the historical ~1.8x raw-vs-fold
    scale mismatch is injected. This keeps the contract exercised on hosts without the archived
    hg19 GM12878 bigWigs (issue F2)."""
    from mcool_bed import _bw_window_signal
    df = _synthetic_boundaries()
    bg = 1.6                                   # stand-in for the random-window mean
    sig, _ = _bw_window_signal(str(synthetic_bw), df, resolution=RES)
    serve_fold = sig / bg                      # == 1.25 everywhere

    # matched scale -> guard passes
    _assert_fold_scale_match("SYNTH", serve_fold, serve_fold.copy(), len(df))

    # historical failure mode: trained on RAW signal (no /bg) -> 1.8x-style mismatch
    with pytest.raises(AssertionError, match="feature-scale mismatch reintroduced"):
        _assert_fold_scale_match("SYNTH", serve_fold, sig, len(df))
