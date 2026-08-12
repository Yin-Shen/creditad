"""Missing-assay encoding for the BCP scorer.

A genuinely ABSENT ChIP assay (column missing OR all-NaN) must yield bcp_score = NaN
(not_applicable, matching the tier engine), NOT a fabricated fold=0.0 that the GBM
reads as 'strongly depleted'. Counts keep fillna(0). When both ChIP assays are present,
score_features returns finite probabilities.

The masking logic is verified deterministically against a tiny fake _model (a fitted
sklearn GBM on 4 toy features) so it does not depend on the shipped deploy artifact.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytestmark = pytest.mark.legacy_bcp

BACKEND = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BACKEND))

from tad_vci import bcp_infer    # noqa: E402

_FEATS = ["votes", "ctcf_fold", "rad21_fold", "loop_anchor_count"]


def _install_fake_model(monkeypatch):
    """Fit a tiny GBM on the 4 production features and install it as the loaded model,
    bypassing _load() / the on-disk artifact so the masking logic is tested in isolation."""
    from sklearn.ensemble import GradientBoostingClassifier

    rng = np.random.default_rng(0)
    n = 400
    votes = rng.integers(1, 6, n).astype(float)
    ctcf = rng.uniform(0.0, 8.0, n)
    rad21 = rng.uniform(0.0, 8.0, n)
    loops = rng.integers(0, 3, n).astype(float)
    X = np.column_stack([votes, ctcf, rad21, loops])
    # credible when supported by votes + CTCF/RAD21 signal above background
    y = ((votes >= 3) & ((ctcf > 1.0) | (rad21 > 1.0))).astype(int)
    gbm = GradientBoostingClassifier(random_state=0).fit(X, y)

    fake = {"gbm": gbm, "features": _FEATS, "apply_platt": False, "platt": None}
    monkeypatch.setattr(bcp_infer, "_model", fake, raising=False)
    monkeypatch.setattr(bcp_infer, "_loaded", True, raising=False)
    return fake


def _present_df(n=20):
    rng = np.random.default_rng(1)
    return pd.DataFrame({
        "votes": rng.integers(1, 6, n),
        "ctcf_fold": rng.uniform(0.5, 6.0, n),
        "rad21_fold": rng.uniform(0.5, 6.0, n),
        "loop_anchor_count": rng.integers(0, 3, n),
    })


def test_both_assays_present_returns_finite(monkeypatch):
    _install_fake_model(monkeypatch)
    df = _present_df()
    p = bcp_infer.score_features(df)
    assert p is not None
    assert len(p) == len(df)
    assert np.isfinite(p).all()
    assert ((p >= 0.0) & (p <= 1.0)).all()


def test_ctcf_column_absent_returns_all_nan(monkeypatch):
    """ctcf_fold column entirely missing -> assay not measured -> all-NaN (not_applicable)."""
    _install_fake_model(monkeypatch)
    df = _present_df().drop(columns=["ctcf_fold"])
    p = bcp_infer.score_features(df)
    assert p is not None
    assert len(p) == len(df)
    assert np.isnan(p).all()


def test_ctcf_all_nan_returns_all_nan(monkeypatch):
    """ctcf_fold present as a column but all-NaN -> still assay-absent -> all-NaN."""
    _install_fake_model(monkeypatch)
    df = _present_df()
    df["ctcf_fold"] = np.nan
    p = bcp_infer.score_features(df)
    assert p is not None
    assert np.isnan(p).all()


def test_absent_assay_not_treated_as_depleted(monkeypatch):
    """The masked (absent-assay) result must be NaN, distinctly different from feeding the
    GBM a fabricated fold=0.0 (which it would read as strongly depleted signal)."""
    _install_fake_model(monkeypatch)
    df = _present_df()
    masked = bcp_infer.score_features(df.drop(columns=["ctcf_fold"]))
    fabricated = bcp_infer.score_features(df.assign(ctcf_fold=0.0))
    assert np.isnan(masked).all()
    # the fabricated-zero path yields finite scores -> proves the two paths differ
    assert np.isfinite(fabricated).all()


def test_stray_nan_in_present_assay_filled_neutral(monkeypatch):
    """A single NaN bin in an otherwise-present assay is filled to 1.0 (at-background),
    so the overall result stays finite (assay is still 'present')."""
    _install_fake_model(monkeypatch)
    df = _present_df()
    df.loc[0, "ctcf_fold"] = np.nan
    p = bcp_infer.score_features(df)
    assert p is not None
    assert np.isfinite(p).all()


def test_counts_keep_fillna_zero(monkeypatch):
    """votes / loop_anchor_count are counts: a NaN there is filled with 0 (real value),
    and the result stays finite because the ChIP assays are present."""
    _install_fake_model(monkeypatch)
    df = _present_df()
    df.loc[0, "loop_anchor_count"] = np.nan
    df.loc[1, "votes"] = np.nan
    p = bcp_infer.score_features(df)
    assert p is not None
    assert np.isfinite(p).all()


def test_deploy_model_missing_assay_if_loadable():
    """If the shipped deploy model loads in this env, the same all-NaN behaviour must hold
    on the real GBM. If it is not loadable, skip gracefully (do not fabricate)."""
    bcp_infer._loaded = False
    bcp_infer._model = None
    if not bcp_infer.available():
        return  # deploy artifact not available here -> nothing to assert
    df = _present_df()
    p_full = bcp_infer.score_features(df)
    p_absent = bcp_infer.score_features(df.drop(columns=["ctcf_fold"]))
    assert p_full is not None and np.isfinite(p_full).all()
    assert p_absent is not None and np.isnan(p_absent).all()
