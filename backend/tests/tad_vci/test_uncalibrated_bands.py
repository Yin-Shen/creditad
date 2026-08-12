"""An uncalibrated cell line must not be graded against another cell line's cutpoints.

Borrowing the GM12878 cutpoints while the fold numerator and background come from the
user's own track mixes two scales. GM12878 has the lowest cutpoints of the calibrated
set, so the error was strictly one-directional. Measured on HepG2 chr1 (n=2,351, same
input, only the declared cell line changed):

    before: BD2 +17.4%, BD3 +31.2%, tier +3.5%, zero downgrades
    after : BD2   1.0% (12 up / 12 down), BD3 0.2%, tier 0.09%

These tests pin the policy, not those exact numbers (which depend on host ChIP files).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import os

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tad_vci.graded_evidence import load_bands  # noqa: E402


def test_unspecified_cell_line_is_never_reported_as_calibrated():
    """The silent path: cell_line None/'' used to return calibrated=True, which both
    suppressed the UI disclosure and let the caller swap in GM12878's stored background."""
    for cell in (None, "", "Other", "HCT116", "NoSuchCell"):
        _c, _r, calibrated, used = load_bands(cell)
        assert calibrated is False, f"{cell!r} must not claim a per-cell calibration"
        assert used == "GM12878"   # the cutpoints it would fall back to, if asked to


def test_known_cell_lines_stay_calibrated():
    for cell in ("GM12878", "K562", "IMR90", "HepG2", "hepg2"):
        _c, _r, calibrated, used = load_bands(cell)
        assert calibrated is True
        assert used.lower() == cell.lower()


def test_in_situ_bands_reproduce_the_packaged_calibration():
    """The in-situ estimator must agree with the packaged bands for a cell we DID
    calibrate — otherwise it is not a valid substitute for an uncalibrated one."""
    import json
    track = Path(os.environ.get("CREDITAD_DATA_ROOT", "__unset__")) / "hepg2/ctcf/ENCFF357NFO.bigWig"
    if not track.is_file():
        pytest.skip("host ChIP track not available")
    from tad_vci.insitu_bands import estimate_bands
    bands, bg = estimate_bands(str(track))
    assert bands is not None and bg and bg > 0
    ref = json.loads((Path(__file__).resolve().parents[2] / "tad_vci" /
                      "calibrated_bands.json").read_text())["HepG2"]["CTCF"]
    assert abs(bands[0][1] - ref["strong_fold"]) / ref["strong_fold"] < 0.05
    assert abs(bands[1][1] - ref["moderate_fold"]) / ref["moderate_fold"] < 0.05


def test_untouchable_track_falls_back_to_not_assessable_not_foreign_bands():
    """If the track cannot be characterised, BD2/BD3 must go not_assessable rather than
    be graded against cutpoints derived from a different experiment."""
    from tad_vci.insitu_bands import estimate_bands
    bands, bg = estimate_bands("/nonexistent/track.bigWig")
    assert bands is None and bg is None
