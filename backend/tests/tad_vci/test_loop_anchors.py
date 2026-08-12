"""Unit tests for the BD4 loop-anchor wiring added to mcool_bed:
_parse_loop_anchors (HiCCUPS-format parse + midpoint + chrom normalisation)
and the per-boundary loop_anchor_count -> grade_loop integration."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from mcool_bed import _parse_loop_anchors  # noqa: E402


def _write_loops(tmp_path):
    p = tmp_path / "loops.txt"
    # HiCCUPS header + 2 loops on chr7, 1 on chrX. midpoints: (x1+x2)//2, (y1+y2)//2
    p.write_text(
        "chr1\tx1\tx2\tchr2\ty1\ty2\to\n"
        "7\t100000\t110000\t7\t300000\t310000\t9\n"      # anchors 105000, 305000
        "chr7\t500000\t510000\tchr7\t800000\t810000\t9\n"  # anchors 505000, 805000
        "X\t10000\t20000\tX\t90000\t100000\t9\n"           # anchors 15000, 95000
    )
    return str(p)


def test_parse_loop_anchors_midpoints_and_norm(tmp_path):
    anchors = _parse_loop_anchors(_write_loops(tmp_path))
    assert set(anchors) == {"7", "X"}                       # 'chr' stripped, both ends kept
    assert sorted(anchors["7"].tolist()) == [105000, 305000, 505000, 805000]
    assert sorted(anchors["X"].tolist()) == [15000, 95000]


def test_parse_loop_anchors_bad_file_returns_empty(tmp_path):
    p = tmp_path / "bad.txt"
    p.write_text("not\ta\tloops\tfile\n1\t2\t3\t4\n")
    assert _parse_loop_anchors(str(p)) == {}


def test_loop_count_grades_bd4(tmp_path):
    """A boundary within 50 kb of >=2 anchors -> strong; 1 -> moderate; 0 -> none."""
    from tad_vci import grade_loop
    anchors = _parse_loop_anchors(_write_loops(tmp_path))["7"]
    tol = 50000

    def count(pos):
        return int(np.sum(np.abs(anchors - pos) <= tol))

    assert grade_loop(count(105000)) == "moderate"   # only 105000 within 50 kb
    assert grade_loop(count(305000)) == "moderate"   # only 305000
    assert grade_loop(count(700000)) == "none"       # nearest anchor 505000/805000 > 50 kb
    # two anchors within 50 kb of a midpoint between 505000 and 545000? craft strong:
    near_two = _parse_loop_anchors(
        _make_two_close(tmp_path))["7"]
    assert grade_loop(int(np.sum(np.abs(near_two - 600000) <= tol))) == "strong"


def _make_two_close(tmp_path):
    p = tmp_path / "loops_close.txt"
    p.write_text(
        "chr1\tx1\tx2\tchr2\ty1\ty2\to\n"
        "7\t580000\t600000\t7\t610000\t630000\t9\n"   # anchors 590000, 620000 (both within 50 kb of 600000)
    )
    return str(p)
