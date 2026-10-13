"""Loop files in the ENCODE / Juicer BEDPE layout must be read, not rejected.

ENCODE distributes HiCCUPS loop calls (for example ENCFF794OIH.bedpe.gz) with the column
names on a comment line that starts "#chr1", a second comment line naming the
juicer_tools version, and 33 tab-separated columns per loop. Release v1.0.4 took the
"#chr1" line for a named header, found no column called "chr1" and returned no anchors,
so `creditad annotate --loops` stopped with "yielded no loop anchors". By the BEDPE
definition the first six columns are chr1 x1 x2 chr2 y1 y2; comment lines are skipped.
"""
import gzip
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from mcool_bed import _parse_loop_anchors, require_loop_anchors  # noqa: E402
from tad_vci.annotate_from_beds import annotate_from_beds  # noqa: E402

# Header line copied from ENCFF794OIH.bedpe.gz.
JUICER_HEADER = (
    "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\tcolor\tobserved\t"
    "expectedBL\texpectedDonut\texpectedH\texpectedV\tfdrBL\tfdrDonut\tfdrH\tfdrV\t"
    "numCollapsed\tcentroid1\tcentroid2\tradius\thighRes_start_1\thighRes_end_1\t"
    "highRes_start_2\thighRes_end_2\tlocalX\tlocalY\tlocalObserved\tlocalPval\tlocalPeakID\n"
    "# juicer_tools version 2.13.07\n"
)
SIX = [
    ("chr7", 100000, 110000, "chr7", 300000, 310000),
    ("chr7", 500000, 510000, "chr7", 800000, 810000),
    ("chr7", 520000, 530000, "chr7", 900000, 910000),
]
EXTRA = (".", ".", ".", ".", "0,255,255", "16.0", "1.0795175", "1.0984638", "2.0739074",
         "1.509245", "0.0", "0.0", "0.0", "0.0", "1", "105000", "305000", "0", "104000",
         "106000", "304000", "306000", "105500", "305500", "4.0", "0.0047", "0")


def _six_col_text():
    return "".join("\t".join(map(str, r)) + "\n" for r in SIX)


def _juicer_text():
    assert len(JUICER_HEADER.split("\n")[0].split("\t")) == 6 + len(EXTRA) == 33
    return JUICER_HEADER + "".join("\t".join(map(str, r + EXTRA)) + "\n" for r in SIX)


def _write(path, text):
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return str(path)


@pytest.mark.parametrize("name", ["loops.bedpe", "loops.bedpe.gz"])
def test_juicer_layout_parses_like_the_six_column_file(tmp_path, name):
    six = _parse_loop_anchors(_write(tmp_path / "six.bedpe", _six_col_text()))
    juicer = _parse_loop_anchors(_write(tmp_path / name, _juicer_text()))
    assert set(juicer) == set(six) == {"7"}
    assert juicer["7"].tolist() == six["7"].tolist() == [105000, 305000, 505000, 525000,
                                                         805000, 905000]


def test_juicer_layout_annotates_identically_to_the_six_column_file(tmp_path):
    primary = tmp_path / "boundaries.bed"
    primary.write_text("chr7\t100000\t110000\nchr7\t500000\t510000\nchr7\t700000\t710000\n")
    panel = tmp_path / "callerA.bed"
    panel.write_text("chr7\t100000\t110000\nchr7\t700000\t710000\n")
    six = _write(tmp_path / "six.bedpe", _six_col_text())
    juicer = _write(tmp_path / "ENCFF_like.bedpe.gz", _juicer_text())
    kw = dict(method_beds=[(str(panel), "callerA")], resolution=10000, hash_tag="t")
    r1 = annotate_from_beds(str(primary), loops_path=six, output_dir=str(tmp_path / "a"), **kw)
    r2 = annotate_from_beds(str(primary), loops_path=juicer, output_dir=str(tmp_path / "b"), **kw)
    with open(r1["TADVCI_annotation"]) as a, open(r2["TADVCI_annotation"]) as b:
        t1, t2 = a.read(), b.read()
    assert t1 == t2
    assert "strong" in t1          # the boundary at 500 kb has three loop ends


def test_juicer_layout_with_no_loops_still_stops_the_run(tmp_path):
    """Header and comment lines alone are a file that parses to zero loops: the clear
    error stays, it is not graded as a measured absence."""
    empty = _write(tmp_path / "empty.bedpe.gz", JUICER_HEADER)
    assert _parse_loop_anchors(empty) == {}
    with pytest.raises(ValueError, match="yielded no loop anchors"):
        require_loop_anchors(empty)


def test_commented_file_with_fewer_than_six_columns_is_rejected(tmp_path):
    bad = _write(tmp_path / "bad.bedpe", "#chr1\tx1\tx2\n# note\nchr7\t1\t2\n")
    with pytest.raises(ValueError, match="yielded no loop anchors"):
        require_loop_anchors(bad)


def test_named_header_without_hash_is_still_read_by_name(tmp_path):
    """The HiCCUPS looplist form (named header, no '#') keeps its by-name reading."""
    text = "chr2\ty1\ty2\tchr1\tx1\tx2\n" + "".join(
        f"{c2}\t{y1}\t{y2}\t{c1}\t{x1}\t{x2}\n" for c1, x1, x2, c2, y1, y2 in SIX)
    a = _parse_loop_anchors(_write(tmp_path / "named.txt", text))
    assert a["7"].tolist() == [105000, 305000, 505000, 525000, 805000, 905000]
