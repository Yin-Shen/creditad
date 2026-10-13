"""A gzipped loop file must grade exactly like the plain file, and a loop file that yields
no anchors must stop the run.

Release v1.0.3 opened `.gz` loop files as text, failed on the first line, returned no
anchors and graded BD4 `none` at every boundary without an error; the repository ships
its loop calls gzipped, so passing them straight to `creditad annotate` produced a table
with a measured absence that was never measured."""
import gzip
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from mcool_bed import _parse_loop_anchors, require_loop_anchors  # noqa: E402
from tad_vci.annotate_from_beds import annotate_from_beds  # noqa: E402
from tad_vci.cli import main as cli_main  # noqa: E402

LOOPS = (
    "chr7\t100000\t110000\tchr7\t300000\t310000\n"
    "chr7\t500000\t510000\tchr7\t800000\t810000\n"
    "chr7\t520000\t530000\tchr7\t900000\t910000\n"
)


def _inputs(tmp_path):
    primary = tmp_path / "boundaries.bed"
    primary.write_text("chr7\t100000\t110000\nchr7\t500000\t510000\nchr7\t700000\t710000\n")
    panel = tmp_path / "callerA.bed"
    panel.write_text("chr7\t100000\t110000\nchr7\t700000\t710000\n")
    plain = tmp_path / "loops.bedpe"
    plain.write_text(LOOPS)
    gz = tmp_path / "loops.bedpe.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(LOOPS)
    return primary, panel, plain, gz


def _table(res):
    with open(res["TADVCI_annotation"]) as fh:
        return fh.read()


def test_gzip_and_plain_parse_identically(tmp_path):
    _, _, plain, gz = _inputs(tmp_path)
    a, b = _parse_loop_anchors(str(plain)), _parse_loop_anchors(str(gz))
    assert set(a) == set(b) == {"7"}
    assert a["7"].tolist() == b["7"].tolist() == [105000, 305000, 505000, 525000, 805000, 905000]


def test_gzip_and_plain_annotate_identically(tmp_path):
    primary, panel, plain, gz = _inputs(tmp_path)
    kw = dict(method_beds=[(str(panel), "callerA")], resolution=10000)
    r1 = annotate_from_beds(str(primary), loops_path=str(plain), output_dir=str(tmp_path / "p"),
                            hash_tag="plain", **kw)
    r2 = annotate_from_beds(str(primary), loops_path=str(gz), output_dir=str(tmp_path / "g"),
                            hash_tag="plain", **kw)
    assert _table(r1) == _table(r2)
    assert "strong" in _table(r1)          # the boundary at 500 kb has three loop ends


def test_unparseable_loop_file_raises(tmp_path):
    primary, panel, _, _ = _inputs(tmp_path)
    bad = tmp_path / "bad.txt"
    bad.write_text("not\ta\tloops\tfile\n")
    with pytest.raises(ValueError, match="yielded no loop anchors"):
        require_loop_anchors(str(bad))
    with pytest.raises(ValueError, match="yielded no loop anchors"):
        annotate_from_beds(str(primary), method_beds=[(str(panel), "callerA")], resolution=10000,
                           loops_path=str(bad), output_dir=str(tmp_path / "o"))


def test_cli_reports_unparseable_loop_file(tmp_path):
    primary, panel, _, _ = _inputs(tmp_path)
    bad = tmp_path / "bad.txt"
    bad.write_text("not\ta\tloops\tfile\n")
    with pytest.raises(SystemExit) as exc:
        cli_main(["annotate", "--boundaries", str(primary), "--panel", f"callerA={panel}",
                  "--loops", str(bad), "--resolution", "10000", "--out", str(tmp_path / "c")])
    assert "yielded no loop anchors" in str(exc.value)
