"""A chromosome that a bigWig header lists but holds no data for was not measured.

bedGraphToBigWig writes every chromosome of the sizes file into the header, so a track
made for some chromosomes lists all of them. Release v1.0.4 read each window on such a
chromosome as zero signal and graded BD2/BD3 `none` at every boundary there: a measured
absence that was never measured. A chromosome missing from the header already gave
`not_assessable`; one listed without data must give the same.

A window without data on a chromosome that has data is still read as zero signal, because
tracks written with zero-valued intervals omitted carry no record there."""
import csv
import os
import sys

import pytest

pyBigWig = pytest.importorskip("pyBigWig")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from tad_vci.annotate_from_beds import annotate_from_beds  # noqa: E402

SIZE = 2_000_000


def _track(path, with_gap=False):
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([("chr7", SIZE), ("chr8", SIZE)])
    starts = list(range(0, SIZE, 1000))
    if with_gap:
        starts = [s for s in starts if not 1_400_000 <= s < 1_600_000]
    bw.addEntries(["chr7"] * len(starts), starts, ends=[s + 1000 for s in starts],
                  values=[3.0 if 480_000 <= s < 520_000 else 1.0 for s in starts])
    bw.close()
    return str(path)


def _rows(tmp_path, bw):
    primary = tmp_path / "boundaries.bed"
    primary.write_text("chr7\t500000\t510000\nchr7\t1500000\t1510000\nchr8\t500000\t510000\n")
    panel = tmp_path / "callerA.bed"
    panel.write_text("chr7\t500000\t510000\nchr8\t500000\t510000\n")
    res = annotate_from_beds(str(primary), method_beds=[(str(panel), "callerA")],
                             resolution=10000, ctcf_bw=bw, rad21_bw=bw, cell_line="GM12878",
                             output_dir=str(tmp_path / "out"))
    with open(res["TADVCI_annotation"]) as fh:
        return {(r["chrom"].replace("chr", ""), int(r["pos"])): r
                for r in csv.DictReader(fh, delimiter="\t")}


def test_listed_chromosome_without_data_is_not_assessable(tmp_path):
    rows = _rows(tmp_path, _track(tmp_path / "chr7_only.bw"))
    chr8 = rows[("8", 500000)]
    assert chr8["BD2_ctcf"] == "not_assessable"
    assert chr8["BD3_rad21"] == "not_assessable"
    assert "BD2" in chr8["not_assessable_criteria"] and "BD3" in chr8["not_assessable_criteria"]
    assert rows[("7", 500000)]["BD2_ctcf"] != "not_assessable"


def test_window_without_data_on_a_measured_chromosome_reads_as_zero(tmp_path):
    rows = _rows(tmp_path, _track(tmp_path / "gap.bw", with_gap=True))
    in_gap = rows[("7", 1500000)]
    assert in_gap["BD2_ctcf"] == "none"
    assert in_gap["BD3_rad21"] == "none"
    assert rows[("8", 500000)]["BD2_ctcf"] == "not_assessable"
