"""BYO (bring-your-own-boundaries) regression suite — closes inventory item A8.

Everything here runs on synthetic inputs written into ``tmp_path``: a few BED lines and,
for the ChIP tests, a ~2 Mb two-chromosome bigWig built on the fly. No ENCODE download,
no demo package, no host-specific path. The point is that the BYO path can be regression
tested on any machine, because until now it was only ever exercised by ad-hoc scripts
in /tmp.

What is pinned here (current intended behaviour, 2026-07-26):

* no user ChIP / no loops  -> BD2, BD3 and BD4 are ``not_assessable``; the demo tracks
  are never silently substituted, and the tier comes from BD1 alone;
* a real panel            -> votes are counted against the panel that could vote on
  *that chromosome* (``n_methods_effective``), and the callers that were silent there
  are named in ``callers_absent_on_chrom``;
* a single-caller upload  -> BD1 is floored to ``weak`` (1/1 is not cross-caller support);
* vote tolerance          -> derived once from the resolution (``resolution * 1`` since
  2026-07-29, was ``* 2``) so the CLI and the API cannot disagree on the same input;
  BD4's proximity window is a separate, fixed 50 kb and does not scale with it;
* uncalibrated cell line + user ChIP -> the cutpoints are estimated in situ from the
  user's own track (``bands_source == "in_situ"``, ``bands_calibrated`` False); the
  GM12878 cutpoints are never borrowed silently.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tad_vci.annotate_from_beds import annotate_from_beds  # noqa: E402

RES = 25_000


def _write_bed(path: Path, chrom: str, starts, width: int = RES) -> Path:
    path.write_text("".join(f"{chrom}\t{s}\t{s + width}\n" for s in starts))
    return path


def _read(out: dict) -> pd.DataFrame:
    return pd.read_csv(out["TADVCI_annotation"], sep="\t")


# --------------------------------------------------------------------------------------
# 1. No user ChIP: the molecular axes must stay unassessed, not be filled from the demo
# --------------------------------------------------------------------------------------

def test_byo_without_chip_leaves_bd2_bd3_bd4_not_assessable(tmp_path: Path):
    primary = _write_bed(tmp_path / "primary.bed", "7", range(100_000, 500_001, 100_000))
    a = _write_bed(tmp_path / "A.bed", "7", [101_000, 201_000, 301_000])
    b = _write_bed(tmp_path / "B.bed", "7", [102_000, 202_000])

    out = annotate_from_beds(
        str(primary), [(str(a), "A"), (str(b), "B")],
        output_dir=str(tmp_path), hash_tag="nochip", resolution=RES,
        cell_line="Other", assembly="hg38",
    )
    df = _read(out)

    assert len(df) == 5
    assert set(df["BD2_ctcf"]) == {"not_assessable"}
    assert set(df["BD3_rad21"]) == {"not_assessable"}
    assert set(df["BD4_loop_anchor"]) == {"not_assessable"}
    # every row must SAY which axes went unassessed, in the served table itself
    assert set(df["not_assessable_criteria"]) == {"BD2,BD3,BD4"}
    # no signal columns were invented
    assert "ctcf" not in df.columns and "rad21" not in df.columns
    assert "ctcf_fold" not in df.columns and "rad21_fold" not in df.columns
    # tier therefore comes from BD1 alone: strong -> D4 (D5 needs a supporting axis)
    assert df.loc[df["votes"] == 2, "D_tier"].tolist() == ["D4", "D4"]
    assert set(df.loc[df["votes"] < 2, "D_tier"]) <= {"D1"}


# --------------------------------------------------------------------------------------
# 2. Panel votes
# --------------------------------------------------------------------------------------

def test_byo_panel_votes_name_their_supporters(tmp_path: Path):
    primary = _write_bed(tmp_path / "primary.bed", "7", [100_000, 300_000, 900_000])
    a = _write_bed(tmp_path / "A.bed", "7", [100_500, 300_500])
    b = _write_bed(tmp_path / "B.bed", "7", [101_000])
    c = _write_bed(tmp_path / "C.bed", "7", [99_000])

    out = annotate_from_beds(
        str(primary), [(str(a), "A"), (str(b), "B"), (str(c), "C")],
        output_dir=str(tmp_path), hash_tag="votes", resolution=RES,
    )
    assert out["n_methods"] == 3
    assert out["bd1_panel"] == ["A", "B", "C"]

    df = _read(out).set_index("pos")
    assert int(df.loc[100_000, "votes"]) == 3
    assert sorted(df.loc[100_000, "supporting_callers"].split(";")) == ["A", "B", "C"]
    assert int(df.loc[300_000, "votes"]) == 1
    assert df.loc[300_000, "supporting_callers"] == "A"
    assert int(df.loc[900_000, "votes"]) == 0
    # n=3: strong >= ceil(0.8*3)=3, moderate >= ceil(0.6*3)=2
    assert df.loc[100_000, "BD1_caller_support"] == "strong"
    assert df.loc[300_000, "BD1_caller_support"] == "weak"
    assert df.loc[900_000, "BD1_caller_support"] == "weak"


def test_byo_single_caller_upload_is_floored_to_weak(tmp_path: Path):
    """One BED and no panel: 1/1 agreement is vacuous, so BD1 must read weak and the
    served table must state the floor in `bd1_rule` rather than leave the reviewer to
    infer it from a `strong` that means nothing."""
    primary = _write_bed(tmp_path / "primary.bed", "7", [100_000, 200_000])

    out = annotate_from_beds(
        str(primary), [], output_dir=str(tmp_path), hash_tag="single", resolution=RES,
    )
    df = _read(out)

    assert out["n_methods"] == 1 and out["bd1_panel"] == []
    assert set(df["votes"]) == {1}
    assert set(df["n_methods_effective"]) == {1}
    assert set(df["BD1_caller_support"]) == {"weak"}
    assert "floored to weak" in df["bd1_rule"].iloc[0]
    # weak BD1 with nothing else assessable is the bottom tier
    assert set(df["D_tier"]) == {"D1"}


# --------------------------------------------------------------------------------------
# 3. Per-chromosome denominator (D13)
# --------------------------------------------------------------------------------------

def test_byo_denominator_is_the_panel_present_on_that_chromosome(tmp_path: Path):
    """A caller that emitted nothing on a chromosome cannot vote there, so it must not
    sit in that chromosome's BD1 denominator. Discriminating case: 2 of 3 callers agree
    on chr8 where the third is silent -> 2/2 = strong; the retired global denominator
    would have scored the same boundary 2/3 = moderate."""
    primary = tmp_path / "primary.bed"
    primary.write_text("7\t100000\t125000\n8\t400000\t425000\n")
    a = tmp_path / "A.bed"
    a.write_text("7\t100500\t125500\n8\t400500\t425500\n")
    b = tmp_path / "B.bed"
    b.write_text("7\t101000\t126000\n8\t401000\t426000\n")
    c = tmp_path / "C.bed"          # silent on chr8
    c.write_text("7\t100200\t125200\n")

    out = annotate_from_beds(
        str(primary), [(str(a), "A"), (str(b), "B"), (str(c), "C")],
        output_dir=str(tmp_path), hash_tag="perchrom", resolution=RES,
    )
    df = _read(out).set_index("chrom")

    assert out["n_methods"] == 3                       # declared panel size unchanged
    assert int(df.loc[7, "n_methods"]) == 3 and int(df.loc[8, "n_methods"]) == 3
    assert int(df.loc[7, "n_methods_effective"]) == 3
    assert int(df.loc[8, "n_methods_effective"]) == 2
    absent7 = df.loc[7, "callers_absent_on_chrom"]
    assert pd.isna(absent7) or str(absent7) == ""
    assert df.loc[8, "callers_absent_on_chrom"] == "C"
    assert int(df.loc[8, "votes"]) == 2
    assert df.loc[8, "BD1_caller_support"] == "strong"
    assert df.loc[7, "BD1_caller_support"] == "strong"


# --------------------------------------------------------------------------------------
# 4. Vote tolerance is derived from the resolution (CLI == API)
# --------------------------------------------------------------------------------------

@pytest.mark.parametrize("resolution", [10_000, 25_000])
def test_byo_default_vote_tolerance_is_one_bin(tmp_path: Path, resolution: int):
    """The default tolerance must be `resolution * 1` and nothing else. A hardcoded 50 kb
    used to make the API count votes the CLI did not, at the same resolution and on the
    same input. Tightened from 2 bins on 2026-07-29 (docs/MEASURED_LIMITATIONS.md §11);
    a call 2 bins away must no longer vote."""
    primary = _write_bed(tmp_path / "primary.bed", "7", [1_000_000], width=resolution)
    near = _write_bed(tmp_path / "near.bed", "7", [1_000_000 + resolution], width=resolution)
    two = _write_bed(tmp_path / "two.bed", "7", [1_000_000 + 2 * resolution], width=resolution)
    far = _write_bed(tmp_path / "far.bed", "7", [1_000_000 + 3 * resolution], width=resolution)

    out = annotate_from_beds(
        str(primary), [(str(near), "Near"), (str(two), "Two"), (str(far), "Far")],
        output_dir=str(tmp_path), hash_tag=f"tol{resolution}", resolution=resolution,
    )
    df = _read(out)
    assert int(df["vote_tol_bp"].iloc[0]) == resolution * 1
    # Near is 1 bin away (inside); Two at 2 bins and Far at 3 bins are outside
    assert int(df["votes"].iloc[0]) == 1
    assert df["supporting_callers"].iloc[0] == "Near"
    assert int(df["supporting_offsets_bp"].iloc[0]) == resolution


def test_byo_explicit_tolerance_still_overrides_the_derived_default(tmp_path: Path):
    primary = _write_bed(tmp_path / "primary.bed", "7", [1_000_000], width=10_000)
    far = _write_bed(tmp_path / "far.bed", "7", [1_030_000], width=10_000)

    default = annotate_from_beds(
        str(primary), [(str(far), "Far")],
        output_dir=str(tmp_path), hash_tag="dflt", resolution=10_000,
    )
    widened = annotate_from_beds(
        str(primary), [(str(far), "Far")],
        output_dir=str(tmp_path), hash_tag="wide", resolution=10_000, tol=50_000,
    )
    assert int(_read(default)["votes"].iloc[0]) == 0
    assert int(_read(default)["vote_tol_bp"].iloc[0]) == 10_000
    assert int(_read(widened)["votes"].iloc[0]) == 1
    assert int(_read(widened)["vote_tol_bp"].iloc[0]) == 50_000


# --------------------------------------------------------------------------------------
# 5. BYO with the user's own ChIP + loops (synthetic bigWig, no ENCODE file)
# --------------------------------------------------------------------------------------

try:                                   # Windows hosts ship winbbi (read-only) instead
    import pyBigWig
except ImportError:                    # pragma: no cover - host dependent
    pyBigWig = None

needs_pybigwig = pytest.mark.skipif(
    pyBigWig is None, reason="pyBigWig (bigWig writer) not installed on this host"
)

# 8 Mb, not 2 Mb: the in-situ bands are percentiles of random windows drawn from THIS
# track, and the v2 strong rung is the p97.5. On a 2 Mb x 2 chromosome fixture the three
# injected peaks occupy ~3.8% of all window positions, so they defined the top 2.5% of
# their own null and the strong cut landed among the peaks themselves. At 8 Mb they are
# ~0.9% of window positions, which is the regime a real genome is in.
CHROM_LEN = 8_000_000
PEAKS = (100_000, 200_000, 300_000)


def _synthetic_bigwig(path: Path, seed: int) -> str:
    """A 2 Mb two-chromosome fold-change-like track with peaks at PEAKS on chr7.

    Chromosome lengths are deliberately not hg19/hg38 chr1 lengths, so
    `validate_bigwig.assembly_of` returns None (undecidable) and the assembly guard
    stays out of the way of the grading assertions below.
    """
    import numpy as np

    chroms = [("chr1", CHROM_LEN), ("chr7", CHROM_LEN)]
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader(chroms)
    rng = np.random.default_rng(seed)
    step = 5_000
    for name, length in chroms:
        n = length // step
        vals = rng.gamma(1.0, 1.0, size=n) + 0.2
        if name == "chr7":
            for pk in PEAKS:
                vals[pk // step] += 40.0
        bw.addEntries(
            [name] * n,
            list(range(0, length, step)),
            ends=list(range(step, length + 1, step)),
            values=[float(v) for v in vals],
        )
    bw.close()
    return str(path)


def _byo_with_tracks(tmp_path: Path, **kwargs) -> pd.DataFrame:
    primary = _write_bed(tmp_path / "primary.bed", "7", range(100_000, 900_001, 100_000))
    a = _write_bed(tmp_path / "A.bed", "7", [p + 1_000 for p in range(100_000, 900_001, 100_000)])
    # B also calls at 800 kb — a position far from every injected peak. It is the only
    # way this fixture can produce BD1=strong with all three supporting axes graded none,
    # which is what exercises the no-support clause below.
    b = _write_bed(tmp_path / "B.bed", "7",
                   [p + 2_000 for p in list(range(100_000, 500_001, 100_000)) + [800_000]])
    loops = tmp_path / "loops.bedpe"
    loops.write_text(
        "chr7\t90000\t110000\tchr7\t190000\t210000\n"
        "chr7\t95000\t105000\tchr7\t290000\t310000\n"
    )
    ctcf = _synthetic_bigwig(tmp_path / "ctcf.bigWig", seed=0)
    rad21 = _synthetic_bigwig(tmp_path / "rad21.bigWig", seed=1)

    out = annotate_from_beds(
        str(primary), [(str(a), "A"), (str(b), "B")],
        output_dir=str(tmp_path), resolution=RES,
        ctcf_bw=ctcf, rad21_bw=rad21, loops_path=str(loops),
        cell_line="Other", **kwargs,
    )
    return _read(out)


@needs_pybigwig
def test_byo_user_chip_and_loops_grade_all_four_axes_with_in_situ_bands(tmp_path: Path):
    """Uncalibrated cell line + the user's own tracks: BD2/BD3/BD4 become assessable, and
    the BD2/BD3 cutpoints come from the user's OWN track (bands_source == in_situ), never
    from the packaged GM12878 calibration."""
    df = _byo_with_tracks(tmp_path, hash_tag="chip")

    assert set(df["bands_source"]) == {"in_situ"}
    assert not bool(df["bands_calibrated"].iloc[0])
    assert set(df["bands_cell"]) == {"in_situ:Other"}
    # the served fold and the cutpoints share one denominator, and it is the user's
    assert float(df["bg_ctcf"].iloc[0]) > 0
    # all four axes assessed -> nothing left in the not-assessable list
    assert set(df["not_assessable_criteria"].fillna("")) == {""}
    assert set(df["BD2_ctcf"]) <= {"strong", "moderate", "weak", "none"}
    assert set(df["BD3_rad21"]) <= {"strong", "moderate", "weak", "none"}
    assert set(df["BD4_loop_anchor"]) <= {"strong", "moderate", "none"}

    peaked = df[df["pos"].isin(PEAKS)]
    flat = df[~df["pos"].isin(PEAKS)]
    assert set(peaked["BD2_ctcf"]) == {"strong"}
    assert set(peaked["BD3_rad21"]) == {"strong"}
    # a 40x spike must not make every unrelated boundary look supported too
    assert "strong" not in set(flat["BD2_ctcf"])
    assert set(peaked["D_tier"]) == {"D5"}
    # BD1 strong while no measured supporting axis reaches its weak band -> v2 no-support clause
    depleted = df[(df["BD1_caller_support"] == "strong")
                  & (df["BD2_ctcf"] == "none") & (df["BD3_rad21"] == "none")
                  & (df["BD4_loop_anchor"] == "none")]
    assert len(depleted) >= 1
    assert set(depleted["D_tier"]) == {"D3"}
    assert set(df["tier_rule_version"]) == {"max_support_v2"}


@needs_pybigwig
def test_byo_not_assessable_mode_refuses_to_grade_an_uncalibrated_track(tmp_path: Path):
    """The conservative opt-out: with `uncalibrated_bands="not_assessable"` the user's
    ChIP is read but NOT graded, because no cutpoints for it exist. It must not fall
    back to GM12878's."""
    df = _byo_with_tracks(tmp_path, hash_tag="na", uncalibrated_bands="not_assessable")

    assert set(df["bands_source"]) == {"not_assessable"}
    assert set(df["BD2_ctcf"]) == {"not_assessable"}
    assert set(df["BD3_rad21"]) == {"not_assessable"}
    assert "ctcf_fold" not in df.columns and "rad21_fold" not in df.columns
    # BD4 does not depend on ChIP calibration, so it stays assessable
    assert set(df["BD4_loop_anchor"]) <= {"strong", "moderate", "none"}
    assert set(df["not_assessable_criteria"]) == {"BD2,BD3"}
