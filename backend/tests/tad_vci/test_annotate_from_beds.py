"""Unit tests for the product annotate-from-BED path (no built-in callers)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from tad_vci.annotate_from_beds import annotate_from_beds


def test_single_caller_panel_no_chip(tmp_path: Path):
    primary = tmp_path / "primary.bed"
    primary.write_text("7\t100000\t125000\n7\t200000\t225000\n")
    out = annotate_from_beds(
        str(primary), [],
        output_dir=str(tmp_path), hash_tag="single", resolution=25000,
        cell_line="GM12878", assembly="hg19",
    )
    assert out["n_methods"] == 1
    df = pd.read_csv(out["TADVCI_annotation"], sep="\t")
    assert len(df) == 2
    assert (df["votes"] == 1).all()
    # Single-method panel cannot express cross-caller agreement: BD1 must floor to
    # "weak", never the vacuous "strong" that grade_votes(1, 1) used to return.
    assert set(df["BD1_caller_support"]) == {"weak"}
    assert set(df["BD2_ctcf"]) == {"not_assessable"}
    assert set(df["BD3_rad21"]) == {"not_assessable"}
    assert set(df["BD4_loop_anchor"]) == {"not_assessable"}
    assert "BD2" in df["not_assessable_criteria"].iloc[0]


def test_zero_support_boundary_not_forced_to_one(tmp_path: Path):
    """A primary boundary that NO panel method supports keeps votes=0 (honest),
    is not fabricated up to votes=1, and still grades BD1 weak (assessable)."""
    primary = tmp_path / "primary.bed"
    primary.write_text("7\t100000\t125000\n7\t900000\t925000\n")
    a = tmp_path / "A.bed"
    a.write_text("7\t100500\t125500\n")   # supports 100000 only
    b = tmp_path / "B.bed"
    b.write_text("7\t100300\t125300\n")   # supports 100000 only
    out = annotate_from_beds(
        str(primary),
        [(str(a), "A"), (str(b), "B")],
        output_dir=str(tmp_path), hash_tag="zero", resolution=25000,
    )
    assert out["n_methods"] == 2
    df = pd.read_csv(out["TADVCI_annotation"], sep="\t")
    row900 = df.loc[df["pos"] == 900000].iloc[0]
    assert int(row900["votes"]) == 0
    sc = row900["supporting_callers"]
    assert pd.isna(sc) or str(sc) == ""
    assert row900["BD1_caller_support"] == "weak"


def test_multi_method_votes_and_n_methods(tmp_path: Path):
    primary = tmp_path / "primary.bed"
    primary.write_text("7\t100000\t125000\n7\t200000\t225000\n")
    a = tmp_path / "A.bed"
    a.write_text("7\t100500\t125500\n")
    b = tmp_path / "B.bed"
    b.write_text("7\t100200\t125200\n7\t200100\t225100\n")
    out = annotate_from_beds(
        str(primary),
        [(str(a), "A"), (str(b), "B")],
        output_dir=str(tmp_path), hash_tag="multi", resolution=25000,
    )
    assert out["n_methods"] == 2
    assert out["bd1_panel"] == ["A", "B"]
    df = pd.read_csv(out["TADVCI_annotation"], sep="\t")
    row100 = df.loc[df["pos"] == 100000].iloc[0]
    row200 = df.loc[df["pos"] == 200000].iloc[0]
    assert int(row100["votes"]) == 2
    assert int(row200["votes"]) == 1
    # n=2: strong if votes>=2, else weak (ceil(0.8*2)=2, ceil(0.6*2)=2)
    assert row100["BD1_caller_support"] == "strong"
    assert row200["BD1_caller_support"] == "weak"


def test_primary_not_double_counted_as_method(tmp_path: Path):
    primary = tmp_path / "primary.bed"
    primary.write_text("7\t100000\t125000\n")
    other = tmp_path / "other.bed"
    other.write_text("7\t100200\t125200\n")
    out = annotate_from_beds(
        str(primary),
        [(str(primary), "PrimaryDup"), (str(other), "Other")],
        output_dir=str(tmp_path), hash_tag="dedup", resolution=25000,
    )
    assert out["n_methods"] == 1
    assert out["bd1_panel"] == ["Other"]
