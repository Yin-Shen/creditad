"""The two demo data layers must stay distinguishable, and a matrix on the wrong
genome build must never be offered for a package.

Shipped layout (example_data/README.md):
  * evidence packages — genome-wide, always shipped
  * browser tracks    — chr7-only subsets, unless the user added full-genome files

The shipped chr7 matrices are hg19 while every current package is hg38, so without
the gate a fresh install would draw hg38 boundary ticks over an hg19 map.
"""
from pathlib import Path

import pytest

from tad_vci import api_routes


def test_shipped_chr7_matrix_is_per_cell_and_hg38(monkeypatch):
    """Without user downloads, each cell line must resolve to ITS OWN hg38 chr7 map.

    Both halves matter. A shared matrix silently draws one cell line's contacts
    under another's boundaries; a GRCh37 matrix silently misplaces every tick.
    """
    monkeypatch.setattr(api_routes, "_shared_raw_roots", lambda: [])
    seen: dict[str, Path] = {}
    for cell in ("GM12878", "IMR90", "HepG2"):
        mcool = api_routes._multicell_track_paths(cell).get("mcool")
        if mcool is None:
            pytest.skip(f"{cell} demo matrix not built on this machine")
        mcool = Path(mcool)
        assert cell.lower() in mcool.name.lower(), (
            f"{cell} resolved to {mcool.name}, which is not that cell line's matrix"
        )
        profile = api_routes._matrix_profile(mcool)
        # _matrix_profile normalises names to the bare form the rest of the app uses.
        assert profile["chroms"] == ["7"], f"{mcool.name} should be chr7-only"
        assert profile["assembly"] == "hg38", (
            f"{mcool.name} is {profile['assembly']}, but every shipped package is hg38"
        )
        seen[cell] = mcool

    assert len({p.resolve() for p in seen.values()}) == len(seen), (
        f"cell lines share a contact matrix: {seen}"
    )


def test_no_cell_agnostic_track_is_reachable(monkeypatch):
    """Every resolved demo track must name the cell line it belongs to.

    `CTCF_chr7.bigWig` / `RAD21_chr7.bigWig` were byte-identical copies of
    GM12878's tracks offered as a fallback for all three cell lines.
    """
    monkeypatch.setattr(api_routes, "_shared_raw_roots", lambda: [])
    for cell in ("GM12878", "IMR90", "HepG2"):
        tracks = api_routes._multicell_track_paths(cell)
        for role in ("ctcf", "rad21"):
            p = tracks.get(role)
            if p is None:
                continue
            assert cell.lower() in Path(p).name.lower(), (
                f"{cell} {role} resolved to {Path(p).name}, which is not that cell line's file"
            )


def test_retired_shared_chr7_matrix_is_not_reachable(monkeypatch):
    """The GRCh37 shared file must not come back as a fallback."""
    monkeypatch.setattr(api_routes, "_shared_raw_roots", lambda: [])
    for cell in ("GM12878", "IMR90", "HepG2"):
        mcool = api_routes._multicell_track_paths(cell).get("mcool")
        if mcool is None:
            continue
        assert Path(mcool).name != "chr7.mcool", (
            "the shared GRCh37 chr7.mcool is reachable again — it is GM12878-only "
            "and on the wrong build"
        )


def test_fallback_signal_tracks_are_chr7_subsets(monkeypatch):
    monkeypatch.setattr(api_routes, "_shared_raw_roots", lambda: [])
    for cell in ("GM12878", "IMR90", "HepG2"):
        tracks = api_routes._multicell_track_paths(cell)
        for key in ("ctcf", "rad21"):
            p = tracks.get(key)
            if p is None:
                continue
            assert "chr7" in Path(p).name.lower(), (
                f"{cell} {key}: shipped signal tracks must be the chr7 subsets, got {p}"
            )


def test_evidence_is_genome_wide_in_every_package():
    for cell in ("GM12878", "IMR90", "HepG2"):
        for res in (25000, 10000):
            ann = api_routes._multicell_combo_dir(cell, res) / "annotation.tsv"
            if not ann.is_file():
                pytest.skip(f"{cell} {res} package not built on this machine")
            chroms = api_routes._evidence_chroms(ann)
            assert len(chroms) >= 23, (
                f"{cell} {res}: evidence table covers {len(chroms)} chromosomes; the "
                "packages are documented as genome-wide"
            )
            assert chroms[0] == "1" and "X" in chroms
