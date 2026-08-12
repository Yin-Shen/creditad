"""One-chromosome smoke test for the five built-in pure-Python callers.

Closes inventory items A10 and H2. Until now the built-in panel (the "Advanced:
built-in 5 pure-Python callers" import path) had no in-repo test at all — the only
coverage was a whole-genome manual run, which nobody re-runs.

The fixture is the checked-in synthetic cooler ``toy_8mb.mcool`` (320 bins x 25 kb =
8 Mb, eight 1 Mb block-diagonal TADs, built by ``_make_toy_mcool.py``). The whole
module runs in well under a second, so it is cheap enough to stay in the default suite.

This is a smoke test, not an accuracy test: it asserts that every caller returns a
non-empty boundary list, sorted, inside the requested region, on the bin grid, with the
declared method name — the failure modes that silently produced an empty BED in
production. The only accuracy-flavoured assertion is that the two callers that operate
on whole domains (network_modularity, corner_ranksum) recover the planted 1 Mb block
edges, which is what makes an all-zero or all-noise regression visible.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_BACKEND = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_BACKEND))

from builtin_callers_v2 import registry  # noqa: E402

FIXTURE = Path(__file__).parent / "toy_8mb.mcool"
RESOLUTION = 25_000
N_BINS = 320
CHROM = "chr_toy"
REGION_END = N_BINS * RESOLUTION            # 8 Mb
PLANTED_EDGES = [t * 1_000_000 for t in range(1, 8)]   # eight 1 Mb blocks -> 7 internal edges

cooler = pytest.importorskip("cooler", reason="cooler not installed on this host")
pytestmark = pytest.mark.skipif(
    not FIXTURE.is_file(), reason=f"toy fixture missing: {FIXTURE}"
)


@pytest.fixture(scope="module")
def toy_matrix() -> np.ndarray:
    clr = cooler.Cooler(f"{FIXTURE}::resolutions/{RESOLUTION}")
    assert clr.chromnames == [CHROM]
    mat = clr.matrix(balance=False).fetch(CHROM)
    assert mat.shape == (N_BINS, N_BINS)
    return mat


def test_registry_exposes_exactly_five_callers():
    assert registry.list_callers() == [
        "contact_contrast", "insulation", "network_modularity",
        "spectral_profile", "topdom_like",
    ]
    assert registry.list_public_callers() == [
        "contact_contrast", "corner_ranksum", "insulation",
        "laplacian_profile", "network_modularity",
    ]


@pytest.mark.parametrize("method_id", sorted(registry.CALLERS))
def test_each_builtin_caller_returns_a_usable_boundary_list(method_id: str, toy_matrix):
    result = registry.get_caller(method_id)(
        toy_matrix, CHROM, RESOLUTION, region_start=0
    )
    b = result.boundaries

    assert result.method_name == method_id
    assert result.chrom == CHROM
    assert result.resolution == RESOLUTION
    assert (result.region_start, result.region_end) == (0, REGION_END)

    assert not b.empty, f"{method_id} produced no boundaries on the toy matrix"
    assert {"chrom", "start", "end"} <= set(b.columns)

    starts = b["start"].astype(int).to_numpy()
    ends = b["end"].astype(int).to_numpy()
    assert (starts == np.sort(starts)).all(), f"{method_id} boundaries are not sorted"
    assert len(set(starts.tolist())) == len(starts), f"{method_id} emitted duplicates"
    assert starts.min() >= 0 and ends.max() <= REGION_END, f"{method_id} left the region"
    assert (starts % RESOLUTION == 0).all(), f"{method_id} is off the bin grid"
    assert (ends - starts == RESOLUTION).all()
    assert set(b["chrom"]) == {CHROM}
    # a caller that "finds" a boundary at every bin is as useless as one that finds none
    assert len(b) < N_BINS // 2


def test_run_all_callers_returns_all_five_in_one_pass(toy_matrix):
    results = registry.run_all_callers(toy_matrix, CHROM, RESOLUTION)
    assert set(results) == set(registry.CALLERS)
    assert all(not r.boundaries.empty for r in results.values())


@pytest.mark.parametrize("method_id", ["network_modularity", "topdom_like"])
def test_domain_scale_callers_recover_the_planted_block_edges(method_id: str, toy_matrix):
    """Guards against a caller degrading to noise while still returning a non-empty list."""
    result = registry.get_caller(method_id)(toy_matrix, CHROM, RESOLUTION, region_start=0)
    called = set(result.boundaries["start"].astype(int))
    hit = [e for e in PLANTED_EDGES if any(abs(e - c) <= 2 * RESOLUTION for c in called)]
    assert len(hit) == len(PLANTED_EDGES), (
        f"{method_id} recovered {len(hit)}/{len(PLANTED_EDGES)} planted 1 Mb edges"
    )


def test_one_chromosome_import_pipeline_end_to_end(tmp_path):
    """The Import path itself (mcool -> 5 callers -> votes -> evidence tier -> BEDs/TSV)
    on a single chromosome. No ChIP and no loops are supplied, so BD2/BD3/BD4 must stay
    not_assessable and the tier must come from BD1 alone."""
    import pandas as pd

    import mcool_bed

    out = mcool_bed.detect_TAD_boundaries_tadvci(
        str(FIXTURE), RESOLUTION, selected_chroms=[CHROM],
        hash_tag="smoke", output_dir=str(tmp_path),
        cell_line="Other", assembly="hg38",
    )

    for method_id in mcool_bed.TADVCI_METHOD_IDS:
        bed = Path(out[method_id])
        assert bed.is_file() and bed.stat().st_size > 0, f"{method_id} wrote an empty BED"

    ann = pd.read_csv(out["TADVCI_annotation"], sep="\t")
    assert len(ann) > 0
    assert set(ann["BD2_ctcf"]) == {"not_assessable"}
    assert set(ann["BD3_rad21"]) == {"not_assessable"}
    assert set(ann["BD4_loop_anchor"]) == {"not_assessable"}
    assert set(ann["D_tier"]) <= {"D1", "D2", "D3", "D4", "D5"}
    assert ann["votes"].max() >= 2, "no boundary agreed on by 2+ built-in callers"
    assert int(ann["votes"].max()) <= 5

    tier_bed = Path(out["TADVCI_tier"])
    assert tier_bed.is_file() and tier_bed.stat().st_size > 0
    tier_rows = [ln.split("\t") for ln in tier_bed.read_text().splitlines() if ln.strip()]
    assert len(tier_rows) == len(ann)
    assert all(int(r[1]) < int(r[2]) <= REGION_END for r in tier_rows)
