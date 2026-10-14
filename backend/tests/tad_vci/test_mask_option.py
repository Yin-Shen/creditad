"""`--mask`: BD2/BD3 not_assessable inside regions where ChIP signal is not a measurement.

A boundary inside an assembly gap or an exclusion-list region used to be graded on CTCF and
RAD21 like any other position, so a low fold there read as a measured absence. With a mask,
such a boundary records BD2/BD3 as not_assessable -- the same record as a run where no ChIP
track was supplied -- and its tier is re-derived by the released rule. Without a mask the
output is unchanged. The bundled hg38 mask reproduces the v1.1 tables from the v1.0 tables.
"""
from __future__ import annotations

import asyncio
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

BACKEND = Path(__file__).resolve().parents[2]
REPO = BACKEND.parent
sys.path.insert(0, str(BACKEND))

from tad_vci.annotate_from_beds import annotate_from_beds  # noqa: E402
from tad_vci.graded_evidence import combine_tier  # noqa: E402
from tad_vci import masking  # noqa: E402

SIZE = 3_000_000
RES = 10_000
# 500000 and 900000 inside the mask; 1500000 outside; 1200000 is the mask END (half-open,
# outside); 2000000 is a mask START (inside).
POSITIONS = [500_000, 900_000, 1_200_000, 1_500_000, 2_000_000]
MASK_TEXT = ("chr7\t400000\t600000\tassembly_gap\n"
             "chr7\t850000\t1200000\tencode_blacklist\n"
             "chr7\t880000\t950000\tassembly_gap\n"
             "chr7\t2000000\t2100000\tencode_blacklist\n")
INSIDE = {500_000, 900_000, 2_000_000}


def _write_track(path: Path) -> str:
    pyBigWig = pytest.importorskip("pyBigWig")
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([("chr7", SIZE)])
    starts = list(range(0, SIZE, 1000))
    vals = [6.0 if any(abs(s - p) < 20_000 for p in POSITIONS) else 0.4 for s in starts]
    bw.addEntries(["chr7"] * len(starts), starts, ends=[s + 1000 for s in starts], values=vals)
    bw.close()
    return str(path)


def _inputs(tmp_path: Path) -> dict:
    primary = tmp_path / "boundaries.bed"
    primary.write_text("".join(f"chr7\t{p}\t{p + RES}\n" for p in POSITIONS))
    panel_a = tmp_path / "A.bed"
    panel_a.write_text("".join(f"chr7\t{p}\t{p + RES}\n" for p in POSITIONS))
    panel_b = tmp_path / "B.bed"
    panel_b.write_text("".join(f"chr7\t{p}\t{p + RES}\n" for p in POSITIONS[:3]))
    loops = tmp_path / "loops.bedpe"
    loops.write_text(f"chr7\t{POSITIONS[0]}\t{POSITIONS[0] + RES}\t"
                     f"chr7\t{POSITIONS[3]}\t{POSITIONS[3] + RES}\n")
    mask = tmp_path / "mask.bed"
    mask.write_text(MASK_TEXT)
    return {"primary": str(primary), "panel": [(str(panel_a), "A"), (str(panel_b), "B")],
            "loops": str(loops), "mask": str(mask)}


def _run(tmp_path: Path, name: str, inp: dict, chip: str | None, **kw) -> dict:
    return annotate_from_beds(inp["primary"], method_beds=inp["panel"], resolution=RES,
                              ctcf_bw=chip, rad21_bw=chip, loops_path=inp["loops"],
                              cell_line="GM12878", assembly="hg38", hash_tag=name,
                              output_dir=str(tmp_path / name), **kw)


def _table(result: dict) -> pd.DataFrame:
    return pd.read_csv(result["TADVCI_annotation"], sep="\t", dtype=str,
                       keep_default_na=False).set_index("pos")


@pytest.fixture
def runs(tmp_path):
    inp = _inputs(tmp_path)
    chip = _write_track(tmp_path / "chip.bw")
    return {
        "inp": inp,
        "plain": _run(tmp_path, "plain", inp, chip),
        "masked": _run(tmp_path, "masked", inp, chip, mask_bed=inp["mask"]),
        "no_chip": _run(tmp_path, "nochip", inp, None),
    }


def test_masked_positions_record_bd2_bd3_not_assessable(runs):
    plain, masked = _table(runs["plain"]), _table(runs["masked"])
    for pos in INSIDE:
        row = masked.loc[str(pos)]
        assert plain.loc[str(pos), "BD2_ctcf"] != "not_assessable"
        assert row["BD2_ctcf"] == "not_assessable"
        assert row["BD3_rad21"] == "not_assessable"
        assert {"BD2", "BD3"} <= set(row["not_assessable_criteria"].split(","))
        tier, meta = combine_tier({"BD1": row["BD1_caller_support"], "BD2": "not_assessable",
                                   "BD3": "not_assessable", "BD4": row["BD4_loop_anchor"]})
        assert row["D_tier"] == tier
        assert row["not_assessable_criteria"] == ",".join(meta["not_assessable"])
    assert runs["masked"]["n_masked"] == len(INSIDE)


def test_bd1_bd4_and_raw_signal_are_kept_on_masked_rows(runs):
    plain, masked = _table(runs["plain"]), _table(runs["masked"])
    changed = {"BD2_ctcf", "BD3_rad21", "D_tier", "not_assessable_criteria"}
    keep = [c for c in plain.columns if c not in changed]
    assert list(masked.columns) == list(plain.columns)
    pd.testing.assert_frame_equal(masked[keep], plain[keep])


def test_rows_outside_the_mask_are_unchanged(runs):
    plain, masked = _table(runs["plain"]), _table(runs["masked"])
    outside = [str(p) for p in POSITIONS if p not in INSIDE]
    pd.testing.assert_frame_equal(masked.loc[outside], plain.loc[outside])


def test_masked_record_equals_a_run_without_chip_tracks(runs):
    """Inside the mask the record says what a run without ChIP input says."""
    masked, no_chip = _table(runs["masked"]), _table(runs["no_chip"])
    for pos in INSIDE:
        for col in ("BD1_caller_support", "BD2_ctcf", "BD3_rad21", "BD4_loop_anchor",
                    "D_tier", "not_assessable_criteria"):
            assert masked.loc[str(pos), col] == no_chip.loc[str(pos), col], (pos, col)


def test_mask_records_file(runs):
    rec = pd.read_csv(runs["masked"]["mask_records"], sep="\t", dtype=str,
                      keep_default_na=False)
    assert list(rec.columns) == masking.RECORD_COLUMNS
    assert sorted(int(p) for p in rec["pos"]) == sorted(INSIDE)
    src = dict(zip(rec["pos"].astype(int), rec["mask_source"]))
    assert src[500_000] == "assembly_gap"
    assert src[900_000] == "assembly_gap+encode_blacklist"
    assert src[2_000_000] == "encode_blacklist"
    plain = _table(runs["plain"])
    for r in rec.itertuples():
        assert r.BD2_ctcf_unmasked == plain.loc[r.pos, "BD2_ctcf"]
        assert r.D_tier_unmasked == plain.loc[r.pos, "D_tier"]


def test_tier_bed_follows_the_masked_tier(runs):
    bed = pd.read_csv(runs["masked"]["TADVCI_tier"], sep="\t", header=None, dtype=str)
    masked = _table(runs["masked"])
    assert dict(zip(bed[1], bed[3])) == dict(zip(masked.index, masked["D_tier"]))


def test_no_mask_leaves_output_and_result_unchanged(tmp_path):
    inp = _inputs(tmp_path)
    chip = _write_track(tmp_path / "chip.bw")
    plain = _run(tmp_path, "plain", inp, chip)
    far = tmp_path / "far.bed"
    far.write_text("chr7\t2900000\t2950000\nchr9\t0\t1000000\n")
    with_far = _run(tmp_path, "far", inp, chip, mask_bed=str(far))
    assert set(plain) == {"TADVCI_annotation", "TADVCI_tier", "n_boundaries", "n_methods",
                          "bd1_panel", "tier_dist", "primary_bed"}
    assert sorted(p.name for p in (tmp_path / "plain").iterdir()) == [
        "tadvci_annotation_beds_plain.tsv", "tadvci_tier_beds_plain.bed"]
    assert with_far["n_masked"] == 0
    assert (Path(plain["TADVCI_annotation"]).read_bytes()
            == Path(with_far["TADVCI_annotation"]).read_bytes())
    assert Path(plain["TADVCI_tier"]).read_bytes() == Path(with_far["TADVCI_tier"]).read_bytes()


def test_flag_positions_is_half_open_and_ignores_the_chr_prefix(tmp_path):
    p = tmp_path / "m.bed.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("# comment\ntrack name=x\n7\t100\t200\nchr8\t0\t50\tlabel\n")
    m = masking.load_mask(str(p))
    chrom = ["7", "7", "chr7", "7", "8", "chr8", "9"]
    pos = [99, 100, 199, 200, 49, 50, 150]
    assert masking.flag_positions(m, chrom, pos).tolist() == [
        False, True, True, False, True, False, False]
    assert masking.mask_sources(m, "chr7", 150) == "m.bed.gz"
    assert masking.mask_sources(m, "8", 0) == "label"


def test_malformed_masks_are_refused(tmp_path):
    bad = tmp_path / "bad.bed"
    bad.write_text("chr7\t200\t100\n")
    with pytest.raises(ValueError, match="invalid interval"):
        masking.load_mask(str(bad))
    empty = tmp_path / "empty.bed"
    empty.write_text("# nothing\n")
    with pytest.raises(ValueError, match="no parseable intervals"):
        masking.load_mask(str(empty))
    with pytest.raises(FileNotFoundError):
        masking.resolve_mask(str(tmp_path / "missing.bed"))


EXPECTED_HG38_LABELS = {"assembly_gap": 603, "cytoband_acen": 48, "cytoband_gvar": 17,
                        "cytoband_stalk": 5, "encode_blacklist": 910}


def test_bundled_hg38_mask():
    path = Path(masking.resolve_mask("hg38"))
    assert path.is_file() and path.parent == masking.MASK_DIR
    assert masking.resolve_mask("HG38", "hg38") == str(path)
    with pytest.raises(ValueError, match="assembly hg19"):
        masking.resolve_mask("hg38", "hg19")
    rows = [line.split("\t") for line in gzip.decompress(path.read_bytes()).decode().splitlines()]
    counts: dict[str, int] = {}
    for r in rows:
        counts[r[3]] = counts.get(r[3], 0) + 1
    assert counts == EXPECTED_HG38_LABELS
    assert {r[0] for r in rows} == {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
    repo_copy = REPO / "data" / "masks" / "hg38_mask.bed.gz"
    if repo_copy.is_file():
        assert (hashlib.sha256(repo_copy.read_bytes()).hexdigest()
                == hashlib.sha256(path.read_bytes()).hexdigest())


PACKAGES = ["GM12878_25kb", "GM12878_10kb", "IMR90_25kb", "IMR90_10kb",
            "HepG2_25kb", "HepG2_10kb"]


@pytest.mark.parametrize("pkg", PACKAGES)
def test_bundled_mask_reproduces_the_v11_tables(pkg):
    """apply_mask(v1.0 table, hg38) == v1.1 table, cell for cell, on every row v1.1 kept.
    v1.1 additionally deleted the rows whose position lies beyond the chromosome end;
    that deletion is not part of the mask and is checked against the v1.1 manifest."""
    v10 = REPO / "data" / "multicell" / pkg / "annotation.tsv.gz"
    v11 = REPO / "data" / "multicell_v1.1" / pkg / "annotation.tsv.gz"
    if not (v10.is_file() and v11.is_file()):
        pytest.skip(f"{pkg}: v1.0 or v1.1 table not in this checkout")
    read = dict(sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    old, ref = pd.read_csv(v10, **read), pd.read_csv(v11, **read)
    man = json.loads((v11.parent / "manifest.json").read_text())["masking"]
    masked, rec = masking.apply_mask(old, masking.load_mask(masking.resolve_mask("hg38")))
    key_ref = set(zip(ref["chrom"], ref["pos"]))
    kept = np.array([k in key_ref for k in zip(masked["chrom"], masked["pos"])])
    assert int((~kept).sum()) == man["n_rows_removed_off_chromosome_end"]
    masked = masked[kept].reset_index(drop=True)
    assert list(masked.columns) == list(ref.columns)
    pd.testing.assert_frame_equal(masked, ref)
    on_kept = rec[[k in key_ref for k in zip(rec["chrom"], rec["pos"])]]
    assert len(on_kept) == man["n_rows_masked"]


# -- command line and HTTP ----------------------------------------------------------------


def test_cli_mask_option(tmp_path, capsys):
    from tad_vci.cli import main
    inp = _inputs(tmp_path)
    args = ["annotate", "--boundaries", inp["primary"], "--resolution", str(RES),
            "--out", str(tmp_path / "cli"), "--tag", "t"]
    for path, name in inp["panel"]:
        args += ["--panel", f"{name}={path}"]
    with pytest.raises(SystemExit, match="mask BED not found"):
        main(args + ["--mask", str(tmp_path / "missing.bed")])
    with pytest.raises(SystemExit, match="assembly hg19"):
        main(args + ["--assembly", "hg19", "--mask", "hg38"])
    capsys.readouterr()
    assert main(args + ["--mask", inp["mask"], "--json"]) == 0
    text = capsys.readouterr().out
    out = json.loads(text[text.index("\n{\n") + 1:])   # after the engine's log lines
    assert out["n_masked"] == len(INSIDE)
    assert Path(out["mask_records"]).is_file()


def _post(body: dict) -> tuple[int, dict]:
    pytest.importorskip("aiohttp")
    from aiohttp.test_utils import make_mocked_request
    from tad_vci import api_routes

    class _Manager:
        def register(self, path, file_type, prediction_path=None):
            return f"tok-{path}"

        def get(self, token):
            return None

    req = make_mocked_request("POST", "/api/tadvci/annotate", app={"manager": _Manager()})

    async def _j():
        return body
    req.json = _j
    resp = asyncio.run(api_routes.handle_tadvci_annotate(req))
    payload = json.loads(resp.body.decode())
    api_routes.IMPORTED_ANNOT.pop(payload.get("session_token"), None)
    return resp.status, payload


def test_http_mask_bed(tmp_path):
    inp = _inputs(tmp_path)
    base = {"primary_bed": inp["primary"], "resolution": RES, "hash_tag": "http",
            "output_dir": str(tmp_path / "http"),
            "method_beds": [{"path": p, "name": n} for p, n in inp["panel"]]}
    status, payload = _post({**base, "mask_bed": str(tmp_path / "missing.bed")})
    assert status == 404 and "mask_bed not found" in payload["error"]
    status, payload = _post({**base, "mask_bed": 7})
    assert status == 400
    status, payload = _post({**base, "mask_bed": "hg38", "assembly": "hg19"})
    assert status == 400 and "hg19" in payload["error"]
    status, payload = _post({**base, "mask_bed": inp["mask"]})
    assert status == 200, payload
    assert payload["n_masked"] == len(INSIDE)
    status, plain = _post(base)
    assert status == 200 and "n_masked" not in plain
