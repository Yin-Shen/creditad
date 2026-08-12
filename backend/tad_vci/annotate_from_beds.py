"""Annotate already-called TAD boundaries from user-supplied BED files.

CrediTAD product path: boundaries are *inputs*, not re-called by built-in
pure-Python detectors. BD1 = multi-caller votes among the provided method panel
(n_methods = panel size). BD2–BD4 only if the user supplies ChIP / loop files —
never silently filled from demo data.
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from tad_vci import annotate
from tad_vci.graded_evidence import load_bands, load_null_quantiles


def _read_boundary_starts(bed_path: str) -> pd.DataFrame:
    """Load chrom/start from a BED-like file (whitespace, # comments)."""
    if not bed_path or not os.path.isfile(bed_path):
        raise FileNotFoundError(f"BED not found: {bed_path}")
    df = pd.read_csv(
        bed_path, sep=r"\s+", header=None, comment="#",
        usecols=[0, 1], names=["chrom", "start"], engine="python",
    )
    df = df[pd.to_numeric(df["start"], errors="coerce").notna()].copy()
    df["chrom"] = df["chrom"].astype(str).str.replace(r"^chr", "", regex=True)
    df["start"] = df["start"].astype(int)
    return df


def _method_positions(method_beds: list[tuple[str, str]]) -> dict[str, dict[str, np.ndarray]]:
    """method_name -> {chrom: sorted start array}."""
    out: dict[str, dict[str, list[int]]] = {}
    for path, name in method_beds:
        label = (name or Path(path).stem).strip() or Path(path).stem
        try:
            df = _read_boundary_starts(path)
        except Exception as exc:
            print(f"[annotate] method BED parse failed {path}: {exc}")
            continue
        bucket = out.setdefault(label, {})
        for cnum, g in df.groupby("chrom"):
            bucket.setdefault(str(cnum), []).extend(g["start"].tolist())
    return {
        m: {c: np.sort(np.array(v, dtype=int)) for c, v in byc.items()}
        for m, byc in out.items()
    }


def _count_votes(
    chrom: str,
    pos: int,
    method_pos: dict[str, dict[str, np.ndarray]],
    tol: int,
) -> tuple[int, list[str]]:
    """(n_votes, supporter names). Thin wrapper over _vote_details."""
    _, supporters, _ = _vote_details(chrom, pos, method_pos, tol)
    return len(supporters), supporters


def _vote_details(
    chrom: str,
    pos: int,
    method_pos: dict[str, dict[str, np.ndarray]],
    tol: int,
) -> tuple[int, list[str], list[int]]:
    """(n_votes, supporter names, SIGNED offset of each supporter's nearest call in bp).

    A method votes when ANY of its calls lies within `tol` of the candidate, so the
    supporting call is generally NOT at the candidate position. The offset is the
    per-supporter provenance of that vote and is written to the annotation
    (`supporting_offsets_bp`) so a card can state which call was borrowed and from how far.
    Consequence worth stating once here: because the candidate set is the union of all
    callers' calls, one call inside a 2-bin window supports every candidate in that window,
    so votes at neighbouring candidates are correlated rather than independent.
    """
    supporters: list[str] = []
    offsets: list[int] = []
    for m, byc in method_pos.items():
        arr = byc.get(chrom)
        if arr is None or len(arr) == 0:
            continue
        i = int(np.searchsorted(arr, pos))
        best: int | None = None
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(arr):
                d = int(arr[j]) - pos
                if abs(d) <= tol and (best is None or abs(d) < abs(best)):
                    best = d
        if best is not None:
            supporters.append(m)
            offsets.append(best)
    return len(supporters), supporters, offsets


def annotate_from_beds(
    primary_bed: str,
    method_beds: Optional[list[tuple[str, str]]] = None,
    *,
    resolution: int = 25000,
    tol: Optional[int] = None,
    loop_tol: int = 50000,
    ctcf_bw: Optional[str] = None,
    rad21_bw: Optional[str] = None,
    loops_path: Optional[str] = None,
    cell_line: Optional[str] = None,
    assembly: Optional[str] = None,
    hash_tag: str = "annotate",
    output_dir: str = "./tad_results",
    origin: str = "user",
    uncalibrated_bands: str = "in_situ",
) -> dict:
    """Grade primary BED candidates with optional multi-method BD1 panel.

    Parameters
    ----------
    primary_bed
        Boundary positions to grade (e.g. the union of several callers, or one caller).
    origin
        "user" for a real upload, "builtin_demo" for a package shipped with the app. Only a
        real upload may be labelled as the reviewer's own boundary in the evidence card.
    method_beds
        Optional list of ``(path, display_name)``. Each file is one BD1 voter.
        When empty, every primary position gets votes=1 and n_methods=1
        (single-caller; cross-caller support is not multi-method informative).
    """
    os.makedirs(output_dir, exist_ok=True)
    t0 = time.time()
    primary = _read_boundary_starts(primary_bed)
    if primary.empty:
        raise ValueError(f"Primary BED has no parseable positions: {primary_bed}")

    # BD1 vote-matching tolerance. A hardcoded 50 kb default made the API and the CLI
    # disagree on the same input: at 10 kb the CLI passed resolution*2 = 20 kb while the API
    # used 50 kb, so the two paths counted different votes and produced different tiers.
    # Derive it from the resolution in ONE place so every entry point matches by construction.
    #
    # 2026-07-29: default tightened from 2 bins to 1 bin. Two bins let one caller call
    # support up to 5 candidates (measured reuse 2.14-2.44 calls per candidate), and
    # 49.9-59.8% of multi-vote candidates had at most one call in their own bin. One bin
    # still absorbs the inter-caller registration offset that is actually documented in
    # these data (DI sits +1 bin from TopDom/SpectralTAD in three of the six packages),
    # which is why this is not taken to zero. Migration: docs/MEASURED_LIMITATIONS.md §11.
    if tol is None:
        tol = resolution * 1

    method_beds = list(method_beds or [])
    # Deduplicate by path while preserving order
    seen_paths: set[str] = set()
    unique_methods: list[tuple[str, str]] = []
    for path, name in method_beds:
        ap = os.path.abspath(path)
        if ap in seen_paths or not os.path.isfile(path):
            continue
        # Skip primary if user also listed it as a method (avoid double-counting)
        if os.path.abspath(path) == os.path.abspath(primary_bed):
            continue
        seen_paths.add(ap)
        unique_methods.append((path, name))

    method_pos = _method_positions(unique_methods)
    panel_names = list(method_pos.keys())
    n_methods = max(1, len(panel_names))

    # BD1 is a FRACTION of the panel, so the denominator has to be the panel that could
    # actually vote on this chromosome. A caller that produced no boundaries at all on a
    # chromosome (e.g. Arrowhead where KR normalisation fails to converge) can never
    # contribute a vote there; counting it in the denominator silently deflates every
    # grade on that chromosome. Measured impact before this fix: 7,321 candidates across
    # 4 chromosomes in the six shipped packages. `callers_absent_on_chrom` records who was
    # excluded so the deflation is visible rather than inferred.
    chrom_panel: dict[str, list[str]] = {}
    for m, byc in method_pos.items():
        for ch, arr in byc.items():
            if arr is not None and len(arr):
                chrom_panel.setdefault(str(ch), []).append(m)

    rows = []
    for _, r in primary.iterrows():
        cnum = str(r["chrom"])
        pos = int(r["start"])
        present = chrom_panel.get(cnum, [])
        absent = [m for m in panel_names if m not in present]
        n_eff = max(1, len(present)) if panel_names else 1
        if panel_names:
            votes, supporters, offsets = _vote_details(cnum, pos, method_pos, tol)
            # Zero panel support stays votes=0 (honest): grade_votes floors at
            # "weak" and never returns "none", so BD1 remains assessable without
            # fabricating a self-vote for a boundary no panel method supports.
        else:
            # No voting panel = single-caller upload. Keep the one primary call as
            # votes=1, but n_methods=1 makes grade_votes floor BD1 to "weak":
            # cross-caller agreement cannot be assessed from a single method.
            votes, supporters, offsets = 1, [], []
        rows.append({
            "chrom": cnum,
            "pos": pos,
            "votes": int(votes),
            # n_methods stays the declared panel size (provenance); n_methods_effective is
            # what BD1 is graded against on this chromosome.
            "n_methods": int(n_methods),
            "n_methods_effective": int(n_eff),
            "callers_absent_on_chrom": ";".join(absent) if absent else "",
            "supporting_callers": ";".join(supporters) if supporters else "",
            # Signed bp offset of the call each supporter voted with, same order as
            # supporting_callers. 0 means the supporter called this exact position.
            "supporting_offsets_bp": ";".join(str(o) for o in offsets) if offsets else "",
            "user_bed": 1 if origin == "user" else 0,
        })
    df = pd.DataFrame(rows)

    # Optional BD4. Its proximity window is NOT the BD1 vote tolerance: the criterion is
    # published as "loop anchors within 50 kb" and the mcool path has always used a fixed
    # 50 kb. Sharing one `tol` silently scaled it with the resolution, so the 10 kb packages
    # graded BD4 at +/-20 kb while the card and BD_SPEC both said 50 kb. Separated
    # 2026-07-29 so the BD1 tolerance can be tuned without moving BD4 underneath it.
    if len(df) and loops_path and os.path.isfile(loops_path):
        from mcool_bed import _parse_loop_anchors
        anchors = _parse_loop_anchors(loops_path)

        def _lc(c, p):
            k = str(c).replace("chr", "")
            return int(np.sum(np.abs(anchors[k] - p) <= loop_tol)) if k in anchors else 0

        df["loop_anchor_count"] = [_lc(c, p) for c, p in zip(df["chrom"], df["pos"])]

    # A track from a different genome build than the boundaries would be sampled at the
    # wrong coordinates and still return numbers. Refuse rather than grade nonsense.
    if assembly:
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from validate_bigwig import check_assembly
        for _p in (ctcf_bw, rad21_bw):
            if _p and os.path.isfile(_p):
                _msg = check_assembly(_p, assembly)
                if _msg:
                    raise ValueError(_msg)

    # Optional BD2/BD3 — only if user supplied paths (never demo auto-fill)
    bg_ctcf = bg_rad21 = None
    ctcf_bands, rad21_bands, bands_calibrated, bands_cell = load_bands(cell_line)
    bands_source = "packaged" if bands_calibrated else uncalibrated_bands
    # Where each observed fold sits in the random-window null. Packaged for a calibrated
    # cell, sampled in situ otherwise; empty when neither is available (the card then
    # reports the fold and the cutpoints without a percentile rather than inventing one).
    null_quantiles: dict = dict(load_null_quantiles(cell_line) or {}) if bands_calibrated else {}
    if len(df):
        if ctcf_bw and os.path.isfile(ctcf_bw):
            from mcool_bed import _bw_window_signal
            df["ctcf"], bg_ctcf = _bw_window_signal(ctcf_bw, df, resolution)
        if rad21_bw and os.path.isfile(rad21_bw):
            from mcool_bed import _bw_window_signal
            df["rad21"], bg_rad21 = _bw_window_signal(rad21_bw, df, resolution)

        # UNCALIBRATED CELL LINE POLICY.
        # Borrowing GM12878 cutpoints while the numerator and background come from the
        # user's own track mixes two scales, and since GM12878 has the lowest cutpoints of
        # the calibrated set the error only ever inflates grades (measured on HepG2 chr1:
        # BD2 +17.4%, BD3 +31.2%, tier +3.5%, zero downgrades). So we do not do that
        # silently. Default: estimate this track's own p75/p90 the same way the packaged
        # bands were produced.
        if not bands_calibrated:
            if uncalibrated_bands == "in_situ":
                from tad_vci.insitu_bands import estimate_bands_and_null
                got_ctcf = got_rad21 = False
                if ctcf_bw and os.path.isfile(ctcf_bw) and "ctcf" in df:
                    est, est_bg, est_q = estimate_bands_and_null(ctcf_bw)
                    if est:
                        # cutpoints AND denominator from the same random-window pass
                        ctcf_bands, bg_ctcf, got_ctcf = est, est_bg, True
                        if est_q:
                            null_quantiles["CTCF"] = est_q
                if rad21_bw and os.path.isfile(rad21_bw) and "rad21" in df:
                    est, est_bg, est_q = estimate_bands_and_null(rad21_bw)
                    if est:
                        rad21_bands, bg_rad21, got_rad21 = est, est_bg, True
                        if est_q:
                            null_quantiles["RAD21"] = est_q
                # A track we could not characterise must not be graded against foreign
                # cutpoints: drop it so BD2/BD3 report not_assessable instead.
                if "ctcf" in df and not got_ctcf:
                    df = df.drop(columns=["ctcf"])
                if "rad21" in df and not got_rad21:
                    df = df.drop(columns=["rad21"])
                bands_source = "in_situ"
                bands_cell = f"in_situ:{cell_line or 'unspecified'}"
            elif uncalibrated_bands == "not_assessable":
                for _c in ("ctcf", "rad21"):
                    if _c in df:
                        df = df.drop(columns=[_c])
                bands_source = "not_assessable"
            # "gm12878_fallback" keeps the historical behaviour; it must be asked for.
        if bands_calibrated and bands_source == "packaged":
            try:
                import json as _json
                from tad_vci.graded_evidence import _BANDS_JSON
                _cal = _json.loads(_BANDS_JSON.read_text())
                if bands_cell in _cal:
                    if bg_ctcf is not None and "ctcf" in df:
                        bg_ctcf = float(_cal[bands_cell]["CTCF"]["mean"])
                    if bg_rad21 is not None and "rad21" in df:
                        bg_rad21 = float(_cal[bands_cell]["RAD21"]["mean"])
            except Exception as _bg_err:
                print(f"[annotate] canonical-background override skipped: {_bg_err}")

        full = annotate(
            df,
            bg_ctcf=bg_ctcf or 1.0,
            bg_rad21=bg_rad21 or 1.0,
            ctcf_bands=ctcf_bands,
            rad21_bands=rad21_bands,
            n_methods=n_methods,
            null_quantiles=null_quantiles,
        )
        full["bands_cell"] = bands_cell
        full["bands_source"] = bands_source
        full["vote_tol_bp"] = int(tol)
        full["loop_window_bp"] = int(loop_tol)
        # BD1 vote model, stated per row: candidates here are the raw supplied positions
        # (typically the union of every caller's calls), so ONE caller call inside the
        # tolerance window supports EVERY candidate in that window. Votes at neighbouring
        # candidates are therefore correlated. The mcool built-in path clusters first
        # (bd1_vote_model=cluster_consensus) and does not have this property.
        full["bd1_vote_model"] = "nearest_within_tol"
        full["bands_calibrated"] = bool(bands_calibrated)
        full["data_cell_line"] = str(cell_line or "")
        full["assembly"] = str(assembly or "")
        full["tier_rule_version"] = "max_support_v2"
        full["chip_window_bp"] = 50000
        full["resolution_bp"] = int(resolution)
        full["upstream"] = "user_supplied_beds"
        full["bd1_panel"] = ",".join(panel_names) if panel_names else Path(primary_bed).stem
        full["bd1_rule"] = (
            f"grade_votes(votes, n_methods={n_methods}): "
            + (
                "single-method panel (n=1) -> BD1 floored to weak "
                "(cross-caller support not assessable)"
                if n_methods < 2
                else "strong>=ceil(0.8*n), moderate>=ceil(0.6*n), else weak"
            )
            + f"; a method votes if any of its calls is within +/-{int(tol)} bp "
              f"(+/-{int(tol) // int(resolution)} bins) of the candidate"
        )
        full["bg_ctcf"] = bg_ctcf if bg_ctcf is not None else np.nan
        full["bg_rad21"] = bg_rad21 if bg_rad21 is not None else np.nan
        if not bands_calibrated and (ctcf_bw or rad21_bw):
            print(
                f"[annotate] BD2/BD3 bands NOT calibrated for '{cell_line}' "
                f"-> bands estimated in situ from the supplied track (bands_source=in_situ)"
            )
    else:
        full = df

    suffix = f"_{hash_tag}" if hash_tag else ""
    ann_tsv = os.path.join(output_dir, f"tadvci_annotation_beds{suffix}.tsv")
    full.to_csv(ann_tsv, sep="\t", index=False)
    tier_bed = os.path.join(output_dir, f"tadvci_tier_beds{suffix}.bed")
    with open(tier_bed, "w") as fh:
        if len(full):
            for r in full.itertuples(index=False):
                fh.write(f"{r.chrom}\t{int(r.pos)}\t{int(r.pos) + int(resolution)}\t{r.D_tier}\n")

    dist = {}
    if len(full) and "D_tier" in full.columns:
        dist = (
            full["D_tier"].value_counts()
            .reindex(["D5", "D4", "D3", "D2", "D1"]).fillna(0).astype(int).to_dict()
        )
    print(
        f"[annotate-from-beds] n={len(full)} primary={primary_bed} "
        f"panel={panel_names or ['(single)']} n_methods={n_methods} "
        f"dist={dist} {time.time() - t0:.1f}s -> {ann_tsv}"
    )
    return {
        "TADVCI_annotation": ann_tsv,
        "TADVCI_tier": tier_bed,
        "n_boundaries": int(len(full)),
        "n_methods": int(n_methods),
        "bd1_panel": panel_names,
        "tier_dist": dist,
        "primary_bed": primary_bed,
    }
