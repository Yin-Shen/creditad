"""TAD-VCI API routes, mounted into the main CrediTAD backend (main.py / :5001).

Exposes evidence-card annotations and application-level versioned decision
records (tad_vci.graded_evidence + tad_vci.curation) for the front end.

Routes:
  GET  /api/tadvci/annotation/{chrom}   -> boundaries + evidence cards for a chrom
  POST /api/tadvci/signout              -> append a versioned decision snapshot
  GET  /api/tadvci/history/{boundary_id}
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

try:
    from aiohttp import web
except ModuleNotFoundError as _exc:  # pragma: no cover - core-only install
    # This module is the HTTP/GUI surface. aiohttp is declared in the `app` extra, not in
    # the core dependencies, so a core `pip install creditad` used to fail here with an
    # opaque "No module named 'aiohttp'". Name the extra instead; the CLI and the library
    # annotation path do not need it (audit 2026-08-12).
    raise ModuleNotFoundError(
        "tad_vci.api_routes needs the optional HTTP/GUI dependencies. Install them with: "
        "pip install 'creditad[app]'   (the `creditad annotate` CLI and "
        "tad_vci.annotate_from_beds do not require them)."
    ) from _exc

from tad_vci.curation import (
    append_snapshot, evidence_card, load_history, make_snapshot, next_version, verify_snapshot,
)
from tad_vci.rebuild_bundle import (
    BundleEligibilityError,
    validate_gm12878_bundle,
)

_REBUILD_DIR = Path(__file__).resolve().parent.parent / "analyses" / "tad_vci_rebuild_20260710"
_CELLS_DIR = _REBUILD_DIR / "cells_ice_v2"
_RECORDS_NAME = "decision_records_v3.jsonl"


def _records_log() -> Path:
    """Where the local review records live.

    The sign-out log is the only audit artefact this tool produces, and it defaulted to
    the source directory — i.e. inside the install tree for a packaged build. Program
    Files is read-only for a normal user (sign-out then fails outright) and the NSIS
    uninstaller removes the install directory wholesale on every upgrade, taking the
    records with it. CREDITAD_RECORDS_DIR lets the launcher point it at a per-user
    directory; the filename never changes, so an existing log keeps working.
    """
    env = (os.environ.get("CREDITAD_RECORDS_DIR") or "").strip()
    if env:
        base = Path(env)
        try:
            base.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            print(f"[records] CREDITAD_RECORDS_DIR={env!r} unusable ({exc}); "
                  f"falling back to the module directory")
            return Path(__file__).resolve().parent / _RECORDS_NAME
        return base / _RECORDS_NAME
    return Path(__file__).resolve().parent / _RECORDS_NAME


_LOG = _records_log()
_ASSEMBLY = "hg19"

# The built-in view is enabled only when the canonical July GM12878 rebuild passes
# the shared bundle verifier. Legacy/raw/KR and incomplete multi-cell artifacts are
# deliberately not served. The mcool/ChIP live on this server (referenced by path,
# not shipped), so the Hi-C heatmap renders live for any chromosome.
# Raw Hi-C / ChIP inputs live on the host (referenced by path, not shipped). Derive the
# data root from the repo layout (.../TAD) or the CREDITAD_DATA_ROOT env var, so no absolute
# host path is baked into the source. Override CREDITAD_DATA_ROOT for a different deployment.
_HERE = Path(__file__).resolve()
# parents[3] is the tree above the repo on a source checkout. A packaging layout that
# puts this file nearer the drive root makes that index raise IndexError AT IMPORT, which
# takes the whole /api/tadvci/* route table down; fall back to the deepest parent instead.
_DEFAULT_DATA_ROOT = _HERE.parents[3] if len(_HERE.parents) > 3 else _HERE.parents[-1]
_DATA_ROOT = Path(os.environ.get("CREDITAD_DATA_ROOT", str(_DEFAULT_DATA_ROOT)))
CELL_REGISTRY: dict[str, dict] = {
    "GM12878": {"annotation": str(_CELLS_DIR / "GM12878_annotation.tsv"),
                "mcool": str(_DATA_ROOT / "GSE63525_GM12878_insitu_DpnII_combined.mcool"),
                "ctcf": str(_DATA_ROOT / "CTCF_ENCFF749HDD.bigWig"),
                "rad21": str(_DATA_ROOT / "RAD21_ENCFF000WCT.bigWig"),
                "assembly": "hg19", "resolution": 25000},
}
DEFAULT_CELL = "GM12878"
_cell_cache: dict[str, dict] = {}   # cell -> {chrom: [cards]}
_cell_cache_hash: dict[str, str] = {}

# Five built-in method tracks shown alongside the evidence tier. Eligible BEDs live in the
# ICE-sidecar output namespace and are written only by build_gm12878_annotations.py.
_TADVCI_METHODS = ["insulation", "topdom_like", "contact_contrast", "spectral_profile", "network_modularity"]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _manifest_identity_matches(metadata: dict, require_hash: bool = True) -> bool:
    try:
        path = Path(metadata["path"])
        stat = path.stat()
        if int(metadata.get("size_bytes", -1)) != int(stat.st_size):
            return False
        if "mtime_ns" in metadata and int(metadata["mtime_ns"]) != int(stat.st_mtime_ns):
            return False
        if require_hash:
            return bool(metadata.get("sha256")) and metadata["sha256"] == _sha256(path)
        return True
    except (KeyError, OSError, TypeError, ValueError):
        return False


def _cell_manifest_status(cell: str, reg: dict) -> tuple[bool, str]:
    """Fail-closed eligibility check for a built-in scientific annotation."""
    if cell != "GM12878":
        return False, "only the canonical GM12878 bundle is eligible"
    try:
        verified = validate_gm12878_bundle(_CELLS_DIR)
        manifest = verified["manifest"]
        if manifest["cell_line"] != cell:
            return False, "verified bundle cell line differs from registry"
        if manifest["assembly"] != reg["assembly"]:
            return False, "verified bundle assembly differs from registry"
        if int(manifest["resolution_bp"]) != int(reg["resolution"]):
            return False, "verified bundle resolution differs from registry"
        if verified["annotation_path"].resolve() != Path(reg["annotation"]).resolve():
            return False, "verified annotation path differs from registry"
    except (BundleEligibilityError, KeyError, OSError, TypeError, ValueError) as exc:
        return False, f"canonical bundle verification failed: {exc}"
    return True, "eligible"


def _method_bed(cell: str, method: str) -> str | None:
    manifest_path = _CELLS_DIR / f"{cell}_run_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
        metadata = manifest["outputs"]["method_beds"][method]
    except (OSError, KeyError, json.JSONDecodeError):
        return None
    return str(metadata["path"]) if _manifest_identity_matches(metadata, True) else None

# Per-import TAD-VCI tier annotations, keyed by the mcool token. Populated by
# main.py when a user runs the 5-method 'tadvci' detection on THEIR uploaded data,
# so the visualization + sign-out read the user's own boundaries, not the static demo.
IMPORTED_ANNOT: dict[str, str] = {}
# Parsed evidence-card cache, keyed by token -> {chrom: [cards]} and (separately) the raw
# DataFrame for triage. The annotation TSV (9k-12k rows) is parsed ONCE per token instead of
# on every chromosome switch (previously ~1.5s each), so chrom navigation is instant.
_imported_cards_cache: dict[str, dict] = {}
_imported_track_cache: dict[str, dict] = {}
_imported_df_cache: dict[str, object] = {}


def set_imported_annotation(token: str, annotation_tsv: str) -> None:
    if token and annotation_tsv:
        IMPORTED_ANNOT[token] = annotation_tsv
        _imported_cards_cache.pop(token, None)   # invalidate stale parse if the token is re-pointed
        _imported_track_cache.pop(token, None)
        _imported_df_cache.pop(token, None)


_KNOWN_ASSEMBLIES = {"hg19", "grch37", "hg38", "grch38"}
_ASSEMBLY_ALIASES = {"grch37": "hg19", "grch38": "hg38"}


def _norm_assembly(val) -> str | None:
    """Normalise a user/registry assembly string to hg19/hg38 (GRCh37/GRCh38 aliases too).
    Returns None for anything unrecognised so callers can fall through to the next source."""
    if not isinstance(val, str):
        return None
    v = val.strip().lower()
    if v in _ASSEMBLY_ALIASES:
        return _ASSEMBLY_ALIASES[v]
    return v if v in _KNOWN_ASSEMBLIES else None


def _assembly_for_boundary(token, chrom: str, bid: str) -> str | None:
    """Read the assembly the user declared at import time from the token's annotation TSV.

    Authoritative provenance for a live upload: prefer the exact boundary's row (so a
    multi-assembly token could not happen, but stay precise), else any non-null assembly value
    in that token's table. Returns None if no token, no annotation, or no assembly column —
    so the sign-out path can fall back to the payload / registry without ever defaulting an
    hg38 upload to hg19."""
    if not token:
        return None
    df = _df_for_token(token)
    if df is None or "assembly" not in df.columns:
        return None
    ck = str(chrom).replace("chr", "")
    try:
        pos = int(bid.split(":")[1])
        row = df[(df["chrom"].astype(str).str.replace("chr", "") == ck) & (df["pos"].astype("Int64") == pos)]
        if len(row):
            a = _norm_assembly(row.iloc[0]["assembly"])
            if a:
                return a
    except (IndexError, ValueError, KeyError):
        pass
    vals = df["assembly"].dropna()
    return _norm_assembly(vals.iloc[0]) if len(vals) else None

# Fields the evidence-tier TRACK needs to paint a tick and to identify what was clicked.
# Everything else on a card is only ever read for ONE boundary at a time, in the panel.
_TRACK_FIELDS = ("boundary_id", "pos", "auto_tier", "votes", "user_bed")


def _track_rows_from_df(df) -> dict:
    """Slim per-chromosome tick lists, built straight from the table.

    Serving a full evidence card per boundary made the tier track cost 21.2 MB and 22.5 s
    for one chromosome of a 10 kb package (5,298 boundaries x 3.9 kB), and held ~250 MB of
    card dicts per token in memory — to draw ticks that need five fields. The full card is
    built on demand for the one boundary a reviewer clicks (`/api/tadvci/card/...`).
    """
    df["chrom"] = df["chrom"].astype(str)
    sub = df.sort_values(["chrom", "pos"])
    chroms = sub["chrom"].str.replace("chr", "", regex=False).to_numpy()
    poss = sub["pos"].astype(int).to_numpy()
    tiers = (sub["D_tier"] if "D_tier" in sub.columns else pd.Series([None] * len(sub))).to_numpy()
    votes = (sub["votes"] if "votes" in sub.columns else pd.Series([None] * len(sub))).to_numpy()
    ub = (sub["user_bed"] if "user_bed" in sub.columns else pd.Series([0] * len(sub))).to_numpy()
    by_chrom: dict[str, list] = {}
    for c, p, t, v, u in zip(chroms, poss, tiers, votes, ub):
        by_chrom.setdefault(str(c), []).append({
            "boundary_id": f"{c}:{int(p)}",
            "pos": int(p),
            "auto_tier": None if t is None or t != t else str(t),
            "votes": None if v is None or v != v else int(v),
            "user_bed": bool(u) if u == u and u is not None else False,
        })
    return by_chrom


def _cards_from_df(df) -> dict:
    """Full evidence cards per chromosome. Only the `view=full` escape hatch uses this."""
    df["chrom"] = df["chrom"].astype(str)
    by_chrom: dict[str, list] = {}
    for _, r in df.sort_values(["chrom", "pos"]).iterrows():
        card = evidence_card(r.to_dict())
        votes = int(r["votes"]) if "votes" in r.index and pd.notna(r["votes"]) else card.get("votes")
        by_chrom.setdefault(str(r["chrom"]).replace("chr", ""), []).append(
            {**card, "pos": int(r["pos"]), "votes": votes, "_row": {k: r[k] for k in r.index}})
    return by_chrom


def _cards_for(chrom: str, cell: str = DEFAULT_CELL) -> list[dict]:
    """Evidence cards for one chromosome of one built-in cell line (cached per cell)."""
    chrom = str(chrom).replace("chr", "")
    cell = cell if cell in CELL_REGISTRY else DEFAULT_CELL
    eligible, _ = _cell_manifest_status(cell, CELL_REGISTRY[cell])
    if not eligible:
        _cell_cache.pop(cell, None)
        _cell_cache_hash.pop(cell, None)
        return []
    manifest = json.loads((_CELLS_DIR / f"{cell}_run_manifest.json").read_text())
    annotation_hash = manifest["outputs"]["annotation"]["sha256"]
    if cell not in _cell_cache or _cell_cache_hash.get(cell) != annotation_hash:
        ann = Path(CELL_REGISTRY[cell]["annotation"])
        _cell_cache[cell] = _cards_from_df(pd.read_csv(ann, sep="\t"))
        _cell_cache_hash[cell] = annotation_hash
    return _cell_cache[cell].get(chrom, [])


def _cards_for_token(token: str, chrom: str) -> list[dict]:
    if token not in _imported_cards_cache:
        path = IMPORTED_ANNOT.get(token)
        if not path or not Path(path).exists():
            return []   # don't cache a miss (token may become valid after re-registration)
        _imported_cards_cache[token] = _cards_from_df(pd.read_csv(path, sep="\t"))
    return _imported_cards_cache[token].get(str(chrom).replace("chr", ""), [])


def _df_for_token(token: str):
    """The token's annotation table, parsed once. Shared by the track, the single-card
    lookup and the assembly lookup so one token holds ONE copy of its data."""
    df = _imported_df_cache.get(token)
    if df is None:
        path = IMPORTED_ANNOT.get(token)
        if not path or not Path(path).exists():
            return None
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception as exc:
            print(f"[tadvci] annotation parse failed for {token}: {exc}")
            return None
        _imported_df_cache[token] = df
    return df


def _track_for_token(token: str, chrom: str) -> list[dict]:
    if token not in _imported_track_cache:
        df = _df_for_token(token)
        if df is None:
            return []   # don't cache a miss (token may become valid after re-registration)
        _imported_track_cache[token] = _track_rows_from_df(df)
    return _imported_track_cache[token].get(str(chrom).replace("chr", ""), [])


def _card_for_position(token: str, chrom: str, pos: int) -> dict | None:
    """One full evidence card, built on demand from the cached table."""
    df = _df_for_token(token)
    if df is None:
        return None
    ck = str(chrom).replace("chr", "")
    rows = df[(df["chrom"].astype(str).str.replace("chr", "", regex=False) == ck)
              & (df["pos"].astype("int64") == int(pos))]
    if not len(rows):
        return None
    r = rows.iloc[0]
    card = evidence_card(r.to_dict())
    votes = int(r["votes"]) if "votes" in r.index and pd.notna(r["votes"]) else card.get("votes")
    return {**card, "pos": int(r["pos"]), "votes": votes}


async def handle_tadvci_annotation(request: web.Request) -> web.Response:
    chrom = request.match_info["chrom"]
    cell = request.query.get("cell", DEFAULT_CELL)
    if cell not in CELL_REGISTRY:
        cell = DEFAULT_CELL
    reg = CELL_REGISTRY[cell]
    eligible, reason = _cell_manifest_status(cell, reg)
    if not eligible:
        return web.json_response(
            {"error": f"built-in annotation is not eligible: {reason}", "cell_line": cell},
            status=503,
        )
    cards = _cards_for(chrom, cell)
    public = [{k: v for k, v in c.items() if k != "_row"} for c in cards]
    return web.json_response({"chrom": str(chrom).replace("chr", ""), "cell_line": cell,
                              "assembly": reg["assembly"], "resolution": reg["resolution"],
                              "n": len(public), "boundaries": public})


async def handle_tadvci_cells(request: web.Request) -> web.Response:
    """List manifest-gated built-in cell views and their current availability."""
    cells = []
    for cell, reg in CELL_REGISTRY.items():
        eligible, reason = _cell_manifest_status(cell, reg)
        mcool_ok = os.path.isfile(reg["mcool"])
        cells.append({"cell_line": cell, "assembly": reg["assembly"],
                      "resolution": reg["resolution"],
                      "annotation_ready": eligible, "heatmap_ready": mcool_ok,
                      "eligibility_reason": reason,
                      "available": eligible and mcool_ok})
    return web.json_response({"cells": cells, "default": DEFAULT_CELL})


async def handle_tadvci_use_cell(request: web.Request) -> web.Response:
    """Register a built-in cell line's Hi-C mcool (+ CTCF/RAD21 BigWigs) so the heatmap +
    signal tracks render live, and tell the frontend which cell annotation to read."""
    cell = request.match_info["cell"]
    reg = CELL_REGISTRY.get(cell)
    if reg is None:
        return web.json_response({"error": f"unknown cell line {cell}"}, status=404)
    if not os.path.isfile(reg["mcool"]):
        return web.json_response(
            {"error": f"Hi-C matrix for {cell} not found on this server ({reg['mcool']})."},
            status=404)
    eligible, reason = _cell_manifest_status(cell, reg)
    if not eligible:
        return web.json_response(
            {"error": f"Annotation for {cell} is not eligible: {reason}."},
            status=503)
    mgr = request.app["manager"]
    hic_token = mgr.register(reg["mcool"], "hic")
    # Point this token at the manifest-gated cell annotation so the imported tier
    # and local review card use the same per-token path as user-provided data.
    set_imported_annotation(hic_token, reg["annotation"])
    bigwig_tokens = []
    for label, path in (("CTCF", reg.get("ctcf")), ("RAD21", reg.get("rad21"))):
        if path and os.path.isfile(path):
            bigwig_tokens.append({"token": mgr.register(path, "bigwig"), "fileName": label})
    # The 5 built-in method tracks (genome-wide per-method boundary BEDs).
    bed_tokens = []
    for method in _TADVCI_METHODS:
        mb = _method_bed(cell, method)
        if mb:
            bed_tokens.append({"token": mgr.register(mb, "bed"), "fileName": method})
    return web.json_response({
        "ok": True, "cell_line": cell, "assembly": reg["assembly"],
        "resolution": reg["resolution"], "hic_token": hic_token, "mcoolPath": reg["mcool"],
        "bigwig_tokens": bigwig_tokens, "bed_tokens": bed_tokens,
    })


# The single-cell TADShop consensus demo endpoint was removed: the UI now ships the
# 3-cell x 2-resolution panel demo below, which uses no consensus set. The archival
# website-native inputs stay on disk under example_data/tadshop_25kb (the upstream
# service is offline, so that data cannot be re-fetched) but are no longer served.

# ---------------------------------------------------------------------------
# Multi-cell demo: 3 cell lines x 2 resolutions. Candidates are the deduplicated
# UNION of the six panel callers' boundaries — no ConsensusTAD is used or shown.
# ---------------------------------------------------------------------------
# Repo layout: .../creditad/backend/tad_vci/api_routes.py → parents[2] = creditad root.
# Packaged Electron layout: resources/backend/tad_vci/... → parents[2] = resources/.
_CREDITAD_ROOT = Path(__file__).resolve().parents[2]
_MULTICELL_DIR = _CREDITAD_ROOT / "example_data" / "multicell"
_BACKEND_DIR = Path(__file__).resolve().parents[1]
_EXTRA_MODE = _BACKEND_DIR / "extra_mode"


def _dev_shared_raw_default() -> Path | None:
    """Development-host fallback, resolved by data_manager (single definition).

    Imported lazily so a data_manager import problem degrades the optional
    full-genome search instead of taking the whole TAD-VCI route table down.
    """
    try:
        from data_manager import dev_shared_raw_default
    except Exception:
        return None
    return dev_shared_raw_default()


def _shared_raw_roots() -> list[Path]:
    """Search order for optional large tracks (mcool / BigWig).

    Portable Windows builds do not ship multi-GB matrices. Set CREDITAD_SHARED_RAW
    to a folder the user populated (see Data Manager catalog), or rely on the
    development default when present.
    """
    roots: list[Path] = []
    env = (os.environ.get("CREDITAD_SHARED_RAW") or "").strip()
    if env:
        roots.append(Path(env))
    # Packaged data/ next to backend or resources
    cands = [_CREDITAD_ROOT / "data" / "raw", _BACKEND_DIR / "data" / "raw"]
    # The development host's tree is appended only when it really is that host
    # (POSIX + directory present). Baking it in unconditionally contradicted the
    # note above and let the packaged build resolve host files through wine's Z:
    # mapping, hiding the out-of-package branch from every smoke test.
    dev = _dev_shared_raw_default()
    if dev is not None:
        cands.append(dev)
    for cand in cands:
        if cand not in roots:
            roots.append(cand)
    return roots


def _first_existing(candidates: list[Path | str]) -> Path | None:
    for c in candidates:
        p = Path(c)
        try:
            if p.is_file() and p.stat().st_size > 0:
                return p
        except OSError:
            continue
    return None


def _multicell_track_paths(cell: str) -> dict:
    """Resolve mcool / CTCF / RAD21 for a multicell demo.

    Layers (first hit wins):
      1. CREDITAD_SHARED_RAW full-genome files (optional user download)
      2. SHIPPED chr7 subsets (Portable / example_data/chr7_tracks + extra_mode)

    Genome-wide *evidence* still comes from example_data/multicell/*/annotation.tsv
    for every cell × resolution. Browser tracks without full-genome downloads
    only have signal on chr7.
    """
    specs = {
        "GM12878": {
            "mcool_name": "4DNFIXP4QG5B.mcool",
            "ctcf_acc": "ENCFF734CUT",
            "rad21_acc": "ENCFF571ZJJ",
        },
        "IMR90": {
            "mcool_name": "4DNFIJTOIGOI.mcool",
            "ctcf_acc": "ENCFF105FHL",
            "rad21_acc": "ENCFF048PZI",
        },
        "HepG2": {
            "mcool_name": "4DNFIS6HAUPP.mcool",
            "ctcf_acc": "ENCFF357NFO",
            "rad21_acc": "ENCFF972ODZ",
        },
    }
    spec = specs[cell]
    key = cell.lower()
    mcool_cands: list[Path] = []
    ctcf_cands: list[Path] = []
    rad21_cands: list[Path] = []

    # (1) Optional full-genome under shared raw
    for root in _shared_raw_roots():
        mcool_cands.append(root / key / "hic" / spec["mcool_name"])
        ctcf_cands.append(root / key / "ctcf" / f"{spec['ctcf_acc']}.bigWig")
        rad21_cands.append(root / key / "rad21" / f"{spec['rad21_acc']}.bigWig")

    # (2) Shipped chr7-only demo tracks (names always contain _chr7)
    chr7_dirs = (
        _EXTRA_MODE / "chr7_tracks",
        _CREDITAD_ROOT / "example_data" / "chr7_tracks",
    )
    for d in chr7_dirs:
        ctcf_cands.append(
            d / f"{cell}_CTCF_{spec['ctcf_acc']}_hg38_chr7.bigWig"
        )
        rad21_cands.append(
            d / f"{cell}_RAD21_{spec['rad21_acc']}_hg38_chr7.bigWig"
        )
    # No cell-agnostic fallback. `CTCF_chr7.bigWig` / `RAD21_chr7.bigWig` used to be
    # appended here for every cell line; they were byte-identical copies of
    # GM12878's tracks, so a missing IMR90 or HepG2 file resolved silently to
    # GM12878's signal — the same failure mode as the shared GRCh37 matrix. A track
    # that is absent must stay absent; the UI reports that honestly.

    # (3) Shipped chr7 heatmap — PER CELL LINE.
    #
    # There used to be a single shared `chr7.mcool` here, used for whichever demo
    # was loaded. It was GM12878 on GRCh37, so the IMR90 and HepG2 demos drew
    # GM12878's contact map under their own boundaries, and all three drew GRCh38
    # ticks over a GRCh37 map. Both are silent errors — the picture looks fine.
    # The shared file is deliberately NOT a fallback: a demo with no map is
    # correct and is explained in the UI, a demo with someone else's map is not.
    for d in chr7_dirs:
        mcool_cands.append(d / f"{cell}_chr7_hg38.mcool")
    return {
        "mcool": _first_existing(mcool_cands),
        "ctcf": _first_existing(ctcf_cands),
        "rad21": _first_existing(rad21_cands),
    }


# Display order in the viewer: the six panel callers only. The union is still the
# candidate set behind annotation.tsv, but it is not shown as its own track.
_MULTICELL_BEDS = [
    "TopDom_boundaries.bed",
    "SpectralTAD_boundaries.bed",
    "OnTAD_boundaries.bed",
    "MSTD_boundaries.bed",
    "arrowhead_boundaries.bed",
    "DI_boundaries.bed",
]

_MULTICELL_CELLS = ("GM12878", "IMR90", "HepG2")


def _multicell_combo_dir(cell: str, resolution: int) -> Path:
    return _MULTICELL_DIR / f"{cell}_{int(resolution) // 1000}kb"


def _annotation_path(combo_dir: Path) -> Path | None:
    """The package's evidence table, plain or gzipped, or None if absent.

    The six tables are 139 MB of tab-separated text and compress ~8.8x. The desktop
    package ships them gzipped: the installer writes its payload twice (NSIS unpacks an
    inner archive, then decompresses it into the install directory), so every megabyte of
    shipped text is paid for twice on disk and twice on the progress bar. The repository
    keeps the plain .tsv, so the analysis scripts are untouched; only the packaged copy is
    compressed. pandas reads either transparently from the suffix."""
    plain = combo_dir / "annotation.tsv"
    if plain.is_file():
        return plain
    gz = combo_dir / "annotation.tsv.gz"
    return gz if gz.is_file() else None


def _open_annotation(path: Path):
    """Text handle for an annotation table, gzipped or not."""
    if path.suffix == ".gz":
        import gzip
        return gzip.open(path, "rt")
    return path.open()


# Lengths that separate the two builds. Any one match identifies the assembly of a
# matrix, which is what the browser-track gate below needs.
_ASSEMBLY_SIGNATURES = {
    "hg19": {"1": 249250621, "7": 159138663, "X": 155270560},
    "hg38": {"1": 248956422, "7": 159345973, "X": 156040895},
}


def _matrix_profile(path) -> dict:
    """Chromosomes and genome build of an .mcool, read from its metadata only.

    Returned assembly is `None` when the build could not be established; callers must
    treat that as "unverified" and refuse to draw the matrix, NOT as "matches". The
    difference is not theoretical: the signature table only knows chr1/chr7/chrX, so an
    hg19 matrix holding chr2 alone profiled as assembly=None and an unverified-equals-fine
    gate served it under an hg38 package. `assembly_source` records how it was decided and
    `unverified_reason` says why it could not be, so neither can be lost downstream.
    """
    out = {"chroms": [], "assembly": None, "assembly_source": None, "unverified_reason": None}
    if path is None:
        out["unverified_reason"] = "no matrix file"
        return out
    try:
        import cooler  # imported lazily: the heatmap is optional

        uris = cooler.fileops.list_coolers(str(path))
        if not uris:
            out["unverified_reason"] = "file contains no cooler group"
            return out
        c = cooler.Cooler(f"{str(path)}::{uris[0]}")
        sizes = {str(k).replace("chr", ""): int(v) for k, v in dict(c.chromsizes).items()}
        declared = _norm_assembly(c.info.get("genome-assembly"))
    except ImportError as exc:
        out["unverified_reason"] = f"cooler is unavailable ({exc})"
        return out
    except Exception as exc:
        out["unverified_reason"] = f"matrix metadata could not be read ({type(exc).__name__}: {exc})"
        return out
    out["chroms"] = sorted(sizes)
    for build, sig in _ASSEMBLY_SIGNATURES.items():
        if any(sizes.get(name) == length for name, length in sig.items()):
            out["assembly"] = build
            out["assembly_source"] = "chromosome_length"
            return out
    # Second source, weaker because it is declared rather than measured: cooler stores the
    # build in its own metadata. Used only when no length in the file identifies a build,
    # so a subset matrix (no chr1/7/X) is not rejected out of hand.
    if declared:
        out["assembly"] = declared
        out["assembly_source"] = "file_metadata"
        return out
    out["unverified_reason"] = (
        f"no chromosome length in this matrix identifies a genome build "
        f"(has {', '.join(out['chroms'][:6]) or 'no chromosomes'}) and it declares none"
    )
    return out


# hash_tag ends up inside output FILE NAMES, so anything that can act as a path
# separator or a parent reference must be rejected before it is concatenated.
_SAFE_TAG_RE = re.compile(r"[A-Za-z0-9._-]{1,64}")


def _resolved_output_dir(raw: str) -> str:
    """Normalise a caller-chosen results directory and make sure it is usable.

    Deliberately NOT confined to a fixed root: this is a local single-user tool and the
    caller (GUI, CLI and the parity test) legitimately names its own output directory.
    What is enforced is that the path is normalised once here and that an unusable one
    fails as a 400 now, rather than as a 500 deep inside the writer.
    """
    if not raw.strip():
        raise ValueError("output_dir must not be empty")
    path = os.path.abspath(os.path.expanduser(raw.strip()))
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as exc:
        raise ValueError(f"output_dir cannot be created: {exc.strerror}") from None
    if not os.access(path, os.W_OK):
        raise ValueError("output_dir is not writable")
    return path


def _parse_resolution(raw, default: int) -> int:
    """Bin size from a request body: absent means default, anything else must be a
    positive whole number.

    `int(raw or default)` treated 0 as "unspecified" and silently produced a 25 kb
    result, and truncated 10000.7 to 10000 — while -5 was correctly rejected, so the
    same field failed loudly or quietly depending on which wrong value was sent.
    Resolution is hashed into every sign-out record, so it may not be guessed.
    """
    if raw in (None, ""):
        return default
    if isinstance(raw, bool):
        raise ValueError("resolution must be an integer number of base pairs")
    if isinstance(raw, float) and not raw.is_integer():
        raise ValueError(f"resolution must be a whole number of base pairs, got {raw}")
    try:
        value = int(raw)
    except (TypeError, ValueError):
        raise ValueError("resolution must be an integer number of base pairs") from None
    if value <= 0:
        raise ValueError("resolution must be positive")
    return value


def _evidence_chroms(annotation_tsv: Path) -> list:
    """Chromosomes present in the evidence table (which is genome-wide)."""
    seen = set()
    try:
        with _open_annotation(annotation_tsv) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = header.index("chrom")
            for line in fh:
                parts = line.split("\t", idx + 1)
                if len(parts) > idx:
                    seen.add(parts[idx].replace("chr", ""))
    except (OSError, ValueError):
        return []

    def key(name):
        return (0, int(name)) if name.isdigit() else (1, name)

    return sorted(seen, key=key)


async def handle_tadvci_multicell_datasets(request: web.Request) -> web.Response:
    """List the available demo datasets (3 cell lines x 2 resolutions)."""
    out = []
    for cell in _MULTICELL_CELLS:
        for res in (25000, 10000):
            d = _multicell_combo_dir(cell, res)
            man = d / "manifest.json"
            if _annotation_path(d) is None or not man.is_file():
                continue
            try:
                m = json.loads(man.read_text())
            except json.JSONDecodeError:
                continue
            out.append({
                "cell_line": cell,
                "resolution": res,
                "assembly": m.get("assembly"),
                "n_candidates": m.get("n_candidates"),
                "bd1_n_methods": m.get("bd1_n_methods"),
                "per_caller_boundary_counts": m.get("per_caller_boundary_counts", {}),
                "tier_counts": m.get("tier_counts", {}),
                "reference_note": m.get("reference_note"),
                "caveats": m.get("caveats", []),
                "consensus_used": bool(m.get("consensus_used", False)),
            })
    # When nothing is found, say WHERE it looked. An installed desktop build that
    # resolves example_data to the wrong place is indistinguishable, from the client, from
    # a backend that has not finished starting — and the user only ever saw "Could not
    # read the shipped packages". The path makes it a one-line diagnosis.
    payload = {"ok": True, "datasets": out}
    if not out:
        payload["ok"] = False
        payload["searched"] = str(_MULTICELL_DIR)
        payload["searched_exists"] = _MULTICELL_DIR.is_dir()
        try:
            payload["searched_entries"] = sorted(
                p.name for p in _MULTICELL_DIR.iterdir())[:12] if _MULTICELL_DIR.is_dir() else []
        except OSError as exc:
            payload["searched_entries"] = [f"<unreadable: {exc}>"]
    return web.json_response(payload)


def _track_layer(path) -> str:
    """Which data layer one browser track came from."""
    return "bundled_chr7" if path is not None and "chr7" in str(path).lower() else "full_genome"


def _signal_tracks(mgr, tracks: dict, bigwig_tokens: list) -> list:
    """Per-assay layer + covered chromosomes for the registered signal tracks."""
    registered = {t["fileName"] for t in bigwig_tokens}
    by_name = {t["fileName"]: t["token"] for t in bigwig_tokens}
    out = []
    for label, key in (("CTCF", "ctcf"), ("RAD21", "rad21")):
        if label not in registered:
            continue
        p = tracks.get(key)
        chroms = []
        r = mgr.get(by_name[label])
        if r is not None:
            try:
                chroms = [c["name"] for c in r.get_all_chroms()]
            except Exception as exc:      # a header we cannot read is not a fatal scope error
                print(f"[scope] {label} chromosome list unavailable: {exc}")
        out.append({"name": label, "source": _track_layer(p), "chroms": chroms})
    return out


def _signal_chroms(mgr, bigwig_tokens: list) -> list:
    """Union of the chromosomes the registered signal tracks cover."""
    seen = set()
    for t in bigwig_tokens:
        r = mgr.get(t["token"])
        if r is None:
            continue
        try:
            seen.update(c["name"] for c in r.get_all_chroms())
        except Exception:
            continue

    def key(name):
        return (0, int(name)) if name.isdigit() else (1, name)

    return sorted(seen, key=key)


def _signal_layer(tracks: dict, bigwig_tokens: list) -> str | None:
    """Session-level layer for the signal tracks: one layer, "mixed", or None."""
    if not bigwig_tokens:
        return None
    registered = {t["fileName"] for t in bigwig_tokens}
    layers = {
        _track_layer(tracks.get(key))
        for label, key in (("CTCF", "ctcf"), ("RAD21", "rad21"))
        if label in registered
    }
    if not layers:
        return None
    return layers.pop() if len(layers) == 1 else "mixed"


async def handle_tadvci_use_multicell_example(request: web.Request) -> web.Response:
    """Load one of the 3-cell x 2-resolution demos (panel callers only, no consensus)."""
    try:
        body = await request.json()
    except Exception:
        body = {}
    # A JSON array or string reached body.get() and raised AttributeError, which the
    # middleware turned into an opaque 500 for what is plainly a bad request.
    if not isinstance(body, dict):
        return web.json_response({"error": "body must be a JSON object"}, status=400)
    cell = str(body.get("cell_line") or body.get("cell") or "GM12878")
    try:
        res = _parse_resolution(body.get("resolution"), 25000)
    except ValueError as exc:
        return web.json_response({"error": str(exc)}, status=400)
    if cell not in _MULTICELL_CELLS:
        return web.json_response({"error": f"unknown cell_line: {cell}"}, status=400)
    if res not in (25000, 10000):
        return web.json_response({"error": f"unsupported resolution: {res}"}, status=400)

    d = _multicell_combo_dir(cell, res)
    ann = _annotation_path(d)
    if ann is None:
        return web.json_response(
            {"error": f"demo package missing ({d}). Run scripts/build_multicell_examples.py."},
            status=404,
        )
    try:
        manifest = json.loads((d / "manifest.json").read_text())
    except (OSError, json.JSONDecodeError):
        manifest = {}

    tracks = _multicell_track_paths(cell)
    mcool = tracks.get("mcool")
    # Heatmap is optional: evidence cards live in annotation.tsv. Portable builds may
    # only ship chr7.mcool or no matrix until the user downloads full-genome files.
    notes = []
    mgr = request.app["manager"]
    hic_token = None
    pkg_assembly = manifest.get("assembly")
    profile = _matrix_profile(mcool)
    map_blocked = None
    if mcool is not None:
        # A matrix on the wrong build must never be drawn under this package's
        # boundaries: chr7 alone differs by 207,310 bp between hg19 and hg38 and the
        # local offset varies along the chromosome, so every tick would sit at the
        # wrong place while still looking plausible.
        #
        # Fail CLOSED on an unverifiable build too. The old test was
        # `profile["assembly"] and ... != pkg_assembly`, so "could not tell" passed as
        # "matches" — an hg19 matrix whose build could not be read was drawn under hg38
        # boundaries with `assembly: null` and no reason given anywhere.
        if profile["assembly"] is None:
            map_blocked = (
                f"the genome build of the available matrix could not be verified "
                f"({profile['unverified_reason']})"
            )
            notes.append(
                f"Contact map not shown: {map_blocked}. Evidence cards are unaffected. "
                "A matrix whose build cannot be confirmed is not drawn, because ticks "
                "from the wrong build look plausible. Add a .mcool that carries chr1, "
                "chr7 or chrX, or that declares its assembly, via Data Manager → "
                "Dataset setup guide."
            )
            mcool = None
        elif pkg_assembly and profile["assembly"] != pkg_assembly:
            map_blocked = (
                f"the available matrix is {profile['assembly']} and this package is "
                f"{pkg_assembly}"
            )
            notes.append(
                f"Contact map not shown: {map_blocked}. Evidence cards are unaffected. "
                "Add a matching full-genome .mcool via Data Manager → Dataset setup guide."
            )
            mcool = None
        else:
            hic_token = mgr.register(str(mcool), "hic")
            if "chr7" in Path(mcool).name.lower():
                notes.append(
                    "Using the packaged chr7 heatmap subset; place a full-genome .mcool "
                    "under CREDITAD_SHARED_RAW for genome-wide Hi-C."
                )
    else:
        notes.append(
            "No Hi-C .mcool found — evidence cards still load; heatmap empty. "
            "See Data Manager → Dataset setup guide."
        )
    if tracks.get("ctcf") and "chr7" in str(tracks["ctcf"]):
        notes.append("CTCF/RAD21 tracks are packaged chr7-only subsets for demo browsing.")

    # GlobalManager.register() de-duplicates by path, and one cell line's 25 kb and 10 kb
    # demos share a single mcool — so both packages would receive the SAME hic_token and
    # set_imported_annotation() would re-point the earlier session's tier table at the other
    # resolution's annotation.tsv. Key the annotation on a per-package token instead; the
    # viewer prefers annotation_token over hic_token when picking the tier source.
    annotation_token = f"multicell:{d.name}"
    set_imported_annotation(annotation_token, str(ann))

    bigwig_tokens = []
    for label, key in (("CTCF", "ctcf"), ("RAD21", "rad21")):
        p = tracks.get(key)
        if p is not None and Path(p).is_file():
            try:
                tok = mgr.register(str(p), "bigwig")
                bigwig_tokens.append({"token": tok, "fileName": label})
            except Exception as reg_err:
                notes.append(f"{label} BigWig register failed: {reg_err}")
        else:
            notes.append(f"{label} BigWig not found (optional).")

    bed_tokens = []
    for name in _MULTICELL_BEDS:
        bp = d / name
        if bp.is_file() and bp.stat().st_size > 0:
            label = name.replace("_boundaries.bed", "")
            bed_tokens.append({"token": mgr.register(str(bp), "bed"), "fileName": label})

    base_note = (
        f"{cell} @ {res // 1000} kb. Candidates are the deduplicated union of six caller "
        "boundary sets; BD1 counts how many of the six support each position. "
        "CrediTAD does not call TADs here, and no consensus set is used."
    )
    if notes:
        base_note = base_note + " " + " ".join(notes)

    return web.json_response({
        "ok": True,
        "demo": d.name,
        "cell_line": cell,
        "assembly": manifest.get("assembly"),
        "resolution": res,
        "hic_token": hic_token,
        "annotation_token": annotation_token,
        "mcoolPath": str(mcool) if mcool else None,
        "bigwig_tokens": bigwig_tokens,
        "bed_tokens": bed_tokens,
        # Two data layers, stated separately so the UI never implies that what the
        # browser can draw is the same as what was annotated:
        #   evidence — always genome-wide, ships with the app
        #   browser tracks — chr7-only subsets unless the user added full-genome files
        "scope": {
            "evidence_chroms": _evidence_chroms(ann),
            "map": {
                "available": hic_token is not None,
                "chroms": profile["chroms"] if hic_token else [],
                "assembly": profile["assembly"],
                # How the build was established: measured chromosome length, the file's
                # own declaration, or nothing. A consumer must not present an unverified
                # matrix as if its build were known.
                "assembly_verified": profile["assembly"] is not None,
                "assembly_source": profile["assembly_source"],
                "blocked_reason": map_blocked,
                "source": (
                    None if hic_token is None
                    else ("bundled_chr7" if "chr7" in Path(mcool).name.lower() else "full_genome")
                ),
            },
            "signal": {
                "available": bool(bigwig_tokens),
                "source": _signal_layer(tracks, bigwig_tokens),
                # Was hardcoded to [] for the full-genome layer, so a consumer reading
                # signal.chroms concluded "no chromosome has signal" while map.chroms
                # listed 24. Report what each track actually covers.
                "chroms": _signal_chroms(mgr, bigwig_tokens),
                # source/chroms above collapse both assays into one answer, which is only
                # right when they came from the same layer; per-track keeps a mixed
                # session (chr7 CTCF + full-genome RAD21) readable.
                "tracks": _signal_tracks(mgr, tracks, bigwig_tokens),
            },
        },
        "n_candidates": int(manifest.get("n_candidates") or 0),
        "bd1_n_methods": int(manifest.get("bd1_n_methods") or 0),
        "per_caller_boundary_counts": manifest.get("per_caller_boundary_counts", {}),
        "tier_counts": manifest.get("tier_counts", {}),
        "reference_note": manifest.get("reference_note"),
        "caveats": manifest.get("caveats", []),
        "tier_rule_version": "max_support_v2",
        "role": "downstream_evidence_annotation",
        "note": base_note,
    })


async def handle_tadvci_annotate(request: web.Request) -> web.Response:
    """Downstream path: grade user-supplied boundary BEDs (no built-in callers).

    Body:
      primary_bed (required): path to candidate boundaries
      method_beds (optional): list of {path, name} multi-caller panel for BD1
      ctcf_bw / rad21_bw / loops_path (optional): BD2–BD4 only when supplied —
        never auto-filled from demo data
      cell_line, assembly, resolution, hash_tag, output_dir
      hic_token (optional): if provided, bind the annotation to that existing
        Hi-C session token so the heatmap + tier share one key
    """
    try:
        body = await request.json()
    except Exception:
        return web.json_response({"error": "Invalid JSON body"}, status=400)
    if not isinstance(body, dict):
        return web.json_response({"error": "body must be a JSON object"}, status=400)

    primary = (body.get("primary_bed") or "").strip()
    if not primary:
        return web.json_response({"error": "primary_bed is required"}, status=400)
    if not os.path.isfile(primary):
        return web.json_response({"error": f"primary_bed not found: {primary}"}, status=404)

    raw_methods = body.get("method_beds") or []
    method_beds: list[tuple[str, str]] = []
    if isinstance(raw_methods, list):
        for item in raw_methods:
            if isinstance(item, str):
                method_beds.append((item, Path(item).stem))
            elif isinstance(item, dict) and item.get("path"):
                method_beds.append((str(item["path"]), str(item.get("name") or Path(item["path"]).stem)))
    # Also accept comma-separated method_bed_paths for simple clients
    for p in (body.get("method_bed_paths") or "").split(",") if isinstance(body.get("method_bed_paths"), str) else []:
        p = p.strip()
        if p:
            method_beds.append((p, Path(p).stem))

    try:
        resolution = _parse_resolution(body.get("resolution"), 25000)
    except ValueError as exc:
        return web.json_response({"error": str(exc)}, status=400)

    # Explicit user-supplied only — never fall back to demo ChIP / loops.
    ctcf_bw = (body.get("ctcf_bw") or "").strip() or None
    rad21_bw = (body.get("rad21_bw") or "").strip() or None
    loops_path = (body.get("loops_path") or "").strip() or None
    for label, path in (("ctcf_bw", ctcf_bw), ("rad21_bw", rad21_bw), ("loops_path", loops_path)):
        if path and not os.path.isfile(path):
            return web.json_response({"error": f"{label} not found: {path}"}, status=404)

    # hash_tag is concatenated into output FILE NAMES, so a separator in it escapes the
    # output directory. Unvalidated it produced a 500 instead of a 400, and the message
    # echoed the resolved server path straight back to the client.
    hash_tag = (body.get("hash_tag") or "annotate").strip() or "annotate"
    if not _SAFE_TAG_RE.fullmatch(hash_tag):
        return web.json_response(
            {"error": "hash_tag may contain only letters, digits, '.', '_' and '-'"},
            status=400)

    output_dir = str(body.get("output_dir") or "./tad_results")
    try:
        output_dir = _resolved_output_dir(output_dir)
    except ValueError as exc:
        return web.json_response({"error": str(exc)}, status=400)
    cell_line = body.get("cell_line")
    assembly = body.get("assembly")

    try:
        import asyncio
        from tad_vci.annotate_from_beds import annotate_from_beds
        result = await asyncio.get_event_loop().run_in_executor(
            None,
            lambda: annotate_from_beds(
                primary_bed=primary,
                method_beds=method_beds,
                resolution=resolution,
                ctcf_bw=ctcf_bw,
                rad21_bw=rad21_bw,
                loops_path=loops_path,
                cell_line=cell_line,
                assembly=assembly,
                hash_tag=hash_tag,
                output_dir=output_dir,
            ),
        )
    except FileNotFoundError as e:
        return web.json_response({"error": f"{e}"}, status=404)
    except Exception as e:
        # The raw exception text carried resolved SERVER paths back to the client. Same
        # rule main.py's /api/register already follows: detail to the server log, a
        # generic message to the caller.
        print(f"[tadvci/annotate] failed for primary_bed={primary!r} "
              f"output_dir={output_dir!r}: {type(e).__name__}: {e}")
        return web.json_response(
            {"error": "annotate failed while grading these boundaries — check the server "
                      "log for details."},
            status=500)

    ann_path = result.get("TADVCI_annotation")
    mgr = request.app["manager"]

    # Session token: reuse hic_token if client already registered mcool, else mint a
    # token for THIS run so /imported/{token} works without a Hi-C matrix.
    #
    # It used to be `mgr.register(primary, "bed")`, and register() de-duplicates by
    # path — so re-running the same candidate BED with a different caller panel handed
    # back the earlier run's token and silently re-pointed that session's tier table at
    # the new annotation. An already-open tab then showed another run's tiers with no
    # indication. Same failure the multicell demos were fixed for above; bind the token
    # to the run instead of to the input file. The BED TRACK tokens below still
    # de-duplicate by path, which is what we want for rendering.
    hic_token = body.get("hic_token")
    if hic_token and mgr.get(hic_token) is not None:
        session_token = hic_token
    else:
        run_key = hashlib.sha1(
            f"{os.path.abspath(primary)}|{output_dir}|{hash_tag}|{ann_path}".encode()
        ).hexdigest()[:12]
        session_token = f"annotate:{hash_tag}:{run_key}"

    set_imported_annotation(session_token, ann_path)

    bed_tokens = []
    # Primary first for track order
    bed_tokens.append({
        "token": mgr.register(primary, "bed"),
        "fileName": body.get("primary_name") or Path(primary).stem.replace("_boundaries", ""),
        "path": primary,
        "role": "primary",
    })
    for path, name in method_beds:
        if os.path.isfile(path):
            bed_tokens.append({
                "token": mgr.register(path, "bed"),
                "fileName": name or Path(path).stem.replace("_boundaries", ""),
                "path": path,
                "role": "method",
            })

    bigwig_tokens = []
    for label, path in (("CTCF", ctcf_bw), ("RAD21", rad21_bw)):
        if path and os.path.isfile(path):
            bigwig_tokens.append({
                "token": mgr.register(path, "bigwig"),
                "fileName": label,
                "path": path,
            })

    return web.json_response({
        "ok": True,
        "role": "downstream_evidence_annotation",
        "session_token": session_token,
        "annotation_token": session_token,
        "hic_token": hic_token if (hic_token and mgr.get(hic_token) is not None) else None,
        "TADVCI_annotation": ann_path,
        "TADVCI_tier": result.get("TADVCI_tier"),
        "n_boundaries": result.get("n_boundaries"),
        "n_methods": result.get("n_methods"),
        "bd1_panel": result.get("bd1_panel"),
        "tier_dist": result.get("tier_dist"),
        "bed_tokens": bed_tokens,
        "bigwig_tokens": bigwig_tokens,
        "cell_line": cell_line,
        "assembly": assembly,
        "resolution": resolution,
        "bd2_assessed": bool(ctcf_bw),
        "bd3_assessed": bool(rad21_bw),
        "bd4_assessed": bool(loops_path),
        "note": (
            "Boundaries graded as inputs. Built-in callers were not run. "
            "BD2–BD4 are assessable only for assays you supplied."
        ),
    })


async def handle_tadvci_imported(request: web.Request) -> web.Response:
    """Per-import TAD-VCI tier (computed on the user's uploaded mcool), keyed by token."""
    token = request.match_info["token"]
    chrom = request.match_info["chrom"]
    # Unknown token -> 404 (consistent with /triage), so the frontend can tell an expired/
    # invalid import apart from a valid import that simply has no boundaries on this chrom.
    annot_path = IMPORTED_ANNOT.get(token)
    if not annot_path or not Path(annot_path).exists():
        return web.json_response(
            {"error": f"no imported annotation for token {token}", "token": token,
             "chrom": str(chrom).replace("chr", ""), "n": 0, "boundaries": []},
            status=404)
    # Default is the slim track payload. `view=full` keeps the historical whole-card dump
    # for anything that still wants it; it is 60x larger and is not what the viewer loads.
    if request.query.get("view") == "full":
        cards = _cards_for_token(token, chrom)
        public = [{k: v for k, v in c.items() if k != "_row"} for c in cards]
    else:
        public = _track_for_token(token, chrom)
    return web.json_response({"token": token, "chrom": str(chrom).replace("chr", ""),
                              "n": len(public), "view": request.query.get("view", "track"),
                              "boundaries": public})


async def handle_tadvci_card(request: web.Request) -> web.Response:
    """One full evidence card. The tier track serves ticks; this serves the panel."""
    token = request.match_info["token"]
    chrom = request.match_info["chrom"]
    try:
        pos = int(request.match_info["pos"])
    except (TypeError, ValueError):
        return web.json_response({"error": "pos must be an integer"}, status=400)
    annot_path = IMPORTED_ANNOT.get(token)
    if not annot_path or not Path(annot_path).exists():
        return web.json_response(
            {"error": f"no imported annotation for token {token}", "token": token}, status=404)
    card = _card_for_position(token, chrom, pos)
    if card is None:
        return web.json_response(
            {"error": f"no boundary at {str(chrom).replace('chr', '')}:{pos} for this token"},
            status=404)
    return web.json_response({"token": token, "boundary": card})


async def handle_tadvci_triage(request: web.Request) -> web.Response:
    """Report that the legacy BCP/conformal triage layer is scientifically disabled."""
    token = request.match_info["token"]
    path = IMPORTED_ANNOT.get(token)
    if not path or not Path(path).exists():
        return web.json_response({"error": "no annotation for token"}, status=404)
    df = _imported_df_cache.get(token)
    if df is None:
        df = pd.read_csv(path, sep="\t")
        _imported_df_cache[token] = df
    chrom_q = request.query.get("chrom")
    scope = "whole sample"
    if chrom_q:
        ck = str(chrom_q).replace("chr", "")
        df = df[df["chrom"].astype(str).str.replace("chr", "") == ck]
        scope = f"chr{ck}"
    n = int(len(df))
    return web.json_response({
        "n": n,
        "scope": scope,
        "scored": 0,
        "not_scored": n,
        "bcp_available": False,
        "auto_confident": 0,
        "routed_to_expert": 0,
        "routed_frac": None,
        "coverage_guarantee": None,
        "concordance_qc": {},
        "major_discordant": 0,
        "bcp_status": "retired_pending_provenance_safe_labels_and_fixed_evaluation",
        "note": "The legacy BCP and conformal routing claims failed the July 2026 data/model "
                "audit and are disabled. Evidence cards and deterministic tiers remain available.",
    })


async def handle_tadvci_signout(request: web.Request) -> web.Response:
    try:
        body = await request.json()
    except Exception:
        return web.json_response({"error": "Invalid JSON body"}, status=400)
    if not isinstance(body, dict):
        return web.json_response({"error": "body must be a JSON object"}, status=400)
    bid = body.get("boundary_id", "")
    rationale = (body.get("rationale") or "").strip()
    final_tier = body.get("final_tier")
    if not rationale:
        return web.json_response({"error": "rationale required"}, status=422)
    if len(rationale) > 10000:
        return web.json_response({"error": "rationale too long (max 10000 chars)"}, status=422)
    if not final_tier:
        return web.json_response({"error": "final_tier required"}, status=422)
    if final_tier not in {"D1", "D2", "D3", "D4", "D5"}:
        return web.json_response({"error": "final_tier must be D1..D5"}, status=422)
    verdict = body.get("verdict", "")
    allowed_verdicts = {
        "supported_for_review",
        "not_supported_for_review",
        "uncertain",
    }
    if verdict not in allowed_verdicts:
        return web.json_response(
            {"error": f"verdict must be one of {sorted(allowed_verdicts)}"}, status=422
        )
    actor_label = str(body.get("actor_label") or "").strip()
    if not actor_label:
        return web.json_response(
            {"error": "actor_label required (label is not authenticated identity)"},
            status=422,
        )
    if len(actor_label) > 200:
        return web.json_response({"error": "actor_label too long (max 200 chars)"}, status=422)
    chrom = bid.split(":")[0] if ":" in bid else "7"
    # Resolve the evidence card from ONE explicitly named source.
    #
    # A boundary_id is only "chrom:pos", and the demo packages share thousands of
    # positions (GM12878_25kb and IMR90_25kb share 825 on chr7 alone, 414 of them at a
    # different tier). The previous order tried the built-in cell FIRST — with an unknown
    # cell_line silently downgraded to GM12878 — and, when no token was given, scanned
    # every registered import and took the first hit. That produced a signed record
    # labelled IMR90 that carried GM12878's evidence and still passed its own checksum.
    tok = body.get("token")
    requested_cell = body.get("cell_line")
    evidence_source = None
    if tok:
        if tok not in IMPORTED_ANNOT:
            return web.json_response(
                {"error": f"unknown token {tok}; load the dataset again before recording"},
                status=404)
        card = next((c for c in _cards_for_token(tok, chrom) if c["boundary_id"] == bid), None)
        evidence_source = f"token:{tok}"
    elif requested_cell in CELL_REGISTRY:
        card = next((c for c in _cards_for(chrom, requested_cell) if c["boundary_id"] == bid), None)
        evidence_source = f"builtin:{requested_cell}"
    else:
        # No silent search. Which dataset a verdict is about is part of the record.
        return web.json_response(
            {"error": "token required (or a built-in cell_line): the dataset a decision "
                      "refers to cannot be inferred from boundary_id alone"},
            status=422)
    if card is None:
        return web.json_response(
            {"error": f"unknown boundary {bid} in {evidence_source}"}, status=404)
    decision = {"final_tier": final_tier, "verdict": verdict, "rationale": rationale}
    v = next_version(bid, _LOG)
    # Context metadata is covered by the v3 per-record checksum.
    try:
        resolution = int(body["resolution"]) if body.get("resolution") not in (None, "") else None
    except (TypeError, ValueError):
        resolution = None
    cell_line = body.get("cell_line") or None
    # Resolve the assembly that goes into the versioned decision snapshot. A live (non-registry)
    # upload must NOT be stamped hg19 just because the served cell defaults to GM12878 — an hg38
    # upload would otherwise enter the sha256 provenance as hg19. Resolution order:
    #   1. the assembly column written into the imported annotation TSV (authoritative; the build
    #      that user declared at import time),  2. an explicit assembly in the sign-out payload
    #      (frontend fallback),  3. the registry assembly of the built-in cell THIS record
    #      was actually read from (never a cell the caller merely mentioned), 4. _ASSEMBLY.
    builtin_cell = requested_cell if evidence_source.startswith("builtin:") else None
    assembly = (_assembly_for_boundary(tok, chrom, bid)
                or _norm_assembly(body.get("assembly"))
                or (CELL_REGISTRY[builtin_cell]["assembly"] if builtin_cell in CELL_REGISTRY else None)
                or _ASSEMBLY)
    snap = make_snapshot(card, assembly, actor_label,
                         datetime.now(timezone.utc).isoformat(timespec="seconds"),
                         decision=decision, prior_version=v,
                         resolution=resolution, cell_line=cell_line)
    append_snapshot(snap, _LOG)
    rec = {**snap.__dict__}
    return web.json_response({
        "ok": True,
        "snapshot": rec,
        "checksum_matches": verify_snapshot(rec),
        "verification_scope": (
            "per-record canonical-field checksum only; no identity authentication, "
            "log-completeness proof, hash chain, or immutable storage"
        ),
    })


async def handle_tadvci_history(request: web.Request) -> web.Response:
    bid = request.match_info["boundary_id"]
    return web.json_response({
        "boundary_id": bid,
        "history": [_public_record(record) for record in load_history(bid, _LOG)],
    })


def _public_record(record: dict) -> dict:
    """Expose neutral v3 names even when reading a legacy v1/v2 record."""
    public = dict(record)
    if "review_decision" not in public and "expert_override" in public:
        public["review_decision"] = public.pop("expert_override")
        public["legacy_field_names_normalized"] = True
    if "actor_label" not in public and "curator" in public:
        public["actor_label"] = public.pop("curator")
        public["legacy_field_names_normalized"] = True
    return public


async def handle_tadvci_records(request: web.Request) -> web.Response:
    """All application-level versioned decision records for the My Labels page."""
    import json
    recs = []
    if _LOG.exists():
        for line in _LOG.read_text(encoding="utf-8").splitlines():
            if line.strip():
                recs.append(_public_record(json.loads(line)))
    recs.sort(key=lambda r: (r.get("boundary_id", ""), r.get("version", 0)))
    return web.json_response({"n": len(recs), "records": recs})


def register_tadvci_routes(app: web.Application) -> None:
    app.add_routes([
        web.get("/api/tadvci/annotation/{chrom}", handle_tadvci_annotation),
        web.get("/api/tadvci/cells", handle_tadvci_cells),
        web.post("/api/tadvci/use_cell/{cell}", handle_tadvci_use_cell),
        web.get("/api/tadvci/multicell_datasets", handle_tadvci_multicell_datasets),
        web.post("/api/tadvci/use_multicell_example", handle_tadvci_use_multicell_example),
        web.post("/api/tadvci/annotate", handle_tadvci_annotate),
        web.get("/api/tadvci/imported/{token}/{chrom}", handle_tadvci_imported),
        web.get("/api/tadvci/card/{token}/{chrom}/{pos}", handle_tadvci_card),
        web.get("/api/tadvci/triage/{token}", handle_tadvci_triage),
        web.post("/api/tadvci/signout", handle_tadvci_signout),
        web.get("/api/tadvci/history/{boundary_id}", handle_tadvci_history),
        web.get("/api/tadvci/records", handle_tadvci_records),
    ])
