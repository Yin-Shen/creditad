"""Download and checksum-verify the active raw GM12878 viewer inputs only."""
from __future__ import annotations

import asyncio
import hashlib
import os
import re
from pathlib import Path

try:
    import aiohttp
    from aiohttp import web
except ModuleNotFoundError as _exc:  # pragma: no cover - core-only install
    # Download/serve helpers for the GUI. aiohttp lives in the `app` extra; the CLI and
    # the library annotation path never import this module (audit 2026-08-12).
    raise ModuleNotFoundError(
        "data_manager needs the optional HTTP/GUI dependencies. Install them with: "
        "pip install 'creditad[app]'   (the `creditad annotate` CLI does not require them)."
    ) from _exc


_BACKEND_DIR = Path(__file__).resolve().parent
_EXTRA_MODE_DIR = _BACKEND_DIR / "extra_mode"
_CHR7_TRACKS_DIR = _EXTRA_MODE_DIR / "chr7_tracks"
# Resolves to creditad/example_data/chr7_tracks in the source tree and to
# resources/example_data/chr7_tracks in the packaged app — the same place
# tad_vci.api_routes looks for the demo matrices.
_EXAMPLE_CHR7_DIR = _BACKEND_DIR.parent / "example_data" / "chr7_tracks"

# ── Two layers (do not mix names) ──────────────────────────────────────────
# 1) SHIPPED DEMO TRACKS (chr7 only) — small, always in Portable / install tree.
#    Filenames MUST contain "_chr7" so they are never confused with full-genome
#    ENCODE downloads.
# 2) OPTIONAL FULL-GENOME — multi-GB; not listed as "example". Users place them
#    under CREDITAD_SHARED_RAW via Dataset setup guide (DEMO_DATA_CATALOG).
#
# Genome-wide *evidence tables* for all 3 cells × 2 resolutions live in
# example_data/multicell/ (annotation.tsv). Browser heatmap/ChIP for demos use
# chr7 subsets unless a full-genome file is available under CREDITAD_SHARED_RAW.

# Every shipped chr7 asset, for all three demo cell lines. The inventory used to
# list GM12878 only (plus the retired shared `chr7.mcool`), so Data Manager showed
# one cell line while the app actually shipped tracks for three, and — after the
# GRCh37 matrix was withdrawn — reported a missing file that no longer exists.
#
# The contact matrices live under example_data/chr7_tracks (one per cell line,
# GRCh38, built by scripts/build_chr7_demo_matrices.py); the ChIP subsets live
# under backend/extra_mode/chr7_tracks. Both directories resolve identically in
# the source tree and in the packaged resources/ layout.
_DEMO_CELLS = (
    {
        "cell": "GM12878",
        "ctcf": "ENCFF734CUT",
        "rad21": "ENCFF571ZJJ",
        "mcool_md5": "9095137002f997939bebedff049e25aa",
        "mcool_mb": 91,
        "ctcf_md5": "5758946f2eb3f031997d0f13dda427b7",
        "ctcf_mb": 5,
        "rad21_md5": "6d34b56e14ea3e30201447cbc1a71021",
        "rad21_mb": 5,
    },
    {
        "cell": "IMR90",
        "ctcf": "ENCFF105FHL",
        "rad21": "ENCFF048PZI",
        "mcool_md5": "05ec7a43e201d9dea62ce054e4e80b13",
        "mcool_mb": 42,
        "ctcf_md5": "b8034daf916b34831816501800f55d8b",
        "ctcf_mb": 5,
        "rad21_md5": "2d379d8e6a1754a09bd9bc096a6f212c",
        "rad21_mb": 6,
    },
    {
        "cell": "HepG2",
        "ctcf": "ENCFF357NFO",
        "rad21": "ENCFF972ODZ",
        "mcool_md5": "2a51cdf6a3327d785bec8757f1dba4cc",
        "mcool_mb": 48,
        "ctcf_md5": "0d4ecb3d7775946df477b694fbfa17e2",
        "ctcf_mb": 4,
        "rad21_md5": "f00592d4d02babc4898f0264d4c42d78",
        "rad21_mb": 5,
    },
)


def _build_data_files() -> list[dict]:
    out: list[dict] = []
    for spec in _DEMO_CELLS:
        cell = spec["cell"]
        out.append(
            {
                "id": f"{cell.lower()}_chr7_mcool",
                "category": "shipped_chr7_demo",
                "cell_line": cell,
                "filename": f"{cell}_chr7_hg38.mcool",
                "dir": str(_EXAMPLE_CHR7_DIR),
                "url": "",
                "md5": spec["mcool_md5"],
                "size_mb": spec["mcool_mb"],
                "required": True,
                "description": (
                    f"SHIPPED demo Hi-C for {cell}: chr7 multi-resolution cooler "
                    "(GRCh38, 10 kb–1 Mb), sliced from that cell line's own matrix. "
                    "Not a full-genome matrix."
                ),
            }
        )
        out.append(
            {
                "id": f"{cell.lower()}_ctcf_chr7_bigwig",
                "category": "shipped_chr7_demo",
                "cell_line": cell,
                "filename": f"{cell}_CTCF_{spec['ctcf']}_hg38_chr7.bigWig",
                "dir": str(_CHR7_TRACKS_DIR),
                "url": "",
                "md5": spec["ctcf_md5"],
                "size_mb": spec["ctcf_mb"],
                "required": True,
                "description": (
                    f"SHIPPED demo CTCF ({cell} {spec['ctcf']} fc/control, GRCh38) — "
                    "chr7 subset only. Full-genome track: Dataset setup guide → ENCODE."
                ),
            }
        )
        out.append(
            {
                "id": f"{cell.lower()}_rad21_chr7_bigwig",
                "category": "shipped_chr7_demo",
                "cell_line": cell,
                "filename": f"{cell}_RAD21_{spec['rad21']}_hg38_chr7.bigWig",
                "dir": str(_CHR7_TRACKS_DIR),
                "url": "",
                "md5": spec["rad21_md5"],
                "size_mb": spec["rad21_mb"],
                "required": True,
                "description": (
                    f"SHIPPED demo RAD21 ({cell} {spec['rad21']} fc/control, GRCh38) — "
                    "chr7 subset only. Full-genome track: Dataset setup guide → ENCODE."
                ),
            }
        )
    return out


DATA_FILES = _build_data_files()

VISIBLE_CATEGORIES = frozenset({"shipped_chr7_demo"})
VISIBLE_FILE_IDS = frozenset(entry["id"] for entry in DATA_FILES)
_download_progress: dict[str, dict] = {}


def _visible_entries() -> tuple[dict, ...]:
    """Return the exact public download inventory, failing closed on bad config."""
    entries = tuple(
        entry
        for entry in DATA_FILES
        if entry.get("id") in VISIBLE_FILE_IDS
        and entry.get("category") in VISIBLE_CATEGORIES
    )
    visible_ids = {entry["id"] for entry in entries}
    if visible_ids != VISIBLE_FILE_IDS or len(entries) != len(VISIBLE_FILE_IDS):
        raise RuntimeError(
            "Data Manager visible inventory is incomplete, duplicated, or mis-categorized"
        )
    return entries


def _visible_entry(file_id: str | None) -> dict | None:
    return next((entry for entry in _visible_entries() if entry["id"] == file_id), None)


def _install_path(entry: dict) -> Path:
    """Canonical install location under the Data Manager inventory (extra_mode/…)."""
    return Path(entry["dir"]) / entry["filename"]


def _file_path(entry: dict) -> Path:
    """Where this entry's bytes live for read/verify.

    Prefer the inventory install path. If missing, reuse an existing shared-raw copy
    of the same ENCODE accession so the user is not forced to re-download ~880 MB.
    Override the shared root with CREDITAD_SHARED_RAW.
    """
    primary = _install_path(entry)
    if primary.is_file() and primary.stat().st_size > 0:
        return primary
    for existing in _shared_copies(entry):
        return existing
    return primary


def dev_shared_raw_default() -> Path | None:
    """The development host's optional full-genome tree, or None.

    Gated on POSIX + the directory actually existing. A packaged Windows build must
    resolve tracks from its own install tree or from CREDITAD_SHARED_RAW; an absolute
    developer path left in the search order resolves through wine's Z: drive during
    the packaging smoke test, so the branch a real user takes was never exercised.
    Override with CREDITAD_DEV_SHARED_RAW.

    The fallback is DERIVED from where this file sits (…/<tree>/creditad/backend/ ->
    <tree>/../_shared_data/raw), never written as a literal. A literal developer home in
    this module is shipped inside the Windows package, where the packaging verifier
    rejects it (`check_no_host_paths`) — correctly, because the string names one
    machine's directory layout to every user who unpacks the installer.
    """
    if os.name != "posix":
        return None
    env = os.environ.get("CREDITAD_DEV_SHARED_RAW")
    if env:
        root = Path(env)
    else:
        here = Path(__file__).resolve()
        if len(here.parents) < 4:
            return None
        root = here.parents[3] / "_shared_data" / "raw"
    return root if root.is_dir() else None


def _shared_copies(entry: dict):
    """Yield existing files under the shared raw tree matching this entry's accession."""
    acc = entry.get("accession") or _accession_of(entry)
    if not acc:
        return
    root = shared_raw_root()
    if not root.is_dir():
        return
    # ENCFF….bigWig may sit two or three levels under the raw root.
    for pattern in (f"*/*/{acc}.bigWig", f"*/*/*/{acc}.bigWig", f"**/{acc}.bigWig"):
        for hit in root.glob(pattern):
            if hit.is_file() and hit.stat().st_size > 0:
                yield hit


def _accession_of(entry: dict) -> str | None:
    m = re.search(r"(ENCFF[A-Z0-9]+)", entry.get("url", "") or entry.get("filename", ""))
    return m.group(1) if m else None


def _file_exists(entry: dict) -> bool:
    path = _file_path(entry)
    return path.is_file() and path.stat().st_size > 0


def _ensure_install_link(entry: dict) -> Path | None:
    """If a shared copy exists, materialise it at the install path (symlink or copy).

    Returns the install path when the file is usable there, else None.
    """
    install = _install_path(entry)
    if install.is_file() and install.stat().st_size > 0:
        return install
    sources = list(_shared_copies(entry))
    if not sources:
        return None
    src = sources[0]
    install.parent.mkdir(parents=True, exist_ok=True)
    try:
        if install.exists() or install.is_symlink():
            install.unlink()
        install.symlink_to(src.resolve())
    except OSError:
        # Cross-device or no symlink permission: hard copy (slow but correct).
        import shutil
        shutil.copy2(src, install)
    if install.is_file() and install.stat().st_size > 0:
        return install
    return None


def _md5_file(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324 - fixed upstream integrity identifier
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


async def handle_data_status(request: web.Request) -> web.Response:
    files = []
    for entry in _visible_entries():
        # Opportunistically link shared copies so Refresh shows Ready without a re-download.
        if not _install_path(entry).is_file():
            _ensure_install_link(entry)
        resolved = _file_path(entry)
        present = resolved.is_file() and resolved.stat().st_size > 0
        info = {
            "id": entry["id"],
            "category": entry["category"],
            # The shipped inventory covers three cell lines; without this the UI
            # cannot group them and the list reads as one undifferentiated pile.
            "cell_line": entry.get("cell_line"),
            "filename": entry["filename"],
            "description": entry["description"],
            "size_mb": entry["size_mb"],
            "required": entry["required"],
            "has_url": bool(entry["url"]),
            "present": present,
            "dir": entry["dir"],
            "resolved_path": str(resolved) if present else None,
            "url": entry["url"] or None,
        }
        if entry["id"] in _download_progress:
            info["download"] = _download_progress[entry["id"]]
        files.append(info)
    return web.json_response({"files": files})


async def handle_data_download(request: web.Request) -> web.Response:
    try:
        body = await request.json()
    except Exception:
        return web.json_response({"error": "Invalid JSON body"}, status=400)
    if not isinstance(body, dict):
        return web.json_response({"error": "Invalid JSON body"}, status=400)
    # Accept both "id" (current UI) and legacy "file_id".
    file_id = body.get("id") or body.get("file_id")
    entry = _visible_entry(file_id)
    if entry is None:
        return web.json_response({"error": f"Unknown file id: {file_id}"}, status=404)
    if _download_progress.get(file_id, {}).get("status") == "downloading":
        return web.json_response({"error": "Download already in progress"}, status=409)

    # Already on disk (install path or linkable shared copy) — finish without network.
    linked = _ensure_install_link(entry)
    if linked is not None:
        size_mb = round(linked.stat().st_size / (1024 * 1024), 1)
        _download_progress[file_id] = {
            "status": "done",
            "percent": 100,
            "downloaded_mb": size_mb,
            "total_mb": size_mb,
            "error": None,
            "source": "local_shared_copy",
        }
        return web.json_response(
            {"status": "done", "id": file_id, "source": "local_shared_copy", "path": str(linked)}
        )

    url = entry["url"]
    if not url:
        return web.json_response(
            {
                "error": (
                    f"No configured download URL for {entry['filename']}; "
                    f"place it manually in {entry['dir']}"
                )
            },
            status=400,
        )
    asyncio.create_task(_download_file(entry))
    return web.json_response({"status": "started", "id": file_id})


async def handle_data_verify(request: web.Request) -> web.Response:
    try:
        body = await request.json()
    except Exception:
        return web.json_response({"error": "Invalid JSON body"}, status=400)
    file_id = body.get("id") if isinstance(body, dict) else None
    entry = _visible_entry(file_id)
    if entry is None:
        return web.json_response({"error": f"Unknown file id: {file_id}"}, status=404)
    path = _file_path(entry)
    if not path.is_file():
        return web.json_response({"verified": False, "reason": "File not found"})
    actual = _md5_file(path)
    expected = entry["md5"]
    verified = actual == expected
    return web.json_response(
        {
            "verified": verified,
            "expected": expected,
            "actual": actual,
            "reason": "MD5 match" if verified else "MD5 mismatch - file may be corrupted",
        }
    )


async def handle_data_download_progress(request: web.Request) -> web.Response:
    file_id = request.match_info["id"]
    if _visible_entry(file_id) is None:
        return web.json_response({"error": f"Unknown file id: {file_id}"}, status=404)
    return web.json_response(_download_progress.get(file_id, {"status": "idle"}))


async def _download_file(entry: dict) -> None:
    """Download into the inventory install path (never overwrite a shared-raw path)."""
    file_id = entry["id"]
    destination = _install_path(entry)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    _download_progress[file_id] = {
        "status": "downloading",
        "percent": 0,
        "downloaded_mb": 0,
        "total_mb": entry["size_mb"],
        "error": None,
    }
    try:
        # Prefer local shared copy one more time (race with concurrent Refresh).
        linked = _ensure_install_link(entry)
        if linked is not None:
            size_mb = round(linked.stat().st_size / (1024 * 1024), 1)
            _download_progress[file_id] = {
                "status": "done",
                "percent": 100,
                "downloaded_mb": size_mb,
                "total_mb": size_mb,
                "error": None,
                "source": "local_shared_copy",
            }
            return

        timeout = aiohttp.ClientTimeout(total=7200, connect=60)
        headers = {
            # ENCODE occasionally rejects empty / bot-like User-Agents.
            "User-Agent": "CrediTAD-DataManager/0.5 (scientific reproducibility download)",
            "Accept": "*/*",
        }
        async with aiohttp.ClientSession(timeout=timeout, headers=headers) as session:
            async with session.get(entry["url"], allow_redirects=True) as response:
                if response.status != 200:
                    raise RuntimeError(
                        f"HTTP {response.status} fetching {entry['url']}: "
                        f"{response.reason or 'download failed'}. "
                        f"If the network blocks ENCODE/S3, place the file manually in "
                        f"{entry['dir']} as {entry['filename']}."
                    )
                total = response.content_length or entry["size_mb"] * 1024 * 1024
                downloaded = 0
                with temporary.open("wb") as handle:
                    async for chunk in response.content.iter_chunked(256 * 1024):
                        handle.write(chunk)
                        downloaded += len(chunk)
                        _download_progress[file_id] = {
                            "status": "downloading",
                            "percent": min(99, int(downloaded * 100 / max(total, 1))),
                            "downloaded_mb": round(downloaded / (1024 * 1024), 1),
                            "total_mb": round(total / (1024 * 1024), 1),
                            "error": None,
                        }
                    handle.flush()
                    os.fsync(handle.fileno())
        actual = _md5_file(temporary)
        if actual != entry["md5"]:
            raise RuntimeError(
                f"MD5 mismatch: expected {entry['md5']}, got {actual}"
            )
        os.replace(temporary, destination)
        size_mb = round(destination.stat().st_size / (1024 * 1024), 1)
        _download_progress[file_id] = {
            "status": "done",
            "percent": 100,
            "downloaded_mb": size_mb,
            "total_mb": size_mb,
            "error": None,
            "source": "network",
        }
    except asyncio.CancelledError:
        _cleanup_tmp(temporary)
        _download_progress[file_id] = {
            "status": "cancelled",
            "error": "Download cancelled",
        }
        raise
    except Exception as exc:
        _cleanup_tmp(temporary)
        msg = str(exc).strip() or repr(exc)
        _download_progress[file_id] = {"status": "error", "error": msg}


def _cleanup_tmp(path: Path) -> None:
    try:
        path.unlink(missing_ok=True)
    except OSError:
        pass


# Public demo-dataset catalog for Windows / offline packaging. These are the
# files the six multicell demos need for heatmaps + BD2/BD3/BD4. Multi-GB mcool
# files are listed for guided manual download (not auto-fetched by Data Manager).
DEMO_DATA_CATALOG = [
    {
        "cell_line": "GM12878",
        "assembly": "hg38",
        "items": [
            {
                "role": "mcool",
                "label": "Hi-C multi-resolution matrix",
                "accession": "4DNFIXP4QG5B",
                "filename": "4DNFIXP4QG5B.mcool",
                "url": "https://data.4dnucleome.org/files-processed/4DNFIXP4QG5B/",
                "size_note": "~20–30 GB (multi-resolution)",
                "place_under": "gm12878/hic/",
                "used_for": "Demo heatmap (Visualization)",
            },
            {
                "role": "ctcf",
                "label": "CTCF ChIP-seq (fold change over control)",
                "accession": "ENCFF734CUT",
                "filename": "ENCFF734CUT.bigWig",
                "url": "https://www.encodeproject.org/files/ENCFF734CUT/@@download/ENCFF734CUT.bigWig",
                "size_note": "~356 MB",
                "place_under": "gm12878/ctcf/",
                "used_for": "Optional full-genome BD2 browser track (not the shipped *_chr7* demo file)",
            },
            {
                "role": "rad21",
                "label": "RAD21 ChIP-seq (fold change over control)",
                "accession": "ENCFF571ZJJ",
                "filename": "ENCFF571ZJJ.bigWig",
                "url": "https://www.encodeproject.org/files/ENCFF571ZJJ/@@download/ENCFF571ZJJ.bigWig",
                "size_note": "~486 MB",
                "place_under": "gm12878/rad21/",
                "used_for": "Optional full-genome BD3 browser track (not the shipped *_chr7* demo file)",
            },
            {
                "role": "loops",
                "label": "HiCCUPS loop calls (BEDPE)",
                "accession": "GM12878_HiCCUPS_loops_hg38",
                "filename": "GM12878_HiCCUPS_loops_hg38.bedpe",
                "url": "",
                "size_note": "~0.4 MB",
                "place_under": "data/loops/",
                "used_for": "BD4",
                "ship_with_package": True,
                "why_no_url": "No single stable public file URL; already in the repo (extra_mode/loops/). Ship with the Windows build.",
            },
        ],
    },
    {
        "cell_line": "IMR90",
        "assembly": "hg38",
        "items": [
            {
                "role": "mcool",
                "label": "Hi-C multi-resolution matrix",
                "accession": "4DNFIJTOIGOI",
                "filename": "4DNFIJTOIGOI.mcool",
                "url": "https://data.4dnucleome.org/files-processed/4DNFIJTOIGOI/",
                "size_note": "~20–30 GB",
                "place_under": "imr90/hic/",
                "used_for": "Demo heatmap",
            },
            {
                "role": "ctcf",
                "label": "CTCF ChIP-seq (fold change over control)",
                "accession": "ENCFF105FHL",
                "filename": "ENCFF105FHL.bigWig",
                "url": "https://www.encodeproject.org/files/ENCFF105FHL/@@download/ENCFF105FHL.bigWig",
                "size_note": "~300–400 MB",
                "place_under": "imr90/ctcf/",
                "used_for": "BD2",
            },
            {
                "role": "rad21",
                "label": "RAD21 ChIP-seq (fold change over control)",
                "accession": "ENCFF048PZI",
                "filename": "ENCFF048PZI.bigWig",
                "url": "https://www.encodeproject.org/files/ENCFF048PZI/@@download/ENCFF048PZI.bigWig",
                "size_note": "~300–400 MB",
                "place_under": "imr90/rad21/",
                "used_for": "BD3",
            },
            {
                "role": "loops",
                "label": "HiCCUPS loop calls (BEDPE)",
                "accession": "IMR90_HiCCUPS_loops_hg38",
                "filename": "IMR90_HiCCUPS_loops_hg38.bedpe",
                "url": "",
                "size_note": "~0.4 MB",
                "place_under": "data/loops/",
                "used_for": "BD4",
                "ship_with_package": True,
                "why_no_url": "No single stable public file URL; already in the repo (extra_mode/loops/). Ship with the Windows build.",
            },
        ],
    },
    {
        "cell_line": "HepG2",
        "assembly": "hg38",
        "items": [
            {
                "role": "mcool",
                "label": "Hi-C multi-resolution matrix",
                "accession": "4DNFIS6HAUPP",
                "filename": "4DNFIS6HAUPP.mcool",
                "url": "https://data.4dnucleome.org/files-processed/4DNFIS6HAUPP/",
                "size_note": "~20–30 GB",
                "place_under": "hepg2/hic/",
                "used_for": "Demo heatmap",
            },
            {
                "role": "ctcf",
                "label": "CTCF ChIP-seq (fold change over control)",
                "accession": "ENCFF357NFO",
                "filename": "ENCFF357NFO.bigWig",
                "url": "https://www.encodeproject.org/files/ENCFF357NFO/@@download/ENCFF357NFO.bigWig",
                "size_note": "~300–400 MB",
                "place_under": "hepg2/ctcf/",
                "used_for": "BD2",
            },
            {
                "role": "rad21",
                "label": "RAD21 ChIP-seq (fold change over control)",
                "accession": "ENCFF972ODZ",
                "filename": "ENCFF972ODZ.bigWig",
                "url": "https://www.encodeproject.org/files/ENCFF972ODZ/@@download/ENCFF972ODZ.bigWig",
                "size_note": "~300–400 MB",
                "place_under": "hepg2/rad21/",
                "used_for": "BD3",
            },
            {
                "role": "loops",
                "label": "HiCCUPS loop calls (BEDPE)",
                "accession": "HepG2_HiCCUPS_loops_hg38",
                "filename": "HepG2_HiCCUPS_loops_hg38.bedpe",
                "url": "",
                "size_note": "~3.9 MB",
                "place_under": "data/loops/",
                "used_for": "BD4",
                "ship_with_package": True,
                "why_no_url": "No single stable public file URL; already in the repo (extra_mode/loops/). Ship with the Windows build.",
            },
        ],
    },
]


def shared_raw_root() -> Path:
    """The directory the optional full-genome files are read from.

    Single source of truth for the UI. It must stay the FIRST root that
    tad_vci.api_routes._multicell_track_paths searches, otherwise the guide would
    tell users to fill a folder the resolver never looks in.
    """
    env = (os.environ.get("CREDITAD_SHARED_RAW") or "").strip()
    if env:
        return Path(env)
    dev = dev_shared_raw_default()
    if dev is not None:
        return dev
    # Packaged builds: the guide must name a folder the resolver actually searches
    # next (tad_vci.api_routes._shared_raw_roots), not a developer's host path.
    return _BACKEND_DIR.parent / "data" / "raw"


def _catalog_with_status() -> list[dict]:
    """The optional-download catalog, annotated with what is actually on disk.

    Without this the guide is a static shopping list: a user who has downloaded
    everything gets the same screen as one who has downloaded nothing, and neither
    can tell whether the application can see the files.
    """
    root = shared_raw_root()
    cells = []
    for cell in DEMO_DATA_CATALOG:
        items = []
        for item in cell["items"]:
            entry = dict(item)
            if item.get("ship_with_package"):
                # Ships in the install tree; nothing for the user to place.
                entry["expected_path"] = None
                entry["present"] = True
            else:
                expected = root / item["place_under"] / item["filename"]
                entry["expected_path"] = str(expected)
                entry["present"] = expected.is_file() and expected.stat().st_size > 0
            items.append(entry)
        cells.append({**cell, "items": items})
    return cells


async def handle_data_catalog(request: web.Request) -> web.Response:
    """Guided inventory for Windows / offline installs (not auto-download)."""
    root = shared_raw_root()
    cells = _catalog_with_status()
    optional = [i for c in cells for i in c["items"] if not i.get("ship_with_package")]
    present = sum(1 for i in optional if i["present"])
    return web.json_response(
        {
            "ok": True,
            "purpose": (
                "Optional full-genome files. The app runs without them: every demo "
                "package ships genome-wide evidence plus a chr7 contact map and chr7 "
                "ChIP tracks. Add a file here and the browser uses it instead of the "
                "chr7 subset, for that cell line and that track only."
            ),
            "shared_raw_root": str(root),
            "shared_raw_exists": root.is_dir(),
            "optional_present": present,
            "optional_total": len(optional),
            "notes": [
                "All demo packages are GRCh38 / hg38.",
                "SHIPPED (small, always present): multicell/* evidence tables (genome-wide) "
                "+ per-cell <CELL>_chr7_hg38.mcool + chr7_tracks/*_chr7.bigWig, for all three cells.",
                "OPTIONAL (large, you place them): the files below, under the folder shown above. "
                "Keep the sub-folder layout exactly — that is where the app looks.",
                "Resolution is per track and per cell line: a full-genome CTCF file is used even "
                "if the matrix for that cell line is still the chr7 subset.",
                "CTCF and RAD21 must be fold-change-over-control (not signal p-value).",
            ],
            "cells": cells,
        }
    )


def register_data_routes(app: web.Application) -> None:
    app.add_routes(
        [
            web.get("/api/data/status", handle_data_status),
            web.get("/api/data/catalog", handle_data_catalog),
            web.post("/api/data/download", handle_data_download),
            web.post("/api/data/verify", handle_data_verify),
            web.get("/api/data/progress/{id}", handle_data_download_progress),
        ]
    )
