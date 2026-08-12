"""Cross-platform BigWig reader adapter.

On Windows uses `winbbi` (pure DLL, the only reliable option on Win).
On Linux / macOS uses `pyBigWig` (libBigWig C bindings).

API matches `winbbi.BigWigReader` so `main.py` can stay unchanged:
    with BigWigReader() as r:
        r.open(path)
        r.get_chromosomes()            -> list[str]
        r.get_chrom_size(name)         -> int
        r.read_zoom_signal(chrom=..., start=..., end=..., num_bins=..., use_closest=True)
                                       -> list[float]
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import List, Optional

_IS_WIN = sys.platform == "win32"

def _chrom_aliases(name: str) -> list[str]:
    """Return candidate spellings: '7' ↔ 'chr7' (cooler vs ENCODE)."""
    n = str(name or "").strip()
    if not n:
        return []
    out = [n]
    if n.lower().startswith("chr"):
        bare = n[3:]
        if bare:
            out.append(bare)
            out.append(f"chr{bare}")
    else:
        out.append(f"chr{n}")
        out.append(f"Chr{n}")
    # de-dupe preserve order
    seen = set()
    uniq = []
    for x in out:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    return uniq


if _IS_WIN:
    # winbbi's Go DLL can corrupt its own handle or hard-fault, which aborts the whole
    # process (Go panics exit 2 and are not catchable from Python). It therefore runs in
    # a child process; see bigwig_helper.py for the two reproduced failures. Everything
    # below is the client side of that one-JSON-object-per-line protocol.
    import json as _json
    import subprocess as _sp
    import threading as _th

    _HELPER = Path(__file__).resolve().parent / "bigwig_helper.py"
    _HEADER_CACHE: dict = {}
    _LOCK = _th.Lock()
    _PROC = None


    def _spawn():
        return _sp.Popen(
            [sys.executable, "-u", str(_HELPER)],
            stdin=_sp.PIPE, stdout=_sp.PIPE, stderr=_sp.DEVNULL,
            text=True, encoding="utf-8",
            cwd=str(Path(__file__).resolve().parent),
        )

    def _ask(payload: dict) -> dict:
        """One request/response. Respawns the helper if it died; never raises Go panics."""
        global _PROC
        with _LOCK:
            for attempt in (1, 2):        # one free retry across a helper restart
                if _PROC is None or _PROC.poll() is not None:
                    _PROC = _spawn()
                try:
                    _PROC.stdin.write(_json.dumps(payload) + "\n")
                    _PROC.stdin.flush()
                    line = _PROC.stdout.readline()
                except (BrokenPipeError, OSError, ValueError):
                    line = ""
                if not line:              # helper exited (crash, or hit its request cap)
                    try:
                        _PROC.kill()
                    except Exception:
                        pass
                    _PROC = None
                    if attempt == 1:
                        continue
                    raise RuntimeError("BigWig helper process exited without a reply")
                try:
                    msg = _json.loads(line)
                except ValueError as exc:
                    raise RuntimeError(f"BigWig helper sent malformed output: {exc}") from None
                if not msg.get("ok"):
                    raise RuntimeError(msg.get("error") or "BigWig helper reported failure")
                return msg.get("result") or {}
        raise RuntimeError("BigWig helper unreachable")

    def _header(path: str) -> dict:
        """{raw_chrom: length}. Cached in THIS process: headers do not change, and
        caching them keeps ordinary browsing to one helper round trip per signal read."""
        key = os.path.abspath(path)
        got = _HEADER_CACHE.get(key)
        if got is None:
            got = _ask({"op": "header", "path": path}).get("chroms") or {}
            _HEADER_CACHE[key] = got
        return got

    class BigWigReader:
        """winbbi-compatible surface backed by the helper process."""

        def __init__(self) -> None:
            self._path: Optional[str] = None

        def __enter__(self) -> "BigWigReader":
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            self.close()

        def open(self, path: str) -> None:
            # Absolute: the helper runs with its own cwd, so a relative path handed
            # straight through would resolve against the wrong directory.
            self._path = os.path.abspath(path)
            _header(self._path)           # surfaces an unreadable file immediately

        def close(self) -> None:
            # Nothing to close: the handle lives (and dies) with the helper process.
            self._path = None

        def get_chromosomes(self) -> List[str]:
            return list(_header(self._path).keys()) if self._path else []

        def chroms(self) -> dict:
            return dict(_header(self._path)) if self._path else {}

        def _resolve_chrom(self, name: str) -> Optional[str]:
            if not self._path:
                return None
            available = _header(self._path)
            for cand in _chrom_aliases(name):
                if cand in available:
                    return cand
            return None

        def get_chrom_size(self, name: str) -> int:
            raw = self._resolve_chrom(name)
            return int(_header(self._path).get(raw, 0)) if raw else 0

        def read_zoom_signal(
            self,
            chrom: str,
            start: int,
            end: int,
            num_bins: int,
            use_closest: bool = True,
        ) -> List[float]:
            if not self._path:
                return [0.0] * num_bins
            raw = self._resolve_chrom(chrom)
            if raw is None:
                return [0.0] * num_bins
            res = _ask({
                "op": "zoom", "path": self._path, "chrom": raw,
                "start": int(start), "end": int(end), "bins": int(num_bins),
                "closest": bool(use_closest),
            })
            return [float(v) for v in (res.get("values") or [])]


else:
    import pyBigWig

    class BigWigReader:
        """pyBigWig-backed reader with winbbi-compatible surface."""

        def __init__(self) -> None:
            self._bw = None

        def __enter__(self) -> "BigWigReader":
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            self.close()

        def open(self, path: str) -> None:
            self._bw = pyBigWig.open(path)
            if self._bw is None:
                raise RuntimeError(f"pyBigWig failed to open: {path}")

        def close(self) -> None:
            if self._bw is not None:
                try:
                    self._bw.close()
                except Exception:
                    pass
                self._bw = None

        def get_chromosomes(self) -> List[str]:
            if self._bw is None:
                return []
            return list(self._bw.chroms().keys())

        def chroms(self) -> dict:
            if self._bw is None:
                return {}
            return dict(self._bw.chroms())

        def _resolve_chrom(self, name: str) -> Optional[str]:
            if self._bw is None:
                return None
            have = self._bw.chroms() or {}
            for cand in _chrom_aliases(name):
                if cand in have:
                    return cand
            lower = {k.lower(): k for k in have}
            for cand in _chrom_aliases(name):
                if cand.lower() in lower:
                    return lower[cand.lower()]
            return None

        def get_chrom_size(self, name: str) -> int:
            if self._bw is None:
                return 0
            raw = self._resolve_chrom(name)
            if raw is None:
                return 0
            return int(self._bw.chroms().get(raw, 0))

        def read_zoom_signal(
            self,
            chrom: str,
            start: int,
            end: int,
            num_bins: int,
            use_closest: bool = True,  # kept for API parity, ignored
        ) -> List[float]:
            if self._bw is None:
                return [0.0] * num_bins

            raw = self._resolve_chrom(chrom)
            if raw is None:
                return [0.0] * num_bins
            chrom_size = int(self._bw.chroms().get(raw, 0))
            if chrom_size == 0:
                return [0.0] * num_bins

            start = max(0, int(start))
            end = min(int(end), chrom_size)
            if end <= start:
                return [0.0] * num_bins

            try:
                vals = self._bw.stats(raw, start, end, type="mean", nBins=num_bins)
            except Exception:
                return [0.0] * num_bins

            return [0.0 if v is None else float(v) for v in vals]
