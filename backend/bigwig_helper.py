"""Out-of-process BigWig reader for Windows.

Why this exists
---------------
`winbbi` ships a Go 1.25 c-shared DLL. Two behaviours make it unsafe to call from the
server process:

* a whole-chromosome zoom read (`start=0, end=chrom_length, num_bins=2000` — the baseline
  probe in `main.py`) leaves the handle corrupted. A following call either returns
  nonsense (`get_chromosomes()` -> `[]` where it had returned `['7']`) or faults:

      unexpected fault address 0xffffffffffffffff
      fatal error: fault  [signal 0xc0000005]
      go-bigwig/gobigwig.CloseBigWig(...)   api.go:126
      go-bigwig/gobigwig.BigWigGetChromosomes(...)  api.go:446

* `BigWigClose()` on a handle that has served a zoom read faults the same way, and
  `winbbi.BigWigReader.__del__` calls `close()`, so garbage collection is enough to
  trigger it.

A Go panic is not a Python exception: it aborts the process with exit code 2, so no
`try/except` in the backend can contain it. On a user's Windows 10 machine this showed up
as a backend that logged "All systems operational" and then vanished, leaving the UI on
"Could not read the shipped packages".

So the DLL runs here, in a child process, and the server talks to it over stdin/stdout.
If it faults, this process dies and the parent raises an ordinary Python error that the
existing handlers already cope with — a blank track instead of a dead application.

Protocol: one JSON object per line in, one JSON object per line out.
    {"op": "header", "path": ...}
    {"op": "zoom", "path": ..., "chrom": ..., "start": ..., "end": ..., "bins": ...,
     "closest": true}
    -> {"ok": true, "result": ...} | {"ok": false, "error": "..."}

A reader is opened fresh per request and never closed (closing is what crashes). The
handles that leaks are reclaimed by exiting after MAX_REQUESTS; the parent respawns.
"""
from __future__ import annotations

import json
import os
import sys

MAX_REQUESTS = int(os.environ.get("CREDITAD_BIGWIG_HELPER_MAX_REQUESTS", "32"))


# Readers are kept alive for the life of this process. `winbbi.BigWigReader.__del__`
# calls `close()`, and closing a handle that served a zoom read faults — so letting one
# be garbage-collected is enough to kill this helper. We both retain a reference AND
# shadow `close` on the instance, so neither collection nor interpreter shutdown can
# reach the crashing call.
_KEEP: list = []


def _reader(path: str):
    import winbbi

    r = winbbi.BigWigReader()
    r.close = lambda: None          # instance attribute shadows the class method
    r.open(path)
    _KEEP.append(r)
    return r


def _header(path: str) -> dict:
    r = _reader(path)
    names = list(r.get_chromosomes())
    return {"chroms": {n: int(r.get_chrom_size(n)) for n in names}}


def _zoom(req: dict) -> dict:
    r = _reader(req["path"])
    vals = r.read_zoom_signal(
        chrom=req["chrom"],
        start=int(req["start"]),
        end=int(req["end"]),
        num_bins=int(req["bins"]),
        use_closest=bool(req.get("closest", True)),
    )
    return {"values": [float(v) for v in vals]}


def main() -> int:
    served = 0
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
            op = req.get("op")
            if op == "header":
                out = {"ok": True, "result": _header(req["path"])}
            elif op == "zoom":
                out = {"ok": True, "result": _zoom(req)}
            elif op == "ping":
                out = {"ok": True, "result": {"pong": True}}
            else:
                out = {"ok": False, "error": f"unknown op: {op!r}"}
        except Exception as exc:                        # noqa: BLE001 - reported to parent
            out = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        sys.stdout.write(json.dumps(out) + "\n")
        sys.stdout.flush()
        served += 1
        if served >= MAX_REQUESTS:
            # Exit cleanly so the leaked (never-closed) handles are reclaimed by the OS.
            # The parent notices EOF and starts a fresh helper for the next request.
            return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
