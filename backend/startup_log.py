"""Structured startup logger for CrediTAD.

Emits machine-parseable lines that the Electron splash window picks up.
Format: ``[STARTUP|<step>|<status>|<message>]``

Steps used:
    env       application environment (Python version, platform, paths)
    deps      package import phase
    db        database init
    port      TCP port reservation
    listen    HTTP server bind
    ready     fully operational

Statuses:
    info      neutral information (always shown)
    active    step is in progress
    done      step finished successfully
    warn      non-fatal issue
    error     fatal issue (renderer shows the retry card)
"""
from __future__ import annotations

import os
import platform
import sys
import time
from contextlib import contextmanager
from typing import Iterator

_T0 = time.monotonic()


def _flush(line: str) -> None:
    print(line, flush=True)
    sys.stdout.flush()


def emit(step: str, status: str, message: str) -> None:
    elapsed = time.monotonic() - _T0
    _flush(f"[STARTUP|{step}|{status}|{message}|{elapsed:.2f}]")


def banner() -> None:
    emit(
        "env",
        "info",
        f"Python {platform.python_version()} ({platform.python_implementation()}) on {platform.platform()}",
    )
    emit("env", "info", f"Executable: {sys.executable}")
    emit("env", "info", f"Working directory: {os.getcwd()}")
    emit(
        "env",
        "info",
        f"Process: PID={os.getpid()} CPU count={os.cpu_count()}",
    )


@contextmanager
def step(name: str, message: str) -> Iterator[None]:
    emit(name, "active", message)
    t0 = time.monotonic()
    try:
        yield
    except Exception as exc:
        emit(name, "error", f"{message} - {type(exc).__name__}: {exc}")
        raise
    else:
        emit(name, "done", f"{message} ({(time.monotonic() - t0) * 1000:.0f} ms)")


def ready(port: int, degraded: dict[str, str] | None = None) -> None:
    """Final startup line. ``degraded`` names subsystems that failed to come up.

    Announcing "All systems operational" after a subsystem failed is how a backend
    whose entire /api/tadvci/* route table is missing still reported a clean start.
    """
    elapsed = time.monotonic() - _T0
    if degraded:
        detail = "; ".join(f"{k}: {v}" for k, v in sorted(degraded.items()))
        emit("ready", "warn",
             f"Started DEGRADED on port {port} (total {elapsed:.2f}s) - {detail}")
        return
    emit("ready", "done", f"All systems operational on port {port} (total {elapsed:.2f}s)")


def fatal(step_name: str, message: str) -> None:
    emit(step_name, "error", message)
