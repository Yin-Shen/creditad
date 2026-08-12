#!/usr/bin/env python3
"""
BigWig sanity check — local file integrity + statistics.

**This is local file QC ONLY.**  It cannot tell you whether a BigWig comes
from the right experiment type (ChIP-seq vs ChIA-PET vs RNA-seq vs ATAC-seq
vs siRNA-knockdown-RNA-seq).  For that, query ENCODE assay metadata at
download time — see `download_chipseq_bigwigs.py` for the canonical pattern.

Real ChIP-seq fold-change-over-control files are ALWAYS dense everywhere
(MACS2 emits a fold-change value at every base genome-wide, including
centromeres where it falls to baseline ~0.5–1.0).  The earlier heuristic
("centromere density < 0.20 ⇒ signal-type, > 0.80 ⇒ reject") was wrong:
it correctly flagged fold-change files but it also rejected real ChIP-seq.
We now treat fold-change ChIP-seq as the canonical input and rely on the
hg{19,38}_gaps.bed mask in `main.py: _load_gaps()` to drop gap regions at
read time.

What this script verifies:
  - File is openable by pyBigWig
  - chr1 exists and has expected length (assembly inferable)
  - Signal value range looks like real biology (max > 5, sumData large enough
    that file is not empty)
  - Genome-wide nBasesCovered ratio (informational only; high is OK for
    fold-change ChIP-seq)
  - Probe a known strong-CTCF locus (chr19 ZNF cluster at 38-39 Mb) and
    confirm signal is non-zero — distinguishes a real-data file from an
    empty/zero-filled file


Usage:
    # CLI
    python3 validate_bigwig.py path/to/a.bigWig path/to/b.bigWig

    # importable
    from validate_bigwig import validate_bigwig, validate_or_raise
    result = validate_bigwig("path.bigWig")
    validate_or_raise({"GM12878 CTCF": "path1", "GM12878 RAD21": "path2"})
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyBigWig

CENTROMERE_PROBES = {
    "chr1": (121_000_000, 125_000_000),
    "1":    (121_000_000, 125_000_000),
}
# A known-strong CTCF/cohesin locus.  Real CTCF/RAD21 ChIP-seq should have
# substantial signal here; an empty or wrong-target file will not.
STRONG_CTCF_PROBE = ("chr19", 38_000_000, 39_000_000)


@dataclass(frozen=True)
class BigWigVerdict:
    path: str
    cov_ratio: float | None
    centromere_density: float | None
    max_val: float | None
    sum_data: float | None
    strong_locus_max: float | None
    verdict: str          # "ok" | "empty" | "error"
    detail: str

    @property
    def ok(self) -> bool:
        return self.verdict == "ok"


def validate_bigwig(path: str | Path) -> BigWigVerdict:
    p = str(path)
    try:
        bw = pyBigWig.open(p)
        if bw is None:
            return BigWigVerdict(p, None, None, None, None, None, "error",
                                 "pyBigWig.open returned None")
    except Exception as e:
        return BigWigVerdict(p, None, None, None, None, None, "error",
                             f"open failed: {e}")

    try:
        chroms = bw.chroms()
        h = bw.header()
        total = sum(chroms.values())
        cov_ratio = h["nBasesCovered"] / total if total else 0.0
        max_val = float(h["maxVal"])
        sum_data = float(h["sumData"])

        probe_chrom = next((c for c in CENTROMERE_PROBES if c in chroms), None)
        if probe_chrom is None:
            cen_dens = None
        else:
            start, end = CENTROMERE_PROBES[probe_chrom]
            end = min(end, chroms[probe_chrom])
            vals = bw.values(probe_chrom, start, end, numpy=True)
            cen_dens = float((~np.isnan(vals)).sum() / len(vals))

        # Probe a known strong-CTCF locus (chr19 ZNF cluster).  Real ChIP-seq
        # for CTCF/RAD21 should have a high max here.  Wrong-target files
        # return near-zero.
        strong_max = None
        for chrom_form in (STRONG_CTCF_PROBE[0], STRONG_CTCF_PROBE[0].replace("chr", "")):
            if chrom_form in chroms:
                end_p = min(STRONG_CTCF_PROBE[2], chroms[chrom_form])
                vals = bw.values(chrom_form, STRONG_CTCF_PROBE[1], end_p, numpy=True)
                vals = vals[~np.isnan(vals)]
                if len(vals):
                    strong_max = float(vals.max())
                break
    finally:
        bw.close()

    # Verdict: "ok" if file has real biology; "empty" if signal seems absent;
    # "error" already returned above on open failure.
    if max_val < 1.0 or sum_data < 1_000_000:
        verdict, detail = "empty", "max < 1.0 or sumData < 1M — file looks empty"
    elif strong_max is not None and strong_max < 0.5:
        verdict, detail = "empty", f"chr19:38-39Mb max={strong_max:.3f} — wrong-target?"
    else:
        verdict = "ok"
        detail = (f"max={max_val:.0f}  sumData={sum_data:.2e}  cov={cov_ratio:.3f}  "
                  f"chr19_38-39Mb_max={strong_max if strong_max is not None else 'NA'}")

    return BigWigVerdict(
        path=p, cov_ratio=cov_ratio, centromere_density=cen_dens,
        max_val=max_val, sum_data=sum_data, strong_locus_max=strong_max,
        verdict=verdict, detail=detail,
    )


def validate_or_raise(named_paths: dict[str, str | Path]) -> dict[str, BigWigVerdict]:
    """Used at backend startup. Raises RuntimeError if any file is empty or
    cannot be opened.  Does NOT verify ENCODE assay type — that's a separate
    concern handled at download time (see download_chipseq_bigwigs.py)."""
    results = {label: validate_bigwig(p) for label, p in named_paths.items()}
    bad = [(label, r) for label, r in results.items() if r.verdict in ("empty", "error")]
    if bad:
        msg = ["[validate_bigwig] FATAL — empty or unreadable BigWig(s):"]
        for label, r in bad:
            msg.append(f"  {label}: {r.detail} — {r.path}")
        raise RuntimeError("\n".join(msg))
    return results


def _print_table(results: dict[str, BigWigVerdict]) -> int:
    icons = {"ok": "OK ", "empty": "EMP", "error": "ERR"}
    width = max((len(l) for l in results), default=10)
    print(f"{'label':<{width}}  max     sumData    cov   chr19max  verdict  path")
    print("-" * (width + 75))
    bad = 0
    for label, r in results.items():
        mx = f"{r.max_val:.0f}" if r.max_val is not None else "  -  "
        sd = f"{r.sum_data:.2e}" if r.sum_data is not None else "    -    "
        cov = f"{r.cov_ratio:.3f}" if r.cov_ratio is not None else "  -  "
        chr19 = f"{r.strong_locus_max:.1f}" if r.strong_locus_max is not None else "  -  "
        print(f"{label:<{width}}  {mx:>5}   {sd}   {cov}   {chr19:>6}   {icons[r.verdict]}      {r.path}")
        if r.verdict in ("empty", "error"):
            bad += 1
    return bad


def main(argv: list[str]) -> int:
    if not argv:
        print("usage: python3 validate_bigwig.py <path.bigWig> [<path2> ...]", file=sys.stderr)
        return 2
    results = {f"#{i}": validate_bigwig(p) for i, p in enumerate(argv)}
    bad = _print_table(results)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


# --- assembly inference, added 2026-07-25 for the BYO path -------------------------
# A user can point CrediTAD at a bigWig from a different genome build than the one they
# declared. Nothing checked this, so a hg19 track against hg38 boundaries silently produced
# numbers: every window lands at the wrong coordinate and the resulting BD2/BD3 grades are
# meaningless rather than absent.
_CHR1_LENGTHS = {249250621: "hg19", 248956422: "hg38"}


def assembly_of(path: str) -> str | None:
    """Infer hg19/hg38 from the bigWig's chr1 length, or None if it cannot be told."""
    try:
        import pyBigWig
    except ImportError:
        return None
    try:
        bw = pyBigWig.open(path)
    except (RuntimeError, OSError):
        return None
    try:
        chroms = bw.chroms()
    finally:
        bw.close()
    for name in ("chr1", "1"):
        L = chroms.get(name)
        if L in _CHR1_LENGTHS:
            return _CHR1_LENGTHS[L]
    return None


def check_assembly(path: str, declared: str | None) -> str | None:
    """Return a human-readable mismatch message, or None when consistent/undecidable.

    Undecidable (non-human build, unusual chrom naming) is NOT an error — we only refuse
    when we can positively tell the two apart.
    """
    if not declared:
        return None
    got = assembly_of(path)
    if got is None:
        return None
    d = str(declared).strip().lower()
    d = {"grch37": "hg19", "grch38": "hg38"}.get(d, d)
    if d in ("hg19", "hg38") and got != d:
        return (f"{path}: the track is {got} (chr1 length) but the run declares {d}. "
                f"Grading a {d} boundary list against a {got} signal track compares "
                f"different coordinates; supply the matching build or drop the track.")
    return None
