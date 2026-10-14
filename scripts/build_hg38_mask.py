#!/usr/bin/env python3
"""build_hg38_mask.py -- build the bundled hg38 mask used by `creditad annotate --mask hg38`.

The mask lists the hg38 regions where a ChIP readout is not a measurement: assembly gaps,
centromeric and heterochromatic cytobands, and the ENCODE exclusion list. It is the region
set the v1.1 evidence tables (data/multicell_v1.1/) were masked with, written as one BED4
file whose fourth column names the source of each interval:

    assembly_gap              rows of the UCSC hg38 `gap` table
    cytoband_acen             UCSC hg38 `cytoBand` rows with gieStain acen
    cytoband_gvar             ... gvar
    cytoband_stalk            ... stalk
    encode_blacklist          ENCODE exclusion list ENCFF356LFX (GRCh38)

Only chr1-22, chrX and chrY are kept, as in the v1.1 masking. Intervals are half-open
(BED) and are written as they appear in the sources, not merged, so the label of every
interval stays readable; a position is masked when it lies in any interval.

Inputs (public, read-only downloads; see data/masks/README.md for dates and checksums):
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
    https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz

Usage:
    python scripts/build_hg38_mask.py --gap gap.txt.gz --cytoband cytoBand.txt.gz \\
        --blacklist ENCFF356LFX.bed.gz --out data/masks/hg38_mask.bed.gz

The output is gzip-compressed with a zero timestamp, so the same inputs give the same bytes.
"""
from __future__ import annotations

import argparse
import gzip
import io
import re
from pathlib import Path

KEEP = re.compile(r"^chr([0-9]+|X|Y)$")
CYTO_STAINS = ("acen", "gvar", "stalk")


def _open_text(path: Path):
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return io.StringIO(raw.decode())


def _chrom_key(chrom: str) -> tuple[int, str]:
    tail = chrom[3:]
    return (int(tail), "") if tail.isdigit() else (100, tail)


def read_gap(path: Path) -> list[tuple[str, int, int, str]]:
    # UCSC gap table: bin chrom chromStart chromEnd ix n size type bridge
    out = []
    for line in _open_text(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and KEEP.match(f[1]):
            out.append((f[1], int(f[2]), int(f[3]), "assembly_gap"))
    return out


def read_cytoband(path: Path) -> list[tuple[str, int, int, str]]:
    # UCSC cytoBand table: chrom chromStart chromEnd name gieStain
    out = []
    for line in _open_text(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5 and KEEP.match(f[0]) and f[4] in CYTO_STAINS:
            out.append((f[0], int(f[1]), int(f[2]), f"cytoband_{f[4]}"))
    return out


def read_blacklist(path: Path) -> list[tuple[str, int, int, str]]:
    out = []
    for line in _open_text(path):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        if KEEP.match(f[0]):
            out.append((f[0], int(f[1]), int(f[2]), "encode_blacklist"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gap", type=Path, required=True, help="UCSC hg38 gap.txt.gz")
    ap.add_argument("--cytoband", type=Path, required=True, help="UCSC hg38 cytoBand.txt.gz")
    ap.add_argument("--blacklist", type=Path, required=True, help="ENCODE ENCFF356LFX.bed.gz")
    ap.add_argument("--out", type=Path, required=True, help="output .bed.gz")
    args = ap.parse_args()

    rows = read_gap(args.gap) + read_cytoband(args.cytoband) + read_blacklist(args.blacklist)
    bad = [r for r in rows if not 0 <= r[1] < r[2]]
    if bad:
        raise SystemExit(f"interval with start >= end: {bad[0]}")
    rows = sorted(set(rows), key=lambda r: (_chrom_key(r[0]), r[1], r[2], r[3]))
    text = "".join(f"{c}\t{s}\t{e}\t{lab}\n" for c, s, e, lab in rows).encode()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0) as gz:
        gz.write(text)
    args.out.write_bytes(buf.getvalue())

    counts: dict[str, int] = {}
    for r in rows:
        counts[r[3]] = counts.get(r[3], 0) + 1
    print(f"{len(rows)} intervals -> {args.out}")
    for lab in sorted(counts):
        print(f"  {lab:18s} {counts[lab]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
