#!/usr/bin/env python3
"""Extract chr7-only BigWig subsets for Portable / example_data shipping.

Full ENCODE fold-change tracks are hundreds of MB each. Demo heatmaps use
chr7.mcool; shipping chr7-only CTCF/RAD21 keeps BD2/BD3 browser tracks working
without multi-GB packages.

Default output: example_data/chr7_tracks/ (+ copy to backend/extra_mode/chr7_tracks/).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pyBigWig

ROOT = Path(__file__).resolve().parents[1]
SHARED = Path(os.environ.get("CREDITAD_DATA_ROOT", ROOT / "data" / "raw"))

SOURCES = [
    {
        "cell": "GM12878",
        "role": "ctcf",
        "src": SHARED / "gm12878/ctcf/ENCFF734CUT.bigWig",
        "out": "GM12878_CTCF_ENCFF734CUT_hg38_chr7.bigWig",
    },
    {
        "cell": "GM12878",
        "role": "rad21",
        "src": SHARED / "gm12878/rad21/ENCFF571ZJJ.bigWig",
        "out": "GM12878_RAD21_ENCFF571ZJJ_hg38_chr7.bigWig",
    },
    {
        "cell": "IMR90",
        "role": "ctcf",
        "src": SHARED / "imr90/ctcf/ENCFF105FHL.bigWig",
        "out": "IMR90_CTCF_ENCFF105FHL_hg38_chr7.bigWig",
    },
    {
        "cell": "IMR90",
        "role": "rad21",
        "src": SHARED / "imr90/rad21/ENCFF048PZI.bigWig",
        "out": "IMR90_RAD21_ENCFF048PZI_hg38_chr7.bigWig",
    },
    {
        "cell": "HepG2",
        "role": "ctcf",
        "src": SHARED / "hepg2/ctcf/ENCFF357NFO.bigWig",
        "out": "HepG2_CTCF_ENCFF357NFO_hg38_chr7.bigWig",
    },
    {
        "cell": "HepG2",
        "role": "rad21",
        "src": SHARED / "hepg2/rad21/ENCFF972ODZ.bigWig",
        "out": "HepG2_RAD21_ENCFF972ODZ_hg38_chr7.bigWig",
    },
]


def _pick_chr7(chroms: dict) -> str:
    if "chr7" in chroms:
        return "chr7"
    if "7" in chroms:
        return "7"
    raise KeyError("no chr7/7 in bigWig chroms")


def extract_chr7(src: Path, dest: Path, bin_size: int = 250) -> None:
    """Write chr7-only bigWig: mean over fixed bins; skip near-zero bins.

    Processes the chromosome in 5 Mb windows so libBigWig never gets multi-million
    bin queries (those are extremely slow).
    """
    print(f"  open {src}", flush=True)
    bw = pyBigWig.open(str(src))
    if bw is None:
        raise RuntimeError(f"cannot open {src}")
    chrom = _pick_chr7(bw.chroms())
    length = int(bw.chroms()[chrom])
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".tmp.bigWig")
    if tmp.exists():
        tmp.unlink()
    out = pyBigWig.open(str(tmp), "w")
    # Cooler .mcool chrom names are bare ("7"); keep both spellings out of the
    # file by writing bare "7" so Hi-C + BigWig share one selector in the UI.
    # bigwig_adapter also aliases chr7 ↔ 7 at read time.
    out_chrom = "7"
    out.addHeader([(out_chrom, length)])

    window = 5_000_000
    starts: list[int] = []
    ends: list[int] = []
    vals: list[float] = []
    n_kept = 0
    for w0 in range(0, length, window):
        w1 = min(length, w0 + window)
        nbin = max(1, (w1 - w0) // bin_size)
        stats = bw.stats(chrom, w0, w1, nBins=nbin, type="mean")
        if not stats:
            continue
        step = (w1 - w0) / nbin
        for i, v in enumerate(stats):
            if v is None:
                continue
            fv = float(v)
            if abs(fv) < 1e-6:
                continue
            s = int(w0 + i * step)
            e = int(w0 + (i + 1) * step)
            if e <= s:
                e = s + 1
            starts.append(s)
            ends.append(e)
            vals.append(fv)
            n_kept += 1
            if len(vals) >= 50_000:
                out.addEntries([out_chrom] * len(vals), starts, ends=ends, values=vals)
                starts, ends, vals = [], [], []
        if w0 and w0 % (window * 5) == 0:
            print(f"  … {w0 / 1e6:.0f} Mb / {length / 1e6:.0f} Mb", flush=True)
    if vals:
        out.addEntries([out_chrom] * len(vals), starts, ends=ends, values=vals)
    out.close()
    bw.close()
    tmp.replace(dest)
    print(
        f"  wrote {dest.name}  {dest.stat().st_size / 1e6:.1f} MB  non-zero bins={n_kept}",
        flush=True,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "example_data" / "chr7_tracks")
    ap.add_argument("--bin-size", type=int, default=250, help="mean-bin width in bp")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    done: dict[tuple[str, str], Path] = {}
    for spec in SOURCES:
        src = Path(spec["src"])
        dest = out_dir / spec["out"]
        if not src.is_file():
            print(f"SKIP missing source: {src}", file=sys.stderr, flush=True)
            continue
        if dest.is_file() and dest.stat().st_size > 1000 and not args.force:
            print(f"exists {dest.name} ({dest.stat().st_size / 1e6:.1f} MB) — skip", flush=True)
            done[(spec["cell"], spec["role"])] = dest
            continue
        print(f"extract {spec['cell']} {spec['role']} …", flush=True)
        extract_chr7(src, dest, bin_size=args.bin_size)
        done[(spec["cell"], spec["role"])] = dest

    if ("GM12878", "ctcf") in done:
        shutil.copy2(done[("GM12878", "ctcf")], out_dir / "CTCF_chr7.bigWig")
    if ("GM12878", "rad21") in done:
        shutil.copy2(done[("GM12878", "rad21")], out_dir / "RAD21_chr7.bigWig")

    extra = ROOT / "backend" / "extra_mode" / "chr7_tracks"
    extra.mkdir(parents=True, exist_ok=True)
    for p in out_dir.glob("*.bigWig"):
        shutil.copy2(p, extra / p.name)
    print(f"copied → {extra}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
