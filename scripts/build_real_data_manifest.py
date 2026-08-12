"""Build the real_data_manifest.tsv for the 5 fixed mcool files.

Probes each file for available resolutions, balance presence, target chrom
availability, and selected regions (Level A / B / C).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

DATASETS = [
    {"dataset_id": "GM12878", "cell_line": "GM12878", "organism": "human",
     "file_path": "$CREDITAD_DATA_ROOT/gm12878/hic/4DNFIXP4QG5B.mcool",
     "role": "primary"},
    {"dataset_id": "K562", "cell_line": "K562", "organism": "human",
     "file_path": "$CREDITAD_DATA_ROOT/k562/hic/4DNFI18UHVRO.mcool",
     "role": "secondary_human"},
    {"dataset_id": "HepG2", "cell_line": "HepG2", "organism": "human",
     "file_path": "$CREDITAD_DATA_ROOT/hepg2/hic/4DNFIS6HAUPP.mcool",
     "role": "secondary_human"},
    {"dataset_id": "IMR90", "cell_line": "IMR90", "organism": "human",
     "file_path": "$CREDITAD_DATA_ROOT/imr90/hic/4DNFIJTOIGOI.mcool",
     "role": "secondary_human"},
    {"dataset_id": "CH12_LX", "cell_line": "CH12-LX", "organism": "mouse",
     "file_path": "$CREDITAD_DATA_ROOT/ch12_lx/hic/4DNFIEAVF99Y.mcool",
     "role": "secondary_mouse"},
]

TARGET_RESOLUTIONS = [10_000, 25_000, 50_000, 100_000]


def probe_file(d: dict) -> dict:
    p = Path(d["file_path"])
    if not p.exists():
        return {**d, "file_exists": False, "file_size_bytes": 0,
                "cooler_open_success": False,
                "available_resolutions": "", "target_resolutions_present": "",
                "target_resolutions_missing": ",".join(str(r) for r in TARGET_RESOLUTIONS),
                "selected_level_a_resolutions": "",
                "selected_level_b_resolution": "",
                "selected_level_c_resolutions": "",
                "balance_available_by_resolution": "",
                "selected_matrix_mode_by_resolution": "",
                "chroms_available": "", "selected_level_a_region": "",
                "selected_level_b_region": "", "selected_level_c_region": "",
                "notes": "FILE MISSING"}

    rec = {**d, "file_exists": True, "file_size_bytes": p.stat().st_size}

    try:
        with h5py.File(p, "r") as f:
            available = sorted(int(k) for k in f["resolutions"].keys())
        rec["available_resolutions"] = ",".join(str(r) for r in available)
        rec["cooler_open_success"] = True
    except Exception as e:
        rec["cooler_open_success"] = False
        rec["available_resolutions"] = ""
        rec["notes"] = f"h5py failed: {e}"
        return rec

    target_present = [r for r in TARGET_RESOLUTIONS if r in available]
    target_missing = [r for r in TARGET_RESOLUTIONS if r not in available]
    rec["target_resolutions_present"] = ",".join(str(r) for r in target_present)
    rec["target_resolutions_missing"] = ",".join(str(r) for r in target_missing)

    # Test cooler open at one resolution and probe chroms
    try:
        clr = cooler.Cooler(f"{p}::resolutions/{target_present[0]}")
        chroms = clr.chromnames
    except Exception as e:
        rec["chroms_available"] = ""
        rec["notes"] = f"cooler open at {target_present[0]} failed: {e}"
        return rec

    rec["chroms_available"] = ",".join(chroms[:30] + (["..."] if len(chroms) > 30 else []))

    # Pick a chr19 (or fallback) for human; for mouse pick chr19 too if present, else chr12
    if d["organism"] == "human":
        candidates = ["chr19", "19", "chr18", "18"]
    else:
        candidates = ["chr19", "19", "chr12", "12", "chr11", "11"]
    selected_chrom = next((c for c in candidates if c in chroms), chroms[0])

    # Get chrom length
    sizes = dict(zip(clr.chromnames, clr.chromsizes))
    chrom_len = int(sizes[selected_chrom])

    # Level A: 10 Mb region centered or first 10 Mb that fits
    if chrom_len >= 40_000_000:
        a_start = 30_000_000
    else:
        a_start = max(0, chrom_len // 2 - 5_000_000)
    a_end = min(chrom_len, a_start + 10_000_000)
    rec["selected_level_a_region"] = f"{selected_chrom}:{a_start}-{a_end}"

    # Level B: >=20 Mb on the same chrom
    if chrom_len >= 20_000_000:
        b_start = 0
        b_end = min(chrom_len, 25_000_000)
    else:
        b_start = 0
        b_end = chrom_len
    rec["selected_level_b_region"] = f"{selected_chrom}:{b_start}-{b_end}"

    # Level C: GM12878 only — full chrom or >=30 Mb
    if d["dataset_id"] == "GM12878":
        rec["selected_level_c_region"] = f"{selected_chrom}:0-{chrom_len}"
    else:
        rec["selected_level_c_region"] = ""

    rec["selected_level_a_resolutions"] = ",".join(str(r) for r in target_present)
    rec["selected_level_b_resolution"] = "25000" if 25000 in target_present else (str(target_present[0]) if target_present else "")
    rec["selected_level_c_resolutions"] = (",".join(str(r) for r in target_present)
                                           if d["dataset_id"] == "GM12878" else "")

    # Balance availability check per resolution
    bal_avail = {}
    matrix_mode = {}
    for r in target_present:
        try:
            clr_r = cooler.Cooler(f"{p}::resolutions/{r}")
            bins_first = clr_r.bins()[:1000]
            has_balance = "weight" in bins_first.columns
            if has_balance:
                # also check it's not all NaN
                w = bins_first["weight"].dropna()
                has_balance = len(w) > 0
            bal_avail[r] = has_balance
            matrix_mode[r] = "balanced" if has_balance else "raw"
        except Exception:
            bal_avail[r] = False
            matrix_mode[r] = "raw"
    rec["balance_available_by_resolution"] = json.dumps({str(k): v for k, v in bal_avail.items()})
    rec["selected_matrix_mode_by_resolution"] = json.dumps({str(k): v for k, v in matrix_mode.items()})

    if not target_missing:
        rec["notes"] = "all 4 target resolutions available"
    else:
        rec["notes"] = f"missing target resolutions: {target_missing}"
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    rows = [probe_file(d) for d in DATASETS]
    df = pd.DataFrame(rows)
    cols = [
        "dataset_id", "cell_line", "organism", "file_path", "file_exists",
        "file_size_bytes", "cooler_open_success", "available_resolutions",
        "target_resolutions_present", "target_resolutions_missing",
        "selected_level_a_resolutions", "selected_level_b_resolution",
        "selected_level_c_resolutions", "balance_available_by_resolution",
        "selected_matrix_mode_by_resolution", "chroms_available",
        "selected_level_a_region", "selected_level_b_region",
        "selected_level_c_region", "role", "notes",
    ]
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    df = df[cols]
    df.to_csv(out / "real_data_manifest.tsv", sep="\t", index=False)
    print(df[["dataset_id", "file_exists", "target_resolutions_present",
              "selected_level_a_region", "selected_matrix_mode_by_resolution"]].to_string(index=False))


if __name__ == "__main__":
    main()
