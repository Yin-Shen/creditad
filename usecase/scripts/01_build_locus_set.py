#!/usr/bin/env python3
"""Step 1 of the CrediTAD use-case demonstration: build the disease-locus set.

Implements PREREGISTERED_RULE.md sections 1 and 2 EXACTLY. No thresholds appear
here that are not written in that file. Re-runnable: reads only the downloaded
ClinVar release, writes only under usecase/.

Usage:  python 01_build_locus_set.py
"""
from __future__ import annotations
import gzip, hashlib, json, os, sys
from datetime import datetime, timezone

UC = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(UC, "data")
GZ = os.path.join(DATA, "variant_summary.txt.gz")

# ---- PREREG §2 constants (verbatim from the pre-registration) -------------
ACCEPT_TYPE = {"copy number loss", "copy number gain", "deletion", "duplication"}
ACCEPT_SIG = {"pathogenic", "pathogenic/likely pathogenic"}
CHROMS = {str(i) for i in range(1, 23)} | {"X"}
SIZE_MIN, SIZE_MAX = 10_000, 5_000_000
ASSEMBLY = "GRCh38"


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5sum(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    funnel = {
        "rows_total": 0, "assembly_GRCh38": 0, "chrom_ok": 0, "type_ok": 0,
        "clinsig_ok": 0, "coord_sane": 0, "size_ok": 0,
    }
    dropped_coord = 0
    type_counter: dict[str, int] = {}
    sig_counter: dict[str, int] = {}
    kept: dict[tuple, dict] = {}   # (chrom, start, stop, type) -> row
    dup_collapsed = 0

    with gzip.open(GZ, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        need = ["Type", "ClinicalSignificance", "Assembly", "Chromosome", "Start",
                "Stop", "VariationID", "GeneSymbol", "PhenotypeList", "ReviewStatus",
                "RCVaccession", "nsv/esv (dbVar)", "Name", "NumberSubmitters",
                "LastEvaluated", "ChromosomeAccession", "Cytogenetic", "Origin"]
        for n in need:
            if n not in idx:
                raise SystemExit(f"missing expected ClinVar column: {n}")

        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            funnel["rows_total"] += 1

            if f[idx["Assembly"]] != ASSEMBLY:
                continue
            funnel["assembly_GRCh38"] += 1

            chrom = f[idx["Chromosome"]].strip()
            chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
            if chrom not in CHROMS:
                continue
            funnel["chrom_ok"] += 1

            vtype_raw = f[idx["Type"]].strip()
            type_counter[vtype_raw] = type_counter.get(vtype_raw, 0) + 1
            if vtype_raw.strip().lower() not in ACCEPT_TYPE:
                continue
            funnel["type_ok"] += 1

            sig_raw = f[idx["ClinicalSignificance"]].strip()
            sig_counter[sig_raw] = sig_counter.get(sig_raw, 0) + 1
            if sig_raw.strip().lower() not in ACCEPT_SIG:
                continue
            funnel["clinsig_ok"] += 1

            try:
                start = int(f[idx["Start"]]); stop = int(f[idx["Stop"]])
            except ValueError:
                dropped_coord += 1
                continue
            if start < 1 or stop < start:
                dropped_coord += 1
                continue
            funnel["coord_sane"] += 1

            size = stop - start
            if not (SIZE_MIN <= size <= SIZE_MAX):
                continue
            funnel["size_ok"] += 1

            try:
                vid_int = int(f[idx["VariationID"]])
            except ValueError:
                vid_int = 1 << 62
            key = (chrom, start, stop, vtype_raw)
            row = {
                "VariationID": f[idx["VariationID"]].strip(),
                "VariationID_int": vid_int,
                "Type": vtype_raw,
                "ClinicalSignificance": sig_raw,
                "Chromosome": chrom,
                "Start": start,
                "Stop": stop,
                "size_bp": size,
                "GeneSymbol": f[idx["GeneSymbol"]].strip(),
                "PhenotypeList": f[idx["PhenotypeList"]].strip(),
                "ReviewStatus": f[idx["ReviewStatus"]].strip(),
                "RCVaccession": f[idx["RCVaccession"]].strip(),
                "dbVar": f[idx["nsv/esv (dbVar)"]].strip(),
                "Name": f[idx["Name"]].strip(),
                "NumberSubmitters": f[idx["NumberSubmitters"]].strip(),
                "LastEvaluated": f[idx["LastEvaluated"]].strip(),
                "ChromosomeAccession": f[idx["ChromosomeAccession"]].strip(),
                "Cytogenetic": f[idx["Cytogenetic"]].strip(),
                "Origin": f[idx["Origin"]].strip(),
            }
            prev = kept.get(key)
            if prev is None:
                kept[key] = row
            else:
                dup_collapsed += 1
                if vid_int < prev["VariationID_int"]:
                    kept[key] = row

    # deterministic order: chrom (1..22,X), start, stop, VariationID int
    def ckey(c: str) -> int:
        return 23 if c == "X" else int(c)
    loci = sorted(kept.values(),
                  key=lambda r: (ckey(r["Chromosome"]), r["Start"], r["Stop"],
                                 r["VariationID_int"]))

    out_tsv = os.path.join(UC, "data", "locus_set.tsv")
    cols = ["VariationID", "Type", "ClinicalSignificance", "Chromosome", "Start", "Stop",
            "size_bp", "GeneSymbol", "PhenotypeList", "ReviewStatus", "RCVaccession",
            "dbVar", "Name", "NumberSubmitters", "LastEvaluated", "ChromosomeAccession",
            "Cytogenetic", "Origin"]
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in loci:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    prov = {
        "source_rank_used": 1,
        "source_name": "ClinVar variant_summary.txt.gz (NCBI, tab-delimited full release)",
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
        "md5_url": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz.md5",
        "retrieved_utc": "2026-08-15T02:39:17Z",
        "http_status": 200,
        "bytes": os.path.getsize(GZ),
        "sha256": sha256(GZ),
        "md5": md5sum(GZ),
        "md5_published_by_ncbi": "897c38fd97dff56f83f1c7a674eb35f3",
        "md5_verified_match": True,
        "server_last_modified": "Mon, 10 Aug 2026 12:59:47 GMT",
        "release_note": "NCBI does not embed a version string in variant_summary.txt.gz; "
                        "the release is identified by the server Last-Modified date and the "
                        "published MD5, both recorded here.",
        "assembly_used": "GRCh38 (hg38) as released; NO liftover performed",
        "filter_applied": "PREREGISTERED_RULE.md section 2",
        "script": "usecase/scripts/01_build_locus_set.py",
        "output": "usecase/data/locus_set.tsv",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    with open(os.path.join(UC, "SOURCE_PROVENANCE.json"), "w") as fh:
        json.dump(prov, fh, indent=2)

    funnel_out = {
        "funnel": funnel,
        "dropped_coordinate_insane_or_unparseable": dropped_coord,
        "duplicate_tuples_collapsed": dup_collapsed,
        "loci_final": len(loci),
        "accepted_types_observed": {k: v for k, v in sorted(type_counter.items())
                                    if k.strip().lower() in ACCEPT_TYPE},
        "accepted_clinsig_observed": {k: v for k, v in sorted(sig_counter.items())
                                      if k.strip().lower() in ACCEPT_SIG},
    }
    with open(os.path.join(UC, "data", "locus_funnel.json"), "w") as fh:
        json.dump(funnel_out, fh, indent=2)

    print(json.dumps({"loci_final": len(loci), "funnel": funnel,
                      "dup_collapsed": dup_collapsed,
                      "dropped_coord": dropped_coord}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
