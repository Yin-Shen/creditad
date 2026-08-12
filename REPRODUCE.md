# Reproducing the CrediTAD results

This repository contains the tool and the scripts. It does **not** contain the data:
the raw ChIP and Hi-C inputs are ENCODE/4DN files that we are not entitled to
redistribute, and they are large. This document says exactly what to fetch, from where,
and what you can then reproduce.

## What is reproducible, and to what standard

`scripts/reproduce_packages.py` regenerates the **six delivered evidence tables**
(GM12878 / IMR90 / HepG2 x 25 kb / 10 kb; 284,744 candidate boundaries in total) from the
delivered candidate BEDs and the molecular inputs, and compares **every column** of the
result against the shipped `annotation.tsv`.

Verified on 2026-08-12, on all six packages: **38 of 38 columns identical, evidence-tier
agreement 1.000000, zero mismatched columns**, and three repeat runs of the same package
produced byte-identical output. The check is exact, not approximate.

What this does **not** do: it does not re-call TADs. The six per-caller boundary BEDs are
*inputs*, produced upstream by six external callers (TopDom, SpectralTAD, OnTAD, MSTD,
Arrowhead, DI). CrediTAD is a downstream evidence-audit workflow and never calls
boundaries itself.

## Data availability — three tiers

### Tier 1 — ships in this repository (~56 MB)

| what | where | form |
|---|---|---|
| six evidence tables (284,744 boundaries) | `data/multicell/<CELL>_<RES>/annotation.tsv.gz` | gzip (17.5 MB; 145 MB raw) |
| candidate BED + six per-caller BEDs + tier BED | `data/multicell/<CELL>_<RES>/*_boundaries.bed.gz` | gzip (5.4 MB) |
| package manifests (accessions, counts, provenance) | `data/multicell/<CELL>_<RES>/manifest.json` | plain |
| HiCCUPS loop calls (the BD4 input) | `data/loops/<CELL>_HiCCUPS_loops_hg38.bedpe.gz` | gzip (1.7 MB) |
| chr7 **preview** ChIP tracks | `data/chr7_preview_tracks/*.bigWig` | 31.6 MB |

`pandas.read_csv` reads the `.gz` files directly; nothing needs unpacking first.

**The chr7 preview tracks are for display, not for analysis.** They are 250 bp re-binned
subsets: measured median window-mean ratio 1.118 (CTCF) and 1.109 (RAD21), Pearson r 0.935
and 0.909, and **0 of 1,528** chr7 windows agree with the genome-wide track to within
1e-6. Substituting them reproduces only **28 of 38** columns of the delivered tables, at
evidence-tier agreement **0.903**. With the genome-wide ENCODE tracks the same command
reproduces **all six packages exactly: 38/38 columns, grade agreement 1.000000 on the
CTCF and RAD21 axes, tier agreement 1.000000**. Both facts are measured, and they are
why this document keeps repeating the distinction.

The matching chr7 `.mcool` contact maps (42–91 MB each) are **not** in the repository —
they exceed GitHub's file-size thresholds — and are in the Zenodo archive.

### Tier 2 — retrieved from ENCODE / 4DN by accession

| cell line | CTCF | RAD21 | Hi-C matrix (4DN) |
|---|---|---|---|
| GM12878 | `ENCFF734CUT` | `ENCFF571ZJJ` | `4DNFIXP4QG5B` |
| IMR90 | `ENCFF105FHL` | `ENCFF048PZI` | `4DNFIJTOIGOI` |
| HepG2 | `ENCFF357NFO` | `ENCFF972ODZ` | `4DNFIS6HAUPP` |

ENCODE files: `https://www.encodeproject.org/files/<ACCESSION>/` — GRCh38, output type
"fold change over control", ~2.3 GB for the six. 4DN matrices:
`https://data.4dnucleome.org/files-processed/<ACCESSION>/`; needed only if you want to
re-run the upstream callers, not to reproduce the evidence tables.

**How they are obtained, stated precisely.** The desktop application (the `app` extra)
ships a Data Manager that resolves every input the tool needs, at two different levels of
automation:

- **Shipped chr7 demo assets (nine files)** — **MD5-verified inside the application, not
  downloaded by it.** These nine entries carry a stored digest but no URL: they ship with
  the package or are resolved from a local copy, and the verify endpoint re-hashes the file
  on disk and reports `MD5 match`. Asking the application to fetch one that is absent
  returns "No configured download URL" and the path where you should place it. (The
  download code path does refuse a file whose MD5 disagrees, but no shipped entry reaches
  it, because none has a URL.)
- **Optional full-genome ENCODE/4DN tracks (the table above)** — listed per accession with
  the file size and the exact destination folder, and **linked for one-click retrieval
  from the source repository**: each ENCODE item renders a Download button pointing at
  `https://www.encodeproject.org/files/<ACC>/@@download/<ACC>.bigWig`, and each 4DN matrix
  renders a "4DN page" link. Once a file is in place the application detects it and the
  row switches to "In use"; until then it shows the path where it is expected.

So you are not hunting for six files: the application tells you which accession, how large
it is, where to put it, gives you the link, and confirms when it can see it. What it does
**not** do is stream those particular files itself — that step is a click through to
ENCODE — and because the full-genome catalog entries carry no stored MD5, **they are not
checksum-verified by the application**. ENCODE publishes an MD5 on each file's own page if
you want to check one. In-app MD5 verification covers the nine shipped demo assets and
nothing else. In short: the application verifies, it does not download.

**The CLI has no download path at all.** `creditad annotate` never references the Data
Manager, and `creditad --help` advertises only `annotate` and `criteria`; in a core
`pip install creditad`, importing the downloader raises a message telling you to install
`creditad[app]`. **A CLI-only user retrieves the six bigWigs from the ENCODE URLs above
themselves.** All three statements above were established by testing the shipped code, not
by reading it.

### Tier 3 — archived to Zenodo

A snapshot of this repository at the tagged release, plus the large files that cannot go
in git: the three chr7 `.mcool` contact maps and the raw per-caller call sets. The Zenodo
DOI is minted at release and recorded in `CITATION.cff` and `.zenodo.json`.

## Running it

```bash
python -m venv .venv && . .venv/bin/activate
pip install '.[test]'

# Tier-1 inputs are already in the repository:
export CREDITAD_PACKAGES_DIR=$PWD/data/multicell
export CREDITAD_LOOPS_DIR=$PWD/data/loops
# Tier-2 inputs you fetch (see the accession table above):
export CREDITAD_DATA_ROOT=/path/to/encode

python scripts/reproduce_packages.py                       # all six
python scripts/reproduce_packages.py --packages GM12878_25kb   # one
```

Exit code 0 iff every column of every checked package matches. If an input is missing the
script says which one and skips that package rather than silently comparing less.
Runtime: about 100 s for all six on one core (see the benchmark in the paper).

## Test suite

```bash
pip install '.[test]'
pytest                     # 290 passed, 26 deselected
```

The 26 deselected tests are quarantined legacy experiments (`legacy_bcp`,
`legacy_reference` markers in `pyproject.toml`); they are not part of the production
contract.

## Environment variables

| variable | meaning | required for |
|---|---|---|
| `CREDITAD_PACKAGES_DIR` | delivered packages | reproduction |
| `CREDITAD_DATA_ROOT` | ENCODE bigWigs (and 4DN mcools) | reproduction |
| `CREDITAD_LOOPS_DIR` | HiCCUPS `.bedpe` files | reproduction |
| `CREDITAD_CALLS_ROOT` | upstream per-caller calls | rebuilding packages from scratch |
| `CREDITAD_SHARED_RAW` | optional large-track search root | two optional tests |

Tests that need data which is not present skip themselves; they do not fail.
