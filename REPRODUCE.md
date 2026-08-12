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

### Tier 1 — ships in this repository (57 MB)

| what | where | form |
|---|---|---|
| six evidence tables (284,744 boundaries) | `data/multicell/<CELL>_<RES>/annotation.tsv.gz` | gzip (17.5 MB; 145 MB raw) |
| candidate BED + six per-caller BEDs + tier BED | `data/multicell/<CELL>_<RES>/*_boundaries.bed.gz` | gzip (5.4 MB) |
| package manifests (accessions, counts, provenance) | `data/multicell/<CELL>_<RES>/manifest.json` | plain |
| HiCCUPS loop calls (the BD4 input), 3 cell lines | `data/loops/<CELL>_HiCCUPS_loops_hg38.bedpe.gz` | gzip (1.6 MB) |
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

**How they are obtained, stated precisely.** This repository ships **no graphical
application**. The released code is the library, the `creditad` CLI, and an optional
aiohttp API server (`pip install '.[server]'`). The Electron desktop client that wraps
this API is developed separately and **is not part of this release** — 0 of the 74
installed files is a frontend asset, and `creditad` is the only console script. What
follows describes the API server, which *is* released, and every statement below was
measured against a live server started from a clean clone of this repository.

**The Data Manager is an HTTP API, and it is reachable without any frontend.** With the
server running (see "Running the API server" below):

```bash
# which inputs the tool can see, and where it expects the ones it cannot
curl -s http://127.0.0.1:5001/api/data/status  | python -m json.tool   # HTTP 200
# the optional full-genome catalog: accession, size, destination folder
curl -s http://127.0.0.1:5001/api/data/catalog | python -m json.tool   # HTTP 200

# re-hash a shipped chr7 asset against its stored MD5
# (requires the one-time links at the end of this section, which is where the
#  server looks for these files; on a bare clone it reports File not found)
curl -s -X POST http://127.0.0.1:5001/api/data/verify \
     -H 'Content-Type: application/json' \
     -d '{"id":"gm12878_ctcf_chr7_bigwig"}'
# -> {"verified": true, "expected": "5758946f…", "actual": "5758946f…", "reason": "MD5 match"}

# asking it to FETCH a file it has no URL for
curl -s -X POST http://127.0.0.1:5001/api/data/download \
     -H 'Content-Type: application/json' -d '{"id":"gm12878_chr7_mcool"}'
# -> HTTP 400 {"error": "No configured download URL for GM12878_chr7_hg38.mcool;
#              place it manually in <dir>"}
```

Two levels of automation, and the difference matters:

- **Shipped chr7 demo assets (nine catalog entries)** — **verified, not downloaded.** These
  entries carry a stored MD5 but no URL. `POST /api/data/verify` re-hashes the file on disk
  and returns `{"verified": true, …, "reason": "MD5 match"}`; for a file that is absent it
  returns `{"verified": false, "reason": "File not found"}` — HTTP 200 either way, so check
  the body, not the status. `POST /api/data/download` on one of them returns **HTTP 400**
  naming the directory to place it in. Six of the nine (the chr7 CTCF and RAD21 bigWigs)
  ship in this repository under `data/chr7_preview_tracks/`; the three chr7 `.mcool` contact
  maps do not (42–91 MB each) and are in the Zenodo archive.
- **Optional full-genome ENCODE/4DN tracks (the table above)** — **listed, never fetched.**
  `GET /api/data/catalog` returns each accession with its size, the exact folder to place it
  in, and the source URL (`https://www.encodeproject.org/files/<ACC>/@@download/<ACC>.bigWig`
  for ENCODE, the 4DN file page for matrices). Once a file is in place the entry reports it
  as present. These entries carry **no stored MD5**, so the server does not and cannot
  checksum them; ENCODE publishes an MD5 on each file's own page.

In short: **the server verifies what ships and points at what does not. It downloads
nothing.** The retrieval step is yours, whether you use the API or not.

**The CLI has no download path at all.** `creditad annotate` never references the Data
Manager, and `creditad --help` advertises only `annotate` and `criteria`; in a core
`pip install creditad`, importing the downloader raises a message telling you to install
the server extra. **A CLI-only user retrieves the six bigWigs from the ENCODE URLs above
themselves** — which is the supported path for reproducing the tables, and the one
`scripts/reproduce_packages.py` expects.

### Running the API server

The server is usable headless. It needs a writable SQLite path, because the built-in
default is a *relative* path (`extra_mode/tad_tokens.db`) that does not exist in this
repository — without `TAD_DB_PATH` it exits during startup with
`sqlite3.OperationalError: unable to open database file`:

```bash
pip install '.[server]'
mkdir -p ~/.creditad
cd backend
TAD_DB_PATH="sqlite:///$HOME/.creditad/tad_tokens.db" CREDITAD_PORT=5001 python main.py
# -> [STARTUP|ready|done|All systems operational on port 5001
#    CrediTAD API ready [Port: 5001]
```

What you get is a **JSON API, not a web page**: `GET /` returns HTTP 404, because no
frontend is served. The routes are `/api/data/*` (above) and `/api/tadvci/*`
(`cells`, `multicell_datasets`, `annotate`, `card/{token}/{chrom}/{pos}`, `triage/{token}`,
`signout`, `history/{boundary_id}`, `records`).

Two of those route groups resolve their inputs from a **development-tree layout**
(`example_data/multicell`, `backend/extra_mode/chr7_tracks`) rather than this repository's
`data/` layout, and there is no environment variable to redirect them. To let the server
see the shipped packages and preview tracks, link them once from the repository root:

```bash
mkdir -p example_data backend/extra_mode
ln -s ../data/multicell            example_data/multicell
ln -s ../data/chr7_preview_tracks  example_data/chr7_tracks
ln -s ../../data/chr7_preview_tracks backend/extra_mode/chr7_tracks
ln -s ../../data/loops               backend/extra_mode/loops
```

Measured effect of those links: `GET /api/tadvci/multicell_datasets` goes from
`{"ok": false, "datasets": []}` to all **six** packages with their candidate counts
(28,259 / 69,479 / 28,476 / 65,157 / 29,513 / 63,860), and `GET /api/data/status` goes from
0 to **6 of 9** entries present — the three missing ones are the chr7 `.mcool` maps that are
not in this repository. This is a packaging wart, disclosed rather than papered over: the
reproduction path (`scripts/reproduce_packages.py`, environment-variable driven) does not
need it, and the CLI does not need it.

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
pytest                     # 239 passed, 39 skipped, 7 deselected
```

The 7 deselected tests are quarantined legacy experiments (`legacy_bcp`,
`legacy_reference` markers in `pyproject.toml`); they are not part of the production
contract. The 39 skips are tests whose data is not present in a clean checkout (archived
hg19 fixtures, host ChIP tracks); they report the missing file and skip rather than fail.

## Environment variables

| variable | meaning | required for |
|---|---|---|
| `CREDITAD_PACKAGES_DIR` | delivered packages | reproduction |
| `CREDITAD_DATA_ROOT` | ENCODE bigWigs (and 4DN mcools) | reproduction |
| `CREDITAD_LOOPS_DIR` | HiCCUPS `.bedpe` files | reproduction |
| `CREDITAD_CALLS_ROOT` | upstream per-caller calls | rebuilding packages from scratch |
| `CREDITAD_SHARED_RAW` | optional large-track search root | two optional tests |

Tests that need data which is not present skip themselves; they do not fail.
