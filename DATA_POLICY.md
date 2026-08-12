# Data policy

Three tiers: what ships here, what the user obtains from ENCODE/4DN, what is archived.
Every size below was measured from the files themselves.

## Tier 1 — ships in this repository (57 MB, 155 files)

| content | path | size |
|---|---|---|
| six evidence tables, 284,744 boundaries | `data/multicell/<CELL>_<RES>/annotation.tsv.gz` | 17.5 MB gz (144.9 MB raw) |
| candidate, per-caller and tier BEDs | `data/multicell/<CELL>_<RES>/*_boundaries.bed.gz` | 5.4 MB gz (24.1 MB raw) |
| package manifests | `data/multicell/<CELL>_<RES>/manifest.json` | 10 kB |
| HiCCUPS loop calls (BD4 input), 3 cell lines | `data/loops/*.bedpe.gz` | 1.6 MB gz (4.9 MB raw) |
| chr7 **preview** ChIP tracks | `data/chr7_preview_tracks/*.bigWig` | 31.6 MB |
| source, tests, scripts, metadata, reproduction report | `backend/`, `scripts/`, `docs/`, root | 1.2 MB |

Tables and BEDs ship gzipped. Raw they are 169 MB, which is not appropriate for git;
gzip is lossless and `pandas.read_csv` reads `.gz` directly, so no downstream code changes.
Largest single file in the repository: 5.6 MB — comfortably under GitHub's thresholds.

**The chr7 preview tracks are display assets, not analysis inputs**, and are labelled as
such in `data/chr7_preview_tracks/README.md`. They are 250 bp re-bins and reproduce only
28/38 columns of the delivered tables (tier agreement 0.903). The genome-wide ENCODE
tracks reproduce all six packages exactly (38/38, tier agreement 1.000000).

## Tier 2 — obtained from ENCODE / 4DN by accession (~2.3 GB + Hi-C)

Not redistributed here: third-party files under the source portals' terms, and large.

**ENCODE ChIP-seq, GRCh38, fold change over control** — `https://www.encodeproject.org/files/<ACC>/`

| cell line | CTCF | RAD21 |
|---|---|---|
| GM12878 | ENCFF734CUT | ENCFF571ZJJ |
| IMR90 | ENCFF105FHL | ENCFF048PZI |
| HepG2 | ENCFF357NFO | ENCFF972ODZ |

**4DN Hi-C matrices** — `https://data.4dnucleome.org/files-processed/<ACC>/`:
GM12878 `4DNFIXP4QG5B`, IMR90 `4DNFIJTOIGOI`, HepG2 `4DNFIS6HAUPP`. Needed only to re-run
the upstream callers.

### What the released code does for you — and what it is not

**This repository ships no graphical application.** The released code is the library, the
`creditad` CLI, and an optional aiohttp **API server** (`pip install '.[server]'`; the old
`[app]` name still resolves as a deprecated alias). Measured on the installed package: 0 of
74 files is a frontend asset, `creditad` is the only console script, and the server returns
**HTTP 404 at `/`** because no UI is served. The Electron desktop client that consumes this
API is developed separately and is not part of this release.

The Data Manager is therefore an HTTP surface, not a page — and it is fully reachable
without a frontend. `REPRODUCE.md` gives the exact `curl` calls and the startup command
(the server needs `TAD_DB_PATH` set; its built-in default is a relative path that does not
exist here). What that surface does:

| | shipped chr7 demo assets (9 entries) | optional full-genome tracks (ENCODE/4DN) |
|---|---|---|
| listed with accession, size, destination folder | yes (`/api/data/status`) | yes (`/api/data/catalog`) |
| retrieved over the network by the server | **no** — 6 of the 9 ship here; `download` returns HTTP 400 naming the folder | **no** — the catalog gives the source URL; retrieval is yours |
| MD5-verified by the server | **yes** — stored digest, file re-hashed on request, returns `MD5 match` | **no** — these entries carry no stored MD5 |
| detected once present on disk | yes | yes |

The nine demo entries carry a stored MD5 but no URL. `POST /api/data/verify` re-hashes the
file on disk against that digest and returns `{"verified": true, …, "reason": "MD5 match"}`;
for an absent file it returns `{"verified": false, "reason": "File not found"}` — HTTP 200 in
both cases, so read the body. `POST /api/data/download` returns **HTTP 400** with the
directory to place the file in. Six of the nine (the chr7 CTCF and RAD21 bigWigs) are in
`data/chr7_preview_tracks/`; the three chr7 `.mcool` maps are not (42–91 MB each) and are in
the Zenodo archive.

For the full-genome tracks the catalog gives each accession, its size, the destination
folder and the source URL
(`https://www.encodeproject.org/files/<ACC>/@@download/<ACC>.bigWig`, or the 4DN file page).
The server does not stream them and **cannot vouch for their integrity**, because those
catalog entries carry no MD5. ENCODE publishes an MD5 on each file's page.

In one sentence: **the server verifies what ships and points at what does not; it downloads
nothing, and it is an API rather than an application.**

**The CLI has no download path.** `creditad --help` offers `annotate` and `criteria`; the
CLI never references the Data Manager, and importing it from a core install raises a message
pointing at the server extra. CLI-only users fetch the six bigWigs from the ENCODE URLs
themselves — the supported reproduction path, and the one `scripts/reproduce_packages.py`
expects.

Every statement in this section was established by driving a live server started from a
clean clone of this repository — not by reading the code.

## Tier 3 — archived to Zenodo

- a snapshot of this repository at the tagged release;
- the three chr7 `.mcool` contact maps (42–91 MB each — above GitHub's thresholds,
  which is why they are not in tier 1);
- the raw per-caller call sets the packages were built from.

The DOI is minted at release and recorded in `CITATION.cff` and `.zenodo.json`.

## Provenance and licence

Boundary BEDs were produced by six external callers run locally; parameter choices follow
published defaults where obtainable. This is a local panel, not a byte-reproduction of any
published pipeline. ENCODE and 4DN files remain under their originating portals' terms;
the code in this repository is MIT.
