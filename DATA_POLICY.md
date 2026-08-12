# Data policy

Three tiers: what ships here, what the user obtains from ENCODE/4DN, what is archived.
Every size below was measured from the files themselves.

## Tier 1 — ships in this repository (~56 MB)

| content | path | size |
|---|---|---|
| six evidence tables, 284,744 boundaries | `data/multicell/<CELL>_<RES>/annotation.tsv.gz` | 17.5 MB gz (144.9 MB raw) |
| candidate, per-caller and tier BEDs | `data/multicell/<CELL>_<RES>/*_boundaries.bed.gz` | 5.4 MB gz (24.1 MB raw) |
| package manifests | `data/multicell/<CELL>_<RES>/manifest.json` | 12 kB |
| HiCCUPS loop calls (BD4 input) | `data/loops/*.bedpe.gz` | 1.7 MB gz (5.2 MB raw) |
| chr7 **preview** ChIP tracks | `data/chr7_preview_tracks/*.bigWig` | 31.6 MB |
| source, tests, scripts, metadata | `backend/`, `tests/`, `scripts/`, root | 1.2 MB |

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

### What the application does for you

The desktop application (`pip install 'creditad[app]'`) ships a Data Manager page that
resolves every input the tool needs. It operates at two levels, and the difference matters:

| | shipped chr7 demo assets (9 files) | optional full-genome tracks (ENCODE/4DN) |
|---|---|---|
| listed with accession, size, destination folder | yes | yes |
| retrieved over the network by the application | no — they ship with the package, or are linked from a local copy | no — one-click link to the source repository |
| MD5-verified by the application | **yes** (stored digest; the file on disk is re-hashed on request) | **no** (no stored MD5 for these entries) |
| detected once present on disk | yes | yes — the row switches to "In use" |

The nine demo entries carry a stored MD5 but no URL: the application resolves them from
the installed package or a local shared copy and re-hashes the file on disk against that
digest, reporting `MD5 match`. Requesting one that is absent returns "No configured
download URL"; the application does not fetch it. The download code path does hash-then-
refuse on mismatch, but no shipped entry reaches it.

For the full-genome tracks each ENCODE item renders a Download button pointing at
`https://www.encodeproject.org/files/<ACC>/@@download/<ACC>.bigWig` and each 4DN matrix a
"4DN page" link, alongside the expected size and the `place_under` folder. The user is
therefore never left to hunt for files — but the application is not streaming those
particular files itself, and it cannot vouch for their integrity, because the catalog
entries for them carry no MD5. ENCODE publishes an MD5 on each file's page.

**The CLI has no download path.** `creditad --help` offers `annotate` and `criteria`; the
CLI never references the Data Manager, and importing it from a core install raises a
message pointing at `creditad[app]`. CLI-only users fetch the six bigWigs from the ENCODE
URLs themselves.

Each statement in this section was established by testing the shipped code — driving the
HTTP handlers and probing a core-only install — not by reading it.

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
