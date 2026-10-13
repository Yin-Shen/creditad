# CrediTAD

**Structured per-boundary evidence records for post-calling review of TAD boundaries.**

CrediTAD is a *downstream* evidence-audit workflow. You give it candidate boundaries that
some caller already produced; it attaches to each one a checkable evidence record with
four named criteria, combines them with a deterministic and inspectable rule, and reports
what it could and could not measure.

CrediTAD **does not call TADs**, does not rank callers, and does not produce a consensus
or a truth label.

## The evidence record

| Criterion | Evidence | Output column | Grades |
|---|---|---|---|
| **BD1** | cross-caller support among the supplied caller panel | `votes` -> `BD1_caller_support` | strong / moderate / weak |
| **BD2** | CTCF ChIP enrichment | `ctcf`, `ctcf_fold` -> `BD2_ctcf` | strong / moderate / weak / none / not_assessable |
| **BD3** | RAD21 (cohesin) ChIP enrichment | `rad21`, `rad21_fold` -> `BD3_rad21` | strong / moderate / weak / none / not_assessable |
| **BD4** | chromatin-loop-anchor overlap | `loop_anchor_count` -> `BD4_loop_anchor` | strong / moderate / none / not_assessable |

In the column list, `measured -> graded`: the left name holds the measured value, the right name holds the grade.

**Paper-to-software correspondence for BD4.** The manuscript calls this criterion
**"loop-anchor engagement (number of loop ends in the window)"**. The software column
is **`loop_anchor_count`** and the `creditad criteria` table prints it as
**"loop-anchor overlap"**. All three name the same quantity: the number of loop ends
falling within the window. The delivered column name is deliberately **not** renamed,
so records already produced stay readable; this table is the mapping.

BD4 grades depend on the density of the loop file supplied, and HiCCUPS call counts
differ several-fold between cell lines, so absolute BD4 rates are not comparable
across cell lines.

These combine into a **D1–D5 evidence tier** by the `max_support_v2` rule, together with a
human-readable derivation of the branch that fired.

Two properties are the point of the tool:

1. **An unmeasured axis is not a negative one.** `not_assessable` (nobody supplied the
   assay) is kept distinct from `none` (the assay was supplied and the position was not
   enriched), in the record *and* in the tier rule.
2. **The aggregation is inspectable.** The rule is a short deterministic function of the
   four grades, printed per boundary, not a fitted model.

### What the tier is not

The tier is **not** a probability, **not** a ground-truth label, and **not** a measure of
caller quality. BD2–BD4 are correlated readouts of the same CTCF/cohesin biology and all
three feed the tier, so enrichment of those signals cannot be used to validate it.

CTCF/cohesin looping does not account for all human TAD boundaries. Absence of CTCF,
cohesin or loop-anchor signal at a position is an absence of *these measured readouts* —
never evidence that the boundary is false.

## Install

```bash
pip install .            # library + `creditad` CLI
pip install '.[test]'    # + test suite
pip install '.[server]'  # + the optional HTTP/JSON API server (no GUI; see REPRODUCE.md)
```

Python >= 3.10.

## Use

```bash
creditad annotate \
    --boundaries candidates.bed \
    --panel TopDom=TopDom.bed --panel OnTAD=OnTAD.bed --panel DI=DI.bed \
    --ctcf CTCF.bigWig --rad21 RAD21.bigWig --loops loops.bedpe \
    --cell-line GM12878 --assembly hg38 --resolution 25000 \
    --out results/
```

Omit `--ctcf`, `--rad21` or `--loops` and that axis is recorded as `not_assessable`
rather than as absent signal. With no `--panel`, BD1 is floored to `weak`: a single
method cannot express cross-caller agreement. BED and loop files may be plain or
gzip-compressed; a loop file that yields no loop anchors stops the run with an error
(release v1.0.3 and earlier read a gzipped loop file as empty without reporting it). A
BEDPE loop file is read from its first six columns (chr1 x1 x2 chr2 y1 y2), and lines
starting with `#` are skipped, so ENCODE / Juicer loop files whose header line starts with
`#chr1` are read as distributed (release v1.0.4 and earlier rejected them).

Quick start, from files that ship in this repository, checked byte for byte against
the expected output:

```bash
python examples/quickstart/run_quickstart.py      # prints PASS and exits 0
```

`examples/reading_records.md` shows how to read the table for a follow-up decision
(caller support against measured molecular support). A container image can be built
from the `Dockerfile`:

```bash
docker build -t creditad:1.0.5 .
docker run --rm --entrypoint python creditad:1.0.5 examples/quickstart/run_quickstart.py
```

Print the rules actually compiled into your installed copy:

```bash
creditad criteria
```

## Repository layout

```
data/multicell/             six evidence tables + candidate/caller BEDs (gzipped), v1.0
data/multicell_v1.1/        the same tables with gap/blacklist positions not_assessable
data/loops/                 HiCCUPS loop calls (BD4 input)
data/chr7_preview_tracks/   chr7 display tracks -- NOT analysis inputs
backend/tad_vci/            evidence engine, tier rule, CLI
backend/builtin_callers_v2/ optional pure-Python detectors (development only)
backend/mcool_bed.py        bigWig / loop-call readers used by the annotation path
backend/validate_bigwig.py  genome-build guard for supplied tracks
backend/tests/              pytest suite
scripts/                    reproduce_packages.py (regenerate the six tables and compare
                            every column) and mask_and_retier.py (derive the v1.1 tables)
```

## Data

The six evidence tables (284,744 boundaries), the candidate and per-caller BEDs, the
HiCCUPS loop calls and chr7 preview tracks ship here under `data/` (56 MB, gzipped
where it matters). The genome-wide ENCODE ChIP tracks and 4DN Hi-C matrices are obtained by
accession. **This repository ships no graphical application source.** The code here is the
library, the `creditad` CLI and an optional aiohttp API server, and a `pip` install gets
no frontend. The Electron desktop client is developed in a separate repository,
https://github.com/Yin-Shen/creditad-desktop, the source of the desktop builds; its
pre-built binaries for Linux and Windows are attached to the GitHub release, each with a
published SHA-256 (see [Desktop application](#desktop-application)). The server's Data Manager
endpoints (`/api/data/status`, `/api/data/catalog`, `/api/data/verify`) report which inputs
are present, re-hash the shipped chr7 assets against stored MD5s, and give the accession,
size and destination folder for each track you still need -- they never fetch the
full-genome ENCODE files. A CLI-only install has no download path at all; you retrieve the
tracks from ENCODE yourself.

**This repository is the single home of the code and of the data that ships with it.**
There is no separate archive deposit and no DOI. Two classes of file are deliberately not
in this tree: the three chr7 `.mcool` contact maps (42.1, 48.3 and 90.8 MB), which exceed
the repository's file-size policy and ship instead inside the desktop builds attached to
the [v1.0.5 release](https://github.com/Yin-Shen/creditad/releases/tag/v1.0.5) under
`resources/example_data/chr7_tracks/`; and the raw per-caller call sets, which are not
redistributed — each package manifest records the call set it was built from by path and
SHA-256, and the per-caller boundary BEDs derived from it ship here under
`data/multicell/<CELL>_<RES>/`, one file per caller, with row counts matching the
manifest's per-caller counts exactly.

Full three-tier breakdown, including exactly what is and is not checksum-verified:
[DATA_POLICY.md](DATA_POLICY.md). Accessions and the reproduction procedure:
[REPRODUCE.md](REPRODUCE.md).

## Reproducibility

`scripts/reproduce_packages.py` regenerates the six delivered evidence tables and compares
every column against the shipped ones. Verified 2026-08-12, again on 2026-10-12 with
release v1.0.4 and on 2026-10-13 with release v1.0.5, across all six packages (284,744
boundaries): **38/38 columns identical, tier agreement 1.000000**, with repeat runs
byte-identical. See [REPRODUCE.md](REPRODUCE.md).

## Availability

This repository, at tagged release **v1.0.5**, is the single home of the code and of the
data that ships with it; there is no separate archive deposit and no DOI. It contains the
package, the test suite, the six evidence tables, the candidate and per-caller BEDs, the
package manifests, the HiCCUPS loop calls and `scripts/reproduce_packages.py`; release v1.1 of
the evidence tables, with gap and blacklist positions recorded as `not_assessable`, ships under
`data/multicell_v1.1/` with `scripts/mask_and_retier.py` and its impact tables. Install
with `pip install .` for the library and command-line interface, or `pip install
'.[server]'` to add the local HTTP/JSON API service. Desktop application builds — an
AppImage for Linux (client 1.0.2), and an installer and a portable archive for Windows
(client 1.0.1) — are attached to the v1.0.5 release; their SHA-256 digests are listed in
`RELEASE_ASSETS.sha256` (see [Desktop application](#desktop-application)). Python >= 3.10.

The ChIP-seq tracks the evidence tables were built from are ENCODE files (GRCh38, fold
change over control), retrieved by accession from
`https://www.encodeproject.org/files/<ACCESSION>/`: GM12878 CTCF `ENCFF734CUT` and RAD21
`ENCFF571ZJJ`; IMR90 CTCF `ENCFF105FHL` and RAD21 `ENCFF048PZI`; HepG2 CTCF `ENCFF357NFO`
and RAD21 `ENCFF972ODZ`. The Hi-C matrices the callers were run on are 4DN files
`4DNFIXP4QG5B` (GM12878), `4DNFIJTOIGOI` (IMR90) and `4DNFIS6HAUPP` (HepG2), needed only
to re-run the upstream callers.

## Citation

See [CITATION.cff](CITATION.cff).

## Licence

MIT — see [LICENSE](LICENSE).

## Desktop application

Pre-built desktop binaries are attached to the [v1.0.5
release](https://github.com/Yin-Shen/creditad/releases/tag/v1.0.5): a Linux AppImage, a
Windows NSIS installer and a Windows portable zip. Their SHA-256 digests are listed in
`RELEASE_ASSETS.sha256`:

| file | client | embedded backend | size | SHA-256 |
|---|---|---|---|---|
| `CrediTAD-Linux.AppImage` | 1.0.2 | v1.0.5 | 540,753,929 B | `81061aef1133d27e758dd86cf23b0fe3b89dfc2cbe6a8bac0312542896bea51e` |
| `CrediTAD-Windows-Setup.exe` | 1.0.1 | v1.0.1 | 452,125,658 B | `a7b162421fc6171f27a52dad658e437215e98df3d037daa2d8ac8e1726688d81` |
| `CrediTAD-Windows-1.0.1.zip` | 1.0.1 | v1.0.1 | 512,613,927 B | `9f5f2c2649a515d9a116ef251fd990e7a9c65c48c9766c7007cd71c247c3c57c` |

The client source of each build is the tag of the same version in the desktop client's
repository, https://github.com/Yin-Shen/creditad-desktop. Each build carries its own Python
runtime, and no separate Python installation is needed. Verify a download against the
SHA-256 in `RELEASE_ASSETS.sha256` before running it.

**Linux AppImage (client 1.0.2).** It starts the backend of release v1.0.5 that it bundles;
its runtime includes `cooler` and `pyBigWig`. Mounted read-only, as the AppImage runtime
mounts it, the bundled service answered its health check, all six demo packages loaded with
their caller, ChIP and contact-map tracks, and the client's three end-to-end scripts
(CLI/API agreement, demo packages, Import page) passed. The AppImage needs FUSE 2
(`libfuse.so.2`); where unprivileged user namespaces are restricted it must be started with
`--no-sandbox`. The v1.0.1 Linux AppImage, still attached to the v1.0.1 and v1.0.2 releases
(SHA-256 `a991dc39b7477b10d28a954e513cd6b7c54e6bdea5ae8c81d000bda604120f4e`), does not
start the backend it bundles (its launcher looks for `runtime/python`, while the runtime
has `runtime/bin/python` and lacks `cooler` and `pyBigWig`); use the 1.0.2 build instead.

**Windows builds (client 1.0.1).** They are the files attached to v1.0.1, carried forward
unchanged; the portable archive keeps its build-time filename. They ship
`runtime/python.exe`, the interpreter their launcher looks for, together with `cooler` and
`winbbi` (the Windows bigWig reader), and were cross-built on Linux and confirmed working
on Windows by the authors; no Windows version matrix was exercised. They embed the backend
of release v1.0.1, which reads a gzip-compressed loop file as empty without reporting it
(every boundary is then graded BD4 `none`): decompress loop files before loading them in
the Windows builds.
