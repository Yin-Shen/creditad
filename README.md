# CrediTAD

**Structured per-boundary evidence records for post-calling review of TAD boundaries.**

CrediTAD is a *downstream* evidence-audit workflow. You give it candidate boundaries that
some caller already produced; it attaches to each one a checkable evidence record with
four named criteria, combines them with a deterministic and inspectable rule, and reports
what it could and could not measure.

CrediTAD **does not call TADs**, does not rank callers, and does not produce a consensus
or a truth label.

## The evidence record

| Criterion | Evidence | Grades |
|---|---|---|
| **BD1** | cross-caller support among the supplied caller panel | strong / moderate / weak |
| **BD2** | CTCF ChIP enrichment | strong / moderate / weak / none / not_assessable |
| **BD3** | RAD21 (cohesin) ChIP enrichment | strong / moderate / weak / none / not_assessable |
| **BD4** | chromatin-loop-anchor overlap | strong / moderate / none / not_assessable |

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
pip install '.[app]'     # + the optional HTTP/GUI surface
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
method cannot express cross-caller agreement.

Print the rules actually compiled into your installed copy:

```bash
creditad criteria
```

## Repository layout

```
data/multicell/             six evidence tables + candidate/caller BEDs (gzipped)
data/loops/                 HiCCUPS loop calls (BD4 input)
data/chr7_preview_tracks/   chr7 display tracks -- NOT analysis inputs
backend/tad_vci/            evidence engine, tier rule, CLI
backend/builtin_callers_v2/ optional pure-Python detectors (development only)
backend/mcool_bed.py        bigWig / loop-call readers used by the annotation path
backend/validate_bigwig.py  genome-build guard for supplied tracks
backend/tests/              pytest suite
scripts/                    package regeneration and reproduction scripts
```

## Data

The six evidence tables (284,744 boundaries), the candidate and per-caller BEDs, the
HiCCUPS loop calls and chr7 preview tracks ship here under `data/` (~56 MB, gzipped where
it matters). The genome-wide ENCODE ChIP tracks and 4DN Hi-C matrices are obtained by
accession: the desktop application's Data Manager lists each one with its size and
destination folder, links it for one-click retrieval from ENCODE/4DN, and detects it once
in place. A CLI-only install has no download path and fetches them from ENCODE directly. Large binaries and the release snapshot are
archived on Zenodo.

Full three-tier breakdown, including exactly what is and is not checksum-verified:
[DATA_POLICY.md](DATA_POLICY.md). Accessions and the reproduction procedure:
[REPRODUCE.md](REPRODUCE.md).

## Reproducibility

`scripts/reproduce_packages.py` regenerates the six delivered evidence tables and compares
every column against the shipped ones. Verified 2026-08-12 across all six packages
(284,744 boundaries): **38/38 columns identical, tier agreement 1.000000**, with repeat
runs byte-identical. See [REPRODUCE.md](REPRODUCE.md).

## Citation

See [CITATION.cff](CITATION.cff).

## Licence

MIT — see [LICENSE](LICENSE).
