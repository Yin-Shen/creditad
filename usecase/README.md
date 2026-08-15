# Pre-registered reading demonstration (ClinVar)

This directory holds the reproduction material for the reading demonstration reported in
the CrediTAD application note (Results R5 and Supplementary Note 19). It shows what a
CrediTAD evidence record adds over a bare coordinate list when a reader looks at one
locus — nothing more. **No accuracy, validation, calling or enrichment claim is made
here, and none can be derived from this material.**

## The discipline: rule first, query second

The selection rule was written, hashed and timestamped **before any locus-level query
touched the candidate tables**, so the example shown could not be chosen after seeing the
answer. The rule is [`PREREGISTERED_RULE.md`](PREREGISTERED_RULE.md); the seal — hash,
byte count, UTC timestamp and append-only amendment history — is
[`PREREG_SEAL.json`](PREREG_SEAL.json).

Check the seal yourself:

```bash
sha256sum PREREGISTERED_RULE.md
# 1be02799aeab779f3e5da9ec038916ca54ede3c5c18d72d7b1f89b54969d58e9
```

That digest is the third entry in the seal's `amendment_history`. The file was amended
twice after sealing, **append-only** — no text above section 7 was ever edited, and each
amendment carries its own dated hash:

| UTC | event | sha256 |
|---|---|---|
| 2026-08-15T02:38:00Z | sealed (12,390 bytes) | `71e963c8…20c0cc0` |
| 2026-08-15T02:51:23Z | amendment 1 — ClinVar inner-bound observation | `05b2f6be…3135978` |
| 2026-08-15T05:18:05Z | amendment 2 — header claim corrected (19,036 bytes) | `1be02799…69d58e9` |

**What the seal does and does not cover.** Amendment 2 exists because the original header
overstated its own scope. The seal predates every locus-level query and predates the
ClinVar download, but it does **not** predate aggregate inspection of the delivered
tables: before sealing, the 38 column names, the six per-package row counts, the `chrom`
domain, the constant window parameters, one package manifest and the GM12878 10 kb
`D_tier` value counts had been read. No pre-seal command returned, filtered or ordered a
single row, so no candidate or variant identity existed when the rule was written. The
consequence, stated plainly: the seal removes post-hoc locus selection and post-hoc
rule-writing, but not knowledge of the corpus marginals. Section 7 of the rule file lists
the pre-seal commands verbatim.

## The source data, and how to get it back

The disease loci come from the ClinVar `variant_summary.txt.gz` **full release**,
retrieved from
`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz` on
**2026-08-15**, GRCh38 coordinates used as released with no liftover.

NCBI does not embed a version string in this file, so the release is identified by the two
values recorded in [`SOURCE_PROVENANCE.json`](SOURCE_PROVENANCE.json):

| field | value |
|---|---|
| MD5 (matches the value NCBI publishes beside the file) | `897c38fd97dff56f83f1c7a674eb35f3` |
| SHA-256 of the retrieved file | `44d47173c6f786a442df0503e2758918e26e9d6c276e4f15f40c71d633596e4b` |
| server `Last-Modified` | `Mon, 10 Aug 2026 12:59:47 GMT` |
| size | 441,836,319 bytes |

ClinVar replaces this file in place. A later download **will** differ; verify the MD5
before assuming you have reproduced the same release.

## What ships here, and what does not

| file | what it is |
|---|---|
| `PREREGISTERED_RULE.md` | the sealed selection rule, with append-only amendments |
| `PREREG_SEAL.json` | hash, timestamps, amendment history, pre-seal inspection list |
| `SOURCE_PROVENANCE.json` | ClinVar release identity, retrieval time, checksums |
| `summary_counts.json` | every count reported in the paper, including the locus funnel |
| `contrast_table.csv` | the coordinate-list-vs-evidence-record contrast, field by field |
| `CONTRAST.md` | the same contrast in prose |
| `selection_record.json` | which record the rule selected, and the two evidence cards |
| `verification.json` | 68 independent re-derivations of the reported values (68/68 match) |
| `locus_funnel.json` | the ClinVar filter funnel, stage by stage |
| `tier_context_hits_vs_corpus.csv` | tier distribution of hits against the whole corpus |
| `data/locus_set.tsv` | the 4,229 filtered loci — the checksummed output of step 1 |
| `scripts/` | the four scripts, in run order |

**`hits_all.csv` (29,762 rows, 23.4 MB) is deliberately not committed here.** It exceeds
the repository's file-size policy. It is regenerated exactly by step 2 below, and it is
also provided as Supplementary Data with the article.

`data/locus_set.tsv` **is** committed (1.0 MB, sha256
`d6851a6711eace36d355308267c3baff1180fe6d9bfc86a3d17716980c3cfd7f`) so that the
verification step runs from a clone without the 441 MB ClinVar download. Step 1
regenerates it from the release identified above.

## Re-running it

The six evidence packages ship in this repository, so steps 3 and 4 need no external
input. Paths default to this repository and are overridable:

| variable | default |
|---|---|
| `CREDITAD_PACKAGES_DIR` | `<repo>/data/multicell` |
| `CREDITAD_BACKEND_DIR` | `<repo>/backend` |

The scripts read the delivered tables either gzipped (as they ship here) or plain.

```bash
# 1. filter ClinVar to the pre-registered locus set  -> data/locus_set.tsv
#    (needs the 441 MB release in data/; the output already ships here)
python scripts/01_build_locus_set.py

# 2. intersect the locus breakpoints against all six packages -> hits_all.csv (23.4 MB)
python scripts/02_intersect.py

# 3. apply the sealed rule; emit the selection and the contrast
python scripts/03_select_and_contrast.py

# 4. re-derive every reported number independently  -> verification.json
python scripts/04_verify.py
```

Step 4 is the check that matters, and it is the one to run first: it needs only
`data/locus_set.tsv`, `selection_record.json` and the shipped packages — not the ClinVar
download and not `hits_all.csv`. It recomputes each reported value from the delivered
tables using the engine's own functions and compares against what is published here.
**68 checks, 68 matching, zero failures**, re-confirmed from this tree on 2026-08-15.

## The numbers this produces

| quantity | value |
|---|---|
| ClinVar rows read | 9,036,351 |
| after the pre-registered filter | 4,381 loci |
| tested (after collapsing 152 duplicate tuples) | 4,229 loci / 8,458 breakpoints |
| hits across all six packages | 29,762 |
| distinct loci with a hit | 4,044 |
| distinct candidate boundaries with a hit | 13,078 |
| corpus rows searched | 284,744 |

The two records read in the article are ClinVar VariationID **541215** (RCV000651427), the
primary example at chr1:1,310,000 (tier D5), and VariationID **58037** (RCV000051780), the
mandatory lowest-tier companion at chr1:830,000 (tier D1). The companion is fixed by the
sealed rule, not chosen afterwards: it exists so the demonstration cannot show only the
flattering end of the tier distribution.
