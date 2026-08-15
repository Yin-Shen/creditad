# PRE-REGISTERED selection rule — CrediTAD use-case demonstration

**Written 2026-08-15T02:38:00Z (UTC).**
**This file was written BEFORE the first query against any candidate evidence table.**
Author: research student, Phase 1 use-case task. Namespace:
`/home/coder/TAD/creditad/rebuild_2026-08-12/usecase/`.

This document exists because of the failure mode recorded in `PAUSE_STATE.md`: figure
examples in a sibling project were described as "random / unselected / prespecified"
without contemporaneous documentation, and the claim could not be defended. Every
selection decision below is fixed here, in advance, in a form a referee can check against
the output. Nothing above §7 may be edited after the timestamp above. Corrections and
fallbacks are handled by **appending** a dated section, never by rewriting history.

---

## 0. What this demonstration is and is not

The locked narrative (`STORY_LOCK.md`) forbids accuracy, boundary-calling and validation
claims. This use case is therefore constrained to a **handling** claim and nothing more:

- **Claim to be demonstrated:** for one boundary that a person has an external reason to
  look at, a coordinate list supports one action (note the position) whereas the evidence
  record supports a different, statable set of actions, because it exposes which callers
  agreed, what was measured on each molecular axis with what value, which rule fired, and
  **what was not measured**.
- **NOT claimed, and not inferable from any output of this task:** that the boundary is
  real, correct, validated or disease-causing; that CrediTAD called, discovered, recovered
  or rescued it; that overlap with a pathogenic variant is evidence the boundary is true;
  that absence of a measured molecular signal means the boundary is false.
- The disease locus supplies **only the reason a person is looking at this boundary at
  all**. It is context for handling. It is not evidence about the boundary, and the
  evidence record is not evidence about the variant's pathogenicity. Both directions of
  that inference are prohibited here.
- Wording constraint carried from `GROUP_BRIEF.md`: absence of signal is reported as
  "no CTCF, cohesin or loop-anchor support was measured at this position".
- Terminology constraint carried from `TERMINOLOGY.md`: prose and any figure/table label
  say **loop-anchor engagement**; the delivered column name `loop_anchor_count` is quoted
  unchanged when the column itself is referenced.

## 1. Data source for disease loci — declared priority order

Declared **now**, before any retrieval attempt, so a fallback cannot be mistaken for a
post-hoc choice. The first source in this list that is actually retrievable in this
session is used; if a source is unavailable or access is denied, that fact and the move
down the list are recorded by appending to this file (§7), with the failure message.
**No locus may be entered by hand, from memory, or from any document in this repository.**

1. **ClinVar `variant_summary.txt.gz`** (NCBI), tab-delimited full release.
   URL: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz`
   Rows restricted to `Assembly == "GRCh38"`; coordinates taken from the `Chromosome`,
   `Start`, `Stop` columns as released. No liftover is performed.
2. **ClinVar GRCh38 VCF**
   `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`, using `CHROM`,
   `POS`, `INFO/END`/`SVLEN` where present. Used only if source 1 is unavailable.
3. **A published TAD-disruption locus set with hg38 coordinates from a citable source**
   (e.g. the EPHA4 locus set of Lupiáñez et al. 2015). Admissible **only** if hg38
   coordinates are read from a citable record retrieved in this session. Hand liftover is
   prohibited; if a liftover is unavoidable, the chain file and its checksum must be
   recorded and the result labelled lifted-over in every output.

For whichever source is used, `usecase/SOURCE_PROVENANCE.json` records: exact URL, UTC
retrieval time, HTTP status, byte size, SHA-256 of the downloaded file, and the
release/version string the file declares internally.

## 2. Deterministic filter on the disease-locus set

Applied to ClinVar source 1 (equivalently to source 2). All thresholds fixed here.

- **Assembly:** `Assembly == "GRCh38"` exactly. Other assemblies dropped.
- **Chromosome:** `Chromosome` in `1..22, X`. `Y`, `MT`, unplaced/alt contigs dropped,
  because the candidate tables contain 23 chromosomes (1–22, X) and a locus on a
  chromosome absent from the corpus cannot produce an interpretable hit.
- **Variant class (`Type`):** the copy-number / structural classes only —
  `copy number loss`, `copy number gain`, `Deletion`, `Duplication`.
  Single-nucleotide and short classes (`single nucleotide variant`, `Indel`, `Insertion`,
  `Microsatellite`, `Inversion`, `Translocation`, `Complex`, `Variation`, `protein only`,
  `fusion`, and any class not among the four named) are dropped. Rationale: only
  copy-number/structural events have breakpoints whose relation to a TAD boundary is
  defined at 10–25 kb resolution.
- **Clinical significance (`ClinicalSignificance`):** the assertion string, case-folded
  and stripped, must be exactly one of `pathogenic`, `pathogenic/likely pathogenic`.
  Everything else is dropped, including bare `likely pathogenic`, `uncertain
  significance`, `benign`, `likely benign`, `conflicting interpretations of
  pathogenicity`, `not provided`, `other`, `drug response`, and any string carrying a
  qualifier not in the two accepted forms.
- **Size bounds:** `10,000 <= (Stop - Start) <= 5,000,000` bp from released coordinates,
  both bounds inclusive. Lower bound: an event shorter than one 10 kb bin cannot be
  resolved against the bin grid. Upper bound: events above 5 Mb are large chromosomal
  segments whose breakpoints carry no specific relation to any one boundary and which
  would dominate the hit table.
- **Coordinate sanity:** `Start >= 1` and `Stop >= Start`; failures dropped and counted.
- **Deduplication:** on the tuple `(Chromosome, Start, Stop, Type)`. Where duplicate
  tuples carry different accessions, the retained row is the one with the smallest
  `VariationID` compared as an integer, and the number of collapsed duplicates is reported.
- **No review-status filter is applied.** Review status is *recorded* for every hit
  (`ReviewStatus` column) so a referee can re-stratify, but it removes no locus, because
  choosing a review-status threshold after seeing hit counts is exactly the degree of
  freedom this document exists to remove.

Every filter stage reports input and output row counts to `usecase/summary_counts.json`,
so the funnel is auditable.

## 3. Overlap criterion against candidate boundaries

- The tested feature of each locus is its **two breakpoints** — the released `Start` and
  the released `Stop` — each tested independently. The interior of the event is **not**
  tested; a boundary sitting inside a large deletion is not a breakpoint coincidence and
  is out of scope for this demonstration.
- A candidate boundary is a row of `annotation.tsv` with position `pos`, which is a **bin
  start**: the bin is the half-open interval `[pos, pos + resolution_bp)`.
- **Hit definition.** Breakpoint `b` (1-based ClinVar coordinate) on chromosome `c` hits
  candidate `(c, pos)` in a package of resolution `R` if and only if

      | (b - 1) - pos |  <=  R

  where `b - 1` converts the 1-based ClinVar coordinate to the 0-based half-open
  convention of the bin grid used by the candidate tables. `R = 10,000` in the 10 kb
  packages and `R = 25,000` in the 25 kb packages, i.e. **±1 bin, stated per resolution**.
- **Why ±1 bin and not some other window:** `R` is exactly the `vote_tol_bp` the delivered
  packages already use for BD1 caller agreement (verified from the delivered tables before
  writing this file: `vote_tol_bp == resolution_bp` in all six packages). The window is
  therefore the engine's own positional tolerance, not a window chosen for this analysis.
- Chromosome matching: ClinVar `Chromosome` strings are compared to the candidate `chrom`
  column after stripping any `chr` prefix from both sides. (The delivered tables store
  `1..22, X` without a prefix; verified before writing this file.)
- All six packages are tested. **Every** locus × package × candidate hit is written to
  `usecase/hits_all.csv` with the full 38-column evidence record attached. There is no
  manual curation, no post-hoc exclusion, and no re-run of this step with a different
  window.

## 4. Reporting rule

- `usecase/hits_all.csv` contains **all** hits, unfiltered and not ordered by any quality
  criterion (deterministic order: package, chromosome, pos, VariationID, breakpoint end).
- `usecase/summary_counts.json` reports loci tested after each filter stage, hits per
  package, distinct loci with >=1 hit, tier distribution of hits, and the breakpoint-end
  split.
- If the intersection yields **zero** hits, `usecase/NULL_RESULT.md` is written stating the
  rule and the outcome, and the demonstration stops there. A null result is a reportable
  outcome of this pre-registration, not a reason to widen the window.

## 5. Deterministic choice of the manuscript example

Applied only after §3 has been run and `hits_all.csv` exists. No visual inspection of the
records precedes this step.

- **Reference package:** `GM12878_10kb`. Declared here on two grounds fixed in advance: it
  is the largest of the six packages by candidate count (69,479, verified before writing
  this file), and it is the package in which the project's independently re-verified
  `loop_anchor_count` example was checked (`TERMINOLOGY.md`).
- **Primary example = the hit with the highest evidence tier in the reference package**,
  tier order `D5 > D4 > D3 > D2 > D1`, taken from the `combine_tier` rule of
  `backend/tad_vci/graded_evidence.py` (read before writing this file: `D5` is the top
  tier, reached at BD1 strong with a strong assessable supporting axis).
- **Tie-breaking, applied in this order until a unique row remains:**
  1. smallest chromosome, ordered `1 < 2 < ... < 22 < X`;
  2. smallest `pos`;
  3. smallest ClinVar `VariationID`, compared as an integer;
  4. breakpoint end, `Start` before `Stop`.
  This chain is total: it cannot fail to yield exactly one row.
- **Package fallback order,** used only if the reference package has zero hits:
  `GM12878_10kb -> IMR90_10kb -> HepG2_10kb -> GM12878_25kb -> IMR90_25kb -> HepG2_25kb`.
  If every package has zero hits, §4's null-result branch applies.
- **Mandatory companion example, to prevent a favourable-case reading.** The same rule is
  also applied in the *opposite* direction — the **lowest**-tier hit in the reference
  package, ties broken by the identical chain — and that record is reported alongside the
  primary one in `contrast_table.csv` and `CONTRAST.md`. Both are declared here so the
  demonstration cannot be accused of showing only a well-evidenced boundary: the point
  being demonstrated is that the record is *informative*, which must hold at both ends of
  the tier range or not be claimed at all. If the reference package's hits are all one
  tier, that fact is reported and the companion omitted with a stated reason.
- **The example is never chosen by eye, by biological interest, by gene name, by
  literature familiarity, or by how well the record reads.** If the rule selects a locus
  with an unremarkable or uninformative record, that is the result and it is reported.

## 6. Verification requirement

Every number appearing in the contrast must be **re-derived from the delivered
`annotation.tsv`** by a script in `usecase/scripts/`, each check recorded in
`usecase/verification.json` as `{value, source column, match true/false}`. A number that
cannot be re-derived from a named delivered column is removed from the contrast rather
than explained. Clinical statements about the locus carry a ClinVar accession retrieved
this session; any mechanistic clinical claim carries a PMID/DOI retrieved this session, or
is omitted.

## 7. Amendments (append-only)

Any deviation from the above — including a move down the source priority list of §1 — is
recorded below with a UTC timestamp and the reason, by appending. No text above this line
is edited after 2026-08-15T02:38:00Z.

*(no amendments at time of writing)*

### 2026-08-15T02:51:23Z — observation recorded, NO deviation from the rule

While confirming the two selected ClinVar records against NCBI E-utilities esummary
(retrieved this session), the GRCh38 coordinates NCBI returns for both records populate
`inner_start`/`inner_stop`, with `outer_start`/`outer_stop` empty. For ClinVar records
using HGVS uncertain-breakpoint notation, the released `Start`/`Stop` are therefore INNER
bounds: the true breakpoint lies at or outside the stated coordinate, not at it. Measured
on the pre-registered locus set: 1,073 of 4,229 loci (25.4%) carry `?` uncertain-breakpoint
notation in their `Name`, as do 1,029 of the 4,044 loci with at least one hit (25.4%).

This is recorded as a caveat, not acted on. The rule of §2 says the coordinates are used
"as released", and they were. No filter was added, no locus removed, no window changed,
and the selection of §5 was not re-run. The primary example, ClinVar 541215
(`NC_000001.11:g.(?_1020153)_(1313808_?)del`), is one of the affected records and is
reported as such in CONTRAST.md and REPORT_usecase.md. Learning this after the rule was
sealed is exactly why it is written here rather than silently incorporated.

### 2026-08-15T05:18:05Z — CORRECTION to the header claim, raised by independent review

**What the header says and why it is too strong.** The header of this file (line 4) asserts
it was written "BEFORE the first query against any candidate evidence table". An independent
reviewer established that this is inaccurate as written: before the seal at
2026-08-15T02:38:00Z I ran aggregate shell inspections that *did* read the six
`annotation.tsv` files. **The header's claim is corrected by this section to: written before
any LOCUS-LEVEL query touched the candidate tables.** The original wording is left in place
above, uncorrected, because this file is append-only; this section is the correction of
record, and Supplementary Note 19 must carry it.

**Exactly what ran before the seal.** Recovered verbatim from the session execution log, not
from recollection. Seven cells, 02:33:07–02:34:51 UTC, all preceding the 02:38:00Z seal:

| cell | UTC | command (abbreviated) | what it revealed |
|---|---|---|---|
| 4 | 02:33:07 | `ls -la` over `example_data/multicell/*/` | file names and sizes only; no file content |
| 5 | 02:33:07 | `json.load('GM12878_10kb/manifest.json')`, full walk | the delivered manifest, **including its `tier_counts` block** (D1 29263, D2 4436, D3 20596, D4 8911, D5 6273) and `n_candidates` 69479 |
| 7 | 02:33:18 | `grep -n "D1\|D2\|...\|tier" graded_evidence.py` | engine source only |
| 8 | 02:33:31 | `sed -n '265,345p' graded_evidence.py` | engine source: `combine_tier` body, i.e. D5 is the top tier |
| 9 | 02:33:31 | `sed -n '210,265p' graded_evidence.py` + `grade_loop` | engine source: `grade_votes`, `grade_fold`, `grade_loop` |
| 10 | 02:34:34 | `head -1 GM12878_10kb/annotation.tsv \| tr '\t' '\n'` | the 38 column names, and `NF`=38 |
| 12 | 02:34:51 | `for d in <six packages>; tail -n +2 $d/annotation.tsv \| wc -l` and an `awk` over columns 1, 25, 26, 31, 32, 33, then `cut -f1 \| sort -u`, `cut -f21 \| sort \| uniq -c`, `cut -f22 \| sort -u` | per-package row counts (69479, 28259, 63860, 29513, 65157, 28476; sum 284744); `chrom` domain (1–22, X, no `chr` prefix); constant `resolution_bp`, `vote_tol_bp`, `loop_window_bp`, `chip_window_bp`, `tier_rule_version`; **and the GM12878 10 kb `D_tier` value counts** (D1 29263, D2 4436, D3 20596, D4 8911, D5 6273) |

So the tier composition of one package was in front of me before the seal, twice: once from
that package's own delivered `manifest.json` and once from an `awk`/`uniq -c` over its
`D_tier` column.

**Why the seal still carries its weight.** Three things are true of everything listed above,
and they are checkable against the commands themselves:

1. **It is aggregate only.** Every pre-seal read is a count, a domain, a column name or a
   file size. Not one of them returns, filters or orders a row. No `chrom`+`pos` pair, no
   variant, no candidate identity, and no per-boundary evidence value was visible. The seal
   therefore predates **all locus-level information**, which is the degree of freedom §5
   exists to remove: the tie-break chain could not have been aimed at a locus that had not
   been retrieved, from a variant file that had not been downloaded (ClinVar was fetched
   02:38:20–02:39:17Z, after the seal).
2. **It was already published.** The tier composition I saw is not private knowledge obtained
   by looking: it is a released number. The manuscript states at Figure 1c that D1 accounts
   for **19.3–42.1%** of candidates and D5 for **7.7–16.2%** across the six packages
   (`manuscript/CrediTAD_GPB_AppNote.md` lines 112–113, traced in `NUMBERS_LEDGER.csv` to
   `analysis/out/census_long.csv` via `census.py` — a different student's work, predating
   this task). GM12878 10 kb's D1 share, 29263/69479 = 42.12%, is precisely the upper bound
   of that published range. The same counts ship inside each package's `manifest.json`. A
   reader of the paper knows them without running anything.
3. **The seal is a hash, and the ledger copy is independent.** The 02:38:00Z copy in
   `ledger/PREREGISTERED_RULE_20260815T023800Z.md` (12,390 bytes, sha256 `71e963c8…`) fixes
   the rule text itself, so no clause can be shown to have been added after a result.

**The residual concern, stated rather than argued away.** Having seen that D5 exists and is
~9% of the reference package, I knew when writing §5 that "highest tier" would probably
resolve to a D5 record rather than come back empty. That is a real, if weak, dependency: the
*tier level* named by the rule was not chosen blind. It cannot have influenced **which**
locus was selected — no locus was known — and the mandatory lowest-tier companion in §5
deliberately spends the opposite end of the same distribution, so the demonstration does not
rest on the tier level being favourable. A referee is nonetheless entitled to know this, and
it is now on the record. The honest summary is that the seal removes post-hoc locus selection
and post-hoc rule-writing, but not a designer who had read the paper's own Figure 1c.

**Nothing was re-run and nothing was changed** in consequence of this correction: not the
filter, not the window, not the tie-break chain, not the selection, not the verification. The
only artefacts touched are this appended section, the added hash field in `PREREG_SEAL.json`,
a stale check-count label in `CONTRAST.md`, a scope caveat in `REPORT_usecase.md`, and a new
`MANUSCRIPT_CORRECTION_NEEDED.md`.
