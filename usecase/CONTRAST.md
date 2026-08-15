# Two-column contrast — what a coordinate list shows vs what the evidence record shows

Generated 2026-08-15T02:52:14Z by
`usecase/scripts/03_select_and_contrast.py`. Machine-readable form:
`usecase/contrast_table.csv`. Every number below is re-derived from the delivered
`annotation.tsv` by `usecase/scripts/04_verify.py`
(68/68 checks, all match) and recorded in
`usecase/verification.json`.

**What this demonstration claims.** That the evidence record changes how a person *handles*
one boundary. Nothing here is evidence that a boundary is correct, real, validated or
disease-causing, that CrediTAD called or recovered it, or that overlap with a pathogenic
variant supports the boundary. The ClinVar locus supplies only the reason someone is
looking at this position at all.

**How these two boundaries were chosen.** By the rule sealed in
`PREREGISTERED_RULE.md` §5 at 2026-08-15T02:38:00Z, before any query touched the candidate
tables: highest-tier hit in GM12878 10 kb (ties: chromosome, then position, then ClinVar
VariationID, then breakpoint end), plus — declared in the same sealed file — the
lowest-tier hit under the identical chain, so the demonstration cannot show only a
well-evidenced boundary. Neither was chosen by eye. Both are 2 of 4,569 hits in that
package; 400 hits sit at D5 and 1,624 at D1.

---

## Example A (primary, rule = highest tier): chr1:1,310,000–1,320,000

Reason for attention: **ClinVar 541215 / VCV000541215 / RCV000651427**, a Pathogenic
deletion at chr1:1,020,153–1,313,808 (293,655 bp), asserted for *Congenital myasthenic
syndrome 8*, review status "criteria provided, single submitter", gene *AGRN*. Its `Stop`
breakpoint lies 3,807 bp from the candidate boundary position — inside the pre-registered
±1 bin (10,000 bp) window.

| | (a) coordinate list shows | (b) evidence record shows |
|---|---|---|
| position | `chr1:1310000-1320000` — a bin start, and nothing else | identical position |
| caller support (BD1) | *not shown* | 5/6 callers: TopDom;SpectralTAD;MSTD;arrowhead;DI; grade strong |
| which callers dissent | *not shown* | OnTAD did not support this position (panel minus supporting callers) |
| caller offsets | *not shown* | 0;0;0;-10000;0 bp from the candidate position |
| CTCF (BD2) | *not shown* | grade strong; fold 2.3671 vs background 0.61100; null percentile 98.51 |
| cohesin (BD3) | *not shown* | grade weak; fold 1.3835 vs background 0.59390; null percentile 89.41 |
| loop-anchor engagement (BD4) | *not shown* | grade none; loop_anchor_count = 0 loop ends in a 50000 bp window |
| evidence tier | *not shown* | D5 (rule max_support_v2) |
| derivation | *not shown* | `BD1 strong + supporting strong → D5 (best supporting: BD2=strong)` |
| what was not measured | *not shown* | all four criteria were assessable at this position |
| threshold provenance | *not shown* | BD2/BD3 bands from GM12878 (packaged), calibrated=True; ChIP window 50000 bp |
| panel definition | *not shown* | 6 callers; `grade_votes` fraction rule, ±1 bin vote tolerance |

**How the record changes handling.** The coordinate list supports one action: note that a
boundary was called near the deletion's distal breakpoint. The record supports a different
set, and each is licensed by a specific field rather than by the tier label:

1. **The tier is D5, but it rests on one axis.** The derivation names its own basis —
   `BD2=strong` — and the record shows the other two measured axes did not agree with it:
   cohesin is `weak` (fold 1.3835, 89.41st percentile of this cell's random-window null)
   and loop-anchor engagement is `none` (`loop_anchor_count` = 0 in a 50 kb window). A
   reader who treats D5 as a summary verdict would miss that. A reader who reads the record
   sees a CTCF-supported position without measured loop engagement, and can decide whether
   that is the kind of support their question needs.
2. **The support is panel-conditional, and the record says by how much.** 5 of 6 callers,
   with OnTAD dissenting, and one supporter (arrowhead) placed at −10,000 bp rather than 0.
   Under the shipped `grade_votes` rule at n=6 (`strong` requires votes ≥ ⌈0.8×6⌉ = 5),
   this boundary sits on the *threshold*: holding its molecular grades fixed and dropping
   one supporter takes BD1 to `moderate` and the tier from D5 to D4; dropping two takes it
   to D3. The list gives no handle on that sensitivity. The record gives the vote count,
   the panel, the tolerance and the rule string, so the sensitivity is computable by the
   reader rather than hidden in the label.
3. **The unmeasured is separated from the unsupported.** All four criteria were assessable
   here, so `none` on loop-anchor engagement means measured-and-not-supported, not
   never-measured. That distinction is what `not_assessable_criteria` exists to carry, and
   it is the distinction that decides whether the next experiment is a new assay or a
   different interpretation.
4. **The thresholds are inspectable.** The grade `strong` is not a bare word: it is
   fold 2.3671 against a GM12878-calibrated cutpoint, with the observed value at the
   98.51st percentile of that cell's own random-window null, in a stated 50 kb window. A
   referee can disagree with the cutpoint and re-grade without re-running the tool.

## Example B (companion, rule = lowest tier): chr1:830,000–840,000

Reason for attention: **ClinVar 58037 / VCV000058037 / RCV000051780**, a Pathogenic copy
number gain at chr1:826,553–4,719,105 (3,892,552 bp), dbVar nsv530301, phenotype recorded
only as "See cases". Its `Start` breakpoint lies 3,448 bp from the candidate position.

| | (a) coordinate list shows | (b) evidence record shows |
|---|---|---|
| position | `chr1:830000-840000` | identical position |
| caller support (BD1) | *not shown* | 1/6 callers: TopDom; grade weak |
| CTCF (BD2) | *not shown* | grade none; fold 0.6297 vs background 0.61100; null percentile 12.10 |
| cohesin (BD3) | *not shown* | grade none; fold 0.5081 vs background 0.59390; null percentile 11.22 |
| loop-anchor engagement (BD4) | *not shown* | grade none; loop_anchor_count = 0 loop ends in a 50000 bp window |
| evidence tier | *not shown* | D1 (rule max_support_v2) |
| derivation | *not shown* | `BD1 weak + no positive supporting → D1 (best supporting: no supporting signal)` |
| what was not measured | *not shown* | all four criteria were assessable at this position |

**How the record changes handling.** In the coordinate list these two rows are
indistinguishable in kind: both are one line, `chr1` and a position, each near the
breakpoint of a Pathogenic variant. The record separates them. Here a single caller
(TopDom) supports the position; CTCF and cohesin were measured and sit near the bottom of
this cell's null distribution (12.10th and 11.22nd percentiles); no loop ends fall in the
window. All four criteria were assessable, so this is a position where the measurements
were taken and none of them returned support — reported as *no CTCF, cohesin or
loop-anchor support was measured at this position*, which is not a statement that the
boundary is false. CTCF-independent boundaries exist, and three correlated readouts of one
extrusion mechanism cannot rule one out. What the record does license is a triage
decision: two boundaries that a coordinate list ranks equally are not equally worth an
experiment, and the record states the grounds for the ordering.

---

## Caveats attached to this contrast

- **ClinVar breakpoint coordinates are inner bounds for a quarter of the set.** NCBI's
  esummary returns the GRCh38 coordinates of both selected records in `inner_start`/
  `inner_stop`, with outer bounds empty. Example A's HGVS name,
  `NC_000001.11:g.(?_1020153)_(1313808_?)del`, uses uncertain-breakpoint notation: the true
  breakpoint is at or beyond 1,313,808, not necessarily at it. 1,073 of 4,229 loci (25.4%)
  and 1,029 of 4,044 hit loci (25.4%) carry this notation. The measured distance of 3,807 bp
  is therefore a distance to a *reported* bound, and the coincidence is at the resolution of
  the variant call, not of the boundary.
- **A 3,807 bp separation is sub-bin, not sub-kilobase agreement.** At 10 kb resolution the
  bin is the unit of position; the window was set to the engine's own vote tolerance
  (`vote_tol_bp` = `resolution_bp`), not to a distance chosen to make this pair look close.
- **Breakpoint proximity is not mechanism.** Nothing here shows the deletion disrupts this
  boundary, or that boundary disruption contributes to the phenotype. Testing that would
  require the patient's own contact map, which is not in this corpus.
- **The clinical assertion is one submitter's,** carried at ClinVar's stated review status,
  and it is reported rather than endorsed. *AGRN* variants are an established cause of
  congenital myasthenic syndrome (PMID 36835142, doi:10.3390/ijms24043730; PMID 38696726,
  doi:10.1093/brain/awae124, both retrieved this session), but no reviewed source retrieved
  here attributes this particular deletion's effect to TAD-boundary disruption.
- **Cell-type mismatch.** *AGRN*-related congenital myasthenic syndrome is a
  neuromuscular-junction disorder; GM12878 is a lymphoblastoid line. The record describes
  chromatin architecture in GM12878, not in the affected tissue. This is a limitation of the
  demonstration, and it is a consequence of the pre-registered rule rather than an oversight:
  the rule fixed the reference package before any locus was known.
