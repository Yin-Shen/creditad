# v1.1 evidence tables (gap and blacklist masking)

Release v1.1 of the six evidence tables records BD2 and BD3 as `not_assessable` for candidates in hg38
assembly gaps or the ENCODE blacklist (ENCFF356LFX), removes the 12 candidates lying beyond their
chromosome end, and re-derives the tier through the released rule. Each package directory holds
`annotation.tsv.gz`, `tier_boundaries.bed.gz`, `union_boundaries.bed.gz`, `masking_records.tsv.gz`
(which rows were masked and why) and `manifest.json`. The per-caller boundary BEDs are unchanged from
v1.0 and are read from `data/multicell/<CELL>_<RES>/`. `scripts/mask_and_retier.py` produces these tables
from the v1.0 tables; the per-package impact tables are in `impact/`. All numbers in the article describe
the v1.0 tables unless stated otherwise.
