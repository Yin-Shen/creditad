# chr7 preview tracks — DISPLAY ONLY, NOT ANALYSIS INPUTS

These bigWigs are 250 bp re-binned subsets of the genome-wide ENCODE tracks,
built for browser display. They are **not** value-identical to the full tracks:
measured median window-mean ratio 1.118 (CTCF) and 1.109 (RAD21), Pearson r 0.935
and 0.909, with **0 of 1,528** chr7 windows agreeing to within 1e-6.

Using them in place of the genome-wide tracks reproduces only 28 of 38 columns of
the delivered tables and gives evidence-tier agreement 0.903. With the genome-wide
ENCODE tracks the same command reproduces all six packages exactly (38/38 columns,
tier agreement 1.000000).

For reproduction, fetch the genome-wide tracks — see ../../REPRODUCE.md.

The matching chr7 .mcool contact maps are not in this repository (42.1, 48.3 and
90.8 MB, above the repository's file-size policy); they ship inside the desktop builds
attached to the v1.0.1 release, under resources/example_data/chr7_tracks/.
