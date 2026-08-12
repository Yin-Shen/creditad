# HiCCUPS loop calls — the BD4 input

One file per cell line for which this repository ships evidence packages. BD4 grades a
candidate boundary by whether it overlaps a loop anchor in its own cell line's calls;
`backend/mcool_bed.py` parses these `.bedpe` files and `data/multicell/*/manifest.json`
records which file each package used.

| cell line | loop calls | gzipped |
|---|---|---|
| GM12878 | 9,445 | 106 kB |
| IMR90 | 8,033 | 93 kB |
| HepG2 | 21,905 | 1360 kB |

`pandas.read_csv` and the parser both read the `.gz` form directly; nothing needs
unpacking.

**Loop-call density differs several-fold between cell lines** (HepG2 carries ~2.7x
GM12878's calls here), and BD4 grading depends on the density of the file supplied.
Absolute BD4 rates are therefore **not comparable across cell lines**, and no such
comparison is made anywhere in this repository or the accompanying paper.

A K562 loop file shipped in v1.0.0 and was removed in v1.0.1: no K562 evidence package
exists, no package manifest referenced it, and no released code path read it. It was a
leftover from development-era cell-line screening, and its presence implied a K562
result that does not exist.

Provenance: HiCCUPS calls on the corresponding Hi-C matrices, lifted to hg38. The
accessions of the matrices are in `../../REPRODUCE.md`.
