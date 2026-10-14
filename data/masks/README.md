# hg38 mask for the ChIP axes

`hg38_mask.bed.gz` lists the hg38 regions where a CTCF or RAD21 ChIP readout is not a
measurement. With `creditad annotate --mask hg38` (or `--mask data/masks/hg38_mask.bed.gz`),
a boundary whose position (BED column 2) lies inside one of these intervals, half-open
`[start, end)`, gets BD2 and BD3 recorded as `not_assessable`, both named in
`not_assessable_criteria`, and its tier re-derived with the released rule. BD1, BD4 and the
raw `ctcf` / `rad21` / fold / null-percentile columns are not changed. Masked rows are listed
in `tadvci_mask_records_<tag>.tsv`. This is the masking of the v1.1 evidence tables in
`data/multicell_v1.1/`; `scripts/mask_and_retier.py --check data/multicell_v1.1` rebuilds
those tables from the v1.0 tables with it and compares them byte for byte.

The same file ships inside the package as `backend/tad_vci/masks/hg38_mask.bed.gz` (the
file `--mask hg38` reads in an installed copy); a test asserts the two are identical.

## Contents

BED4, sorted, gzip-compressed with a zero timestamp; 1,583 intervals on chr1–22, chrX and
chrY; the fourth column names the source of each interval. Intervals are written as in the
sources, not merged (merged, they cover 266,135,336 bp in 1,039 intervals).

| label | source | intervals |
|---|---|---|
| `assembly_gap` | UCSC hg38 `gap` table | 603 |
| `cytoband_acen` | UCSC hg38 `cytoBand` table, gieStain `acen` (centromeres) | 48 |
| `cytoband_gvar` | UCSC hg38 `cytoBand` table, gieStain `gvar` (variable heterochromatin) | 17 |
| `cytoband_stalk` | UCSC hg38 `cytoBand` table, gieStain `stalk` (acrocentric stalks) | 5 |
| `encode_blacklist` | ENCODE exclusion list ENCFF356LFX (GRCh38) | 910 |

`hg38_mask.bed.gz`: 14,547 B, md5 `02fbf9ec52b172298f690b8c9b316cc0`, SHA-256
`053e6a188dcf943855069ffb5bbd8ff9fe993fb7cd7f774e79744d82578ff829`.

`hg38.chrom.sizes` holds the lengths of chr1–22, chrX, chrY and chrM from the UCSC file of
the same name (md5 `75afe0a7c5316fe36c565fd9c77e8b27`). `scripts/mask_and_retier.py` uses it
to remove the candidates lying at or beyond the end of their chromosome, as the v1.1 tables
do; `--mask` itself does not remove rows.

## Sources

| file | URL | retrieved | md5 | SHA-256 |
|---|---|---|---|---|
| `gap.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz | 2026-10-13 | `da739b566b6a0db5eb0dd6c8c206adb6` | `c466135bf8553f007960b3bd7d747343d014143a0415c127d1452bbe484bf1c1` |
| `cytoBand.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz | 2026-10-13 | `6677dd50caf19b0c3626b99ce67ca073` | `e514e48b0a12a516fd5231b44241a20a282f290bdd6131cadcddc6276ee26082` |
| `ENCFF356LFX.bed.gz` | https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz | local copy dated 2026-04-17; md5 equal to the ENCODE record, checked 2026-10-13 | `393688b4f06c9ce26165d47433dd8c37` | `a9d086ce90ca67f933b29adfbff56ea768bbce2ec4b0d51e2d405e5a1d61bb56` |
| `hg38.chrom.sizes` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes | 2026-10-13 | `c42e9f75fff906c6a7143cb2fac86602` | `e1d9152418038457a959e949d99c9caf7ae3e4f87cfbb4ffe7d8b9a54ba1202b` |

The UCSC tables carry no version number; the checksums above identify the copies used.
ENCFF356LFX is the ENCODE annotation ENCSR636HFF ("exclusion list regions", GRCh38, Kundaje
lab), described in Amemiya, Kundaje and Boyle, *Sci Rep* 9, 9354 (2019).

Rebuild (same bytes from the same inputs):

```bash
python scripts/build_hg38_mask.py --gap gap.txt.gz --cytoband cytoBand.txt.gz \
    --blacklist ENCFF356LFX.bed.gz --out data/masks/hg38_mask.bed.gz
```

## Relation to the v1.1 tables

The v1.1 tables were masked with two files: an hg38 gap file of 243 intervals (md5
`e6071892d973f45d3d50419765c1f24f`) and ENCFF356LFX. That gap file is, interval for
interval, the merge (overlapping and abutting intervals joined) of the UCSC `gap` rows on
chr1–22, chrX and chrY with the UCSC `cytoBand` acen, gvar and stalk bands of the copies
listed above: it covers every base of the UCSC gap table on those chromosomes (150,610,728
bp) and 112,126,770 bp of centromeric and heterochromatic bands besides. `hg38_mask.bed.gz`
keeps the same sources unmerged, so it masks exactly the same positions while its labels say
which source each interval comes from. In the v1.1 masking records, `assembly_gap` covers
both the `assembly_gap` and the `cytoband_*` labels.

## Licence

UCSC: "All the files and tables in this directory are freely usable for any purpose"
(README.txt of https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/, retrieved
2026-10-13). ENCODE: "External data users may freely download, analyze and publish results
based on any ENCODE data without restrictions" (ENCODE data use policy,
https://www.encodeproject.org/help/citing-encode/, retrieved 2026-10-13). The files in this
directory are derived from these sources and are redistributed under those terms; please cite
the ENCODE blacklist paper above when using the mask.
