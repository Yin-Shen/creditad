# CrediTAD v1.0.5

Bug fixes for loop-file and ChIP-track input, the demo packages and read-only installs, an
optional mask that withholds the ChIP grades in gap and exclusion-list regions, a widened
pandas range, and corrected desktop documentation. The tier rule, the grading and the
released evidence tables are unchanged: `scripts/reproduce_packages.py` regenerates all six
packages exactly with this release (284,744 rows, 38 of 38 columns identical), and without
`--mask` every output is byte-identical to that of the previous code.

## Fixed
- **ENCODE / Juicer loop files.** Loop files in the BEDPE layout that ENCODE distributes
  (for example ENCFF794OIH.bedpe.gz: a header line starting `#chr1`, a
  `# juicer_tools version` line, 33 columns) were rejected with "yielded no loop anchors":
  v1.0.4 read the `#chr1` line as a named header and found no column called `chr1`. Lines
  starting with `#` are now skipped and the first six columns (chr1 x1 x2 chr2 y1 y2) are
  read. On ENCFF794OIH.bedpe.gz this gives the same 99,462 anchors on 23 chromosomes as a
  six-column copy of the file. A loop file that parses to no loops still stops the run with
  the same error.
- **ChIP tracks that list a chromosome without data.** bedGraphToBigWig writes every
  chromosome of the sizes file into the bigWig header, so a track made for some
  chromosomes lists all of them. Releases up to v1.0.4 read every window on a listed
  chromosome without data as zero signal and graded BD2/BD3 `none` at every boundary there.
  Such a chromosome is now treated like one missing from the header: BD2/BD3 are
  `not_assessable` and named in `not_assessable_criteria`. A window without data on a
  chromosome that has data is still read as zero signal, because tracks written with
  zero-valued intervals omitted carry no record there. The six released packages are
  unaffected (their tracks have data on every chromosome they list).
- **Demo caller tracks.** `POST /api/tadvci/use_multicell_example` registered only
  uncompressed `<caller>_boundaries.bed`, but this repository ships them as `.bed.gz`, so a
  demo served from the repository tree showed no caller tracks. The gzipped files are now
  used; an uncompressed file is still preferred when present.
- **Read-only installs.** When the backend runs from a read-only location (the Linux
  AppImage mounts its contents read-only), the Import page's default relative output
  directory `./tad_results` was answered with HTTP 400 "output_dir cannot be created:
  Read-only file system". A relative output directory now resolves under the per-user data
  directory that the token database already falls back to (`$XDG_DATA_HOME/CrediTAD` or
  `~/.local/share/CrediTAD`; `%LOCALAPPDATA%\CrediTAD` on Windows) when the working
  directory is not writable. The same applies to `POST /api/tad/run`. With a writable
  working directory nothing changes. Writability is judged with `os.access`, so a location
  that refuses writes without being reported as read-only is not detected.

## Changed
- **pandas.** The declared range is now `pandas>=1.5,<4` (`requirements.txt`:
  `pandas>=2.0.0,<4`). The desktop runtimes ship pandas 3.0.5. Under that runtime
  (Python 3.11.15, pandas 3.0.5, numpy 2.4.6) the test suite passes and
  `scripts/reproduce_packages.py` regenerates the six packages exactly; the six regenerated
  tables are byte-identical to those produced under Python 3.13.12 with pandas 2.3.1.
  `requirements.lock.txt` still records the Python 3.13 test environment (pandas 2.3.1).
- **Desktop documentation.** The README states which backend each desktop build embeds and
  records that the v1.0.1 Linux AppImage does not start the backend it bundles (its
  launcher looks for `runtime/python`; the bundled runtime has `runtime/bin/python` and
  lacks `cooler` and `pyBigWig`); the 1.0.2 AppImage attached here fixes this. The README
  and REPRODUCE.md now point to the desktop source repository,
  https://github.com/Yin-Shen/creditad-desktop.
- A development-tree path was removed from REPRODUCE.md and from a comment in
  `backend/sql.py`.
- **v1.1 masking scripts.** `scripts/mask_and_retier.py` and `scripts/mask_common.py` no
  longer contain machine-specific paths. They call the package's masking code (the code
  behind `--mask`) and take `--packages-dir`, `--mask`, `--chrom-sizes` and `--out`;
  `--check data/multicell_v1.1` compares every file written with the released v1.1 tables.
  All annotation tables, BEDs, masking records and summary tables come out byte-identical;
  in the manifests, the fields naming file locations and the generator differ, as do two
  fields copied from the v1.0 manifests (`hic_matrix` and the loop-file path in
  `loop_calls`). `data/multicell_v1.1/README.md` now states that the gap set includes the
  centromeric and heterochromatic cytobands.

## Desktop builds attached to this release
- `CrediTAD-Linux.AppImage`: desktop client 1.0.2 with the backend of this release
  (540,753,908 B, SHA-256 `63d2cfbd637ab1e1083b7a2d313aac23697fe451ed3f347836e69a715c938ba3`).
  It starts the backend it bundles: the client finds the interpreter at
  `runtime/bin/python`, and the runtime now includes `cooler` and `pyBigWig`. Mounted
  read-only, its Import page writes results under the per-user data directory (see
  "Read-only installs" above). Client source: tag v1.0.2 of
  https://github.com/Yin-Shen/creditad-desktop.
- `CrediTAD-Windows-Setup.exe` and `CrediTAD-Windows-1.0.1.zip`: the Windows builds of
  v1.0.1, carried forward unchanged (client 1.0.1, backend v1.0.1). They read a
  gzip-compressed loop file as empty, so decompress loop files before loading them there,
  and they grade BD2/BD3 `none` on a chromosome that a ChIP track lists without data.
- `RELEASE_ASSETS.sha256` lists the three files.

## Added
- **`--mask` for the ChIP axes.** `creditad annotate --mask BED` (and `mask_bed` in the body
  of `POST /api/tadvci/annotate`) records BD2 and BD3 as `not_assessable`, and names them in
  `not_assessable_criteria`, for every boundary whose position lies inside an interval of
  the BED (half-open, chromosome names compared without a leading `chr`); the tier of those
  boundaries is re-derived with the released rule. BD1, BD4 and the raw `ctcf` / `rad21` /
  fold / null-percentile columns are not changed, the 38-column table keeps its layout, and
  the masked rows are listed in `tadvci_mask_records_<tag>.tsv`. A masked boundary carries
  the same BD2, BD3, tier and `not_assessable_criteria` as a run without ChIP input.
  `--mask hg38` selects the mask that ships with the package (`data/masks/hg38_mask.bed.gz`,
  also installed as `tad_vci/masks/hg38_mask.bed.gz`): the UCSC hg38 gap table, the UCSC
  cytoBand acen, gvar and stalk bands and the ENCODE exclusion list ENCFF356LFX, chr1–22, X
  and Y, 1,583 intervals; sources, retrieval dates, checksums and licence terms are in
  `data/masks/README.md`, and `scripts/build_hg38_mask.py` rebuilds it byte for byte. It is
  the mask of the v1.1 tables: annotating the six packages with `--mask hg38` gives the
  annotation tables of `data/multicell_v1.1/` byte for byte on the 284,732 rows v1.1 keeps
  (4,474 masked rows); v1.1 also removes 12 candidates lying beyond their chromosome end,
  which `--mask` does not do. The bundled mask is refused for a run declared as another
  assembly.
- `data/masks/hg38.chrom.sizes` (chr1–22, X, Y, M; UCSC), used by `scripts/mask_and_retier.py`.
- Tests: `backend/tests/tad_vci/test_mask_option.py`,
  `backend/tests/tad_vci/test_bigwig_chromosome_without_data.py`,
  `backend/tests/tad_vci/test_loops_encode_juicer_layout.py`,
  `backend/tests/tad_vci/test_multicell_demo_gz_beds.py`,
  `backend/tests/tad_vci/test_output_dir_read_only_install.py`.
