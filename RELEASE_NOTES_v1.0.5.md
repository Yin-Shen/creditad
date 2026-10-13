# CrediTAD v1.0.5

Bug fixes for loop-file and ChIP-track input, the demo packages and read-only installs, a
widened pandas range, and corrected desktop documentation. The tier rule, the grading and the released
evidence tables are unchanged: `scripts/reproduce_packages.py` regenerates all six packages
exactly with this release (284,744 rows, 38 of 38 columns identical).

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

## Desktop builds attached to this release
- `CrediTAD-Linux.AppImage`: desktop client 1.0.2 with the backend of this release
  (540,753,929 B, SHA-256 `81061aef1133d27e758dd86cf23b0fe3b89dfc2cbe6a8bac0312542896bea51e`).
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
- Tests: `backend/tests/tad_vci/test_bigwig_chromosome_without_data.py`,
  `backend/tests/tad_vci/test_loops_encode_juicer_layout.py`,
  `backend/tests/tad_vci/test_multicell_demo_gz_beds.py`,
  `backend/tests/tad_vci/test_output_dir_read_only_install.py`.
