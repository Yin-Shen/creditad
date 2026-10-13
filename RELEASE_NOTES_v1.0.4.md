# CrediTAD v1.0.4

Bug fix, reproducibility material for the article, and packaging. The tier rule, the grading
and the released evidence tables are unchanged: `scripts/reproduce_packages.py` regenerates
all six packages exactly with this release (284,744 rows, 38 of 38 columns identical).

## Fixed
- **Gzipped loop files.** v1.0.3 opened a `.gz` loop file as text, failed on its first line
  and treated it as containing no loops, so BD4 was graded `none` at every boundary without an
  error. The repository ships its loop calls gzipped, so passing them straight to
  `creditad annotate` produced this. Loop files are now read plain or gzipped, and a supplied
  loop file that yields no loop anchors stops the run with an error instead of grading BD4.
  The command-line interface reports the error and exits with status 1; the local HTTP
  service answers the request with an error (HTTP 500, detail in the server log). Results
  computed from uncompressed loop files, including every released table and every number in
  the article, are unaffected.

## Changed
- The descriptive text of the v1.1 masking metadata (`data/multicell_v1.1/*/manifest.json`,
  `scripts/mask_and_retier.py`, `scripts/mask_common.py`) was simplified; no value, table or
  predicate changed, and the masking script regenerates the v1.1 tables byte for byte.

## Added
- `examples/quickstart/`: a run on files that ship in the repository, checked byte for byte
  against the expected output.
- `examples/reading_records.md`: reading the table for a follow-up decision.
- `Dockerfile`: the command-line interface, the test suite and the quick start in a clean Python environment.
- Tests for the loop-file fix (`backend/tests/tad_vci/test_loops_gzip_and_unparseable.py`).
