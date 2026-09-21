## CrediTAD v1.0.3

Per-boundary evidence records for reviewing called TAD boundary lists. This tag is the repository state
cited by the article. Code paths and the six v1.0 evidence tables are unchanged from v1.0.2.

* Every package manifest (`data/multicell/<CELL>_<RES>/manifest.json`) now records the HiCCUPS loop file
  it was graded against: shipped path, MD5 of the gzipped and of the uncompressed file, the upstream
  GEO / ENCODE accession and the liftOver applied.
* `data/multicell_v1.1/`: release v1.1 of the six evidence tables (BD2 and BD3 recorded as
  `not_assessable` in hg38 assembly gaps and the ENCODE blacklist; 12 out-of-range rows removed), with
  masking records, per-package impact tables and `scripts/mask_and_retier.py`.
* `pyproject.toml` version 1.0.3, so `pip install` reports the tagged version.

The desktop binaries (Linux AppImage, Windows installer and portable archive) remain attached to the
v1.0.2 release, unchanged since v1.0.1; their SHA-256 digests are listed in `RELEASE_ASSETS.sha256`.
