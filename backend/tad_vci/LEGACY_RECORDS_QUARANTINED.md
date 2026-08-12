# Legacy decision records are quarantined

The former `curation_records.jsonl` and browser-demo files predate the July 2026 audit. They contain
demo actor labels, stale evidence definitions and claims about conformal routing or expert
confirmation. They were moved intact to
`backend/analyses/tad_vci_20260529/_quarantined_tad_vci_demo/` and are no longer served, packaged as
active examples, or collected by the production test suite. `demo/README.md` is only a retirement
pointer.

New records are appended to `decision_records_v3.jsonl` using neutral `review_decision` and
`actor_label` fields plus a full SHA-256 checksum over the canonical record fields. The actor label is
not authenticated. This remains an application-level convention, not a signature, hash chain,
completeness proof or immutable storage.
