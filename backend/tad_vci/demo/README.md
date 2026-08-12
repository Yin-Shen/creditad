# Retired legacy browser demo

The pre-audit browser demo was moved intact to
`backend/analyses/tad_vci_20260529/_quarantined_tad_vci_demo/` on 2026-07-10.
It depended on retired BCP/conformal routing, BD5, an unauthenticated demo
actor record, and false immutable/signed wording.  It is not an active test,
application entry point, scientific artifact, or submission example.

Current behavior is covered by `backend/tests/tad_vci/`, the Vue frontend, and
the manifest-gated July rebuild.  Do not restore the legacy HTML, screenshots,
JSONL records or Playwright scripts to this active namespace.
