# Legacy class-API cleanroom_v2 tests (quarantined 2026-05-13)

These 8 test files were written against the old class-based API
(`InsulationScoreBoundaryCaller().call(...)`) and against the 5-method
set including `DIZeroCrossingBoundaryCaller` (DI).

On 2026-05-13 the production set was reduced to 4 methods
(`insulation`, `topdom_like`, `contact_contrast`, `spectral_profile`)
and the module was switched to a function-based API
(`call_boundaries(matrix, chrom, resolution, region_start=0, params=None)`).
See `../../../../../TADSuite/evidence/tad_boundary_run_20260513_002918/DECISION_DROP_DI.md`.

These tests are preserved here for reference. They are NOT picked up
by pytest because the parent directory uses a `test_` prefix convention
and pytest does not recurse into directories whose names do not match
the conftest discovery pattern. To re-enable a test, port it to the
new function API and move it back into `backend/tests/`.

A green test suite against the new API should be written as a separate
task (not part of the 2026-05-13 method-replace change).
