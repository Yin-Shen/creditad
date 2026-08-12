# Tests not included in the public repository

Five test modules from the development tree are **not** shipped here, because the code
they exercise is not shipped here either. They test scripts in the project's retired
`backend/analyses/` namespace (superseded hg19 / `max_support_v1` analysis generation,
~1.8 GB), which is excluded from the public repository by design.

| test module | imports | lives in |
| `test_bcp_config_consistency.py` | `bcp_config` | `backend/analyses/tad_vci_20260529/` |
| `test_bcp_conformal.py` | `conformal_student` | `backend/analyses/tad_vci_20260529/` |
| `test_gm_annotation_builder_gates.py` | `build_gm12878_annotations` | `backend/analyses/tad_vci_rebuild_20260710/` |
| `test_gm_calibration_provenance.py` | `generate_gm12878_calibration` | `backend/analyses/tad_vci_rebuild_20260710/calibration/` |
| `test_ice_sidecar_manifest.py` | `compute_ice_sidecar` | `backend/analyses/tad_vci_rebuild_20260710/` |

Removing them changes no production coverage: none of them import `tad_vci`,
`mcool_bed`, or the CLI. The shipped suite covers the evidence engine, the tier rule, the
CLI, the annotation path, the curation/review record and the HTTP input guards.

The development tree retains all of them.
