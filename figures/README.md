# Figure source data

Source data for the published Hi-C contact-map panel of Figure 1, with the script that
exported it. This is the numeric input behind the drawn panel: **no value here is
rescaled, logged or coloured** — these are the cooler's own ICE-balanced contact values.

| file | what it is |
|---|---|
| `source_data/fig1_hic_contact_matrix.tsv` | the exported contact values |
| `source_data/fig1_hic_contact_matrix_meta.json` | export parameters and provenance |
| `source_data/fig1_hic_contact_matrix.sha256` | checksum of the `.tsv` |
| `scripts/export_hic_matrix.py` | the exporter that produced all three |

Verify:

```bash
cd figures/source_data && sha256sum -c fig1_hic_contact_matrix.sha256
```

## Provenance

| field | value |
|---|---|
| region | `chr7:86375000-87475000` |
| resolution | 10,000 bp |
| balancing | ICE (cooler 'weight' column) |
| geometry | rotated_triangle |
| cells exported | 2,448 |
| `.tsv` sha256 | `26b99da04f4d303bda6c5ead403ec195e04cccda06b77909a6c34932d0239209` |
| exporter sha256 | `9a9555b14363af3db6d0ae4cf8e8ccc5dec6e7cb3660625c251229872a366d66` |

`meta.json` is the authoritative record and carries the full parameter set, including the
colour-scale rule the application uses.

## Re-running the export

The exporter reads the GM12878 chr7 multi-resolution cooler
(`GM12878_chr7_hg38.mcool`, 90.8 MB). **That file is not in this repository** — it exceeds
the file-size policy and ships inside the desktop builds attached to the
[v1.0.2 release](https://github.com/Yin-Shen/creditad/releases/tag/v1.0.2), under
`resources/example_data/chr7_tracks/`.

`export_hic_matrix.py` is shipped **verbatim as it ran**, so its checksum above is
meaningful; it therefore still carries the absolute `MCOOL` path from the machine it ran
on. To re-run, point that constant (near the top of the file) at your own copy of the
cooler:

```python
MCOOL = "/path/to/resources/example_data/chr7_tracks/GM12878_chr7_hg38.mcool"
```

Requires `cooler` and `numpy`. The script writes all three files above into
`figures/source_data/`.
