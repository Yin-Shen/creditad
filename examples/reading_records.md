# Reading the records: a follow-up filter

The evidence table has one row per boundary. Two of its columns answer the question
"which boundaries in this list have both kinds of support?": `BD1_caller_support`
(caller agreement) and the three molecular grades `BD2_ctcf`, `BD3_rad21`,
`BD4_loop_anchor`. The filter below is the reading used for Figure 3D of the article. It
is a recipe on top of the table, not a feature of the software, and it adds no threshold
to those of the criteria.

```python
import pandas as pd

t = pd.read_csv("tadvci_annotation_beds_<tag>.tsv", sep="\t")
rank = {"none": 0, "weak": 1, "moderate": 2, "strong": 3}
mol = t[["BD2_ctcf", "BD3_rad21", "BD4_loop_anchor"]]

measured = mol.where(mol != "not_assessable")
best = measured.apply(lambda c: c.map(rank)).max(axis=1)        # NaN if nothing measured

t["caller_support"] = t["BD1_caller_support"].isin(["strong", "moderate"]).map(
    {True: ">=4 of 6", False: "<=3 of 6"})                       # for a six-caller panel
t["molecular_support"] = best.map({0: "none", 1: "weak only",
                                   2: "moderate or strong", 3: "moderate or strong"}
                                  ).fillna("not measured")

print(pd.crosstab(t["caller_support"], t["molecular_support"]))

both = t[(t["caller_support"] == ">=4 of 6") & (t["molecular_support"] == "moderate or strong")]
callers_only = t[(t["caller_support"] == ">=4 of 6") & (t["molecular_support"] == "none")]
needs_data = t[t["not_assessable_criteria"].fillna("") != ""]
```

- `both`: caller agreement and an external readout coincide at the boundary.
- `callers_only`: strong caller agreement with no measured CTCF, cohesin or loop-anchor
  support; these are recorded at D3. They may still insulate in the contact map; the
  record says only that no external readout coincides with them.
- `needs_data`: a criterion could not be assessed because its input was not supplied;
  the column lists which. An unmeasured criterion is never counted as negative.

BD1 grades are fractions of the panel (strong >= ceil(0.8 n), moderate >= ceil(0.6 n)), so
"at least four of six" is BD1 moderate or strong only for a six-caller panel.
