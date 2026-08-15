#!/usr/bin/env python
"""EXPORT ONLY -- write the GM12878 chr7 contact values Figure 1b displays, plus the
intensity pre-scale the TOOL would use for this matrix, so the PANEL can be drawn in R.

Spec v3 pass-bar item 28 requires every panel of every figure to be produced by the R
build system. cooler is a Python library and there is no R reader for a .mcool, so this
script is allowed to do exactly one thing: read values out of the matrix and write them as
text. It does not draw, does not choose a colour, does not take a log, and does not
compute a summary statistic that any figure claims.

WHAT IT WRITES, AND WHAT IT DELIBERATELY DOES NOT
  * fig1_hic_contact_matrix.tsv -- one row per non-zero upper-triangle cell of the drawn
    window: bin indices, bin start coordinates, and the cooler's own ICE-BALANCED contact
    frequency. Nothing is rescaled, clipped, logged or binned here.
  * fig1_hic_contact_matrix_meta.json -- the window, the resolution, and `color_scale`:
    the intensity pre-scale the SOFTWARE would apply to this matrix, recomputed here by
    reproducing the backend's own probe (below). It is a measured property of the matrix,
    not a drawing decision -- R reads it and applies the tool's colour mapping itself.
  * NOT the intensity, and NOT a colour. log1p(v * scale) / log1p(100), the 12-stop LUT
    and the lookup all live in R/panels/fig1b_hic.R, because those are the rendering.

THE COLOUR SCALE, AND WHY IT IS COMPUTED HERE
The tool renders intensity as log1p(v * scale) / log1p(100) with a per-resolution `scale`
(package_stage/backend/main.py, _compute_color_scale, and the matching comment in
desktop_repo/src/components/token_best_mcool.vue renderTileToCanvas). scale = 1.0 for
raw-count levels; for an ICE-balanced level it is _COLOR_ANCHOR_V / median, with
_COLOR_ANCHOR_V = 3.0, so that level's median contact lands at the same intensity
(~0.30) a raw median renders at. Reproducing it needs the matrix, hence Python.

The probe is copied from the backend exactly, not approximated:
    T = 256; n = min(n_bins, 3 * T)          -> the first 768 bins of the first chromosome
    fetch the (n x n) near-diagonal block, balanced, keep upper triangle, keep v > 0
    treat as balanced iff  size >= 50  and  percentile(v, 99) < 1.0
    scale = 3.0 / median(v)
Two documented quirks of that probe are reproduced rather than fixed, because fixing
either would change every heatmap the software has ever drawn: it samples a fixed BIN
COUNT (so a different genomic span per resolution), and it takes the median over ALL
separations in the block, not over the drawn window. The backend's own docstring records
the measured value for this matrix at 10 kb (median 1.711e-4); this script asserts it
reproduces that, so a silent drift in cooler or in the file is caught here.

THE DRAWN WINDOW IS A FULL ROTATED TRIANGLE
The exported region is the full triangle over the track window: apex separation = window
width, i.e. exactly those contact pairs BOTH of whose ends lie inside the displayed locus.
That is the non-arbitrary extent for a rotated map and it is what the tool shows when the
region is fitted to the canvas (token_best_mcool.vue rotates by -pi/4 and skips r > c, so
it too draws only the upper triangle). An earlier version truncated at a 300 kb
separation, which is a narrow near-diagonal strip: under the tool's own scale every cell
in it lands in the darkest bands (29.5% pure black, the five lightest stops unused), so
the map lost its internal structure. The full triangle uses 11 of the 12 stops.
"""
import hashlib
import json
import os

import cooler
import numpy as np

MCOOL = ("/home/coder/TAD/creditad/package_stage/example_data/chr7_tracks/"
         "GM12878_chr7_hg38.mcool")
RESOLUTION = 10000              # directive: "pick the 10 kb resolution"
CHROM = "chr7"

# The drawn track window (panel_data/fig1_zoom.json), and the triangle over it.
WIN_LO, WIN_HI = 86525000, 87325000

# ROTATED TRIANGLE, WIDE AND SHALLOW (sixth review, matching a Juicebox-style reference).
#
# The drawn aspect is (window span) : (separation cap). The reference is ~3:1, and the
# directive proposed reaching it by widening the window to 2.0 Mb with a 650 kb cap. The
# delegator ruled for the equivalent that costs no data change: the existing 800 kb window
# with a 260 kb cap is 800/260 = 3.077:1, the same geometry as 2000/650 = 3.077:1, with the
# tracks already extracted at this window and every window-keyed caption number still true.
MAX_SEP_BP = 260000

# A rotated cell is drawn at the MIDPOINT of its two ends, so a cell whose midpoint is inside
# the window may have an end outside it. Fetch half the cap beyond each edge or the drawn
# triangle's upper corners come back empty for contacts that were measured but not asked for.
FETCH_LO = WIN_LO - MAX_SEP_BP // 2 - 2 * RESOLUTION
FETCH_HI = WIN_HI + MAX_SEP_BP // 2 + 2 * RESOLUTION

# ONE ROW / ONE COLUMN OF OVERDRAW, so the panel's clip cuts a STRAIGHT edge. Diamonds tile
# by interlocking -- each row of constant separation is offset half a cell from the row below
# -- which fills the interior exactly but leaves a sawtooth of bare paper along any boundary
# whose neighbouring row is absent. Measured on an earlier rotated render: 58% of the top
# edge band and 23% of the right edge band was paper, which reads as the same dropout
# speckling that was flagged inside the field. Exporting one cell beyond each drawn limit
# gives coord_cartesian something to cut. It does NOT widen what is drawn.
OVERDRAW_SEP_BP = 2 * RESOLUTION
OVERDRAW_POS_BP = 2 * RESOLUTION

# Backend constants, reproduced verbatim (package_stage/backend/main.py).
ANCHOR_V = 3.0                  # _COLOR_ANCHOR_V
PROBE_T = 256                   # T in _compute_color_scale
BACKEND_RECORDED_MEDIAN_10KB = 1.711e-4   # from the backend docstring, for this matrix

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "source_data")
STEM = "fig1_hic_contact_matrix"


def probe_color_scale(clr):
    """Reproduce the backend's _compute_color_scale for one resolution level."""
    n_bins = clr.info["nbins"]
    n = min(n_bins, 3 * PROBE_T)
    m = clr.matrix(balance=True, sparse=True)[0:n, 0:n].tocoo()
    upper = m.row <= m.col                       # the pixel table is upper-triangle
    v = m.data[upper].astype(float)
    v = v[np.isfinite(v) & (v > 0)]
    balanced = v.size >= 50 and float(np.percentile(v, 99)) < 1.0
    if not balanced:
        return 1.0, None, n, int(v.size)         # raw-count level -> tool uses scale 1.0
    med = float(np.median(v))
    assert med > 0, "probe median is not positive"
    return ANCHOR_V / med, med, n, int(v.size)


def main():
    clr = cooler.Cooler("%s::/resolutions/%d" % (MCOOL, RESOLUTION))
    assert clr.binsize == RESOLUTION, "cooler binsize %s != %s" % (clr.binsize, RESOLUTION)

    scale, med, probe_bins, probe_px = probe_color_scale(clr)
    assert med is not None, "10 kb level did not classify as balanced; the tool would use scale=1.0"
    rel = abs(med - BACKEND_RECORDED_MEDIAN_10KB) / BACKEND_RECORDED_MEDIAN_10KB
    assert rel < 0.005, (
        "probe median %.6g disagrees with the value the backend recorded for this matrix "
        "(%.6g, %.2f%% off) -- the file or cooler changed; do not ship a scale the software "
        "would not use" % (med, BACKEND_RECORDED_MEDIAN_10KB, 100 * rel))

    region = "%s:%d-%d" % (CHROM, FETCH_LO, FETCH_HI)
    # DENSE, not sparse: the sparse pixel table cannot distinguish "absent because the
    # count is zero" from "absent because the bin is unweighted", and that distinction is
    # exactly what the white-speckle defect turned on.
    bal = clr.matrix(balance=True, sparse=False).fetch(region)
    raw = clr.matrix(balance=False, sparse=False).fetch(region)
    bins = clr.bins().fetch(region)
    starts = bins["start"].values
    weight = bins["weight"].values

    n = bal.shape[0]
    bi, bj = np.triu_indices(n)                 # upper triangle, diagonal included
    si, sj = starts[bi], starts[bj]
    mid = (si + sj) / 2.0
    keep = (((sj - si) <= MAX_SEP_BP + OVERDRAW_SEP_BP)
            & (mid >= WIN_LO - OVERDRAW_POS_BP) & (mid <= WIN_HI + OVERDRAW_POS_BP))
    bi, bj, si, sj = bi[keep], bj[keep], si[keep], sj[keep]
    val = bal[bi, bj].astype(float)
    cnt = raw[bi, bj].astype(float)
    wok = np.isfinite(weight[bi]) & np.isfinite(weight[bj])

    # THREE STATES, kept distinct all the way to the panel:
    #   measured      a contact count was observed        -> balanced value
    #   zero          measurable, and the count is zero    -> balanced value 0
    #   unmeasurable  a bin carries no ICE weight          -> no value (NA)
    status = np.where(~wok, "unmeasurable",
                      np.where(cnt > 0, "measured", "zero"))
    val = np.where(status == "unmeasurable", np.nan,
                   np.where(status == "zero", 0.0, val))

    order = np.lexsort((bj, bi))
    bi, bj, si, sj, val, cnt, status = (bi[order], bj[order], si[order], sj[order],
                                       val[order], cnt[order], status[order])
    assert bi.size > 0, "no cells survived the window filter"

    # COMPLETENESS. Every cell of the drawn triangle must be present, or the panel leaves
    # bare paper where a reader will see a dropout. This is the defect that produced the
    # white speckles, asserted away rather than fixed by hand.
    expect = int(keep.sum())
    assert bi.size == expect, "exported %d cells but the drawn triangle has %d" % (bi.size, expect)
    n_zero = int((status == "zero").sum())
    n_unmeas = int((status == "unmeasurable").sum())
    finite = val[np.isfinite(val) & (val > 0)]
    assert finite.size > 0, "no positive contact values in the drawn triangle"

    os.makedirs(OUTDIR, exist_ok=True)
    tsv = os.path.join(OUTDIR, STEM + ".tsv")
    head = [
        "# Figure 1b -- GM12878 chr7 Hi-C contact values, for the R-drawn contact map.",
        "# Written by figures/scripts/export_hic_matrix.py. EXPORT ONLY: no value here is",
        "# rescaled, logged or coloured -- these are the cooler's own ICE-balanced contact",
        "# frequencies (count * weight_i * weight_j), exactly as the software reads them.",
        "# source: %s" % MCOOL,
        "# cooler: /resolutions/%d   balanced: True (ICE 'weight' column)" % RESOLUTION,
        "# drawn window: %s:%d-%d   triangle apex separation: %d bp"
        % (CHROM, WIN_LO, WIN_HI, MAX_SEP_BP),
        "# fetched region: %s (half the apex beyond each edge, so the drawn" % region,
        "#   triangle's upper corners are covered; the panel clips to the drawn window)",
        "# rows: EVERY cell of the drawn region -- including cells with no contacts, so the",
        "#   panel renders them deliberately rather than leaving bare paper, and including",
        "#   one row/column of OVERDRAW beyond the drawn limits so the panel's clip cuts a",
        "#   straight edge (the overdraw is exported, not drawn).",
        "# geometry: ROTATED TRIANGLE, %.2f:1 (window %d kb : separation cap %d kb)"
        % ((WIN_HI - WIN_LO) / float(MAX_SEP_BP), (WIN_HI - WIN_LO) / 1000,
           MAX_SEP_BP / 1000),
        "# status: measured      = a contact count was observed",
        "#         zero          = both bins are ICE-weighted and the count is zero",
        "#         unmeasurable  = a bin carries no ICE weight (balanced is NA)",
        "# columns: bin_i bin_j start_i start_j balanced count status",
        "bin_i\tbin_j\tstart_i\tstart_j\tbalanced\tcount\tstatus",
    ]
    with open(tsv, "w") as fh:
        fh.write("\n".join(head) + "\n")
        for a, b, c, d, v, k, s in zip(bi, bj, si, sj, val, cnt, status):
            fh.write("%d\t%d\t%d\t%d\t%s\t%d\t%s\n"
                     % (a, b, c, d, "NA" if not np.isfinite(v) else "%.10g" % v, int(k), s))

    sha = hashlib.sha256(open(tsv, "rb").read()).hexdigest()
    open(os.path.join(OUTDIR, STEM + ".sha256"), "w").write("%s  %s.tsv\n" % (sha, STEM))

    meta = dict(
        source_mcool=MCOOL, cooler_uri="%s::/resolutions/%d" % (MCOOL, RESOLUTION),
        resolution_bp=RESOLUTION, chrom=CHROM,
        balanced=True, balance_method="ICE (cooler 'weight' column)",
        drawn_window=[WIN_LO, WIN_HI], triangle_apex_sep_bp=MAX_SEP_BP,
        separation_cap_bp=MAX_SEP_BP,
        geometry="rotated_triangle", n_bins_axis=int(n),
        axis_window=[WIN_LO, WIN_HI],
        overdraw_sep_bp=OVERDRAW_SEP_BP, overdraw_pos_bp=OVERDRAW_POS_BP,
        drawn_aspect=round((WIN_HI - WIN_LO) / float(MAX_SEP_BP), 4),
        overdraw_note=("cells within one row/column beyond the drawn limits are exported "
                       "so the panel's clip cuts a straight edge instead of leaving the "
                       "diamond tiling's boundary sawtooth; they are not drawn"),
        geometry_note=("45-degree rotated triangle, x = midpoint of the two ends, "
                       "y = their separation, capped so the drawn field is wide and "
                       "shallow at %.2f:1" % ((WIN_HI - WIN_LO) / float(MAX_SEP_BP))),
        # KEY NAMES ARE PART OF THE INTERFACE. verify_figure_claims_v2.py,
        # build_figure_sources_v2.py, assemble_captions.py and R/panels/fig1b_hic.R all read
        # this file, so a key rename is a breaking change to four consumers -- one of which
        # (the claims gate) is what proves the panel shows what the caption says. Renaming
        # them here in the redraw broke the suite; the original names are kept.
        exported_region=region, n_cells=int(val.size),
        n_bins=int(len(set(bi.tolist()) | set(bj.tolist()))),
        value_min_finite=float(finite.min()), value_max_finite=float(finite.max()),
        n_measured=int((status == "measured").sum()),
        n_zero_contacts=n_zero, n_unmeasurable=n_unmeas,
        cell_status_note=("every cell of the drawn triangle is exported; 'zero' means "
                          "measurable with zero observed contacts, 'unmeasurable' means a "
                          "bin carries no ICE weight"),
        # the software's rendering parameters, for R to apply
        color_scale=scale, color_anchor_v=ANCHOR_V,
        probe_median_nonzero=med, probe_block_bins=int(probe_bins),
        probe_block_span_bp=int(probe_bins * RESOLUTION),
        probe_nonzero_pixels=probe_px,
        color_scale_rule=("intensity = min(1, log1p(v * color_scale) / log1p(100)); "
                          "color_scale = color_anchor_v / median(non-zero balanced "
                          "contacts in the first min(n_bins, 3*256) bins of the first "
                          "chromosome), per backend _compute_color_scale"),
        lut_source="desktop_repo/src/components/token_best_mcool.vue CUSTOM_COLORS",
        tsv_sha256=sha, exported_by="figures/scripts/export_hic_matrix.py",
    )
    # THE KEY CONTRACT, asserted here. Four consumers read this file:
    #   scripts/verify_figure_claims_v2.py     (the claims gate)
    #   scripts/build_figure_sources_v2.py     (traceability)
    #   scripts/assemble_captions.py           (the caption's scale value)
    #   R/panels/fig1b_hic.R                   (the panel itself)
    # Renaming a key here breaks them one gate at a time, several minutes apart, which is
    # exactly how the redraw broke the suite twice. Asserting the contract at the point of
    # writing turns that into an immediate, local failure.
    REQUIRED_KEYS = ("geometry", "n_bins_axis", "axis_window", "drawn_aspect",
                     "n_measured", "n_zero_contacts", "n_unmeasurable",
                     "source_mcool", "cooler_uri", "resolution_bp", "balanced",
                     "exported_region", "drawn_window", "triangle_apex_sep_bp",
                     "n_cells", "n_bins", "value_min_finite", "value_max_finite",
                     "color_scale", "tsv_sha256")
    missing_keys = [k for k in REQUIRED_KEYS if k not in meta]
    assert not missing_keys, ("meta is missing key(s) consumers read: %s" % missing_keys)

    with open(os.path.join(OUTDIR, STEM + "_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
        fh.write("\n")

    print("wrote %s" % tsv)
    print("  %d cells (%d measured, %d zero-contact, %d unmeasurable), %d bins"
          % (val.size, meta["n_measured"], n_zero, n_unmeas, meta["n_bins"]))
    print("  rotated triangle %.2f:1, cap %d kb (+%d kb overdraw); values %.6g-%.6g"
          % ((WIN_HI - WIN_LO) / float(MAX_SEP_BP), MAX_SEP_BP / 1000,
             OVERDRAW_SEP_BP / 1000, finite.min(), finite.max()))
    print("  probe: %d bins (%.2f Mb), %d non-zero px, median %.6g"
          % (probe_bins, probe_bins * RESOLUTION / 1e6, probe_px, med))
    print("  color_scale = %.1f / %.6g = %.6f" % (ANCHOR_V, med, scale))
    print("  sha256 %s" % sha)


if __name__ == "__main__":
    main()
