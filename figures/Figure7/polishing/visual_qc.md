# Figure 7 print-size visual QC

Status: **PASS** (2026-08-02).

The final composite and companion Supplementary Figures 4 and 6 were inspected
at their exported sizes. Figure 7 is 7.1 x 10.645 inches and 2130 x 3193 pixels
at 300 DPI.

- Panel identity and order: A--L match `panel_map.csv` and the enforced
  A=7A; B=7C; C=SI4A; D=SI4B; E=SI4C; F=SI4E; G=SI7B; H=7B-ECDF; I=7F;
  J=QC-filtered copy-number heatmap; K=7D; L=7E contract.
- C--E: all three UMAPs are distinct at print size; C retains nine noncolliding
  cluster labels, including separate 4c and 10 labels.
- G and I: pathway labels, dendrograms, interval markers, and color scales are
  legible without clipping.
- H: the overall, 4N-origin, and 2N-origin vehicle-versus-gemcitabine ECDF
  curves are separately discernible on their common scale. No compressed lower
  strip remains.
- J: the heatmap contains 9,832 rows and 22 chromosome columns. Injected-origin,
  dose, and mouse annotation bars are visible; sparse chromosome labels
  (1, 5, 9, 13, 17, 22) do not collide; the copy-number color scale remains
  visible. The caption, rather than a space-consuming 16-mouse key, defines the
  annotation bars.
- F: positive-enrichment asterisks remain visible above the correct bars.
- K and L: all eight mouse labels and both statistics boxes are legible and do
  not obscure points.
- Supplementary Figure 4I: the pointwise-supported 0.296--0.486 region,
  family-wise-supported 0.414--0.426 region, raw-excess peak at 0.452, and
  strongest standardized evidence at 0.420 remain distinguishable.
- Supplementary Figure 6: the former copy-number heatmap is absent and the
  remaining A--E layout has no empty labeled slot.
- Readability, typography, and clipping: panel letters, axes, tick labels, pathway labels,
  color bars, and data marks remain within the devices.

The final hash comparison against the independently generated Manager output is
recorded in `figure_byte_identity_report.tsv`.
