# Figure 7 print-size visual QC

Status: **PASS** (2026-08-01).

The final composite was inspected directly at its intended full-page size of
7.1 x 10.645 inches. The PNG is 2130 x 3193 pixels at 300 DPI; the vector PDF is
a single 511 x 766 point page (7.10 x 10.64 inches after PDF-point rounding).

- Panel identity and order: A--K match `panel_map.csv` and the enforced
  A=7A; B=7C; C=SI4A; D=SI4B; E=SI4C; F=SI4E; G=SI7B; H=7B; I=7F; J=7D;
  K=7E contract.
- Titles and prose: no panel-internal title, subtitle, caption, formula, or
  methods paragraph remains. Panel letters are uniform, unobstructed, and in
  row-major reading order.
- Clipping: axes, tick labels, pathway labels, color bars, dendrograms, panel
  letters, and data marks remain inside the device. The shortened panel-A y
  label and outer gutter prevent the former rotated-label clipping.
- UMAP labels: C displays all nine cluster identifiers without collisions;
  in particular, 4c and 10 remain separate and legible.
- Heatmaps: G and I have dedicated enlarged rows for labels. G retains both
  row and column dendrograms and an explicit NES color scale; I retains its
  collection strips, interval markers, and mean-gene-z-score color scale.
- Pseudotime redistribution: H is 2.295 inches high at the target size. The
  overall, 4N-origin, and 2N-origin vehicle-versus-gemcitabine ECDF curves are
  separately discernible on their unchanged common 0--1.05 scale. The lower
  strip clearly shows the equal-mouse density contrast, gray simultaneous null
  envelope, pale pointwise-support region (0.296--0.486), dark family-wise
  support (0.414--0.426), raw-excess peak at 0.452, and strongest standardized
  evidence at 0.420 without implying that the entire rounded 0.30--0.49 state
  window has simultaneous support.
- Composition: F's positive-enrichment asterisks are bold and visible above
  the correct bars at print size.
- Mouse-level panels: J and K each display all eight mouse identifiers. J's
  statistics box is in the unoccupied upper-left corner; K's remains in the
  upper-right. Neither box obscures a point or identifier.
- Readability, typography, and spacing: reader-facing axis/facet/pathway labels, consistent
  sans-serif type, margins, gutters, and color-bar typography are legible at
  the target size. Compatible graphical keys are described once in the
  manuscript legend rather than repeated inside small panels.

Independent production comparison:

- The polishing PNG and the Figure-7-only Manager PNG are byte-identical
  (SHA-256 `eb2b522b74b49c391a98624463bb32c3e18fa2ce18a4f2ed4e035903c2b83e03`).
- Cairo embeds a creation timestamp, so separately exported PDFs are not
  byte-identical. Rasterizing both vector PDFs independently at 150 DPI with
  `pdftocairo -png -singlefile -r 150` produced byte-identical 1065 x 1596
  renders (SHA-256
  `31a7c4a67950388bbb993a60852b7b6a79796ad81e2dd6cdf7bfa4359f4fea58`).
