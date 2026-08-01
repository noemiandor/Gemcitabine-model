# Figure 7 print-size visual QC

Status: **PASS** (2026-08-01).

The final composite was inspected directly at its intended full-page size of
7.1 x 9.7 inches. The PNG is 2130 x 2910 pixels at 300 DPI; the vector PDF is a
single 511 x 698 point page (7.10 x 9.69 inches after PDF-point rounding).

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
  (SHA-256 `cff7eeebc2bf0828ebc121753adac0a487e6fb1df48ea8aeac9c8e5c4f11d034`).
- Cairo embeds a creation timestamp, so separately exported PDFs are not
  byte-identical. Rasterizing both vector PDFs independently at 150 DPI with
  `pdftocairo -png -singlefile -r 150` produced byte-identical 1065 x 1455
  renders (SHA-256
  `c253198bd7d214357055c53716d528a45d6540e3e220fabdf1b3f211d39c4104`).
