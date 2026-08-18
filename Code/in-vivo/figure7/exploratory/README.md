# Exploratory Figure 7I: treated-versus-vehicle expression within the state interval

This directory contains a candidate replacement analysis for main Figure 7I.
It is intentionally isolated from the reviewed Figure 7 workflow and does not
modify the canonical composite or frozen panel-7F reference.

## Question and relation to the existing panel

The current Figure 7I source panel does **not** test 2N versus 4N expression.
It fits pseudotime-dependent expression across the full CellCycle trajectory
and contrasts the accumulated-state interval with its two flanks. Dose and
injected origin are adjustment terms in that model.

The exploratory analysis here asks a different question: among tumor cells in
the exact inclusive pseudotime interval 0.296--0.486, which human-tumor genes
and pathways differ between gemcitabine-treated and vehicle tumors?

## Statistical design

- The input cell universe is the reviewed 2,881-cell CellCycle table.
- Only exact `GRCh38-` RNA features are retained. Exact `GRCm39-` features are
  removed before expression filtering, normalization, and modeling;
  unrecognized prefixes fail closed.
- The interval contains 417 cells from all 16 tumors.
- Counts are summed once per mouse within the interval. The 16 mouse
  pseudobulks, not the 417 cells, are the independent observations.
- Genes must have CPM > 1 in at least four mice, matching the size of the
  smallest dose group.
- The model uses edgeR TMM normalization and limma-voom with robust empirical
  Bayes moderation. It adjusts for injected 2N/4N origin and each mouse's mean
  pseudotime within the interval.
- The primary contrast is
  `0.5 * (30 mg/kg) + 0.5 * (120 mg/kg) - vehicle`. Because both treated dose
  groups contain four mice, this is also an equal-treated-mouse average.
  Dose-specific contrasts are written as diagnostics.
- A positive log2 fold change means higher expression in treated tumors.
- The same equal-dose comparison is also fit separately within the 2N-origin
  and 4N-origin tumors. Each stratum contains eight mice: four vehicle, two at
  30 mg/kg, and two at 120 mg/kg. Because injected origin is constant within a
  stratum, these models omit that term but retain adjustment for mean
  within-interval pseudotime. Their expression filter is CPM > 1 in at least
  two mice, matching the smallest within-origin dose group.
- A formal treatment-by-origin model is fit to all 16 mouse pseudobulks using
  six origin-by-dose group means plus mean within-interval pseudotime. Its
  interaction contrast is
  `[equal-dose treated - vehicle]_4N - [equal-dose treated - vehicle]_2N`.
  Positive interaction effects are more positive in 4N; negative effects are
  more positive in 2N. "More positive" does not necessarily mean a larger
  absolute response: it can also represent weaker depletion.
- The 2N and 4N treatment contrasts displayed for interaction-selected
  pathways are estimated from that same joint model and expressed-gene
  universe. This keeps the paired estimates directly aligned with the formal
  interaction test.
- Counts of FDR-significant pathways are retained as descriptive audit values,
  not compared as a measure of response strength. Such counts depend on
  standard errors, gene-set overlap, separate expression filters, and the
  significance threshold as well as effect magnitude.
- A separate threshold-free comparison directly estimates global response
  magnitude on the joint model's single expressed-gene universe. For each
  origin it averages the squared equal-dose treated-minus-vehicle log2 fold
  changes across genes and subtracts the corresponding moderated sampling
  variance. The formal estimand is the resulting 2N-minus-4N difference in
  noise-corrected mean squared log2 fold change. A positive value denotes a
  larger total transcriptional displacement in 2N; it does not imply that the
  same genes or pathways move in the same direction.
- Uncertainty for that magnitude difference uses an HC2 whole-mouse
  Rademacher wild bootstrap with fixed voom weights. Every gene within a mouse
  receives the same bootstrap multiplier, preserving cross-gene dependence;
  the default is 4,999 replicates. This is an exploratory, small-sample global
  test and is reported with its interval-selected-analysis caveat.
- For a pathway-level Figure 7I candidate, the primary model's moderated
  gene-level t statistics are passed to the same GSEA implementation used by
  the existing Figure 7I workflow. The analysis uses the pinned Homo sapiens
  MSigDB Hallmark, Reactome, and GO biological-process collections, adaptive
  `fgseaMultilevel` settings, collection-wide BH correction, and display rule
  (FDR <= 0.05; up to three pathways per direction and collection; no
  nonsignificant backfill) from the reviewed workflow. The pooled, 2N-origin,
  and 4N-origin panels all use this same rule. Positive normalized enrichment
  scores indicate enrichment in treated tumors.
- The Figure 7I interaction panel contains only GSEA of the formal interaction.
  It displays the 10 most negative and 10 most positive interaction NES across
  all three collections, ranked by NES within sign. A negative interaction NES
  denotes a more positive treated-minus-vehicle response in 2N-origin tumors;
  a positive interaction NES denotes a more positive response in 4N-origin
  tumors. Each pathway is represented once by its interaction NES, collection
  is encoded by color, and point size represents the interaction FDR. The
  separate-origin results remain audit outputs but are not displayed and are
  not used as evidence of an interaction.
- A sensitivity model omits only the mean within-interval pseudotime
  adjustment; injected origin remains an adjustment term and the gene universe
  and treatment contrast remain unchanged.

The interval itself was localized from a treated-versus-vehicle density
comparison. Consequently, these differential-expression results are
exploratory and conditional on that data-selected interval; they should not be
presented as independent confirmatory treatment tests. The origin-stratified
analyses are additionally lower-powered because each model contains only eight
independent mouse pseudobulks. Claims of origin-dependent treatment response
must be based on the direct interaction result, not on differing within-origin
significance calls.

## Run

With the reviewed raw Seurat object already present at its standard location:

```bash
Rscript Code/in-vivo/figure7/exploratory/treated_vs_vehicle_interval_de.R
```

To run only the threshold-free global magnitude audit without loading a
pathway database or regenerating pathway panels:

```bash
Rscript Code/in-vivo/figure7/exploratory/treated_vs_vehicle_interval_de.R \
  --global-magnitude-only=true
```

The isolated default output is:

```text
Results/in-vivo/figure7/exploratory/treated_vs_vehicle_interval_de/
```

It contains pooled, 2N-origin, and 4N-origin pathway-enrichment candidate
PDF/PNG files; a shared-plus-formal-interaction comparison PDF/PNG; the
pooled gene volcano as an audit figure; complete primary and dose-specific DE
tables; the joint interaction design and exact contrasts; complete interaction
DE and GSEA tables; paired joint-model origin effects; complete and selected
pooled/stratified GSEA tables; leading-edge genes; the exact
gene-set contract and membership; mouse coverage; design matrices and contrast
definitions; species and expression-filter audits; portable input provenance;
global response-magnitude estimates, their whole-mouse bootstrap audit;
package versions; and session information.

Run the focused synthetic tests from the repository root with:

```bash
Rscript -e 'testthat::test_file("Code/in-vivo/figure7/exploratory/tests/test_treated_vs_vehicle_interval_de.R")'
```
