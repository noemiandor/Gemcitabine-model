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
are differentially expressed in gemcitabine-treated versus vehicle tumors?

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
- A sensitivity model omits only the mean within-interval pseudotime
  adjustment; injected origin remains an adjustment term and the gene universe
  and treatment contrast remain unchanged.

The interval itself was localized from a treated-versus-vehicle density
comparison. Consequently, this differential-expression result is exploratory
and conditional on that data-selected interval; it should not be presented as
an independent confirmatory treatment test.

## Run

With the reviewed raw Seurat object already present at its standard location:

```bash
Rscript Code/in-vivo/figure7/exploratory/treated_vs_vehicle_interval_de.R
```

The isolated default output is:

```text
Results/in-vivo/figure7/exploratory/treated_vs_vehicle_interval_de/
```

It contains the candidate PDF/PNG, complete primary and dose-specific DE
tables, mouse coverage, the design matrix and contrast definitions, species
and expression-filter audits, portable input provenance, package versions,
and session information.

Run the focused synthetic tests from the repository root with:

```bash
Rscript -e 'testthat::test_file("Code/in-vivo/figure7/exploratory/tests/test_treated_vs_vehicle_interval_de.R")'
```
