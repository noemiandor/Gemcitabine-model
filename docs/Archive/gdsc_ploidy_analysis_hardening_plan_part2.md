## Revised Codex implementation plan: harden GDSC–ploidy analysis without prematurely changing canonical outputs

### Guiding principle

Do **not** change the canonical manuscript-facing result yet. Keep the current `Z_SCORE`-based pipeline and current enrichment logic reproducible. Add audits, full correlation statistics, and supplemental multi-metric outputs so that output-changing decisions can be reviewed explicitly.

The implementation should proceed in two phases:

1. **Phase 1: audit + supplemental outputs only**

   * Preserve current canonical `Z_SCORE` pipeline.
   * Add cell-line matching collision audits.
   * Keep current raw-name matching for canonical outputs unless a regression report explicitly accepts a normalized-key output change.
   * Add normalized-key matching as supplemental audit/statistics first.
   * Add full correlation statistics for canonical `Z_SCORE`.
   * Add supplemental multi-metric outputs for `Z_SCORE`, `LN_IC50`, and `AUC`.
   * Add drug annotation audit tables.
   * Do **not** yet replace PubChem/legacy categories with manual overrides for enrichment.

2. **Phase 2: curated category switch**

   * Only after reviewing the annotation audit, create an approved curated drug-class mapping.
   * Switch enrichment to `final_category` only when that curated file exists and has explicit approval/status columns.

---

# 0. Dependency policy for this pass

This plan intentionally uses `data.table` for tabular joins, grouped summaries, `fwrite()`, `:=`, `uniqueN()`, and `frank()`. This supersedes the earlier "do not introduce new dependencies beyond `openxlsx`" constraint: add `data.table` as an explicit dependency in `dependencies.R` and README dependency documentation before using it in the analysis.

Required dependency changes:

```r
require_packages(c("EnrichIntersect", "xlsx", "plyr", "openxlsx", "data.table"))
```

Do not use invalid base-pipe/data.table placeholder syntax such as:

```r
dt |>
  .[, .N, by = group]
```

Use one of these valid styles instead:

```r
dt <- data.table::as.data.table(dt)
dt[, .N, by = group]
```

or, when avoiding attachment:

```r
dt <- data.table::as.data.table(dt)
dt[, list(n = .N), by = "group"]
```

If adding `data.table` is not acceptable, rewrite the examples in this document using base R before implementation. Do not mix `data.table` syntax into production code without an explicit dependency check.

Use the active script layout:

```text
Code/gdsc_ploidy_analysis/src/
Code/gdsc_ploidy_analysis/data/raw/
Code/gdsc_ploidy_analysis/data/derived/
Code/gdsc_ploidy_analysis/data/manual/
Code/gdsc_ploidy_analysis/output/metadata/
Code/gdsc_ploidy_analysis/output/tables/
```

Do not introduce a new `output/logs/session_info.txt` path unless the current metadata layout is intentionally migrated. Runtime metadata should continue to go under `output/metadata/`.

---

# 1. Cell-line normalization: add collision audit before matching

## Problem

A normalized cell-line key can map multiple original names to the same key. Therefore, do **not** simply run:

```r
rownames(appCL) <- appCL$CELL_LINE_KEY
```

without collision handling.

## Required implementation

Add this helper in `analysis_helpers.R`:

```r
normalize_cell_line_name <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9]", "", x)
  x
}
```

Add a collision-audit helper:

```r
audit_normalized_key_collisions <- function(dt, raw_col, key_col, source_name) {
  tmp <- unique(as.data.frame(dt)[, c(raw_col, key_col), drop = FALSE])
  names(tmp) <- c("raw_name", "normalized_key")

  out <- data.table::as.data.table(tmp)
  out <- out[, .(
    n_raw_names = data.table::uniqueN(raw_name),
    raw_names = paste(sort(unique(raw_name)), collapse = " | ")
  ), by = normalized_key]
  out <- out[n_raw_names > 1]

  out[, source := source_name]
  out[]
}
```

Apply to both GDSC and ploidy tables:

```r
dr$CELL_LINE_KEY <- normalize_cell_line_name(dr$CELL_LINE_NAME)
appCL$CELL_LINE_KEY <- normalize_cell_line_name(appCL$`Cell iname`)
```

Write collision tables:

```r
gdsc_collisions <- audit_normalized_key_collisions(
  dr, raw_col = "CELL_LINE_NAME", key_col = "CELL_LINE_KEY", source_name = "GDSC"
)

ploidy_collisions <- audit_normalized_key_collisions(
  appCL, raw_col = "Cell iname", key_col = "CELL_LINE_KEY", source_name = "CellPassports_ploidy"
)

data.table::fwrite(gdsc_collisions, file.path(output_dir, "tables", "cell_line_key_collisions_gdsc.tsv"), sep = "\t")
data.table::fwrite(ploidy_collisions, file.path(output_dir, "tables", "cell_line_key_collisions_ploidy.tsv"), sep = "\t")
```

## Collision handling rule

For **ploidy data**, fail if a normalized key maps to multiple rows with materially different ploidy values.

Suggested rule:

```r
ploidy_collision_tolerance <- 0.05

appCL_dt <- data.table::as.data.table(appCL)

ploidy_key_summary <- appCL_dt[, .(
  n_rows = .N,
  n_raw_names = data.table::uniqueN(`Cell iname`),
  ploidy_min = min(ploidy, na.rm = TRUE),
  ploidy_max = max(ploidy, na.rm = TRUE),
  ploidy_range = max(ploidy, na.rm = TRUE) - min(ploidy, na.rm = TRUE),
  raw_names = paste(sort(unique(`Cell iname`)), collapse = " | ")
), by = CELL_LINE_KEY]

bad_ploidy_collisions <- ploidy_key_summary[n_raw_names > 1 & ploidy_range > ploidy_collision_tolerance]

if (nrow(bad_ploidy_collisions) > 0) {
  data.table::fwrite(
    bad_ploidy_collisions,
    file.path(output_dir, "tables", "cell_line_key_collisions_ploidy_FAIL.tsv"),
    sep = "\t"
  )
  stop("Normalized ploidy cell-line keys are not unique and have discordant ploidy values. Review cell_line_key_collisions_ploidy_FAIL.tsv.")
}
```

If duplicates have identical or near-identical ploidy, collapse them explicitly and write a record of the collapse.

Do **not** use rownames until normalized keys are confirmed unique.

Preferred approach for supplemental normalized-key analyses: use explicit joins by `CELL_LINE_KEY`, not rownames.

## Canonical-output protection

Changing the matching key can change the canonical `Z_SCORE` output by adding or removing matched cell lines. Therefore, Phase 1 must be audit-first:

1. Keep the existing canonical analysis matching behavior for manuscript-facing legacy outputs.
2. Write normalized-key matching QC tables and collision tables.
3. Generate normalized-key correlation statistics as clearly named supplemental outputs, for example:

   ```text
   output/tables/drug_ploidy_correlations_by_cancer_Z_SCORE_normalized_key_SUPPLEMENTAL.tsv
   output/tables/drug_ploidy_correlations_by_cancer_Z_SCORE_normalized_key_SUPPLEMENTAL.xlsx
   ```

4. Compare canonical raw-name matching against supplemental normalized-key matching and write:

   ```text
   output/tables/cell_line_matching_delta_raw_vs_normalized.tsv
   output/tables/correlation_delta_raw_vs_normalized_Z_SCORE.tsv
   ```

5. Do not switch canonical outputs to normalized-key matching unless a milestone regression report explicitly accepts the result change and identifies the first changed cell lines, correlations, and enrichment outputs.

Log `ploidy_collision_tolerance` in `output/metadata/run_parameters.tsv`.

---

# 2. Preserve canonical `Z_SCORE` pipeline

## Problem

`LN_IC50` and `AUC` are present, so adding multi-metric analysis will materially expand and potentially change interpretation. Do not silently replace the current analysis.

## Required implementation

Set the current canonical metric explicitly:

```r
canonical_metric <- "Z_SCORE"
supplemental_metrics <- intersect(c("Z_SCORE", "LN_IC50", "AUC"), names(dr))
```

All current manuscript-facing enrichment outputs should continue to use `canonical_metric`.

Add a run-parameter output. Prefer appending to the existing `output/metadata/run_config.tsv` convention if that file remains the canonical metadata table; otherwise write the same fields to `output/metadata/run_parameters.tsv`.

```r
run_params <- data.frame(
  parameter = c(
    "canonical_metric",
    "supplemental_metrics",
    "canonical_matching_mode",
    "supplemental_matching_mode",
    "ploidy_collision_tolerance"
  ),
  value = c(
    canonical_metric,
    paste(supplemental_metrics, collapse = ";"),
    "legacy_raw_name_matching",
    "normalized_cell_line_key",
    as.character(ploidy_collision_tolerance)
  )
)

data.table::fwrite(
  run_params,
  file.path(output_dir, "metadata", "run_parameters.tsv"),
  sep = "\t"
)
```

## Output convention

Canonical outputs should be named clearly:

```text
drug_ploidy_correlations_by_cancer_Z_SCORE.xlsx
drug_ploidy_correlations_by_cancer_Z_SCORE.tsv
```

Supplemental outputs should be named:

```text
drug_ploidy_correlations_by_cancer_all_metrics.xlsx
drug_ploidy_correlations_by_cancer_all_metrics.tsv
```

Do **not** change manuscript figures from `Z_SCORE` to `LN_IC50` or `AUC` until directionality is reviewed.

---

# 3. Add full correlation statistics for canonical Z_SCORE and supplemental metrics

## Required outputs

For the canonical `Z_SCORE` analysis, save:

```text
output/tables/drug_ploidy_correlations_by_cancer_Z_SCORE.xlsx
```

This workbook must have **one sheet per cancer type**.

Each sheet should include:

```text
cancer_type
drug
metric
n
pearson_r
pearson_p
pearson_ci_low
pearson_ci_high
pearson_fdr
spearman_rho
spearman_p
spearman_ci_low
spearman_ci_high
spearman_ci_method
spearman_fdr
rank_by_pearson_desc
rank_by_pearson_asc
rank_by_spearman_desc
rank_by_spearman_asc
response_direction_note
```

Also save the same table as TSV.

For supplemental metrics, save:

```text
output/tables/drug_ploidy_correlations_by_cancer_all_metrics.xlsx
output/tables/drug_ploidy_correlations_by_cancer_all_metrics.tsv
```

Each supplemental cancer-type sheet should contain all available metrics with `metric` as a column.

## Correlation helper

Add to `analysis_helpers.R`:

```r
compute_drug_ploidy_correlation <- function(response, ploidy, min_n = 10) {
  ok <- is.finite(response) & is.finite(ploidy)
  n <- sum(ok)

  if (n < min_n) {
    return(data.frame(
      n = n,
      pearson_r = NA_real_,
      pearson_p = NA_real_,
      pearson_ci_low = NA_real_,
      pearson_ci_high = NA_real_,
      spearman_rho = NA_real_,
      spearman_p = NA_real_,
    spearman_ci_low = NA_real_,
    spearman_ci_high = NA_real_,
    spearman_ci_method = "not_estimated_n_below_minimum"
    ))
  }

  pearson <- suppressWarnings(cor.test(response[ok], ploidy[ok], method = "pearson"))
  spearman <- suppressWarnings(cor.test(response[ok], ploidy[ok], method = "spearman", exact = FALSE))

  pearson_ci <- pearson$conf.int

  rho <- unname(spearman$estimate)
  if (is.finite(rho) && abs(rho) < 1 && n > 3) {
    z <- atanh(rho)
    se <- 1 / sqrt(n - 3)
    spearman_ci <- tanh(c(z - 1.96 * se, z + 1.96 * se))
    spearman_ci_method <- "approximate_fisher_transform"
  } else {
    spearman_ci <- c(NA_real_, NA_real_)
    spearman_ci_method <- "not_estimated"
  }

  data.frame(
    n = n,
    pearson_r = unname(pearson$estimate),
    pearson_p = pearson$p.value,
    pearson_ci_low = pearson_ci[1],
    pearson_ci_high = pearson_ci[2],
    spearman_rho = rho,
    spearman_p = spearman$p.value,
    spearman_ci_low = spearman_ci[1],
    spearman_ci_high = spearman_ci[2],
    spearman_ci_method = spearman_ci_method
  )
}
```

The Spearman confidence interval is an approximation, especially in the presence of ties. Every output containing `spearman_ci_low` and `spearman_ci_high` must also include `spearman_ci_method`, and run metadata must state that `approximate_fisher_transform` is not an exact Spearman confidence interval. Do not describe these columns generically as exact 95% Spearman confidence intervals. If approximate intervals are considered too easy to misuse during implementation, omit them or replace them with fixed-seed bootstrap intervals in a separately named output.

Use BH FDR within each cancer type and metric:

```r
res <- data.table::as.data.table(res)
res[, pearson_fdr := p.adjust(pearson_p, method = "BH"), by = .(cancer_type, metric)]
res[, spearman_fdr := p.adjust(spearman_p, method = "BH"), by = .(cancer_type, metric)]
```

Use neutral ranking names until response directionality is settled:

```r
res[, rank_by_pearson_desc := data.table::frank(-pearson_r, ties.method = "min", na.last = "keep"), by = .(cancer_type, metric)]
res[, rank_by_pearson_asc := data.table::frank(pearson_r, ties.method = "min", na.last = "keep"), by = .(cancer_type, metric)]
res[, rank_by_spearman_desc := data.table::frank(-spearman_rho, ties.method = "min", na.last = "keep"), by = .(cancer_type, metric)]
res[, rank_by_spearman_asc := data.table::frank(spearman_rho, ties.method = "min", na.last = "keep"), by = .(cancer_type, metric)]
```

## Excel writing

Use `openxlsx`, not `xlsx`.

```r
safe_sheet_name <- function(x) {
  x <- gsub("[\\[\\]\\:\\*\\?\\/\\\\]", "_", as.character(x))
  substr(x, 1, 31)
}
```

```r
write_correlations_xlsx <- function(cor_dt, path) {
  wb <- openxlsx::createWorkbook()

  for (ct in sort(unique(cor_dt$cancer_type))) {
    sheet <- safe_sheet_name(ct)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, cor_dt[cancer_type == ct])
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = 1:ncol(cor_dt), widths = "auto")
  }

  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}
```

---

# 4. Drug annotation audit first; do not yet switch enrichment categories

## Problem

`custom_set_candidate.tsv` was reconstructed as a missing-category patch, not as an authoritative curated category table. Applying it universally could change many category assignments.

## Required implementation

Create an audit table, but do **not** yet change the enrichment analysis to use manual overrides universally.

Write:

```text
output/tables/drug_annotation_audit.tsv
```

Columns:

```text
drug
drug_key
pubchem_category
legacy_group_used_for_enrichment
custom_set_candidate_category
would_manual_override_change_legacy_group
would_manual_override_change_pubchem_group
all_pubchem_categories
annotation_source
annotation_status
notes
```

Where:

```r
would_manual_override_change_legacy_group =
  !is.na(custom_set_candidate_category) &
  custom_set_candidate_category != legacy_group_used_for_enrichment
```

Use this audit table to decide later whether `custom_set_candidate.tsv` should become an authoritative curation table.

---

# 5. Do not remove “last listed category” yet; label it legacy

## Problem

The current pipeline may use the last listed category from a multi-category annotation. That is not ideal, but removing it without a curated replacement makes the enrichment categories under-specified.

## Required implementation

Keep the current category used for enrichment as:

```text
legacy_group_used_for_enrichment
```

Do not call it `final_category`.

If a drug has multiple categories, keep them in:

```text
all_pubchem_categories
```

Add a flag:

```text
multi_category_flag
```

Suggested logic:

```r
annotation_audit[, multi_category_flag := grepl(";", all_pubchem_categories, fixed = TRUE)]
```

For now, enrichment should continue to use `legacy_group_used_for_enrichment`.

Add a warning:

```r
warning("Enrichment still uses legacy_group_used_for_enrichment. Review drug_annotation_audit.tsv before switching to curated final_category.")
```

---

# 6. Define the future curated mapping schema, but do not require it yet

Create a template file only:

```text
output/tables/drug_class_final_curated_TEMPLATE.tsv
```

Columns:

```text
drug
drug_key
legacy_group_used_for_enrichment
all_pubchem_categories
custom_set_candidate_category
proposed_final_category
approved_final_category
curation_status
curator
curation_notes
```

Allowed `curation_status` values:

```text
unreviewed
approved
needs_review
exclude_from_enrichment
```

Phase 2 should only switch enrichment to `approved_final_category` when:

1. The curated file exists.
2. Every drug used in enrichment has `curation_status == "approved"` or `"exclude_from_enrichment"`.
3. There are no blank `approved_final_category` values among included drugs.

Until then, keep legacy enrichment unchanged.

---

# 7. Enrichment outputs: preserve legacy, add optional hardened version later

For now, keep current enrichment output as legacy:

```text
output/tables/drugsVsPloidyCorr_legacy_Z_SCORE.xlsx
```

Use explicit enrichment run modes instead of hard-coding a single expensive permutation count:

```r
analysis_mode <- "dev" # allowed: "dev", "manuscript"
enrichment_permute_n <- if (analysis_mode == "manuscript") 10000L else 300L
```

Recommended defaults:

```text
dev:        enrichment_permute_n = 300
manuscript: enrichment_permute_n = 10000
```

The `dev` mode is for local regression and smoke testing. The `manuscript` mode is for final tables only and should be run deliberately because it can be much slower, especially if repeated across metrics. Log `analysis_mode`, `enrichment_permute_n`, selected-drug thresholds, and metric set in `output/metadata/run_parameters.tsv`.

If enrichment code is touched, make only non-controversial improvements:

* document the p-value cutoff;
* document the permutation count;
* save the set of selected drugs used for each enrichment test;
* save enrichment metadata.

Do not yet:

* switch to `final_category`;
* change thresholds;
* change the canonical permutation count;
* change selected-drug definitions;

unless the goal is to create a clearly named supplemental/hardened output, not replace the canonical one.

Recommended new outputs:

```text
output/tables/class_enrichment_legacy_Z_SCORE.tsv
output/tables/class_enrichment_legacy_metadata.tsv
output/tables/class_enrichment_selected_drugs_legacy_Z_SCORE.tsv
```

If a hardened enrichment is added later, name it separately:

```text
output/tables/class_enrichment_curated_Z_SCORE.tsv
```

---

# 8. Add gemcitabine rank summary without changing main conclusions

Create:

```text
output/tables/gemcitabine_rank_summary_Z_SCORE.tsv
output/tables/gemcitabine_rank_summary_all_metrics.tsv
```

Use normalized drug names and aliases:

```r
gem_aliases <- c("GEMCITABINE", "GEMZAR")
```

Report:

```text
cancer_type
metric
drug
n
pearson_r
pearson_p
pearson_fdr
pearson_ci_low
pearson_ci_high
spearman_rho
spearman_p
spearman_fdr
spearman_ci_low
spearman_ci_high
spearman_ci_method
rank_by_pearson_desc
rank_by_pearson_asc
rank_by_spearman_desc
rank_by_spearman_asc
response_direction_note
```

Do not decide in code whether gemcitabine is “low-ploidy selective.” Use neutral rank fields and let the manuscript interpret after confirming metric direction.

Do not name the table, columns, or notes as if they prove low-ploidy selectivity. Acceptable names include `gemcitabine_rank_summary_*`, `gemcitabine_ploidy_association_*`, or `gemcitabine_correlation_summary_*`. Avoid names such as `gemcitabine_low_ploidy_support_*` in code or generated outputs.

---

# 9. Add tissue-adjusted all-cancer model as supplemental only

Do not replace the existing `allcancers` result. Add a supplemental tissue-adjusted model:

```r
lm(response ~ ploidy + TCGA_DESC)
```

Use explicit fit eligibility thresholds:

```r
tissue_model_min_n <- 20L
tissue_model_min_tissues <- 3L
tissue_model_min_rows_per_tissue <- 2L
```

For each drug and metric:

1. Drop rows with non-finite response or ploidy and missing `TCGA_DESC`.
2. Collapse unused tissue levels before fitting.
3. Skip the model with `model_status = "skipped_insufficient_n"` if `n < tissue_model_min_n`.
4. Skip the model with `model_status = "skipped_insufficient_tissue_count"` if fewer than `tissue_model_min_tissues` tissues remain after filtering.
5. Skip or flag sparse tissues with `model_status = "skipped_sparse_tissue_levels"` if any retained tissue has fewer than `tissue_model_min_rows_per_tissue` rows, unless the implementation explicitly collapses rare tissues into an `OTHER` group and records that decision.
6. Fit with `tryCatch()`, not `try()`, and record the error message if fitting fails.
7. After fitting, check rank deficiency. If `fit$rank < length(stats::coef(fit))` or the ploidy coefficient is absent/`NA`, record `model_status = "rank_deficient_no_ploidy_coef"` and do not report a ploidy p-value as interpretable.
8. If the fit is rank-deficient but the ploidy coefficient is estimable, record `model_status = "rank_deficient_ploidy_estimable"` and keep the ploidy coefficient with a note.
9. Only use `model_status = "ok"` when the model is full-rank or the implementation has explicitly determined the ploidy coefficient is estimable.

Output:

```text
output/tables/all_cancers_tissue_adjusted_ploidy_models_Z_SCORE.tsv
output/tables/all_cancers_tissue_adjusted_ploidy_models_all_metrics.tsv
```

Columns:

```text
drug
metric
n
n_cancer_types
n_tissues
beta_ploidy
se_ploidy
t_ploidy
p_ploidy
fdr_ploidy
model_formula
model_status
rank_deficient
skipped_reason
notes
```

Log `tissue_model_min_n`, `tissue_model_min_tissues`, and `tissue_model_min_rows_per_tissue` in `output/metadata/run_parameters.tsv`. This should be described as a confounding check, not the primary result.

---

# 10. PubChem refresh and duplicate-selection safety fixes

These are clear correctness fixes and should be implemented before optional scientific switches.

## PubChem subset refresh safety

The refresh entrypoint under `Code/gdsc_ploidy_analysis/src/` must not accidentally replace the full annotation cache with only the refreshed subset.

Required behavior:

1. Default refresh mode reads the existing cache, refreshes only requested drugs, merges refreshed rows back into the full cache, and writes the full merged cache.
2. A subset-only output is allowed only when the user provides an explicit argument such as `--subset-output`.
3. The refresh should fail if the requested subset contains drugs that cannot be resolved and no existing cached annotation is available, unless an explicit `--allow-missing` style flag is supplied.
4. The refresh should write a metadata table under `output/metadata/` recording requested drug count, refreshed drug count, preserved cached row count, failed drug count, and whether subset-only mode was used.

Regression test after implementation:

```text
1. Start from a full cached annotation table.
2. Refresh a small subset, for example 3 drugs.
3. Assert the default output still contains all original cached drugs plus refreshed values for the requested subset.
4. Assert subset-only output is produced only when the explicit subset-output option is supplied.
```

## Duplicate-selection safety

The duplicate handling strategy `lowest_rmse` must require an `RMSE` column and must not silently fall back to first-row selection.

Required behavior:

1. If `duplicate_strategy == "lowest_rmse"` and `RMSE` is absent, stop with a clear error.
2. If `duplicate_strategy == "first_after_sort"` is used, document the sort columns and write duplicate-selection metadata.
3. Always write a duplicate audit table showing all duplicate candidate rows, the selected row, and the selection reason.
4. Log `duplicate_strategy`, required columns, duplicate counts, and selected-row counts in `output/metadata/run_parameters.tsv` or a companion metadata table.

Regression test after implementation:

```text
1. Run with current canonical duplicate strategy and confirm canonical outputs are unchanged.
2. Run a fixture with duplicate rows and no RMSE while requesting lowest_rmse.
3. Assert the run fails with a clear missing-RMSE error.
4. Run a fixture with duplicate rows and RMSE.
5. Assert the selected row is the row with the lowest finite RMSE.
```

---

# 11. Regression gates for result-changing switches

Split implementation milestones into two classes.

Audit-only additions:

```text
cell-line collision tables
drug annotation audit table
curated mapping template
additional neutral correlation/rank tables
metadata tables
session metadata under output/metadata/
```

Expected behavior for audit-only additions: canonical legacy outputs should remain byte-identical where feasible, or semantically identical when XLSX/PDF metadata makes byte identity unstable.

Scientific-result-changing switches:

```text
switching canonical matching from raw names to normalized keys
switching enrichment from legacy categories to curated categories
changing enrichment selected-drug thresholds
changing canonical enrichment permutation count
changing canonical metric from Z_SCORE to LN_IC50 or AUC
changing duplicate-selection strategy for canonical results
changing PubChem annotations used by canonical enrichment
```

Required behavior for every scientific-result-changing switch:

1. Implement the switch in its own milestone.
2. Run the baseline-vs-candidate regression comparison immediately after the switch.
3. Write a regression report under `output/metadata/` with changed file hashes, changed table dimensions, changed cell-line matches, changed drug annotations, changed correlations, changed selected enrichment drugs, and changed enrichment summaries as applicable.
4. Identify the first divergent intermediate table, not just the final PDF/XLSX.
5. Explicitly state whether the change is accepted for canonical outputs or retained as supplemental only.

Do not combine multiple scientific-result-changing switches in one milestone. If two changes are made together, the regression report cannot explain which change caused the final result delta.

---

# 12. Revised acceptance criteria

The implementation is complete when:

1. The current canonical `Z_SCORE` pipeline still runs.
2. Existing canonical enrichment outputs are either unchanged or reproduced under a clearly labeled `legacy_Z_SCORE` name.
3. Cell-line normalized-key collision audits are written.
4. Ploidy key collisions with discordant ploidy values fail loudly.
5. Full canonical `Z_SCORE` correlation statistics are saved to TSV and XLSX.
6. The canonical XLSX has one sheet per cancer type.
7. Supplemental all-metric correlation statistics are saved separately.
8. Drug annotation audit is written, but enrichment still uses the legacy category unless an approved curated mapping is supplied.
9. A curated drug-class mapping template is generated.
10. Gemcitabine rank summary is generated for canonical `Z_SCORE` and supplemental metrics.
11. All-cancer tissue-adjusted models are generated as supplemental checks.
12. Enrichment mode and permutation count are explicitly logged.
13. Tissue-adjusted model eligibility thresholds and model statuses are explicitly logged.
14. PubChem subset refresh preserves unrequested cached rows by default.
15. `lowest_rmse` duplicate selection fails clearly when `RMSE` is unavailable.
16. Session information and run parameters are saved under `output/metadata/`.
17. Any scientific-result-changing switch has its own regression report before acceptance.

---

# 13. Explicit non-goals for this implementation pass

Do not change the main manuscript figure yet.

Do not switch the main result from `Z_SCORE` to `LN_IC50` or `AUC`.

Do not let `custom_set_candidate.tsv` supersede PubChem or legacy categories yet.

Do not remove the legacy “last listed category” rule until an approved curated replacement exists.

Do not interpret correlation direction in code beyond neutral rank labels.

Do not silently collapse cell-line key collisions unless they are exact or near-exact duplicates and the collapse is recorded.

Do not switch canonical outputs to normalized matching, curated categories, different permutation counts, or a different duplicate-selection rule without a dedicated regression report accepting that result change.
