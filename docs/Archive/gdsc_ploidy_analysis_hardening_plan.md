# GDSC Ploidy Analysis Hardening Plan

## Objective

Make `Code/gdsc_ploidy_analysis` reproducible, auditable, and safe to rerun without relying on live PubChem state, silent error swallowing, or implicit local package state. The end state should preserve the current scientific workflow while making each external dependency, annotation decision, and analysis failure explicit.

## Scope

This plan covers the active pipeline under `Code/gdsc_ploidy_analysis`, especially:

- `run_gdsc_ploidy_analysis.R`
- `annotate_from_pubchem.R`
- `README.md`
- `custom_set_candidate.tsv`
- `data/`
- `output/`
- vendored `.Rlibs/`
- legacy/reference scripts

The plan does not attempt to redesign the biological analysis. Changes should first reproduce the current committed outputs as closely as possible, then make intentional behavior changes visible through documented output diffs.

## Current Problems

1. The main pipeline calls PubChem during every run through `ChemmineR::pubchemName2CID()` and `textreadr::read_html()`, so offline reruns can fail and online reruns can change as PubChem content or formatting changes.
2. PubChem JSON is parsed as text with line windows and regexes, making classification dependent on response formatting rather than response structure.
3. `try()` suppresses enrichment and plotting failures, allowing partial or invalid outputs to look successful.
4. The final plot label says `IC50`, but the script correlates ploidy with `Z_SCORE`.
5. Dependency handling is implicit. The repo contains `.Rlibs/textreadr`, but the active script does not add that path to `.libPaths()`, and package versions are not pinned.
6. `SIGNALING` and `CYTOTOXIC` are assumed to exist during matrix ordering, but that assumption is not validated after annotation normalization.
7. Duplicate drug/cell-line records are resolved by keeping the first row, which is order-dependent.
8. Source files, raw data, generated output, vendored packages, and legacy scripts are mixed in one directory.

## Target Design

The main run should be deterministic by default:

```text
raw inputs + cached annotation table + explicit config -> analysis outputs
```

Network access should happen only in an explicit refresh command:

```text
raw drug list -> PubChem refresh script -> reviewed cached annotation table
```

The active pipeline should:

- Use cached annotations by default.
- Fail fast with clear messages on missing dependencies, missing required groups, malformed inputs, or enrichment failures.
- Parse PubChem responses with structured JSON parsing when refreshing annotations.
- Record enough metadata to audit which inputs, package versions, and options produced an output.
- Keep generated outputs separate from source-controlled inputs and scripts.

## Milestone Regression Protocol

Every phase below must include a regression test before it is considered complete. The test must answer two questions:

1. Did the final scientific outputs change?
2. If they changed, which implementation change caused the difference?

### Baseline Before Any Refactor

Before changing behavior, create a reproducible baseline from the current committed pipeline:

```text
Code/gdsc_ploidy_analysis/baseline/
  outputs/
  intermediate/
  metadata/
  checksums.tsv
  result_summary.tsv
```

The baseline should contain:

- Final output checksums for `drugsVsPloidyCorr.RData`, `coxIn.RData`, `drugsVsPloidyCorr.xlsx`, `drugsVsPloidyCorr.pdf`, and `ploidyVsDrugSensitivity.pdf`.
- Intermediate tables exported in text form, including raw correlations, normalized annotations, filtered annotations, low-ploidy enrichment matrix, and high-ploidy enrichment matrix.
- Input file checksums.
- R package versions and `sessionInfo()`.
- Runtime configuration, including metric, random seed, permutation count, annotation source, duplicate handling rule, and category filtering threshold.

Binary artifacts such as PDFs and Excel files should not be the only comparison targets. Each run must also write machine-readable TSV or RDS intermediates so differences can be attributed without manually inspecting binary files.

### Required Comparison Outputs

After each milestone, run the full pipeline into a milestone-specific directory:

```text
Code/gdsc_ploidy_analysis/output_milestones/<phase_name>/
```

Then compare that directory against the baseline and write:

```text
Code/gdsc_ploidy_analysis/output_milestones/<phase_name>/regression_report.md
Code/gdsc_ploidy_analysis/output_milestones/<phase_name>/regression_diff_summary.tsv
```

The report must include:

- Whether final outputs changed.
- Which final outputs changed.
- Which intermediate table first diverged.
- The row, column, drug, cancer type, group, or plot label responsible for the first meaningful divergence.
- Whether the difference is expected or unexpected.
- The implementation reason for expected differences.
- The follow-up action for unexpected differences.

### Suggested Comparison Script

Add a reusable comparison script early in the work:

```text
Code/gdsc_ploidy_analysis/tests/compare_analysis_outputs.R
```

The script should compare:

- Input checksums.
- Package/version metadata.
- Annotation cache rows and normalized groups.
- Duplicate-resolution summaries.
- Raw drug/ploidy correlation vectors by cancer type.
- Enrichment matrices, with numeric tolerance for permutation-based p-values.
- Workbook sheet values after reading them into data frames.
- PDF text labels where feasible, or plot metadata exported alongside PDFs.

The script should emit a non-zero exit code for unexpected differences. Expected differences must be declared in a small milestone config file, for example:

```text
Code/gdsc_ploidy_analysis/tests/expected_differences/<phase_name>.yml
```

This prevents broad "known differences" from hiding unrelated regressions.

## Phase 1: Add Cached Annotation Inputs

### Implementation Tasks

1. Create a new cached annotation file, for example:

   ```text
   Code/gdsc_ploidy_analysis/data/derived/pubchem_drug_annotations.tsv
   ```

2. Include at least these columns:

   ```text
   drug
   drugName
   pubchem_cid
   drugCategory_Pubchem
   annotation_source
   annotation_status
   retrieved_at
   pubchem_url
   notes
   ```

3. Populate the initial cache from the currently committed `output/coxIn.RData` or by running the existing annotation path once in a controlled refresh step.
4. Treat `custom_set_candidate.tsv` as a manual override input, not as an emergency patch inside the pipeline.
5. Add a function such as `load_drug_annotations(annotation_file, custom_file)` that:

   - Reads the cached annotation file.
   - Validates required columns.
   - Applies manual overrides deterministically.
   - Reports drugs that remain unclassified.
   - Returns the same columns currently expected by the downstream analysis.

6. Update `run_gdsc_ploidy_analysis.R` to use the cached annotation file by default.
7. Add a command-line option or config variable for the annotation file path.

### Acceptance Criteria

- Running `run_gdsc_ploidy_analysis.R` does not make network requests.
- The run succeeds without internet access if required R packages are installed.
- The cached annotation file is sufficient to recreate `coxIn` groups used by enrichment.
- Missing cached annotations fail with a clear message unless an explicit `allow_unclassified = TRUE` option is set.

### Milestone Testing

1. Run the original pipeline once to create or refresh the baseline artifacts.
2. Run the cached-annotation pipeline with network disabled or blocked.
3. Compare normalized annotations, `coxIn.RData`, enrichment matrices, workbook sheets, and final PDFs against baseline.
4. Expected result: no scientific output changes if the cache was generated from the same annotations used by the baseline.
5. If outputs change, the regression report must identify whether the first divergence is in the cached annotation table, manual override application, category normalization, or downstream enrichment.
6. Do not proceed until every annotation difference is either corrected or explicitly accepted with a drug-level explanation.

## Phase 2: Split PubChem Refresh From Main Analysis

### Implementation Tasks

1. Replace the active role of `annotate_from_pubchem.R` with a refresh-only script or function, for example:

   ```text
   Code/gdsc_ploidy_analysis/src/refresh_pubchem_annotations.R
   ```

2. The refresh script should accept:

   - Input drug list.
   - Existing annotation cache.
   - Output cache path.
   - Optional `--force-refresh` flag.
   - Optional `--limit` flag for testing.

3. Refresh behavior should:

   - Reuse existing cache rows unless `--force-refresh` is set.
   - Query PubChem only for missing or explicitly refreshed drugs.
   - Rate-limit requests.
   - Retry transient HTTP failures with bounded retries.
   - Record failures as `annotation_status = "failed"` rather than silently returning `NA`.
   - Write a log file listing refreshed, reused, failed, and manually overridden drugs.

4. Keep network code out of `run_gdsc_ploidy_analysis.R`.

### Acceptance Criteria

- There are two distinct entrypoints: one for analysis, one for annotation refresh.
- The analysis entrypoint can run offline.
- The refresh entrypoint clearly documents that it requires network access.
- Refresh output is reviewable as a TSV diff before being used by the main analysis.

### Milestone Testing

1. Run the main analysis without invoking the refresh entrypoint.
2. Compare all outputs to the Phase 1 milestone output.
3. Expected result: no changes, because this phase should only split network refresh mechanics from analysis execution.
4. Run the refresh entrypoint in dry-run or `--limit` mode against a small drug subset and verify that it does not modify the committed cache unless an explicit output path is provided.
5. If final outputs change after this phase, the regression report must show whether the analysis accidentally consumed refreshed annotations or whether path/config defaults changed.

## Phase 3: Parse PubChem JSON Structurally

### Implementation Tasks

1. Replace `textreadr::read_html(url)` and line-based parsing with `jsonlite::fromJSON()`.
2. Fetch the PubChem PUG View JSON with a standard HTTP client such as `httr2`, `httr`, or base R `utils::download.file()` plus `jsonlite`.
3. Implement a helper such as:

   ```r
   extract_pubchem_drug_classes <- function(pubchem_json) {
     # Traverse Section nodes by TOCHeading and return structured String values.
   }
   ```

4. Traverse nested `Section` records recursively to find sections where `TOCHeading == "Drug Classes"`.
5. Extract values from structured `Information` fields, not from formatted text lines.
6. Implement fallback extraction from structured `Mechanism of Action` sections only when `Drug Classes` is unavailable.
7. Preserve the current keyword heuristics for inflammatory, cytotoxic, signaling, and metabolic categories only as a clearly named fallback.
8. Add unit-like fixture tests using saved minimal PubChem JSON examples:

   - Drug with a `Drug Classes` section.
   - Drug without `Drug Classes` but with `Mechanism of Action`.
   - Drug with no usable category.
   - Malformed or empty response.

### Acceptance Criteria

- PubChem classification no longer depends on line position or response pretty-printing.
- Category extraction can be tested from local JSON fixtures without network access.
- Malformed JSON or missing expected sections produces a structured failed status and an explanatory note.

### Milestone Testing

1. Run fixture tests for structured JSON parsing before running the full pipeline.
2. Refresh annotations for a controlled subset into a temporary cache, not the canonical cache.
3. Diff old parser output against structured parser output at drug level.
4. Run the main pipeline with the unchanged canonical cache and compare outputs to Phase 2.
5. Expected full-pipeline result: no final output changes when the canonical cache is unchanged.
6. If adopting a newly refreshed structured cache, run a separate intentional-change regression and explain every drug whose category changed, including the PubChem section or fallback heuristic that caused the new category.

## Phase 4: Replace Silent `try()` With Explicit Error Handling

### Implementation Tasks

1. Add a helper for enrichment:

   ```r
   run_enrichment_or_stop <- function(values, annotations, cancer, direction, permute_n, pvalue_cutoff) {
     tryCatch(
       enrichment(values, annotations, permute.n = permute_n, normalize = FALSE, pvalue.cutoff = pvalue_cutoff)$pvalue,
       error = function(e) {
         stop(sprintf("Enrichment failed for cancer=%s direction=%s: %s", cancer, direction, conditionMessage(e)), call. = FALSE)
       }
     )
   }
   ```

2. Record context before each enrichment call:

   - Cancer type.
   - Direction, such as `low_ploidy_sensitive` or `high_ploidy_sensitive`.
   - Number of drugs.
   - Number of annotation groups.

3. Replace `try(barplot(...))` with either:

   - A fail-fast plotting helper, if plots are required outputs.
   - A warning-and-skip helper, if missing plots are acceptable.

4. If warning-and-skip is chosen for plotting, write a machine-readable plot warning log under `output/logs/`.
5. Ensure graphics devices are closed on error with `on.exit(dev.off(), add = TRUE)` around each PDF block.

### Acceptance Criteria

- Enrichment failures stop the pipeline with cancer and direction in the error message.
- Plotting failures are either fatal or logged explicitly, not silently suppressed.
- PDF devices are not left open after errors.
- A failed run does not leave outputs that appear complete without an accompanying failure log.

### Milestone Testing

1. Run the full pipeline on normal inputs and compare outputs to Phase 3.
2. Expected result: no scientific output changes, because error handling should not alter successful computations.
3. Add negative tests using small fixtures that force enrichment failure and plotting failure.
4. Verify those failures include cancer type, direction, drug count, and the original error message.
5. If normal-run outputs change, the regression report must identify whether a previously swallowed failure is now exposed. If so, classify the previous output as invalid or partial rather than treating the change as a numerical regression.

## Phase 5: Correct Metric Naming and Scientific Labels

### Implementation Tasks

1. Introduce a single metric configuration variable:

   ```r
   metric <- "Z_SCORE"
   metric_label <- "GDSC Z-score"
   ```

2. Use `metric_label` in all plot labels, output metadata, and README text.
3. If the intended metric is actually `LN_IC50`, make that an explicit configuration change and rerun output comparison.
4. Add validation that `metric` exists in the GDSC input table.
5. Add a short README note explaining what positive and negative correlations mean for the selected metric.

### Acceptance Criteria

- No output produced from `Z_SCORE` is labeled as `IC50`.
- The selected metric appears in output metadata.
- The README defines the metric and the interpretation of correlation direction.

### Milestone Testing

1. Run the full pipeline with `metric = "Z_SCORE"` and compare numeric intermediates to Phase 4.
2. Expected result: numeric outputs do not change; plot text, metadata, and documentation should change.
3. The regression report should classify PDF checksum changes as label-only changes if exported plot metadata confirms identical plotted values and changed axis text.
4. If any correlation or enrichment value changes, the report must identify whether the configured metric accidentally changed from `Z_SCORE` to another column.

## Phase 6: Make Dependencies Explicit and Reproducible

### Implementation Tasks

1. Choose one dependency strategy:

   - Preferred: use `renv` with `renv.lock`.
   - Minimal: add an installation script and package-version report.

2. If using `renv`:

   - Initialize `Code/gdsc_ploidy_analysis/renv.lock`.
   - Record versions for `EnrichIntersect`, `xlsx`, `plyr`, `ChemmineR`, `jsonlite`, and the selected HTTP package.
   - Document `renv::restore()` in the README.

3. If keeping `.Rlibs`, add this near the top of entrypoint scripts:

   ```r
   local_lib <- file.path(base_dir, ".Rlibs")
   if (dir.exists(local_lib)) {
     .libPaths(c(local_lib, .libPaths()))
   }
   ```

4. Do not vendor only one package unless the dependency model explains why that package is special.
5. Add startup checks that report missing packages with installation guidance.
6. Write session metadata at the end of a successful run:

   ```text
   output/metadata/session_info.txt
   output/metadata/input_manifest.tsv
   ```

### Acceptance Criteria

- A fresh environment has documented steps to install the required R packages.
- Package versions used for a run are recorded.
- `.Rlibs` is either intentionally used and documented or removed from source control.

### Milestone Testing

1. Run the pipeline in the existing environment and compare outputs to Phase 5.
2. Run the pipeline in a restored clean environment if feasible.
3. Compare package metadata between baseline and restored environment.
4. Expected result: no scientific output changes when package versions match.
5. If outputs change because dependency versions changed, the regression report must name the package/version difference and the first divergent intermediate table.
6. If exact reproduction requires package pinning, do not accept this phase until the dependency lockfile or install script reproduces the prior result or documents why reproduction is impossible.

## Phase 7: Validate Required Groups Before Ordering

### Implementation Tasks

1. After category normalization and frequency filtering, validate:

   ```r
   required_groups <- c("SIGNALING", "CYTOTOXIC")
   missing_groups <- setdiff(required_groups, unique(coxIn$group))
   if (length(missing_groups) > 0) {
     stop("Missing required annotation groups after filtering: ", paste(missing_groups, collapse = ", "), call. = FALSE)
   }
   ```

2. Write a category summary table before filtering and after filtering:

   ```text
   output/tables/drug_category_counts_before_filter.tsv
   output/tables/drug_category_counts_after_filter.tsv
   ```

3. Document why singleton categories are removed.
4. Consider making `required_groups` and minimum category frequency configurable.

### Acceptance Criteria

- Missing required groups fail early with a clear message.
- Category filtering decisions are visible in output tables.
- The README documents why `SIGNALING` and `CYTOTOXIC` are used for matrix ordering.

### Milestone Testing

1. Run the full pipeline and compare outputs to Phase 6.
2. Expected result: no scientific output changes for valid inputs, except for new category-count QC tables and metadata.
3. Add a fixture or temporary annotation cache that removes `SIGNALING` or `CYTOTOXIC`.
4. Verify the pipeline fails before enrichment ordering and reports the missing group.
5. If normal outputs change, the regression report must identify whether the validation was added before or after category filtering and whether filtering semantics changed.

## Phase 8: Handle Duplicate Drug/Cell-Line Records Explicitly

### Implementation Tasks

1. Add a duplicate audit before correlations:

   ```r
   duplicate_rows <- dr_sub[duplicated(dr_sub[, c("DRUG_NAME", "CELL_LINE_NAME")]) |
                            duplicated(dr_sub[, c("DRUG_NAME", "CELL_LINE_NAME")], fromLast = TRUE), ]
   ```

2. Write duplicate summaries to:

   ```text
   output/qc/duplicate_drug_cell_line_records.tsv
   output/qc/duplicate_resolution_summary.tsv
   ```

3. Choose and document a deterministic resolution strategy:

   - Use the mean `Z_SCORE` across duplicates.
   - Use the row with lowest `RMSE`.
   - Use a specific release/result identifier ordering.

4. Implement the selected strategy in a helper such as:

   ```r
   resolve_duplicate_drug_cell_lines <- function(drug_table, metric, strategy = "lowest_rmse") {
     # Return one row per DRUG_NAME and CELL_LINE_NAME with audit columns.
   }
   ```

5. Validate that the resulting table has exactly one row per `DRUG_NAME` and `CELL_LINE_NAME`.

### Acceptance Criteria

- Duplicate handling is deterministic and documented.
- The selected resolution strategy is recorded in output metadata.
- The pipeline no longer depends on input row order for duplicate records.

### Milestone Testing

1. First implement duplicate auditing without changing the existing first-row behavior.
2. Run the full pipeline and compare outputs to Phase 7.
3. Expected result for audit-only implementation: no scientific output changes.
4. Add a fixture with duplicated drug/cell-line rows in shuffled order and verify the audit reports them.
5. When switching to a new duplicate-resolution strategy, run an intentional-change regression.
6. The intentional-change report must identify every drug/cancer correlation that changed because of duplicate resolution and summarize how many downstream enrichment values changed.
7. Do not accept a new duplicate strategy without a written scientific rationale for why it is preferable to first-row selection.

## Phase 9: Clean Directory Layout and Repository Hygiene

### Proposed Layout

```text
Code/gdsc_ploidy_analysis/
  README.md
  renv.lock
  src/
    run_gdsc_ploidy_analysis.R
    refresh_pubchem_annotations.R
    annotations.R
    pubchem_client.R
    enrichment_helpers.R
    plotting.R
    validation.R
  data/
    raw/
      GDSC2_fitted_dose_response_24Jul22.txt
      ploidyAcrossCellLines_V1.txt
      small_molecule_20200407234909.csv
    manual/
      custom_set_candidate.tsv
    derived/
      pubchem_drug_annotations.tsv
  references/
    GDSC_analysis_original.R
    GDSC_analysis_repro_reference.R
  output/
    .gitkeep
```

### Implementation Tasks

1. Move active code into `src/`.
2. Move legacy scripts into `references/` and mark them as non-runnable historical references.
3. Move generated outputs out of version control or into a release artifact location.
4. Add `.gitignore` entries for generated files:

   ```text
   Code/gdsc_ploidy_analysis/output/*
   !Code/gdsc_ploidy_analysis/output/.gitkeep
   Code/gdsc_ploidy_analysis/.Rlibs/
   ```

5. Keep raw input data only if licensing and repository size are acceptable.
6. If raw data should remain tracked, add source/provenance documentation for each input file.
7. Update all paths in scripts and README after moves.

### Acceptance Criteria

- Source code, raw inputs, manual inputs, derived annotation cache, references, and generated outputs are clearly separated.
- Generated binary outputs are not committed by default.
- Historical scripts are preserved without being confused for maintained entrypoints.

### Milestone Testing

1. Move files without changing analysis logic.
2. Run the full pipeline from the new entrypoint and compare outputs to Phase 8.
3. Expected result: no scientific output changes.
4. Compare input manifests before and after the move to prove the same raw files and annotation cache were consumed.
5. If outputs change, the regression report must identify the first path-dependent difference, such as loading the wrong cache, missing `.Rlibs`, writing to a different output directory, or failing to source the intended helper.

## Phase 10: Verification and Regression Testing

### Implementation Tasks

1. Add a smoke test script:

   ```text
   Code/gdsc_ploidy_analysis/tests/smoke_test.R
   ```

2. The smoke test should run on a small fixture dataset and verify:

   - Cached annotations are loaded.
   - Duplicate resolution is deterministic.
   - Required groups are validated.
   - Enrichment errors include context.
   - Plot labels use the configured metric label.

3. Add a full-run validation checklist:

   - Run analysis offline.
   - Compare output row/column dimensions to current committed outputs.
   - Compare category counts before and after filtering.
   - Compare enrichment workbook sheets.
   - Inspect PDF labels.
   - Confirm session metadata and input manifests are written.

4. If exact numeric reproduction is expected, add checksums for key intermediate tables.
5. If exact numeric reproduction is not expected after dependency or annotation changes, document tolerated differences and why.
6. Add the baseline-generation script:

   ```text
   Code/gdsc_ploidy_analysis/tests/create_analysis_baseline.R
   ```

7. Add the milestone comparison script:

   ```text
   Code/gdsc_ploidy_analysis/tests/compare_analysis_outputs.R
   ```

8. Add a result-attribution template:

   ```text
   Code/gdsc_ploidy_analysis/tests/templates/regression_report_template.md
   ```

9. Require each milestone PR or commit to include a completed regression report.

### Acceptance Criteria

- A developer can validate the pipeline without making PubChem requests.
- Full-run output differences are either absent or documented.
- The README contains the exact commands for smoke and full validation.
- Every milestone has a regression report that states whether final results changed and why.
- Unexpected final-result changes block the milestone until fixed or explicitly accepted.

### Milestone Testing

1. Run all smoke tests.
2. Run baseline creation from a clean checkout of the pre-hardening commit if available.
3. Run milestone comparison for every completed phase.
4. Verify that known intentional changes are narrow and declared in the relevant expected-differences file.
5. Verify that the final hardening output can be traced back through phase-by-phase reports to explain all differences from the original baseline.

## Suggested Implementation Order

1. Add baseline export and output comparison scripts.
2. Create the original baseline and commit or archive its metadata.
3. Add cached annotation loading while preserving current output paths, then run milestone regression.
4. Add network refresh as a separate script using the current parsing logic temporarily, then run milestone regression.
5. Replace PubChem parsing with structured JSON and fixtures, then run fixture tests and milestone regression.
6. Replace `try()` calls and add device cleanup, then run success-path and failure-path tests.
7. Fix metric labels and required-group validation, then run numeric and label-difference tests.
8. Add duplicate audit first, run no-change regression, then decide whether to change duplicate resolution.
9. Add dependency management and session metadata, then run environment reproduction tests.
10. Reorganize directory layout once behavior is stable, then run path-equivalence regression.
11. Update README with exact validation commands and all accepted result changes.

This order minimizes scientific-output churn early. It first removes the most serious reproducibility issue, then tightens correctness and maintainability.

## Deliverables

- Cached annotation TSV used by default.
- Refresh-only PubChem annotation script with structured JSON parsing.
- Refactored main analysis script with fail-fast enrichment behavior.
- Corrected metric labels and output metadata.
- Dependency restoration instructions or `renv.lock`.
- Duplicate handling audit and deterministic resolution.
- Updated README with commands and interpretation notes.
- Optional repository cleanup commit moving code/data/output into clearer subdirectories.

## Risks and Decisions Needed

1. Decide whether the canonical annotation source is the current committed `coxIn.RData`, a one-time PubChem refresh, or manually reviewed categories.
2. Decide whether generated outputs should remain tracked for manuscript reproducibility or be regenerated as artifacts.
3. Decide whether duplicate drug/cell-line records should use lowest `RMSE`, mean `Z_SCORE`, or another scientific rule.
4. Decide whether exact reproduction of current PDFs/workbooks is required or whether improved annotations/dependencies may intentionally change outputs.
5. Decide whether to adopt `renv` or keep a simpler package installation script.

## Definition of Done

The hardening work is complete when a fresh checkout can:

1. Restore documented R dependencies.
2. Run the main GDSC ploidy analysis without network access.
3. Produce outputs with explicit metric labels, metadata, and QC tables.
4. Fail with actionable messages for missing annotations, missing required groups, enrichment errors, malformed inputs, or missing packages.
5. Refresh PubChem annotations only through a separate explicit command.
6. Explain all generated output differences from the pre-hardening commit.
