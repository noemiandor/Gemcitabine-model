# GDSC Drug-Class Simplification Plan

## Purpose

Simplify the Figure 1 GDSC drug-class enrichment workflow so that category assignment is a direct, auditable input instead of a derived workflow with PubChem caches, legacy fallbacks, proposed categories, manual override reconstruction, and generated curation tables.

The new source of truth is the downloaded workbook:

```text
/Users/4470246/Downloads/drug_class_final_used_with_primary_secondary_corrected.xlsx
```

Initial inspection shows:

- Sheet `drug_class_final_used_with_prim` has 286 drug rows and 44 columns.
- Sheet `Drug counts` summarizes 19 primary anticancer classes.
- `primary_anticancer_class` has 34 drugs in `Other`, down from 58 in the previous `Other reviewed` class.
- The workbook covers the same 286-drug universe as the current correlation-eligible drug set.
- SHA-256 checksum of the downloaded workbook, computed on 2026-06-30:
  ```text
  58f1d714b3821e1f9bc2bf280141b28398b6ae267f95d7a8d29c70c1f503af72
  ```
- Excel readers may coerce 13 numeric-only `drug_key` values unless explicitly normalized:
  ```text
  123138, 123829, 150412, 50869, 615590, 630600, 667880,
  720427, 729189, 741909, 743380, 765771, 776928
  ```
- The workbook still contains old-schema columns from the previous 10-class workflow. Those columns are provenance only and must not drive enrichment.
- `secondary_suggested_class` is not a controlled ontology. It should be treated as a secondary mechanism note, not as a second enrichment class.

## End State

The GDSC analysis should:

1. Read one reviewed category workbook.
2. Use `primary_anticancer_class` as the only enrichment grouping variable for Figure 1.
3. Preserve `secondary_suggested_class`, `classification_note`, `revision_note`, `gdsc_putative_target`, `gdsc_pathway_name`, `source_annotation_used`, and URL columns as metadata only.
4. Fail fast if the workbook is incomplete, duplicated, has missing primary classes, or does not exactly cover the correlation-eligible drug universe.
5. Stop generating or consuming legacy category inputs such as:
   - `drug_class_final_curated.tsv`
   - `drug_class_category_schema.tsv`
   - `custom_set_candidate.tsv`
   - `pubchem_drug_annotations.tsv`
   - `wrong_assignment_examples.tsv`
6. Remove the `legacy` and `proposal` category modes after the workbook-driven outputs have been reviewed and accepted.
7. Remove deprecated generated tables after sign-off, including at least:
   - `drug_annotation_audit.tsv`
   - `drug_class_final_curated_TEMPLATE.tsv`
   - `drug_class_curation_evidence_by_drug_id.tsv`
   - `drug_class_curation_input.tsv`
   - `drug_class_final_used.tsv`
   - `drug_class_assignment_diff.tsv`
   - `drug_class_unassigned_failures.tsv`
   - `drug_class_reviewed_exclusions.tsv`
   - `drug_class_low_count_actions.tsv`
8. Keep Figure 1 materialization and Manager provenance aligned with the workbook-driven run, so promoted assets and `metadata/input_manifest.tsv` point to the reviewed workbook rather than legacy curation files.

## Non-Goals

- Do not reintroduce PubChem lookup or DrugBank scraping in the main analysis.
- Do not infer missing classes from GDSC fields during the main analysis.
- Do not silently assign unclassified drugs to `Other`.
- Do not keep duplicate "used" copies of the workbook as generated outputs unless needed for a temporary validation report.
- Do not change the drug-response metric, duplicate drug-cell-line resolution policy, or ploidy matching in this refactor.

## Proposed Input Contract

Canonical input file:

```text
Code/gdsc_ploidy_analysis/data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx
```

This keeps the change local to the existing module conventions. A later repository-wide data-layout cleanup can move it under `Data/public_data/gdsc_ploidy_analysis/manual/` if desired.

The implementation must copy the downloaded workbook into the repo before any code points to it. The code must not read from `/Users/4470246/Downloads`.

Expected workbook identity:

```text
source_path=/Users/4470246/Downloads/drug_class_final_used_with_primary_secondary_corrected.xlsx
repo_path=Code/gdsc_ploidy_analysis/data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx
sha256=58f1d714b3821e1f9bc2bf280141b28398b6ae267f95d7a8d29c70c1f503af72
import_date=2026-06-30
```

Required sheet:

```text
drug_class_final_used_with_prim
```

Required columns for analysis:

```text
drug
drug_key
primary_anticancer_class
secondary_suggested_class
include_in_enrichment
```

Interpretation of these fields:

- `primary_anticancer_class` is the only class used for enrichment, heatmap columns, and Figure 1 drug-class labels.
- `secondary_suggested_class` should be renamed or handled in code as `secondary_mechanism_note` where possible. It is descriptive metadata only.
- `include_in_enrichment` must be parsed as a boolean. The current workbook has all drugs included; future workbooks should still be validated.

Required columns for provenance/diagnostics:

```text
gdsc_putative_target
gdsc_pathway_name
classification_note
revision_note
source_annotation_used
revision_source_basis
pubchem_query_url
drugbank_query_url
```

Optional ordering sheet:

```text
Drug counts
```

If present, its class counts must exactly match the computed included counts. Use its `primary_anticancer_class` order for heatmap columns and legends. If absent, order classes by decreasing count, then alphabetically, with `Other` last.

Ignored legacy/provenance columns:

```text
approved_final_category_id
approved_final_display_label
category_mode
category_source
legacy_group_used_for_enrichment
schema_*
curated_mapping_*
previous_*
superseded_*
```

These columns may remain in the workbook as audit history, but they must never drive `coxIn`, heatmap labels, category IDs, class ordering, or manuscript-facing output.

Key normalization rules:

- Convert all `drug`, `drug_key`, and GDSC `DRUG_NAME` keys to trimmed character strings before matching.
- Preserve numeric-looking drug keys as exact integer-like strings, e.g. `123138`, not `123138.0`.
- Apply the existing uppercase key normalization after character coercion.

Inclusion rules:

- Use only rows with `include_in_enrichment == TRUE`.
- Fail if any included row has missing `primary_anticancer_class`.
- Fail if the correlation-eligible drug universe has a key missing from the workbook.
- Fail if the workbook has duplicate `drug_key` values after normalization.
- Fail if the workbook has extra keys not present in the correlation-eligible universe, unless an explicit `--allow-extra-class-rows` development option is introduced and defaults to false.

Low-count class policy:

- Retain all primary classes in the first candidate run so the reviewed workbook is represented faithfully.
- Flag classes with five or fewer included drugs as exploratory in review outputs and interpretation. Initial low-count classes are expected to include `Hedgehog pathway inhibitors` (n=1), `Proteasome inhibitors` (n=2), `Hormone therapy` (n=4), and `p53/MDM2 pathway` (n=5).
- Do not collapse or exclude low-count classes after seeing p-values. If collapsing is desired later, define the rule before re-running and document it as a separate analysis decision.

Per-run assignment output:

Write a direct audit table from the workbook, e.g.

```text
drug_class_primary_secondary_used.tsv
```

Required columns:

```text
drug
drug_key
normalized_drug_key
primary_anticancer_class
primary_anticancer_class_slug
secondary_mechanism_note
include_in_enrichment
source_workbook
source_workbook_sha256
```

This table replaces the old generated curation/audit tables for active runs.

## Implementation Phases

### Phase 1: Add Workbook-Driven Category Loading Beside Existing Path

Goal: produce candidate results from the new workbook without deleting legacy paths yet.

Tasks:

1. Copy the workbook into the module:
   ```text
   Code/gdsc_ploidy_analysis/data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx
   ```
2. Add a narrow reader/validator in a new `src/drug_class_workbook.R`:
   - `read_primary_secondary_drug_classes(path, sheet = "drug_class_final_used_with_prim")`
   - `normalize_excel_drug_key(x)`
   - `validate_primary_secondary_class_coverage(correlation_drugs, class_rows)`
   - `make_enrichment_class_table(class_rows)`
   - `read_primary_secondary_class_order(path, sheet = "Drug counts")`
   Prefer `openxlsx` because the module already uses it; avoid adding a new Excel dependency unless necessary.
3. Add a command-line argument:
   ```text
   --drug-class-workbook=Code/gdsc_ploidy_analysis/data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx
   ```
4. Add a temporary category mode:
   ```text
   --category-mode=primary_secondary
   ```
   This mode should bypass PubChem, `custom_set_candidate.tsv`, schema validation, curated mapping loading, proposal generation, and legacy normalization.
5. Refactor the runner so category assignment is built after correlation-eligible drugs are known:
   - Split "compute drug-ploidy correlations" from "build category assignments".
   - Add a single `build_category_assignment(...)` helper returning `coxIn`, `category_suffix`, `category_source_label`, `group_order`, validation rows, and metadata rows.
   - Move legacy input checks, schema loading, PubChem loading, audit generation, and curation table generation inside the legacy-only branch during the transition.
   - Allow `--analysis-mode=manuscript --category-mode=primary_secondary` during Phase 1 validation.
6. For `primary_secondary`, build `coxIn` directly:
   ```text
   drug = normalized drug_key
   group = primary_anticancer_class
   category_id = normalized primary_anticancer_class or a stable slug
   secondary_mechanism_note = metadata only
   ```
7. Write new outputs with distinct names:
   ```text
   drugsVsPloidyCorr_primary_secondary_Z_SCORE.xlsx
   class_enrichment_primary_secondary_Z_SCORE.tsv
   class_enrichment_primary_secondary_metadata.tsv
   class_enrichment_selected_drugs_primary_secondary_Z_SCORE.tsv
   drug_class_primary_secondary_counts.tsv
   drug_class_workbook_validation.tsv
   drug_class_primary_secondary_used.tsv
   ```
8. Keep `drugsVsPloidyCorr.xlsx` as the compatibility workbook copied from the active category mode, because plotting scripts currently consume this generic name.
9. Update Manager in Phase 1, not after cleanup:
   - Add `--gdsc-drug-class-workbook`.
   - Pass the workbook to the GDSC runner.
   - Include the workbook in `input_paths_for_module gdsc` and in `metadata/input_manifest.tsv`.
   - Keep `--gdsc-category-mode primary_secondary` as a temporary transition option until Phase 3 removes category modes.
10. Update `Code/tools/materialize_figure_assets.py` and `figures/Figure1/manifest.tsv` expectations when new `primary_secondary` filenames are introduced, so Manager-promoted Figure 1 assets cannot silently point to old curated outputs.
11. Update plot scripts for primary-class labels:
   - Preserve acronyms and symbols such as `WNT`, `p53/MDM2`, and `DNA`.
   - Wrap or reduce font size for long labels such as `Chaperone/protein-homeostasis inhibitors`.
   - Validate labels in both heatmaps and enrichment panels.
12. Add tests before deleting any legacy code:
   - workbook reader succeeds on the real workbook
   - missing required columns fail
   - duplicate normalized keys fail
   - missing/extra keys fail
   - all 13 numeric-looking keys are preserved exactly
   - `Drug counts` equals computed included counts and supplies column order
   - ignored old-schema columns do not affect `coxIn`
   - `secondary_suggested_class` is carried only as metadata
   - figure-label formatting preserves all 19 primary labels
   - a temp-output integration run for `primary_secondary` does not require legacy category files

Validation:

```sh
Rscript Code/gdsc_ploidy_analysis/tests/smoke_test.R

Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R \
  --analysis-mode=manuscript \
  --category-mode=primary_secondary \
  --drug-class-workbook=Code/gdsc_ploidy_analysis/data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx \
  --output-dir=Results/_validation/public_data/gdsc_ploidy_analysis/runs/primary_secondary_check
```

Inspect:

- Class counts, especially `Other`.
- Missing/extra/duplicate key validation.
- Figure 1 enrichment panels and clustered heatmaps.
- Whether any class label is too long for heatmap axes.
- Whether the new enrichment output has expected dimensions: cancer types by primary anticancer classes.
- `lowpIsSens` and `highpIsSens` sheets are cancer contexts by 19 primary classes, have identical row/column order, contain numeric p-values in `[0,1]`, and do not have unexpected `NA`.
- `metadata/input_manifest.tsv` contains the workbook path and SHA-256 checksum.

### Phase 2: Review Candidate Results Before Removing Legacy Code

Goal: confirm that the simplified primary-class results are scientifically acceptable.

Tasks:

1. Compare old curated outputs versus new primary-secondary outputs:
   - category counts
   - enrichment p-values
   - selected drugs by direction
   - Figure 1 heatmaps
   - `ploidyVsDrugSensitivity.pdf` readability and category coloring
2. Generate a compact review table:
   ```text
   tables/drug_class_primary_secondary_review_summary.tsv
   ```
   Suggested columns:
   ```text
   primary_anticancer_class
   n_drugs
   n_low_ploidy_selected
   n_high_ploidy_selected
   min_low_pvalue
   min_high_pvalue
   ```
3. Decide whether `secondary_suggested_class` should appear anywhere:
   - Default recommendation: keep it out of enrichment and figures.
   - Use it only in supplemental review tables or manuscript notes if needed.
4. Review low-count and `Other` behavior before interpreting class enrichment:
   - Treat n<5 classes as exploratory unless a prespecified collapsing rule is added.
   - Do not interpret `Other` as a biological mechanism class.
5. Update manuscript wording for Figure 1:
   - Avoid "PubChem drug classes."
   - Use "reviewed anticancer mechanism classes based on GDSC target/pathway annotations with manual primary/secondary review."
   - Recompute or soften any old broad-class claims such as cytotoxic/signaling patterns until the new 19-class results are inspected.
   - If broad parent families are desired for narrative interpretation, define them as a separate display/interpretation mapping, not as the enrichment grouping variable.

Validation:

```sh
bash Manager.sh --mode check-only --modules gdsc --run-id primary_secondary_check

bash Manager.sh \
  --mode standard \
  --modules gdsc \
  --gdsc-category-mode primary_secondary \
  --run-id primary_secondary_candidate \
  --overwrite \
  --no-update-latest
```

Acceptance criteria:

- The pipeline completes with no missing category assignments.
- `Other` is reduced relative to the old curated result and contains only intentionally retained drugs.
- Figure 1 panels are readable and use the new primary class names.
- The exact workbook path and checksum are recorded in `metadata/input_manifest.tsv` and run metadata.
- The user/scientific reviewer confirms the biological class ontology is acceptable.

### Phase 3: Make Workbook-Driven Classes The Only Supported Path

Goal: remove category modes and fallbacks after the candidate results are accepted.

Tasks:

1. Remove `--category-mode` from `run_gdsc_ploidy_analysis.R`; if keeping the CLI for compatibility, accept only `primary_secondary` and emit an error for old values.
2. Remove Manager options:
   - `--gdsc-category-mode`
   - `gdsc_category_mode`
   - references to `curated|legacy|proposal`
3. Update Manager GDSC input checks and `docs/manuscript_figure_module_registry.tsv` to require only:
   - `Code/gdsc_ploidy_analysis/data/raw/GDSC2_fitted_dose_response_24Jul22.txt`
   - `Code/gdsc_ploidy_analysis/data/raw/ploidyAcrossCellLines_V1.txt`
   - the new primary/secondary workbook
4. Remove main-analysis dependencies on:
   - `small_molecule_20200407234909.csv`, unless still needed for non-category logic
   - `custom_set_candidate.tsv`
   - `drug_class_category_schema.tsv`
   - `drug_class_final_curated.tsv`
   - `pubchem_drug_annotations.tsv`
   - `wrong_assignment_examples.tsv`
5. Remove `legacy` and `proposal` branches from `src/run_gdsc_ploidy_analysis.R`.
6. Simplify `src/annotations.R` to keep only workbook reading, validation, class-count summaries, and small helper functions still used by the run.
7. Remove or archive no-longer-called files:
   - `src/pubchem_client.R`
   - `src/refresh_pubchem_annotations.R`
   - top-level `refresh_pubchem_annotations.R`
   - `references/annotate_from_pubchem_legacy.R`
   - `tests/test_pubchem_client.R`
   - obsolete portions of `tests/test_annotations.R`
8. Update `README.md`, `docs/FigureCodeMap.md`, `docs/manuscript_figure_module_registry.tsv`, `docs/manuscript_figure_output_standardization_plan.md`, `figures/Figure1/manifest.tsv`, and relevant manuscript Methods/caption text to describe the new source of truth.

Validation:

```sh
Rscript Code/gdsc_ploidy_analysis/tests/smoke_test.R

Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R \
  --analysis-mode=manuscript \
  --output-dir=Results/_validation/public_data/gdsc_ploidy_analysis/runs/primary_secondary_final

bash Manager.sh --mode check-only --modules gdsc --run-id primary_secondary_final_check

bash Manager.sh \
  --mode standard \
  --modules gdsc \
  --run-id primary_secondary_final \
  --overwrite \
  --no-update-latest
```

### Phase 4: Remove Deprecated Inputs, Outputs, And Baselines

Goal: clean the repository after the simplified workflow is validated.

Tasks:

1. Delete deprecated input tables:
   ```text
   Code/gdsc_ploidy_analysis/data/manual/custom_set_candidate.tsv
   Code/gdsc_ploidy_analysis/data/manual/drug_class_category_schema.tsv
   Code/gdsc_ploidy_analysis/data/manual/drug_class_final_curated.tsv
   Code/gdsc_ploidy_analysis/data/derived/pubchem_drug_annotations.tsv
   Code/gdsc_ploidy_analysis/tests/fixtures/wrong_assignment_examples.tsv
   ```
2. Delete deprecated generated outputs from tracked/current result directories only after the new results have been promoted:
   ```text
   drug_annotation_audit.tsv
   drug_class_final_curated_TEMPLATE.tsv
   drug_class_curation_evidence_by_drug_id.tsv
   drug_class_curation_input.tsv
   drug_class_final_used.tsv
   drug_class_assignment_diff.tsv
   drug_class_unassigned_failures.tsv
   drug_class_reviewed_exclusions.tsv
   drug_class_low_count_actions.tsv
   class_enrichment_legacy_*
   drugsVsPloidyCorr_legacy_*
   ```
3. Replace old baselines or remove milestone regression artifacts that exist only to validate the old annotation workflow:
   ```text
   Code/gdsc_ploidy_analysis/baseline/
   Code/gdsc_ploidy_analysis/output_milestones/
   ```
   Do this only if the repo no longer needs historical regression comparisons.
4. Commit cleanup separately from the implementation commit so provenance remains clear.

Validation:

```sh
rg -n "pubchem|legacy|proposal|drug_class_final_curated|custom_set_candidate|wrong_assignment|drug_class_category_schema" \
  Code/gdsc_ploidy_analysis Code/tools Manager.sh docs figures GemcitabinePaper.tex

Rscript Code/gdsc_ploidy_analysis/tests/smoke_test.R
bash Manager.sh --mode check-only --modules gdsc --run-id post_cleanup_check
```

Any remaining hits must be either:

- historical documentation intentionally retained under `docs/Archive` or
- manuscript explanatory text that no longer implies active code usage.

## Output Naming After Cleanup

Recommended final output names:

```text
drugsVsPloidyCorr_primary_secondary_Z_SCORE.xlsx
drugsVsPloidyCorr.xlsx                    # compatibility copy
class_enrichment_primary_secondary_Z_SCORE.tsv
class_enrichment_primary_secondary_metadata.tsv
class_enrichment_selected_drugs_primary_secondary_Z_SCORE.tsv
drug_class_primary_secondary_counts.tsv
drug_class_workbook_validation.tsv
ploidy_enrichment_panels_primary_secondary_ABC.png
ploidy_enrichment_panels_primary_secondary_ABC.pdf
ploidy_enrichment_clustered_primary_secondary_lowpIsSens_clustermap.png
ploidy_enrichment_clustered_primary_secondary_highpIsSens_clustermap.png
ploidy_enrichment_clustered_primary_secondary_shared_order.png
```

If stable manuscript-facing paths are already handled by `figures/Figure1/manifest.tsv`, the internal run filenames can change as long as the manifest records the new source files.

## Risks And Controls

| Risk | Control |
|---|---|
| Excel reads numeric drug keys as numbers | Normalize all keys through a string-preserving helper and test the 13 numeric-only keys explicitly. |
| New class names make heatmap labels unreadable | Inspect Figure 1 and clustered heatmaps; add wrapping rules to plotting helpers if needed. |
| Removing old tables breaks downstream scripts | Use `rg` before deletion and run Manager check-only plus standard GDSC run. |
| Review workbook contains extra/non-analysis rows | Fail by default on extras; only allow extras with explicit development flag. |
| `secondary_suggested_class` is mistaken for enrichment class | Keep enrichment strictly on `primary_anticancer_class`; write secondary only to metadata/review summaries. |
| Manuscript still says PubChem categories | Update captions/methods to say reviewed primary anticancer mechanism classes from GDSC target/pathway annotations and manual review. |
| Old 10-class columns in the workbook accidentally drive output | Add tests that mutate ignored legacy columns and confirm `coxIn` and class counts are unchanged. |
| Low-count classes create overinterpreted single-drug signals | Flag n<5 classes as exploratory in review outputs and manuscript interpretation unless a prespecified collapse rule is adopted. |
| Manager promotes stale curated files | Update `materialize_figure_assets.py`, `figures/Figure1/manifest.tsv`, and Manager input manifests in Phase 1. |

## Definition Of Done

- The GDSC run no longer reads PubChem caches or generated curation TSVs.
- The only drug-class assignment input is the reviewed primary/secondary workbook.
- Figure 1 source panels are regenerated from workbook-driven primary classes.
- The input manifest records the workbook path and checksum.
- Deprecated category modes, refresh code, legacy tests, and deprecated generated tables are removed or archived.
- `Rscript Code/gdsc_ploidy_analysis/tests/smoke_test.R` passes.
- `bash Manager.sh --mode check-only --modules gdsc` passes.
- A standard GDSC manager run regenerates Figure 1 source outputs without requiring any deprecated category files.
