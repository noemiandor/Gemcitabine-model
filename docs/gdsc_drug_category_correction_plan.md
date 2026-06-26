# GDSC Drug-Category Correction Plan

## Goal

Correct the drug-class assignments used by the GDSC ploidy-enrichment analysis without introducing new silent misassignments and without leaving any enrichment-eligible drug unassigned.

This plan addresses two sources of error:

1. Legacy normalization logic in `Code/gdsc_ploidy_analysis/src/annotations.R`.
2. Incorrect legacy labels already present in `Code/gdsc_ploidy_analysis/data/derived/pubchem_drug_annotations.tsv`, which was seeded from `baseline_coxIn_RData`.

## Current Failure Mode

The current enrichment pipeline loads cached annotations, assigns `coxIn$group`, then normalizes those groups immediately before enrichment:

- `Code/gdsc_ploidy_analysis/src/run_gdsc_ploidy_analysis.R:358`: loads cached/manual annotations.
- `Code/gdsc_ploidy_analysis/src/run_gdsc_ploidy_analysis.R:378`: applies `normalize_legacy_drug_group()`.
- `Code/gdsc_ploidy_analysis/src/run_gdsc_ploidy_analysis.R:448-464`: passes the already-normalized `coxIn$group` into enrichment.

The problematic normalization rules are:

- `Code/gdsc_ploidy_analysis/src/annotations.R:95`: keeps only the last semicolon-separated PubChem category.
- `Code/gdsc_ploidy_analysis/src/annotations.R:110`: maps `(Anti-)Inflammatory` to `IMMUNOSUPPRESSIVE AGENTS`.
- `Code/gdsc_ploidy_analysis/src/annotations.R:111-119`: maps `PARP INHIBITORS`, `SIGNAL TRANSDUCTION INHIBITORS`, `ENZYME INHIBITORS`, and `TARGETED THERAPIES` to `SIGNALING`.

These rules cause biologically incorrect assignments such as PARP inhibitors to `SIGNALING`, JQ1 to `IMMUNOSUPPRESSIVE AGENTS`, and broad antineoplastic labels to override more informative target biology. Some errors are not introduced by normalization but are already present in the cached table, for example `Entinostat`, `Elesclomol`, and `Leflunomide` entering as `Signaling`. Refametinib should not be treated as a wrong `SIGNALING` example unless the reviewed schema intentionally separates it into a more specific ERK/MAPK signaling subclass.

## Design Principles

- Treat PubChem categories as raw evidence, not as the final enrichment group.
- Do not infer final biology from the last PubChem category in a multi-category string.
- Do not allow broad administrative categories such as `Antineoplastic Agents`, `Targeted therapies`, `Enzyme Inhibitors`, or `Signal Transduction Inhibitors` to override more specific target/pathway information.
- Use the GDSC `PUTATIVE_TARGET` and `PATHWAY_NAME` fields as primary evidence where available.
- Keep a reviewed final drug-to-category table under version control.
- Fail the analysis if any pre-filter correlation-eligible drug lacks an approved final category or an approved exclusion.
- Make every category change auditable by writing old category, evidence fields, final category, and reason.
- Anchor evidence at GDSC `DRUG_ID` where available, then roll up to the current `drug_key`-level analysis only after conflicts among aliases, salts, probes, or duplicate IDs have been resolved explicitly.
- Freeze category granularity, low-count handling, and parent-category collapse rules before inspecting revised enrichment results.

## Definitions

Use these definitions consistently in code, tests, and documentation:

- **Correlation-eligible drug universe**: the unique drug set created from `R` immediately after drug-ploidy correlations are computed and before annotation/category filtering. In current code this begins at `Code/gdsc_ploidy_analysis/src/run_gdsc_ploidy_analysis.R:357`. This is the required curation universe. Implement this as a concrete helper, `derive_correlation_eligible_drugs(R)`, which returns the normalized drug key, display drug name, available GDSC IDs, and a deterministic ordering before any category joins or filters.
- **Approved curation universe**: every correlation-eligible drug after joining to the reviewed curated mapping and category schema. This universe may include reviewed exclusions, but it must not contain missing, blank, unapproved, or schema-invalid categories.
- **Enrichment-used drug universe**: the subset of the approved curation universe with `include_in_enrichment = TRUE`, after any explicitly documented minimum-count or plotting filters are applied.
- **Reviewed exclusion**: a drug intentionally not used for enrichment because it lacks a meaningful category under the approved schema or is unsuitable for another documented reason. This is not the same as unassigned. Reviewed exclusions require `include_in_enrichment = FALSE`, a non-empty `exclusion_reason`, and `curation_status = approved`.

Validation must occur on the correlation-eligible drug universe before any `NA`, empty-category, singleton-category, or display filtering is applied. A successful manuscript run must prove that no drug disappeared because of missing category assignment.

The current analysis remains `drug_key`-level unless the full enrichment pipeline is refactored to `DRUG_ID`-level identities. Therefore the final curated mapping must contain exactly one resolved row per correlation-eligible `drug_key`. If multiple GDSC `DRUG_ID` values map to the same `drug_key`, the evidence can be stored at `DRUG_ID` resolution, but the final mapping must either resolve the key to one approved category or mark the key as an approved exclusion.

## Implementation Steps

### 1. Define The Category Schema

Add a version-controlled schema file:

`Code/gdsc_ploidy_analysis/data/manual/drug_class_category_schema.tsv`

Required columns:

- `schema_version`: stable version string for the reviewed schema.
- `category_id`: stable machine-readable ID, for example `SIGNALING`, `EPIGENETICS_TRANSCRIPTION`, or `GENOME_INTEGRITY_DNA_DAMAGE`.
- `display_label`: label shown in Figure 1 and exported tables.
- `parent_family`: optional broader family used for grouping or summaries.
- `plot_order`: integer order used by heatmaps and bar plots.
- `include_in_enrichment_default`: default logical value.
- `manuscript_allowed`: logical value; manuscript mode must reject categories where this is false.
- `min_count_policy`: how to handle categories with too few drugs, for example `retain`, `collapse_to_parent`, or `exclude_with_report`.
- `min_count_threshold`: prespecified threshold used by `min_count_policy`.
- `collapse_parent_category_id`: required when `min_count_policy = "collapse_to_parent"`.
- `allowed_aliases`: semicolon-separated accepted aliases from old labels or external sources.
- `deprecated_or_broad_terms`: semicolon-separated terms that should not be accepted as final categories without explicit curator approval.
- `is_deprecated_final_class`: logical value; deprecated final classes must not be used in manuscript mode.
- `notes`

The schema should be reviewed before curation begins and before revised enrichment results are inspected. It must resolve the granularity problem explicitly by listing actual approved category IDs, parent/child relationships, and target/pathway decision rules. For example, decide whether PARP/topoisomerase/DNA-replication agents are all `CYTOTOXIC`, a new `GENOME_INTEGRITY_DNA_DAMAGE` class, or a parent/child pair. Similarly, decide whether MEK/ERK, PI3K/AKT/mTOR, and other kinase-pathway drugs share `SIGNALING` or use more specific signaling subclasses. The same decision must be reflected in Figure 1 labels and interpretation.

Broad administrative terms such as `ANTINEOPLASTIC AGENTS`, `TARGETED THERAPIES`, `ENZYME INHIBITORS`, and `SIGNAL TRANSDUCTION INHIBITORS` must be deprecated or marked `manuscript_allowed = FALSE`. They can be retained as raw evidence or legacy comparison labels, but they should not be interpretable final enrichment classes. If broad-only evidence cannot be resolved to a specific approved category, the drug should become `OTHER_REVIEWED` or an approved exclusion.

Schema validation must assert unique `category_id` values, unique non-missing `plot_order` values among manuscript-visible categories, valid boolean fields, no duplicate aliases assigned to multiple categories, valid `collapse_parent_category_id` values, and consistency between `approved_final_category_id` and `display_label`. Prefer deriving the final display label from the schema rather than duplicating free text in downstream mapping files.

The plotting code must stop hard-coding assumptions that `SIGNALING` and `CYTOTOXIC` are required ordering anchors. Heatmap row/column order should come from the schema where applicable, and validation should check the approved schema rather than a fixed two-group list.

### 2. Create A Curation Input Table

Add a script or function that builds two evidence outputs:

`Code/gdsc_ploidy_analysis/output/tables/drug_class_curation_evidence_by_drug_id.tsv`

This table should contain one row per GDSC `DRUG_ID`/drug-name identity where `DRUG_ID` is available. It is the evidence layer used to detect aliases, salts, duplicate IDs, and target/pathway conflicts before rolling up to the analysis key.

Also build:

`Code/gdsc_ploidy_analysis/output/tables/drug_class_curation_input.tsv`

This table should contain one row per correlation-eligible `drug_key` and include:

- `drug`
- `drug_key`
- `drug_id_list`: semicolon-separated unique `DRUG_ID` values from GDSC.
- current legacy enrichment group
- raw `drugCategory_Pubchem`
- parsed PubChem category tokens
- `gdsc_putative_target_values`: semicolon-separated distinct `PUTATIVE_TARGET` values from `GDSC2_fitted_dose_response_24Jul22.txt`
- `gdsc_pathway_name_values`: semicolon-separated distinct `PATHWAY_NAME` values from `GDSC2_fitted_dose_response_24Jul22.txt`
- `gdsc_target_pathway_conflict_flag`: true when one drug key maps to conflicting target/pathway evidence requiring manual review.
- `rollup_conflict_status`: for example `none`, `resolved_same_category`, `requires_manual_resolution`, or `exclude`.
- `rollup_resolution_reason`
- number of GDSC rows supporting each target/pathway value, either as separate count columns or a structured summary string
- current manual override, if any
- whether the drug appears in the maintained wrong-assignment fixture
- proposed final category
- evidence reason

This is the table scientists should review. It must make the bad assignments visible before any correction is applied.

### 3. Add A Reviewed Final Mapping

Create a version-controlled file:

`Code/gdsc_ploidy_analysis/data/manual/drug_class_final_curated.tsv`

Required columns:

- `drug`
- `drug_key`
- `drug_id_list`
- `approved_final_category_id`
- `approved_final_display_label`
- `include_in_enrichment`
- `exclusion_reason_code`
- `exclusion_reason`
- `rollup_resolution_status`
- `rollup_resolution_reason`
- `curation_status`
- `primary_evidence_source`
- `evidence_summary`
- `curator`
- `curation_date`
- `curation_notes`

Rules:

- `approved_final_category_id` must be non-empty for every correlation-eligible drug unless the row is an approved exclusion with `include_in_enrichment = FALSE`.
- `approved_final_category_id` must exist in `drug_class_category_schema.tsv`.
- `curation_status` must be `approved` before the drug is used in enrichment.
- Drugs that genuinely do not fit a specific class must either receive an explicit reviewed category such as `OTHER_REVIEWED` or be marked as an approved exclusion. They must not remain `NA`, blank, or silently dropped.
- `OTHER_REVIEWED` requires an explicit `include_in_enrichment` decision. If included, it must be interpreted as a heterogeneous reviewed category; if excluded, it must be reported as a reviewed exclusion, not hidden.
- Broad categories such as `ANTINEOPLASTIC AGENTS`, `TARGETED THERAPIES`, `ENZYME INHIBITORS`, and `SIGNAL TRANSDUCTION INHIBITORS` are not valid manuscript final classes. Broad-only evidence should lead to `OTHER_REVIEWED` or an approved exclusion.
- The final curated mapping must have one row per correlation-eligible `drug_key`. Duplicate final `drug_key` rows are invalid. If distinct `DRUG_ID` values share a `drug_key`, the final mapping must document a valid `rollup_resolution_status` and reason.
- `exclusion_reason_code` must come from a controlled vocabulary, for example `AMBIGUOUS_MECHANISM`, `NO_SPECIFIC_TARGET_EVIDENCE`, `BROAD_ONLY_EVIDENCE`, `DUPLICATE_ALIAS_COLLAPSED`, `SALT_OR_PRODRUG_COLLAPSED`, `CONTROL_OR_NON_THERAPEUTIC`, `LOW_COUNT_CATEGORY_EXCLUDED`, `CONFLICTING_GDSC_EVIDENCE`, or `OUT_OF_SCOPE_FOR_ENRICHMENT`.

### 4. Replace Last-Category Normalization With Schema-Checked Resolution

Refactor `normalize_legacy_drug_group()` into smaller functions:

- `derive_correlation_eligible_drugs(R)`: extracts the pre-filter drug universe directly after correlations are computed.
- `parse_pubchem_categories(category_string)`: returns all semicolon-separated tokens.
- `normalize_category_token(token)`: standardizes spelling/case only.
- `build_drug_evidence_table(drug_key, gdsc_rows, pubchem_rows, manual_rows)`: aggregates `DRUG_ID`, GDSC target/pathway, PubChem tokens, and legacy labels.
- `infer_candidate_categories(pubchem_tokens, gdsc_targets, gdsc_pathways, schema)`: proposes categories from all evidence fields for reviewer convenience only.
- `resolve_final_drug_category(drug, evidence, curated_mapping, schema, mode)`: returns the approved final category or fails.

The final enrichment group should be chosen by explicit mode:

1. `category_mode = "curated"`: manuscript mode. Require `drug_class_final_curated.tsv` with `curation_status = approved`, schema-valid `approved_final_category_id`, and valid reviewed exclusions for every correlation-eligible drug.
2. `category_mode = "legacy"`: legacy comparison mode only. Use legacy labels and write outputs with `legacy` in the filename and metadata.
3. `category_mode = "proposal"`: development diagnostics only. Allow inferred proposals, write them as `assignment_source = inferred_proposal`, and never name the outputs curated or use them for manuscript Figure 1.

The production/manuscript path should require curated categories or approved exclusions for every correlation-eligible drug. It must fail if the curated mapping or schema is absent and must not silently use fallback inference.

### 5. Correct The Problematic Rule Families

Remove the following behaviors from final assignment:

- Do not choose only `parts[length(parts)]`.
- Do not automatically convert `(Anti-)Inflammatory` to `IMMUNOSUPPRESSIVE AGENTS`.
- Do not collapse `PARP INHIBITORS`, `SIGNAL TRANSDUCTION INHIBITORS`, `ENZYME INHIBITORS`, or `TARGETED THERAPIES` to `SIGNALING`.

Replacement behavior:

- Use all parsed PubChem categories as evidence.
- Prefer specific mechanisms over broad umbrella terms.
- Use GDSC target/pathway when it is more specific than PubChem.
- Examples from `wrongass.txt` should be resolved explicitly:
  - HDAC/BET/DNMT/chromatin targets -> an approved epigenetics/transcription category.
  - PARP/TOP/DNA-replication targets -> the approved genome-integrity/DNA-damage or cytotoxic category defined by the schema.
  - Proteasome targets -> the approved protein-degradation/proteasome category, or an explicitly justified schema category.
  - PI3K/MEK pathway drugs -> `SIGNALING` or a more specific approved signaling subclass from the schema. Refametinib belongs in this family unless the schema defines a more specific ERK/MAPK subclass.
  - Antioxidant or DHODH/metabolic drugs -> reviewed metabolic/other categories, or approved exclusion using a controlled reason code.

### 6. Save The Final Used Mapping

The pipeline should always write the exact mapping used for enrichment:

`Code/gdsc_ploidy_analysis/output/tables/drug_class_final_used.tsv`

This file must include every correlation-eligible drug, not only the included subset. Required columns:

- `drug`
- `drug_key`
- `drug_id_list`
- `approved_final_category_id`
- `approved_final_display_label`
- `legacy_group_used_for_enrichment`
- `pubchem_category`
- `gdsc_putative_target`
- `gdsc_pathway_name`
- `assignment_source`
- `include_in_enrichment`
- `included_in_enrichment`
- `exclusion_reason`
- `exclusion_reason_code`
- `low_count_filter_status`
- `low_count_filter_reason`
- `category_source`
- `category_mode`
- `analysis_metric`
- `schema_version`
- `schema_checksum`
- `curated_mapping_version`
- `curated_mapping_checksum`
- `curation_status`
- `evidence_summary`

This file is the canonical answer to “what drug categories informed Figure 1?”

Also write curated-vs-legacy comparison outputs:

- `drug_class_assignment_diff.tsv`: old group vs curated group, with evidence and reason.
- `class_enrichment_selected_drugs_curated_Z_SCORE.tsv`: selected drugs and curated groups actually used.
- `class_enrichment_curated_Z_SCORE.tsv`: curated enrichment p-values.
- `drugsVsPloidyCorr_curated_Z_SCORE.xlsx`: curated enrichment workbook used by Figure 1 plotting.
- Curated figure outputs with explicit names, for example `ploidy_enrichment_panels_curated_ABC.*` and `ploidy_enrichment_clustered_curated_*`.
- Keep legacy files only as comparison artifacts. Files with `legacy` in the name should not be cited as curated Figure 1 provenance.

Update enrichment metadata so `category_source` is `curated_mapping` for curated outputs. Do not continue labeling curated outputs as `legacy_group_used_for_enrichment`. If generic filenames such as `drugsVsPloidyCorr.xlsx` remain for backward compatibility, the pipeline must also write metadata/checksums proving which category mode generated them.

### 7. Add Validation Gates

Add hard failures before enrichment:

- No `NA`, empty string, or whitespace-only final categories in included curated rows.
- No correlation-eligible drug missing from `drug_class_final_curated.tsv`.
- No drug with `curation_status != "approved"` in manuscript mode.
- No final category outside an allowed category schema.
- No reviewed exclusion without `exclusion_reason`.
- No `include_in_enrichment = TRUE` row with an excluded or manuscript-disallowed schema category.
- No duplicate `drug_key` without an explicit resolved conflict.
- Plot-ordering code uses schema order and does not require hard-coded `SIGNALING`/`CYTOTOXIC` anchors.
- No duplicate aliases assigned to multiple schema categories.
- No invalid boolean fields in schema or curated mapping.
- No missing or duplicated `plot_order` among manuscript-visible categories.
- No `collapse_to_parent` policy without a valid `collapse_parent_category_id`.
- No manuscript run with `category_mode` other than `curated`.
- No broad administrative category used as a manuscript final class.

Add audit outputs:

- `drug_class_assignment_diff.tsv`: old group vs final group, with reason.
- `drug_class_unassigned_failures.tsv`: should have zero rows on successful runs.
- `drug_class_category_counts.tsv`: final category counts before and after low-frequency filtering.
- `drug_class_reviewed_exclusions.tsv`: all approved exclusions with reasons.
- `drug_class_low_count_actions.tsv`: category-level low-count decisions, thresholds, parent collapses, and affected drugs.

### 8. Add Tests To Prevent New Errors

Add tests under `Code/gdsc_ploidy_analysis/tests/` that assert:

- The wrong-assignment fixture drugs resolve to exact expected approved category IDs, not merely “not the old bad label.”
- Multi-category PubChem strings do not resolve by last token alone.
- Multi-category PubChem parsing handles trailing semicolons, empty tokens, and extra whitespace without creating fake categories.
- `(Anti-)Inflammatory` is not automatically mapped to `IMMUNOSUPPRESSIVE AGENTS`.
- PARP inhibitors are not mapped to `SIGNALING`.
- Broad labels such as `Targeted therapies`, `Enzyme Inhibitors`, and `Antineoplastic Agents` do not override specific GDSC target/pathway evidence.
- Every correlation-eligible drug has either a schema-valid approved final category or an approved exclusion.
- Duplicate curated keys fail unless an explicit conflict-resolution field is present and valid.
- Missing curated keys fail before category filtering.
- Unapproved curation statuses fail in manuscript mode.
- Invalid schema categories fail.
- Whitespace-only categories fail.
- Broad administrative categories fail when used as manuscript final classes.
- Missing curated keys write a failure table with the expected drug rows.
- Known validation failures emit stable error messages that tests can assert.
- Rerunning the pipeline writes `drug_class_final_used.tsv`, `class_enrichment_selected_drugs_curated_Z_SCORE.tsv`, and legacy comparison outputs.

The wrong-assignment examples should become a permanent repo fixture, not a local file dependency:

`Code/gdsc_ploidy_analysis/tests/fixtures/wrong_assignment_examples.tsv`

Required fixture columns:

- `drug_key`
- `legacy_bad_group`
- `expected_category_id`
- `expected_include_in_enrichment`
- `expected_exclusion_reason_code`
- `evidence_summary`

Update existing tests that lock in legacy behavior. In particular, tests around `normalize_legacy_drug_group()` should be replaced or re-scoped so they verify raw token parsing and schema-checked resolution rather than the old bad mappings.

### 9. Update Manifests, README, And Regression Baselines

Add the new schema, curated mapping, and wrong-assignment fixture to:

- input manifests written by the pipeline
- `Code/gdsc_ploidy_analysis/README.md`
- baseline creation and comparison scripts
- checksum/regression reports

The regression report must identify the first changed stage: curation input, final category mapping, selected drug table, enrichment matrix, or heatmap rendering.

### 10. Review Expected Output Changes

Changing final drug categories will change enrichment p-values and Figure 1 heatmaps. That is expected. The review should compare:

- old vs new category assignments
- old vs new category counts
- old vs new enrichment matrices
- old vs new Figure 1 heatmaps
- included vs reviewed-excluded drug counts
- category schema and display labels used in Figure 1

Do not approve the change solely because the original listed errors disappear. Approval requires confirming that all new category changes are supported by target/pathway or curated evidence.

## Acceptance Criteria

- `run_gdsc_ploidy_analysis.R` fails if any correlation-eligible drug lacks an approved final category or an approved exclusion.
- `drug_class_category_schema.tsv` exists and controls accepted category IDs, display labels, manuscript eligibility, plot order, and inclusion policy.
- `drug_class_final_curated.tsv` covers the full pre-filter correlation-eligible drug universe.
- `drug_class_final_curated.tsv` has exactly one resolved row per correlation-eligible `drug_key`.
- `drug_class_final_used.tsv` exists after every successful run, contains every correlation-eligible drug, and fully explains included drugs, reviewed exclusions, low-count filtering, schema version, mapping version, and category mode.
- All drugs listed in the wrong-assignment fixture have exact expected approved category IDs or explicit approved exclusions.
- No final assignment depends on taking only the last PubChem category.
- No final assignment is produced only from broad PubChem labels when more specific GDSC target/pathway evidence is available.
- No curated output is labeled as legacy provenance.
- Curated workbook and figure artifact names are distinct from legacy comparison artifacts, or generic compatibility outputs include category-mode metadata and checksums.
- No drug can disappear because of missing category assignment before enrichment.
- Figure 1 heatmaps are regenerated from the corrected enrichment workbook.
