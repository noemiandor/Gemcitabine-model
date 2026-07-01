# Manuscript Figure Output Standardization Plan

## Goal

Make this repository predictable to rerun and audit by enforcing one convention:

- `Code/` contains code and static code-adjacent fixtures only.
- `Data/` contains raw, manual, and reusable processed inputs.
- `Results/` contains generated analysis outputs from reproducible runs.
- `figures/` contains manuscript-facing figure assets and per-figure manifests.

The end state should include a root-level manager entrypoint, preferably `Manager.sh`, that can regenerate every figure panel this repository is currently able to generate. The manager should support multiple rerun levels, including plotting from saved fit summaries without refitting, full refitting where supported, and check-only/dry-run modes.

This plan is intentionally not an implementation. It defines the target layout, dependency order, migration steps, and validation gates.

## End-State User View

After implementation, a root-level `Manager.sh` command should make the current
locally reproducible manuscript panels visible through one interface. The table
below is the intended user-facing contract for what the manager can generate or
materialize into `figures/`.

| Figure | Panels generated | Content | Source module | Default mode | Optional flags | External/manual caveats |
| --- | --- | --- | --- | --- | --- | --- |
| Figure 1 | GDSC enrichment panels and source tables | Low- and high-ploidy reviewed primary drug-class enrichment summaries, including clustered heatmap-style panels where enabled | `Code/gdsc_ploidy_analysis/` | `standard` | `--gdsc-analysis-mode`, `--gdsc-enrichment-permute-n`, `--gdsc-drug-class-workbook` | Final panel lettering/composition may remain manuscript-facing assembly under `figures/Figure1/` |
| Figure 2 | Public-data ploidy support panels | GDSC support outputs and local CCLE ploidy-support output for Fig. 2B when configured | `Code/gdsc_ploidy_analysis/`, `Code/ccle_ploidy_analysis/` | `standard` | `--modules gdsc,ccle`; CCLE metric/source options once finalized | Fig. 2B metric policy must be fixed before freezing; remaining Figure 2 panels may be external/manual |
| Figure 3 | Fig. 3H drug-response/ploidy panels; Fig. 3I LCI overlays when source analysis is provided | Gemcitabine dose-response fits, AUC/EC50/IC50/ploidy association outputs, and live-cell imaging overlay/time-course panels | `Code/in-vitro/drug_response/`, `Code/lci_overlays/` | `standard` | `--refresh-cloneid-ploidy`, `--lci-analysis-dir`, `--lci-render`, `--lci-panel-only` | LCI rendering requires an external analysis directory; earlier Figure 3 panels remain external/manual unless source code is added |
| Figure 4 | Fig. 4B PKPD-derived dFdCTP driver component | Baseline-subtracted PK-derived dFdCTP signal driver plots used to support the PKPD model | `Code/in-vitro/pkpd_live_dead_model/` | `saved-fit` or `standard` | `--pkpd-saved-fit`, `--modules pkpd` | Fig. 4C is external/resolved; other Figure 4 panels may remain external/manual |
| Figure 5 | Fig. 5B-D model panels and optional Fig. 5F metabolomics support | Cohort live/dead fits, dFdCTP signal curves, effective beta/Hill-corrected signal plot, ploidy parameter fold change, dose-response comparison, and dCMP/metabolomics support where available | `Code/in-vitro/pkpd_live_dead_model/`, `Code/Gemcitabine_Metabolomics_Heatmap/` | `saved-fit` for model plots; `standard` for metabolomics | `--pkpd-refit`, `--pkpd-smoke-fit`, `--pkpd-fit-output`, `--metabolomics-input` | Full model refitting is optional and expensive; immunoblot/schematic panels remain external/manual |
| Figure 6 | Metabolomics panels | PCA, response/category summaries, pathway enrichment heatmaps, ordered metabolite heatmaps, and volcano-style outputs | `Code/Gemcitabine_Metabolomics_Heatmap/` | `standard` | `--metabolomics-input`, `--modules metabolomics` | The manager should materialize source panels, while final composite assembly may remain manual |
| Supplementary | GDSC and PKPD supplementary panels | Supplementary enrichment summaries and cohort model-fit/source plots such as Supp. Fig. 1 | GDSC and PKPD modules | `standard` or `saved-fit` | Same module-specific flags as above | Only locally reproducible supplementary panels are in scope |
| In vivo pending | Optional pseudotime/TGI panels | Pseudotime-shift, TGI, and ploidy association plots/tables | `Code/in-vivo/` | Not run by default | `--include-in-vivo` | Results are pending; manuscript use should remain gated by explicit opt-in |

The manager should support `check-only`, `panels-only`, `saved-fit`,
`standard`, and `full-refit` modes. It should also allow scoped module runs with
`--modules`, custom run/output roots, dry runs, overwrite protection, and
explicit control over whether `latest.txt` pointers are updated.

## Non-Goals

- Do not attempt to reproduce manually assembled or external panels unless code and source inputs already exist locally.
- Do not delete old outputs during the initial migration. First generate canonical outputs under `Results/`, update manifests, and only then decide what old `Figs/`, `Code/*/output/`, or `Data/*/invitro_fitting_outputs/` paths can be archived or removed.
- Do not require full PKPD refitting for ordinary manuscript panel regeneration. Saved-summary plotting must remain a first-class workflow.
- Do not make symlinks mandatory. A text pointer such as `Results/.../latest.txt` is more portable than a filesystem symlink.
- Do not switch module default output locations during the first migration steps. Initially, modules should accept explicit canonical output directories while preserving current defaults for backward compatibility. Default changes belong in the final compatibility phase.
- Do not move static module input fixtures until all output contracts are stable. For example, `Code/*/data/` can remain in place initially if code already treats those files as module fixtures; later phases can decide whether they should move under `Data/`.

## Target Layout

```text
Data/
  <domain>/<module>/
    raw/
    manual/
    processed/

Results/
  <domain>/<module>/
    latest.txt
    runs/<run_id>/
      figures/
      tables/
      models/
      metadata/
      logs/
      qc/

figures/
  Figure1/
    manifest.tsv
    panel_1A.png
    panel_1B.png
  Figure2/
    manifest.tsv
  Figure3/
    manifest.tsv
  Figure4/
    manifest.tsv
  Figure5/
    manifest.tsv
  Figure6/
    manifest.tsv
  Supplementary/
    manifest.tsv
```

`run_id` should be deterministic when requested and timestamped otherwise. Recommended default format:

```text
YYYYMMDDTHHMMSS_<short_label>
```

Examples:

```text
20260629T143000_full_pkpd_beta_hill_confluence
20260629T143000_saved_pkpd_alsogoodfit
20260629T143000_gdsc_primary_secondary_zscore
```

Run ID rules:

- `manager_run_id` identifies one root-manager invocation.
- `module_run_id` should normally be `<manager_run_id>_<module_slug>` when a manager invocation launches multiple modules.
- Standalone module runs may use the provided `--run-id` directly.
- Allowed characters should be `[A-Za-z0-9._-]`; spaces and shell-special characters should be rejected.
- Timestamps should be generated in local time unless the manager records `timezone` in `run_config.tsv`; UTC is also acceptable if stated consistently.
- If a run directory already exists, the command must fail unless `--overwrite` is explicitly provided.

## Latest Pointer Contract

Each module may maintain:

```text
Results/<domain>/<module>/latest.txt
```

`latest.txt` should be a small tab-separated file, not a free-form path. Required fields:

```text
key	value
run_id	<module_run_id>
module	<module_name>
mode	<mode_that_created_the_run>
run_path	<repo-relative Results/.../runs/<run_id> path>
command_sha256	<sha256 of normalized command string>
output_manifest_sha256	<sha256 of metadata/output_manifest.tsv>
updated_at	<ISO-8601 timestamp>
```

Rules:

- Update `latest.txt` atomically only after the module run and validation have succeeded.
- Never update it for failed, interrupted, dry-run, or check-only runs.
- Add a manager-level `--no-update-latest` option for validation runs.
- Resolve `latest.txt` to a concrete run path before writing any `figures/FigureX/manifest.tsv`; figure manifests should not point to mutable `latest.txt`.
- Consider adding `approved.txt` later for manuscript-approved frozen runs. `latest.txt` should mean only "last successful validated run", not "approved for the manuscript".

## Per-Run Required Files

Each run directory should contain:

- `metadata/run_config.tsv`: command, script version, module, key arguments, run_id, output root.
- `metadata/input_manifest.tsv`: every required input path, role, source kind, checksum, byte size, and whether it is raw/manual/processed.
- `metadata/output_manifest.tsv`: every generated output path, role, source kind, figure-panel mapping if applicable, checksum, and byte size.
- `logs/stdout.log` and `logs/stderr.log` for manager-launched commands.

Existing module-specific `metadata/`, `tables/`, and `qc/` conventions should be preserved where possible, but the root output directory should move under `Results/`.

Minimum manifest columns:

```text
path
repo_relative_path
absolute_path
role
source_kind
module
generated_by
command_id
sha256
checksum_unavailable_reason
byte_size
mtime_utc
figure
panel
notes
```

`sha256` should be filled for local files whenever feasible. If a checksum is impossible or inappropriate, `checksum_unavailable_reason` must explain why. Absolute paths may be included for audit trails, but repo-relative paths should be the primary stable reference whenever the file lives inside the repository.

Add a manifest validator before treating a run as successful. Validation should check required headers, nonempty required fields, path existence for local files, checksum consistency where checksums are available, and that output paths are inside the intended result directory.

## Manuscript Figure Manifest

Each `figures/FigureX/manifest.tsv` should describe what the manuscript-facing assets are and where they came from.

Required columns:

```text
figure
panel
asset_path
source_file
source_kind
generated_by
command
input_data
result_run_dir
run_id
caption_role
asset_status
not_regenerated_reason
local_provenance_path
citation_or_uri
notes
```

Recommended `source_kind` values:

- `generated_panel`: copied or exported from a reproducible `Results/.../runs/<run_id>/figures/` file.
- `generated_table`: table used directly by a panel or caption.
- `manual_composite`: assembled outside code, but tracked as the manuscript asset.
- `external`: cited or externally assembled and not regenerated by this repo.
- `placeholder`: temporary manuscript asset that should not be considered final.

The manifest should point first to the manuscript-facing `figures/FigureX/...` asset, then to the canonical source under `Results/.../runs/<run_id>/...`.

Manual or external panels should be first-class manifest rows. They may leave `result_run_dir` and `run_id` empty, but must include `not_regenerated_reason`. If there is local provenance, such as `Figs/GemcitabinePaper_Figures.pptx`, record it in `local_provenance_path` and include a checksum for any local manuscript-facing asset. If the panel is external, record the citation key or URI in `citation_or_uri`.

## Dependency Graph By Module

### GDSC Ploidy Analysis: Figure 1, Figure 2A Support, GDSC Supplement

Current entrypoint:

```bash
Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R
```

Current default output:

```text
Code/gdsc_ploidy_analysis/output/
```

Target output:

```text
Results/public_data/gdsc_ploidy_analysis/runs/<run_id>/
```

Important dependencies:

1. `Code/gdsc_ploidy_analysis/src/run_gdsc_ploidy_analysis.R` reads raw GDSC, ploidy, and the reviewed primary/secondary drug-class workbook from `Code/gdsc_ploidy_analysis/data/`.
2. The R pipeline writes enrichment workbooks such as `drugsVsPloidyCorr_<category>_<metric>.xlsx`.
3. The same R pipeline immediately calls Python plotting scripts:
   - `Code/gdsc_ploidy_analysis/src/plot_ploidy_enrichment_panels.py`
   - `Code/gdsc_ploidy_analysis/src/plot_ploidy_enrichment_clustered_heatmaps.py`
4. Those Python scripts consume the workbook written in step 2 and generate Figure 1 source panels and exploratory clustered heatmaps.

Implementation implications:

- During migration, the manager should pass `--output-dir=Results/public_data/gdsc_ploidy_analysis/runs/<run_id>/`. The script's current default can remain `Code/gdsc_ploidy_analysis/output/` until the final compatibility phase.
- The R-to-Python internal call must pass the same canonical result directory.
- The compatibility workbook `drugsVsPloidyCorr.xlsx` can remain inside the run directory, but the canonical workbook should be the category/metric-specific filename.
- The manuscript-oriented manager command uses `--analysis-mode=manuscript` and the reviewed primary/secondary drug-class workbook. Current code uses `Z_SCORE` as the canonical metric; any Figure 1/Figure 2A manuscript wording that refers to IC50 must be reconciled before freezing the manuscript-facing assets.
- The manager should copy selected assets into:
  - `figures/Figure1/`
  - `figures/Supplementary/`
- The manager should not use `Code/gdsc_ploidy_analysis/output/` once this migration is complete.

Recommended manager modes:

- `--gdsc-only`: rerun GDSC and regenerate Figure 1/SI source panels.
- `--skip-gdsc`: use `Results/public_data/gdsc_ploidy_analysis/latest.txt`.
- `--drug-class-workbook <path>` for the reviewed primary/secondary class workbook.
- `--analysis-mode dev|manuscript`; manuscript figure regeneration should use `manuscript`.
- `--gdsc-enrichment-permute-n N`, mapped to the current `--enrichment-permute-n` script option.

### CCLE Ploidy Analysis: Figure 2B Local Support

Current entrypoint:

```bash
Rscript Code/ccle_ploidy_analysis/run_ccle_ploidy_analysis.R
```

Current default output:

```text
Code/ccle_ploidy_analysis/output/
```

Target output:

```text
Results/public_data/ccle_ploidy_analysis/runs/<run_id>/
```

Important dependencies:

1. The CCLE script reads raw and manual inputs under `Code/ccle_ploidy_analysis/data/`.
2. It writes correlation tables, a barplot, and metadata.
3. No other local script currently consumes CCLE outputs as input, but Figure 2B local support and documentation do.

Implementation implications:

- Keep `Code/ccle_ploidy_analysis/data/` for now if those files are static module fixtures. A later cleanup can decide whether they should move to `Data/public_data/ccle_ploidy_analysis/`.
- During migration, the manager should pass the canonical `--output-dir` explicitly while preserving the current module default until the final compatibility phase.
- Decide which Figure 2B support mode is intended:
  - `Z_SCORE` support reproduces the current local default.
  - `--metric=IC50 --metric-source=grbrowser` is the relevant mode if the manuscript intends true IC50 support.
  - If Figure 2B remains external-only, the local output should be recorded as support/provenance but not copied as a manuscript panel.
- The manager should copy or reference `ccle_drug_ploidy_correlations_<metric>.pdf` into `figures/Figure2/` if the manuscript continues using the local support plot; otherwise manifest it as local support while marking the panel as external.

### Figure 3H Drug-Response/Ploidy Analysis

Current entrypoints:

```bash
Rscript Code/in-vitro/drug_response/query_cloneid_fig3h_ploidy.R
Rscript Code/in-vitro/drug_response/plot_gemcitabine_ploidy_auc_association.R
```

Current outputs:

```text
Data/in-vitro/drug_response/fig3h_cloneid_lineage_map.tsv
Data/in-vitro/drug_response/fig3h_cloneid_ploidy.tsv
Data/in-vitro/drug_response/fig3h_cloneid_ploidy_profile_details.tsv
Data/in-vitro/drug_response/fig3h_cloneid_id_mapping.tsv
Figs/gemcitabine_*.png
Figs/gemcitabine_*.pdf
Figs/gemcitabine_*.tsv
```

Target outputs:

```text
Data/in-vitro/drug_response/manual/fig3h_cloneid_lineage_map.tsv
Data/in-vitro/drug_response/processed/fig3h_cloneid_ploidy.tsv
Data/in-vitro/drug_response/processed/fig3h_cloneid_ploidy_profile_details.tsv
Data/in-vitro/drug_response/processed/fig3h_cloneid_id_mapping.tsv
Results/in-vitro/drug_response/runs/<run_id>/{figures,tables,metadata,logs}/
```

Important dependencies:

1. `query_cloneid_fig3h_ploidy.R` depends on `cloneid`, DB access, Java support, and `fig3h_cloneid_lineage_map.tsv`. It currently writes three explicit output files via `--output`, `--details-output`, and `--mapping-output`.
2. `plot_gemcitabine_ploidy_auc_association.R` consumes:
   - `Data/in-vitro/drug_response/Gemcitabine.txt`
   - `fig3h_cloneid_ploidy.tsv`
3. The plotting script writes dose-response fits, AUC/EC50/IC50 correlations, and source plots.

Implementation implications:

- Split this workflow into two manager stages:
  - `prepare_fig3h_ploidy`: optional, requires external CLONEID infrastructure.
  - `plot_fig3h_drug_response`: local, consumes an existing ploidy table.
- Do not treat an unavailable CLONEID/DB environment as a successful validation. If unavailable, the manager should explicitly say it is reusing the existing processed ploidy table.
- Either preserve the query script's current three explicit output flags or add a new `--output-dir` that writes all three processed outputs consistently. Until such an option exists, manager commands should pass all three output paths explicitly.
- Add `--output-dir` named option support to `plot_gemcitabine_ploidy_auc_association.R`; keep positional arguments temporarily for backward compatibility.
- Move plotting outputs from top-level `Figs/` to `Results/in-vitro/drug_response/runs/<run_id>/figures/` and `tables/`.
- Copy only selected manuscript-facing panels into `figures/Figure3/`, with manifest rows for each.

Recommended manager modes:

- `--refresh-cloneid-ploidy`: rerun `query_cloneid_fig3h_ploidy.R`.
- Default: reuse existing processed ploidy table and rerun plotting.
- `--skip-fig3h`: use latest existing run.

### Live-Cell Imaging Overlays: Figure 3I

Current entrypoint:

```bash
python3 Code/lci_overlays/generate_timecourse_overlays.py
```

Current default output:

```text
<analysis-dir>/Overlays/<well>_timecourse/
```

Target output:

```text
Results/in-vitro/lci_overlays/runs/<run_id>/figures/
```

Important dependencies:

1. `generate_timecourse_overlays.py` calls `generate_class_colored_overlay.py` once per requested time point.
2. It then assembles the generated per-time overlays into a time-course panel.
3. The input `analysis-dir` may be a large external or local image-analysis folder. This dependency needs to be explicit in `input_manifest.tsv`.

Implementation implications:

- Manager should require an explicit `--lci-analysis-dir` unless a checked-in default is truly valid.
- Generated overlays and assembled panels should land under `Results/in-vitro/lci_overlays/runs/<run_id>/figures/`.
- The Figure 3 manifest should distinguish individual overlays from the final assembled source panel.
- `--skip-render` should remain available for rebuilding the assembled panel from already rendered overlays.
- Do not rely on the current module-local `Code/lci_overlays/Manager.sh` unless it is repaired first; it refers to `scripts/generate_timecourse_overlays.py`, while the active script is `Code/lci_overlays/generate_timecourse_overlays.py`. The root manager should call the active Python script directly until the module-local wrapper is fixed or retired.

Recommended manager modes:

- `--lci-render`: render all overlays.
- `--lci-panel-only`: use existing overlays and rebuild the time-course panel.
- Default: skip LCI unless the required analysis directory exists and is explicitly provided.

### PKPD Live/Dead Model: Figure 4B dFdCTP Component, Figure 5B-D, Supplementary Figure 1

Current entrypoints:

```bash
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py
```

Current inputs and outputs:

```text
Data/GemDelayKillTerm/processed/counts_by_well_time_wellAggregated.parquet
Data/in-vitro/pkpd_live_dead_model/raw/Gemcitabine_PlateMap_20240111.xlsx
Data/in-vitro/pkpd_live_dead_model/raw/drugKinetics/GemcitabineExposure_PKPD.xlsx
Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/<fit_id>/
```

Target outputs:

```text
Data/GemDelayKillTerm/processed/
Results/in-vitro/pkpd_live_dead_model/runs/<run_id>/
  figures/
  tables/
  models/
  metadata/
  logs/
```

Important dependencies:

1. `Data/GemDelayKillTerm/portScript.R` currently generates the count parquet files consumed by the fitter/plotter. This is code living under `Data/`, which violates the target convention.
2. `run_invitro_fit.py` consumes processed counts, platemap, and PKPD workbook. It writes `joint_fit_summary.tsv` and `optimizer_attempts.tsv`.
3. `plot_invitro_fit_outputs.py` consumes a fit-output folder containing `joint_fit_summary.tsv`. It reconstructs the model configuration and regenerates downstream plots.
4. Saved-summary plotting must work both for archived saved fits and for newly generated full refits.

Implementation implications:

- Move or wrap `Data/GemDelayKillTerm/portScript.R` into an appropriate code location, such as `Code/in-vitro/pkpd_live_dead_model/preprocess_counts.R`, while keeping generated parquets under `Data/GemDelayKillTerm/processed/`.
- During migration, manager commands should pass `--output-dir` to `run_invitro_fit.py`; the current default under `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/<timestamp>` can remain until the final compatibility phase.
- Make PKPD `--check-inputs` write its validation manifest under the requested output directory, not directly under `Data/in-vitro/pkpd_live_dead_model/`, before relying on it in manager validation.
- Keep compatibility copies of `joint_fit_summary.tsv` and `optimizer_attempts.tsv` at the run root until the plotter, tests, and validation commands explicitly support `models/joint_fit_summary.tsv` and `models/optimizer_attempts.tsv`. Moving summaries only into `models/` before that would break saved-summary plotting.
- Add an explicit input/output split to `plot_invitro_fit_outputs.py` before using it from the root manager:
  - `--fit-output <path-or-id>` or `--fit-summary-dir <path>` for the folder containing `joint_fit_summary.tsv`.
  - `--output-dir <Results/.../runs/<run_id>>` for generated figures and tables.
  - Legacy positional behavior should remain valid during migration and should continue writing back into the input fit folder only when `--output-dir` is omitted.
- Preserve a mode where the plotter reads saved summaries from legacy paths under `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/` during migration.
- Record the saved summary path and checksum in the PKPD run's `input_manifest.tsv` when plotting from an existing saved fit.

Recommended manager modes:

- `--pkpd-saved-fit <path-or-id>`: regenerate plots from an existing `joint_fit_summary.tsv` without fitting.
- `--pkpd-refit`: run full fitting, then feed the resulting summary into the plotter.
- `--pkpd-smoke-fit`: run the existing smoke test and validate compatibility, but do not use it for manuscript panels.
- `--pkpd-check-inputs`: validate counts, platemap, PKPD workbook, surfaces, and live/dead alignment.

Minimum acceptance criteria:

- Saved-summary mode reproduces:
  - `cohort_joint_fit_2n.png`
  - `cohort_joint_fit_4n.png`
  - `dfdctp_signal_curve_2n.png`
  - `dfdctp_signal_curve_4n.png`
  - `dfdctp_signal_curve_combined_ploidy.png`
  - `effective_dfdctp_signal_curve_combined_ploidy.png`
  - `ploidy_parameter_log2_fold_change.png`
  - `dose_response_ploidy_comparison.png`
- Full-refit mode produces a new fit summary and then regenerates the same downstream plot set from that summary.

### Metabolomics: Figure 5F Support And Figure 6

Current entrypoints:

```bash
python3 Code/Gemcitabine_Metabolomics_Heatmap/run_full_2fold_metabolomics_analysis.py
python3 Code/Gemcitabine_Metabolomics_Heatmap/Pathway_enrichment_heatmap_Gemcitabine_2fold_corrected.py
python3 Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_Heatmap_Gemcitabine_zscore.py
```

Current output behavior:

- `run_full_2fold_metabolomics_analysis.py` defaults to `metabolomics_2fold_full_package/`.
- `Pathway_enrichment_heatmap_Gemcitabine_2fold_corrected.py` requires `--outdir`.
- `Metabolomics_Heatmap_Gemcitabine_zscore.py` defaults to `ordered_response_heatmap_output/`.

Target output:

```text
Results/in-vitro/metabolomics/runs/<run_id>/
  figures/
  tables/
  metadata/
  logs/
```

Important dependencies:

1. These scripts appear to consume the same metabolomics workbook independently. The workbook currently lives at `Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_2N_4N_Full.xlsm`; the target input location should be `Data/in-vitro/metabolomics/raw/Metabolomics_2N_4N_Full.xlsm`.
2. The pathway-enrichment and z-score heatmap scripts do not currently consume tables generated by `run_full_2fold_metabolomics_analysis.py`.
3. If future refactoring shares preprocessed tables between them, the manager should make that dependency explicit. For the first migration, keep them independent and simply standardize their output roots.

Implementation implications:

- Add or standardize `--output-dir`/`--outdir` behavior across the three scripts.
- Add a manager-level `--metabolomics-input` option, defaulting to the canonical `Data/in-vitro/metabolomics/raw/Metabolomics_2N_4N_Full.xlsm` after the workbook is copied or moved there. During migration, the manager may explicitly point to the current workbook path.
- Write all tables to `Results/in-vitro/metabolomics/runs/<run_id>/tables/`.
- Write all source plots to `Results/in-vitro/metabolomics/runs/<run_id>/figures/`.
- Copy selected Figure 6 source panels into `figures/Figure6/`.
- If a dCMP-specific Figure 5F panel is later added, write it under the same metabolomics run and add a `figures/Figure5/manifest.tsv` row.

### In-Vivo Pseudotime/TGI Analyses: Pending In-Vivo Manuscript Sections

Current entrypoints:

```bash
Rscript Code/in-vivo/pseudotimeAssociations.R
Rscript Code/in-vivo/ploidy_vs_TGI
```

Current output:

```text
Figs/CellCycleCells_*.png
Figs/CellCycleCells_*.pdf
Figs/pseudotimeAssociations_reported_stats.csv
Figs/ploidy_vs_TGI_Day24_p90_treated_by_dose.png
Figs/ploidy_vs_TGI_Day24_p90_treated_by_dose.pdf
```

Target output:

```text
Results/in-vivo/pseudotime_associations/runs/<run_id>/{figures,tables,metadata,logs}/
```

Important dependencies:

1. Both scripts consume `Data/in-vivo/CellCycelCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv`.
2. They do not currently feed downstream figure-generation scripts, but they generate statistics and plots that may be cited in pending in-vivo sections.

Implementation implications:

- Add named `--input` and `--output-dir` support.
- Move statistics CSV files to `tables/`, not `figures/`.
- Write plots to `figures/`.
- Keep these modules optional in the root manager until the in-vivo manuscript figures are finalized.

## Root Manager Design

Create a root-level `Manager.sh` as the primary orchestrator because the repository uses both R and Python. A thin `Manager.R` should stay out of scope unless an R-only orchestration need appears.

Important sequencing rule: do not make the root manager authoritative until module-level output contracts have been stabilized and independently validated. The root manager should consume explicit, tested module CLIs; it should not discover or normalize legacy behavior while those contracts are still changing.

Recommended command surface:

```bash
bash Manager.sh --check-inputs
bash Manager.sh --dry-run
bash Manager.sh --run-id 20260629T143000_manuscript_panels --mode saved-fit
bash Manager.sh --run-id 20260629T143000_full_refit --mode full-refit --pkpd-refit --jobs 6
bash Manager.sh --modules gdsc,ccle,pkpd,metabolomics
```

Recommended manager options:

```text
--run-id <id>
--mode check-only|saved-fit|standard|full-refit|panels-only
--modules <comma-separated modules>
--output-root Results
--figure-root figures
--overwrite
--jobs <N>
--dry-run
--verbose
--no-update-latest

--gdsc-drug-class-workbook <path>
--gdsc-analysis-mode dev|manuscript
--gdsc-enrichment-permute-n <N>

--pkpd-saved-fit <path-or-id>
--pkpd-refit
--pkpd-smoke-fit
--pkpd-fit-output <path>

--refresh-cloneid-ploidy
--lci-analysis-dir <path>
--lci-render
--lci-panel-only

--include-in-vivo
--metabolomics-input <path>
```

Recommended mode semantics:

- `check-only`: validate required inputs and print planned commands. No generated outputs except optional temporary logs.
- `panels-only`: copy/link selected outputs from `Results/.../latest.txt` into `figures/` and write figure manifests. No analysis rerun.
- `saved-fit`: rerun deterministic/local analyses and regenerate PKPD panels from an existing saved fit summary. No full PKPD optimization.
- `standard`: rerun local deterministic analyses, reuse external/cache-dependent inputs, and skip long full refits unless explicitly requested.
- `full-refit`: include long-running PKPD full fitting and feed the new `joint_fit_summary.tsv` back into the plotter.

Mode precedence and constraints:

- `--modules` narrows the set of modules to execute; it should not change a module's mode semantics.
- `--pkpd-refit` is valid only in `full-refit` or an explicit PKPD-only validation mode. It should be rejected in `saved-fit` and `panels-only`.
- `--pkpd-saved-fit` is required for saved-fit PKPD regeneration unless a concrete validated run is resolved from `latest.txt`.
- `--lci-render` and `--lci-panel-only` are mutually exclusive.
- `panels-only` must fail if the required latest/approved result pointer is missing or if the pointed run lacks an output manifest.
- `check-only` and `dry-run` must not update `latest.txt`.
- Missing optional infrastructure, such as CLONEID DB access or an LCI analysis directory, must be reported as skipped or blocked, not silently treated as success.

Module mode matrix:

| Module | check-only | panels-only | saved-fit | standard | full-refit |
|---|---|---|---|---|---|
| GDSC | validate inputs and command | materialize from latest GDSC run | rerun primary-secondary manuscript GDSC | rerun primary-secondary manuscript GDSC | same as standard |
| CCLE | validate inputs and command | materialize/support from latest CCLE run | rerun selected CCLE support mode | rerun selected CCLE support mode | same as standard |
| Figure 3H | validate Gemcitabine table and ploidy table | materialize from latest drug-response run | rerun plotting from existing ploidy table | rerun plotting from existing ploidy table | same as standard unless `--refresh-cloneid-ploidy` is explicitly requested |
| LCI overlays | validate explicit analysis dir if provided | materialize from latest LCI run | skip unless requested | skip unless `--lci-analysis-dir` is provided | same as standard |
| PKPD | validate counts, platemap, PKPD workbook, and fit-summary path | materialize from latest PKPD plot run | plot from saved fit summary; no optimizer | plot from saved fit summary; no optimizer | run full fit, then plot from new summary |
| Metabolomics | validate workbook and commands | materialize from latest metabolomics run | rerun source plots from workbook | rerun source plots from workbook | same as standard |
| In-vivo | validate CSV if `--include-in-vivo` | materialize from latest in-vivo run if included | skip unless included | skip unless included | same as standard |

The manager should write a top-level run manifest:

```text
Results/manager/runs/<run_id>/metadata/module_runs.tsv
```

Suggested columns:

```text
module
status
mode
command
result_run_dir
latest_pointer
started_at
finished_at
stdout_log
stderr_log
notes
```

## Figure Materialization Rules

The manager should materialize manuscript-facing assets after module runs complete.

For each figure:

1. Create `figures/FigureX/`.
2. Copy selected source outputs from `Results/.../runs/<run_id>/figures/`.
3. Preserve stable manuscript-facing filenames such as `panel_1A_gdsc_workflow.png`.
4. Write or update `figures/FigureX/manifest.tsv`.
5. Leave manual/external panels as manifest rows with `source_kind=manual_composite` or `source_kind=external`, and set `asset_status` separately.

Do not make `figures/` a general dumping ground. It should contain only files intended to be cited or assembled into manuscript figures.

`figures/` may contain both source panels and final Overleaf composites, but filenames should make the role explicit. For example:

```text
figures/Figure5/panel_5B_cohort_joint_fit_2n.png
figures/Figure5/panel_5C_effective_dfdctp_signal.png
figures/Figure5/Figure5_v4_Overleaf.png
```

The manifest should distinguish source panels from final composites using `source_kind` and `caption_role`.

## Update `FigureCodeMap.md`

After the output migration, update `docs/FigureCodeMap.md` so each panel follows this order:

1. Manuscript-facing asset under `figures/FigureX/`.
2. Canonical generated source under `Results/.../runs/<run_id>/`.
3. Generating command.
4. Input data under `Data/...`.
5. Remaining manual/external caveats.

The “Output found” column should not point first to `Code/.../output/`, `Figs/`, or legacy `Data/.../invitro_fitting_outputs/` once canonical outputs exist.

## Implementation Phases

### Phase 1: Define Output Contract And Compatibility Helpers

Deliverables:

- Define the canonical result directory helper, run ID rules, latest-pointer writer, and manifest schemas.
- Add or document a manifest validator before requiring every module to use it.
- Create a module registry document or machine-readable table listing module name, domain, current entrypoint, canonical output root, required inputs, expected outputs, and whether the module is optional.
- Do not create the root `Manager.sh` as an authoritative runner yet.
- Do not switch module defaults yet.

Validation:

```bash
git diff --check -- docs/manuscript_figure_output_standardization_plan.md
```

If helper code is added in this phase, validate it with a temporary output root such as `/tmp/gemcitabine_results_contract_check` or `Results/_validation/<run_id>/`; do not mutate real `latest.txt` pointers.

### Phase 2: Add Explicit Canonical Output Options

Update modules so they can write to canonical `Results/.../runs/<run_id>/` directories when explicitly requested. Preserve current defaults until the final compatibility phase.

Suggested order:

1. GDSC: already accepts `--output-dir`; validate canonical output with explicit `--analysis-mode=manuscript` and the reviewed primary/secondary workbook.
2. CCLE: already accepts `--output-dir`; decide whether Figure 2B support uses legacy Z-score, GR-browser IC50, or external-only provenance.
3. PKPD plotter: add separate saved-fit input and output destination, for example `--fit-output <path-or-id>` plus `--output-dir <Results/.../runs/<run_id>>`; preserve current positional folder behavior.
4. PKPD fitter: keep root-level `joint_fit_summary.tsv` and `optimizer_attempts.tsv` compatibility while allowing explicit `--output-dir` under `Results/`.
5. Figure 3H query: either keep explicit `--output`, `--details-output`, and `--mapping-output` or add `--output-dir`; record the lineage map as a manual input.
6. Figure 3H plotter: add named `--input`, `--ploidy-file`, and `--output-dir`; preserve current positional arguments.
7. In-vivo scripts: add named `--input` and `--output-dir`; split statistics into `tables/` and plots into `figures/`.
8. Metabolomics scripts: standardize `--output-dir`/`--outdir` and output subfolders; support an explicit `--metabolomics-input` path in the manager plan.
9. LCI overlays: keep `--output-dir`; repair or retire the stale module-local `Manager.sh`, or have the future root manager call `generate_timecourse_overlays.py` directly.

Validation should use a disposable run root and `--no-update-latest` once that flag exists. Representative commands after the relevant CLI changes:

```bash
Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R \
  --analysis-mode=manuscript \
  --drug-class-workbook=Code/gdsc_ploidy_analysis/data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx \
  --output-dir=Results/_validation/public_data/gdsc_ploidy_analysis/runs/test_gdsc

Rscript Code/ccle_ploidy_analysis/run_ccle_ploidy_analysis.R \
  --metric=IC50 \
  --metric-source=grbrowser \
  --output-dir=Results/_validation/public_data/ccle_ploidy_analysis/runs/test_ccle_ic50

python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py \
  --check-inputs \
  --output-dir Results/_validation/in-vitro/pkpd_live_dead_model/runs/test_check_inputs

python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py \
  --fit-output Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/alsoGoodFit_20260514T093906 \
  --output-dir Results/_validation/in-vitro/pkpd_live_dead_model/runs/test_saved_fit

Rscript Code/in-vitro/drug_response/query_cloneid_fig3h_ploidy.R --check-config

Rscript Code/in-vitro/drug_response/plot_gemcitabine_ploidy_auc_association.R \
  Data/in-vitro/drug_response/Gemcitabine.txt \
  Results/_validation/in-vitro/drug_response/runs/test_fig3h \
  Data/in-vitro/drug_response/fig3h_cloneid_ploidy.tsv

Rscript Code/in-vivo/pseudotimeAssociations.R \
  --input Data/in-vivo/CellCycelCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv \
  --output-dir Results/_validation/in-vivo/pseudotime_associations/runs/test_pseudotime

Rscript Code/in-vivo/ploidy_vs_TGI \
  --input Data/in-vivo/CellCycelCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv \
  --output-dir Results/_validation/in-vivo/pseudotime_associations/runs/test_ploidy_tgi

python3 Code/Gemcitabine_Metabolomics_Heatmap/run_full_2fold_metabolomics_analysis.py \
  --input Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_2N_4N_Full.xlsm \
  --out Results/_validation/in-vitro/metabolomics/runs/test_full

python3 Code/Gemcitabine_Metabolomics_Heatmap/Pathway_enrichment_heatmap_Gemcitabine_2fold_corrected.py \
  --input Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_2N_4N_Full.xlsm \
  --outdir Results/_validation/in-vitro/metabolomics/runs/test_pathway

python3 Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_Heatmap_Gemcitabine_zscore.py \
  --input Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_2N_4N_Full.xlsm \
  --outdir Results/_validation/in-vitro/metabolomics/runs/test_zscore

python3 Code/lci_overlays/generate_timecourse_overlays.py \
  --analysis-dir <analysis-dir> \
  --output-dir Results/_validation/in-vitro/lci_overlays/runs/test_lci/figures \
  --skip-render
```

The in-vivo named-argument commands, PKPD `--fit-output`/`--output-dir` split, and any future Figure 3H query `--output-dir` are target commands; validation should use current positional or explicit-file behavior until those CLIs are implemented. If `query_cloneid_fig3h_ploidy.R --check-config` fails because CLONEID/DB/Java infrastructure is unavailable, record that explicitly and validate downstream plotting with the existing processed ploidy table instead.

### Phase 3: Per-Module Canonical Runs

Run each module independently into `Results/_validation/...` or a named `Results/.../runs/<run_id>/` directory using explicit output paths. Do not copy anything into `figures/` yet.

For each module:

- Confirm expected figures and tables exist.
- Confirm no generated output was written to `Code/*/output/` or top-level `Figs/`.
- Confirm the module can be rerun with `--overwrite` or fails clearly without it.
- Confirm optional external dependencies are reported explicitly.

This phase should pass module by module before any aggregate manager command is considered.

### Phase 4: Module Output Manifests And Validation

- Add `metadata/output_manifest.tsv` for each module run.
- Add a reusable manifest validator.
- Validate both input and output manifests for schema, local path existence, checksums where available, and intended output-root containment.
- Do not create manuscript-facing `figures/FigureX/manifest.tsv` yet.

Validation:

```bash
python3 Code/tools/validate_manifest.py Results/<domain>/<module>/runs/<run_id>/metadata/input_manifest.tsv
python3 Code/tools/validate_manifest.py Results/<domain>/<module>/runs/<run_id>/metadata/output_manifest.tsv
```

The validator path is illustrative; implement it wherever shared repository utilities are placed.

### Phase 5: Figure Materialization And Figure Manifests

- Create `figures/FigureX/` directories.
- Copy selected source panels from validated `Results/.../runs/<run_id>/figures/`.
- Write `figures/FigureX/manifest.tsv`.
- Add rows for manual/external panels with `not_regenerated_reason`, local provenance path or citation, and checksums for local assets.
- Resolve `latest.txt` to concrete run paths before writing manifests; do not store mutable latest pointers as source paths.

Validation:

```bash
python3 Code/tools/validate_figure_manifest.py figures/Figure1/manifest.tsv
python3 Code/tools/validate_figure_manifest.py figures/Figure5/manifest.tsv
```

The validator should check required columns and verify that generated-panel rows point to existing `Results/` files and manuscript-facing assets.

### Phase 6: Root Manager Skeleton

Only after Phases 2-5 have stable contracts, create root `Manager.sh`.

First manager deliverables:

- Parse `--run-id`, `--mode`, `--modules`, `--output-root`, `--figure-root`, `--dry-run`, `--check-inputs`, `--overwrite`, and `--no-update-latest`.
- Read a module registry rather than hard-coding every command inline.
- Implement dry-run and check-only behavior.
- Run one or two low-risk modules first, such as CCLE and GDSC, before adding PKPD/metabolomics/LCI.

Validation:

```bash
bash Manager.sh --dry-run --mode check-only --modules gdsc,ccle --run-id 20260629T143000_check
bash Manager.sh --mode check-only --modules gdsc,ccle --run-id 20260629T143000_check --no-update-latest
```

### Phase 7: Manager Standard Mode

Implement aggregate `standard` mode only after module-by-module canonical runs pass.

Standard mode should:

- Run GDSC with reviewed primary-secondary manuscript settings.
- Run the selected CCLE support mode or manifest Figure 2B as external-only.
- Run Figure 3H plotting from the existing processed ploidy table.
- Run PKPD saved-summary plotting from an explicit saved fit.
- Run metabolomics source plots.
- Skip LCI unless `--lci-analysis-dir` is provided.
- Skip in-vivo unless `--include-in-vivo` is provided.

Validation:

```bash
bash Manager.sh --mode standard --run-id 20260629T143000_standard --overwrite --no-update-latest
```

Expected outputs:

- `Results/public_data/gdsc_ploidy_analysis/runs/20260629T143000_standard_gdsc/`
- `Results/public_data/ccle_ploidy_analysis/runs/20260629T143000_standard_ccle/`
- `Results/in-vitro/drug_response/runs/20260629T143000_standard_drug_response/`
- `Results/in-vitro/pkpd_live_dead_model/runs/20260629T143000_standard_pkpd_saved_fit/`
- `Results/in-vitro/metabolomics/runs/20260629T143000_standard_metabolomics/`
- `Results/manager/runs/20260629T143000_standard/metadata/module_runs.tsv`
- validated `figures/FigureX/manifest.tsv` files for regenerated panels.

### Phase 8: Full-Refit Mode

Implement manager `full-refit` mode for PKPD:

1. Run `run_invitro_fit.py` into `Results/in-vitro/pkpd_live_dead_model/runs/<manager_run_id>_pkpd_fit/`.
2. Validate root-level compatibility copies of `joint_fit_summary.tsv` and `optimizer_attempts.tsv`, or validate the new `models/` layout only after the plotter supports it.
3. Feed that concrete fit-output path into `plot_invitro_fit_outputs.py`.
4. Materialize Figure 5 and Supplementary Figure 1 assets from the refit outputs.

Validation:

```bash
bash Manager.sh --mode full-refit --pkpd-refit --jobs 6 --run-id 20260629T143000_full_refit
```

If full refitting cannot run because of compute time, missing dependencies, or infrastructure constraints, the manager must report that explicitly and should not mark full-refit validation as passed. Standard manuscript-panel regeneration should remain valid without full refitting.

### Phase 9: Switch Defaults And Retire Legacy Output Locations

After the manager can regenerate current panels:

- Change module defaults from `Code/*/output` and `Figs` to canonical `Results/.../runs/<run_id>/`.
- Keep backward-compatible `--output-dir` overrides.
- Stop committing generated output under `Code/`.
- Decide whether old `Figs/` should be archived or replaced by lowercase `figures/`.
- Move static module input fixtures into `Data/` only where doing so is low-risk and improves clarity.
- Update `GemcitabinePaper.tex` only after manuscript-facing `figures/` assets are stable.

Validation:

```bash
git grep -n "Code/.*/output\\|Figs/" -- Code docs GemcitabinePaper.tex
```

Any remaining references should be intentional legacy/archive notes or explicit manuscript-asset references pending final LaTeX migration.

## Acceptance Criteria

The migration is complete when:

1. A clean run of `bash Manager.sh --mode standard --run-id <id>` regenerates all locally reproducible source panels.
2. A clean run of `bash Manager.sh --mode panels-only --run-id <id>` repopulates `figures/` from existing `Results/` without rerunning analyses.
3. A clean run of `bash Manager.sh --mode full-refit --pkpd-refit --jobs 6 --run-id <id>` produces a new PKPD `joint_fit_summary.tsv`, regenerates downstream PKPD plots, and updates Figure 5/Supplementary manifests.
4. No default module writes generated outputs under `Code/`.
5. No default module writes analysis outputs directly into top-level `Figs/`.
6. `FigureCodeMap.md` points to `figures/FigureX/` and `Results/.../runs/<run_id>/` paths for every locally reproducible panel.
7. Every generated manuscript-facing asset has a manifest row documenting its source command, inputs, run_id, and status.
