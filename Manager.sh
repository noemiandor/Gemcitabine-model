#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${repo_root}"

run_id="$(date +"%Y%m%dT%H%M%S_manuscript")"
source_run_id=""
mode="full-refit"
modules="gdsc,ccle,drug_response,pkpd,metabolomics,si_figure4,in_vivo_figure7"
output_root="Results"
figure_root="figures"
module_registry="docs/manuscript_figure_module_registry.tsv"
overwrite=false
dry_run=false
no_update_latest=false
jobs=1

gdsc_analysis_mode="manuscript"
gdsc_enrichment_permute_n="1000"
gdsc_drug_class_workbook="Code/gdsc_ploidy_analysis/data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx"

pkpd_saved_fit="Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/alsoGoodFit_20260514T093906"
pkpd_fit_model_preset="beta_hill_baseline_confluence"
pkpd_refit=false
pkpd_smoke_fit=false
pkpd_fit_output=""

refresh_cloneid_ploidy=false
lci_analysis_dir=""
lci_render=false
lci_panel_only=false
include_in_vivo=false
metabolomics_input="Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_2N_4N_Full.xlsm"
figure7_full_analysis=false
figure7_refresh_inputs="auto"
figure7_intermediate_dir=""
figure7_python=""
figure7_overwrite_intermediates=false
figure7_cell_ploidy_input="Data/in-vivo/all_ploidy.tsv"
figure7_sample_info_input="Data/in-vivo/sample_info.xlsx"
figure7_growth_curve_input="Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx"
figure7_panels_ae_only=false
figure7_tgi_day="17"
figure7_figure_name="Figure7"
figure7_seurat_rds=""
figure7_gene_set_artifact=""
figure7_reference_id="taoli_state_pathway_etp2_24_day17_v1"
figure7_state_pathway_results_root=""
figure7_canonical_reference_root="${FIGURE7_CANONICAL_REFERENCE_ROOT:-Data/in-vivo/figure7/saved_state_pathway/${figure7_reference_id}}"
figure7_reference_root="${figure7_canonical_reference_root}"
figure7_data_root="${FIGURE7_DATA_ROOT:-Data/in-vivo}"
si_figure4_seurat_metadata=""
si_figure4_scvelo_metrics=""
si_figure4_seurat_metadata_provenance=""
si_figure4_input_source=""
si_figure4_prep_run_dir=""
skip_analysis_loop=false

usage() {
  cat <<'EOF'
Usage: bash Manager.sh [options]

Core options:
  --run-id ID
  --source-run-id ID              Required only for panels-only; immutable source manager run
  --mode check-only|saved-fit|standard|full-refit|panels-only
                                  Default: full-refit
  --modules comma,separated,names
  --output-root DIR
  --figure-root DIR
  --overwrite
  --jobs N
  --dry-run
  --check-inputs                 Alias for --mode check-only
  --no-update-latest

Module options:
  --gdsc-analysis-mode dev|manuscript
  --gdsc-enrichment-permute-n N
  --gdsc-drug-class-workbook PATH
  --pkpd-saved-fit PATH
  --pkpd-refit
  --pkpd-smoke-fit
  --pkpd-fit-output PATH
  --refresh-cloneid-ploidy
  --lci-analysis-dir PATH
  --lci-render
  --lci-panel-only
  --include-in-vivo
  --metabolomics-input PATH
  Module name: si_figure4          Uses published Figure 7 scVelo/Seurat inputs, then run intermediates
  --figure7-full-analysis          Opt-in full pathway recomputation; requires both paths below
  --figure7-refresh-inputs         Run the end-to-end input workflow and publish validated CSVs to Data
  --figure7-intermediate-dir PATH  Isolated full-workflow intermediate directory
  --figure7-python PATH            Python executable containing scVelo dependencies
  --figure7-cell-ploidy-input PATH Cell-level ploidy input for the refresh workflow
  --figure7-sample-info-input PATH Sample metadata workbook for the refresh workflow
  --figure7-growth-curve-input PATH
                                  Tumor growth workbook for the refresh workflow
  --figure7-overwrite-intermediates
                                  Permit replacement of partial full-workflow intermediates
  --figure7-panels-ae-only         Generate/materialize 7A-7E while canonical panel 7F is unavailable
  --figure7-tgi-day DAY            TGI endpoint day (default: 17)
  --figure7-figure-name NAME       Materialization folder under --figure-root (default: Figure7)
  --figure7-seurat-rds ABSOLUTE_PATH
  --figure7-gene-set-artifact PATH  Pinned, versioned local gene-set artifact
  --figure7-state-pathway-results-root PATH
                                  Export panel-7F reference from this completed state-pathway result tree
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-id) run_id="$2"; shift 2 ;;
    --source-run-id) source_run_id="$2"; shift 2 ;;
    --mode) mode="$2"; shift 2 ;;
    --modules) modules="$2"; shift 2 ;;
    --output-root) output_root="$2"; shift 2 ;;
    --figure-root) figure_root="$2"; shift 2 ;;
    --overwrite) overwrite=true; shift ;;
    --jobs) jobs="$2"; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    --check-inputs) mode="check-only"; shift ;;
    --no-update-latest) no_update_latest=true; shift ;;
    --gdsc-analysis-mode) gdsc_analysis_mode="$2"; shift 2 ;;
    --gdsc-enrichment-permute-n) gdsc_enrichment_permute_n="$2"; shift 2 ;;
    --gdsc-drug-class-workbook) gdsc_drug_class_workbook="$2"; shift 2 ;;
    --pkpd-saved-fit) pkpd_saved_fit="$2"; shift 2 ;;
    --pkpd-refit) pkpd_refit=true; shift ;;
    --pkpd-smoke-fit) pkpd_smoke_fit=true; shift ;;
    --pkpd-fit-output) pkpd_fit_output="$2"; shift 2 ;;
    --refresh-cloneid-ploidy) refresh_cloneid_ploidy=true; shift ;;
    --lci-analysis-dir) lci_analysis_dir="$2"; shift 2 ;;
    --lci-render) lci_render=true; shift ;;
    --lci-panel-only) lci_panel_only=true; shift ;;
    --include-in-vivo) include_in_vivo=true; shift ;;
    --metabolomics-input) metabolomics_input="$2"; shift 2 ;;
    --figure7-full-analysis) figure7_full_analysis=true; shift ;;
    --figure7-refresh-inputs) figure7_refresh_inputs=true; shift ;;
    --figure7-intermediate-dir) figure7_intermediate_dir="$2"; shift 2 ;;
    --figure7-python) figure7_python="$2"; shift 2 ;;
    --figure7-cell-ploidy-input) figure7_cell_ploidy_input="$2"; shift 2 ;;
    --figure7-sample-info-input) figure7_sample_info_input="$2"; shift 2 ;;
    --figure7-growth-curve-input) figure7_growth_curve_input="$2"; shift 2 ;;
    --figure7-overwrite-intermediates) figure7_overwrite_intermediates=true; shift ;;
    --figure7-panels-ae-only) figure7_panels_ae_only=true; shift ;;
    --figure7-tgi-day) figure7_tgi_day="$2"; shift 2 ;;
    --figure7-figure-name) figure7_figure_name="$2"; shift 2 ;;
    --figure7-seurat-rds) figure7_seurat_rds="$2"; shift 2 ;;
    --figure7-gene-set-artifact) figure7_gene_set_artifact="$2"; shift 2 ;;
    --figure7-state-pathway-results-root) figure7_state_pathway_results_root="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "${mode}" in
  check-only|saved-fit|standard|full-refit|panels-only) ;;
  *) echo "Invalid --mode: ${mode}" >&2; exit 2 ;;
esac

if [[ "${figure7_refresh_inputs}" == "auto" ]]; then
  if [[ "${mode}" == "full-refit" && "${figure7_full_analysis}" != true && ",${modules}," == *",in_vivo_figure7,"* ]]; then
    figure7_refresh_inputs=true
  else
    figure7_refresh_inputs=false
  fi
fi

if [[ ! "${run_id}" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "Invalid --run-id; use only letters, numbers, dots, underscores, and hyphens" >&2
  exit 2
fi
if [[ -n "${source_run_id}" && ! "${source_run_id}" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "Invalid --source-run-id; use only letters, numbers, dots, underscores, and hyphens" >&2
  exit 2
fi

if [[ "${mode}" == "panels-only" && -z "${source_run_id}" ]]; then
  echo "--source-run-id is required with --mode panels-only" >&2
  exit 2
fi
if [[ "${mode}" != "panels-only" && -n "${source_run_id}" ]]; then
  echo "--source-run-id is valid only with --mode panels-only" >&2
  exit 2
fi
if [[ "${figure7_full_analysis}" == true || -n "${figure7_seurat_rds}" || -n "${figure7_gene_set_artifact}" ]]; then
  if [[ ( "${mode}" != "full-refit" && "${mode}" != "check-only" ) || "${figure7_full_analysis}" != true || -z "${figure7_seurat_rds}" || -z "${figure7_gene_set_artifact}" ]]; then
    echo "--figure7-full-analysis, --figure7-seurat-rds, and --figure7-gene-set-artifact must be supplied together with --mode full-refit or check-only" >&2
    exit 2
  fi
  if [[ "${figure7_seurat_rds}" != /* ]]; then
    echo "--figure7-seurat-rds must be an absolute path" >&2
    exit 2
  fi
fi
if [[ "${figure7_panels_ae_only}" == true && "${figure7_full_analysis}" == true ]]; then
  echo "--figure7-panels-ae-only and --figure7-full-analysis are mutually exclusive" >&2
  exit 2
fi
if [[ "${figure7_refresh_inputs}" == true ]]; then
  if [[ "${mode}" != "full-refit" && "${mode}" != "check-only" ]]; then
    echo "--figure7-refresh-inputs requires --mode full-refit or check-only" >&2
    exit 2
  fi
  if [[ "${figure7_full_analysis}" == true ]]; then
    echo "--figure7-refresh-inputs and --figure7-full-analysis are mutually exclusive" >&2
    exit 2
  fi
  if [[ -z "${figure7_intermediate_dir}" ]]; then
    figure7_intermediate_dir="${output_root}/in-vivo/figure7/intermediates/${run_id}"
  fi
  figure7_intermediate_dir="$(
    python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${figure7_intermediate_dir}"
  )"
fi
if [[ ! "${figure7_tgi_day}" =~ ^[0-9]+$ ]]; then
  echo "--figure7-tgi-day must be a non-negative integer" >&2
  exit 2
fi
if [[ ! "${figure7_figure_name}" =~ ^Figure7([._-][A-Za-z0-9._-]+)?$ ]]; then
  echo "--figure7-figure-name must be Figure7 or a Figure7-prefixed folder name" >&2
  exit 2
fi

if [[ "${pkpd_refit}" == true && "${mode}" != "full-refit" && "${modules}" != "pkpd" ]]; then
  echo "--pkpd-refit is valid only with --mode full-refit or --modules pkpd" >&2
  exit 2
fi
if [[ "${lci_render}" == true && "${lci_panel_only}" == true ]]; then
  echo "--lci-render and --lci-panel-only are mutually exclusive" >&2
  exit 2
fi

IFS=',' read -r -a module_list <<< "${modules}"
if [[ "${include_in_vivo}" == true && ",${modules}," != *",in_vivo,"* ]]; then
  module_list+=("in_vivo")
fi
if [[ "${figure7_refresh_inputs}" == true && ",${modules}," != *",in_vivo_figure7,"* ]]; then
  echo "--figure7-refresh-inputs requires --modules to include in_vivo_figure7" >&2
  exit 2
fi
si_figure4_selected=false
for module in "${module_list[@]}"; do
  if [[ "${module}" == "si_figure4" ]]; then
    si_figure4_selected=true
    break
  fi
done
if [[ "${si_figure4_selected}" == true ]]; then
  if [[ -z "${figure7_intermediate_dir}" ]]; then
    figure7_intermediate_dir="${output_root}/in-vivo/figure7/intermediates/${run_id}"
  fi
  figure7_intermediate_dir="$(
    python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${figure7_intermediate_dir}"
  )"
fi
if [[ -n "${figure7_state_pathway_results_root}" ]]; then
  if [[ "${mode}" == "panels-only" ]]; then
    echo "--figure7-state-pathway-results-root cannot be used with --mode panels-only" >&2
    exit 2
  fi
  if [[ "${figure7_panels_ae_only}" == true ]]; then
    echo "--figure7-state-pathway-results-root cannot be used with --figure7-panels-ae-only" >&2
    exit 2
  fi
  if [[ "${figure7_refresh_inputs}" == true ]]; then
    echo "--figure7-state-pathway-results-root cannot be used with --figure7-refresh-inputs" >&2
    exit 2
  fi
  figure7_module_selected=false
  for module in "${module_list[@]}"; do
    if [[ "${module}" == "in_vivo_figure7" ]]; then
      figure7_module_selected=true
      break
    fi
  done
  if [[ "${figure7_module_selected}" != true ]]; then
    echo "--figure7-state-pathway-results-root requires --modules to include in_vivo_figure7" >&2
    exit 2
  fi
  if [[ ! -d "${figure7_state_pathway_results_root}" ]]; then
    echo "Missing Figure 7 state-pathway results directory: ${figure7_state_pathway_results_root}" >&2
    exit 1
  fi
  figure7_state_pathway_results_root="$(
    python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${figure7_state_pathway_results_root}"
  )"
fi

quote_args() {
  local out=""
  local arg
  for arg in "$@"; do
    printf -v quoted "%q" "${arg}"
    out+="${quoted} "
  done
  printf "%s" "${out% }"
}

require_file() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    return 1
  fi
}

require_dir() {
  local path="$1"
  if [[ ! -d "${path}" ]]; then
    echo "Missing required directory: ${path}" >&2
    return 1
  fi
}

registry_module_name() {
  case "$1" in
    pkpd|pkpd_fit) printf "pkpd_live_dead_model" ;;
    *) printf "%s" "$1" ;;
  esac
}

validate_module_registry() {
  local module
  local registry_name
  require_file "${module_registry}"
  for module in "$@"; do
    registry_name="$(registry_module_name "${module}")"
    if ! awk -F '\t' -v target="${registry_name}" 'NR > 1 && $1 == target { found = 1 } END { exit found ? 0 : 1 }' "${module_registry}"; then
      echo "Module is not listed in ${module_registry}: ${module}" >&2
      return 1
    fi
  done
}

module_run_dir() {
  local module="$1"
  local selected_run_id="${2:-${run_id}}"
  case "${module}" in
    gdsc) printf "%s/public_data/gdsc_ploidy_analysis/runs/%s_gdsc" "${output_root}" "${selected_run_id}" ;;
    ccle) printf "%s/public_data/ccle_ploidy_analysis/runs/%s_ccle" "${output_root}" "${selected_run_id}" ;;
    drug_response) printf "%s/in-vitro/drug_response/runs/%s_drug_response" "${output_root}" "${selected_run_id}" ;;
    lci_overlays) printf "%s/in-vitro/lci_overlays/runs/%s_lci_overlays" "${output_root}" "${selected_run_id}" ;;
    pkpd) printf "%s/in-vitro/pkpd_live_dead_model/runs/%s_pkpd_saved_fit" "${output_root}" "${selected_run_id}" ;;
    pkpd_fit) printf "%s/in-vitro/pkpd_live_dead_model/runs/%s_pkpd_fit" "${output_root}" "${selected_run_id}" ;;
    metabolomics) printf "%s/in-vitro/metabolomics/runs/%s_metabolomics" "${output_root}" "${selected_run_id}" ;;
    metabolomics_pathway) printf "%s/in-vitro/metabolomics/runs/%s_metabolomics_pathway" "${output_root}" "${selected_run_id}" ;;
    metabolomics_zscore) printf "%s/in-vitro/metabolomics/runs/%s_metabolomics_zscore" "${output_root}" "${selected_run_id}" ;;
    in_vivo) printf "%s/in-vivo/pseudotime_associations/runs/%s_in_vivo" "${output_root}" "${selected_run_id}" ;;
    si_figure4) printf "%s/in-vivo/SI_figure4/runs/%s_si_figure4" "${output_root}" "${selected_run_id}" ;;
    in_vivo_figure7) printf "%s/in-vivo/figure7/runs/%s_figure7" "${output_root}" "${selected_run_id}" ;;
    *) echo "Unknown module: ${module}" >&2; return 1 ;;
  esac
}

input_paths_for_module() {
  local module="$1"
  case "${module}" in
    gdsc)
      printf "%s\n" \
        Code/gdsc_ploidy_analysis/data/raw/GDSC2_fitted_dose_response_24Jul22.txt \
        Code/gdsc_ploidy_analysis/data/raw/ploidyAcrossCellLines_V1.txt \
        "${gdsc_drug_class_workbook}"
      ;;
    ccle)
      printf "%s\n" \
        Code/ccle_ploidy_analysis/data/raw/Cell_app_export.txt \
        Code/ccle_ploidy_analysis/data/manual/breast_ccle_ploidy.tsv \
        Code/ccle_ploidy_analysis/data/raw/DrugAliases.txt \
        Code/ccle_ploidy_analysis/data/derived/ccle_expression_columns.tsv
      printf "%s\n" Code/ccle_ploidy_analysis/data/raw/breast_cancer_drug_sensitivity/*.tsv
      ;;
    drug_response)
      printf "%s\n" \
        Data/in-vitro/drug_response/Gemcitabine.txt \
        Data/in-vitro/drug_response/fig3h_cloneid_ploidy.tsv
      ;;
    pkpd|pkpd_fit)
      printf "%s\n" \
        Data/GemDelayKillTerm/processed/counts_by_well_time_wellAggregated.parquet \
        Data/in-vitro/pkpd_live_dead_model/raw/Gemcitabine_PlateMap_20240111.xlsx \
        Data/in-vitro/pkpd_live_dead_model/raw/drugKinetics/GemcitabineExposure_PKPD.xlsx
      if [[ "${module}" == "pkpd" ]]; then
        printf "%s\n" "${pkpd_saved_fit}/joint_fit_summary.tsv"
      fi
      ;;
    metabolomics|metabolomics_pathway|metabolomics_zscore)
      printf "%s\n" "${metabolomics_input}" ;;
    lci_overlays)
      printf "%s\n" "${lci_analysis_dir}" ;;
    in_vivo)
      printf "%s\n" Data/in-vivo/CellCycelCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv ;;
    si_figure4)
      printf "%s\n" \
        Code/in-vivo/figure7/figure7_config.yaml \
        Code/in-vivo/figure7/zenodo_upstream_analysis_documents.tsv \
        "${si_figure4_seurat_metadata}" \
        "${si_figure4_scvelo_metrics}" \
        "${si_figure4_seurat_metadata_provenance}" ;;
    in_vivo_figure7)
      printf "%s\n" Code/in-vivo/figure7/zenodo_upstream_analysis_documents.tsv
      if [[ "${figure7_refresh_inputs}" == true ]]; then
        printf "%s\n" \
          Code/in-vivo/figure7/figure7_config.yaml \
          Code/in-vivo/figure7/zenodo_required_files.tsv \
          "${figure7_cell_ploidy_input}" \
          "${figure7_sample_info_input}" \
          "${figure7_growth_curve_input}"
      else
        printf "%s\n" \
          Code/in-vivo/figure7/figure7_config.yaml \
          Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv \
          Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
        if [[ "${figure7_panels_ae_only}" != true && ( -z "${figure7_state_pathway_results_root}" || -d "${figure7_reference_root}" ) ]]; then
          printf "%s\n" \
            "${figure7_reference_root}/panel_7F_pathway_activity_plot_data.tsv" \
            "${figure7_reference_root}/panel_7F_selected_pathway_gsea.tsv" \
            "${figure7_reference_root}/panel_7F_leading_edge_genes.tsv" \
            "${figure7_reference_root}/state_pathway_gene_ranking_complete.tsv" \
            "${figure7_reference_root}/state_pathway_gsea_complete.tsv" \
            "${figure7_reference_root}/state_pathway_sample_bin_coverage.tsv" \
            "${figure7_reference_root}/state_pathway_design_qc.tsv" \
            "${figure7_reference_root}/state_pathway_provenance.tsv"
        fi
        if [[ "${figure7_full_analysis}" == true ]]; then
          printf "%s\n" "${figure7_seurat_rds}" "${figure7_gene_set_artifact}"
        fi
      fi
      ;;
  esac
}

check_module_inputs() {
  local module="$1"
  local path
  if [[ "${module}" == "si_figure4" && "${si_figure4_input_source}" == "planned_figure7_scvelo_generation" ]]; then
    require_file Code/in-vivo/SI_figure4/generate_supplementary_figure4.R
    require_file Code/in-vivo/figure7/run_figure7.R
    require_file Code/in-vivo/figure7/figure7_config.yaml
    require_file Code/in-vivo/figure7/zenodo_required_files.tsv
    return 0
  fi
  if [[ "${module}" == "in_vivo_figure7" && -n "${figure7_state_pathway_results_root}" ]]; then
    require_dir "${figure7_state_pathway_results_root}"
    require_file Code/in-vivo/figure7/export_state_pathway_reference.R
  fi
  while IFS= read -r path; do
    [[ -z "${path}" ]] && continue
    if [[ "${module}" == "lci_overlays" ]]; then
      require_dir "${path}"
    else
      require_file "${path}"
    fi
  done < <(input_paths_for_module "${module}")
}

command_for_module() {
  local module="$1"
  local run_dir="$2"
  case "${module}" in
    gdsc)
      quote_args Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R \
        "--analysis-mode=${gdsc_analysis_mode}" \
        "--drug-class-workbook=${gdsc_drug_class_workbook}" \
        "--enrichment-permute-n=${gdsc_enrichment_permute_n}" \
        "--output-dir=${run_dir}"
      ;;
    ccle)
      quote_args Rscript Code/ccle_ploidy_analysis/run_ccle_ploidy_analysis.R \
        "--metric=Z_SCORE" \
        "--metric-source=legacy" \
        "--output-dir=${run_dir}"
      ;;
    drug_response)
      quote_args Rscript Code/in-vitro/drug_response/plot_gemcitabine_ploidy_auc_association.R \
        --input Data/in-vitro/drug_response/Gemcitabine.txt \
        --output-dir "${run_dir}" \
        --ploidy-file Data/in-vitro/drug_response/fig3h_cloneid_ploidy.tsv
      ;;
    lci_overlays)
      local extra=()
      if [[ "${lci_panel_only}" == true ]]; then
        extra+=(--skip-render)
      fi
      quote_args python3 Code/lci_overlays/generate_timecourse_overlays.py \
        --analysis-dir "${lci_analysis_dir}" \
        --output-dir "${run_dir}/figures" \
        "${extra[@]}"
      ;;
    pkpd)
      quote_args python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py \
        --fit-output "${pkpd_saved_fit}" \
        --output-dir "${run_dir}"
      ;;
    pkpd_fit)
      local fit_args=(
        python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py
        --output-dir "${run_dir}"
        --n-jobs "${jobs}"
        --model-preset "${pkpd_fit_model_preset}"
      )
      if [[ "${overwrite}" == true ]]; then
        fit_args+=(--overwrite)
      fi
      if [[ "${pkpd_smoke_fit}" == true ]]; then
        fit_args+=(--smoke-test)
      fi
      quote_args "${fit_args[@]}"
      ;;
    metabolomics)
      quote_args python3 Code/Gemcitabine_Metabolomics_Heatmap/run_full_2fold_metabolomics_analysis.py \
        --input "${metabolomics_input}" \
        --output-dir "${run_dir}"
      ;;
    metabolomics_pathway)
      quote_args python3 Code/Gemcitabine_Metabolomics_Heatmap/Pathway_enrichment_heatmap_Gemcitabine_2fold_corrected.py \
        --input "${metabolomics_input}" \
        --output-dir "${run_dir}"
      ;;
    metabolomics_zscore)
      quote_args python3 Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_Heatmap_Gemcitabine_zscore.py \
        --input "${metabolomics_input}" \
        --output-dir "${run_dir}"
      ;;
    in_vivo)
      quote_args Rscript Code/in-vivo/pseudotimeAssociations.R \
        --input Data/in-vivo/CellCycelCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv \
        --output-dir "${run_dir}"
      ;;
    si_figure4)
      quote_args Rscript Code/in-vivo/SI_figure4/generate_supplementary_figure4.R \
        --seurat-metadata "${si_figure4_seurat_metadata}" \
        --scvelo-metrics "${si_figure4_scvelo_metrics}" \
        --seurat-metadata-provenance "${si_figure4_seurat_metadata_provenance}" \
        --config Code/in-vivo/figure7/figure7_config.yaml \
        --output-dir "${run_dir}"
      ;;
    in_vivo_figure7)
      local figure7_mode="standard"
      local figure7_args=()
      if [[ "${figure7_refresh_inputs}" == true ]]; then
        figure7_mode="full-workflow"
        figure7_args=(
          Rscript Code/in-vivo/figure7/run_figure7.R
          "--mode=${figure7_mode}"
          --config=Code/in-vivo/figure7/figure7_config.yaml
          "--intermediate-dir=${figure7_intermediate_dir}"
          "--seurat-metadata-output=${figure7_intermediate_dir}/seurat_metadata.csv"
          "--seurat-metadata-provenance-output=${figure7_intermediate_dir}/seurat_metadata_provenance.tsv"
          "--cell-ploidy-input=${figure7_cell_ploidy_input}"
          "--sample-info-input=${figure7_sample_info_input}"
          "--growth-curve-input=${figure7_growth_curve_input}"
          "--tgi-day=${figure7_tgi_day}"
          "--output-dir=${run_dir}"
        )
        if [[ -n "${figure7_python}" ]]; then
          figure7_args+=("--python=${figure7_python}")
        fi
        if [[ "${figure7_overwrite_intermediates}" == true ]]; then
          figure7_args+=(--overwrite-intermediates=true)
        fi
      else
        figure7_args=(
          Rscript Code/in-vivo/figure7/run_figure7.R
          "--mode=${figure7_mode}"
          --config=Code/in-vivo/figure7/figure7_config.yaml
          --cellcycle-input=Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
          --non-cellcycle-input=Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
          "--tgi-day=${figure7_tgi_day}"
          "--output-dir=${run_dir}"
        )
      fi
      if [[ "${figure7_panels_ae_only}" == true ]]; then
        figure7_args+=(--panel-set=a-e)
      elif [[ "${figure7_refresh_inputs}" != true ]]; then
        figure7_args+=("--saved-state-pathway-dir=${figure7_reference_root}")
        if [[ -n "${figure7_state_pathway_results_root}" ]]; then
          figure7_args+=("--state-pathway-results-root=${figure7_state_pathway_results_root}")
        fi
      fi
      if [[ "${figure7_full_analysis}" == true ]]; then
        figure7_args[2]="--mode=full-analysis"
        figure7_args+=(
          "--seurat-rds=${figure7_seurat_rds}"
          "--gene-set-artifact=${figure7_gene_set_artifact}"
        )
      fi
      quote_args "${figure7_args[@]}"
      ;;
  esac
}

write_latest() {
  local module="$1"
  local run_dir="$2"
  local command_string="$3"
  local latest_path
  latest_path="$(dirname "$(dirname "${run_dir}")")/latest.txt"
  latest_path="$(python3 -c 'import os,sys; print(os.path.normpath(sys.argv[1]))' "${latest_path}")"
  local manifest_path="${run_dir}/metadata/output_manifest.tsv"
  local command_sha output_sha
  command_sha="$(printf "%s" "${command_string}" | shasum -a 256 | awk '{print $1}')"
  output_sha="$(shasum -a 256 "${manifest_path}" | awk '{print $1}')"
  mkdir -p "$(dirname "${latest_path}")"
  {
    printf "key\tvalue\n"
    printf "run_id\t%s\n" "$(basename "${run_dir}")"
    printf "module\t%s\n" "${module}"
    printf "mode\t%s\n" "${mode}"
    printf "run_path\t%s\n" "${run_dir}"
    printf "command_sha256\t%s\n" "${command_sha}"
    printf "output_manifest_sha256\t%s\n" "${output_sha}"
    printf "updated_at\t%s\n" "$(date -Iseconds)"
  } > "${latest_path}.tmp"
  mv "${latest_path}.tmp" "${latest_path}"
}

record_module_run() {
  local module="$1" status="$2" command_string="$3" run_dir="$4" stdout_log="$5" stderr_log="$6" notes="$7" started_at="$8" finished_at="$9"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${module}" "${status}" "${mode}" "${command_string}" "${run_dir}" "" "${started_at}" "${finished_at}" \
    "${stdout_log}" "${stderr_log}" "${notes}" >> "${module_runs_file}"
}

figure7_reference_export_command() {
  quote_args Rscript Code/in-vivo/figure7/export_state_pathway_reference.R \
    "--results-root=${figure7_state_pathway_results_root}" \
    "--output-dir=${figure7_reference_root}"
}

record_figure7_reference_export() {
  local status="$1" command_string="$2" stdout_log="$3" stderr_log="$4" started_at="$5" finished_at="$6"
  local metadata_path="${manager_run_dir}/metadata/figure7_state_pathway_export.tsv"
  mkdir -p "$(dirname "${metadata_path}")"
  {
    printf "key\tvalue\n"
    printf "status\t%s\n" "${status}"
    printf "source_results_root\t%s\n" "${figure7_state_pathway_results_root}"
    printf "source_report_html\t%s\n" "${figure7_state_pathway_results_root}/report/pseudotime_state_pathways_report.html"
    printf "canonical_reference_id\t%s\n" "${figure7_reference_id}"
    printf "exported_reference_dir\t%s\n" "${figure7_reference_root}"
    printf "canonical_data_reference_dir\t%s\n" "${figure7_canonical_reference_root}"
    printf "canonical_data_materialization_metadata\t%s\n" "${manager_run_dir}/metadata/figure7_state_pathway_materialization.tsv"
    printf "exporter_script\t%s\n" "Code/in-vivo/figure7/export_state_pathway_reference.R"
    printf "command\t%s\n" "${command_string}"
    printf "stdout_log\t%s\n" "${stdout_log}"
    printf "stderr_log\t%s\n" "${stderr_log}"
    printf "started_at\t%s\n" "${started_at}"
    printf "finished_at\t%s\n" "${finished_at}"
  } > "${metadata_path}"
}

prepare_figure7_reference() {
  local command_string stdout_log stderr_log started_at finished_at status
  command_string="$(figure7_reference_export_command)"
  stdout_log="${manager_run_dir}/logs/figure7_state_pathway_export.stdout.log"
  stderr_log="${manager_run_dir}/logs/figure7_state_pathway_export.stderr.log"
  mkdir -p "${manager_run_dir}/logs" "$(dirname "${figure7_reference_root}")"
  started_at="$(date -Iseconds)"
  status=0
  bash -c "${command_string}" >"${stdout_log}" 2>"${stderr_log}" || status=$?
  finished_at="$(date -Iseconds)"
  if [[ "${status}" -eq 0 ]]; then
    record_figure7_reference_export "ok" "${command_string}" "${stdout_log}" "${stderr_log}" "${started_at}" "${finished_at}"
  else
    record_figure7_reference_export "failed" "${command_string}" "${stdout_log}" "${stderr_log}" "${started_at}" "${finished_at}"
  fi
  return "${status}"
}

figure7_reference_filenames() {
  printf "%s\n" \
    panel_7F_pathway_activity_plot_data.tsv \
    panel_7F_selected_pathway_gsea.tsv \
    panel_7F_leading_edge_genes.tsv \
    state_pathway_gene_ranking_complete.tsv \
    state_pathway_gsea_complete.tsv \
    state_pathway_sample_bin_coverage.tsv \
    state_pathway_design_qc.tsv \
    state_pathway_provenance.tsv
}

record_figure7_reference_materialization() {
  local status="$1" started_at="$2" finished_at="$3"
  local metadata_path="${manager_run_dir}/metadata/figure7_state_pathway_materialization.tsv"
  {
    printf "key\tvalue\n"
    printf "status\t%s\n" "${status}"
    printf "canonical_reference_id\t%s\n" "${figure7_reference_id}"
    printf "source_reference_dir\t%s\n" "${figure7_reference_root}"
    printf "target_data_reference_dir\t%s\n" "${figure7_canonical_reference_root}"
    printf "materialized_file_count\t8\n"
    printf "started_at\t%s\n" "${started_at}"
    printf "finished_at\t%s\n" "${finished_at}"
  } > "${metadata_path}"
}

materialize_figure7_reference_to_data() {
  local filename source_path target_path temporary_path
  local materialized_count=0

  while IFS= read -r filename; do
    [[ -z "${filename}" ]] && continue
    require_file "${figure7_reference_root}/${filename}"
  done < <(figure7_reference_filenames)

  mkdir -p "${figure7_canonical_reference_root}"
  while IFS= read -r filename; do
    [[ -z "${filename}" ]] && continue
    source_path="${figure7_reference_root}/${filename}"
    target_path="${figure7_canonical_reference_root}/${filename}"
    temporary_path="${figure7_canonical_reference_root}/.${filename}.tmp.${run_id}"
    cp "${source_path}" "${temporary_path}"
    mv "${temporary_path}" "${target_path}"
    if ! cmp -s "${source_path}" "${target_path}"; then
      echo "Materialized Figure 7 reference does not match exported artifact: ${target_path}" >&2
      return 1
    fi
    materialized_count=$((materialized_count + 1))
  done < <(figure7_reference_filenames)

  if [[ "${materialized_count}" -ne 8 ]]; then
    echo "Expected to materialize 8 Figure 7 reference files; found ${materialized_count}" >&2
    return 1
  fi
}

figure7_input_materialized_count=0
figure7_input_source_scvelo=""
figure7_input_source_seurat_metadata=""
figure7_input_source_seurat_metadata_provenance=""
figure7_input_source_cellcycle=""
figure7_input_source_noncellcycle=""
figure7_input_target_scvelo=""
figure7_input_target_seurat_metadata=""
figure7_input_target_seurat_metadata_provenance=""
figure7_input_target_cellcycle=""
figure7_input_target_noncellcycle=""
figure7_input_sha_scvelo=""
figure7_input_sha_seurat_metadata=""
figure7_input_sha_seurat_metadata_provenance=""
figure7_input_sha_cellcycle=""
figure7_input_sha_noncellcycle=""

si_figure4_validate_input_pair() {
  python3 Code/tools/validate_si_figure4_inputs.py \
    --seurat-metadata "$1" \
    --scvelo-metrics "$2" \
    --seurat-metadata-provenance "$3" \
    --config Code/in-vivo/figure7/figure7_config.yaml
}

record_si_figure4_input_materialization() {
  local status="$1"
  local source_seurat="$2"
  local source_scvelo="$3"
  local source_provenance="$4"
  local target_seurat="$5"
  local target_scvelo="$6"
  local target_provenance="$7"
  local metadata_path="${manager_run_dir}/metadata/si_figure4_input_materialization.tsv"
  {
    printf "key\tvalue\n"
    printf "status\t%s\n" "${status}"
    printf "input_source\t%s\n" "${si_figure4_input_source}"
    printf "source_seurat_metadata\t%s\n" "${source_seurat}"
    printf "source_scvelo_metrics\t%s\n" "${source_scvelo}"
    printf "source_seurat_metadata_provenance\t%s\n" "${source_provenance}"
    printf "target_seurat_metadata\t%s\n" "${target_seurat}"
    printf "target_scvelo_metrics\t%s\n" "${target_scvelo}"
    printf "target_seurat_metadata_provenance\t%s\n" "${target_provenance}"
    printf "seurat_metadata_sha256\t%s\n" "$(shasum -a 256 "${target_seurat}" | awk '{print $1}')"
    printf "scvelo_metrics_sha256\t%s\n" "$(shasum -a 256 "${target_scvelo}" | awk '{print $1}')"
    printf "seurat_metadata_provenance_sha256\t%s\n" "$(shasum -a 256 "${target_provenance}" | awk '{print $1}')"
    printf "materialized_at\t%s\n" "$(date -Iseconds)"
  } > "${metadata_path}"
}

publish_si_figure4_input_pair() {
  local source_seurat="$1"
  local source_scvelo="$2"
  local source_provenance="$3"
  local data_root_real target_seurat target_scvelo target_provenance
  local temp_seurat temp_scvelo temp_provenance
  si_figure4_validate_input_pair "${source_seurat}" "${source_scvelo}" "${source_provenance}"
  data_root_real="$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${figure7_data_root}")"
  target_seurat="${data_root_real}/seurat_metadata.csv"
  target_scvelo="${data_root_real}/scvelo_cell_metrics.csv"
  target_provenance="${data_root_real}/seurat_metadata_provenance.tsv"
  mkdir -p "${data_root_real}"
  temp_seurat="${data_root_real}/.seurat_metadata.csv.tmp.${run_id}"
  temp_scvelo="${data_root_real}/.scvelo_cell_metrics.csv.tmp.${run_id}"
  temp_provenance="${data_root_real}/.seurat_metadata_provenance.tsv.tmp.${run_id}"
  cp "${source_seurat}" "${temp_seurat}"
  cp "${source_scvelo}" "${temp_scvelo}"
  cp "${source_provenance}" "${temp_provenance}"
  cmp -s "${source_seurat}" "${temp_seurat}"
  cmp -s "${source_scvelo}" "${temp_scvelo}"
  cmp -s "${source_provenance}" "${temp_provenance}"
  mv "${temp_seurat}" "${target_seurat}"
  mv "${temp_scvelo}" "${target_scvelo}"
  mv "${temp_provenance}" "${target_provenance}"
  cmp -s "${source_seurat}" "${target_seurat}"
  cmp -s "${source_scvelo}" "${target_scvelo}"
  cmp -s "${source_provenance}" "${target_provenance}"
  record_si_figure4_input_materialization \
    "ok" "${source_seurat}" "${source_scvelo}" "${source_provenance}" \
    "${target_seurat}" "${target_scvelo}" "${target_provenance}"
  si_figure4_seurat_metadata="${target_seurat}"
  si_figure4_scvelo_metrics="${target_scvelo}"
  si_figure4_seurat_metadata_provenance="${target_provenance}"
}

si_figure4_input_prep_command() {
  local intermediate_seurat="$1"
  local intermediate_scvelo="$2"
  local intermediate_provenance="$3"
  local -a args=(
    Rscript Code/in-vivo/figure7/run_figure7.R
    --mode=prepare-scvelo-inputs
    --config=Code/in-vivo/figure7/figure7_config.yaml
    "--intermediate-dir=${figure7_intermediate_dir}"
    "--seurat-metadata-output=${intermediate_seurat}"
    "--scvelo-metrics=${intermediate_scvelo}"
    "--seurat-metadata-provenance-output=${intermediate_provenance}"
    "--output-dir=${si_figure4_prep_run_dir}"
  )
  if [[ -n "${figure7_python}" ]]; then
    args+=("--python=${figure7_python}")
  fi
  if [[ "${figure7_overwrite_intermediates}" == true ]]; then
    args+=(--overwrite-intermediates=true)
  fi
  quote_args "${args[@]}"
}

resolve_si_figure4_inputs() {
  local data_root_real published_seurat published_scvelo published_provenance
  local intermediate_seurat intermediate_scvelo intermediate_provenance
  local published_count=0 intermediate_count=0 prep_command
  data_root_real="$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${figure7_data_root}")"
  published_seurat="${data_root_real}/seurat_metadata.csv"
  published_scvelo="${data_root_real}/scvelo_cell_metrics.csv"
  published_provenance="${data_root_real}/seurat_metadata_provenance.tsv"
  intermediate_seurat="${figure7_intermediate_dir}/seurat_metadata.csv"
  intermediate_scvelo="${figure7_intermediate_dir}/scvelo_cell_metrics.csv"
  intermediate_provenance="${figure7_intermediate_dir}/seurat_metadata_provenance.tsv"
  [[ -f "${published_seurat}" ]] && published_count=$((published_count + 1))
  [[ -f "${published_scvelo}" ]] && published_count=$((published_count + 1))
  [[ -f "${published_provenance}" ]] && published_count=$((published_count + 1))
  [[ -f "${intermediate_seurat}" ]] && intermediate_count=$((intermediate_count + 1))
  [[ -f "${intermediate_scvelo}" ]] && intermediate_count=$((intermediate_count + 1))
  [[ -f "${intermediate_provenance}" ]] && intermediate_count=$((intermediate_count + 1))

  if [[ "${published_count}" -eq 3 ]]; then
    si_figure4_validate_input_pair "${published_seurat}" "${published_scvelo}" "${published_provenance}"
    si_figure4_seurat_metadata="${published_seurat}"
    si_figure4_scvelo_metrics="${published_scvelo}"
    si_figure4_seurat_metadata_provenance="${published_provenance}"
    si_figure4_input_source="published_data"
    return 0
  fi
  if [[ "${intermediate_count}" -eq 3 ]]; then
    si_figure4_input_source="existing_figure7_intermediate"
    if [[ "${dry_run}" == true || "${mode}" == "check-only" ]]; then
      si_figure4_validate_input_pair "${intermediate_seurat}" "${intermediate_scvelo}" "${intermediate_provenance}"
      si_figure4_seurat_metadata="${intermediate_seurat}"
      si_figure4_scvelo_metrics="${intermediate_scvelo}"
      si_figure4_seurat_metadata_provenance="${intermediate_provenance}"
      printf "[si_figure4_input_publish] %s -> %s\n" "${figure7_intermediate_dir}" "${data_root_real}"
    else
      publish_si_figure4_input_pair "${intermediate_seurat}" "${intermediate_scvelo}" "${intermediate_provenance}"
    fi
    return 0
  fi

  si_figure4_input_source="planned_figure7_scvelo_generation"
  si_figure4_prep_run_dir="${output_root}/in-vivo/figure7/input-prep/runs/${run_id}_scvelo_inputs"
  si_figure4_seurat_metadata="${intermediate_seurat}"
  si_figure4_scvelo_metrics="${intermediate_scvelo}"
  si_figure4_seurat_metadata_provenance="${intermediate_provenance}"
  if [[ "${published_count}" -gt 0 || "${intermediate_count}" -gt 0 ]]; then
    figure7_overwrite_intermediates=true
  fi
  prep_command="$(si_figure4_input_prep_command \
    "${intermediate_seurat}" "${intermediate_scvelo}" "${intermediate_provenance}")"
  if [[ "${dry_run}" == true || "${mode}" == "check-only" ]]; then
    printf "[si_figure4_input_prep] %s\n" "${prep_command}"
    return 0
  fi
  if [[ -e "${si_figure4_prep_run_dir}" && "${overwrite}" != true ]]; then
    echo "SI Figure 4 input-prep run exists; use --overwrite: ${si_figure4_prep_run_dir}" >&2
    return 1
  fi
  if [[ -e "${si_figure4_prep_run_dir}" && "${overwrite}" == true ]]; then
    rm -rf "${si_figure4_prep_run_dir}"
  fi
  bash -c "${prep_command}"
  si_figure4_validate_input_pair "${intermediate_seurat}" "${intermediate_scvelo}" "${intermediate_provenance}"
  si_figure4_input_source="generated_figure7_scvelo_intermediate"
  publish_si_figure4_input_pair "${intermediate_seurat}" "${intermediate_scvelo}" "${intermediate_provenance}"
}

figure7_run_config_value() {
  local run_config="$1"
  local key="$2"
  awk -F '\t' -v target="${key}" '
    NR > 1 && $1 == target { count += 1; value = $2 }
    END {
      if (count != 1) exit 1
      print value
    }
  ' "${run_config}"
}

figure7_validate_publish_csv() {
  local artifact="$1"
  local path="$2"
  local header=""
  local required_column
  IFS= read -r header < "${path}" || true
  header="${header%$'\r'}"
  if [[ -z "${header}" ]]; then
    echo "Cannot publish empty Figure 7 CSV: ${path}" >&2
    return 1
  fi
  if [[ "${artifact}" == "seurat_metadata_provenance" ]]; then
    if [[ "${header}" != $'key\tvalue' ]]; then
      echo "Unexpected Seurat metadata provenance header; refusing publication: ${path}" >&2
      return 1
    fi
    return 0
  fi
  if [[ "${artifact}" == "scvelo" ]]; then
    if [[ "${header}" != cell,velocity_cell,* ]]; then
      echo "Unexpected scVelo metrics header; refusing publication: ${path}" >&2
      return 1
    fi
    return 0
  fi
  if [[ "${artifact}" == "seurat_metadata" ]]; then
    for required_column in cell UMAP_1 UMAP_2 Dose; do
      if [[ ",${header}," != *",${required_column},"* ]]; then
        echo "Seurat metadata CSV is missing ${required_column}; refusing publication: ${path}" >&2
        return 1
      fi
    done
    if [[ ",${header}," != *",cluster_final,"* && ",${header}," != *",clusters,"* ]]; then
      echo "Seurat metadata CSV is missing cluster_final/clusters; refusing publication: ${path}" >&2
      return 1
    fi
    return 0
  fi
  for required_column in \
    cell_id sample_id initial_ploidy gemcitabine_dose gemcitabine_dose_mg_per_kg \
    pseudotime cell_ploidy; do
    if [[ ",${header}," != *",${required_column},"* ]]; then
      echo "${artifact} CSV is missing ${required_column}; refusing publication: ${path}" >&2
      return 1
    fi
  done
}

record_figure7_input_materialization() {
  local status="$1"
  local started_at="$2"
  local finished_at="$3"
  local metadata_path="${manager_run_dir}/metadata/figure7_input_materialization.tsv"
  {
    printf "key\tvalue\n"
    printf "status\t%s\n" "${status}"
    printf "source_run_config\t%s\n" "${current_figure7_run_config:-}"
    printf "materialized_file_count\t%s\n" "${figure7_input_materialized_count}"
    printf "source_scvelo_metrics\t%s\n" "${figure7_input_source_scvelo}"
    printf "target_scvelo_metrics\t%s\n" "${figure7_input_target_scvelo}"
    printf "scvelo_sha256\t%s\n" "${figure7_input_sha_scvelo}"
    printf "source_seurat_metadata\t%s\n" "${figure7_input_source_seurat_metadata}"
    printf "target_seurat_metadata\t%s\n" "${figure7_input_target_seurat_metadata}"
    printf "seurat_metadata_sha256\t%s\n" "${figure7_input_sha_seurat_metadata}"
    printf "source_seurat_metadata_provenance\t%s\n" "${figure7_input_source_seurat_metadata_provenance}"
    printf "target_seurat_metadata_provenance\t%s\n" "${figure7_input_target_seurat_metadata_provenance}"
    printf "seurat_metadata_provenance_sha256\t%s\n" "${figure7_input_sha_seurat_metadata_provenance}"
    printf "source_cellcycle\t%s\n" "${figure7_input_source_cellcycle}"
    printf "target_cellcycle\t%s\n" "${figure7_input_target_cellcycle}"
    printf "cellcycle_sha256\t%s\n" "${figure7_input_sha_cellcycle}"
    printf "source_noncellcycle\t%s\n" "${figure7_input_source_noncellcycle}"
    printf "target_noncellcycle\t%s\n" "${figure7_input_target_noncellcycle}"
    printf "noncellcycle_sha256\t%s\n" "${figure7_input_sha_noncellcycle}"
    printf "started_at\t%s\n" "${started_at}"
    printf "finished_at\t%s\n" "${finished_at}"
  } > "${metadata_path}"
}

materialize_figure7_generated_inputs() {
  local run_dir="$1"
  local run_config="${run_dir}/metadata/run_config.tsv"
  local workflow_mode executed_stages data_root_real
  local artifact source expected_source target expected_hash observed_hash source_real expected_real
  local temporary_path
  local -a artifacts=()
  local -a sources=()
  local -a targets=()
  local -a expected_hashes=()
  local -a temporary_paths=()
  local i

  require_file "${run_config}" || return 1
  current_figure7_run_config="${run_config}"
  workflow_mode="$(figure7_run_config_value "${run_config}" mode)" || {
    echo "Figure 7 run_config.tsv does not contain one mode value" >&2
    return 1
  }
  if [[ "${workflow_mode}" != "full-workflow" ]]; then
    echo "Refusing input publication from non-full-workflow Figure 7 run: ${workflow_mode}" >&2
    return 1
  fi
  executed_stages="$(figure7_run_config_value "${run_config}" workflow_executed_stages)" || {
    echo "Figure 7 run_config.tsv does not contain workflow_executed_stages" >&2
    return 1
  }
  data_root_real="$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${figure7_data_root}")"

  if [[ ",${executed_stages}," == *",scvelo_metrics,"* ]]; then
    artifacts+=("scvelo" "seurat_metadata" "seurat_metadata_provenance")
    sources+=(
      "$(figure7_run_config_value "${run_config}" workflow_scvelo_metrics)"
      "$(figure7_run_config_value "${run_config}" workflow_seurat_metadata)"
      "$(figure7_run_config_value "${run_config}" workflow_seurat_metadata_provenance)"
    )
    targets+=(
      "${data_root_real}/scvelo_cell_metrics.csv"
      "${data_root_real}/seurat_metadata.csv"
      "${data_root_real}/seurat_metadata_provenance.tsv"
    )
    expected_hashes+=(
      "$(figure7_run_config_value "${run_config}" workflow_scvelo_sha256)"
      "$(figure7_run_config_value "${run_config}" workflow_seurat_metadata_sha256)"
      "$(figure7_run_config_value "${run_config}" workflow_seurat_metadata_provenance_sha256)"
    )
  fi
  if [[ ",${executed_stages}," == *",celllevel_inputs,"* ]]; then
    artifacts+=("cellcycle" "noncellcycle")
    sources+=(
      "$(figure7_run_config_value "${run_config}" workflow_cellcycle_input)"
      "$(figure7_run_config_value "${run_config}" workflow_noncellcycle_input)"
    )
    targets+=(
      "${data_root_real}/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
      "${data_root_real}/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    )
    expected_hashes+=(
      "$(figure7_run_config_value "${run_config}" workflow_cellcycle_sha256)"
      "$(figure7_run_config_value "${run_config}" workflow_noncellcycle_sha256)"
    )
  fi

  if [[ ",${executed_stages}," == *",scvelo_metrics,"* ]]; then
    si_figure4_validate_input_pair "${sources[1]}" "${sources[0]}" "${sources[2]}" || return 1
  fi

  figure7_input_materialized_count=0
  if [[ "${#artifacts[@]}" -eq 0 ]]; then
    return 0
  fi

  for i in "${!artifacts[@]}"; do
    artifact="${artifacts[$i]}"
    source="${sources[$i]}"
    target="${targets[$i]}"
    expected_hash="${expected_hashes[$i]}"
    expected_source="${figure7_intermediate_dir}/$(basename "${target}")"
    require_file "${source}" || return 1
    source_real="$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${source}")"
    expected_real="$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${expected_source}")"
    if [[ "${source_real}" != "${expected_real}" ]]; then
      echo "Figure 7 ${artifact} source is outside the configured intermediate directory: ${source_real}" >&2
      return 1
    fi
    if [[ ! "${expected_hash}" =~ ^[0-9a-f]{64}$ ]]; then
      echo "Invalid recorded SHA-256 for Figure 7 ${artifact}: ${expected_hash}" >&2
      return 1
    fi
    observed_hash="$(shasum -a 256 "${source}" | awk '{print $1}')"
    if [[ "${observed_hash}" != "${expected_hash}" ]]; then
      echo "Figure 7 ${artifact} changed after validation; refusing publication" >&2
      return 1
    fi
    figure7_validate_publish_csv "${artifact}" "${source}" || return 1
    mkdir -p "$(dirname "${target}")"
    temporary_path="$(dirname "${target}")/.$(basename "${target}").tmp.${run_id}"
    if ! cp "${source}" "${temporary_path}" || ! cmp -s "${source}" "${temporary_path}"; then
      echo "Could not stage validated Figure 7 ${artifact} for publication: ${target}" >&2
      return 1
    fi
    temporary_paths+=("${temporary_path}")
  done

  for i in "${!artifacts[@]}"; do
    artifact="${artifacts[$i]}"
    source="${sources[$i]}"
    target="${targets[$i]}"
    temporary_path="${temporary_paths[$i]}"
    if ! mv "${temporary_path}" "${target}" || ! cmp -s "${source}" "${target}"; then
      echo "Published Figure 7 ${artifact} does not match its validated source: ${target}" >&2
      return 1
    fi
    figure7_input_materialized_count=$((figure7_input_materialized_count + 1))
    case "${artifact}" in
      scvelo)
        figure7_input_source_scvelo="${source}"
        figure7_input_target_scvelo="${target}"
        figure7_input_sha_scvelo="${expected_hashes[$i]}"
        ;;
      seurat_metadata)
        figure7_input_source_seurat_metadata="${source}"
        figure7_input_target_seurat_metadata="${target}"
        figure7_input_sha_seurat_metadata="${expected_hashes[$i]}"
        ;;
      seurat_metadata_provenance)
        figure7_input_source_seurat_metadata_provenance="${source}"
        figure7_input_target_seurat_metadata_provenance="${target}"
        figure7_input_sha_seurat_metadata_provenance="${expected_hashes[$i]}"
        ;;
      cellcycle)
        figure7_input_source_cellcycle="${source}"
        figure7_input_target_cellcycle="${target}"
        figure7_input_sha_cellcycle="${expected_hashes[$i]}"
        ;;
      noncellcycle)
        figure7_input_source_noncellcycle="${source}"
        figure7_input_target_noncellcycle="${target}"
        figure7_input_sha_noncellcycle="${expected_hashes[$i]}"
        ;;
    esac
  done
  if [[ -n "${figure7_input_target_scvelo}" &&
        -n "${figure7_input_target_seurat_metadata}" &&
        -n "${figure7_input_target_seurat_metadata_provenance}" ]]; then
    si_figure4_validate_input_pair \
      "${figure7_input_target_seurat_metadata}" \
      "${figure7_input_target_scvelo}" \
      "${figure7_input_target_seurat_metadata_provenance}"
  fi
}

run_module() {
  local module="$1"
  local run_dir="$2"
  local command_string="$3"
  local module_notes=""
  if [[ "${module}" == "in_vivo_figure7" ]]; then
    module_notes="tgi_day=${figure7_tgi_day};figure_name=${figure7_figure_name}"
  elif [[ "${module}" == "si_figure4" ]]; then
    module_notes="input_source=${si_figure4_input_source}"
  fi

  check_module_inputs "${module}"

  if [[ "${dry_run}" == true || "${mode}" == "check-only" ]]; then
    if [[ "${module}" == "in_vivo_figure7" && -n "${figure7_state_pathway_results_root}" ]]; then
      printf "[in_vivo_figure7_export] %s\n" "$(figure7_reference_export_command)"
    fi
    printf "[%s] %s\n" "${module}" "${command_string}"
    return 0
  fi

  if [[ -e "${run_dir}" && "${overwrite}" != true ]]; then
    echo "Run directory exists; use --overwrite to replace: ${run_dir}" >&2
    return 1
  fi
  if [[ -e "${run_dir}" && "${overwrite}" == true ]]; then
    rm -rf "${run_dir}"
  fi
  mkdir -p "${run_dir}/metadata" "${run_dir}/logs"

  if [[ "${module}" == "in_vivo_figure7" && -n "${figure7_state_pathway_results_root}" ]]; then
    if ! prepare_figure7_reference; then
      local export_finished_at
      export_finished_at="$(date -Iseconds)"
      record_module_run \
        "${module}" "failed" "${command_string}" "${run_dir}" \
        "${manager_run_dir}/logs/figure7_state_pathway_export.stdout.log" \
        "${manager_run_dir}/logs/figure7_state_pathway_export.stderr.log" \
        "state_pathway_export_failed;state_pathway_source_results_root=${figure7_state_pathway_results_root}" \
        "${export_finished_at}" "${export_finished_at}"
      echo "Figure 7 state-pathway export failed; see ${manager_run_dir}/logs/figure7_state_pathway_export.stderr.log" >&2
      return 1
    fi
    check_module_inputs "${module}"
    module_notes="${module_notes};state_pathway_source_results_root=${figure7_state_pathway_results_root};state_pathway_reference_dir=${figure7_reference_root}"
  fi

  local stdout_log="${run_dir}/logs/stdout.log"
  local stderr_log="${run_dir}/logs/stderr.log"
  local status=0
  local started_at finished_at
  started_at="$(date -Iseconds)"
  bash -c "${command_string}" >"${stdout_log}" 2>"${stderr_log}" || status=$?
  finished_at="$(date -Iseconds)"
  if [[ "${status}" -ne 0 ]]; then
    local failure_notes="exit_status=${status}"
    if [[ -n "${module_notes}" ]]; then
      failure_notes="${failure_notes};${module_notes}"
    fi
    record_module_run "${module}" "failed" "${command_string}" "${run_dir}" "${stdout_log}" "${stderr_log}" "${failure_notes}" "${started_at}" "${finished_at}"
    echo "Module failed: ${module}; see ${stderr_log}" >&2
    return "${status}"
  fi

  if [[ "${module}" == "si_figure4" || "${module}" == "in_vivo_figure7" ]]; then
    python3 Code/tools/publish_figure7_upstream_analysis.py \
      --module "${module}" \
      --source-manifest Code/in-vivo/figure7/zenodo_upstream_analysis_documents.tsv \
      --output-dir "${run_dir}/metadata/upstream_analysis" \
      >>"${stdout_log}" 2>>"${stderr_log}" || status=$?
    finished_at="$(date -Iseconds)"
    if [[ "${status}" -ne 0 ]]; then
      record_module_run \
        "${module}" "failed" "${command_string}" "${run_dir}" "${stdout_log}" "${stderr_log}" \
        "upstream_analysis_documentation_publication_failed;exit_status=${status};${module_notes}" \
        "${started_at}" "${finished_at}"
      echo "Module failed while publishing upstream-analysis documentation: ${module}" >&2
      return "${status}"
    fi
  fi

  local input_args=()
  local input_path
  while IFS= read -r input_path; do
    [[ -z "${input_path}" ]] && continue
    if [[ -f "${input_path}" ]]; then
      input_args+=(--path "${input_path}")
    fi
  done < <(input_paths_for_module "${module}")
  if [[ "${#input_args[@]}" -gt 0 ]]; then
    python3 Code/tools/write_file_manifest.py \
      --manifest-type input \
      --module "${module}" \
      --output "${run_dir}/metadata/input_manifest.tsv" \
      --generated-by Manager.sh \
      --command-id "${run_id}" \
      "${input_args[@]}"
    python3 Code/tools/validate_manifest.py "${run_dir}/metadata/input_manifest.tsv" --repo-root "${repo_root}"
  fi
  python3 Code/tools/write_file_manifest.py \
    --manifest-type output \
    --module "${module}" \
    --output "${run_dir}/metadata/output_manifest.tsv" \
    --generated-by Manager.sh \
    --command-id "${run_id}" \
    --scan-dir "${run_dir}"
  python3 Code/tools/validate_manifest.py "${run_dir}/metadata/output_manifest.tsv" \
    --output-root "${run_dir}" --repo-root "${repo_root}"

  if [[ "${module}" == "in_vivo_figure7" && "${figure7_refresh_inputs}" == true ]]; then
    local input_materialization_started_at input_materialization_finished_at input_materialization_status
    input_materialization_started_at="$(date -Iseconds)"
    if ! materialize_figure7_generated_inputs "${run_dir}"; then
      input_materialization_finished_at="$(date -Iseconds)"
      record_figure7_input_materialization "failed" "${input_materialization_started_at}" "${input_materialization_finished_at}"
      finished_at="${input_materialization_finished_at}"
      record_module_run \
        "${module}" "failed" "${command_string}" "${run_dir}" "${stdout_log}" "${stderr_log}" \
        "input_materialization_failed;${module_notes};figure7_data_root=${figure7_data_root}" \
        "${started_at}" "${finished_at}"
      echo "Figure 7 validated input publication failed: ${figure7_data_root}" >&2
      return 1
    fi
    input_materialization_finished_at="$(date -Iseconds)"
    input_materialization_status="ok"
    if [[ "${figure7_input_materialized_count}" -eq 0 ]]; then
      input_materialization_status="skipped_no_generated_inputs"
    fi
    record_figure7_input_materialization "${input_materialization_status}" "${input_materialization_started_at}" "${input_materialization_finished_at}"
    module_notes="${module_notes};input_materialization_status=${input_materialization_status};input_materialized_file_count=${figure7_input_materialized_count};figure7_data_root=${figure7_data_root}"
    finished_at="${input_materialization_finished_at}"
  fi

  if [[ "${module}" == "in_vivo_figure7" && -n "${figure7_state_pathway_results_root}" ]]; then
    local materialization_started_at materialization_finished_at
    materialization_started_at="$(date -Iseconds)"
    if ! materialize_figure7_reference_to_data; then
      materialization_finished_at="$(date -Iseconds)"
      record_figure7_reference_materialization "failed" "${materialization_started_at}" "${materialization_finished_at}"
      finished_at="${materialization_finished_at}"
      record_module_run \
        "${module}" "failed" "${command_string}" "${run_dir}" "${stdout_log}" "${stderr_log}" \
        "state_pathway_materialization_failed;${module_notes};state_pathway_canonical_data_dir=${figure7_canonical_reference_root}" \
        "${started_at}" "${finished_at}"
      echo "Figure 7 canonical Data reference materialization failed: ${figure7_canonical_reference_root}" >&2
      return 1
    fi
    materialization_finished_at="$(date -Iseconds)"
    record_figure7_reference_materialization "ok" "${materialization_started_at}" "${materialization_finished_at}"
    module_notes="${module_notes};state_pathway_canonical_data_dir=${figure7_canonical_reference_root}"
    finished_at="${materialization_finished_at}"
  fi

  if [[ "${no_update_latest}" != true ]]; then
    write_latest "${module}" "${run_dir}" "${command_string}"
  fi
  record_module_run "${module}" "ok" "${command_string}" "${run_dir}" "${stdout_log}" "${stderr_log}" "${module_notes}" "${started_at}" "${finished_at}"
}

manager_run_dir="${output_root}/manager/runs/${run_id}"
module_runs_file="${manager_run_dir}/metadata/module_runs.tsv"
if [[ -n "${figure7_state_pathway_results_root}" ]]; then
  figure7_reference_root="${manager_run_dir}/artifacts/figure7_state_pathway_reference/${figure7_reference_id}"
fi

validate_module_registry "${module_list[@]}"

if [[ "${mode}" == "check-only" || "${dry_run}" == true ]]; then
  echo "Manager ${mode} for run_id=${run_id}"
else
  if [[ -e "${manager_run_dir}" && "${overwrite}" != true ]]; then
    echo "Manager run directory exists; use --overwrite to replace: ${manager_run_dir}" >&2
    exit 1
  fi
  if [[ -e "${manager_run_dir}" && "${overwrite}" == true ]]; then
    rm -rf "${manager_run_dir}"
  fi
  mkdir -p "${manager_run_dir}/metadata"
  printf "module\tstatus\tmode\tcommand\tresult_run_dir\tlatest_pointer\tstarted_at\tfinished_at\tstdout_log\tstderr_log\tnotes\n" \
    > "${module_runs_file}"
fi

completed_modules=()
completed_run_dirs=()

add_completed_run() {
  completed_modules+=("$1")
  completed_run_dirs+=("$2")
}

if [[ "${mode}" == "panels-only" ]]; then
  for module in "${module_list[@]}"; do
    case "${module}" in
      gdsc|ccle|drug_response|pkpd|in_vivo|si_figure4|in_vivo_figure7)
        run_dir="$(module_run_dir "${module}" "${source_run_id}")"
        [[ -d "${run_dir}" ]] || { echo "Missing panels-only source run: ${run_dir}" >&2; exit 1; }
        add_completed_run "${module}" "${run_dir}"
        ;;
      metabolomics)
        for submodule in metabolomics metabolomics_pathway metabolomics_zscore; do
          run_dir="$(module_run_dir "${submodule}" "${source_run_id}")"
          [[ -d "${run_dir}" ]] || { echo "Missing panels-only source run: ${run_dir}" >&2; exit 1; }
          add_completed_run "${submodule}" "${run_dir}"
        done
        ;;
      lci_overlays)
        run_dir="$(module_run_dir "${module}" "${source_run_id}")"
        [[ -d "${run_dir}" ]] || { echo "Missing panels-only source run: ${run_dir}" >&2; exit 1; }
        add_completed_run "${module}" "${run_dir}"
        ;;
      *) echo "Unknown module in --modules: ${module}" >&2; exit 2 ;;
    esac
  done
  dry_run=false
  no_update_latest=true
  skip_analysis_loop=true
  for i in "${!completed_modules[@]}"; do
    now="$(date -Iseconds)"
    record_module_run \
      "${completed_modules[$i]}" "source_selected" "" "${completed_run_dirs[$i]}" "" "" \
      "source_run_id=${source_run_id};materialization_operation_id=${run_id}" "${now}" "${now}"
  done
fi

if [[ "${skip_analysis_loop}" != true ]]; then
  for module in "${module_list[@]}"; do
    case "${module}" in
      gdsc|ccle|drug_response|pkpd|metabolomics|lci_overlays|in_vivo|si_figure4|in_vivo_figure7) ;;
      *) echo "Unknown module in --modules: ${module}" >&2; exit 2 ;;
    esac
    if [[ "${module}" == "lci_overlays" && -z "${lci_analysis_dir}" ]]; then
      echo "Skipping lci_overlays: --lci-analysis-dir was not provided."
      continue
    fi
    if [[ "${module}" == "in_vivo" && "${include_in_vivo}" != true ]]; then
      echo "Skipping in_vivo: use --include-in-vivo to include pending in-vivo outputs."
      continue
    fi

    if [[ "${module}" == "si_figure4" ]]; then
      resolve_si_figure4_inputs
    fi

    if [[ "${module}" == "pkpd" && ( "${mode}" == "full-refit" || "${pkpd_refit}" == true ) ]]; then
      fit_run_dir="${pkpd_fit_output:-$(module_run_dir pkpd_fit)}"
      fit_command="$(command_for_module pkpd_fit "${fit_run_dir}")"
      run_module pkpd_fit "${fit_run_dir}" "${fit_command}"
      add_completed_run pkpd "${fit_run_dir}"
      plot_command="$(quote_args python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py --fit-output "${fit_run_dir}" --output-dir "${fit_run_dir}")"
      if [[ "${dry_run}" == true || "${mode}" == "check-only" ]]; then
        printf "[pkpd_plot] %s\n" "${plot_command}"
      else
        bash -c "${plot_command}" >>"${fit_run_dir}/logs/stdout.log" 2>>"${fit_run_dir}/logs/stderr.log"
        python3 Code/tools/write_file_manifest.py \
          --manifest-type output \
          --module pkpd \
          --output "${fit_run_dir}/metadata/output_manifest.tsv" \
          --generated-by Manager.sh \
          --command-id "${run_id}" \
          --scan-dir "${fit_run_dir}"
        python3 Code/tools/validate_manifest.py "${fit_run_dir}/metadata/output_manifest.tsv" \
          --output-root "${fit_run_dir}" --repo-root "${repo_root}"
      fi
      continue
    fi

    if [[ "${module}" == "metabolomics" ]]; then
      full_dir="$(module_run_dir metabolomics)"
      pathway_dir="$(module_run_dir metabolomics_pathway)"
      zscore_dir="$(module_run_dir metabolomics_zscore)"
      run_module metabolomics "${full_dir}" "$(command_for_module metabolomics "${full_dir}")"
      run_module metabolomics_pathway "${pathway_dir}" "$(command_for_module metabolomics_pathway "${pathway_dir}")"
      run_module metabolomics_zscore "${zscore_dir}" "$(command_for_module metabolomics_zscore "${zscore_dir}")"
      add_completed_run metabolomics "${full_dir}"
      add_completed_run metabolomics_pathway "${pathway_dir}"
      add_completed_run metabolomics_zscore "${zscore_dir}"
      continue
    fi

    run_dir="$(module_run_dir "${module}")"
    run_module "${module}" "${run_dir}" "$(command_for_module "${module}" "${run_dir}")"
    add_completed_run "${module}" "${run_dir}"
  done
fi

if [[ "${mode}" != "check-only" && "${dry_run}" != true ]]; then
  materialize_source_run_id="${run_id}"
  if [[ "${mode}" == "panels-only" ]]; then
    materialize_source_run_id="${source_run_id}"
  fi
  materialize_args=(
    python3 Code/tools/materialize_figure_assets.py
    --figure-root "${figure_root}"
    --figure7-tgi-day "${figure7_tgi_day}"
    --figure7-figure-name "${figure7_figure_name}"
    --source-run-id "${materialize_source_run_id}"
    --operation-id "${run_id}"
    --overwrite
  )
  for i in "${!completed_modules[@]}"; do
    materialize_args+=(--module-run "${completed_modules[$i]}=${completed_run_dirs[$i]}")
  done
  "${materialize_args[@]}"
  for manifest in "${figure_root}"/Figure*/manifest.tsv "${figure_root}"/Supplementary/manifest.tsv; do
    [[ -f "${manifest}" ]] || continue
    python3 Code/tools/validate_figure_manifest.py "${manifest}" --repo-root "${repo_root}"
  done
fi

echo "Manager completed: ${run_id}"
