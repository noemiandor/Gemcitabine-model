#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${repo_root}"

run_id="$(date +"%Y%m%dT%H%M%S_manuscript")"
source_run_id=""
mode="standard"
modules="gdsc,ccle,drug_response,pkpd,metabolomics,in_vivo_figure7"
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
figure7_panels_ae_only=false
figure7_seurat_rds=""
figure7_gene_set_artifact=""
figure7_reference_id="taoli_04i_etp2_24_day17_v1"
figure7_state_pathway_results_root=""
figure7_canonical_reference_root="${FIGURE7_CANONICAL_REFERENCE_ROOT:-Data/in-vivo/figure7/saved_state_pathway/${figure7_reference_id}}"
figure7_reference_root="${figure7_canonical_reference_root}"
skip_analysis_loop=false

usage() {
  cat <<'EOF'
Usage: bash Manager.sh [options]

Core options:
  --run-id ID
  --source-run-id ID              Required only for panels-only; immutable source manager run
  --mode check-only|saved-fit|standard|full-refit|panels-only
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
  --figure7-full-analysis          Opt-in full pathway recomputation; requires both paths below
  --figure7-panels-ae-only         Generate/materialize 7A-7E while canonical panel 7F is unavailable
  --figure7-seurat-rds ABSOLUTE_PATH
  --figure7-gene-set-artifact PATH  Pinned, versioned local gene-set artifact
  --figure7-state-pathway-results-root PATH
                                  Export panel-7F reference from this completed 04i result tree
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
    --figure7-panels-ae-only) figure7_panels_ae_only=true; shift ;;
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
if [[ -n "${figure7_state_pathway_results_root}" ]]; then
  if [[ "${mode}" == "panels-only" ]]; then
    echo "--figure7-state-pathway-results-root cannot be used with --mode panels-only" >&2
    exit 2
  fi
  if [[ "${figure7_panels_ae_only}" == true ]]; then
    echo "--figure7-state-pathway-results-root cannot be used with --figure7-panels-ae-only" >&2
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
    in_vivo_figure7)
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
      ;;
  esac
}

check_module_inputs() {
  local module="$1"
  local path
  if [[ "${module}" == "in_vivo_figure7" && -n "${figure7_state_pathway_results_root}" ]]; then
    require_dir "${figure7_state_pathway_results_root}"
    require_file Code/in-vivo/figure7/export_04i_state_pathway_reference.R
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
    in_vivo_figure7)
      local figure7_mode="standard"
      local figure7_args=(
        Rscript Code/in-vivo/figure7/run_figure7.R
        "--mode=${figure7_mode}"
        --config=Code/in-vivo/figure7/figure7_config.yaml
        --cellcycle-input=Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
        --non-cellcycle-input=Data/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
        "--output-dir=${run_dir}"
      )
      if [[ "${figure7_panels_ae_only}" == true ]]; then
        figure7_args+=(--panel-set=a-e)
      else
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
  quote_args Rscript Code/in-vivo/figure7/export_04i_state_pathway_reference.R \
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
    printf "source_report_html\t%s\n" "${figure7_state_pathway_results_root}/report/04i_pseudotime_state_pathways_report.html"
    printf "canonical_reference_id\t%s\n" "${figure7_reference_id}"
    printf "exported_reference_dir\t%s\n" "${figure7_reference_root}"
    printf "canonical_data_reference_dir\t%s\n" "${figure7_canonical_reference_root}"
    printf "canonical_data_materialization_metadata\t%s\n" "${manager_run_dir}/metadata/figure7_state_pathway_materialization.tsv"
    printf "exporter_script\t%s\n" "Code/in-vivo/figure7/export_04i_state_pathway_reference.R"
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

run_module() {
  local module="$1"
  local run_dir="$2"
  local command_string="$3"
  local module_notes=""

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
    module_notes="state_pathway_source_results_root=${figure7_state_pathway_results_root};state_pathway_reference_dir=${figure7_reference_root}"
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
      gdsc|ccle|drug_response|pkpd|in_vivo|in_vivo_figure7)
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
      gdsc|ccle|drug_response|pkpd|metabolomics|lci_overlays|in_vivo|in_vivo_figure7) ;;
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
