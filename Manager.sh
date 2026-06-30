#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${repo_root}"

run_id="$(date +"%Y%m%dT%H%M%S_manuscript")"
mode="standard"
modules="gdsc,ccle,drug_response,pkpd,metabolomics"
output_root="Results"
figure_root="figures"
module_registry="docs/manuscript_figure_module_registry.tsv"
overwrite=false
dry_run=false
no_update_latest=false
jobs=1

gdsc_category_mode="curated"
gdsc_analysis_mode="manuscript"
gdsc_enrichment_permute_n="300"

pkpd_saved_fit="Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/alsoGoodFit_20260514T093906"
pkpd_refit=false
pkpd_smoke_fit=false
pkpd_fit_output=""

refresh_cloneid_ploidy=false
lci_analysis_dir=""
lci_render=false
lci_panel_only=false
include_in_vivo=false
metabolomics_input="Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_2N_4N_Full.xlsm"
skip_analysis_loop=false

usage() {
  cat <<'EOF'
Usage: bash Manager.sh [options]

Core options:
  --run-id ID
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
  --gdsc-category-mode curated|legacy|proposal
  --gdsc-analysis-mode dev|manuscript
  --gdsc-enrichment-permute-n N
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
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-id) run_id="$2"; shift 2 ;;
    --mode) mode="$2"; shift 2 ;;
    --modules) modules="$2"; shift 2 ;;
    --output-root) output_root="$2"; shift 2 ;;
    --figure-root) figure_root="$2"; shift 2 ;;
    --overwrite) overwrite=true; shift ;;
    --jobs) jobs="$2"; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    --check-inputs) mode="check-only"; shift ;;
    --no-update-latest) no_update_latest=true; shift ;;
    --gdsc-category-mode) gdsc_category_mode="$2"; shift 2 ;;
    --gdsc-analysis-mode) gdsc_analysis_mode="$2"; shift 2 ;;
    --gdsc-enrichment-permute-n) gdsc_enrichment_permute_n="$2"; shift 2 ;;
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
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "${mode}" in
  check-only|saved-fit|standard|full-refit|panels-only) ;;
  *) echo "Invalid --mode: ${mode}" >&2; exit 2 ;;
esac

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
  case "${module}" in
    gdsc) printf "%s/public_data/gdsc_ploidy_analysis/runs/%s_gdsc" "${output_root}" "${run_id}" ;;
    ccle) printf "%s/public_data/ccle_ploidy_analysis/runs/%s_ccle" "${output_root}" "${run_id}" ;;
    drug_response) printf "%s/in-vitro/drug_response/runs/%s_drug_response" "${output_root}" "${run_id}" ;;
    lci_overlays) printf "%s/in-vitro/lci_overlays/runs/%s_lci_overlays" "${output_root}" "${run_id}" ;;
    pkpd) printf "%s/in-vitro/pkpd_live_dead_model/runs/%s_pkpd_saved_fit" "${output_root}" "${run_id}" ;;
    pkpd_fit) printf "%s/in-vitro/pkpd_live_dead_model/runs/%s_pkpd_fit" "${output_root}" "${run_id}" ;;
    metabolomics) printf "%s/in-vitro/metabolomics/runs/%s_metabolomics" "${output_root}" "${run_id}" ;;
    metabolomics_pathway) printf "%s/in-vitro/metabolomics/runs/%s_metabolomics_pathway" "${output_root}" "${run_id}" ;;
    metabolomics_zscore) printf "%s/in-vitro/metabolomics/runs/%s_metabolomics_zscore" "${output_root}" "${run_id}" ;;
    in_vivo) printf "%s/in-vivo/pseudotime_associations/runs/%s_in_vivo" "${output_root}" "${run_id}" ;;
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
        Code/gdsc_ploidy_analysis/data/raw/small_molecule_20200407234909.csv \
        Code/gdsc_ploidy_analysis/data/manual/drug_class_final_curated.tsv \
        Code/gdsc_ploidy_analysis/data/derived/pubchem_drug_annotations.tsv
      ;;
    ccle)
      printf "%s\n" \
        Code/ccle_ploidy_analysis/data/raw/Cell_app_export.txt \
        Code/ccle_ploidy_analysis/data/manual/breast_ccle_ploidy.tsv \
        Code/ccle_ploidy_analysis/data/raw/DrugAliases.txt
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
  esac
}

check_module_inputs() {
  local module="$1"
  local path
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
        "--category-mode=${gdsc_category_mode}" \
        "--enrichment-permute-n=${gdsc_enrichment_permute_n}" \
        "--output-dir=${run_dir}"
      ;;
    ccle)
      quote_args Rscript Code/ccle_ploidy_analysis/run_ccle_ploidy_analysis.R \
        "--metric=IC50" \
        "--metric-source=grbrowser" \
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
      local fit_args=(python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --output-dir "${run_dir}" --n-jobs "${jobs}")
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

run_module() {
  local module="$1"
  local run_dir="$2"
  local command_string="$3"

  check_module_inputs "${module}"

  if [[ "${dry_run}" == true || "${mode}" == "check-only" ]]; then
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

  local stdout_log="${run_dir}/logs/stdout.log"
  local stderr_log="${run_dir}/logs/stderr.log"
  local status=0
  local started_at finished_at
  started_at="$(date -Iseconds)"
  bash -c "${command_string}" >"${stdout_log}" 2>"${stderr_log}" || status=$?
  finished_at="$(date -Iseconds)"
  if [[ "${status}" -ne 0 ]]; then
    record_module_run "${module}" "failed" "${command_string}" "${run_dir}" "${stdout_log}" "${stderr_log}" "exit_status=${status}" "${started_at}" "${finished_at}"
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
    python3 Code/tools/validate_manifest.py "${run_dir}/metadata/input_manifest.tsv"
  fi
  python3 Code/tools/write_file_manifest.py \
    --manifest-type output \
    --module "${module}" \
    --output "${run_dir}/metadata/output_manifest.tsv" \
    --generated-by Manager.sh \
    --command-id "${run_id}" \
    --scan-dir "${run_dir}"
  python3 Code/tools/validate_manifest.py "${run_dir}/metadata/output_manifest.tsv" --output-root "${run_dir}"

  if [[ "${no_update_latest}" != true ]]; then
    write_latest "${module}" "${run_dir}" "${command_string}"
  fi
  record_module_run "${module}" "ok" "${command_string}" "${run_dir}" "${stdout_log}" "${stderr_log}" "" "${started_at}" "${finished_at}"
}

manager_run_dir="${output_root}/manager/runs/${run_id}"
module_runs_file="${manager_run_dir}/metadata/module_runs.tsv"

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
      gdsc|ccle|drug_response|pkpd)
        run_dir="$(module_run_dir "${module}")"
        [[ -d "${run_dir}" ]] || { echo "Missing panels-only source run: ${run_dir}" >&2; exit 1; }
        add_completed_run "${module}" "${run_dir}"
        ;;
      metabolomics)
        for submodule in metabolomics metabolomics_pathway metabolomics_zscore; do
          run_dir="$(module_run_dir "${submodule}")"
          [[ -d "${run_dir}" ]] || { echo "Missing panels-only source run: ${run_dir}" >&2; exit 1; }
          add_completed_run "${submodule}" "${run_dir}"
        done
        ;;
      lci_overlays)
        run_dir="$(module_run_dir "${module}")"
        [[ -d "${run_dir}" ]] || { echo "Missing panels-only source run: ${run_dir}" >&2; exit 1; }
        add_completed_run "${module}" "${run_dir}"
        ;;
      in_vivo)
        run_dir="$(module_run_dir "${module}")"
        [[ -d "${run_dir}" ]] || { echo "Missing panels-only source run: ${run_dir}" >&2; exit 1; }
        add_completed_run "${module}" "${run_dir}"
        ;;
      *) echo "Unknown module in --modules: ${module}" >&2; exit 2 ;;
    esac
  done
  dry_run=false
  no_update_latest=true
  skip_analysis_loop=true
fi

if [[ "${skip_analysis_loop}" != true ]]; then
  for module in "${module_list[@]}"; do
    case "${module}" in
      gdsc|ccle|drug_response|pkpd|metabolomics|lci_overlays|in_vivo) ;;
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
        python3 Code/tools/validate_manifest.py "${fit_run_dir}/metadata/output_manifest.tsv" --output-root "${fit_run_dir}"
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
  materialize_args=(python3 Code/tools/materialize_figure_assets.py --figure-root "${figure_root}" --run-id "${run_id}" --overwrite)
  for i in "${!completed_modules[@]}"; do
    materialize_args+=(--module-run "${completed_modules[$i]}=${completed_run_dirs[$i]}")
  done
  "${materialize_args[@]}"
  for manifest in "${figure_root}"/Figure*/manifest.tsv "${figure_root}"/Supplementary/manifest.tsv; do
    [[ -f "${manifest}" ]] || continue
    python3 Code/tools/validate_figure_manifest.py "${manifest}"
  done
fi

echo "Manager completed: ${run_id}"
