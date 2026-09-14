#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_PATH="${1:-${SCRIPT_DIR}/05_mouse_cell_analysis_config.yaml}"
RESULTS_ROOT="${MOUSE05_RESULTS_ROOT:-/share/lab_crd/taoli/Project/BreastCancerOrthotopicModels/Results/05_mouse_cell_analysis}"
SIF_PATH="${MOUSE05_SIF:-/share/lab_crd/taoli/Docker/gemcitabine-model_full.sif}"
LOG_DIR="${RESULTS_ROOT}/logs"
STATUS_FILE="${LOG_DIR}/pipeline_status.tsv"

mkdir -p "${LOG_DIR}"

# The runner always re-executes itself inside the approved full SIF. This keeps
# every R stage in one immutable environment even when the outer nohup shell is
# started from a host R module session.
if [[ "${MOUSE05_INSIDE_SIF:-0}" != "1" ]]; then
  if [[ ! -f "${SIF_PATH}" ]]; then
    printf 'SIF does not exist: %s\n' "${SIF_PATH}" >&2
    exit 2
  fi
  if ! command -v apptainer >/dev/null 2>&1; then
    printf 'apptainer is not available on host PATH.\n' >&2
    exit 2
  fi
  exec apptainer exec \
    --cleanenv \
    --bind /share/lab_crd:/share/lab_crd \
    --pwd "${PROJECT_ROOT}" \
    --env MOUSE05_INSIDE_SIF=1 \
    --env MOUSE05_RESULTS_ROOT="${RESULTS_ROOT}" \
    --env MOUSE05_SIF="${SIF_PATH}" \
    "${SIF_PATH}" \
    bash "${SCRIPT_DIR}/05_run_mouse_pipeline.sh" "${CONFIG_PATH}"
fi

if [[ ! -f "${STATUS_FILE}" ]]; then
  printf 'stage\tstatus\tstarted_at\tfinished_at\texit_code\tlog\n' > "${STATUS_FILE}"
fi

run_stage() {
  local stage="$1"
  local script="$2"
  local stage_log="${LOG_DIR}/${stage}.log"
  local done_file="${LOG_DIR}/${stage}.done"
  local started_at
  local finished_at

  if [[ -f "${done_file}" ]]; then
    printf '%s\tSKIPPED_ALREADY_DONE\t%s\t%s\t0\t%s\n' "${stage}" "$(date --iso-8601=seconds)" "$(date --iso-8601=seconds)" "${stage_log}" >> "${STATUS_FILE}"
    return 0
  fi

  started_at="$(date --iso-8601=seconds)"
  printf '%s\tRUNNING\t%s\tNA\tNA\t%s\n' "${stage}" "${started_at}" "${stage_log}" >> "${STATUS_FILE}"
  if Rscript "${SCRIPT_DIR}/${script}" "${CONFIG_PATH}" >> "${stage_log}" 2>&1; then
    finished_at="$(date --iso-8601=seconds)"
    touch "${done_file}"
    printf '%s\tCOMPLETED\t%s\t%s\t0\t%s\n' "${stage}" "${started_at}" "${finished_at}" "${stage_log}" >> "${STATUS_FILE}"
  else
    local rc=$?
    finished_at="$(date --iso-8601=seconds)"
    printf '%s\tFAILED\t%s\t%s\t%s\t%s\n' "${stage}" "${started_at}" "${finished_at}" "${rc}" "${stage_log}" >> "${STATUS_FILE}"
    return "${rc}"
  fi
}

printf 'pipeline_pid\t%s\n' "$$" > "${LOG_DIR}/pipeline.pid"
printf 'pipeline_started_at\t%s\n' "$(date --iso-8601=seconds)" > "${LOG_DIR}/pipeline_runtime.tsv"
printf 'host\t%s\n' "$(hostname)" >> "${LOG_DIR}/pipeline_runtime.tsv"
printf 'r_version\t%s\n' "$(R --version | head -n 1)" >> "${LOG_DIR}/pipeline_runtime.tsv"
printf 'sif_path\t%s\n' "${SIF_PATH}" >> "${LOG_DIR}/pipeline_runtime.tsv"
printf 'sif_sha256\t%s\n' "$(sha256sum "${SIF_PATH}" | awk '{print $1}')" >> "${LOG_DIR}/pipeline_runtime.tsv"
printf 'inside_sif\t%s\n' "${MOUSE05_INSIDE_SIF}" >> "${LOG_DIR}/pipeline_runtime.tsv"

run_stage "05a" "05a_extract_mouse_cells.R"
run_stage "05b" "05b_mouse_cell_qc.R"
run_stage "05c" "05c_mouse_sct_cca_integration.R"
run_stage "05d" "05d_mouse_clustering.R"
run_stage "05e" "05e_mouse_broad_annotation.R"
run_stage "05f" "05f_mouse_immune_annotation.R"
run_stage "05g" "05g_mouse_annotation_finalize.R"

printf 'pipeline_finished_at\t%s\n' "$(date --iso-8601=seconds)" >> "${LOG_DIR}/pipeline_runtime.tsv"
printf 'pipeline_status\tCOMPLETED\n' >> "${LOG_DIR}/pipeline_runtime.tsv"
