#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 INPUT_DIR OUTPUT_DIR" >&2
  exit 64
fi

input_dir="$1"
output_dir="$2"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
phase="${CLUSTER_STANDALONE_PHASE:-orchestrate}"

if [[ "${phase}" == "orchestrate" ]]; then
  pre_filter_sif="${CLUSTER_STANDALONE_PRE_FILTER_SIF:-}"
  post_filter_sif="${CLUSTER_STANDALONE_POST_FILTER_SIF:-}"
  if [[ -z "${pre_filter_sif}" || -z "${post_filter_sif}" ]]; then
    echo "CLUSTER_STANDALONE_PRE_FILTER_SIF and CLUSTER_STANDALONE_POST_FILTER_SIF are required." >&2
    exit 64
  fi
  if [[ ! -f "${pre_filter_sif}" ]]; then
    echo "Pre-filter SIF does not exist: ${pre_filter_sif}" >&2
    exit 66
  fi
  if [[ ! -f "${post_filter_sif}" ]]; then
    echo "Post-filter SIF does not exist: ${post_filter_sif}" >&2
    exit 66
  fi
  if ! command -v apptainer >/dev/null 2>&1; then
    echo "apptainer is required to orchestrate the two runtime phases." >&2
    exit 69
  fi

  mkdir -p "${output_dir}"
  apptainer exec --cleanenv \
    --env CLUSTER_STANDALONE_PHASE=pre_filter \
    --env CLUSTER_STANDALONE_ACTIVE_SIF="${pre_filter_sif}" \
    "${pre_filter_sif}" \
    bash "${script_dir}/run_cluster_standalone.sh" "${input_dir}" "${output_dir}"

  apptainer exec --cleanenv \
    --env CLUSTER_STANDALONE_PHASE=post_filter \
    --env CLUSTER_STANDALONE_ACTIVE_SIF="${post_filter_sif}" \
    "${post_filter_sif}" \
    bash "${script_dir}/run_cluster_standalone.sh" "${input_dir}" "${output_dir}"
  exit 0
fi

if [[ "${phase}" != "pre_filter" && "${phase}" != "post_filter" ]]; then
  echo "Invalid CLUSTER_STANDALONE_PHASE: ${phase}" >&2
  exit 64
fi

mkdir -p "${output_dir}"
if [[ "${phase}" == "post_filter" ]]; then
  post_filter_empty_library="${output_dir}/post_filter_empty_r_library"
  mkdir -p "${post_filter_empty_library}"
  export R_LIBS_USER="${post_filter_empty_library}"
  unset CLUSTER_STANDALONE_R_LIBRARY R_MAKEVARS_USER MAKEFLAGS || true
  exec Rscript --vanilla "${script_dir}/cluster_pipeline_standalone.R" \
    "${input_dir}" \
    "${output_dir}" \
    "${phase}"
fi

runtime_id="$(
  Rscript --vanilla -e \
    'cat(sprintf("R-%s-%s", getRversion(), R.version$platform))'
)"
if [[ ! "${runtime_id}" =~ ^R-[A-Za-z0-9._-]+$ ]]; then
  echo "Cannot derive a safe R runtime identifier: ${runtime_id}" >&2
  exit 70
fi
private_r_library="${output_dir}/cluster_standalone_r_library/${runtime_id}"
mkdir -p "${private_r_library}"
# On the target shared filesystem, R's access(2)-based writability check runs
# through the container user mapping and requires the "other" write bit. The
# parent result directory remains access-restricted.
chmod 0777 "${private_r_library}"

export R_LIBS_USER="${private_r_library}"
export CLUSTER_STANDALONE_R_LIBRARY="${private_r_library}"
available_cpus="$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')"
if [[ "${available_cpus}" =~ ^[0-9]+$ ]] && (( available_cpus > 2 )); then
  build_jobs="$((available_cpus - 2))"
else
  build_jobs="1"
fi
export MAKEFLAGS="-j${build_jobs}"
unset R_MAKEVARS_USER || true

Rscript --vanilla "${script_dir}/bootstrap_dependencies.R" \
  "${private_r_library}" \
  "${script_dir}/vendor"

exec Rscript --vanilla "${script_dir}/cluster_pipeline_standalone.R" \
  "${input_dir}" \
  "${output_dir}" \
  "${phase}"
