#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

crosswalk=""
paired_sensitivity=""
workspace=""
expected_input_manifest=""
fcs_dir=""
output_dir=""
expected_samples="16"
min_human_cells="1000"

usage() {
  cat <<'EOF'
Usage:
  bash Code/in-vivo/flow_cytometry/run_endpoint_flow.sh \
    --crosswalk FILE --paired-sensitivity FILE --workspace FILE \
    --expected-input-manifest FILE --fcs-dir DIR --output-dir DIR [options]

Runs the frozen FlowJo-count extractor first, then reconstructs the exact
workspace gates over the raw FCS events and renders Supplementary Figure 9.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --crosswalk) crosswalk="$2"; shift 2 ;;
    --paired-sensitivity) paired_sensitivity="$2"; shift 2 ;;
    --workspace) workspace="$2"; shift 2 ;;
    --expected-input-manifest) expected_input_manifest="$2"; shift 2 ;;
    --fcs-dir) fcs_dir="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --expected-samples) expected_samples="$2"; shift 2 ;;
    --min-human-cells) min_human_cells="$2"; shift 2 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 64 ;;
  esac
done

for required_value in crosswalk paired_sensitivity workspace expected_input_manifest fcs_dir output_dir; do
  if [[ -z "${!required_value}" ]]; then
    printf 'Missing required option: --%s\n' "${required_value//_/-}" >&2
    exit 64
  fi
done

cd "${repo_root}"

python3 Code/in-vivo/flow_cytometry/extract_endpoint_flow.py \
  --crosswalk "${crosswalk}" \
  --paired-sensitivity "${paired_sensitivity}" \
  --workspace "${workspace}" \
  --expected-input-manifest "${expected_input_manifest}" \
  --fcs-dir "${fcs_dir}" \
  --expected-samples "${expected_samples}" \
  --min-human-cells "${min_human_cells}" \
  --output-dir "${output_dir}"

scripts/agentRrunner.sh Code/in-vivo/flow_cytometry/reconstruct_endpoint_flow.R \
  --crosswalk "${crosswalk}" \
  --workspace "${workspace}" \
  --fcs-dir "${fcs_dir}" \
  --frozen-table "${output_dir}/endpoint_flow_per_sample.tsv" \
  --expected-samples "${expected_samples}" \
  --min-human-cells "${min_human_cells}" \
  --output-dir "${output_dir}"
