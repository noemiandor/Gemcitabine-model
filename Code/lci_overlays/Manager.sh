#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "${script_dir}/generate_timecourse_overlays.py" --well D6_2 --crop-fraction 0.2 --columns 6 "$@"
python3 "${script_dir}/generate_timecourse_overlays.py" --well E6_2 --crop-fraction 0.2 --columns 6 "$@"
