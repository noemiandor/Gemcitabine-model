#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

python3 scripts/generate_timecourse_overlays.py --well D6_2 --crop-fraction 0.2 --columns 6
python3 scripts/generate_timecourse_overlays.py --well E6_2 --crop-fraction 0.2 --columns 6
