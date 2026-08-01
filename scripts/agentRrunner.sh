#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
rscript_bin="${RSCRIPT_BIN:-Rscript}"

usage() {
  cat <<'EOF'
Usage:
  scripts/agentRrunner.sh --check
  scripts/agentRrunner.sh <script.R> [arguments ...]
  scripts/agentRrunner.sh -e <expression>

The command runs Rscript --vanilla from the repository root. It preserves the
caller's PATH and R library environment and exports GEMCITABINE_MODEL_REPO_ROOT.
Set RSCRIPT_BIN to select a non-default Rscript executable.
EOF
}

if [[ $# -eq 0 ]]; then
  usage >&2
  exit 64
fi

if [[ "${1}" == "-h" || "${1}" == "--help" ]]; then
  usage
  exit 0
fi

if ! command -v "${rscript_bin}" >/dev/null 2>&1; then
  printf 'error: Rscript executable not found: %s\n' "${rscript_bin}" >&2
  exit 127
fi

export GEMCITABINE_MODEL_REPO_ROOT="${repo_root}"

if [[ "${1}" == "--check" ]]; then
  printf 'Repository root: %s\n' "${repo_root}"
  "${rscript_bin}" --version
  "${rscript_bin}" --vanilla -e \
    'cat("R runner smoke test: OK\n"); cat("Library paths:\n"); cat(paste0("- ", .libPaths(), collapse = "\n"), "\n", sep = "")'
  exit 0
fi

cd "${repo_root}"
exec "${rscript_bin}" --vanilla "$@"
