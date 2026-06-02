#!/usr/bin/env bash
set -euo pipefail

# Submit prepared trajectory scVelo tasks as one Slurm array.

PROJECT="${PROJECT:-/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels}"
CODE_DIR="${CODE_DIR:-$PROJECT/Code/in-vivo}"
RESULTS="${RESULTS:-$PROJECT/Results}"
ENV_NAME="${ENV_NAME:-rna_velocity_py310}"
TRAJECTORY_SCRIPT="${TRAJECTORY_SCRIPT:-04_trajectory.R}"
TRAJECTORY_RESULT_NAME="${TRAJECTORY_RESULT_NAME:-04_trajectory}"
ENV_DIR="${ENV_DIR:-}"
CONDA_MODULE="${CONDA_MODULE:-Miniconda3/23.9.0-0}"
R_MODULE="${R_MODULE:-R/4.4}"
PARTITION="${PARTITION:-}"
QOS="${QOS:-xlarge}"
TASK_TIME="${TASK_TIME:-04:00:00}"

SCVELO_CPUS="${SCVELO_CPUS:-4}"
SCVELO_N_JOBS="${SCVELO_N_JOBS:-$SCVELO_CPUS}"
SKIP_PREP="${SKIP_PREP:-FALSE}"
SCVELO_MEM="${SCVELO_MEM:-256G}"

TRAJECTORY_ROOT="$RESULTS/$TRAJECTORY_RESULT_NAME"
SUMMARY_DIR="$TRAJECTORY_ROOT/02_summary"
SUMMARY_FILE="$SUMMARY_DIR/scvelo_run_summary.csv"
TASK_LIST="$SUMMARY_DIR/scvelo_tasks.txt"
LOOM_ROOT="$TRAJECTORY_ROOT/00_inputs/velocyto_loom"
SBATCH_WORKER="$SUMMARY_DIR/submit_scvelo_array_worker.sbatch"

find_conda_env_dir() {
  local env_name="$1"
  conda env list | awk -v env_name="$env_name" '
    $1 == env_name { print $NF; exit }
    $NF ~ "/" env_name "$" { print $NF; exit }
  '
}

mkdir -p "$SUMMARY_DIR"

module load "$CONDA_MODULE"
module load "$R_MODULE"
source "$(conda info --base)/etc/profile.d/conda.sh"

if [ -z "$ENV_DIR" ]; then
  if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/python" ]; then
    ENV_DIR="$CONDA_PREFIX"
  else
    ENV_DIR="$PROJECT/envs/$ENV_NAME"
  fi
fi

if [ ! -x "$ENV_DIR/bin/python" ]; then
  detected_env_dir="$(find_conda_env_dir "$ENV_NAME" || true)"
  if [ -n "$detected_env_dir" ] && [ -x "$detected_env_dir/bin/python" ]; then
    ENV_DIR="$detected_env_dir"
  fi
fi

set +e
conda activate "$ENV_DIR" >/dev/null 2>&1
activate_status=$?
set -e
if [ "$activate_status" -ne 0 ]; then
  set +e
  conda activate "$ENV_NAME" >/dev/null 2>&1
  activate_status=$?
  set -e
fi
if [ "$activate_status" -eq 0 ] && [ -n "${CONDA_PREFIX:-}" ]; then
  ENV_DIR="$CONDA_PREFIX"
fi

if [ -x "$ENV_DIR/bin/python" ]; then
  PYTHON_BIN="$ENV_DIR/bin/python"
else
  PYTHON_BIN="$(command -v python || true)"
fi

if [ -z "$PYTHON_BIN" ] || [ ! -x "$PYTHON_BIN" ]; then
  echo "ERROR: Cannot find executable python after activating ENV_DIR=$ENV_DIR" >&2
  echo "Run these commands to find the real environment path, then pass ENV_DIR explicitly:" >&2
  echo "  module load $CONDA_MODULE" >&2
  echo "  source \"\$(conda info --base)/etc/profile.d/conda.sh\"" >&2
  echo "  conda activate $ENV_NAME" >&2
  echo "  echo \\\$CONDA_PREFIX" >&2
  echo "  ENV_DIR=\\\$CONDA_PREFIX bash $CODE_DIR/submit_04_trajectory_scvelo_array.sh" >&2
  exit 1
fi

"$PYTHON_BIN" -c "import scvelo, scanpy, anndata, loompy; print('Using python:', '$PYTHON_BIN'); print('scvelo:', scvelo.__version__)"

RSCRIPT_BIN="$(command -v Rscript || true)"
if [ -z "$RSCRIPT_BIN" ] || [ ! -x "$RSCRIPT_BIN" ]; then
  echo "ERROR: Rscript is not available after module load $R_MODULE" >&2
  echo "Try: module load $R_MODULE && which Rscript" >&2
  exit 1
fi

echo "scVelo array submit settings:"
echo "  PROJECT=$PROJECT"
echo "  TRAJECTORY_SCRIPT=$TRAJECTORY_SCRIPT"
echo "  TRAJECTORY_RESULT_NAME=$TRAJECTORY_RESULT_NAME"
echo "  ENV_DIR=$ENV_DIR"
echo "  PYTHON_BIN=$PYTHON_BIN"
echo "  RSCRIPT_BIN=$RSCRIPT_BIN"
echo "  SCVELO_N_JOBS=$SCVELO_N_JOBS"
echo "  SCVELO_MEM=$SCVELO_MEM"
echo "  QOS=$QOS"
echo "  LOOM_ROOT=$LOOM_ROOT"
echo "  TASK_LIST=$TASK_LIST"
echo "  SKIP_PREP=$SKIP_PREP"

if [ "$SKIP_PREP" != "TRUE" ]; then
  echo "Preparing scVelo task scripts. This does not run scVelo; it writes run_scvelo.sh files."
  cd "$PROJECT"

  RUN_SCVELO=FALSE \
  RUN_VELOCYTO=FALSE \
  SCVELO_N_JOBS="$SCVELO_N_JOBS" \
  SCVELO_PYTHON="$PYTHON_BIN" \
  VELOCYTO_OUTPUT_ROOT="$LOOM_ROOT" \
  VELOCYTO_WORK_ROOT="$LOOM_ROOT" \
  "$RSCRIPT_BIN" "$CODE_DIR/$TRAJECTORY_SCRIPT"

  if [ ! -f "$SUMMARY_FILE" ]; then
    echo "ERROR: Missing summary file after preparation: $SUMMARY_FILE" >&2
    exit 1
  fi

  echo "Writing scVelo task list for all prepared trajectory analyses."
  SUMMARY_FILE="$SUMMARY_FILE" TASK_LIST="$TASK_LIST" "$RSCRIPT_BIN" -e '
  x <- read.csv(Sys.getenv("SUMMARY_FILE"), stringsAsFactors = FALSE)
  x <- x[x$input_exists == TRUE & x$status %in% c("prepared", "prepared_only"), ]
  writeLines(unique(x$run_script), Sys.getenv("TASK_LIST"))
  cat("Tasks written:\n")
  keep <- intersect(c("analysis_group", "tn_scope", "ploidy_scope", "analysis_type", "n_cells", "run_script"), names(x))
  print(x[, keep, drop = FALSE], row.names = FALSE)
  '
else
  echo "Skipping preparation and submitting existing task list."
fi

cat > "$SBATCH_WORKER" <<'EOF'
#!/usr/bin/env bash
#SBATCH -J scvelo_gem
#SBATCH -o slurm_scvelo_%A_%a.out
#SBATCH -e slurm_scvelo_%A_%a.err

set -euo pipefail

: "${TASK_LIST:?TASK_LIST is required}"
: "${ENV_DIR:?ENV_DIR is required}"
: "${CONDA_MODULE:?CONDA_MODULE is required}"

module load "$CONDA_MODULE"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

RUN_SCRIPT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")
if [ -z "$RUN_SCRIPT" ]; then
  echo "ERROR: No run script found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID in $TASK_LIST" >&2
  exit 2
fi
if [ ! -f "$RUN_SCRIPT" ]; then
  echo "ERROR: run script does not exist: $RUN_SCRIPT" >&2
  exit 2
fi

echo "JOB_ID=$SLURM_JOB_ID"
echo "ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "HOST=$(hostname)"
echo "RUN_SCRIPT=$RUN_SCRIPT"
echo "START=$(date)"

bash "$RUN_SCRIPT"

echo "END=$(date)"
EOF

submit_array() {
  local task_list="$1"
  local cpus="$2"
  local mem="$3"
  local label="$4"
  local n_tasks

  n_tasks=$(grep -c . "$task_list" || true)
  if [ "$n_tasks" -eq 0 ]; then
    echo "No $label tasks to submit."
    return 0
  fi

  echo "Submitting $label: n=$n_tasks cpus=$cpus mem=$mem partition=${PARTITION:-default} time=$TASK_TIME"
  local sbatch_args=()
  if [ -n "$PARTITION" ]; then
    sbatch_args+=(--partition="$PARTITION")
  fi
  if [ -n "$QOS" ]; then
    sbatch_args+=(--qos="$QOS")
  fi
  sbatch \
    "${sbatch_args[@]}" \
    --array="1-${n_tasks}" \
    --cpus-per-task="$cpus" \
    --mem="$mem" \
    --time="$TASK_TIME" \
    --export=ALL,TASK_LIST="$task_list",ENV_DIR="$ENV_DIR",CONDA_MODULE="$CONDA_MODULE" \
    "$SBATCH_WORKER"
}

if [ ! -f "$TASK_LIST" ]; then
  echo "ERROR: Missing task list: $TASK_LIST" >&2
  echo "Run without SKIP_PREP=TRUE first, or create the task list from $SUMMARY_FILE." >&2
  exit 1
fi

submit_array "$TASK_LIST" "$SCVELO_CPUS" "$SCVELO_MEM" "velocity"

echo "Submitted available scVelo task array."
echo "Check queue: squeue -u \$USER"
