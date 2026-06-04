#!/usr/bin/env bash
set -euo pipefail

# Submit 04_trajectory.R as ROOT/END/method/group-level Slurm tasks.

PROJECT="${PROJECT:-/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels}"
CODE_DIR="${CODE_DIR:-$PROJECT/Code/in-vivo}"
RESULTS="${RESULTS:-$PROJECT/Results}"
ENV_NAME="${ENV_NAME:-rna_velocity_py310}"
ENV_DIR="${ENV_DIR:-}"
CONDA_MODULE="${CONDA_MODULE:-Miniconda3/23.9.0-0}"
R_MODULE="${R_MODULE:-R/4.4}"
PARTITION="${PARTITION:-}"
QOS="${QOS:-xlarge}"
TASK_TIME="${TASK_TIME:-04:00:00}"
TASK_CPUS="${TASK_CPUS:-4}"
TASK_MEM="${TASK_MEM:-256G}"
SCVELO_N_JOBS="${SCVELO_N_JOBS:-$TASK_CPUS}"
MAX_ARRAY_JOBS="${MAX_ARRAY_JOBS:-}"
DRY_RUN="${DRY_RUN:-FALSE}"
TRAJECTORY_METHODS="${TRAJECTORY_METHODS:-scvelo,paga}"

TRAJECTORY_ROOT="$RESULTS/04_trajectory"
TASK_DIR="${TASK_DIR:-$TRAJECTORY_ROOT/00_inputs/hpc_tasks}"
MANIFEST="${MANIFEST:-$TASK_DIR/04_trajectory_hpc_tasks.tsv}"
SBATCH_WORKER="${SBATCH_WORKER:-$TASK_DIR/04_trajectory_hpc_worker.sbatch}"
LOG_DIR="${LOG_DIR:-$TASK_DIR/slurm_logs}"

find_conda_env_dir() {
  local env_name="$1"
  conda env list | awk -v env_name="$env_name" '
    $1 == env_name { print $NF; exit }
    $NF ~ "/" env_name "$" { print $NF; exit }
  '
}

mkdir -p "$TASK_DIR" "$LOG_DIR"

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

PYTHON_BIN="$ENV_DIR/bin/python"
if [ ! -x "$PYTHON_BIN" ]; then
  echo "ERROR: Cannot find executable python in ENV_DIR=$ENV_DIR" >&2
  exit 1
fi

"$PYTHON_BIN" -c "import scvelo, scanpy, anndata, loompy, igraph, leidenalg; print('Using python:', '$PYTHON_BIN')"

RSCRIPT_BIN="$(command -v Rscript || true)"
if [ -z "$RSCRIPT_BIN" ] || [ ! -x "$RSCRIPT_BIN" ]; then
  echo "ERROR: Rscript is not available after module load $R_MODULE" >&2
  exit 1
fi

echo "Writing 04_trajectory HPC task manifest: $MANIFEST"
CODE_DIR="$CODE_DIR" MANIFEST="$MANIFEST" TRAJECTORY_METHODS="$TRAJECTORY_METHODS" "$RSCRIPT_BIN" - <<'RSCRIPT'
code_dir <- Sys.getenv("CODE_DIR")
manifest <- Sys.getenv("MANIFEST")
methods <- tolower(unlist(strsplit(Sys.getenv("TRAJECTORY_METHODS", unset = "scvelo,paga"), "[,;|]+")))
methods <- trimws(methods)
methods <- methods[nzchar(methods)]
allowed_methods <- c("scvelo", "paga")
bad_methods <- setdiff(methods, allowed_methods)
if (length(bad_methods) > 0) {
  stop("Unsupported TRAJECTORY_METHODS: ", paste(bad_methods, collapse = ", "), call. = FALSE)
}

source(file.path(code_dir, "Utils.R"))
config <- load_in_vivo_config(file.path(code_dir, "in_vivo_config.yaml"))

split_regex <- paste0("[,;|", intToUtf8(c(65292, 65307, 12289)), "]+")
strip_quotes <- function(value) {
  quote_chars <- paste0("\"'", intToUtf8(c(8220, 8221, 8216, 8217)))
  gsub(paste0("^[", quote_chars, "]+|[", quote_chars, "]+$"), "", value)
}
parse_clusters <- function(value, default = character(0)) {
  if (is.null(value) || length(value) == 0) value <- default
  if (is.list(value) && !is.data.frame(value)) value <- unlist(value, recursive = TRUE, use.names = FALSE)
  value <- as.character(value)
  value <- value[!is.na(value)]
  parts <- unlist(strsplit(value, split_regex))
  parts <- strip_quotes(trimws(parts))
  parts <- gsub("\\s+", "", parts)
  parts <- parts[nzchar(parts) & !(toupper(parts) %in% c("NA", "NAN", "NULL", "NONE"))]
  unique(parts)
}
cfg_clusters <- function(config, name, env_name = NULL, default = character(0)) {
  if (!is.null(env_name)) {
    env_value <- trimws(Sys.getenv(env_name, unset = ""))
    if (nzchar(env_value)) return(parse_clusters(env_value))
  }
  parsed <- parse_clusters(config[[name]])
  if (length(parsed) > 0) return(parsed)
  parse_clusters(default)
}
cfg_value_local <- function(config, name, default = NULL) {
  value <- config[[name]]
  if (!is.null(value) && length(value) > 0 && nzchar(trimws(as.character(value[1])))) {
    return(as.character(value[1]))
  }
  default
}

roots <- cfg_clusters(
  config,
  "Trajectory_root_clusters",
  "TRAJECTORY_ROOT_CLUSTERS",
  default = cfg_value_local(config, "PseudoTrajectory_root_cluster", "14")
)
ends <- cfg_clusters(config, "Trajectory_end_clusters", "TRAJECTORY_END_CLUSTERS", default = "13")
if (length(roots) == 0) stop("No Trajectory_root_clusters resolved.", call. = FALSE)

groups <- data.frame(
  trajectory_group = c(
    "All_cells",
    "CellLine_ploidy_all", "CellLine_2N", "CellLine_4N",
    "Tumor_ploidy_all", "Tumor_2N", "Tumor_4N"
  ),
  tn_scope = c("All", "CellLine", "CellLine", "CellLine", "Tumor", "Tumor", "Tumor"),
  ploidy_scope = c("ploidy_all", "ploidy_all", "2N", "4N", "ploidy_all", "2N", "4N"),
  group_safe = c(
    "All_cells",
    "CellLine/ploidy_all", "CellLine/2N", "CellLine/4N",
    "Tumor/ploidy_all", "Tumor/2N", "Tumor/4N"
  ),
  stringsAsFactors = FALSE
)

rows <- list()
task_id <- 0L
for (root in roots) {
  scvelo_ends <- c("NULL", ends)
  if ("scvelo" %in% methods) {
    for (end in scvelo_ends) {
      end_label <- if (identical(end, "NULL")) "NULL" else end
      branch <- paste0("ROOT_", sanitize_path_component(root), "_END_", sanitize_path_component(end_label))
      for (i in seq_len(nrow(groups))) {
        task_id <- task_id + 1L
        rows[[length(rows) + 1L]] <- cbind(
          data.frame(task_id = task_id, method = "scvelo", root_cluster = root, end_cluster = end_label, stringsAsFactors = FALSE),
          groups[i, , drop = FALSE],
          data.frame(trajectory_branch = branch, stringsAsFactors = FALSE)
        )
      }
    }
  }
  if ("paga" %in% methods) {
    branch <- paste0("ROOT_", sanitize_path_component(root), "_END_NULL")
    for (i in seq_len(nrow(groups))) {
      task_id <- task_id + 1L
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(task_id = task_id, method = "paga", root_cluster = root, end_cluster = "NULL", stringsAsFactors = FALSE),
        groups[i, , drop = FALSE],
        data.frame(trajectory_branch = branch, stringsAsFactors = FALSE)
      )
    }
  }
}

manifest_df <- do.call(rbind, rows)
dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)
write.table(manifest_df, manifest, sep = "\t", quote = FALSE, row.names = FALSE)
print(manifest_df, row.names = FALSE)
cat("Wrote ", nrow(manifest_df), " tasks to ", manifest, "\n", sep = "")
RSCRIPT

cat > "$SBATCH_WORKER" <<'EOF'
#!/usr/bin/env bash
#SBATCH -J traj04
#SBATCH -o 04_trajectory_%A_%a.out
#SBATCH -e 04_trajectory_%A_%a.err

set -euo pipefail

: "${MANIFEST:?MANIFEST is required}"
: "${CODE_DIR:?CODE_DIR is required}"
: "${PROJECT:?PROJECT is required}"
: "${ENV_DIR:?ENV_DIR is required}"
: "${CONDA_MODULE:?CONDA_MODULE is required}"
: "${R_MODULE:?R_MODULE is required}"
: "${SCVELO_N_JOBS:?SCVELO_N_JOBS is required}"

module load "$CONDA_MODULE"
module load "$R_MODULE"
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda activate "$ENV_DIR" >/dev/null 2>&1; then
  if [ -f "$ENV_DIR/bin/activate" ]; then
    source "$ENV_DIR/bin/activate"
  elif [ -x "$ENV_DIR/bin/python" ]; then
    export PATH="$ENV_DIR/bin:$PATH"
  else
    echo "ERROR: Cannot activate ENV_DIR=$ENV_DIR and no executable python was found there." >&2
    exit 1
  fi
fi

line=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR == n + 1 { print; exit }' "$MANIFEST")
if [ -z "$line" ]; then
  echo "ERROR: No manifest row for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

IFS=$'\t' read -r task_id method root_cluster end_cluster trajectory_group tn_scope ploidy_scope group_safe trajectory_branch <<< "$line"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export SCVELO_PYTHON="$ENV_DIR/bin/python"
export SCVELO_N_JOBS="$SCVELO_N_JOBS"
export RUN_VELOCYTO="${RUN_VELOCYTO:-FALSE}"
export TRAJECTORY_HPC_TASK=TRUE
export TRAJECTORY_TASK_METHOD="$method"
export TRAJECTORY_TASK_ROOT_CLUSTER="$root_cluster"
export TRAJECTORY_TASK_END_CLUSTER="$end_cluster"
export TRAJECTORY_TASK_GROUP="$trajectory_group"
export TRAJECTORY_TASK_TN_SCOPE="$tn_scope"
export TRAJECTORY_TASK_PLOIDY_SCOPE="$ploidy_scope"
export TRAJECTORY_TASK_GROUP_SAFE="$group_safe"

cd "$PROJECT"

echo "JOB_ID=$SLURM_JOB_ID"
echo "ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "HOST=$(hostname)"
echo "TASK_ID=$task_id"
echo "METHOD=$method"
echo "ROOT=$root_cluster"
echo "END=$end_cluster"
echo "GROUP_SAFE=$group_safe"
echo "BRANCH=$trajectory_branch"
echo "START=$(date)"

Rscript "$CODE_DIR/04_trajectory_hpc_task.R"

echo "END=$(date)"
EOF

n_tasks=$(awk 'END { if (NR > 0) print NR - 1; else print 0 }' "$MANIFEST")
if [ "$n_tasks" -le 0 ]; then
  echo "No tasks were written to $MANIFEST" >&2
  exit 1
fi

echo "04_trajectory HPC submit settings:"
echo "  PROJECT=$PROJECT"
echo "  CODE_DIR=$CODE_DIR"
echo "  ENV_DIR=$ENV_DIR"
echo "  MANIFEST=$MANIFEST"
echo "  SBATCH_WORKER=$SBATCH_WORKER"
echo "  LOG_DIR=$LOG_DIR"
echo "  TASKS=$n_tasks"
echo "  TASK_CPUS=$TASK_CPUS"
echo "  TASK_MEM=$TASK_MEM"
echo "  TASK_TIME=$TASK_TIME"
echo "  TRAJECTORY_METHODS=$TRAJECTORY_METHODS"
echo "  DRY_RUN=$DRY_RUN"

if [ "$DRY_RUN" = "TRUE" ]; then
  echo "DRY_RUN=TRUE; manifest and worker were written but sbatch was not called."
  exit 0
fi

array_spec="1-${n_tasks}"
if [ -n "$MAX_ARRAY_JOBS" ]; then
  array_spec="${array_spec}%${MAX_ARRAY_JOBS}"
fi

sbatch_args=()
if [ -n "$PARTITION" ]; then
  sbatch_args+=(--partition="$PARTITION")
fi
if [ -n "$QOS" ]; then
  sbatch_args+=(--qos="$QOS")
fi

sbatch \
  "${sbatch_args[@]}" \
  --array="$array_spec" \
  --cpus-per-task="$TASK_CPUS" \
  --mem="$TASK_MEM" \
  --time="$TASK_TIME" \
  --chdir="$LOG_DIR" \
  --export=ALL,MANIFEST="$MANIFEST",PROJECT="$PROJECT",CODE_DIR="$CODE_DIR",ENV_DIR="$ENV_DIR",CONDA_MODULE="$CONDA_MODULE",R_MODULE="$R_MODULE",SCVELO_N_JOBS="$SCVELO_N_JOBS" \
  "$SBATCH_WORKER"

echo "Submitted 04_trajectory task array."
