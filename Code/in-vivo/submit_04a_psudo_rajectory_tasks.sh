#!/usr/bin/env bash
set -euo pipefail

# Submit 04a_psudo_rajectory.R as ROOT/group-level Monocle3 Slurm tasks.

PROJECT="${PROJECT:-/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels}"
CODE_DIR="${CODE_DIR:-$PROJECT/Code/in-vivo}"
RESULTS="${RESULTS:-$PROJECT/Results}"
R_MODULE="${R_MODULE:-R/4.4}"
PARTITION="${PARTITION:-}"
QOS="${QOS:-xlarge}"
TASK_TIME="${TASK_TIME:-04:00:00}"
TASK_CPUS="${TASK_CPUS:-4}"
TASK_MEM="${TASK_MEM:-128G}"
MAX_ARRAY_JOBS="${MAX_ARRAY_JOBS:-}"
DRY_RUN="${DRY_RUN:-FALSE}"

PSEUDOTRAJ_ROOT="$RESULTS/04a_psudo_rajectory"
TASK_DIR="${TASK_DIR:-$PSEUDOTRAJ_ROOT/hpc_tasks}"
MANIFEST="${MANIFEST:-$TASK_DIR/04a_psudo_rajectory_hpc_tasks.tsv}"
SBATCH_WORKER="${SBATCH_WORKER:-$TASK_DIR/04a_psudo_rajectory_hpc_worker.sbatch}"
LOG_DIR="${LOG_DIR:-$TASK_DIR/slurm_logs}"

mkdir -p "$TASK_DIR" "$LOG_DIR"

module load "$R_MODULE"

RSCRIPT_BIN="$(command -v Rscript || true)"
if [ -z "$RSCRIPT_BIN" ] || [ ! -x "$RSCRIPT_BIN" ]; then
  echo "ERROR: Rscript is not available after module load $R_MODULE" >&2
  exit 1
fi

echo "Writing 04a_psudo_rajectory HPC task manifest: $MANIFEST"
CODE_DIR="$CODE_DIR" MANIFEST="$MANIFEST" "$RSCRIPT_BIN" - <<'RSCRIPT'
code_dir <- Sys.getenv("CODE_DIR")
manifest <- Sys.getenv("MANIFEST")

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
  branch <- paste0("ROOT_", sanitize_path_component(root))
  for (i in seq_len(nrow(groups))) {
    task_id <- task_id + 1L
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(task_id = task_id, method = "monocle3", root_cluster = root, stringsAsFactors = FALSE),
      groups[i, , drop = FALSE],
      data.frame(trajectory_branch = branch, stringsAsFactors = FALSE)
    )
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
#SBATCH -J mono04a
#SBATCH -o 04a_psudo_%A_%a.out
#SBATCH -e 04a_psudo_%A_%a.err

set -euo pipefail

: "${MANIFEST:?MANIFEST is required}"
: "${CODE_DIR:?CODE_DIR is required}"
: "${PROJECT:?PROJECT is required}"
: "${R_MODULE:?R_MODULE is required}"

module load "$R_MODULE"

line=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR == n + 1 { print; exit }' "$MANIFEST")
if [ -z "$line" ]; then
  echo "ERROR: No manifest row for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

IFS=$'\t' read -r task_id method root_cluster trajectory_group tn_scope ploidy_scope group_safe trajectory_branch <<< "$line"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export TRAJECTORY_HPC_TASK=TRUE
export TRAJECTORY_TASK_METHOD="$method"
export TRAJECTORY_TASK_ROOT_CLUSTER="$root_cluster"
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
echo "GROUP_SAFE=$group_safe"
echo "BRANCH=$trajectory_branch"
echo "START=$(date)"

Rscript "$CODE_DIR/04a_psudo_rajectory_hpc_task.R"

echo "END=$(date)"
EOF

n_tasks=$(awk 'END { if (NR > 0) print NR - 1; else print 0 }' "$MANIFEST")
if [ "$n_tasks" -le 0 ]; then
  echo "No tasks were written to $MANIFEST" >&2
  exit 1
fi

echo "04a_psudo_rajectory HPC submit settings:"
echo "  PROJECT=$PROJECT"
echo "  CODE_DIR=$CODE_DIR"
echo "  MANIFEST=$MANIFEST"
echo "  SBATCH_WORKER=$SBATCH_WORKER"
echo "  LOG_DIR=$LOG_DIR"
echo "  TASKS=$n_tasks"
echo "  TASK_CPUS=$TASK_CPUS"
echo "  TASK_MEM=$TASK_MEM"
echo "  TASK_TIME=$TASK_TIME"
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
  --export=ALL,MANIFEST="$MANIFEST",PROJECT="$PROJECT",CODE_DIR="$CODE_DIR",R_MODULE="$R_MODULE" \
  "$SBATCH_WORKER"

echo "Submitted 04a_psudo_rajectory task array."
