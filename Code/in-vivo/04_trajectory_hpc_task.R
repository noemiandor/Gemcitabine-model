#!/usr/bin/env Rscript

resolve_current_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_match <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_match) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_match[1]), mustWork = TRUE)))
  }
  getwd()
}

require_env <- function(name) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    stop(name, " is required.", call. = FALSE)
  }
  value
}

script_dir <- resolve_current_script_dir()

method <- tolower(require_env("TRAJECTORY_TASK_METHOD"))
if (!(method %in% c("scvelo", "paga"))) {
  stop("TRAJECTORY_TASK_METHOD must be 'scvelo' or 'paga'.", call. = FALSE)
}
require_env("TRAJECTORY_TASK_ROOT_CLUSTER")
require_env("TRAJECTORY_TASK_GROUP_SAFE")
if (identical(method, "scvelo")) {
  require_env("TRAJECTORY_TASK_END_CLUSTER")
} else {
  Sys.setenv(TRAJECTORY_TASK_END_CLUSTER = "NULL")
}

Sys.setenv(TRAJECTORY_HPC_TASK = "TRUE")

source(file.path(script_dir, "04_trajectory.R"))
