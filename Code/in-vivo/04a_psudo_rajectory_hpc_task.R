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

method <- tolower(Sys.getenv("TRAJECTORY_TASK_METHOD", unset = "monocle3"))
if (!(method %in% c("monocle3", "monocle"))) {
  stop("TRAJECTORY_TASK_METHOD must be 'monocle3'.", call. = FALSE)
}
require_env("TRAJECTORY_TASK_ROOT_CLUSTER")
require_env("TRAJECTORY_TASK_GROUP_SAFE")

Sys.setenv(
  TRAJECTORY_HPC_TASK = "TRUE",
  TRAJECTORY_TASK_METHOD = "monocle3"
)

source(file.path(script_dir, "04a_psudo_rajectory.R"))
