#!/usr/bin/env Rscript

file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_argument) > 0L) {
  normalizePath(sub("^--file=", "", file_argument[[1L]]), mustWork = FALSE)
} else {
  normalizePath("Code/in-vivo/04h_pseudotime_TGI_essential.R", mustWork = FALSE)
}
script_dir <- dirname(script_path)
utility_path <- file.path(script_dir, "04h_pseudotime_TGI_essential_util.R")
if (!file.exists(utility_path)) {
  stop("Missing essential utility script: ", utility_path, call. = FALSE)
}
sys.source(utility_path, envir = .GlobalEnv)

font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

arguments <- essential_parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
input_root <- normalizePath(
  essential_arg(arguments, "input_root", file.path(repo_root, "Data", "in-vivo")),
  mustWork = FALSE
)
cellcycle_input <- normalizePath(
  essential_arg(
    arguments,
    "cellcycle_input",
    file.path(
      input_root,
      "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    )
  ),
  mustWork = FALSE
)
noncellcycle_input <- normalizePath(
  essential_arg(
    arguments,
    "noncellcycle_input",
    file.path(
      input_root,
      "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    )
  ),
  mustWork = FALSE
)
output_root <- normalizePath(
  essential_arg(
    arguments,
    "output_root",
    "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04h_pseudotime_TGI_essential"
  ),
  mustWork = FALSE
)
seed <- as.integer(essential_arg(arguments, "seed", "1"))
n_perm <- as.integer(essential_arg(arguments, "n_perm", "10000"))
n_boot <- as.integer(essential_arg(arguments, "n_boot", "1000"))
workers <- as.integer(essential_arg(arguments, "workers", "4"))
overwrite <- essential_bool(essential_arg(arguments, "overwrite", "FALSE"))
internal_method <- essential_arg(arguments, "internal_method", NULL)

numeric_parameters <- c(seed = seed, n_perm = n_perm, n_boot = n_boot, workers = workers)
if (anyNA(numeric_parameters) || seed < 0L || n_perm < 1L || n_boot < 1L || workers < 1L) {
  stop("seed, n_perm, n_boot, and workers must be valid positive integers", call. = FALSE)
}
if (!file.exists(cellcycle_input)) stop("Missing CellCycle input: ", cellcycle_input, call. = FALSE)
if (!file.exists(noncellcycle_input)) stop("Missing NonCellCycle input: ", noncellcycle_input, call. = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required", call. = FALSE)

method_specs <- essential_method_specs()

if (!is.null(internal_method)) {
  if (!internal_method %in% names(method_specs)) {
    stop("Unknown internal method: ", internal_method, call. = FALSE)
  }
  essential_run_method(
    cellcycle_path = cellcycle_input,
    noncellcycle_path = noncellcycle_input,
    output_root = output_root,
    spec = method_specs[[internal_method]],
    seed = seed,
    n_perm = n_perm,
    n_boot = n_boot
  )
  message("Completed essential method: ", internal_method)
  quit(save = "no", status = 0L)
}

requested_methods <- essential_arg(arguments, "methods", "all")
if (tolower(requested_methods) == "all") {
  selected_methods <- names(method_specs)
} else {
  selected_methods <- trimws(strsplit(requested_methods, ",", fixed = TRUE)[[1L]])
  unknown <- setdiff(selected_methods, names(method_specs))
  if (length(unknown) > 0L) {
    stop("Unknown method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
}
selected_methods <- unique(selected_methods)
workers <- min(workers, length(selected_methods))

existing_files <- if (dir.exists(output_root)) {
  list.files(output_root, all.files = TRUE, no.. = TRUE)
} else {
  character()
}
if (length(existing_files) > 0L && !overwrite) {
  stop(
    "Output directory is not empty. Use --overwrite=TRUE to replace essential outputs: ",
    output_root,
    call. = FALSE
  )
}
if (overwrite && dir.exists(output_root)) {
  replace_targets <- file.path(output_root, c("Figures", "stats", "plot_data", "README.md"))
  existing_targets <- replace_targets[file.exists(replace_targets) | dir.exists(replace_targets)]
  if (length(existing_targets) > 0L) unlink(existing_targets, recursive = TRUE, force = TRUE)
}
essential_ensure_dir(output_root)
essential_ensure_dir(file.path(output_root, "Figures"))
essential_ensure_dir(file.path(output_root, "stats"))
essential_ensure_dir(file.path(output_root, "plot_data"))

rscript <- file.path(R.home("bin"), "Rscript")
run_child <- function(method, rscript, script_path, cellcycle_input, noncellcycle_input, output_root, seed, n_perm, n_boot) {
  child_arguments <- c(
    shQuote(script_path),
    shQuote(paste0("--internal_method=", method)),
    shQuote(paste0("--cellcycle_input=", cellcycle_input)),
    shQuote(paste0("--noncellcycle_input=", noncellcycle_input)),
    shQuote(paste0("--output_root=", output_root)),
    paste0("--seed=", seed),
    paste0("--n_perm=", n_perm),
    paste0("--n_boot=", n_boot)
  )
  status <- system2(rscript, args = child_arguments)
  as.integer(status)
}

message(
  "Running ", length(selected_methods), " essential method(s) with ", workers,
  " isolated worker(s)"
)
if (workers > 1L && length(selected_methods) > 1L) {
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  statuses <- unlist(parallel::parLapplyLB(
    cluster,
    selected_methods,
    run_child,
    rscript = rscript,
    script_path = script_path,
    cellcycle_input = cellcycle_input,
    noncellcycle_input = noncellcycle_input,
    output_root = output_root,
    seed = seed,
    n_perm = n_perm,
    n_boot = n_boot
  ))
  parallel::stopCluster(cluster)
  on.exit(NULL, add = FALSE)
} else {
  statuses <- vapply(selected_methods, function(method) {
    run_child(
      method,
      rscript,
      script_path,
      cellcycle_input,
      noncellcycle_input,
      output_root,
      seed,
      n_perm,
      n_boot
    )
  }, integer(1L))
}
names(statuses) <- selected_methods
if (any(statuses != 0L)) {
  failed <- names(statuses)[statuses != 0L]
  stop("Essential method worker(s) failed: ", paste(failed, collapse = ", "), call. = FALSE)
}

essential_write_readme(output_root)
inventory <- essential_validate_inventory(output_root, selected_methods)
message("Completed essential pseudotime-TGI analysis")
message("Output root: ", output_root)
message(
  "Inventory: ", inventory$figures, " PDF; ", inventory$stats,
  " statistics CSV; ", inventory$plot_data, " plot-data CSV"
)
