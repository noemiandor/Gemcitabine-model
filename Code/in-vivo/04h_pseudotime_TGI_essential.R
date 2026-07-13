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
cellcycle_input_provided <- !is.null(arguments[["cellcycle_input"]]) &&
  length(arguments[["cellcycle_input"]]) == 1L && nzchar(arguments[["cellcycle_input"]])
noncellcycle_input_provided <- !is.null(arguments[["noncellcycle_input"]]) &&
  length(arguments[["noncellcycle_input"]]) == 1L && nzchar(arguments[["noncellcycle_input"]])
generated_input_root <- normalizePath(
  essential_arg(arguments, "generated_input_root", input_root),
  mustWork = FALSE
)
cellcycle_input <- normalizePath(
  essential_arg(
    arguments,
    "cellcycle_input",
    file.path(
      generated_input_root,
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
      generated_input_root,
      "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    )
  ),
  mustWork = FALSE
)
scvelo_metrics_input <- normalizePath(
  essential_arg(arguments, "scvelo_metrics_input", file.path(input_root, "scvelo_cell_metrics.csv")),
  mustWork = FALSE
)
cell_ploidy_input <- normalizePath(
  essential_arg(arguments, "cell_ploidy_input", file.path(input_root, "all_ploidy.tsv")),
  mustWork = FALSE
)
sample_info_input <- normalizePath(
  essential_arg(arguments, "sample_info_input", file.path(input_root, "sample_info.xlsx")),
  mustWork = FALSE
)
growth_curve_input <- normalizePath(
  essential_arg(
    arguments,
    "growth_curve_input",
    file.path(input_root, "dt_Gem_VT_20241223_v4.xlsx")
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
figures_only <- essential_bool(essential_arg(arguments, "figures_only", "FALSE"))
tables_root <- normalizePath(
  essential_arg(arguments, "tables_root", output_root),
  mustWork = FALSE
)
seed <- as.integer(essential_arg(arguments, "seed", "1"))
n_perm <- as.integer(essential_arg(arguments, "n_perm", "10000"))
n_boot <- as.integer(essential_arg(arguments, "n_boot", "1000"))
workers <- as.integer(essential_arg(arguments, "workers", "4"))
overwrite <- essential_bool(essential_arg(arguments, "overwrite", "FALSE"))
internal_method <- essential_arg(arguments, "internal_method", NULL)
internal_figures_method <- essential_arg(arguments, "internal_figures_method", NULL)

if (!is.null(internal_method) && !is.null(internal_figures_method)) {
  stop("Only one internal method mode can be selected", call. = FALSE)
}
if (figures_only) {
  numeric_parameters <- c(seed = seed, workers = workers)
  if (anyNA(numeric_parameters) || seed < 0L || workers < 1L) {
    stop("seed and workers must be valid positive integers", call. = FALSE)
  }
  if (!dir.exists(tables_root)) stop("Missing figures-only tables root: ", tables_root, call. = FALSE)
} else {
  numeric_parameters <- c(seed = seed, n_perm = n_perm, n_boot = n_boot, workers = workers)
  if (anyNA(numeric_parameters) || seed < 0L || n_perm < 1L || n_boot < 1L || workers < 1L) {
    stop("seed, n_perm, n_boot, and workers must be valid positive integers", call. = FALSE)
  }
  if (!cellcycle_input_provided || !noncellcycle_input_provided) {
    generated <- essential_generate_analysis_inputs(
      scvelo_path = scvelo_metrics_input,
      cell_ploidy_path = cell_ploidy_input,
      sample_info_path = sample_info_input,
      growth_curve_path = growth_curve_input,
      cellcycle_output = if (!cellcycle_input_provided) cellcycle_input else NULL,
      noncellcycle_output = if (!noncellcycle_input_provided) noncellcycle_input else NULL
    )
    generated_names <- names(generated)
    message(
      "Generated missing cell-level analysis input(s) from Data/in-vivo sources: ",
      paste(generated_names, collapse = ", ")
    )
  }
  if (!file.exists(cellcycle_input)) stop("Missing CellCycle input: ", cellcycle_input, call. = FALSE)
  if (!file.exists(noncellcycle_input)) stop("Missing NonCellCycle input: ", noncellcycle_input, call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required", call. = FALSE)

method_specs <- essential_method_specs()

if (!is.null(internal_figures_method)) {
  if (!figures_only) stop("Internal figures method requires --figures_only=TRUE", call. = FALSE)
  if (!internal_figures_method %in% names(method_specs)) {
    stop("Unknown internal figures method: ", internal_figures_method, call. = FALSE)
  }
  essential_regenerate_figures_from_tables(
    tables_root = tables_root,
    output_root = output_root,
    spec = method_specs[[internal_figures_method]]
  )
  message("Regenerated essential figures: ", internal_figures_method)
  quit(save = "no", status = 0L)
}

if (!is.null(internal_method)) {
  if (figures_only) stop("Internal analysis method cannot run in figures-only mode", call. = FALSE)
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

if (figures_only) {
  target_dirs <- file.path(output_root, "Figures", selected_methods)
  nonempty_targets <- vapply(target_dirs, function(path) {
    dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0L
  }, logical(1L))
  if (any(nonempty_targets) && !overwrite) {
    stop(
      "Figure output already exists. Use --overwrite=TRUE to replace selected method figures: ",
      paste(selected_methods[nonempty_targets], collapse = ", "),
      call. = FALSE
    )
  }
  if (overwrite) {
    existing_targets <- target_dirs[dir.exists(target_dirs) | file.exists(target_dirs)]
    if (length(existing_targets) > 0L) unlink(existing_targets, recursive = TRUE, force = TRUE)
  }
  essential_ensure_dir(file.path(output_root, "Figures"))

  rscript <- file.path(R.home("bin"), "Rscript")
  run_figures_child <- function(method, rscript, script_path, tables_root, output_root) {
    child_arguments <- c(
      shQuote(script_path),
      "--figures_only=TRUE",
      shQuote(paste0("--internal_figures_method=", method)),
      shQuote(paste0("--tables_root=", tables_root)),
      shQuote(paste0("--output_root=", output_root))
    )
    status <- system2(rscript, args = child_arguments)
    as.integer(status)
  }

  message(
    "Regenerating figures for ", length(selected_methods), " method(s) with ", workers,
    " isolated worker(s); no statistical tests will be run"
  )
  if (workers > 1L && length(selected_methods) > 1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    statuses <- unlist(parallel::parLapplyLB(
      cluster,
      selected_methods,
      run_figures_child,
      rscript = rscript,
      script_path = script_path,
      tables_root = tables_root,
      output_root = output_root
    ))
    parallel::stopCluster(cluster)
    on.exit(NULL, add = FALSE)
  } else {
    statuses <- vapply(selected_methods, function(method) {
      run_figures_child(method, rscript, script_path, tables_root, output_root)
    }, integer(1L))
  }
  names(statuses) <- selected_methods
  if (any(statuses != 0L)) {
    failed <- names(statuses)[statuses != 0L]
    stop("Figure-only worker(s) failed: ", paste(failed, collapse = ", "), call. = FALSE)
  }
  inventory <- essential_validate_figures_only_inventory(output_root, selected_methods)
  message("Completed figures-only regeneration from existing tables")
  message("Tables root: ", tables_root)
  message("Output root: ", output_root)
  message("Inventory: ", inventory$figures, " PDF; no statistics were recomputed or written")
  quit(save = "no", status = 0L)
}

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
