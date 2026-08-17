#!/usr/bin/env Rscript

parse_args <- function(tokens) {
  out <- list()
  for (token in tokens) {
    if (!startsWith(token, "--") || !grepl("=", token, fixed = TRUE)) {
      stop("Arguments must use --name=value syntax: ", token, call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
    out[[pieces[[1L]]]] <- paste(pieces[-1L], collapse = "=")
  }
  out
}

arg_value <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value) || length(value) != 1L || !nzchar(value)) value <- default
  if (isTRUE(required) && (is.null(value) || !nzchar(value))) {
    stop("Missing required argument --", name, call. = FALSE)
  }
  value
}

script_path <- function() {
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(hit)) return(NA_character_)
  normalizePath(sub("^--file=", "", hit[[1L]]), mustWork = TRUE)
}

resolve_path <- function(path, root, must_work = FALSE) {
  candidate <- if (grepl("^/", path)) path else file.path(root, path)
  normalizePath(candidate, mustWork = must_work)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  own_path <- script_path()
  module_dir <- dirname(own_path)
  repo_root <- normalizePath(file.path(module_dir, "..", "..", ".."), mustWork = TRUE)
  source(file.path(module_dir, "velocity_pseudotime_panel.R"), local = TRUE)
  embedding_path <- resolve_path(
    arg_value(args, "embedding-input", required = TRUE),
    repo_root,
    must_work = TRUE
  )
  embedding_lineage_path <- resolve_path(
    arg_value(args, "embedding-lineage", required = TRUE),
    repo_root,
    must_work = TRUE
  )
  embedding_inventory_path <- resolve_path(
    arg_value(args, "embedding-inventory", required = TRUE),
    repo_root,
    must_work = TRUE
  )
  cellcycle_path <- resolve_path(
    arg_value(
      args,
      "cellcycle-input",
      "Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    ),
    repo_root,
    must_work = TRUE
  )
  si_metadata_path <- resolve_path(
    arg_value(
      args,
      "si-metadata-input",
      "Data/in-vivo/SIfigures/si_figures_cell_metadata.csv"
    ),
    repo_root,
    must_work = TRUE
  )
  output_path <- resolve_path(
    arg_value(
      args,
      "output",
      "Data/in-vivo/figure7/processed/cellcycle_velocity_pseudotime_umap.tsv"
    ),
    repo_root
  )
  provenance_path <- resolve_path(
    arg_value(
      args,
      "provenance-output",
      if (grepl("[.]tsv$", output_path)) {
        sub("[.]tsv$", ".provenance.tsv", output_path)
      } else {
        paste0(output_path, ".provenance.tsv")
      }
    ),
    repo_root
  )
  root_cluster <- arg_value(args, "root-cluster", "6")
  table <- velocity_build_panel_table(
    embedding_path,
    cellcycle_path,
    si_metadata_path,
    root_cluster = root_cluster
  )
  velocity_write_panel_bundle(
    table,
    output_path,
    provenance_path,
    embedding_path,
    embedding_lineage_path,
    embedding_inventory_path,
    cellcycle_path,
    si_metadata_path,
    builder_path = own_path,
    panel_logic_path = file.path(module_dir, "velocity_pseudotime_panel.R"),
    scvelo_generator_path = file.path(
      repo_root,
      "Code/in-vivo/figure7/generate_scvelo_cell_metrics.R"
    ),
    root_cluster = root_cluster
  )
  message("Wrote CellCycle velocity panel table: ", output_path)
  message("Wrote CellCycle velocity provenance: ", provenance_path)
}

if (identical(environment(), globalenv())) main()
