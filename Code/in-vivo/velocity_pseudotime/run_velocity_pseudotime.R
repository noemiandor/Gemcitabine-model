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

arg_flag <- function(args, name, default = FALSE) {
  normalized <- tolower(arg_value(args, name, if (default) "true" else "false"))
  if (!normalized %in% c("true", "false", "1", "0", "yes", "no")) {
    stop("--", name, " must be true or false", call. = FALSE)
  }
  normalized %in% c("true", "1", "yes")
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

run_r_stage <- function(label, script, values, log_path) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  cli <- paste0(
    "--", names(values), "=",
    vapply(values, as.character, character(1L))
  )
  message("[velocity-pseudotime] ", label)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), vapply(cli, shQuote, character(1L))),
    stdout = log_path,
    stderr = log_path
  )
  if (!identical(status, 0L)) {
    detail <- if (file.exists(log_path)) {
      paste(tail(readLines(log_path, warn = FALSE), 80L), collapse = "\n")
    } else {
      "No stage log was written."
    }
    stop(label, " failed; log: ", log_path, "\n", detail, call. = FALSE)
  }
  invisible(log_path)
}

bundle_is_valid <- function(table_path, provenance_path) {
  if (!file.exists(table_path) || !file.exists(provenance_path)) return(FALSE)
  tryCatch(
    {
      velocity_validate_panel_table(velocity_read_tsv(table_path))
      velocity_validate_provenance(provenance_path, table_path)
      TRUE
    },
    error = function(error) FALSE
  )
}

write_run_config <- function(
  output_dir,
  repo_root,
  mode,
  source_kind,
  table_path,
  provenance_path
) {
  table_locator <- velocity_portable_locator(table_path, repo_root)
  provenance_locator <- velocity_portable_locator(provenance_path, repo_root)
  canonical_publication_allowed <- identical(
    source_kind,
    "reviewed_frozen_table"
  ) && velocity_is_canonical_bundle(table_path, provenance_path, repo_root)
  data <- data.frame(
    key = c(
      "module", "mode", "source_kind", "source_table",
      "source_table_sha256", "source_provenance",
      "source_provenance_sha256", "root_cluster", "n_cells",
      "scvelo_mode", "canonical_publication_allowed"
    ),
    value = c(
      "in_vivo_velocity_pseudotime", mode, source_kind,
      table_locator, velocity_sha256(table_path), provenance_locator,
      velocity_sha256(provenance_path), "6", "2881", "stochastic",
      tolower(as.character(canonical_publication_allowed))
    ),
    stringsAsFactors = FALSE
  )
  velocity_write_tsv(data, file.path(output_dir, "metadata", "run_config.tsv"))
  canonical_publication_allowed
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  own_path <- script_path()
  module_dir <- dirname(own_path)
  repo_root <- normalizePath(file.path(module_dir, "..", "..", ".."), mustWork = TRUE)
  source(file.path(module_dir, "velocity_pseudotime_panel.R"), local = TRUE)
  mode <- arg_value(args, "mode", "plot-only")
  if (!mode %in% c("plot-only", "full-workflow")) {
    stop("--mode must be plot-only or full-workflow", call. = FALSE)
  }
  output_dir <- resolve_path(
    arg_value(args, "output-dir", required = TRUE),
    repo_root
  )
  intermediate_dir <- resolve_path(
    arg_value(
      args,
      "intermediate-dir",
      "Results/in-vivo/velocity_pseudotime/intermediates"
    ),
    repo_root
  )
  frozen_table <- resolve_path(
    arg_value(
      args,
      "frozen-table",
      "Data/in-vivo/figure7/processed/cellcycle_velocity_pseudotime_umap.tsv"
    ),
    repo_root
  )
  frozen_provenance <- resolve_path(
    arg_value(
      args,
      "frozen-provenance",
      "Data/in-vivo/figure7/processed/cellcycle_velocity_pseudotime_umap.provenance.tsv"
    ),
    repo_root
  )
  generated_table <- file.path(
    intermediate_dir,
    "cellcycle_velocity_pseudotime_umap.tsv"
  )
  generated_provenance <- file.path(
    intermediate_dir,
    "cellcycle_velocity_pseudotime_umap.provenance.tsv"
  )
  embedding_path <- file.path(intermediate_dir, "scvelo_velocity_umap.tsv")
  embedding_lineage_path <- file.path(
    intermediate_dir,
    "scvelo_velocity_umap.lineage.tsv"
  )
  embedding_inventory_path <- file.path(
    intermediate_dir,
    "scvelo_velocity_umap.inputs.tsv"
  )
  metrics_path <- file.path(intermediate_dir, "scvelo_cell_metrics.csv")
  builder_path <- file.path(module_dir, "build_velocity_pseudotime_table.R")
  panel_logic_path <- file.path(module_dir, "velocity_pseudotime_panel.R")
  scvelo_generator_path <- file.path(
    repo_root,
    "Code/in-vivo/figure7/generate_scvelo_cell_metrics.R"
  )
  config_path <- file.path(repo_root, "Code/in-vivo/figure7/figure7_config.yaml")
  environment_lock_path <- file.path(
    repo_root,
    "Code/in-vivo/figure7/environment_lock.tsv"
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
  source_kind <- ""
  table_path <- ""
  provenance_path <- ""
  if (bundle_is_valid(frozen_table, frozen_provenance)) {
    source_kind <- if (velocity_is_canonical_bundle(
      frozen_table,
      frozen_provenance,
      repo_root
    )) {
      "reviewed_frozen_table"
    } else {
      "noncanonical_frozen_override"
    }
    table_path <- frozen_table
    provenance_path <- frozen_provenance
  } else if (identical(mode, "plot-only")) {
    stop(
      "The reviewed frozen velocity panel bundle is absent or invalid. ",
      "Use --mode=full-workflow with the pinned loom/RDS inputs to build it.",
      call. = FALSE
    )
  } else {
    dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
    raw_data_dir <- resolve_path(
      arg_value(
        args,
        "raw-data-dir",
        "Results/in-vivo/figure7/raw/zenodo_21463392"
      ),
      repo_root
    )
    explicit_seurat <- arg_value(args, "seurat-rds", "")
    explicit_loom <- arg_value(args, "loom-root", "")
    seurat_rds <- resolve_path(
      if (nzchar(explicit_seurat)) {
        explicit_seurat
      } else {
        file.path(raw_data_dir, "integrated_sct_cca_seurat_final_reclustered.rds")
      },
      repo_root
    )
    loom_root <- resolve_path(
      if (nzchar(explicit_loom)) {
        explicit_loom
      } else {
        file.path(raw_data_dir, "velocyto_loom")
      },
      repo_root
    )
    loom_count <- if (dir.exists(loom_root)) {
      length(list.files(
        loom_root,
        pattern = "[.]loom$",
        recursive = TRUE,
        full.names = TRUE
      ))
    } else {
      0L
    }
    if ((!file.exists(seurat_rds) || loom_count != 18L) &&
        !nzchar(explicit_seurat) && !nzchar(explicit_loom) &&
        arg_flag(args, "download-missing-raw", TRUE)) {
      run_r_stage(
        "Download and validate pinned Figure 7 velocity inputs",
        file.path(repo_root, "Code/in-vivo/figure7/download_figure7_raw_data.R"),
        list(
          "raw-data-dir" = raw_data_dir,
          roles = "all",
          "allow-download" = "true"
        ),
        file.path(intermediate_dir, "logs", "00_download_raw.log")
      )
      loom_count <- length(list.files(
        loom_root,
        pattern = "[.]loom$",
        recursive = TRUE,
        full.names = TRUE
      ))
    }
    if (!file.exists(seurat_rds) || !dir.exists(loom_root) || loom_count != 18L) {
      stop(
        "Full velocity regeneration requires the selected Seurat RDS and ",
        "exactly 18 pinned loom files. Observed loom files: ", loom_count,
        call. = FALSE
      )
    }
    overwrite_intermediates <- arg_flag(args, "overwrite-intermediates", FALSE)
    generated_cache_current <- FALSE
    if (!overwrite_intermediates &&
        bundle_is_valid(generated_table, generated_provenance)) {
      generated_cache_current <- tryCatch(
        {
          velocity_generated_bundle_matches_dependencies(
            generated_table,
            generated_provenance,
            embedding_path,
            embedding_lineage_path,
            embedding_inventory_path,
            cellcycle_path,
            si_metadata_path,
            builder_path,
            panel_logic_path,
            scvelo_generator_path,
            seurat_rds,
            loom_root,
            config_path,
            environment_lock_path,
            repo_root
          )
          TRUE
        },
        error = function(error) {
          message(
            "Ignoring stale generated velocity table cache: ",
            conditionMessage(error)
          )
          FALSE
        }
      )
    }
    if (generated_cache_current) {
      source_kind <- "validated_generated_table_cache"
      table_path <- generated_table
      provenance_path <- generated_provenance
    } else {
      embedding_cache_current <- FALSE
      if (!overwrite_intermediates) {
        embedding_cache_current <- tryCatch(
          {
            velocity_validate_embedding_lineage(
              embedding_path,
              embedding_lineage_path,
              embedding_inventory_path,
              seurat_rds,
              loom_root,
              scvelo_generator_path,
              config_path,
              environment_lock_path,
              repo_root
            )
            TRUE
          },
          error = function(error) {
            message(
              "Ignoring absent or stale scVelo embedding cache: ",
              conditionMessage(error)
            )
            FALSE
          }
        )
      }
      if (!embedding_cache_current) {
        run_r_stage(
          "Compute stochastic scVelo pseudotime and UMAP velocity vectors",
          scvelo_generator_path,
          list(
            seurat_rds = seurat_rds,
            loom_root = loom_root,
            output = metrics_path,
            velocity_embedding_output = embedding_path,
            config = config_path,
            python = arg_value(args, "python", Sys.which("python3")),
            n_jobs = arg_value(args, "jobs", "1"),
            work_dir = file.path(intermediate_dir, "scvelo_work"),
            root_clusters = "6"
          ),
          file.path(intermediate_dir, "logs", "01_scvelo.log")
        )
        velocity_write_embedding_lineage(
          embedding_path,
          embedding_lineage_path,
          embedding_inventory_path,
          seurat_rds,
          loom_root,
          scvelo_generator_path,
          config_path,
          environment_lock_path,
          repo_root
        )
      }
      run_r_stage(
        "Build the exact 2,881-cell plot-facing table",
        builder_path,
        list(
          "embedding-input" = embedding_path,
          "embedding-lineage" = embedding_lineage_path,
          "embedding-inventory" = embedding_inventory_path,
          "cellcycle-input" = cellcycle_path,
          "si-metadata-input" = si_metadata_path,
          output = generated_table,
          "provenance-output" = generated_provenance,
          "root-cluster" = "6"
        ),
        file.path(intermediate_dir, "logs", "02_panel_table.log")
      )
      if (!bundle_is_valid(generated_table, generated_provenance)) {
        stop("New velocity panel table failed its provenance contract", call. = FALSE)
      }
      source_kind <- if (embedding_cache_current) {
        "validated_embedding_rebuilt_unreviewed"
      } else {
        "raw_regenerated_unreviewed"
      }
      table_path <- generated_table
      provenance_path <- generated_provenance
    }
  }
  velocity_render_panel(table_path, provenance_path, output_dir)
  canonical_publication_allowed <- write_run_config(
    output_dir,
    repo_root,
    mode,
    source_kind,
    table_path,
    provenance_path
  )
  message(
    "Velocity/pseudotime panel complete; source=", source_kind,
    "; canonical_publication_allowed=", canonical_publication_allowed
  )
}

if (identical(environment(), globalenv())) main()
