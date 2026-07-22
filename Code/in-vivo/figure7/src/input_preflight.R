# End-to-end input preparation and preflight for the Figure 7 full workflow.

figure7_flag <- function(args, name, default = FALSE) {
  value <- figure7_arg(args, name, if (isTRUE(default)) "true" else "false")
  normalized <- tolower(trimws(as.character(value)))
  if (!normalized %in% c("true", "false", "1", "0", "yes", "no", "y", "n")) {
    figure7_stop("--", gsub("_", "-", name), " must be true or false")
  }
  normalized %in% c("true", "1", "yes", "y")
}

figure7_workflow_path <- function(path, repo_root) {
  if (is.null(path) || !nzchar(path)) return("")
  if (!grepl("^/", path)) path <- file.path(repo_root, path)
  normalizePath(path, mustWork = FALSE)
}

figure7_workflow_integer_arg <- function(args, name, default, maximum = 16L) {
  raw <- figure7_arg(args, name, as.character(default))
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || !is.finite(value) || value < 1 ||
      value > maximum || value != floor(value)) {
    figure7_stop("--", name, " must be an integer from 1 to ", maximum)
  }
  as.integer(value)
}

figure7_workflow_paths <- function(args, repo_root, output_dir, config) {
  intermediate_dir <- figure7_workflow_path(
    figure7_arg(args, "intermediate-dir", file.path("Results", "in-vivo", "figure7", "intermediates")),
    repo_root
  )
  run_key <- basename(output_dir)
  reference_id <- as.character(config$state_pathways$reference_id)
  raw_data_dir <- figure7_workflow_path(
    figure7_arg(args, "raw-data-dir", as.character(config$raw_data$default_root)),
    repo_root
  )
  loom_root_arg <- figure7_arg(args, "loom-root", "")
  seurat_rds_arg <- figure7_arg(args, "seurat-rds", "")
  list(
    intermediate_dir = intermediate_dir,
    raw_data_dir = raw_data_dir,
    raw_manifest = figure7_workflow_path(as.character(config$raw_data$required_manifest), repo_root),
    download_missing_raw = figure7_flag(args, "download-missing-raw", TRUE),
    download_workers = figure7_workflow_integer_arg(args, "download-workers", 4L),
    download_connections_per_file = figure7_workflow_integer_arg(
      args, "download-connections-per-file", 2L
    ),
    loom_root_explicit = nzchar(loom_root_arg),
    seurat_rds_explicit = nzchar(seurat_rds_arg),
    scvelo_metrics = figure7_workflow_path(
      figure7_arg(args, "scvelo-metrics", file.path(intermediate_dir, "scvelo_cell_metrics.csv")),
      repo_root
    ),
    seurat_metadata = figure7_workflow_path(
      figure7_arg(args, "seurat-metadata-output", file.path("Data", "in-vivo", "seurat_metadata.csv")),
      repo_root
    ),
    cellcycle = figure7_workflow_path(
      figure7_arg(
        args,
        "cellcycle-input",
        file.path(intermediate_dir, "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv")
      ),
      repo_root
    ),
    noncellcycle = figure7_workflow_path(
      figure7_arg(
        args,
        "non-cellcycle-input",
        file.path(intermediate_dir, "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv")
      ),
      repo_root
    ),
    loom_root = figure7_workflow_path(
      if (nzchar(loom_root_arg)) loom_root_arg else file.path(raw_data_dir, as.character(config$raw_data$loom_subdir)),
      repo_root
    ),
    seurat_rds = figure7_workflow_path(
      if (nzchar(seurat_rds_arg)) seurat_rds_arg else file.path(raw_data_dir, as.character(config$raw_data$seurat_filename)),
      repo_root
    ),
    cell_ploidy = figure7_workflow_path(
      figure7_arg(args, "cell-ploidy-input", file.path("Data", "in-vivo", "all_ploidy.tsv")),
      repo_root
    ),
    sample_info = figure7_workflow_path(
      figure7_arg(args, "sample-info-input", file.path("Data", "in-vivo", "sample_info.xlsx")),
      repo_root
    ),
    growth_curve = figure7_workflow_path(
      figure7_arg(args, "growth-curve-input", file.path("Data", "in-vivo", "dt_Gem_VT_20241223_v4.xlsx")),
      repo_root
    ),
    python = figure7_arg(args, "python", Sys.which("python")),
    scvelo_work_dir = file.path(intermediate_dir, "scvelo_work"),
    state_pathway_root = figure7_workflow_path(
      figure7_arg(args, "state-pathway-output-root", file.path(intermediate_dir, "state_pathways", run_key)),
      repo_root
    ),
    saved_reference = figure7_workflow_path(
      figure7_arg(
        args,
        "saved-state-pathway-dir",
        file.path(intermediate_dir, "saved_state_pathway", run_key, reference_id)
      ),
      repo_root
    ),
    log_dir = file.path(intermediate_dir, "logs", run_key)
  )
}

figure7_require_workflow_file <- function(path, argument) {
  if (!nzchar(path) || !file.exists(path)) {
    figure7_stop("Missing required full-workflow input --", argument, ": ", if (nzchar(path)) path else "not supplied")
  }
  invisible(path)
}

figure7_require_workflow_dir <- function(path, argument) {
  if (!nzchar(path) || !dir.exists(path)) {
    figure7_stop("Missing required full-workflow input --", argument, ": ", if (nzchar(path)) path else "not supplied")
  }
  invisible(path)
}

figure7_preflight_workflow <- function(paths, include_panel_f = TRUE, overwrite_intermediates = FALSE) {
  has_cellcycle <- file.exists(paths$cellcycle)
  has_noncellcycle <- file.exists(paths$noncellcycle)
  has_scvelo <- file.exists(paths$scvelo_metrics)

  if (xor(has_cellcycle, has_noncellcycle) && !isTRUE(overwrite_intermediates)) {
    present <- if (has_cellcycle) paths$cellcycle else paths$noncellcycle
    missing <- if (has_cellcycle) paths$noncellcycle else paths$cellcycle
    figure7_stop(
      "Incomplete cell-level intermediate pair. Present: ", present,
      "; missing: ", missing,
      ". Supply both files or rerun with --overwrite-intermediates=true."
    )
  }

  cell_pair_ready <- has_cellcycle && has_noncellcycle
  state <- if (cell_pair_ready) {
    "cell_tables_ready"
  } else if (has_scvelo) {
    "scvelo_ready"
  } else {
    "raw_inputs_ready"
  }

  if (!cell_pair_ready) {
    figure7_require_workflow_file(paths$cell_ploidy, "cell-ploidy-input")
    figure7_require_workflow_file(paths$sample_info, "sample-info-input")
    figure7_require_workflow_file(paths$growth_curve, "growth-curve-input")
  }
  needs_loom <- identical(state, "raw_inputs_ready")
  needs_seurat <- needs_loom || isTRUE(include_panel_f)
  missing_loom <- needs_loom && !dir.exists(paths$loom_root)
  missing_seurat <- needs_seurat && !file.exists(paths$seurat_rds)

  if (missing_loom && isTRUE(paths$loom_root_explicit)) {
    figure7_stop("Explicit --loom-root does not exist: ", paths$loom_root)
  }
  if (missing_seurat && isTRUE(paths$seurat_rds_explicit)) {
    figure7_stop("Explicit --seurat-rds does not exist: ", paths$seurat_rds)
  }
  raw_download_roles <- c(if (missing_loom) "loom", if (missing_seurat) "seurat_rds")
  if (length(raw_download_roles) && !isTRUE(paths$download_missing_raw)) {
    figure7_stop(
      "Required raw Figure 7 data are missing and automatic download is disabled: ",
      paste(raw_download_roles, collapse = ", "),
      ". Supply valid --loom-root/--seurat-rds or use --download-missing-raw=true."
    )
  }
  if (length(raw_download_roles)) figure7_require_workflow_file(paths$raw_manifest, "raw-data-manifest")
  raw_validation_roles <- c(
    if (needs_loom && !isTRUE(paths$loom_root_explicit)) "loom",
    if (needs_seurat && !isTRUE(paths$seurat_rds_explicit)) "seurat_rds"
  )
  if (length(raw_validation_roles)) figure7_require_workflow_file(paths$raw_manifest, "raw-data-manifest")

  if (needs_loom) {
    if (!nzchar(paths$python) || !file.exists(paths$python)) {
      figure7_stop("Missing required scVelo Python executable: ", if (nzchar(paths$python)) paths$python else "not supplied")
    }
  }

  raw_data_status <- if (length(raw_download_roles)) {
    "planned_zenodo_download"
  } else if (length(raw_validation_roles)) {
    "planned_zenodo_validation"
  } else if (needs_loom || needs_seurat) {
    if (isTRUE(paths$loom_root_explicit) || isTRUE(paths$seurat_rds_explicit)) "explicit_local_inputs" else "reused_zenodo_cache"
  } else if (cell_pair_ready) {
    "not_required_existing_cell_tables"
  } else {
    "not_required_existing_scvelo"
  }

  list(
    state = state,
    cell_pair_ready = cell_pair_ready,
    scvelo_ready = has_scvelo,
    needs_raw_download = length(raw_download_roles) > 0L,
    raw_download_roles = raw_download_roles,
    needs_raw_stage = length(raw_validation_roles) > 0L,
    raw_validation_roles = raw_validation_roles,
    raw_data_status = raw_data_status,
    needs_loom = needs_loom,
    needs_seurat = needs_seurat,
    needs_scvelo = identical(state, "raw_inputs_ready"),
    needs_cell_tables = !cell_pair_ready
  )
}

figure7_stage_cli_args <- function(values) {
  keep <- !vapply(values, function(value) is.null(value) || length(value) == 0L || !nzchar(as.character(value)), logical(1L))
  values <- values[keep]
  paste0("--", names(values), "=", vapply(values, as.character, character(1L)))
}

figure7_run_stage <- function(label, script, values, log_path) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  rscript <- file.path(R.home("bin"), "Rscript")
  cli <- figure7_stage_cli_args(values)
  command_args <- c(shQuote(script), vapply(cli, shQuote, character(1L)))
  message("[full-workflow] ", label)
  status <- system2(rscript, args = command_args, stdout = log_path, stderr = log_path)
  if (!identical(status, 0L)) {
    log_lines <- if (file.exists(log_path)) readLines(log_path, warn = FALSE) else character()
    detail <- if (length(log_lines)) paste(tail(log_lines, 60L), collapse = "\n") else "No stage log was written."
    figure7_stop(label, " failed; log: ", log_path, "\n", detail)
  }
  invisible(log_path)
}

figure7_prepare_full_workflow <- function(
  args,
  paths,
  preflight,
  script_dir,
  config_path,
  config,
  include_panel_f,
  overwrite_intermediates = FALSE
) {
  dir.create(paths$intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$log_dir, recursive = TRUE, showWarnings = FALSE)
  executed <- character()
  raw_data_status <- preflight$raw_data_status

  if (isTRUE(preflight$needs_raw_stage)) {
    figure7_run_stage(
      "Download and verify Figure 7 raw data from Zenodo",
      file.path(script_dir, "download_figure7_raw_data.R"),
      list(
        raw_data_dir = paths$raw_data_dir,
        manifest = paths$raw_manifest,
        roles = paste(preflight$raw_validation_roles, collapse = ","),
        "download-workers" = paths$download_workers,
        "download-connections-per-file" = paths$download_connections_per_file,
        allow_download = if (isTRUE(paths$download_missing_raw)) "TRUE" else "FALSE"
      ),
      file.path(paths$log_dir, "00_raw_data_download.log")
    )
    audit_path <- file.path(paths$raw_data_dir, "provenance", "downloaded_files_checksums.tsv")
    figure7_require_workflow_file(audit_path, "raw-data-download-audit")
    audit <- figure7_read_tsv(audit_path, c("role", "filename", "action", "md5", "sha256"))
    current_audit <- audit[audit$role %in% preflight$raw_validation_roles, , drop = FALSE]
    if (!nrow(current_audit)) figure7_stop("Raw-data download audit does not cover requested roles")
    downloaded <- any(current_audit$action == "downloaded")
    executed <- c(executed, if (downloaded) "raw_data_download" else "raw_data_validation")
    raw_data_status <- if (downloaded) "downloaded_from_zenodo" else "reused_zenodo_cache"
  }
  if (isTRUE(preflight$needs_scvelo)) {
    figure7_require_workflow_dir(paths$loom_root, "loom-root")
    figure7_require_workflow_file(paths$seurat_rds, "seurat-rds")
  } else if (isTRUE(include_panel_f)) {
    figure7_require_workflow_file(paths$seurat_rds, "seurat-rds")
  }
  seurat_rds_sha256 <- ""
  if (isTRUE(preflight$needs_seurat)) {
    seurat_rds_sha256 <- if (isTRUE(paths$seurat_rds_explicit)) {
      figure7_verify_checksum(paths$seurat_rds, config$raw_data$seurat_rds_sha256, "Figure 7 Seurat RDS")
    } else {
      as.character(config$raw_data$seurat_rds_sha256)
    }
  }

  if (isTRUE(preflight$needs_scvelo)) {
    figure7_run_stage(
      "Generate scVelo cell metrics",
      file.path(script_dir, "generate_scvelo_cell_metrics.R"),
      list(
        seurat_rds = paths$seurat_rds,
        loom_root = paths$loom_root,
        output = paths$scvelo_metrics,
        seurat_metadata_output = paths$seurat_metadata,
        python = paths$python,
        work_dir = paths$scvelo_work_dir
      ),
      file.path(paths$log_dir, "01_scvelo_metrics.log")
    )
    executed <- c(executed, "scvelo_metrics")
    figure7_require_workflow_file(paths$scvelo_metrics, "scvelo-metrics")
    figure7_require_workflow_file(paths$seurat_metadata, "seurat-metadata-output")
  }

  if (isTRUE(preflight$needs_cell_tables)) {
    figure7_run_stage(
      "Generate CellCycle and NonCellCycle cell-level inputs",
      file.path(script_dir, "generate_pseudotime_distribution_with_ploidy_dose_tgi.R"),
      list(
        scvelo_metrics_input = paths$scvelo_metrics,
        cell_ploidy_input = paths$cell_ploidy,
        sample_info_input = paths$sample_info,
        growth_curve_input = paths$growth_curve,
        cellcycle_output = paths$cellcycle,
        noncellcycle_output = paths$noncellcycle
      ),
      file.path(paths$log_dir, "02_celllevel_inputs.log")
    )
    executed <- c(executed, "celllevel_inputs")
  }
  figure7_require_workflow_file(paths$cellcycle, "cellcycle-input")
  figure7_require_workflow_file(paths$noncellcycle, "non-cellcycle-input")

  saved_reference <- ""
  if (isTRUE(include_panel_f)) {
    if (dir.exists(paths$state_pathway_root) && !isTRUE(overwrite_intermediates)) {
      figure7_stop(
        "State-pathway output already exists: ", paths$state_pathway_root,
        ". Use a new --output-dir/--state-pathway-output-root or --overwrite-intermediates=true."
      )
    }
    if (dir.exists(paths$saved_reference)) {
      if (!isTRUE(overwrite_intermediates)) {
        figure7_stop(
          "Saved state-pathway reference already exists: ", paths$saved_reference,
          ". Use a new destination or --overwrite-intermediates=true."
        )
      }
      unlink(paths$saved_reference, recursive = TRUE, force = TRUE)
    }
    figure7_run_stage(
      "Generate state-pathway support results",
      file.path(script_dir, "generate_pseudotime_state_pathways_support.R"),
      list(
        cell_metadata = paths$cellcycle,
        noncell_metadata = paths$noncellcycle,
        seurat_rds = paths$seurat_rds,
        config = config_path,
        output_root = paths$state_pathway_root,
        overwrite = if (isTRUE(overwrite_intermediates)) "TRUE" else "FALSE"
      ),
      file.path(paths$log_dir, "03_state_pathway_support.log")
    )
    executed <- c(executed, "state_pathway_support")
    figure7_run_stage(
      "Export compact state-pathway reference",
      file.path(script_dir, "export_state_pathway_reference.R"),
      list(
        results_root = paths$state_pathway_root,
        output_dir = paths$saved_reference,
        verify_source_checksums = "FALSE"
      ),
      file.path(paths$log_dir, "04_state_pathway_export.log")
    )
    executed <- c(executed, "state_pathway_export")
    if (!dir.exists(paths$saved_reference)) {
      figure7_stop("State-pathway exporter did not create the expected reference directory: ", paths$saved_reference)
    }
    saved_reference <- paths$saved_reference
  }

  list(
    initial_state = preflight$state,
    executed_stages = executed,
    cellcycle = paths$cellcycle,
    noncellcycle = paths$noncellcycle,
    scvelo_metrics = paths$scvelo_metrics,
    seurat_metadata = if (file.exists(paths$seurat_metadata)) paths$seurat_metadata else "",
    raw_data_dir = paths$raw_data_dir,
    raw_manifest = paths$raw_manifest,
    raw_data_status = raw_data_status,
    raw_download_roles = preflight$raw_download_roles,
    raw_validation_roles = preflight$raw_validation_roles,
    download_workers = paths$download_workers,
    download_connections_per_file = paths$download_connections_per_file,
    loom_root = paths$loom_root,
    seurat_rds = paths$seurat_rds,
    seurat_rds_sha256 = seurat_rds_sha256,
    state_pathway_root = if (isTRUE(include_panel_f)) paths$state_pathway_root else "",
    saved_reference = saved_reference,
    log_dir = paths$log_dir
  )
}
