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
  reference_id <- as.character(
    config$state_pathways$generated_reference_id
  )
  canonical_cellcycle <- file.path(
    "Data", "in-vivo", "figure7", "processed",
    paste0(
      "CellCycleCells_pseudotime_distribution_per_sample_cell_level_",
      "with_ploidy_dose_tgi.csv"
    )
  )
  canonical_noncellcycle <- file.path(
    "Data", "in-vivo", "figure7", "processed",
    paste0(
      "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_",
      "with_ploidy_dose_tgi.csv"
    )
  )
  cellcycle_arg <- figure7_arg(args, "cellcycle-input", "")
  noncellcycle_arg <- figure7_arg(args, "non-cellcycle-input", "")
  generated_cellcycle <- file.path(
    intermediate_dir,
    basename(canonical_cellcycle)
  )
  generated_noncellcycle <- file.path(
    intermediate_dir,
    basename(canonical_noncellcycle)
  )
  raw_data_dir <- figure7_workflow_path(
    figure7_arg(args, "raw-data-dir", as.character(config$raw_data$default_root)),
    repo_root
  )
  loom_root_arg <- figure7_arg(args, "loom-root", "")
  seurat_rds_arg <- figure7_arg(args, "seurat-rds", "")
  seurat_upstream_dir_arg <- figure7_arg(
    args,
    "seurat-upstream-dir",
    file.path(intermediate_dir, "seurat_upstream")
  )
  cellranger_root_arg <- figure7_arg(args, "cellranger-root", "")
  deposited_seurat_rds <- file.path(
    raw_data_dir,
    as.character(config$raw_data$seurat_filename)
  )
  state_pathway_root <- figure7_workflow_path(
    figure7_arg(
      args,
      "state-pathway-output-root",
      file.path(intermediate_dir, "state_pathways", reference_id)
    ),
    repo_root
  )
  saved_reference <- figure7_workflow_path(
    figure7_arg(
      args,
      "saved-state-pathway-dir",
      file.path(intermediate_dir, "saved_state_pathway", reference_id)
    ),
    repo_root
  )
  list(
    intermediate_dir = intermediate_dir,
    raw_data_dir = raw_data_dir,
    raw_manifest = figure7_workflow_path(as.character(config$raw_data$required_manifest), repo_root),
    download_missing_raw = figure7_flag(args, "download-missing-raw", TRUE),
    download_workers = figure7_workflow_integer_arg(args, "download-workers", 4L),
    download_connections_per_file = figure7_workflow_integer_arg(
      args, "download-connections-per-file", 2L
    ),
    analysis_jobs = figure7_workflow_integer_arg(
      args, "jobs", 16L, maximum = 64L
    ),
    loom_root_explicit = nzchar(loom_root_arg),
    seurat_rds_explicit = nzchar(seurat_rds_arg),
    explicit_seurat_rds = figure7_workflow_path(
      seurat_rds_arg,
      repo_root
    ),
    deposited_seurat_rds = figure7_workflow_path(
      deposited_seurat_rds,
      repo_root
    ),
    seurat_upstream_dir = figure7_workflow_path(
      seurat_upstream_dir_arg,
      repo_root
    ),
    cellranger_root = figure7_workflow_path(
      cellranger_root_arg,
      repo_root
    ),
    cellranger_root_explicit = nzchar(cellranger_root_arg),
    cell_pair_explicit = nzchar(cellcycle_arg) || nzchar(noncellcycle_arg),
    cellcycle_explicit = nzchar(cellcycle_arg),
    noncellcycle_explicit = nzchar(noncellcycle_arg),
    frozen_cellcycle = figure7_workflow_path(
      figure7_arg(
        args,
        "frozen-cellcycle-input",
        canonical_cellcycle
      ),
      repo_root
    ),
    frozen_noncellcycle = figure7_workflow_path(
      figure7_arg(
        args,
        "frozen-non-cellcycle-input",
        canonical_noncellcycle
      ),
      repo_root
    ),
    frozen_cell_pair_explicit =
      !is.null(args[["frozen_cellcycle_input"]]) ||
      !is.null(args[["frozen_non_cellcycle_input"]]),
    scvelo_bundle_explicit = !is.null(args[["scvelo_metrics"]]),
    scvelo_metrics = figure7_workflow_path(
      figure7_arg(args, "scvelo-metrics", file.path(intermediate_dir, "scvelo_cell_metrics.csv")),
      repo_root
    ),
    scvelo_stage_manifest = file.path(
      intermediate_dir,
      "scvelo_stage_manifest.tsv"
    ),
    cellcycle = figure7_workflow_path(
      if (nzchar(cellcycle_arg)) cellcycle_arg else generated_cellcycle,
      repo_root
    ),
    noncellcycle = figure7_workflow_path(
      if (nzchar(noncellcycle_arg)) noncellcycle_arg else generated_noncellcycle,
      repo_root
    ),
    generated_cellcycle = figure7_workflow_path(
      generated_cellcycle,
      repo_root
    ),
    generated_noncellcycle = figure7_workflow_path(
      generated_noncellcycle,
      repo_root
    ),
    cell_table_provenance = file.path(
      intermediate_dir,
      "cell_table_provenance.tsv"
    ),
    loom_root = figure7_workflow_path(
      if (nzchar(loom_root_arg)) loom_root_arg else file.path(raw_data_dir, as.character(config$raw_data$loom_subdir)),
      repo_root
    ),
    seurat_rds = figure7_workflow_path(
      if (nzchar(seurat_rds_arg)) seurat_rds_arg else deposited_seurat_rds,
      repo_root
    ),
    # The default deposited location is only a fallback candidate.  It must
    # not be treated as the active lineage input while validating a deeper
    # cache that may have been authored from a reconstructed Seurat object.
    seurat_rds_lineage_active = nzchar(seurat_rds_arg),
    cell_ploidy = figure7_workflow_path(
      figure7_arg(
        args, "cell-ploidy-input",
        as.character(
          config$versioned_source_artifacts$endpoint_ploidy$default_path
        )
      ),
      repo_root
    ),
    sample_info = figure7_workflow_path(
      figure7_arg(
        args, "sample-info-input",
        as.character(
          config$versioned_source_artifacts$sample_info$default_path
        )
      ),
      repo_root
    ),
    growth_curve = figure7_workflow_path(
      figure7_arg(
        args, "growth-curve-input",
        as.character(
          config$versioned_source_artifacts$growth_curve$default_path
        )
      ),
      repo_root
    ),
    python = figure7_arg(args, "python", Sys.which("python")),
    scvelo_work_dir = file.path(intermediate_dir, "scvelo_work"),
    state_pathway_root = state_pathway_root,
    state_pathway_root_explicit =
      !is.null(args[["state_pathway_output_root"]]),
    saved_reference = saved_reference,
    saved_reference_explicit =
      !is.null(args[["saved_state_pathway_dir"]]),
    frozen_reference = figure7_workflow_path(
      figure7_arg(
        args,
        "frozen-state-pathway-dir",
        file.path(
          as.character(config$state_pathways$reference_root),
          as.character(config$state_pathways$reference_id)
        )
      ),
      repo_root
    ),
    frozen_reference_explicit =
      !is.null(args[["frozen_state_pathway_dir"]]),
    state_stage_manifest = file.path(
      state_pathway_root,
      "00_manifest",
      "figure7_stage_manifest.tsv"
    ),
    reference_stage_manifest = file.path(
      dirname(saved_reference),
      paste0(basename(saved_reference), ".stage_manifest.tsv")
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

figure7_quarantine_existing <- function(path, label) {
  if (!file.exists(path) && !dir.exists(path)) return("")
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
  quarantine <- paste0(path, ".stale.", stamp, ".", Sys.getpid())
  if (!file.rename(path, quarantine)) {
    figure7_stop(
      "Cannot preserve stale ",
      label,
      " before regeneration: ",
      path
    )
  }
  message("[full-workflow] Preserved stale ", label, " at ", quarantine)
  quarantine
}

figure7_path_is_within <- function(path, root) {
  normalized_path <- normalizePath(path, mustWork = FALSE)
  normalized_root <- normalizePath(root, mustWork = FALSE)
  identical(normalized_path, normalized_root) ||
    startsWith(normalized_path, paste0(normalized_root, .Platform$file.sep))
}

figure7_quarantine_generated_cache <- function(
  path,
  label,
  intermediate_dir,
  explicit
) {
  if (isTRUE(explicit) ||
      !figure7_path_is_within(path, intermediate_dir)) {
    figure7_stop(
      "Explicit/external ",
      label,
      " failed validation and will not be modified: ",
      path,
      ". Choose a fresh path under --intermediate-dir."
    )
  }
  figure7_quarantine_existing(path, label)
}

figure7_state_result_files <- function(root) {
  file.path(
    root,
    c(
      "00_manifest/frozen_interval_definition.csv",
      "00_manifest/feature_species_audit.csv",
      "00_manifest/gene_set_contract.csv",
      "00_manifest/gene_set_membership.csv",
      "00_manifest/input_checksums.csv",
      "00_manifest/package_versions.csv",
      "binning/00_manifest/analysis_parameters.csv",
      "binning/01_qc/primary_coverage_check.csv",
      "binning/02_pseudobulk/sample_bin_metadata.csv",
      "binning/ETP_reference_balanced_threshold_2_24/00_manifest/model_parameters.csv",
      "binning/ETP_reference_balanced_threshold_2_24/01_qc/model_design_rank_audit.csv",
      "binning/ETP_reference_balanced_threshold_2_24/03_gene_models/gene_primary_adjacent_state_contrast.csv",
      "binning/ETP_reference_balanced_threshold_2_24/03_gene_models/gene_symbol_resolution.csv",
      "binning/ETP_reference_balanced_threshold_2_24/04_gsea/all_collections_primary_adjacent_state_gsea.csv",
      "binning/ETP_reference_balanced_threshold_2_24/04_gsea/all_collections_leading_edge_genes.csv",
      "binning/ETP_reference_balanced_threshold_2_24/04_gsea/pathway_activity_over_pseudotime.csv"
    )
  )
}

figure7_state_results_complete <- function(root) {
  dir.exists(root) && all(file.exists(figure7_state_result_files(root)))
}

figure7_saved_reference_complete <- function(path) {
  dir.exists(path) &&
    all(file.exists(file.path(path, figure7_generated_state_required_files())))
}

figure7_saved_reference_match_inputs <- function(
  path,
  paths,
  config,
  config_path
) {
  if (!figure7_saved_reference_complete(path)) return(FALSE)
  reference <- tryCatch(
    figure7_validate_generated_state_reference(
      path,
      config,
      expected_inputs = list(
        cellcycle = paths$cellcycle,
        noncellcycle = paths$noncellcycle,
        seurat_rds = figure7_active_lineage_seurat(paths)
      ),
      config_path = config_path
    ),
    error = function(error) NULL
  )
  if (is.null(reference)) return(FALSE)
  paths$saved_reference <- path
  isTRUE(tryCatch(
    figure7_reference_stage_manifest_matches(
      paths,
      config_path,
      config
    ),
    error = function(error) FALSE
  ))
}

figure7_state_results_match_inputs <- function(paths, config_path, config) {
  if (!figure7_state_results_complete(paths$state_pathway_root)) return(FALSE)
  expected <- figure7_state_dependency_values(paths, config_path, config)
  figure7_stage_manifest_matches(
    paths$state_stage_manifest,
    expected = expected,
    outputs = figure7_state_output_hashes(paths),
    required_keys = c(
      "support_script_sha256",
      "feature_species_policy_code_sha256",
      "r_environment_contract_sha256",
      "cellcycle_sha256",
      "noncellcycle_sha256",
      "seurat_rds_sha256",
      "feature_species_policy",
      "gene_set_provider",
      "gene_set_package_version",
      "gene_set_release",
      "gene_set_collections",
      "analysis_seed",
      "config_contract_sha256",
      "parameter_contract_sha256",
      "package_versions_sha256",
      "analysis_parameters_sha256",
      "model_parameters_sha256"
    )
  )
}

figure7_read_key_value_file <- function(path) {
  if (!file.exists(path)) return(character())
  table <- tryCatch(
    figure7_read_tsv(path, c("key", "value")),
    error = function(error) NULL
  )
  if (is.null(table) || anyDuplicated(table$key)) return(character())
  stats::setNames(as.character(table$value), table$key)
}

figure7_sha256_text <- function(value) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    figure7_stop("R package 'digest' is required")
  }
  unname(digest::digest(
    paste(as.character(value), collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  ))
}

figure7_stage_manifest_values <- function(values) {
  value_names <- names(values)
  values <- as.character(values)
  names(values) <- value_names
  values <- values[order(names(values))]
  fingerprint_values <- values[!startsWith(names(values), "audit_")]
  c(
    values,
    stage_fingerprint = figure7_sha256_text(
      paste(
        names(fingerprint_values),
        fingerprint_values,
        sep = "=",
        collapse = "\n"
      )
    )
  )
}

figure7_write_stage_manifest <- function(values, path) {
  values <- figure7_stage_manifest_values(values)
  figure7_write_tsv(
    data.frame(
      key = names(values),
      value = unname(values),
      stringsAsFactors = FALSE
    ),
    path
  )
  invisible(path)
}

figure7_stage_manifest_matches <- function(
  path,
  expected = character(),
  outputs = character(),
  required_keys = character()
) {
  observed <- figure7_read_key_value_file(path)
  if (!length(observed) ||
      !"stage_fingerprint" %in% names(observed) ||
      any(!required_keys %in% names(observed))) {
    return(FALSE)
  }
  body <- observed[names(observed) != "stage_fingerprint"]
  body <- body[order(names(body))]
  fingerprint_body <- body[!startsWith(names(body), "audit_")]
  expected_fingerprint <- figure7_sha256_text(
    paste(
      names(fingerprint_body),
      fingerprint_body,
      sep = "=",
      collapse = "\n"
    )
  )
  if (!identical(
    unname(observed[["stage_fingerprint"]]),
    expected_fingerprint
  )) {
    return(FALSE)
  }
  if (length(expected) &&
      {
        contract_expected <- expected[
          !startsWith(names(expected), "audit_")
        ]
        !all(names(contract_expected) %in% names(observed)) ||
       !identical(
         unname(observed[names(contract_expected)]),
         unname(as.character(contract_expected))
       )
      }) {
    return(FALSE)
  }
  if (length(outputs)) {
    if (is.null(names(outputs)) || any(!nzchar(names(outputs))) ||
        any(!file.exists(outputs)) ||
        !all(names(outputs) %in% names(observed))) {
      return(FALSE)
    }
    output_hashes <- vapply(outputs, figure7_sha256, character(1L))
    if (!identical(
      unname(observed[names(outputs)]),
      unname(output_hashes)
    )) {
      return(FALSE)
    }
  }
  TRUE
}

figure7_script_hash <- function(config_path, filename) {
  figure7_sha256(file.path(dirname(config_path), filename))
}

figure7_scvelo_config_contract <- function(config) {
  si <- config$si_figures
  figure7_contract_sha256(c(
    umap_reduction = as.character(si$umap_reduction),
    pca_reduction = as.character(si$pca_reduction),
    cluster_id_field = as.character(si$cluster_id_field),
    ploidy_field = as.character(si$ploidy_field),
    context_field = as.character(si$context_field)
  ))
}

figure7_scvelo_parameter_contract <- function(jobs = 16L) {
  paste(
    "mode=stochastic",
    "min_shared_counts=20",
    "n_top_genes=2000",
    "n_pcs=30",
    "n_neighbors=30",
    paste0("n_jobs=", as.integer(jobs)),
    "cell_match_contract=all_seurat_cells_exactly_once",
    "root_clusters=6",
    "end_clusters=",
    "use_metadata_pca=true",
    "pca_prefix=PCA_",
    sep = ";"
  )
}

figure7_cell_table_parameter_contract <- function(config) {
  paste(
    paste0("configured_day=", figure7_tgi_day(config)),
    "matched_control=mean_by_initial_ploidy",
    "cellcycle_clusters=4c,6,10",
    "cellcycle_mapping=CellCycle/NonCellCycle",
    sep = ";"
  )
}

figure7_active_lineage_seurat <- function(paths) {
  active <- if (is.null(paths$seurat_rds_lineage_active)) {
    # Direct callers historically pass a selected lineage path without the
    # explicit planner flag.
    TRUE
  } else {
    isTRUE(paths$seurat_rds_lineage_active)
  }
  path <- as.character(paths$seurat_rds)
  if (!active || length(path) != 1L || !nzchar(path) ||
      !file.exists(path)) {
    return("")
  }
  normalizePath(path, mustWork = TRUE)
}

figure7_cell_table_dependency_values <- function(paths, config, config_path) {
  environment_lock <- figure7_environment_lock_path(config_path, config)
  values <- c(
    schema_version = "1",
    artifact = "figure7_cell_level_intermediate_pair",
    generator_script_sha256 = figure7_script_hash(
      config_path,
      "generate_pseudotime_distribution_with_ploidy_dose_tgi.R"
    ),
    r_environment_contract_sha256 =
      figure7_environment_stage_contract_sha256(
        environment_lock,
        "r",
        "celllevel"
      ),
    audit_environment_lock_sha256 = figure7_sha256(environment_lock),
    parameter_contract = figure7_cell_table_parameter_contract(config)
  )
  dependencies <- c(
    scvelo_metrics_sha256 = paths$scvelo_metrics,
    cell_ploidy_sha256 = paths$cell_ploidy,
    sample_info_sha256 = paths$sample_info,
    growth_curve_sha256 = paths$growth_curve
  )
  dependencies <- dependencies[file.exists(dependencies)]
  c(
    values,
    vapply(dependencies, figure7_sha256, character(1L))
  )
}

figure7_r_package_versions <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing)) {
    figure7_stop(
      "Missing required R package(s) before raw-data access: ",
      paste(missing, collapse = ", ")
    )
  }
  paste(
    paste(
      packages,
      vapply(
        packages,
        function(package) as.character(utils::packageVersion(package)),
        character(1L)
      ),
      sep = "="
    ),
    collapse = ";"
  )
}

figure7_python_package_versions <- function(python) {
  if (!nzchar(python) || !file.exists(python)) {
    figure7_stop(
      "Missing required scVelo Python executable before raw-data access: ",
      if (nzchar(python)) python else "not supplied"
    )
  }
  packages <- c("numpy", "pandas", "anndata", "scanpy", "scvelo")
  code <- paste0(
    "import importlib, importlib.metadata\n",
    "names=[",
    paste(sprintf("'%s'", packages), collapse = ","),
    "]\n",
    "out=[]\n",
    "for name in names:\n",
    " importlib.import_module(name)\n",
    " out.append(name+'='+importlib.metadata.version(name))\n",
    "print(';'.join(out))"
  )
  output <- suppressWarnings(system2(
    python,
    c("-c", shQuote(code)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    figure7_stop(
      "Python dependency preflight failed before raw-data access: ",
      paste(output, collapse = "\n")
    )
  }
  value <- tail(output[nzchar(output)], 1L)
  if (length(value) != 1L ||
      !all(vapply(packages, function(package) {
        grepl(paste0("(^|;)", package, "="), value)
      }, logical(1L)))) {
    figure7_stop(
      "Python dependency preflight did not report numpy/pandas/anndata/",
      "scanpy/scvelo versions"
    )
  }
  value
}

figure7_require_download_runtime <- function(environment_lock) {
  figure7_validate_r_environment(environment_lock, "download")
  figure7_r_package_versions("digest")
  downloaders <- c(
    Sys.which("aria2c"),
    Sys.which("wget"),
    Sys.which("curl")
  )
  if (!any(nzchar(downloaders)) && !capabilities("libcurl")) {
    figure7_stop(
      "A resumable downloader (aria2c, wget, curl, or R libcurl) is ",
      "required before raw-data access"
    )
  }
  invisible(TRUE)
}

figure7_require_workflow_dependencies <- function(
  needs_scvelo,
  needs_cell_tables,
  needs_state,
  needs_export,
  needs_upstream,
  needs_download,
  python,
  config_path,
  config
) {
  needs_computation <- isTRUE(needs_scvelo) ||
    isTRUE(needs_cell_tables) ||
    isTRUE(needs_state) ||
    isTRUE(needs_export) ||
    isTRUE(needs_upstream)
  needs_lock <- needs_computation || isTRUE(needs_download)
  environment_lock <- if (needs_lock) {
    figure7_environment_lock_path(config_path, config)
  } else {
    ""
  }
  if (needs_lock) {
    figure7_read_environment_lock(environment_lock)
  }
  if (isTRUE(needs_download)) {
    figure7_require_download_runtime(environment_lock)
  }
  if (isTRUE(needs_scvelo)) {
    figure7_r_package_versions(c("digest", "yaml", "Seurat", "readr"))
    figure7_validate_r_environment(environment_lock, "scvelo")
    figure7_validate_python_environment(
      environment_lock,
      python,
      stage = "scvelo"
    )
  }
  if (isTRUE(needs_cell_tables)) {
    figure7_r_package_versions(c("digest", "readr", "readxl"))
    figure7_validate_r_environment(environment_lock, "celllevel")
  }
  if (isTRUE(needs_state)) {
    figure7_r_package_versions(c(
      "Seurat", "SeuratObject", "Matrix", "yaml", "digest", "edgeR", "limma",
      "fgsea", "msigdbr", "readr"
    ))
    figure7_validate_r_environment(environment_lock, "state")
  }
  if (isTRUE(needs_export)) {
    figure7_r_package_versions(c("readr", "digest"))
    figure7_validate_r_environment(environment_lock, "state_export")
  }
  if (isTRUE(needs_upstream)) {
    figure7_validate_r_environment(
      environment_lock,
      "seurat_upstream"
    )
  }
  invisible(TRUE)
}

figure7_loom_manifest_table <- function(paths) {
  manifest <- utils::read.delim(
    paths$raw_manifest,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )
  required <- c("role", "filename", "size_bytes", "md5")
  if (!all(required %in% names(manifest))) {
    figure7_stop("Raw Figure 7 manifest lacks role/filename/size_bytes/md5")
  }
  loom <- manifest[manifest$role == "loom", , drop = FALSE]
  loom$size_bytes <- suppressWarnings(as.numeric(loom$size_bytes))
  if (nrow(loom) != 18L ||
      anyDuplicated(loom$filename) ||
      any(!is.finite(loom$size_bytes) | loom$size_bytes <= 0) ||
      any(!grepl("^[0-9a-f]{32}$", loom$md5))) {
    figure7_stop("Raw Figure 7 manifest must pin exactly 18 loom files")
  }
  loom
}

figure7_loom_manifest_dependencies <- function(paths) {
  loom <- figure7_loom_manifest_table(paths)
  stats::setNames(
    as.character(loom$md5),
    paste0("loom_md5:", loom$filename)
  )
}

figure7_observed_loom_sha_dependencies <- function(paths) {
  loom_manifest <- figure7_loom_manifest_table(paths)
  manifest_names <- sub(
    "^loom_md5:",
    "",
    names(figure7_loom_manifest_dependencies(paths))
  )
  audit_path <- file.path(
    paths$raw_data_dir,
    "provenance",
    "downloaded_files_checksums.tsv"
  )
  if (!isTRUE(paths$loom_root_explicit) && file.exists(audit_path)) {
    audit <- tryCatch(
      utils::read.delim(
        audit_path,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(error) NULL
    )
    if (!is.null(audit) &&
        all(c("role", "filename", "sha256") %in% names(audit))) {
      audit <- audit[
        audit$role == "loom" & audit$filename %in% manifest_names,
        ,
        drop = FALSE
      ]
      audit <- audit[match(manifest_names, audit$filename), , drop = FALSE]
      if (nrow(audit) == length(manifest_names) &&
          !anyNA(audit$filename) &&
          all(grepl("^[0-9a-f]{64}$", audit$sha256))) {
        return(stats::setNames(
          as.character(audit$sha256),
          paste0("loom_sha256:", audit$filename)
        ))
      }
    }
  }
  candidates <- list.files(
    paths$loom_root,
    pattern = "[.]loom$",
    recursive = TRUE,
    full.names = TRUE
  )
  matched <- vapply(
    manifest_names,
    function(filename) {
      hit <- candidates[basename(candidates) == filename]
      if (length(hit) != 1L) {
        figure7_stop(
          "Expected exactly one deposited loom named ",
          filename,
          " below ",
          paths$loom_root
        )
      }
      hit
    },
    character(1L)
  )
  matched <- matched[match(loom_manifest$filename, basename(matched))]
  observed_sizes <- unname(file.info(matched)$size)
  observed_md5 <- unname(tools::md5sum(matched))
  if (!identical(
        as.numeric(observed_sizes),
        as.numeric(loom_manifest$size_bytes)
      ) ||
      !identical(
        as.character(observed_md5),
        as.character(loom_manifest$md5)
      )) {
    figure7_stop(
      "Explicit/local loom inputs differ from the pinned Zenodo ",
      "size/MD5 contract"
    )
  }
  stats::setNames(
    vapply(matched, figure7_sha256, character(1L)),
    paste0("loom_sha256:", manifest_names)
  )
}

figure7_scvelo_dependency_values <- function(paths, config_path, config) {
  environment_lock <- figure7_environment_lock_path(config_path, config)
  seurat_rds <- figure7_active_lineage_seurat(paths)
  c(
    schema_version = "1",
    artifact = "figure7_scvelo_bundle",
    if (nzchar(seurat_rds)) c(
      seurat_rds_sha256 = figure7_sha256(seurat_rds)
    ),
    audit_config_sha256 = figure7_sha256(config_path),
    config_contract_sha256 = figure7_scvelo_config_contract(config),
    raw_manifest_sha256 = figure7_sha256(paths$raw_manifest),
    generator_script_sha256 = figure7_script_hash(
      config_path,
      "generate_scvelo_cell_metrics.R"
    ),
    r_environment_contract_sha256 =
      figure7_environment_stage_contract_sha256(
        environment_lock,
        "r",
        "scvelo"
      ),
    python_environment_contract_sha256 =
      figure7_environment_stage_contract_sha256(
        environment_lock,
        "python",
        "scvelo"
      ),
    audit_environment_lock_sha256 = figure7_sha256(environment_lock),
    parameter_contract = figure7_scvelo_parameter_contract(paths$analysis_jobs),
    figure7_loom_manifest_dependencies(paths)
  )
}

figure7_state_output_hashes <- function(paths) {
  files <- figure7_state_result_files(paths$state_pathway_root)
  stats::setNames(
    files,
    paste0(
      "output_sha256:",
      substring(
        files,
        nchar(paths$state_pathway_root) + 2L
      )
    )
  )
}

figure7_reference_output_hashes <- function(paths) {
  files <- file.path(
    paths$saved_reference,
    figure7_generated_state_required_files()
  )
  stats::setNames(
    files,
    paste0("output_sha256:", basename(files))
  )
}

figure7_environment_lock_path <- function(config_path, config) {
  configured <- as.character(config$raw_data$environment_lock)
  if (length(configured) != 1L || !nzchar(configured)) {
    figure7_stop("Figure 7 config must declare raw_data.environment_lock")
  }
  path <- if (grepl("^/", configured)) {
    configured
  } else {
    repo_root <- normalizePath(
      file.path(dirname(config_path), "..", "..", ".."),
      mustWork = FALSE
    )
    file.path(repo_root, configured)
  }
  figure7_require_workflow_file(path, "environment-lock")
  normalizePath(path, mustWork = TRUE)
}

figure7_reference_dependency_values <- function(
  paths,
  config_path,
  config
) {
  state_manifest_available <- file.exists(paths$state_stage_manifest)
  state_manifest <- character()
  if (state_manifest_available) {
    if (!figure7_stage_manifest_matches(paths$state_stage_manifest)) {
      figure7_stop(
        "Available state-stage manifest fails its self-authenticating ",
        "fingerprint"
      )
    }
    state_manifest <- figure7_read_key_value_file(
      paths$state_stage_manifest
    )
  }
  environment_lock <- figure7_environment_lock_path(config_path, config)
  state_dependencies <- figure7_state_dependency_values(
    paths,
    config_path,
    config
  )
  c(
    schema_version = "1",
    artifact = "figure7_generated_state_reference",
    generated_reference_id =
      as.character(config$state_pathways$generated_reference_id),
    state_dependencies[intersect(
      c(
        "cellcycle_sha256",
        "noncellcycle_sha256",
        "seurat_rds_sha256"
      ),
      names(state_dependencies)
    )],
    audit_figure7_config_sha256 = figure7_sha256(config_path),
    config_contract_sha256 =
      state_dependencies[["config_contract_sha256"]],
    support_script_sha256 =
      state_dependencies[["support_script_sha256"]],
    feature_species_policy_code_sha256 =
      state_dependencies[["feature_species_policy_code_sha256"]],
    exporter_script_sha256 = figure7_script_hash(
      config_path,
      "export_state_pathway_reference.R"
    ),
    exporter_common_io_sha256 = figure7_sha256(file.path(
      dirname(config_path),
      "src",
      "common_io.R"
    )),
    r_environment_contract_sha256 =
      figure7_environment_stage_contract_sha256(
        environment_lock,
        "r",
        "state_export"
      ),
    audit_environment_lock_sha256 = figure7_sha256(environment_lock),
    if (state_manifest_available) c(
      audit_state_stage_manifest_sha256 =
        figure7_sha256(paths$state_stage_manifest),
      state_stage_fingerprint =
        unname(state_manifest[["stage_fingerprint"]])
    ),
    gene_set_release = as.character(config$gene_sets$database_release),
    analysis_seed = as.character(config$statistics$seed)
  )
}

figure7_reference_stage_manifest_matches <- function(
  paths,
  config_path,
  config
) {
  expected <- figure7_reference_dependency_values(
    paths,
    config_path,
    config
  )
  outputs <- figure7_reference_output_hashes(paths)
  figure7_stage_manifest_matches(
    paths$reference_stage_manifest,
    expected = expected,
    outputs = outputs,
    required_keys = c(
      "generated_reference_id",
      "cellcycle_sha256",
      "noncellcycle_sha256",
      "seurat_rds_sha256",
      "audit_figure7_config_sha256",
      "config_contract_sha256",
      "support_script_sha256",
      "feature_species_policy_code_sha256",
      "exporter_script_sha256",
      "exporter_common_io_sha256",
      "r_environment_contract_sha256",
      "audit_state_stage_manifest_sha256",
      "state_stage_fingerprint",
      "gene_set_release",
      "analysis_seed",
      names(outputs)
    )
  )
}

figure7_state_dependency_values <- function(paths, config_path, config) {
  environment_lock <- figure7_environment_lock_path(config_path, config)
  seurat_rds <- figure7_active_lineage_seurat(paths)
  c(
    schema_version = "1",
    artifact = "figure7_state_pathway_results",
    if (file.exists(paths$cellcycle)) c(
      cellcycle_sha256 = figure7_sha256(paths$cellcycle)
    ),
    if (file.exists(paths$noncellcycle)) c(
      noncellcycle_sha256 = figure7_sha256(paths$noncellcycle)
    ),
    if (nzchar(seurat_rds)) c(
      seurat_rds_sha256 = figure7_sha256(seurat_rds)
    ),
    audit_config_sha256 = figure7_sha256(config_path),
    support_script_sha256 = figure7_script_hash(
      config_path,
      "generate_pseudotime_state_pathways_support.R"
    ),
    feature_species_policy_code_sha256 = figure7_sha256(file.path(
      dirname(config_path),
      "src",
      "feature_species_policy.R"
    )),
    r_environment_contract_sha256 =
      figure7_environment_stage_contract_sha256(
        environment_lock,
        "r",
        "state"
      ),
    audit_environment_lock_sha256 = figure7_sha256(environment_lock),
    feature_species_policy = as.character(
      config$state_pathways$generated_feature_species_policy
    ),
    gene_set_provider = as.character(config$gene_sets$provider),
    gene_set_package_version = as.character(
      config$gene_sets$package_version
    ),
    gene_set_release = as.character(config$gene_sets$database_release),
    gene_set_collections = paste(
      as.character(unlist(config$state_pathways$collections)),
      collapse = ","
    ),
    analysis_seed = as.character(config$statistics$seed),
    model = as.character(config$state_pathways$model),
    config_contract_sha256 =
      figure7_state_config_contract_sha256(config),
    parameter_contract_sha256 = figure7_sha256_text(c(
      unname(figure7_feature_species_contract_values(config)),
      config$state_pathways$assay,
      config$state_pathways$counts_layer,
      config$state_pathways$pseudotime_bins,
      config$state_pathways$minimum_cells_per_sample_bin,
      config$state_pathways$expression_filter,
      config$state_pathways$normalization,
      config$state_pathways$observation_model,
      config$state_pathways$mouse_block,
      config$state_pathways$spline_df,
      config$state_pathways$grid_size,
      config$state_pathways$contrast,
      config$state_pathways$gsea_rank_statistic
    ))
  )
}

figure7_cell_pair_valid <- function(paths, config, config_path) {
  if (!file.exists(paths$cellcycle) || !file.exists(paths$noncellcycle)) {
    return(FALSE)
  }
  observed_hashes <- c(
    cellcycle_sha256 = figure7_sha256(paths$cellcycle),
    noncellcycle_sha256 = figure7_sha256(paths$noncellcycle)
  )
  canonical_hashes <- c(
    cellcycle_sha256 = as.character(config$inputs$cellcycle_sha256),
    noncellcycle_sha256 = as.character(config$inputs$noncellcycle_sha256)
  )
  # The reviewed A-E inputs are a byte contract, not a checkout-location
  # contract.  Exact copies remain reusable after relocation and do not need a
  # generated-stage manifest.
  if (identical(unname(observed_hashes), unname(canonical_hashes))) {
    return(TRUE)
  }
  expected <- figure7_cell_table_dependency_values(
    paths,
    config,
    config_path
  )
  figure7_stage_manifest_matches(
    paths$cell_table_provenance,
    expected = expected,
    outputs = c(
      cellcycle_sha256 = paths$cellcycle,
      noncellcycle_sha256 = paths$noncellcycle
    ),
    required_keys = c(
      "generator_script_sha256",
      "r_environment_contract_sha256",
      "parameter_contract",
      "scvelo_metrics_sha256",
      "cell_ploidy_sha256",
      "sample_info_sha256",
      "growth_curve_sha256"
    )
  )
}

figure7_scvelo_bundle_valid <- function(paths, config_path, config) {
  if (!file.exists(paths$scvelo_metrics)) return(FALSE)
  stage_expected <- figure7_scvelo_dependency_values(
    paths,
    config_path,
    config
  )
  loom_keys <- sub(
    "^loom_md5:",
    "loom_sha256:",
    names(figure7_loom_manifest_dependencies(paths))
  )
  if (isTRUE(paths$loom_root_explicit)) {
    stage_expected <- c(
      stage_expected,
      figure7_observed_loom_sha_dependencies(paths)
    )
  }
  figure7_stage_manifest_matches(
    paths$scvelo_stage_manifest,
    expected = stage_expected,
    outputs = c(
      scvelo_metrics_sha256 = paths$scvelo_metrics
    ),
    required_keys = c(
      "seurat_rds_sha256",
      "r_environment_contract_sha256",
      "python_environment_contract_sha256",
      "python_executable",
      "python_package_versions",
      "r_package_versions",
      loom_keys
    )
  )
}

figure7_write_scvelo_stage_manifest <- function(
  paths,
  config_path,
  config
) {
  values <- c(
    figure7_scvelo_dependency_values(paths, config_path, config),
    figure7_observed_loom_sha_dependencies(paths),
    python_executable = basename(paths$python),
    python_package_versions =
      figure7_python_package_versions(paths$python),
    r_package_versions = figure7_r_package_versions(
      c("digest", "yaml", "Seurat", "readr")
    ),
    figure7_r_runtime_provenance(),
    scvelo_metrics_sha256 = figure7_sha256(paths$scvelo_metrics)
  )
  figure7_write_stage_manifest(values, paths$scvelo_stage_manifest)
  if (!figure7_scvelo_bundle_valid(paths, config_path, config)) {
    figure7_stop(
      "New scVelo metrics failed its complete stage fingerprint"
    )
  }
  invisible(paths$scvelo_stage_manifest)
}

figure7_write_cell_pair_stage_manifest <- function(
  paths,
  config,
  config_path
) {
  values <- c(
    figure7_cell_table_dependency_values(paths, config, config_path),
    figure7_r_runtime_provenance(),
    cellcycle_sha256 = figure7_sha256(paths$cellcycle),
    noncellcycle_sha256 = figure7_sha256(paths$noncellcycle)
  )
  figure7_write_stage_manifest(values, paths$cell_table_provenance)
  if (!figure7_cell_pair_valid(paths, config, config_path)) {
    figure7_stop(
      "New CellCycle/NonCellCycle pair failed its complete stage fingerprint"
    )
  }
  invisible(paths$cell_table_provenance)
}

figure7_write_state_stage_manifest <- function(
  paths,
  config_path,
  config
) {
  package_versions <- file.path(
    paths$state_pathway_root,
    "00_manifest",
    "package_versions.csv"
  )
  analysis_parameters <- file.path(
    paths$state_pathway_root,
    "binning",
    "00_manifest",
    "analysis_parameters.csv"
  )
  model_parameters <- file.path(
    paths$state_pathway_root,
    "binning",
    as.character(config$state_pathways$model),
    "00_manifest",
    "model_parameters.csv"
  )
  values <- c(
    figure7_state_dependency_values(paths, config_path, config),
    figure7_r_runtime_provenance(),
    package_versions_sha256 = figure7_sha256(package_versions),
    analysis_parameters_sha256 = figure7_sha256(analysis_parameters),
    model_parameters_sha256 = figure7_sha256(model_parameters),
    vapply(
      figure7_state_output_hashes(paths),
      figure7_sha256,
      character(1L)
    )
  )
  figure7_write_stage_manifest(values, paths$state_stage_manifest)
  if (!figure7_state_results_match_inputs(paths, config_path, config)) {
    figure7_stop(
      "New state-pathway result tree failed its complete stage fingerprint"
    )
  }
  invisible(paths$state_stage_manifest)
}

figure7_write_reference_stage_manifest <- function(
  paths,
  config_path,
  config
) {
  values <- c(
    figure7_reference_dependency_values(paths, config_path, config),
    vapply(
      figure7_reference_output_hashes(paths),
      figure7_sha256,
      character(1L)
    )
  )
  figure7_write_stage_manifest(values, paths$reference_stage_manifest)
  if (!figure7_saved_reference_match_inputs(
    paths$saved_reference,
    paths,
    config,
    config_path
  )) {
    figure7_stop(
      "New generated state reference failed its complete stage fingerprint"
    )
  }
  invisible(paths$reference_stage_manifest)
}

figure7_paths_with_cell_pair <- function(paths, cellcycle, noncellcycle) {
  paths$cellcycle <- cellcycle
  paths$noncellcycle <- noncellcycle
  paths
}

figure7_cell_pair_has_canonical_hashes <- function(paths, config) {
  file.exists(paths$cellcycle) &&
    file.exists(paths$noncellcycle) &&
    identical(
      c(
        figure7_sha256(paths$cellcycle),
        figure7_sha256(paths$noncellcycle)
      ),
      c(
        as.character(config$inputs$cellcycle_sha256),
        as.character(config$inputs$noncellcycle_sha256)
      )
    )
}

figure7_frozen_reference_valid <- function(path, config) {
  !is.null(tryCatch(
    figure7_validate_state_reference(
      path,
      config,
      verify_checksums = TRUE
    ),
    error = function(error) NULL
  ))
}

figure7_select_cell_pair <- function(
  paths,
  config_path,
  config,
  overwrite_intermediates = FALSE
) {
  if (xor(paths$cellcycle_explicit, paths$noncellcycle_explicit)) {
    figure7_stop(
      "--cellcycle-input and --non-cellcycle-input must be supplied together"
    )
  }

  if (isTRUE(paths$cell_pair_explicit)) {
    selected <- figure7_paths_with_cell_pair(
      paths,
      paths$cellcycle,
      paths$noncellcycle
    )
    if (!figure7_cell_pair_valid(selected, config, config_path)) {
      figure7_stop(
        "Explicit Figure 7 A-E input pair failed the reviewed-byte or ",
        "generated-stage contract; explicit inputs are read-only"
      )
    }
    return(list(
      paths = selected,
      ready = TRUE,
      source = if (figure7_cell_pair_has_canonical_hashes(selected, config)) {
        "reviewed_frozen"
      } else {
        "generated_cache"
      }
    ))
  }

  frozen <- figure7_paths_with_cell_pair(
    paths,
    paths$frozen_cellcycle,
    paths$frozen_noncellcycle
  )
  if (figure7_cell_pair_has_canonical_hashes(frozen, config)) {
    return(list(
      paths = frozen,
      ready = TRUE,
      source = "reviewed_frozen"
    ))
  }
  if (isTRUE(paths$frozen_cell_pair_explicit)) {
    figure7_stop(
      "Explicit frozen Figure 7 A-E pair failed its byte-pinned contract; ",
      "the files will not be modified"
    )
  }

  generated <- figure7_paths_with_cell_pair(
    paths,
    paths$generated_cellcycle,
    paths$generated_noncellcycle
  )
  generated_valid <- !isTRUE(overwrite_intermediates) &&
    figure7_cell_pair_valid(generated, config, config_path)
  list(
    paths = generated,
    ready = generated_valid,
    source = if (generated_valid) "generated_cache" else "raw_rebuild"
  )
}

figure7_preflight_workflow <- function(
  paths,
  config_path,
  config,
  include_panel_f = TRUE,
  overwrite_intermediates = FALSE
) {
  selection <- figure7_select_cell_pair(
    paths,
    config_path,
    config,
    overwrite_intermediates
  )
  paths <- selection$paths
  cell_pair_ready <- isTRUE(selection$ready)

  historical_reference_ready <- isTRUE(include_panel_f) &&
    figure7_frozen_reference_valid(paths$frozen_reference, config)
  if (isTRUE(include_panel_f) &&
      isTRUE(paths$frozen_reference_explicit)) {
    figure7_stop(
      "Full-workflow panel 7F cannot use an explicit historical mixed ",
      "reference; omit --frozen-state-pathway-dir so a generated ",
      "GRCh-only v2 reference is built or reused"
    )
  }
  # Full-workflow is the corrected recomputation path. The byte-pinned v1
  # reference remains renderable for historical audit, but it must never
  # short-circuit GRCh-only model generation.
  frozen_reference_ready <- FALSE
  scvelo_bundle_paths <- c(
    paths$scvelo_metrics,
    paths$scvelo_stage_manifest
  )
  scvelo_bundle_exists <- file.exists(scvelo_bundle_paths)
  has_scvelo <- !cell_pair_ready &&
    !isTRUE(overwrite_intermediates) &&
    figure7_scvelo_bundle_valid(paths, config_path, config)

  generated_cell_exists <- file.exists(c(
    paths$generated_cellcycle,
    paths$generated_noncellcycle
  ))
  if (!cell_pair_ready && any(generated_cell_exists)) {
    message(
      "[full-workflow] Existing generated cell-table cache is partial/stale; ",
      "dependent stages will be regenerated"
    )
  }

  if (!cell_pair_ready && any(scvelo_bundle_exists) &&
      !all(scvelo_bundle_exists)) {
    message(
      "[full-workflow] Existing scVelo cache is partial; it will be ",
      "regenerated"
    )
  }
  if (!cell_pair_ready && all(scvelo_bundle_exists) && !has_scvelo &&
      !isTRUE(overwrite_intermediates)) {
    message(
      "[full-workflow] Existing scVelo cache is stale; it will be regenerated"
    )
  }

  ae_state <- if (cell_pair_ready) {
    "cell_tables_ready"
  } else if (has_scvelo) {
    "scvelo_ready"
  } else {
      "raw_inputs_ready"
  }

  saved_reference_ready <- isTRUE(include_panel_f) &&
    !frozen_reference_ready &&
    !isTRUE(overwrite_intermediates) &&
    figure7_saved_reference_match_inputs(
      paths$saved_reference,
      paths,
      config,
      config_path
    )
  state_results_ready <- isTRUE(include_panel_f) &&
    !frozen_reference_ready &&
    !saved_reference_ready &&
    !isTRUE(overwrite_intermediates) &&
    figure7_state_results_match_inputs(paths, config_path, config)
  needs_scvelo <- !cell_pair_ready && !has_scvelo
  needs_cell_tables <- !cell_pair_ready
  needs_state <- isTRUE(include_panel_f) &&
    !frozen_reference_ready &&
    !saved_reference_ready &&
    !state_results_ready
  needs_export <- isTRUE(include_panel_f) &&
    !frozen_reference_ready &&
    !saved_reference_ready
  needs_loom <- needs_scvelo
  needs_seurat <- needs_scvelo || needs_state
  environment_lock <- if (needs_seurat) {
    figure7_environment_lock_path(config_path, config)
  } else {
    ""
  }
  seurat_selection <- if (needs_seurat) {
    figure7_select_seurat_source(
      needs_seurat = TRUE,
      explicit_rds = paths$explicit_seurat_rds,
      deposited_rds = paths$deposited_seurat_rds,
      deposited_sha256 = as.character(config$raw_data$seurat_rds_sha256),
      upstream_dir = paths$seurat_upstream_dir,
      cellranger_root = paths$cellranger_root,
      module_dir = dirname(config_path),
      environment_lock = environment_lock,
      config = config,
      all_ploidy = paths$cell_ploidy,
      sample_info = paths$sample_info
    )
  } else {
    list(
      required = FALSE,
      action = "not_required",
      source = "not_required_existing_validated_caches",
      rds = "",
      rds_sha256 = "",
      final_stage_manifest = "",
      final_stage_manifest_sha256 = "",
      reconstruction_manifest = "",
      reconstruction_manifest_sha256 = "",
      transitive_dependencies = character(),
      scientific_code_contracts = character(),
      upstream_common_contract = character(),
      resumable_stage = "",
      resumable_stage_manifest = "",
      resumable_stage_manifest_sha256 = "",
      upstream_dir = paths$seurat_upstream_dir,
      cellranger_root = paths$cellranger_root
    )
  }
  if (isTRUE(seurat_selection$required)) {
    paths$seurat_rds <- seurat_selection$rds
    paths$seurat_rds_lineage_active <- TRUE
  }
  needs_upstream <- needs_seurat &&
    identical(seurat_selection$action, "build_upstream")

  if (!cell_pair_ready) {
    figure7_require_workflow_file(paths$cell_ploidy, "cell-ploidy-input")
    figure7_require_workflow_file(paths$sample_info, "sample-info-input")
    figure7_require_workflow_file(paths$growth_curve, "growth-curve-input")
    figure7_verify_checksum(
      paths$cell_ploidy,
      config$versioned_source_artifacts$endpoint_ploidy$sha256,
      "Figure 7 endpoint-ploidy input"
    )
    figure7_verify_checksum(
      paths$sample_info,
      config$versioned_source_artifacts$sample_info$sha256,
      "Figure 7 sample-info input"
    )
    figure7_verify_checksum(
      paths$growth_curve,
      config$versioned_source_artifacts$growth_curve$sha256,
      "Figure 7 growth-curve input"
    )
  } else if (needs_upstream) {
    figure7_require_workflow_file(paths$cell_ploidy, "cell-ploidy-input")
    figure7_require_workflow_file(paths$sample_info, "sample-info-input")
    figure7_verify_checksum(
      paths$cell_ploidy,
      config$versioned_source_artifacts$endpoint_ploidy$sha256,
      "Figure 7 endpoint-ploidy input"
    )
    figure7_verify_checksum(
      paths$sample_info,
      config$versioned_source_artifacts$sample_info$sha256,
      "Figure 7 sample-info input"
    )
  }
  missing_loom <- needs_loom && !dir.exists(paths$loom_root)
  missing_deposited_seurat <- needs_seurat &&
    identical(seurat_selection$action, "download_deposited")

  if (missing_loom && isTRUE(paths$loom_root_explicit)) {
    figure7_stop("Explicit --loom-root does not exist: ", paths$loom_root)
  }
  raw_download_roles <- c(
    if (missing_loom) "loom",
    if (missing_deposited_seurat) "seurat_rds"
  )
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
    if (needs_seurat &&
        seurat_selection$action %in%
          c("validate_deposited", "download_deposited")) {
      "seurat_rds"
    }
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
  } else if (needs_upstream) {
    "planned_cellranger_reconstruction"
  } else if (needs_loom || needs_seurat) {
    seurat_selection$source
  } else if (cell_pair_ready) {
    "not_required_existing_validated_caches"
  } else {
    "not_required_existing_scvelo"
  }

  figure7_require_workflow_dependencies(
    needs_scvelo = needs_scvelo,
    needs_cell_tables = needs_cell_tables,
    needs_state = needs_state,
    needs_export = needs_export,
    needs_upstream = needs_upstream,
    needs_download = length(raw_download_roles) > 0L,
    python = paths$python,
    config_path = config_path,
    config = config
  )

  f_state <- if (!isTRUE(include_panel_f)) {
    ""
  } else if (frozen_reference_ready) {
    "frozen_reference_ready"
  } else if (saved_reference_ready) {
    "generated_reference_ready"
  } else if (state_results_ready) {
    "state_results_ready"
  } else {
    "state_raw_inputs_ready"
  }

  list(
    state = if (isTRUE(include_panel_f)) {
      paste(ae_state, f_state, sep = "__")
    } else {
      ae_state
    },
    ae_state = ae_state,
    f_state = f_state,
    paths = paths,
    cell_pair_source = selection$source,
    cell_pair_ready = cell_pair_ready,
    scvelo_ready = has_scvelo,
    needs_raw_download = length(raw_download_roles) > 0L,
    raw_download_roles = raw_download_roles,
    needs_raw_stage = length(raw_validation_roles) > 0L,
    raw_validation_roles = raw_validation_roles,
    raw_data_status = raw_data_status,
    needs_loom = needs_loom,
    needs_seurat = needs_seurat,
    needs_upstream = needs_upstream,
    seurat_selection = seurat_selection,
    needs_scvelo = needs_scvelo,
    needs_cell_tables = needs_cell_tables,
    frozen_reference_ready = frozen_reference_ready,
    historical_reference_ready = historical_reference_ready,
    saved_reference_ready = saved_reference_ready,
    state_results_ready = state_results_ready,
    needs_state = needs_state,
    needs_export = needs_export
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
  paths <- preflight$paths
  dir.create(paths$intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$log_dir, recursive = TRUE, showWarnings = FALSE)
  executed <- character()
  raw_data_status <- preflight$raw_data_status

  if (isTRUE(preflight$needs_raw_stage)) {
    figure7_run_stage(
      "Download and verify Figure 7 raw data from Zenodo",
      file.path(script_dir, "download_figure7_raw_data.R"),
      list(
        "raw-data-dir" = paths$raw_data_dir,
        manifest = paths$raw_manifest,
        roles = paste(preflight$raw_validation_roles, collapse = ","),
        "download-workers" = paths$download_workers,
        "download-connections-per-file" = paths$download_connections_per_file,
        "allow-download" = if (isTRUE(paths$download_missing_raw)) "TRUE" else "FALSE"
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

  seurat_selection <- preflight$seurat_selection
  if (isTRUE(preflight$needs_upstream)) {
    figure7_run_stage(
      if (nzchar(seurat_selection$resumable_stage)) {
        paste0(
          "Resume Seurat reconstruction from validated ",
          seurat_selection$resumable_stage,
          " cache"
        )
      } else {
        "Reconstruct final Seurat object from Cell Ranger matrices"
      },
      file.path(
        script_dir,
        "generate_final_seurat_from_cellranger.R"
      ),
      figure7_seurat_upstream_generator_args(
        seurat_selection,
        figure7_environment_lock_path(config_path, config),
        config_path,
        paths$cell_ploidy,
        paths$sample_info,
        paths$analysis_jobs
      ),
      file.path(paths$log_dir, "00b_seurat_upstream.log")
    )
    validation <- figure7_validate_seurat_upstream_artifact(
      output_root = seurat_selection$upstream_dir,
      module_dir = script_dir,
      environment_lock =
        figure7_environment_lock_path(config_path, config),
      config = config,
      all_ploidy = paths$cell_ploidy,
      sample_info = paths$sample_info,
      expected_rds = seurat_selection$rds,
      cellranger_root = if (
        dir.exists(seurat_selection$cellranger_root)
      ) {
        seurat_selection$cellranger_root
      } else {
        ""
      }
    )
    paths$seurat_rds <- validation$rds
    paths$seurat_rds_lineage_active <- TRUE
    seurat_selection$action <- "reuse"
    seurat_selection$source <- if (
      nzchar(seurat_selection$resumable_stage)
    ) {
      paste0(
        "resumed_reconstructed_rds_from_",
        seurat_selection$resumable_stage
      )
    } else {
      "reconstructed_rds_from_cellranger"
    }
    seurat_selection$rds <- validation$rds
    seurat_selection$rds_sha256 <- validation$rds_sha256
    seurat_selection$final_stage_manifest <-
      validation$final_stage_manifest
    seurat_selection$final_stage_manifest_sha256 <-
      validation$final_stage_manifest_sha256
    seurat_selection$reconstruction_manifest <-
      validation$reconstruction_manifest
    seurat_selection$reconstruction_manifest_sha256 <-
      validation$reconstruction_manifest_sha256
    seurat_selection$transitive_dependencies <-
      validation$dependencies
    seurat_selection$scientific_code_contracts <-
      validation$scientific_code_contracts
    seurat_selection$upstream_common_contract <-
      validation$upstream_common_contract
    executed <- c(executed, "seurat_upstream")
    raw_data_status <- seurat_selection$source
  }

  if (isTRUE(preflight$needs_scvelo)) {
    figure7_require_workflow_dir(paths$loom_root, "loom-root")
    figure7_require_workflow_file(paths$seurat_rds, "seurat-rds")
  } else if (isTRUE(preflight$needs_seurat)) {
    figure7_require_workflow_file(paths$seurat_rds, "seurat-rds")
  }
  seurat_rds_sha256 <- ""
  if (isTRUE(preflight$needs_seurat)) {
    seurat_rds_sha256 <- if (
      nzchar(seurat_selection$reconstruction_manifest)
    ) {
      observed <- figure7_sha256(paths$seurat_rds)
      if (!identical(observed, seurat_selection$rds_sha256)) {
        figure7_stop(
          "Validated reconstructed Seurat RDS changed before computation"
        )
      }
      observed
    } else {
      figure7_verify_checksum(
        paths$seurat_rds,
        config$raw_data$seurat_rds_sha256,
        "Figure 7 deposited Seurat RDS"
      )
    }
  }

  if (isTRUE(preflight$needs_scvelo)) {
    for (stale in c(
      paths$scvelo_metrics,
      paths$scvelo_stage_manifest,
      paths$scvelo_work_dir
    )) {
      if (file.exists(stale) || dir.exists(stale)) {
        figure7_quarantine_generated_cache(
          stale,
          "generated Figure 7 scVelo cache",
          paths$intermediate_dir,
          paths$scvelo_bundle_explicit
        )
      }
    }
    figure7_run_stage(
      "Generate scVelo cell metrics",
      file.path(script_dir, "generate_scvelo_cell_metrics.R"),
      list(
        seurat_rds = paths$seurat_rds,
        loom_root = paths$loom_root,
        output = paths$scvelo_metrics,
        config = config_path,
        python = paths$python,
        n_jobs = paths$analysis_jobs,
        work_dir = paths$scvelo_work_dir
      ),
      file.path(paths$log_dir, "01_scvelo_metrics.log")
    )
    executed <- c(executed, "scvelo_metrics")
    figure7_require_workflow_file(paths$scvelo_metrics, "scvelo-metrics")
    figure7_write_scvelo_stage_manifest(paths, config_path, config)
  }

  if (isTRUE(preflight$needs_cell_tables)) {
    for (stale in c(
      paths$cellcycle,
      paths$noncellcycle,
      paths$cell_table_provenance
    )) {
      if (file.exists(stale)) {
        figure7_quarantine_generated_cache(
          stale,
          "generated Figure 7 cell-table cache",
          paths$intermediate_dir,
          paths$cell_pair_explicit
        )
      }
    }
    figure7_run_stage(
      "Generate CellCycle and NonCellCycle cell-level inputs",
      file.path(script_dir, "generate_pseudotime_distribution_with_ploidy_dose_tgi.R"),
      list(
        scvelo_metrics_input = paths$scvelo_metrics,
        cell_ploidy_input = paths$cell_ploidy,
        sample_info_input = paths$sample_info,
        growth_curve_input = paths$growth_curve,
        cellcycle_output = paths$cellcycle,
        noncellcycle_output = paths$noncellcycle,
        tgi_day = figure7_tgi_day(config)
      ),
      file.path(paths$log_dir, "02_celllevel_inputs.log")
    )
    executed <- c(executed, "celllevel_inputs")
    figure7_write_cell_pair_stage_manifest(paths, config, config_path)
  }
  figure7_require_workflow_file(paths$cellcycle, "cellcycle-input")
  figure7_require_workflow_file(paths$noncellcycle, "non-cellcycle-input")

  saved_reference <- ""
  reference_id <- ""
  reference_kind <- ""
  canonical_publication_allowed <- FALSE
  if (isTRUE(include_panel_f)) {
    if (isTRUE(preflight$saved_reference_ready)) {
      figure7_validate_generated_state_reference(
        paths$saved_reference,
        config,
        expected_inputs = list(
          cellcycle = paths$cellcycle,
          noncellcycle = paths$noncellcycle,
          seurat_rds = paths$seurat_rds
        ),
        config_path = config_path
      )
      saved_reference <- paths$saved_reference
      reference_id <- as.character(
        config$state_pathways$generated_reference_id
      )
      reference_kind <- as.character(
        config$state_pathways$generated_reference_kind
      )
    } else {
      if (dir.exists(paths$saved_reference)) {
        figure7_quarantine_generated_cache(
          paths$saved_reference,
          "generated state-pathway reference",
          paths$intermediate_dir,
          paths$saved_reference_explicit
        )
      }

      state_results_ready <- isTRUE(preflight$state_results_ready)
      if (!state_results_ready) {
        if (dir.exists(paths$state_pathway_root)) {
          figure7_quarantine_generated_cache(
            paths$state_pathway_root,
            "state-pathway result tree",
            paths$intermediate_dir,
            paths$state_pathway_root_explicit
          )
        }
        figure7_run_stage(
          "Generate state-pathway support results",
          file.path(
            script_dir,
            "generate_pseudotime_state_pathways_support.R"
          ),
          list(
            cell_metadata = paths$cellcycle,
            noncell_metadata = paths$noncellcycle,
            seurat_rds = paths$seurat_rds,
            config = config_path,
            output_root = paths$state_pathway_root,
            overwrite = "FALSE"
          ),
          file.path(paths$log_dir, "03_state_pathway_support.log")
        )
        executed <- c(executed, "state_pathway_support")
        figure7_write_state_stage_manifest(paths, config_path, config)
      }

      figure7_run_stage(
        "Export compact state-pathway reference",
        file.path(script_dir, "export_state_pathway_reference.R"),
        list(
          results_root = paths$state_pathway_root,
          output_dir = paths$saved_reference,
          config = config_path
        ),
        file.path(paths$log_dir, "04_state_pathway_export.log")
      )
      executed <- c(executed, "state_pathway_export")
      figure7_write_reference_stage_manifest(paths, config_path, config)
      figure7_validate_generated_state_reference(
        paths$saved_reference,
        config,
        expected_inputs = list(
          cellcycle = paths$cellcycle,
          noncellcycle = paths$noncellcycle,
          seurat_rds = paths$seurat_rds
        ),
        config_path = config_path
      )
      saved_reference <- paths$saved_reference
      reference_id <- as.character(
        config$state_pathways$generated_reference_id
      )
      reference_kind <- as.character(
        config$state_pathways$generated_reference_kind
      )
    }
  }

  list(
    initial_state = preflight$state,
    executed_stages = executed,
    cellcycle = paths$cellcycle,
    noncellcycle = paths$noncellcycle,
    scvelo_metrics = paths$scvelo_metrics,
    raw_data_dir = paths$raw_data_dir,
    raw_manifest = paths$raw_manifest,
    raw_data_status = raw_data_status,
    raw_download_roles = preflight$raw_download_roles,
    raw_validation_roles = preflight$raw_validation_roles,
    download_workers = paths$download_workers,
    download_connections_per_file = paths$download_connections_per_file,
    analysis_jobs = paths$analysis_jobs,
    loom_root = paths$loom_root,
    seurat_rds = paths$seurat_rds,
    seurat_rds_sha256 = seurat_rds_sha256,
    seurat_source = seurat_selection$source,
    seurat_reconstruction_manifest =
      seurat_selection$reconstruction_manifest,
    seurat_reconstruction_manifest_sha256 =
      seurat_selection$reconstruction_manifest_sha256,
    seurat_transitive_dependencies =
      seurat_selection$transitive_dependencies,
    seurat_final_stage_manifest =
      seurat_selection$final_stage_manifest,
    seurat_final_stage_manifest_sha256 =
      seurat_selection$final_stage_manifest_sha256,
    seurat_scientific_code_contracts =
      seurat_selection$scientific_code_contracts,
    seurat_upstream_common_contract =
      seurat_selection$upstream_common_contract,
    state_pathway_root = if (isTRUE(include_panel_f)) {
      paths$state_pathway_root
    } else {
      ""
    },
    saved_reference = saved_reference,
    state_pathway_reference_id = reference_id,
    state_pathway_reference_kind = reference_kind,
    canonical_publication_allowed = canonical_publication_allowed,
    log_dir = paths$log_dir
  )
}
