if (!exists("figure7_preflight_workflow", mode = "function")) {
  source(
    file.path(module_dir, "src", "generated_state_pathway_reference.R"),
    local = FALSE
  )
  source(file.path(module_dir, "src", "input_preflight.R"), local = FALSE)
}
if (!exists("figure7_select_seurat_source", mode = "function")) {
  source(
    file.path(module_dir, "src", "seurat_upstream_selection.R"),
    local = FALSE
  )
}

testthat::test_that("full-workflow reuses the exact canonical A-E pair without raw data", {
  output_dir <- tempfile("figure7_cached_ae_output_")
  intermediate_dir <- tempfile("figure7_cached_ae_intermediates_")
  cellcycle <- file.path(
    repo_root,
    "Data/in-vivo/figure7/processed",
    "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  )
  noncellcycle <- file.path(
    repo_root,
    "Data/in-vivo/figure7/processed",
    "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  )
  result <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      file.path(module_dir, "run_figure7.R"),
      "--mode=full-workflow",
      "--panel-set=a-e",
      "--preflight-only=true",
      "--download-missing-raw=false",
      paste0("--config=", file.path(module_dir, "figure7_config.yaml")),
      paste0("--cellcycle-input=", cellcycle),
      paste0("--non-cellcycle-input=", noncellcycle),
      paste0("--intermediate-dir=", intermediate_dir),
      paste0("--output-dir=", output_dir)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  testthat::expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
  testthat::expect_true(any(grepl(
    "workflow_state\\tcell_tables_ready",
    result
  )))
  testthat::expect_true(any(grepl(
    "raw_data_status\\tnot_required_existing_validated_caches",
    result
  )))
  testthat::expect_false(dir.exists(output_dir))
})

testthat::test_that("workflow lineage metadata uses explicit portable sentinels", {
  deposited <- list(
    seurat_source = "deposited_rds",
    seurat_reconstruction_manifest = "",
    seurat_reconstruction_manifest_sha256 = "",
    seurat_final_stage_manifest = "",
    seurat_final_stage_manifest_sha256 = ""
  )
  testthat::expect_identical(
    figure7_workflow_scalar(
      deposited,
      "seurat_reconstruction_manifest_sha256",
      "not_available"
    ),
    "not_available"
  )
  testthat::expect_identical(
    figure7_workflow_scalar(
      deposited,
      "seurat_final_stage_manifest_sha256",
      "not_available"
    ),
    "not_available"
  )

  reconstructed <- list(
    seurat_source = "reused_reconstructed_rds",
    seurat_reconstruction_manifest = "seurat_upstream/reconstruction_manifest.tsv",
    seurat_reconstruction_manifest_sha256 = paste(rep("a", 64L), collapse = ""),
    seurat_final_stage_manifest = "seurat_upstream/final/manifest.tsv",
    seurat_final_stage_manifest_sha256 = paste(rep("b", 64L), collapse = "")
  )
  for (name in names(reconstructed)) {
    testthat::expect_identical(
      figure7_workflow_scalar(reconstructed, name, "not_available"),
      reconstructed[[name]]
    )
  }
})

testthat::test_that("full-workflow persists no-lineage metadata without blank hashes", {
  root <- tempfile("figure7_cached_ae_metadata_")
  output_dir <- file.path(root, "output")
  intermediate_dir <- file.path(root, "intermediates")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  cellcycle <- file.path(
    repo_root,
    "Data/in-vivo/figure7/processed",
    "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  )
  noncellcycle <- file.path(
    repo_root,
    "Data/in-vivo/figure7/processed",
    "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  )
  result <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      file.path(module_dir, "run_figure7.R"),
      "--mode=full-workflow",
      "--panel-set=a-e",
      "--download-missing-raw=false",
      paste0("--config=", file.path(module_dir, "figure7_config.yaml")),
      paste0("--cellcycle-input=", cellcycle),
      paste0("--non-cellcycle-input=", noncellcycle),
      paste0("--intermediate-dir=", intermediate_dir),
      paste0("--output-dir=", output_dir)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  testthat::expect_null(
    attr(result, "status"),
    info = paste(result, collapse = "\n")
  )
  run_config <- figure7_read_tsv(
    file.path(output_dir, "metadata", "run_config.tsv"),
    c("key", "value")
  )
  testthat::expect_identical(anyDuplicated(run_config$key), 0L)
  values <- stats::setNames(run_config$value, run_config$key)
  testthat::expect_identical(
    values[["workflow_seurat_source"]],
    "not_required_existing_validated_caches"
  )
  testthat::expect_identical(
    values[["workflow_seurat_reconstruction_manifest"]],
    "not_applicable"
  )
  testthat::expect_identical(
    values[["workflow_seurat_reconstruction_manifest_sha256"]],
    "not_available"
  )
  testthat::expect_identical(
    values[["workflow_seurat_final_stage_manifest"]],
    "not_applicable"
  )
  testthat::expect_identical(
    values[["workflow_seurat_final_stage_manifest_sha256"]],
    "not_available"
  )
})

testthat::test_that("missing raw inputs fail before creating output when download is disabled", {
  output_dir <- tempfile("figure7_no_raw_output_")
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_read_config(config_path)
  paths <- figure7_workflow_paths(
    list(
      intermediate_dir = tempfile("figure7_no_raw_intermediates_"),
      raw_data_dir = tempfile("figure7_no_raw_cache_"),
      download_missing_raw = "false",
      python = Sys.which("python3")
    ),
    repo_root,
    output_dir,
    config
  )
  # Simulate a relocated checkout with neither reviewed frozen branch nor a
  # generated cache. Explicit invalid inputs are tested separately because
  # those are intentionally read-only and must not trigger regeneration.
  paths$frozen_cellcycle <- tempfile("figure7_missing_cellcycle_")
  paths$frozen_noncellcycle <- tempfile("figure7_missing_noncellcycle_")
  paths$frozen_reference <- tempfile("figure7_missing_state_reference_")
  testthat::expect_error(
    figure7_preflight_workflow(
      paths,
      config_path,
      config,
      include_panel_f = TRUE
    ),
    "automatic download is disabled: loom, seurat_rds"
  )
  testthat::expect_false(dir.exists(output_dir))
})

testthat::test_that("generated cell-table stage is reusable and detects tampering", {
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_read_config(config_path)
  root <- tempfile("figure7_cell_stage_")
  dir.create(root)
  paths <- list(
    cellcycle = file.path(root, "CellCycleCells.csv"),
    noncellcycle = file.path(root, "NonCellCycleCells.csv"),
    scvelo_metrics = file.path(root, "scvelo_cell_metrics.csv"),
    cell_ploidy = file.path(
      repo_root,
      config$versioned_source_artifacts$endpoint_ploidy$default_path
    ),
    sample_info = file.path(
      repo_root,
      config$versioned_source_artifacts$sample_info$default_path
    ),
    growth_curve = file.path(
      repo_root,
      config$versioned_source_artifacts$growth_curve$default_path
    ),
    cell_table_provenance = file.path(root, "cell_table_provenance.tsv")
  )
  testthat::expect_true(file.copy(
    file.path(
      repo_root,
      "Data/in-vivo/figure7/processed",
      "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    ),
    paths$cellcycle
  ))
  testthat::expect_true(file.copy(
    file.path(
      repo_root,
      "Data/in-vivo/figure7/processed",
      "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    ),
    paths$noncellcycle
  ))
  writeLines("synthetic upstream metric", paths$scvelo_metrics)

  figure7_write_cell_pair_stage_manifest(paths, config, config_path)
  testthat::expect_true(figure7_cell_pair_valid(paths, config, config_path))
  manifest <- figure7_read_tsv(
    paths$cell_table_provenance,
    c("key", "value")
  )
  testthat::expect_false(any(grepl("^/", manifest$value)))

  write("tamper", paths$cellcycle, append = TRUE)
  testthat::expect_false(figure7_cell_pair_valid(paths, config, config_path))
})

testthat::test_that("explicit state cache manifests follow their selected roots", {
  config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
  root <- tempfile("figure7_explicit_cache_roots_")
  state_root <- file.path(root, "custom-state")
  reference_root <- file.path(root, "custom-reference")
  paths <- figure7_workflow_paths(
    list(
      intermediate_dir = file.path(root, "intermediates"),
      state_pathway_output_root = state_root,
      saved_state_pathway_dir = reference_root
    ),
    repo_root,
    file.path(root, "run"),
    config
  )
  testthat::expect_identical(
    paths$state_stage_manifest,
    file.path(
      normalizePath(state_root, mustWork = FALSE),
      "00_manifest",
      "figure7_stage_manifest.tsv"
    )
  )
  testthat::expect_identical(
    paths$reference_stage_manifest,
    paste0(normalizePath(reference_root, mustWork = FALSE), ".stage_manifest.tsv")
  )
})

testthat::test_that("Zenodo manifest pins the complete raw Figure 7 deposit", {
  manifest <- utils::read.delim(
    file.path(module_dir, "zenodo_required_files.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  testthat::expect_equal(nrow(manifest), 19L)
  testthat::expect_equal(sum(manifest$role == "loom"), 18L)
  testthat::expect_equal(sum(manifest$role == "seurat_rds"), 1L)
  testthat::expect_equal(
    sum(as.numeric(manifest$size_bytes)),
    11395115098
  )
  testthat::expect_true(all(grepl("^[0-9a-f]{32}$", manifest$md5)))
  testthat::expect_true(all(basename(manifest$filename) == manifest$filename))
})

testthat::test_that("raw-stage caches bind stage-scoped environment pins", {
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_read_config(config_path)
  lock_path <- file.path(module_dir, "environment_lock.tsv")
  lock <- figure7_read_environment_lock(lock_path)
  testthat::expect_equal(nrow(lock), 71L)
  testthat::expect_identical(
    figure7_normalize_r_version("1.7-3"),
    "1.7.3"
  )
  testthat::expect_identical(
    figure7_normalize_r_version("1.7.3"),
    "1.7.3"
  )
  testthat::expect_equal(
    lock$version[
      lock$ecosystem == "python" &
        lock$stage == "scvelo" &
        lock$package == "scvelo"
    ],
    "0.3.4"
  )
  testthat::expect_identical(
    stats::setNames(
      lock$version[
        lock$ecosystem == "python" &
          lock$stage == "scvelo" &
          lock$package %in% c(
            "pynndescent", "llvmlite", "joblib", "threadpoolctl",
            "statsmodels", "networkx", "numpy-groupies", "matplotlib"
          )
      ],
      lock$package[
        lock$ecosystem == "python" &
          lock$stage == "scvelo" &
          lock$package %in% c(
            "pynndescent", "llvmlite", "joblib", "threadpoolctl",
            "statsmodels", "networkx", "numpy-groupies", "matplotlib"
          )
      ]
    )[c(
      "pynndescent", "llvmlite", "joblib", "threadpoolctl",
      "statsmodels", "networkx", "numpy-groupies", "matplotlib"
    )],
    c(
      pynndescent = "0.6.0",
      llvmlite = "0.43.0",
      joblib = "1.5.3",
      threadpoolctl = "3.6.0",
      statsmodels = "0.14.6",
      networkx = "3.4.2",
      `numpy-groupies` = "0.11.3",
      matplotlib = "3.8.4"
    )
  )
  testthat::expect_equal(
    lock$version[
      lock$ecosystem == "r" &
        lock$stage == "seurat_upstream" &
        lock$package == "R"
    ],
    "4.5.0"
  )
  testthat::expect_error(
    figure7_validate_r_environment(
      lock_path,
      "celllevel",
      "yaml"
    ),
    "does not explicitly pin R package(s): yaml",
    fixed = TRUE
  )
  testthat::expect_identical(
    unname(figure7_environment_stage_contract(
      lock_path,
      "r",
      "si_raw"
    )[["r/si_raw/future.apply"]]),
    "1.20.2"
  )
  runtime <- figure7_r_runtime_provenance()
  testthat::expect_true(all(nzchar(runtime)))
  python_pins <- lock[
    lock$ecosystem == "python" & lock$stage == "scvelo",
    ,
    drop = FALSE
  ]
  fake_python <- file.path(tempdir(), "figure7_fake_locked_python")
  writeLines(
    c(
      "#!/bin/sh",
      paste0(
        "printf '%s\\n' ",
        paste(
          shQuote(paste0(python_pins$package, "=", python_pins$version)),
          collapse = " "
        )
      )
    ),
    fake_python,
    useBytes = TRUE
  )
  Sys.chmod(fake_python, "0755")
  testthat::expect_silent(figure7_validate_python_environment(
    lock_path,
    fake_python,
    stage = "scvelo"
  ))
  bad_python <- sub(
    "scvelo=0.3.4",
    "scvelo=0.0.0",
    readLines(fake_python, warn = FALSE),
    fixed = TRUE
  )
  writeLines(bad_python, fake_python, useBytes = TRUE)
  testthat::expect_error(
    figure7_validate_python_environment(
      lock_path,
      fake_python,
      stage = "scvelo"
    ),
    "does not match the exact raw-computation lock"
  )

  root <- tempfile("figure7_lock_fingerprints_")
  dir.create(root)
  paths <- list(
    seurat_rds = file.path(root, "seurat.rds"),
    raw_manifest = file.path(module_dir, "zenodo_required_files.tsv"),
    analysis_jobs = 1L,
    scvelo_metrics = file.path(root, "scvelo.csv"),
    cell_ploidy = file.path(root, "ploidy.tsv"),
    sample_info = file.path(root, "sample.xlsx"),
    growth_curve = file.path(root, "growth.xlsx"),
    cellcycle = file.path(root, "cellcycle.csv"),
    noncellcycle = file.path(root, "noncellcycle.csv")
  )
  for (path in paths[c(
    "seurat_rds", "scvelo_metrics", "cell_ploidy", "sample_info",
    "growth_curve", "cellcycle", "noncellcycle"
  )]) {
    writeLines(basename(path), path, useBytes = TRUE)
  }
  expected_lock_hash <- figure7_sha256(lock_path)
  scvelo_dependencies <- figure7_scvelo_dependency_values(
    paths,
    config_path,
    config
  )
  celllevel_dependencies <- figure7_cell_table_dependency_values(
    paths,
    config,
    config_path
  )
  state_dependencies <- figure7_state_dependency_values(
    paths,
    config_path,
    config
  )
  testthat::expect_identical(
    unname(scvelo_dependencies[["audit_environment_lock_sha256"]]),
    expected_lock_hash
  )
  testthat::expect_identical(
    unname(celllevel_dependencies[["audit_environment_lock_sha256"]]),
    expected_lock_hash
  )
  testthat::expect_identical(
    unname(state_dependencies[["audit_environment_lock_sha256"]]),
    expected_lock_hash
  )
  testthat::expect_identical(
    unname(scvelo_dependencies[["r_environment_contract_sha256"]]),
    figure7_environment_stage_contract_sha256(
      lock_path,
      "r",
      "scvelo"
    )
  )
  testthat::expect_identical(
    unname(
      scvelo_dependencies[["python_environment_contract_sha256"]]
    ),
    figure7_environment_stage_contract_sha256(
      lock_path,
      "python",
      "scvelo"
    )
  )

  copied_lock <- file.path(root, "environment_lock.tsv")
  testthat::expect_true(file.copy(lock_path, copied_lock))
  copied_config <- config
  copied_config$raw_data$environment_lock <- copied_lock
  before <- figure7_cell_table_dependency_values(
    paths,
    copied_config,
    config_path
  )
  write(
    "r\taudit_only\tunconsumed-package\t1.0.0",
    copied_lock,
    append = TRUE
  )
  after <- figure7_cell_table_dependency_values(
    paths,
    copied_config,
    config_path
  )
  testthat::expect_false(identical(
    unname(before[["audit_environment_lock_sha256"]]),
    unname(after[["audit_environment_lock_sha256"]])
  ))
  testthat::expect_identical(
    unname(before[["r_environment_contract_sha256"]]),
    unname(after[["r_environment_contract_sha256"]])
  )
  corrupted <- lock
  corrupted$version[
    corrupted$ecosystem == "r" &
      corrupted$stage == "seurat_upstream" &
      corrupted$package == "yaml"
  ] <- "0.0.0"
  utils::write.table(
    corrupted,
    copied_lock,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  testthat::expect_error(
    figure7_read_environment_lock(copied_lock),
    "changes a reviewed exact pin"
  )
})

testthat::test_that("scVelo cache contract binds every consumed metadata field", {
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_read_config(config_path)
  baseline <- figure7_scvelo_config_contract(config)

  changed_ploidy <- config
  changed_ploidy$si_figures$ploidy_field <- "different_ploidy"
  testthat::expect_false(identical(
    baseline,
    figure7_scvelo_config_contract(changed_ploidy)
  ))

  changed_context <- config
  changed_context$si_figures$context_field <- "different_context"
  testthat::expect_false(identical(
    baseline,
    figure7_scvelo_config_contract(changed_context)
  ))

  unrelated <- config
  unrelated$tgi$day <- unrelated$tgi$day + 1L
  testthat::expect_identical(
    baseline,
    figure7_scvelo_config_contract(unrelated)
  )
})

testthat::test_that("validated A-E caches never invoke Seurat selection", {
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_read_config(config_path)
  paths <- figure7_workflow_paths(
    list(
      intermediate_dir = tempfile("figure7_no_touch_intermediates_"),
      raw_data_dir = tempfile("figure7_no_touch_raw_"),
      download_missing_raw = "false"
    ),
    repo_root,
    tempfile("figure7_no_touch_output_"),
    config
  )
  result <- local({
    original <- figure7_select_seurat_source
    on.exit(assign(
      "figure7_select_seurat_source",
      original,
      envir = .GlobalEnv
    ))
    assign(
      "figure7_select_seurat_source",
      function(...) stop("Seurat selection was touched"),
      envir = .GlobalEnv
    )
    figure7_preflight_workflow(
      paths,
      config_path,
      config,
      include_panel_f = FALSE
    )
  })
  testthat::expect_false(result$needs_seurat)
  testthat::expect_identical(
    result$seurat_selection$source,
    "not_required_existing_validated_caches"
  )
})

testthat::test_that("scVelo and state caches attest absent Seurat ancestors", {
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_read_config(config_path)
  root <- tempfile("figure7_attested_ancestors_")
  dir.create(root)

  scvelo <- list(
    seurat_rds = file.path(root, "selected_seurat.rds"),
    seurat_rds_lineage_active = TRUE,
    raw_manifest = file.path(module_dir, "zenodo_required_files.tsv"),
    raw_data_dir = file.path(root, "raw"),
    loom_root = file.path(root, "missing_loom"),
    loom_root_explicit = FALSE,
    analysis_jobs = 1L,
    scvelo_metrics = file.path(root, "scvelo.csv"),
    scvelo_stage_manifest = file.path(root, "scvelo_manifest.tsv")
  )
  writeLines("reviewed reconstructed RDS", scvelo$seurat_rds)
  writeLines("synthetic scVelo metrics", scvelo$scvelo_metrics)
  loom_sha_keys <- sub(
    "^loom_md5:",
    "loom_sha256:",
    names(figure7_loom_manifest_dependencies(scvelo))
  )
  scvelo_values <- c(
    figure7_scvelo_dependency_values(scvelo, config_path, config),
    stats::setNames(
      rep(paste(rep("a", 64L), collapse = ""), length(loom_sha_keys)),
      loom_sha_keys
    ),
    python_executable = "python",
    python_package_versions = "synthetic",
    r_package_versions = "synthetic",
    scvelo_metrics_sha256 = figure7_sha256(scvelo$scvelo_metrics)
  )
  figure7_write_stage_manifest(
    scvelo_values,
    scvelo$scvelo_stage_manifest
  )
  testthat::expect_true(
    figure7_scvelo_bundle_valid(scvelo, config_path, config)
  )
  unlink(scvelo$seurat_rds)
  testthat::expect_true(
    figure7_scvelo_bundle_valid(scvelo, config_path, config)
  )
  writeLines("changed present RDS", scvelo$seurat_rds)
  testthat::expect_false(
    figure7_scvelo_bundle_valid(scvelo, config_path, config)
  )

  state_root <- file.path(root, "state")
  state <- list(
    state_pathway_root = state_root,
    state_stage_manifest = file.path(
      state_root,
      "00_manifest",
      "figure7_stage_manifest.tsv"
    ),
    cellcycle = file.path(root, "cellcycle.csv"),
    noncellcycle = file.path(root, "noncellcycle.csv"),
    seurat_rds = file.path(root, "state_seurat.rds"),
    seurat_rds_lineage_active = TRUE
  )
  for (path in figure7_state_result_files(state_root)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(paste("state output", basename(path)), path)
  }
  writeLines("cellcycle", state$cellcycle)
  writeLines("noncellcycle", state$noncellcycle)
  writeLines("state reconstructed RDS", state$seurat_rds)
  figure7_write_state_stage_manifest(state, config_path, config)
  testthat::expect_true(
    figure7_state_results_match_inputs(state, config_path, config)
  )

  complete <- figure7_read_key_value_file(state$state_stage_manifest)
  complete <- complete[names(complete) != "stage_fingerprint"]
  testthat::expect_identical(
    unname(complete[["state_result_contract_sha256"]]),
    figure7_state_result_contract_sha256()
  )
  species_audit_key <-
    "output_sha256:00_manifest/feature_species_audit.csv"
  testthat::expect_identical(
    unname(complete[[species_audit_key]]),
    figure7_sha256(file.path(
      state_root,
      "00_manifest",
      "feature_species_audit.csv"
    ))
  )
  species_audit_path <- file.path(
    state_root,
    "00_manifest",
    "feature_species_audit.csv"
  )
  species_audit_contents <- readLines(species_audit_path, warn = FALSE)
  writeLines("changed species audit", species_audit_path)
  testthat::expect_false(
    figure7_state_results_match_inputs(state, config_path, config)
  )
  writeLines(species_audit_contents, species_audit_path)
  testthat::expect_true(
    figure7_state_results_match_inputs(state, config_path, config)
  )
  incomplete <- complete[names(complete) != "seurat_rds_sha256"]
  figure7_write_stage_manifest(incomplete, state$state_stage_manifest)
  testthat::expect_false(
    figure7_state_results_match_inputs(state, config_path, config)
  )
  figure7_write_stage_manifest(complete, state$state_stage_manifest)
  legacy <- complete[names(complete) != "state_result_contract_sha256"]
  figure7_write_stage_manifest(legacy, state$state_stage_manifest)
  testthat::expect_false(
    figure7_state_results_match_inputs(state, config_path, config)
  )
  figure7_write_stage_manifest(complete, state$state_stage_manifest)

  unlink(c(state$cellcycle, state$noncellcycle, state$seurat_rds))
  testthat::expect_true(
    figure7_state_results_match_inputs(state, config_path, config)
  )
  writeLines("changed present state RDS", state$seurat_rds)
  testthat::expect_false(
    figure7_state_results_match_inputs(state, config_path, config)
  )
})

testthat::test_that("generated references attest a removed state tree", {
  config_path <- file.path(module_dir, "figure7_config.yaml")
  config <- figure7_read_config(config_path)
  reviewed <- figure7_test_state_reference()
  root <- tempfile("figure7_generated_reference_")
  dir.create(root)
  reference <- file.path(
    root,
    as.character(config$state_pathways$generated_reference_id)
  )
  dir.create(reference)
  compact_files <- setdiff(
    figure7_generated_state_required_files(),
    "state_pathway_provenance.tsv"
  )
  testthat::expect_true(all(file.copy(
    file.path(reviewed$path, compact_files),
    file.path(reference, compact_files)
  )))
  ranking_path <- file.path(
    reference,
    "state_pathway_gene_ranking_complete.tsv"
  )
  ranking_rows <- seq_len(10001L)
  figure7_write_tsv(
    data.frame(
      rank = ranking_rows,
      gene_id = paste0("GRCh38-GENE", ranking_rows),
      gene_symbol = paste0("GENE", ranking_rows),
      moderated_t = seq(5, -5, length.out = length(ranking_rows)),
      stringsAsFactors = FALSE
    ),
    ranking_path
  )

  cellcycle <- file.path(root, "cellcycle.csv")
  noncellcycle <- file.path(root, "noncellcycle.csv")
  seurat_rds <- file.path(root, "reconstructed.rds")
  writeLines("cellcycle lineage", cellcycle)
  writeLines("noncellcycle lineage", noncellcycle)
  writeLines("reconstructed Seurat lineage", seurat_rds)
  compact_hashes <- vapply(
    file.path(reference, compact_files),
    figure7_sha256,
    character(1L)
  )
  compact_hash_keys <- paste0(
    sub("[.]tsv$", "", compact_files),
    "_sha256"
  )
  provenance <- c(
    reference_kind =
      as.character(config$state_pathways$generated_reference_kind),
    canonical_publication_allowed = "false",
    generated_reference_id =
      as.character(config$state_pathways$generated_reference_id),
    source_code_revision =
      paste0("sha256:", paste(rep("a", 64L), collapse = "")),
    seurat_rds_sha256 = figure7_sha256(seurat_rds),
    cellcycle_metadata_sha256 = figure7_sha256(cellcycle),
    noncellcycle_metadata_sha256 = figure7_sha256(noncellcycle),
    assay = as.character(config$state_pathways$assay),
    counts_layer = as.character(config$state_pathways$counts_layer),
    etp_method = as.character(config$etp$method),
    etp_threshold = as.character(config$etp$threshold),
    spline_df = as.character(config$state_pathways$spline_df),
    pseudotime_bins = as.character(config$state_pathways$pseudotime_bins),
    minimum_cells_per_sample_bin = as.character(
      config$state_pathways$minimum_cells_per_sample_bin
    ),
    grid_size = as.character(config$state_pathways$grid_size),
    seed = as.character(config$statistics$seed),
    gene_set_source = "msigdbr",
    gene_set_species = "Homo sapiens",
    gene_set_collections = paste(
      as.character(unlist(config$state_pathways$collections)),
      collapse = ","
    ),
    gene_set_release = as.character(config$gene_sets$database_release),
    msigdbr_package_version =
      as.character(config$gene_sets$package_version),
    gene_set_membership_sha256 = paste(rep("b", 64L), collapse = ""),
    figure7_config_sha256 = figure7_sha256(config_path),
    figure7_config_contract_sha256 =
      figure7_state_config_contract_sha256(config),
    feature_species_policy_id =
      as.character(config$feature_species$policy_id),
    feature_species_policy =
      as.character(config$state_pathways$generated_feature_species_policy),
    human_feature_prefix =
      as.character(config$feature_species$human_prefix),
    mouse_feature_prefix =
      as.character(config$feature_species$mouse_prefix),
    unknown_feature_policy =
      as.character(config$feature_species$unknown_feature_policy),
    n_input_features = "13001",
    n_human_features_retained = "10001",
    n_mouse_features_excluded = "3000",
    n_ambiguous_features = "0",
    feature_species_audit_sha256 = paste(rep("c", 64L), collapse = ""),
    feature_species_policy_code_sha256 = figure7_sha256(file.path(
      module_dir,
      "src",
      "feature_species_policy.R"
    )),
    pathway_selection_rule =
      as.character(config$state_pathways$pathway_selector),
    activity_table_sha256 = figure7_sha256(file.path(
      reference,
      "panel_7F_pathway_activity_plot_data.tsv"
    )),
    stats::setNames(compact_hashes, compact_hash_keys)
  )
  figure7_write_tsv(
    data.frame(
      key = names(provenance),
      value = unname(provenance),
      stringsAsFactors = FALSE
    ),
    file.path(reference, "state_pathway_provenance.tsv")
  )

  state_root <- file.path(root, "state")
  paths <- list(
    saved_reference = reference,
    reference_stage_manifest = file.path(
      root,
      "generated_reference.stage_manifest.tsv"
    ),
    state_pathway_root = state_root,
    state_stage_manifest = file.path(
      state_root,
      "00_manifest",
      "figure7_stage_manifest.tsv"
    ),
    cellcycle = cellcycle,
    noncellcycle = noncellcycle,
    seurat_rds = seurat_rds,
    seurat_rds_lineage_active = TRUE
  )
  dir.create(dirname(paths$state_stage_manifest), recursive = TRUE)
  figure7_write_stage_manifest(
    figure7_state_dependency_values(paths, config_path, config),
    paths$state_stage_manifest
  )
  testthat::expect_silent(figure7_validate_generated_state_reference(
    reference,
    config,
    expected_inputs = list(
      cellcycle = cellcycle,
      noncellcycle = noncellcycle,
      seurat_rds = seurat_rds
    ),
    config_path = config_path
  ))
  figure7_write_reference_stage_manifest(paths, config_path, config)
  testthat::expect_true(figure7_saved_reference_match_inputs(
    reference,
    paths,
    config,
    config_path
  ))
  complete_reference <- figure7_read_key_value_file(
    paths$reference_stage_manifest
  )
  complete_reference <- complete_reference[
    names(complete_reference) != "stage_fingerprint"
  ]
  testthat::expect_identical(
    unname(complete_reference[["state_result_contract_sha256"]]),
    figure7_state_result_contract_sha256()
  )
  legacy_reference <- complete_reference[
    names(complete_reference) != "state_result_contract_sha256"
  ]
  figure7_write_stage_manifest(
    legacy_reference,
    paths$reference_stage_manifest
  )
  testthat::expect_false(figure7_saved_reference_match_inputs(
    reference,
    paths,
    config,
    config_path
  ))
  figure7_write_stage_manifest(
    complete_reference,
    paths$reference_stage_manifest
  )

  unlink(c(
    paths$state_stage_manifest,
    paths$cellcycle,
    paths$noncellcycle,
    paths$seurat_rds
  ))
  testthat::expect_true(figure7_saved_reference_match_inputs(
    reference,
    paths,
    config,
    config_path
  ))
  unrelated <- config
  unrelated$si_figures$plot_shuffle_seed <-
    as.integer(unrelated$si_figures$plot_shuffle_seed) + 1L
  testthat::expect_true(figure7_saved_reference_match_inputs(
    reference,
    paths,
    unrelated,
    config_path
  ))
  writeLines("changed present reconstructed RDS", paths$seurat_rds)
  testthat::expect_false(figure7_saved_reference_match_inputs(
    reference,
    paths,
    config,
    config_path
  ))
})
