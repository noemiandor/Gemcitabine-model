workflow_paths_fixture <- function(root) {
  paths <- list(
    intermediate_dir = root,
    raw_data_dir = file.path(root, "raw"),
    raw_manifest = file.path(module_dir, "zenodo_required_files.tsv"),
    download_missing_raw = TRUE,
    download_workers = 4L,
    download_connections_per_file = 2L,
    loom_root_explicit = TRUE,
    seurat_rds_explicit = TRUE,
    scvelo_metrics = file.path(root, "scvelo_cell_metrics.csv"),
    seurat_metadata = file.path(root, "seurat_metadata.csv"),
    seurat_metadata_provenance = file.path(root, "seurat_metadata_provenance.tsv"),
    cellcycle = file.path(root, "CellCycleCells.csv"),
    noncellcycle = file.path(root, "NonCellCycleCells.csv"),
    loom_root = file.path(root, "velocyto_loom"),
    seurat_rds = file.path(root, "seurat.rds"),
    cell_ploidy = file.path(root, "all_ploidy.tsv"),
    sample_info = file.path(root, "sample_info.xlsx"),
    growth_curve = file.path(root, "growth.xlsx"),
    python = file.path(R.home("bin"), "Rscript"),
    scvelo_work_dir = file.path(root, "scvelo_work"),
    state_pathway_root = file.path(root, "state_pathways"),
    saved_reference = file.path(root, "saved_reference"),
    log_dir = file.path(root, "logs")
  )
  for (path in c(paths$seurat_rds, paths$cell_ploidy, paths$sample_info, paths$growth_curve)) {
    writeLines("fixture", path)
  }
  dir.create(paths$loom_root)
  paths
}

testthat::test_that("scVelo input-bundle preflight reuses a complete bundle", {
  root <- tempfile("figure7_scvelo_pair_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  writeLines("metrics", paths$scvelo_metrics)
  writeLines("metadata", paths$seurat_metadata)
  writeLines("provenance", paths$seurat_metadata_provenance)
  state <- figure7_preflight_scvelo_inputs(paths)
  testthat::expect_identical(state$state, "scvelo_input_bundle_ready")
  testthat::expect_true(state$pair_ready)
  testthat::expect_false(state$needs_generate)
})

testthat::test_that("scVelo input-bundle preflight rejects a partial bundle", {
  root <- tempfile("figure7_scvelo_partial_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  writeLines("metrics", paths$scvelo_metrics)
  testthat::expect_error(
    figure7_preflight_scvelo_inputs(paths),
    "Incomplete scVelo/Seurat metadata input bundle"
  )
})

testthat::test_that("scVelo input-bundle preflight schedules generation from explicit raw inputs", {
  root <- tempfile("figure7_scvelo_generate_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  state <- figure7_preflight_scvelo_inputs(paths)
  testthat::expect_identical(state$state, "scvelo_input_bundle_generation_required")
  testthat::expect_false(state$pair_ready)
  testthat::expect_true(state$needs_generate)
  testthat::expect_false(state$needs_raw_stage)
  testthat::expect_identical(state$raw_data_status, "explicit_local_inputs")
})

testthat::test_that("full-workflow preflight starts from an existing cell-table pair", {
  root <- tempfile("figure7_workflow_pair_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  writeLines("cellcycle", paths$cellcycle)
  writeLines("noncellcycle", paths$noncellcycle)
  state <- figure7_preflight_workflow(paths, include_panel_f = TRUE)
  testthat::expect_identical(state$state, "cell_tables_ready")
  testthat::expect_false(state$needs_scvelo)
  testthat::expect_false(state$needs_cell_tables)
})

testthat::test_that("full-workflow preflight starts from existing scVelo metrics", {
  root <- tempfile("figure7_workflow_scvelo_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  writeLines("metrics", paths$scvelo_metrics)
  state <- figure7_preflight_workflow(paths, include_panel_f = TRUE)
  testthat::expect_identical(state$state, "scvelo_ready")
  testthat::expect_false(state$needs_scvelo)
  testthat::expect_true(state$needs_cell_tables)
})

testthat::test_that("full-workflow preflight starts from loom and Seurat when intermediates are absent", {
  root <- tempfile("figure7_workflow_raw_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  state <- figure7_preflight_workflow(paths, include_panel_f = TRUE)
  testthat::expect_identical(state$state, "raw_inputs_ready")
  testthat::expect_true(state$needs_scvelo)
  testthat::expect_true(state$needs_cell_tables)
  testthat::expect_false(state$needs_raw_stage)
  testthat::expect_identical(state$raw_data_status, "explicit_local_inputs")
})

testthat::test_that("missing default raw inputs schedule a pinned Zenodo download", {
  root <- tempfile("figure7_workflow_zenodo_legacy_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  unlink(paths$loom_root, recursive = TRUE)
  unlink(paths$seurat_rds)
  paths$loom_root_explicit <- FALSE
  paths$seurat_rds_explicit <- FALSE
  paths$loom_root <- file.path(paths$raw_data_dir, "velocyto_loom")
  paths$seurat_rds <- file.path(paths$raw_data_dir, "integrated_sct_cca_seurat_final_reclustered.rds")
  state <- figure7_preflight_workflow(paths, include_panel_f = TRUE)
  testthat::expect_true(state$needs_raw_stage)
  testthat::expect_true(state$needs_raw_download)
  testthat::expect_identical(state$raw_download_roles, c("loom", "seurat_rds"))
  testthat::expect_identical(state$raw_data_status, "planned_zenodo_download")

  paths$download_missing_raw <- FALSE
  testthat::expect_error(
    figure7_preflight_workflow(paths, include_panel_f = TRUE),
    "automatic download is disabled"
  )
})

testthat::test_that("an explicit missing raw path fails instead of silently downloading", {
  root <- tempfile("figure7_workflow_explicit_missing_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  unlink(paths$seurat_rds)
  testthat::expect_error(
    figure7_preflight_workflow(paths, include_panel_f = TRUE),
    "Explicit --seurat-rds does not exist"
  )
})

testthat::test_that("full-workflow rejects a partial cell-table pair", {
  root <- tempfile("figure7_workflow_partial_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  writeLines("cellcycle", paths$cellcycle)
  testthat::expect_error(
    figure7_preflight_workflow(paths, include_panel_f = TRUE),
    "Incomplete cell-level intermediate pair"
  )
})

testthat::test_that("run_figure7 exposes preflight-only without creating a run", {
  input <- figure7_test_inputs()
  output_dir <- tempfile("figure7_preflight_output_")
  cellcycle <- file.path(
    repo_root, "Data", "in-vivo", "figure7", "processed",
    "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  )
  noncellcycle <- file.path(
    repo_root, "Data", "in-vivo", "figure7", "processed",
    "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      file.path(module_dir, "run_figure7.R"),
      "--mode=full-workflow", "--panel-set=a-e", "--preflight-only=true",
      paste0("--config=", file.path(module_dir, "figure7_config.yaml")),
      paste0("--cellcycle-input=", cellcycle),
      paste0("--non-cellcycle-input=", noncellcycle),
      paste0("--output-dir=", output_dir)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  testthat::expect_null(attr(status, "status"))
  testthat::expect_true(any(grepl("workflow_state\\tcell_tables_ready", status)))
  testthat::expect_false(dir.exists(output_dir))
})

testthat::test_that("prepare-scvelo-inputs reuses an existing bundle and records provenance", {
  root <- tempfile("figure7_prepare_scvelo_"); dir.create(root)
  intermediate <- file.path(root, "intermediates"); dir.create(intermediate)
  scvelo <- file.path(intermediate, "scvelo_cell_metrics.csv")
  metadata <- file.path(intermediate, "seurat_metadata.csv")
  metadata_provenance <- file.path(intermediate, "seurat_metadata_provenance.tsv")
  output_dir <- file.path(root, "prep_run")
  writeLines(c("cell,TN,clusters,Ploidy,Dose", "c1,Tumor,6,2N,0mg/kg"), scvelo)
  writeLines(c("cell,UMAP_1,UMAP_2,Dose,clusters,sample", "c1,0,0,0mg/kg,6,s1"), metadata)
  writeLines(c("key\tvalue", "fixture\ttrue"), metadata_provenance)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      file.path(module_dir, "run_figure7.R"),
      "--mode=prepare-scvelo-inputs",
      paste0("--config=", file.path(module_dir, "figure7_config.yaml")),
      paste0("--intermediate-dir=", intermediate),
      paste0("--scvelo-metrics=", scvelo),
      paste0("--seurat-metadata-output=", metadata),
      paste0("--seurat-metadata-provenance-output=", metadata_provenance),
      paste0("--output-dir=", output_dir)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  testthat::expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))
  run_config <- figure7_read_tsv(
    file.path(output_dir, "metadata", "run_config.tsv"), c("key", "value")
  )
  values <- stats::setNames(run_config$value, run_config$key)
  testthat::expect_identical(unname(values[["module"]]), "in_vivo_figure7_input_prep")
  testthat::expect_identical(unname(values[["workflow_initial_state"]]), "scvelo_input_bundle_ready")
  testthat::expect_identical(unname(values[["workflow_executed_stages"]]), "none")
  testthat::expect_identical(unname(values[["workflow_scvelo_metrics"]]), normalizePath(scvelo))
  testthat::expect_identical(unname(values[["workflow_seurat_metadata"]]), normalizePath(metadata))
  testthat::expect_identical(
    unname(values[["workflow_seurat_metadata_provenance"]]),
    normalizePath(metadata_provenance)
  )
  testthat::expect_false(dir.exists(file.path(output_dir, "figures")))
})

testthat::test_that("raw full-workflow dispatches scVelo before cell-level generation", {
  root <- tempfile("figure7_workflow_dispatch_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  preflight <- figure7_preflight_workflow(paths, include_panel_f = FALSE)
  calls <- character()
  original_runner <- get("figure7_run_stage", envir = .GlobalEnv)
  on.exit(assign("figure7_run_stage", original_runner, envir = .GlobalEnv), add = TRUE)
  assign(
    "figure7_run_stage",
    function(label, script, values, log_path) {
      calls <<- c(calls, label)
      if (grepl("scVelo", label, fixed = TRUE)) {
        writeLines("metrics", paths$scvelo_metrics)
        writeLines("metadata", paths$seurat_metadata)
        writeLines("provenance", paths$seurat_metadata_provenance)
      }
      if (grepl("CellCycle", label, fixed = TRUE)) {
        writeLines("cellcycle", paths$cellcycle)
        writeLines("noncellcycle", paths$noncellcycle)
      }
      invisible(log_path)
    },
    envir = .GlobalEnv
  )
  test_config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
  test_config$raw_data$seurat_rds_sha256 <- figure7_sha256(paths$seurat_rds)
  result <- figure7_prepare_full_workflow(
    args = list(), paths = paths, preflight = preflight,
    script_dir = module_dir, config_path = file.path(module_dir, "figure7_config.yaml"),
    config = test_config,
    include_panel_f = FALSE
  )
  testthat::expect_identical(
    calls,
    c("Generate scVelo cell metrics", "Generate CellCycle and NonCellCycle cell-level inputs")
  )
  testthat::expect_identical(result$executed_stages, c("scvelo_metrics", "celllevel_inputs"))
  testthat::expect_identical(result$seurat_metadata, paths$seurat_metadata)
  testthat::expect_identical(
    result$seurat_metadata_provenance,
    paths$seurat_metadata_provenance
  )
})

testthat::test_that("automatic Zenodo materialization runs before scVelo", {
  root <- tempfile("figure7_workflow_download_dispatch_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  unlink(paths$loom_root, recursive = TRUE)
  unlink(paths$seurat_rds)
  paths$loom_root_explicit <- FALSE
  paths$seurat_rds_explicit <- FALSE
  paths$loom_root <- file.path(paths$raw_data_dir, "velocyto_loom")
  paths$seurat_rds <- file.path(paths$raw_data_dir, "integrated_sct_cca_seurat_final_reclustered.rds")
  preflight <- figure7_preflight_workflow(paths, include_panel_f = FALSE)
  calls <- character()
  original_runner <- get("figure7_run_stage", envir = .GlobalEnv)
  on.exit(assign("figure7_run_stage", original_runner, envir = .GlobalEnv), add = TRUE)
  assign(
    "figure7_run_stage",
    function(label, script, values, log_path) {
      calls <<- c(calls, label)
      if (grepl("raw data", label, fixed = TRUE)) {
        dir.create(paths$loom_root, recursive = TRUE)
        writeLines("loom", file.path(paths$loom_root, "fixture.loom"))
        writeLines("rds", paths$seurat_rds)
        provenance <- file.path(paths$raw_data_dir, "provenance")
        dir.create(provenance, recursive = TRUE)
        figure7_write_tsv(
          data.frame(
            role = c("loom", "seurat_rds"), filename = c("fixture.loom", basename(paths$seurat_rds)),
            action = "downloaded", md5 = paste(rep("a", 32), collapse = ""),
            sha256 = paste(rep("b", 64), collapse = ""), stringsAsFactors = FALSE
          ),
          file.path(provenance, "downloaded_files_checksums.tsv")
        )
      }
      if (grepl("scVelo", label, fixed = TRUE)) {
        writeLines("metrics", paths$scvelo_metrics)
        writeLines("metadata", paths$seurat_metadata)
        writeLines("provenance", paths$seurat_metadata_provenance)
      }
      if (grepl("CellCycle", label, fixed = TRUE)) {
        writeLines("cellcycle", paths$cellcycle)
        writeLines("noncellcycle", paths$noncellcycle)
      }
      invisible(log_path)
    },
    envir = .GlobalEnv
  )
  result <- figure7_prepare_full_workflow(
    args = list(), paths = paths, preflight = preflight,
    script_dir = module_dir, config_path = file.path(module_dir, "figure7_config.yaml"),
    config = figure7_read_config(file.path(module_dir, "figure7_config.yaml")),
    include_panel_f = FALSE
  )
  testthat::expect_identical(
    calls,
    c(
      "Download and verify Figure 7 raw data from Zenodo",
      "Generate scVelo cell metrics",
      "Generate CellCycle and NonCellCycle cell-level inputs"
    )
  )
  testthat::expect_identical(result$executed_stages, c("raw_data_download", "scvelo_metrics", "celllevel_inputs"))
  testthat::expect_identical(result$raw_data_status, "downloaded_from_zenodo")
})

testthat::test_that("full-workflow panel F exports generated sources without frozen byte checks", {
  root <- tempfile("figure7_workflow_panel_f_"); dir.create(root)
  paths <- workflow_paths_fixture(root)
  writeLines("cellcycle", paths$cellcycle)
  writeLines("noncellcycle", paths$noncellcycle)
  preflight <- figure7_preflight_workflow(paths, include_panel_f = TRUE)
  calls <- list()
  original_runner <- get("figure7_run_stage", envir = .GlobalEnv)
  on.exit(assign("figure7_run_stage", original_runner, envir = .GlobalEnv), add = TRUE)
  assign(
    "figure7_run_stage",
    function(label, script, values, log_path) {
      calls[[length(calls) + 1L]] <<- list(label = label, script = script, values = values)
      if (grepl("Export compact", label, fixed = TRUE)) dir.create(paths$saved_reference, recursive = TRUE)
      invisible(log_path)
    },
    envir = .GlobalEnv
  )
  test_config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
  test_config$raw_data$seurat_rds_sha256 <- figure7_sha256(paths$seurat_rds)
  result <- figure7_prepare_full_workflow(
    args = list(), paths = paths, preflight = preflight,
    script_dir = module_dir, config_path = file.path(module_dir, "figure7_config.yaml"),
    config = test_config,
    include_panel_f = TRUE
  )
  testthat::expect_identical(
    vapply(calls, `[[`, character(1L), "label"),
    c("Generate state-pathway support results", "Export compact state-pathway reference")
  )
  testthat::expect_identical(calls[[2L]]$values$verify_source_checksums, "FALSE")
  testthat::expect_identical(result$executed_stages, c("state_pathway_support", "state_pathway_export"))
})

testthat::test_that("pinned Zenodo manifest contains the complete raw-data contract", {
  manifest <- utils::read.delim(
    file.path(module_dir, "zenodo_required_files.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character"
  )
  testthat::expect_equal(nrow(manifest), 19L)
  testthat::expect_equal(sum(manifest$role == "loom"), 18L)
  testthat::expect_equal(sum(as.numeric(manifest$size_bytes)), 11395115098)
  rds <- manifest[manifest$role == "seurat_rds", , drop = FALSE]
  testthat::expect_identical(
    rds$sha256,
    "727b8a5e5868da911c3b0873838fb1b0023377ed21ea5498dc6493acbbef6d98"
  )
})

testthat::test_that("aria2c input dispatches parallel files and split connections", {
  env <- new.env(parent = globalenv())
  expressions <- parse(file.path(module_dir, "download_figure7_raw_data.R"))
  for (expression in expressions) eval(expression, envir = env)
  root <- tempfile("figure7_fake_aria2_"); dir.create(root)
  args_log <- file.path(root, "aria2_args.txt")
  fake_aria2 <- file.path(root, "aria2c")
  writeLines(
    c(
      "#!/bin/sh",
      paste0("args_log=", shQuote(args_log)),
      "printf '%s\\n' \"$@\" > \"$args_log\"",
      "input_file=",
      "for arg in \"$@\"; do case \"$arg\" in --input-file=*) input_file=${arg#*=};; esac; done",
      "dir=",
      "while IFS= read -r line; do",
      "  case \"$line\" in",
      "    '  dir='*) dir=${line#  dir=};;",
      "    '  out='*) out=${line#  out=}; mkdir -p \"$dir\"; printf fixture > \"$dir/$out\";;",
      "  esac",
      "done < \"$input_file\""
    ),
    fake_aria2
  )
  Sys.chmod(fake_aria2, mode = "0755")
  selected <- data.frame(
    url = c("https://example.invalid/a", "https://example.invalid/b"),
    path = c(file.path(root, "raw", "a.loom"), file.path(root, "raw", "b.loom")),
    filename = c("a.loom", "b.loom"),
    size_bytes = c(7, 7), stringsAsFactors = FALSE
  )
  env$download_parallel_aria2(
    selected, 1:2, fake_aria2, download_workers = 4L,
    connections_per_file = 2L, raw_data_dir = file.path(root, "raw")
  )
  testthat::expect_true(all(file.exists(paste0(selected$path, ".part"))))
  observed_args <- readLines(args_log, warn = FALSE)
  testthat::expect_true("--max-concurrent-downloads=4" %in% observed_args)
  testthat::expect_true("--split=2" %in% observed_args)
  testthat::expect_true("--max-connection-per-server=2" %in% observed_args)
  testthat::expect_error(env$parse_positive_integer("0", "download-workers", 16L), "integer from 1 to 16")
  testthat::expect_error(env$parse_positive_integer("2.5", "download-workers", 16L), "integer from 1 to 16")
})

testthat::test_that("scVelo multiprocessing uses a short local socket path", {
  env <- new.env(parent = globalenv())
  expressions <- parse(file.path(module_dir, "generate_scvelo_cell_metrics.R"))
  for (expression in head(expressions, -1L)) eval(expression, envir = env)
  path <- env$create_scvelo_multiprocessing_tmpdir()
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)
  testthat::expect_identical(dirname(path), normalizePath("/tmp", mustWork = TRUE))
  testthat::expect_true(startsWith(basename(path), "f7mp_"))
  testthat::expect_lt(nchar(path), 80L)
  testthat::expect_equal(unname(file.access(path, mode = 2L)), 0L)
})

testthat::test_that("raw-data downloader materializes, reuses, and rejects corrupt offline fixtures", {
  root <- tempfile("figure7_download_fixture_"); dir.create(root)
  source_dir <- file.path(root, "source"); dir.create(source_dir)
  raw_dir <- file.path(root, "raw")
  loom_source <- file.path(source_dir, "sample-one.loom")
  rds_source <- file.path(source_dir, "fixture.rds")
  writeBin(charToRaw("loom fixture bytes"), loom_source)
  writeBin(charToRaw("rds fixture bytes"), rds_source)
  fixture_paths <- c(loom_source, rds_source)
  manifest <- data.frame(
    role = c("loom", "seurat_rds"),
    filename = basename(fixture_paths),
    size_bytes = as.character(file.info(fixture_paths)$size),
    md5 = unname(tools::md5sum(fixture_paths)),
    sha256 = c("", vapply(rds_source, figure7_sha256, character(1L))),
    url = paste0("file://", normalizePath(fixture_paths)),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(root, "manifest.tsv")
  utils::write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  downloader <- file.path(module_dir, "download_figure7_raw_data.R")
  command <- c(
    downloader,
    paste0("--raw-data-dir=", raw_dir),
    paste0("--manifest=", manifest_path),
    "--roles=all"
  )
  first <- system2(file.path(R.home("bin"), "Rscript"), command, stdout = TRUE, stderr = TRUE)
  testthat::expect_null(attr(first, "status"), info = paste(first, collapse = "\n"))
  testthat::expect_true(file.exists(file.path(raw_dir, "velocyto_loom", basename(loom_source))))
  testthat::expect_true(file.exists(file.path(raw_dir, basename(rds_source))))
  audit <- figure7_read_tsv(
    file.path(raw_dir, "provenance", "downloaded_files_checksums.tsv"),
    c("filename", "action", "md5", "sha256", "download_client")
  )
  testthat::expect_true(all(audit$action == "downloaded"))
  testthat::expect_true(all(audit$download_client %in% c("curl", "R_libcurl")))

  second <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(command, "--allow-download=false"), stdout = TRUE, stderr = TRUE
  )
  testthat::expect_null(attr(second, "status"), info = paste(second, collapse = "\n"))
  audit <- figure7_read_tsv(file.path(raw_dir, "provenance", "downloaded_files_checksums.tsv"))
  testthat::expect_true(all(audit$action == "reused"))

  writeBin(charToRaw("corrupt"), file.path(raw_dir, basename(rds_source)))
  corrupt <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(command, "--roles=seurat_rds", "--allow-download=false"), stdout = TRUE, stderr = TRUE
  ))
  testthat::expect_true(!is.null(attr(corrupt, "status")) && attr(corrupt, "status") != 0L)
  testthat::expect_match(paste(corrupt, collapse = "\n"), "missing or invalid")
})

testthat::test_that("scVelo loom resolver accepts flat Zenodo and legacy nested layouts", {
  env <- new.env(parent = globalenv())
  expressions <- parse(file.path(module_dir, "generate_scvelo_cell_metrics.R"))
  for (i in seq_len(length(expressions) - 1L)) eval(expressions[[i]], envir = env)
  root <- tempfile("figure7_loom_layout_"); dir.create(root)
  metadata <- data.frame(sample_folder = c("flat", "nested"), stringsAsFactors = FALSE)
  writeLines("flat", file.path(root, "flat.loom"))
  dir.create(file.path(root, "nested"))
  writeLines("nested", file.path(root, "nested", "nested.loom"))
  observed <- env$resolve_loom_files(root, metadata, "sample_folder")
  testthat::expect_identical(
    basename(observed),
    c("flat.loom", "nested.loom")
  )
})

testthat::test_that("persistent Seurat metadata export includes UMAP reduction coordinates", {
  testthat::skip_if_not_installed("Seurat")
  env <- new.env(parent = globalenv())
  expressions <- parse(file.path(module_dir, "generate_scvelo_cell_metrics.R"))
  for (i in seq_len(length(expressions) - 1L)) eval(expressions[[i]], envir = env)
  cells <- c("cell_c", "cell_a", "cell_b")
  counts <- matrix(1, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), cells))
  metadata <- data.frame(
    Dose = c("120", "0", "30"), IDs = c("4N-Tumor", "2N-Tumor", "4N-Tumor"),
    sample_folder = c("s3", "s1", "s2"), cluster_final = c("6", "4c", "10"),
    cluster_cell_cycle_annotation = c(
      "cell_cycle_candidate", "not_cell_cycle_candidate", "cell_cycle_candidate"
    ),
    integrated_snn_res.0.6 = c("6", "4", "10"),
    custom_extra = c("z", "x", "y"), row.names = cells, stringsAsFactors = FALSE
  )
  object <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata)
  umap <- matrix(
    c(3, 30, 1, 10, 2, 20), ncol = 2, byrow = TRUE,
    dimnames = list(cells, c("UMAP_1", "UMAP_2"))
  )
  pca <- matrix(
    c(3, 4, 1, 2, 2, 3), ncol = 2, byrow = TRUE,
    dimnames = list(cells, c("PC_1", "PC_2"))
  )
  object[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = "RNA")
  object[["pca"]] <- Seurat::CreateDimReducObject(embeddings = pca, key = "PC_", assay = "RNA")
  rds <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".csv")
  saveRDS(object, rds)
  invisible(env$build_all_cells_metadata(rds, n_pcs = 2L, seurat_metadata_output = output))
  observed <- readr::read_csv(output, show_col_types = FALSE)
  testthat::expect_identical(observed$cell, sort(cells))
  testthat::expect_identical(observed$UMAP_1, c(1, 2, 3))
  testthat::expect_identical(observed$UMAP_2, c(10, 20, 30))
  testthat::expect_true("custom_extra" %in% names(observed))
})
