scvelo_env <- new.env(parent = globalenv())
sys.source(
  file.path(module_dir, "generate_scvelo_cell_metrics.R"),
  envir = scvelo_env
)

download_env <- new.env(parent = globalenv())
sys.source(
  file.path(module_dir, "download_figure7_raw_data.R"),
  envir = download_env
)

if (!exists("figure7_observed_loom_sha_dependencies", mode = "function")) {
  source(
    file.path(module_dir, "src", "generated_state_pathway_reference.R"),
    local = FALSE
  )
  source(file.path(module_dir, "src", "input_preflight.R"), local = FALSE)
}

testthat::test_that("scVelo output covers every deposited Seurat cell exactly once", {
  metadata <- data.frame(
    cell = c("cell-a", "cell-b", "cell-c"),
    Dose = c("0mg/kg", "30mg/kg", "120mg/kg"),
    Ploidy = c("2N", "4N", "4N"),
    TN = c("tumor", "tumor", "cell-culture"),
    clusters = c("0", "2", "4c"),
    sample = c("A1", "A2", "A5"),
    stringsAsFactors = FALSE
  )
  metrics <- metadata[c(3L, 1L, 2L), , drop = FALSE]
  metrics$velocity_pseudotime <- c(1, 0, 0.5)
  metrics <- metrics[
    ,
    c(
      "cell", "Dose", "Ploidy", "TN", "clusters", "sample",
      "velocity_pseudotime"
    )
  ]
  path <- tempfile("scvelo_metrics_", fileext = ".csv")
  utils::write.csv(metrics, path, row.names = FALSE)
  validated <- scvelo_env$validate_scvelo_metrics_output(
    path,
    metadata,
    "clusters"
  )
  testthat::expect_identical(
    as.character(validated$cell),
    metadata$cell
  )

  duplicate <- metrics
  duplicate$cell[[3L]] <- duplicate$cell[[2L]]
  utils::write.csv(duplicate, path, row.names = FALSE)
  testthat::expect_error(
    scvelo_env$validate_scvelo_metrics_output(
      path,
      metadata,
      "clusters"
    ),
    "exactly one row for every deposited Seurat cell"
  )

  missing <- metrics[-1L, , drop = FALSE]
  utils::write.csv(missing, path, row.names = FALSE)
  testthat::expect_error(
    scvelo_env$validate_scvelo_metrics_output(
      path,
      metadata,
      "clusters"
    ),
    "exactly one row for every deposited Seurat cell"
  )

  nonfinite <- metrics
  nonfinite$velocity_pseudotime[[1L]] <- Inf
  utils::write.csv(nonfinite, path, row.names = FALSE)
  testthat::expect_error(
    scvelo_env$validate_scvelo_metrics_output(
      path,
      metadata,
      "clusters"
    ),
    "finite and within \\[0,1\\]"
  )

  mismatched <- metrics
  mismatched$Ploidy[mismatched$cell == "cell-a"] <- "4N"
  utils::write.csv(mismatched, path, row.names = FALSE)
  testthat::expect_error(
    scvelo_env$validate_scvelo_metrics_output(
      path,
      metadata,
      "clusters"
    ),
    "differs from deposited Seurat metadata: Ploidy"
  )

  worker <- scvelo_env$embedded_scvelo_worker()
  testthat::expect_match(
    worker,
    'matched_meta.shape[0] != meta.shape[0]',
    fixed = TRUE
  )
  testthat::expect_match(
    worker,
    'matched_meta["velocity_obs"].astype(str).duplicated().any()',
    fixed = TRUE
  )
})

testthat::test_that("local loom inputs must satisfy all 18 pinned size and MD5 rows", {
  root <- tempfile("loom_contract_")
  loom_root <- file.path(root, "velocyto_loom")
  dir.create(loom_root, recursive = TRUE)
  filenames <- sprintf("sample_%02d.loom", seq_len(18L))
  contents <- sprintf("loom-%02d", seq_len(18L))
  paths <- file.path(loom_root, filenames)
  for (i in seq_along(paths)) {
    writeChar(contents[[i]], paths[[i]], eos = NULL)
  }
  manifest <- data.frame(
    role = "loom",
    filename = filenames,
    size_bytes = as.numeric(file.info(paths)$size),
    md5 = unname(tools::md5sum(paths)),
    sha256 = vapply(paths, figure7_sha256, character(1L)),
    url = paste0("file://", paths),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(root, "manifest.tsv")
  utils::write.table(
    manifest,
    manifest_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  workflow_paths <- list(
    raw_manifest = manifest_path,
    raw_data_dir = root,
    loom_root = loom_root,
    loom_root_explicit = TRUE
  )
  observed <- figure7_observed_loom_sha_dependencies(workflow_paths)
  testthat::expect_length(observed, 18L)
  testthat::expect_setequal(
    unname(observed),
    manifest$sha256
  )

  writeChar("tamper!", paths[[1L]], eos = NULL)
  testthat::expect_error(
    figure7_observed_loom_sha_dependencies(workflow_paths),
    "differ from the pinned Zenodo size/MD5 contract"
  )
})

testthat::test_that("verified local download atomically replaces corrupt cache", {
  root <- tempfile("download_contract_")
  dir.create(root)
  source <- file.path(root, "source.rds")
  writeChar("reviewed-bytes", source, eos = NULL)
  raw_root <- file.path(root, "raw")
  dir.create(raw_root)
  target <- file.path(raw_root, "fixture.rds")
  writeChar("corrupt-bytes!", target, eos = NULL)
  manifest <- data.frame(
    role = "seurat_rds",
    filename = "fixture.rds",
    size_bytes = as.numeric(file.info(source)$size),
    md5 = unname(tools::md5sum(source)),
    sha256 = figure7_sha256(source),
    url = paste0("file://", normalizePath(source, mustWork = TRUE)),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(root, "manifest.tsv")
  utils::write.table(
    manifest,
    manifest_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  result <- download_env$download_figure7_raw_data(
    raw_data_dir = raw_root,
    manifest_path = manifest_path,
    roles = "seurat_rds",
    aria2_bin = "",
    wget_bin = "",
    curl_bin = "",
    download_workers = 1L,
    connections_per_file = 1L,
    allow_download = TRUE
  )
  testthat::expect_identical(result$downloaded, 1L)
  testthat::expect_identical(readBin(target, "raw", n = 100L), readBin(
    source,
    "raw",
    n = 100L
  ))
  testthat::expect_identical(result$audit$path, "fixture.rds")
  testthat::expect_false(grepl("^/", result$audit$path))
  testthat::expect_false(file.exists(paste0(target, ".part")))

  reused <- download_env$download_figure7_raw_data(
    raw_data_dir = raw_root,
    manifest_path = manifest_path,
    roles = "seurat_rds",
    aria2_bin = "",
    wget_bin = "",
    curl_bin = "",
    download_workers = 1L,
    connections_per_file = 1L,
    allow_download = FALSE
  )
  testthat::expect_identical(reused$reused, 1L)
  testthat::expect_identical(reused$audit$action, "reused")
})
