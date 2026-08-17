#!/usr/bin/env Rscript

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("R package 'testthat' is required", call. = FALSE)
}

cmd <- commandArgs(FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
test_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
module_dir <- normalizePath(file.path(dirname(test_path), ".."), mustWork = TRUE)
source(file.path(module_dir, "velocity_pseudotime_panel.R"), local = FALSE)

synthetic_velocity_table <- function() {
  counts <- velocity_expected_cluster_counts()
  cluster <- rep(names(counts), counts)
  n <- length(cluster)
  index <- seq_len(n)
  phase <- 2 * pi * (index - 1) / n
  cluster_offset <- match(cluster, names(counts)) - 2
  data.frame(
    cell_id = sprintf("cell-%05d", index),
    sample_id = sprintf("mouse-%02d", ((index - 1L) %% 16L) + 1L),
    cluster = cluster,
    UMAP_1 = cos(phase) * 4 + cluster_offset * 1.5,
    UMAP_2 = sin(phase) * 3 + cluster_offset * 0.8,
    velocity_UMAP_1 = -sin(phase) * 0.12 + 0.01,
    velocity_UMAP_2 = cos(phase) * 0.12 + 0.01,
    velocity_pseudotime = (index - 1) / (n - 1),
    is_root = cluster == "6",
    stringsAsFactors = FALSE
  )
}

testthat::test_that("the panel table fixes the root and required vector schema", {
  table <- synthetic_velocity_table()
  testthat::expect_silent(velocity_validate_panel_table(table))

  missing_vector <- table
  missing_vector$velocity_UMAP_2 <- NULL
  testthat::expect_error(
    velocity_validate_panel_table(missing_vector),
    "columns must be exactly"
  )

  wrong_root <- table
  wrong_root$is_root[wrong_root$cluster == "10"] <- TRUE
  testthat::expect_error(
    velocity_validate_panel_table(wrong_root),
    "every and only cluster 6"
  )

  no_root <- table
  no_root$is_root <- FALSE
  testthat::expect_error(
    velocity_validate_panel_table(no_root),
    "every and only cluster 6"
  )
})

testthat::test_that("a checksum-bound bundle renders all required visual elements", {
  table <- synthetic_velocity_table()
  root <- tempfile("velocity_bundle_")
  dir.create(root, recursive = TRUE)
  table_path <- file.path(root, "velocity.tsv")
  provenance_path <- file.path(root, "velocity.provenance.tsv")
  source_paths <- file.path(
    root,
    c(
      "embedding.tsv", "embedding.lineage.tsv", "embedding.inventory.tsv",
      "cellcycle.csv", "si.csv"
    )
  )
  for (path in source_paths) writeLines(basename(path), path)
  velocity_write_panel_bundle(
    table,
    table_path,
    provenance_path,
    source_paths[[1L]],
    source_paths[[2L]],
    source_paths[[3L]],
    source_paths[[4L]],
    source_paths[[5L]],
    builder_path = file.path(module_dir, "build_velocity_pseudotime_table.R"),
    panel_logic_path = file.path(module_dir, "velocity_pseudotime_panel.R"),
    scvelo_generator_path = file.path(
      module_dir,
      "..",
      "figure7",
      "generate_scvelo_cell_metrics.R"
    )
  )
  testthat::expect_silent(
    velocity_validate_provenance(provenance_path, table_path)
  )
  grid <- velocity_build_grid(table)
  testthat::expect_gt(nrow(grid), 12L)
  testthat::expect_true(all(c("UMAP_1_to", "UMAP_2_to") %in% names(grid)))

  rendered <- velocity_render_panel(
    table_path,
    provenance_path,
    file.path(root, "run"),
    width = 4,
    height = 3.6,
    dpi = 90L
  )
  testthat::expect_true(file.exists(rendered$pdf))
  testthat::expect_true(file.exists(rendered$png))
  testthat::expect_identical(
    rendered$plot$scales$get_scales("colour")$name,
    "Velocity\npseudotime"
  )
  testthat::expect_match(rendered$plot$labels$subtitle, "2,881 cells")
  labels <- velocity_cluster_labels(table)
  testthat::expect_setequal(labels$cluster, c("4c", "6", "10"))
  testthat::expect_identical(labels$label[labels$cluster == "6"], "6  ROOT")

  writeLines("changed SI dependency", source_paths[[5L]])
  testthat::expect_error(
    velocity_generated_bundle_matches_dependencies(
      table_path,
      provenance_path,
      source_paths[[1L]],
      source_paths[[2L]],
      source_paths[[3L]],
      source_paths[[4L]],
      source_paths[[5L]],
      file.path(module_dir, "build_velocity_pseudotime_table.R"),
      file.path(module_dir, "velocity_pseudotime_panel.R"),
      file.path(module_dir, "..", "figure7", "generate_scvelo_cell_metrics.R"),
      tempfile("unused_rds_"),
      tempfile("unused_looms_"),
      tempfile("unused_config_"),
      tempfile("unused_lock_"),
      root
    ),
    "stale input or code hashes"
  )

  tampered <- velocity_read_tsv(provenance_path)
  tampered$value[tampered$key == "root_cluster"] <- "10"
  velocity_write_tsv(tampered, provenance_path)
  testthat::expect_error(
    velocity_validate_provenance(provenance_path, table_path),
    "reviewed contract"
  )
})

testthat::test_that("embedding lineage binds the RDS, 18 looms, and code", {
  root <- tempfile("velocity_lineage_")
  repo_root <- file.path(root, "repo")
  loom_root <- file.path(repo_root, "raw", "looms")
  dir.create(loom_root, recursive = TRUE)
  embedding <- file.path(repo_root, "cache", "embedding.tsv")
  dir.create(dirname(embedding), recursive = TRUE)
  velocity_write_tsv(
    data.frame(cell = sprintf("cell-%05d", seq_len(35513L))),
    embedding
  )
  for (index in seq_len(18L)) {
    writeLines(
      paste("loom", index),
      file.path(loom_root, sprintf("sample-%02d.loom", index))
    )
  }
  seurat_rds <- file.path(repo_root, "raw", "final.rds")
  generator <- file.path(repo_root, "code", "generator.R")
  config <- file.path(repo_root, "code", "config.yaml")
  lock <- file.path(repo_root, "code", "environment_lock.tsv")
  dir.create(dirname(generator), recursive = TRUE)
  for (path in c(seurat_rds, generator, config, lock)) writeLines(path, path)
  lineage <- file.path(repo_root, "cache", "embedding.lineage.tsv")
  inventory <- file.path(repo_root, "cache", "embedding.inputs.tsv")
  velocity_write_embedding_lineage(
    embedding,
    lineage,
    inventory,
    seurat_rds,
    loom_root,
    generator,
    config,
    lock,
    repo_root
  )
  testthat::expect_silent(velocity_validate_embedding_lineage(
    embedding,
    lineage,
    inventory,
    seurat_rds,
    loom_root,
    generator,
    config,
    lock,
    repo_root
  ))
  writeLines("mutated", file.path(loom_root, "sample-07.loom"))
  testthat::expect_error(
    velocity_validate_embedding_lineage(
      embedding,
      lineage,
      inventory,
      seurat_rds,
      loom_root,
      generator,
      config,
      lock,
      repo_root
    ),
    "differs from the current 18 loom files"
  )
})

testthat::test_that("portable locators reject absolute and traversal reuse", {
  repo_root <- tempfile("velocity_repo_")
  dir.create(file.path(repo_root, "Data"), recursive = TRUE)
  source <- file.path(repo_root, "Data", "table.tsv")
  writeLines("x", source)
  testthat::expect_identical(
    velocity_portable_locator(source, repo_root),
    "Data/table.tsv"
  )
  testthat::expect_identical(
    velocity_resolve_repo_locator("Data/table.tsv", repo_root),
    normalizePath(source, mustWork = TRUE)
  )
  testthat::expect_error(
    velocity_resolve_repo_locator(source, repo_root),
    "Unsafe or non-repository"
  )
  testthat::expect_error(
    velocity_resolve_repo_locator("../table.tsv", repo_root),
    "Unsafe or non-repository"
  )
  external <- tempfile("velocity_external_")
  writeLines("x", external)
  testthat::expect_match(
    velocity_portable_locator(external, repo_root),
    "^external:",
    perl = TRUE
  )
})

testthat::test_that("the source stage exports the reviewed UMAP velocity embedding", {
  worker_path <- file.path(
    module_dir,
    "..",
    "figure7",
    "generate_scvelo_cell_metrics.R"
  )
  env <- new.env(parent = globalenv())
  sys.source(worker_path, envir = env)
  worker <- env$embedded_scvelo_worker()
  testthat::expect_match(
    worker,
    'scv.tl.velocity_embedding(adata, basis="umap")',
    fixed = TRUE
  )
  testthat::expect_match(worker, '"velocity_UMAP_1"', fixed = TRUE)
  testthat::expect_match(worker, '"is_root"', fixed = TRUE)
  testthat::expect_match(worker, 'args.root_key_name', fixed = TRUE)
})
