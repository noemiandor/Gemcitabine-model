species_env <- new.env(parent = globalenv())
sys.source(
  file.path(module_dir, "src", "feature_species_policy.R"),
  envir = species_env
)

testthat::test_that("feature classifier accepts only exact configured prefixes", {
  annotation <- species_env$figure7_classify_feature_species(c(
    "GRCh38-TP53",
    "GRCm39-Trp53",
    "GRCh38-DUP.1",
    "GRCh38-DUP",
    "grch38-lower",
    "GRCh37-old",
    "GRCh38_wrong",
    "hg38-alias",
    "UNPREFIXED",
    " GRCh38-spaced"
  ))
  testthat::expect_identical(
    annotation$species,
    c(
      "human", "mouse", "human", "human",
      "ambiguous", "ambiguous", "ambiguous", "ambiguous", "ambiguous",
      "ambiguous"
    )
  )
  testthat::expect_identical(
    annotation$symbol[seq_len(4L)],
    c("TP53", "Trp53", "DUP", "DUP")
  )
  testthat::expect_identical(
    annotation$feature[annotation$species == "human"],
    c("GRCh38-TP53", "GRCh38-DUP.1", "GRCh38-DUP")
  )
})

testthat::test_that("human boundary is fail-closed and preserves duplicates/order", {
  testthat::expect_error(
    species_env$figure7_select_human_features(
      c("GRCh38-TP53", "UNPREFIXED"),
      "fixture"
    ),
    "ambiguous or unprefixed"
  )
  selected <- species_env$figure7_select_human_features(
    c("GRCh38-Zeta", "GRCm39-Huge", "GRCh38-Alpha", "GRCh38-Alpha.1"),
    "fixture"
  )
  testthat::expect_identical(
    selected$features,
    c("GRCh38-Zeta", "GRCh38-Alpha", "GRCh38-Alpha.1")
  )
  testthat::expect_identical(
    selected$symbols,
    c("Zeta", "Alpha", "Alpha")
  )
  testthat::expect_identical(
    selected$audit$n_mouse_features_excluded,
    1L
  )
})

testthat::test_that("mouse load is removed before expression summaries", {
  counts <- matrix(
    c(
      10, 20,
      1000000, 1,
      5, 5
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(
      c("GRCh38-A", "GRCm39-MouseLoad", "GRCh38-B"),
      c("cell_1", "cell_2")
    )
  )
  result <- species_env$figure7_filter_human_feature_matrix(
    counts,
    "fixture"
  )
  testthat::expect_identical(
    rownames(result$matrix),
    c("GRCh38-A", "GRCh38-B")
  )
  testthat::expect_equal(colSums(result$matrix), c(cell_1 = 15, cell_2 = 25))
  testthat::expect_false(any(grepl("^GRCm39-", rownames(result$matrix))))
})

testthat::test_that("panel 7F joins every Assay5 counts layer before filtering", {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  if (!"JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    testthat::skip("Installed SeuratObject has no Assay5 JoinLayers support")
  }

  counts <- Matrix::Matrix(
    matrix(
      c(
        10, 20, 30, 40,
        1000, 2000, 3000, 4000,
        5, 6, 7, 8
      ),
      nrow = 3L,
      byrow = TRUE,
      dimnames = list(
        c("GRCh38-A", "GRCm39-MouseLoad", "GRCh38-B"),
        paste0("cell_", seq_len(4L))
      )
    ),
    sparse = TRUE
  )
  object <- Seurat::CreateSeuratObject(counts)
  object$batch <- c("batch_1", "batch_1", "batch_2", "batch_2")
  object[["RNA"]] <- split(object[["RNA"]], f = object$batch)
  testthat::expect_identical(
    SeuratObject::Layers(object[["RNA"]]),
    c("counts.batch_1", "counts.batch_2")
  )

  rds <- tempfile(fileext = ".rds")
  saveRDS(object, rds)
  on.exit(unlink(rds), add = TRUE)
  support_env <- new.env(parent = globalenv())
  sys.source(
    file.path(
      module_dir,
      "generate_pseudotime_state_pathways_support.R"
    ),
    envir = support_env
  )
  loaded <- support_env$load_counts(
    rds,
    assay = "RNA",
    counts_layer = "counts"
  )

  testthat::expect_identical(
    colnames(loaded$counts),
    paste0("cell_", seq_len(4L))
  )
  testthat::expect_identical(
    rownames(loaded$counts),
    c("GRCh38-A", "GRCh38-B")
  )
  testthat::expect_equal(
    as.matrix(loaded$counts),
    as.matrix(counts[c("GRCh38-A", "GRCh38-B"), , drop = FALSE])
  )
  testthat::expect_identical(
    loaded$species_audit$n_mouse_features_excluded,
    1L
  )
})
