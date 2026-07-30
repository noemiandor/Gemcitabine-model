support_env <- new.env(parent = globalenv())
sys.source(
  file.path(
    module_dir,
    "generate_pseudotime_state_pathways_support.R"
  ),
  envir = support_env
)

mock_fgsea_frame <- function(pathways, values) {
  keys <- names(pathways)
  out <- data.frame(
    pathway = keys,
    pval = vapply(keys, function(key) values[[key]]$pval, numeric(1L)),
    padj = vapply(keys, function(key) values[[key]]$padj, numeric(1L)),
    log2err = vapply(
      keys,
      function(key) values[[key]]$log2err,
      numeric(1L)
    ),
    ES = vapply(keys, function(key) values[[key]]$ES, numeric(1L)),
    NES = vapply(keys, function(key) values[[key]]$NES, numeric(1L)),
    size = vapply(keys, function(key) values[[key]]$size, numeric(1L)),
    stringsAsFactors = FALSE
  )
  out$leadingEdge <- lapply(keys, function(key) paste0("GENE_", key))
  out
}

testthat::test_that(
  "adaptive GSEA retries only unresolved keys and recomputes collection BH",
  {
    stats <- stats::setNames(6:1, paste0("G", seq_len(6L)))
    pathways <- list(
      A = c("G1", "G2"),
      B = c("G2", "G3"),
      C = c("G3", "G4"),
      D = c("G4", "G5")
    )
    collection <- list(
      collection = "H",
      label = "hallmark",
      sets = pathways
    )
    calls <- list()
    runner <- function(
      pathways,
      stats,
      min_size,
      max_size,
      nperm_simple,
      seed
    ) {
      calls[[length(calls) + 1L]] <<- list(
        keys = names(pathways),
        budget = nperm_simple,
        seed = seed
      )
      finite <- function(
        pval,
        nes,
        log2err = 0.1,
        subset_padj = 0.999
      ) {
        list(
          pval = pval,
          padj = subset_padj,
          log2err = log2err,
          ES = nes / 2,
          NES = nes,
          size = 2
        )
      }
      unresolved <- list(
        pval = NA_real_,
        padj = NA_real_,
        log2err = NA_real_,
        ES = NA_real_,
        NES = NA_real_,
        size = 2
      )
      values <- list(
        A = finite(
          0.02,
          1.5,
          log2err = NA_real_,
          subset_padj = NA_real_
        ),
        B = if (nperm_simple >= 100L) finite(0.001, -2) else unresolved,
        C = finite(0.04, 0.5),
        D = if (nperm_simple >= 1000L) finite(0.003, -3) else unresolved
      )
      mock_fgsea_frame(pathways, values)
    }

    observed <- support_env$run_fgsea_collection(
      stats = stats,
      collection_obj = collection,
      ranking_id = "fixture",
      min_size = 1L,
      max_size = 10L,
      nperm_simple = 10L,
      seed = 7L,
      nperm_simple_max = 1000L,
      nperm_simple_multiplier = 10L,
      fgsea_runner = runner
    )

    testthat::expect_identical(
      lapply(calls, `[[`, "keys"),
      list(c("A", "B", "C", "D"), c("B", "D"), "D")
    )
    testthat::expect_identical(
      vapply(calls, `[[`, integer(1L), "budget"),
      c(10L, 100L, 1000L)
    )
    testthat::expect_identical(
      vapply(calls, `[[`, integer(1L), "seed"),
      rep(7L, 3L)
    )

    by_pathway <- observed[match(names(pathways), observed$pathway), ]
    testthat::expect_equal(by_pathway$pval, c(0.02, 0.001, 0.04, 0.003))
    testthat::expect_equal(
      by_pathway$padj,
      stats::p.adjust(by_pathway$pval, method = "BH")
    )
    testthat::expect_identical(
      as.integer(by_pathway$nPermSimple),
      c(10L, 100L, 10L, 1000L)
    )
    testthat::expect_identical(
      as.integer(by_pathway$retry_round),
      c(0L, 1L, 0L, 2L)
    )
    testthat::expect_identical(
      as.character(by_pathway$direction),
      c("positive", "negative", "positive", "negative")
    )
    testthat::expect_true(is.na(by_pathway$log2err[[1L]]))
    testthat::expect_true(all(support_env$fgsea_core_finite(observed)))
  }
)

testthat::test_that("adaptive GSEA fails closed at its configured cap", {
  stats <- stats::setNames(4:1, paste0("G", seq_len(4L)))
  collection <- list(
    collection = "C2:CP:REACTOME",
    label = "reactome",
    sets = list(A = c("G1", "G2"))
  )
  unresolved_runner <- function(
    pathways,
    stats,
    min_size,
    max_size,
    nperm_simple,
    seed
  ) {
    mock_fgsea_frame(
      pathways,
      list(A = list(
        pval = NA_real_,
        padj = NA_real_,
        log2err = NA_real_,
        ES = NA_real_,
        NES = NA_real_,
        size = 2
      ))
    )
  }
  testthat::expect_error(
    support_env$run_fgsea_collection(
      stats,
      collection,
      "fixture",
      1L,
      10L,
      10L,
      3L,
      100L,
      10L,
      unresolved_runner
    ),
    "exhausted nPermSimple=100.*A"
  )
})

testthat::test_that("adaptive GSEA rejects incomplete retry keys", {
  stats <- stats::setNames(4:1, paste0("G", seq_len(4L)))
  collection <- list(
    collection = "C5:GO:BP",
    label = "go_bp",
    sets = list(A = c("G1", "G2"), B = c("G2", "G3"))
  )
  incomplete_runner <- function(
    pathways,
    stats,
    min_size,
    max_size,
    nperm_simple,
    seed
  ) {
    values <- list(
      A = list(
        pval = 0.1,
        padj = 0.2,
        log2err = 0.1,
        ES = 0.5,
        NES = 1,
        size = 2
      ),
      B = list(
        pval = 0.2,
        padj = 0.2,
        log2err = 0.1,
        ES = -0.5,
        NES = -1,
        size = 2
      )
    )
    mock_fgsea_frame(pathways[names(pathways) == "A"], values)
  }
  testthat::expect_error(
    support_env$run_fgsea_collection(
      stats,
      collection,
      "fixture",
      1L,
      10L,
      10L,
      4L,
      100L,
      10L,
      incomplete_runner
    ),
    "incomplete or duplicate pathway keys"
  )
})
