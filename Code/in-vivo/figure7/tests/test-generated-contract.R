testthat::test_that(
  "generated panel 7F selects only FDR-significant pathways without backfill",
  {
    config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
    collections <- as.character(unlist(config$state_pathways$collections))
    rows <- list()
    for (collection in collections) {
      positive_padj <- if (identical(collection, "H")) {
        c(0.05, 0.0500001, 0.20, 0.60, 0.90)
      } else {
        c(0.001, 0.002, 0.003, 0.004, 0.005)
      }
      negative_padj <- c(0.0015, 0.0025, 0.0035, 0.0045, 0.0055)
      rows[[length(rows) + 1L]] <- data.frame(
        collection_id = collection,
        pathway_id = paste0(collection, "_positive_", seq_along(positive_padj)),
        padj = positive_padj,
        NES = seq(5, 1, length.out = length(positive_padj)),
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- data.frame(
        collection_id = collection,
        pathway_id = paste0(collection, "_negative_", seq_along(negative_padj)),
        padj = negative_padj,
        NES = -seq(5, 1, length.out = length(negative_padj)),
        stringsAsFactors = FALSE
      )
    }
    complete <- do.call(rbind, rows)
    selected <- figure7_select_generated_pathways(complete, config)
    counts <- table(factor(
      selected$collection_id,
      levels = collections
    ))

    testthat::expect_equal(nrow(selected), 21L)
    testthat::expect_identical(as.integer(counts), c(5L, 8L, 8L))
    testthat::expect_true(all(
      selected$padj <= figure7_generated_pathway_fdr_threshold()
    ))
    testthat::expect_true("H_positive_1" %in% selected$pathway_id)
    testthat::expect_false(any(
      paste0("H_positive_", 2:5) %in% selected$pathway_id
    ))
    testthat::expect_false(any(grepl("_positive_5$", selected$pathway_id)))
    testthat::expect_false(any(grepl("_negative_5$", selected$pathway_id)))
    testthat::expect_equal(
      sum(
        selected$collection_id == "H" &
          selected$selected_direction == "positive"
      ),
      1L
    )

    no_go <- complete
    no_go$padj[no_go$collection_id == "C5:GO:BP"] <- 0.5
    testthat::expect_error(
      figure7_select_generated_pathways(no_go, config),
      "no FDR-significant pathway"
    )
    testthat::expect_error(
      figure7_select_generated_pathways(
        rbind(complete, complete[1L, , drop = FALSE]),
        config
      ),
      "unique keys"
    )
    invalid_cap <- config
    invalid_cap$state_pathways$top_positive_per_collection <- 3L
    testthat::expect_error(
      figure7_select_generated_pathways(complete, invalid_cap),
      "exactly four candidates per sign"
    )
  }
)

testthat::test_that(
  "state-pathway v3 adjusts for injected initial ploidy and excludes endpoint CN score",
  {
    support_env <- new.env(parent = globalenv())
    sys.source(
      file.path(
        module_dir,
        "generate_pseudotime_state_pathways_support.R"
      ),
      envir = support_env
    )
    metadata <- data.frame(
      bin_midpoint = seq(0.05, 0.95, length.out = 12L),
      dose_mg = rep(c(0, 30, 120), each = 4L),
      initial_ploidy = rep(c("2N", "4N"), 6L),
      stringsAsFactors = FALSE
    )
    spec <- support_env$model_spec_initial_ploidy()
    design <- support_env$prepare_design(
      metadata,
      spline_df = 5L,
      model_spec = spec,
      include_dose = TRUE
    )

    testthat::expect_identical(
      spec$model_id,
      "initial_ploidy_adjusted_grch_human_only_v3"
    )
    testthat::expect_identical(spec$covariate_mode, "initial_ploidy")
    testthat::expect_identical(spec$covariate_terms, "initial_ploidy_factor")
    testthat::expect_true("initial_ploidy_factor4N" %in% colnames(design$design))
    testthat::expect_false(any(grepl(
      "ETP|endpoint|copy.number|cn_score",
      colnames(design$design),
      ignore.case = TRUE
    )))
    testthat::expect_identical(
      as.character(unlist(
        figure7_read_config(
          file.path(module_dir, "figure7_config.yaml")
        )$state_pathways$nuisance_terms
      )),
      c("dose_mg_factor", "initial_ploidy_factor")
    )
  }
)

testthat::test_that("generated panel 7F renderer supports 5/8/8 pathway facets", {
  testthat::skip_if_not_installed("ggplot2")
  config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
  collections <- as.character(unlist(config$state_pathways$collections))
  collection_labels <- c(
    "H" = "Hallmark",
    "C2:CP:REACTOME" = "Reactome",
    "C5:GO:BP" = "GO biological process"
  )
  counts <- c(5L, 8L, 8L)
  pathway_metadata <- do.call(rbind, lapply(
    seq_along(collections),
    function(collection_index) {
      collection <- collections[[collection_index]]
      data.frame(
        collection_id = collection,
        collection_label = unname(collection_labels[[collection]]),
        collection_display_order = collection_index,
        pathway_id = paste0(collection, "_pathway_", seq_len(counts[[collection_index]])),
        pathway_label = paste("Pathway", seq_len(counts[[collection_index]])),
        stringsAsFactors = FALSE
      )
    }
  ))
  pathway_metadata$pathway_display_order <- seq_len(nrow(pathway_metadata))
  grid <- seq(0, 1, length.out = 501L)
  activity <- do.call(rbind, lapply(
    seq_len(nrow(pathway_metadata)),
    function(index) {
      data.frame(
        pathway_metadata[
          rep(index, length(grid)),
          ,
          drop = FALSE
        ],
        pseudotime = grid,
        standardized_activity = sin(2 * pi * grid) + index / 100,
        row.names = NULL
      )
    }
  ))

  built <- ggplot2::ggplot_build(figure7_panel_f_plot(activity, config))
  testthat::expect_equal(nrow(built$data[[1L]]), 21L * 501L)
  testthat::expect_equal(
    length(unique(built$layout$layout$PANEL)),
    3L
  )
})
