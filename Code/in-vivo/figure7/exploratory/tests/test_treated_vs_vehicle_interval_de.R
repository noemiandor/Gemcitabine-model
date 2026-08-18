locate_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, ".git")) &&
        file.exists(file.path(current, "Manager.sh"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Cannot locate repository root", call. = FALSE)
    current <- parent
  }
}

configured_root <- Sys.getenv("GEMCITABINE_REPO_ROOT", "")
repo_root <- if (nzchar(configured_root)) {
  normalizePath(configured_root, mustWork = TRUE)
} else {
  locate_repo_root()
}
analysis_env <- new.env(parent = globalenv())
sys.source(
  file.path(
    repo_root,
    "Code/in-vivo/figure7/exploratory/treated_vs_vehicle_interval_de.R"
  ),
  envir = analysis_env
)

testthat::test_that("the exact reviewed interval includes both endpoints", {
  meta <- data.frame(
    cell_id = paste0("cell_", seq_len(6L)),
    sample_id = rep(c("mouse_a", "mouse_b"), each = 3L),
    initial_ploidy = rep(c("2N", "4N"), each = 3L),
    gemcitabine_dose_mg_per_kg = rep(c(0, 30), each = 3L),
    pseudotime = c(0.295999, 0.296, 0.300, 0.480, 0.486, 0.486001),
    stringsAsFactors = FALSE
  )
  interval <- list(
    start = 0.296,
    end = 0.486,
    include_start = TRUE,
    include_end = TRUE
  )
  selected <- analysis_env$select_interval_cells(meta, interval)
  testthat::expect_identical(
    selected$cell_id,
    c("cell_2", "cell_3", "cell_4", "cell_5")
  )
})

testthat::test_that("pseudobulk aggregation produces one library per mouse", {
  meta <- data.frame(
    cell_id = paste0("cell_", seq_len(6L)),
    sample_id = c("mouse_a", "mouse_a", "mouse_a", "mouse_b", "mouse_b", "mouse_b"),
    initial_ploidy = c(rep("2N", 3L), rep("4N", 3L)),
    gemcitabine_dose_mg_per_kg = c(rep(0, 3L), rep(30, 3L)),
    pseudotime = seq(0.31, 0.46, length.out = 6L),
    stringsAsFactors = FALSE
  )
  counts <- Matrix::Matrix(
    matrix(
      seq_len(18L),
      nrow = 3L,
      dimnames = list(paste0("GRCh38-gene_", seq_len(3L)), meta$cell_id)
    ),
    sparse = TRUE
  )
  result <- analysis_env$construct_mouse_pseudobulk(meta, counts)
  testthat::expect_identical(dim(result$counts), c(3L, 2L))
  testthat::expect_identical(colnames(result$counts), c("mouse_a", "mouse_b"))
  testthat::expect_equal(
    as.numeric(result$counts[, "mouse_a"]),
    as.numeric(Matrix::rowSums(counts[, seq_len(3L), drop = FALSE]))
  )
  testthat::expect_equal(
    as.numeric(result$counts[, "mouse_b"]),
    as.numeric(Matrix::rowSums(counts[, 4:6, drop = FALSE]))
  )
  testthat::expect_identical(result$metadata$n_cells, c(3L, 3L))
})

testthat::test_that("mixed mouse-level treatment metadata fails closed", {
  meta <- data.frame(
    cell_id = c("cell_a", "cell_b"),
    sample_id = c("mouse_a", "mouse_a"),
    initial_ploidy = c("2N", "2N"),
    gemcitabine_dose_mg_per_kg = c(0, 30),
    pseudotime = c(0.35, 0.36),
    stringsAsFactors = FALSE
  )
  testthat::expect_error(
    analysis_env$summarize_interval_mice(meta),
    "exactly one nonmissing value per mouse"
  )
})

testthat::test_that("the primary contrast is equal-dose treated minus vehicle", {
  sample_id <- paste0("mouse_", seq_len(16L))
  mouse_meta <- data.frame(
    sample_id = sample_id,
    initial_ploidy = rep(c("2N", "4N"), times = c(8L, 8L)),
    dose_mg = c(0, 0, 0, 0, 30, 30, 120, 120, 0, 0, 0, 0, 30, 30, 120, 120),
    treatment = c(rep("vehicle", 4L), rep("treated", 4L), rep("vehicle", 4L), rep("treated", 4L)),
    n_cells = rep(10L, 16L),
    mean_pseudotime = seq(0.36, 0.42, length.out = 16L),
    stringsAsFactors = FALSE,
    row.names = sample_id
  )
  design <- analysis_env$prepare_interval_design(mouse_meta)
  contrasts <- analysis_env$treatment_contrasts(design$design)
  primary <- contrasts$treated_equal_dose_minus_vehicle
  testthat::expect_identical(design$rank, ncol(design$design))
  testthat::expect_equal(primary[["dose_groupvehicle"]], -1)
  testthat::expect_equal(primary[["dose_groupdose_30"]], 0.5)
  testthat::expect_equal(primary[["dose_groupdose_120"]], 0.5)
  nuisance <- setdiff(names(primary), c(
    "dose_groupvehicle", "dose_groupdose_30", "dose_groupdose_120"
  ))
  testthat::expect_true(all(primary[nuisance] == 0))
})

testthat::test_that("the formal interaction is the 4N-minus-2N treatment contrast", {
  sample_id <- paste0("mouse_", seq_len(16L))
  mouse_meta <- data.frame(
    sample_id = sample_id,
    initial_ploidy = rep(c("2N", "4N"), each = 8L),
    dose_mg = rep(c(0, 0, 0, 0, 30, 30, 120, 120), 2L),
    treatment = rep(c(rep("vehicle", 4L), rep("treated", 4L)), 2L),
    n_cells = rep(10L, 16L),
    mean_pseudotime = seq(0.35, 0.43, length.out = 16L),
    stringsAsFactors = FALSE,
    row.names = sample_id
  )
  design <- analysis_env$prepare_treatment_origin_interaction_design(mouse_meta)
  contrasts <- analysis_env$treatment_origin_interaction_contrasts(design$design)
  interaction <- contrasts$treatment_by_origin_4N_minus_2N
  coefficient <- function(origin, dose) {
    paste0("initial_ploidy_factor", origin, ":dose_group", dose)
  }
  expected <- setNames(
    c(1, -0.5, -0.5, -1, 0.5, 0.5),
    c(
      coefficient("2N", c("vehicle", "dose_30", "dose_120")),
      coefficient("4N", c("vehicle", "dose_30", "dose_120"))
    )
  )
  testthat::expect_identical(design$rank, ncol(design$design))
  testthat::expect_identical(ncol(design$design), 7L)
  testthat::expect_equal(unname(interaction[names(expected)]), unname(expected))
  testthat::expect_equal(interaction[["mean_pseudotime_z"]], 0)

  group_means <- setNames(rep(0, ncol(design$design)), colnames(design$design))
  group_means[coefficient("2N", c("vehicle", "dose_30", "dose_120"))] <- c(10, 12, 14)
  group_means[coefficient("4N", c("vehicle", "dose_30", "dose_120"))] <- c(20, 26, 28)
  testthat::expect_equal(sum(group_means * interaction), 4)
})

testthat::test_that("interaction signs identify the origin with the stronger treatment effect", {
  set.seed(713L)
  sample_id <- paste0("mouse_", seq_len(16L))
  origin <- rep(c("2N", "4N"), each = 8L)
  dose <- rep(c(0, 0, 0, 0, 30, 30, 120, 120), 2L)
  mouse_meta <- data.frame(
    sample_id = sample_id,
    initial_ploidy = origin,
    dose_mg = dose,
    treatment = ifelse(dose == 0, "vehicle", "treated"),
    n_cells = rep(10L, 16L),
    mean_pseudotime = rep(seq(0.35, 0.42, length.out = 8L), 2L),
    stringsAsFactors = FALSE,
    row.names = sample_id
  )
  counts <- matrix(
    stats::rpois(120L * 16L, lambda = 100),
    nrow = 120L,
    dimnames = list(paste0("GRCh38-gene_", seq_len(120L)), sample_id)
  )
  counts[1L, origin == "4N" & dose > 0] <- 700
  counts[2L, origin == "2N" & dose > 0] <- 700
  design <- analysis_env$prepare_treatment_origin_interaction_design(mouse_meta)
  contrasts <- analysis_env$treatment_origin_interaction_contrasts(design$design)
  model <- analysis_env$fit_interval_voom(
    Matrix::Matrix(counts, sparse = TRUE),
    design,
    cpm_threshold = 0,
    minimum_mice = 4L
  )
  table <- analysis_env$extract_contrast_table(
    model,
    contrasts$treatment_by_origin_4N_minus_2N,
    "treatment_by_origin_4N_minus_2N",
    function(feature) sub("^GRCh38-", "", feature),
    positive_direction = "treatment_effect_more_positive_in_4N",
    negative_direction = "treatment_effect_more_positive_in_2N"
  )
  four_n <- table[table$feature_id == "GRCh38-gene_1", , drop = FALSE]
  two_n <- table[table$feature_id == "GRCh38-gene_2", , drop = FALSE]
  testthat::expect_gt(four_n$log2_fold_change, 0)
  testthat::expect_identical(
    four_n$direction,
    "treatment_effect_more_positive_in_4N"
  )
  testthat::expect_lt(two_n$log2_fold_change, 0)
  testthat::expect_identical(
    two_n$direction,
    "treatment_effect_more_positive_in_2N"
  )
})

testthat::test_that("origin-stratified designs omit the constant origin term", {
  sample_id <- paste0("mouse_2N_", seq_len(8L))
  mouse_meta <- data.frame(
    sample_id = sample_id,
    initial_ploidy = rep("2N", 8L),
    dose_mg = c(0, 0, 0, 0, 30, 30, 120, 120),
    treatment = c(rep("vehicle", 4L), rep("treated", 4L)),
    n_cells = rep(10L, 8L),
    mean_pseudotime = seq(0.35, 0.42, length.out = 8L),
    stringsAsFactors = FALSE,
    row.names = sample_id
  )
  testthat::expect_invisible(
    analysis_env$validate_origin_stratum(mouse_meta, "2N")
  )
  design <- analysis_env$prepare_interval_design(
    mouse_meta,
    adjust_mean_pseudotime = TRUE,
    adjust_initial_ploidy = FALSE
  )
  contrasts <- analysis_env$treatment_contrasts(design$design)
  primary <- contrasts$treated_equal_dose_minus_vehicle
  testthat::expect_identical(design$rank, ncol(design$design))
  testthat::expect_identical(ncol(design$design), 4L)
  testthat::expect_false(any(grepl("initial_ploidy", colnames(design$design))))
  testthat::expect_equal(primary[["dose_groupvehicle"]], -1)
  testthat::expect_equal(primary[["dose_groupdose_30"]], 0.5)
  testthat::expect_equal(primary[["dose_groupdose_120"]], 0.5)
  testthat::expect_equal(primary[["mean_pseudotime_z"]], 0)

  invalid <- mouse_meta[-1L, , drop = FALSE]
  testthat::expect_error(
    analysis_env$validate_origin_stratum(invalid, "2N"),
    "does not contain eight mice"
  )
})

testthat::test_that("positive fold change means higher expression in treated mice", {
  set.seed(712L)
  sample_id <- paste0("mouse_", seq_len(16L))
  dose <- c(0, 0, 0, 0, 30, 30, 120, 120, 0, 0, 0, 0, 30, 30, 120, 120)
  mouse_meta <- data.frame(
    sample_id = sample_id,
    initial_ploidy = rep(c("2N", "4N"), each = 8L),
    dose_mg = dose,
    treatment = ifelse(dose == 0, "vehicle", "treated"),
    n_cells = rep(10L, 16L),
    mean_pseudotime = rep(seq(0.35, 0.42, length.out = 8L), 2L),
    stringsAsFactors = FALSE,
    row.names = sample_id
  )
  counts <- matrix(
    stats::rpois(120L * 16L, lambda = 100),
    nrow = 120L,
    dimnames = list(paste0("GRCh38-gene_", seq_len(120L)), sample_id)
  )
  treated <- dose > 0
  counts[1L, treated] <- stats::rpois(sum(treated), lambda = 500)
  counts[1L, !treated] <- stats::rpois(sum(!treated), lambda = 30)
  counts[2L, treated] <- stats::rpois(sum(treated), lambda = 25)
  counts[2L, !treated] <- stats::rpois(sum(!treated), lambda = 450)
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  design <- analysis_env$prepare_interval_design(mouse_meta)
  contrasts <- analysis_env$treatment_contrasts(design$design)
  model <- analysis_env$fit_interval_voom(
    counts,
    design,
    cpm_threshold = 0,
    minimum_mice = 4L
  )
  table <- analysis_env$extract_contrast_table(
    model,
    contrasts$treated_equal_dose_minus_vehicle,
    "treated_equal_dose_minus_vehicle",
    function(feature) sub("^GRCh38-", "", feature)
  )
  gene_up <- table[table$feature_id == "GRCh38-gene_1", , drop = FALSE]
  gene_down <- table[table$feature_id == "GRCh38-gene_2", , drop = FALSE]
  testthat::expect_gt(gene_up$log2_fold_change, 0)
  testthat::expect_identical(gene_up$direction, "higher_in_treated")
  testthat::expect_lt(gene_down$log2_fold_change, 0)
  testthat::expect_identical(gene_down$direction, "lower_in_treated")
})

testthat::test_that("pathway ranking reuses moderated t and resolves duplicate symbols", {
  support <- analysis_env$load_figure7_support(repo_root)
  table <- data.frame(
    feature_id = c("GRCh38-a1", "GRCh38-a2", "GRCh38-b"),
    gene_symbol = c("GENEA", "GENEA", "GENEB"),
    contrast_id = rep("treated_equal_dose_minus_vehicle", 3L),
    moderated_t = c(2, -5, 1),
    p_value = c(0.04, 0.001, 0.2),
    fdr = c(0.08, 0.01, 0.3),
    stringsAsFactors = FALSE
  )
  ranking <- analysis_env$prepare_pathway_ranking(table, support)
  testthat::expect_identical(names(ranking$stats), c("GENEB", "GENEA"))
  testthat::expect_equal(unname(ranking$stats), c(1, -5))
  retained <- ranking$symbol_resolution[
    ranking$symbol_resolution$retained_for_gsea,
    ,
    drop = FALSE
  ]
  testthat::expect_identical(retained$gene, c("GRCh38-a2", "GRCh38-b"))
})

testthat::test_that("pathway display caps the reviewed selector at three per direction", {
  support <- analysis_env$load_figure7_support(repo_root)
  config <- support$read_config(file.path(
    repo_root,
    "Code/in-vivo/figure7/figure7_config.yaml"
  ))
  collections <- as.character(unlist(config$state_pathways$collections))
  gsea <- do.call(rbind, lapply(collections, function(collection) {
    data.frame(
      collection = collection,
      collection_label = collection,
      pathway = paste0(gsub("[^A-Za-z]", "", collection), "_", seq_len(13L)),
      pathway_label = paste("Pathway", seq_len(13L)),
      padj = c(rep(0.01, 12L), 0.2),
      NES = c(seq(2.5, 1.5, length.out = 6L), seq(-1.5, -2.5, length.out = 6L), 3),
      stringsAsFactors = FALSE
    )
  }))
  selected <- analysis_env$select_pathway_results(gsea, support, config)
  counts <- table(selected$collection, selected$selected_direction)
  testthat::expect_true(all(counts[, "positive"] == 3L))
  testthat::expect_true(all(counts[, "negative"] == 3L))
  testthat::expect_true(all(selected$padj <= 0.05))
  testthat::expect_false(any(grepl("_13$", selected$pathway)))
})

testthat::test_that("formal interaction panel selects the ten NES extremes per sign", {
  collections <- c("H", "C2:CP:REACTOME", "C5:GO:BP")
  make_rows <- function(pathways, nes, padj) {
    data.frame(
      collection = rep(collections, length.out = length(pathways)),
      collection_label = rep(c("hallmark", "reactome", "go_bp"), length.out = length(pathways)),
      pathway = pathways,
      pathway_label = gsub("_", " ", pathways, fixed = TRUE),
      NES = nes,
      padj = padj,
      stringsAsFactors = FALSE
    )
  }
  shared_2n <- make_rows(
    paste0("shared_", seq_len(5L)),
    c(-2.2, -1.8, -1.5, 1.6, 2.0),
    seq(0.001, 0.005, length.out = 5L)
  )
  shared_4n <- shared_2n
  shared_4n$NES <- shared_4n$NES * 0.9
  interaction <- make_rows(
    paste0("formal_interaction_", seq_len(24L)),
    c(seq(-1.3, -2.4, length.out = 12L), seq(1.3, 2.4, length.out = 12L)),
    rep(c(0.4, 0.03, 0.01, 0.2), 6L)
  )
  interaction_2n <- interaction
  interaction_2n$NES <- c(rep(1, 12L), rep(-1, 12L))
  interaction_2n$padj <- 0.5
  interaction_4n <- interaction
  interaction_4n$NES <- 0
  interaction_4n$padj <- 0.5
  comparison <- analysis_env$prepare_origin_interaction_pathways(
    shared_2n,
    shared_4n,
    interaction,
    origin_effect_gsea_2n = interaction_2n,
    origin_effect_gsea_4n = interaction_4n,
    fdr_threshold = 0.05,
    maximum_per_direction = 10L
  )
  testthat::expect_identical(nrow(comparison$shared), 5L)
  testthat::expect_identical(nrow(comparison$selected_interactions), 20L)
  testthat::expect_true(all(grepl(
    "^formal_interaction_",
    comparison$selected_interactions$pathway
  )))
  testthat::expect_setequal(
    unique(comparison$selected_interactions$panel_id),
    c("negative_interaction", "positive_interaction")
  )
  testthat::expect_true(all(comparison$selected_2N$NES < 0))
  testthat::expect_true(all(comparison$selected_4N$NES > 0))
  testthat::expect_true(all(comparison$selected_2N$selected_direction == "negative"))
  testthat::expect_true(all(comparison$selected_4N$selected_direction == "positive"))
  testthat::expect_equal(comparison$selected_2N$NES, sort(interaction$NES)[1:10])
  testthat::expect_equal(
    comparison$selected_4N$NES,
    sort(interaction$NES, decreasing = TRUE)[1:10]
  )
  testthat::expect_true(all(
    comparison$selected_interactions$origin_effect_nes_order_concordant
  ))
  testthat::expect_true(all(comparison$selected_2N$NES_2N == 1))
  testthat::expect_true(all(comparison$selected_2N$NES_4N == 0))
  testthat::expect_true(all(comparison$selected_4N$NES_2N == -1))
  testthat::expect_true(all(comparison$selected_4N$NES_4N == 0))
  testthat::expect_identical(
    comparison$selected_2N$selected_rank_within_direction,
    seq_len(10L)
  )
  testthat::expect_identical(
    comparison$selected_4N$selected_rank_within_direction,
    seq_len(10L)
  )
})
