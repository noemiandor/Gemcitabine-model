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
