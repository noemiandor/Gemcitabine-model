#!/usr/bin/env Rscript

# Exploratory candidate for main Figure 7I.
#
# This analysis deliberately does not alter the reviewed Figure 7 workflow or
# its frozen outputs. It reuses that workflow's exact CellCycle metadata,
# fail-closed GRCh38-only count boundary, and edgeR/limma-voom framework, then
# asks a different question: within the prespecified pseudotime interval, which
# genes and pathways differ between gemcitabine-treated and vehicle tumors?

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

exploratory_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) return(NA_character_)
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

exploratory_repo_root <- function(script = exploratory_script_path()) {
  if (is.na(script)) return(normalizePath(getwd(), mustWork = TRUE))
  normalizePath(file.path(dirname(script), "..", "..", "..", ".."), mustWork = TRUE)
}

parse_scalar <- function(value, template) {
  if (is.logical(template)) {
    normalized <- tolower(as.character(value))
    if (!normalized %in% c("true", "false", "t", "f", "1", "0", "yes", "no", "y", "n")) {
      stop("Expected a logical value; received: ", value, call. = FALSE)
    }
    return(normalized %in% c("true", "t", "1", "yes", "y"))
  }
  if (is.integer(template)) return(as.integer(value))
  if (is.numeric(template)) return(as.numeric(value))
  as.character(value)
}

parse_cli_args <- function(argv, defaults) {
  out <- defaults
  index <- 1L
  while (index <= length(argv)) {
    token <- argv[[index]]
    if (!grepl("^--", token)) stop("Unexpected positional argument: ", token, call. = FALSE)
    key_value <- sub("^--", "", token)
    if (grepl("=", key_value, fixed = TRUE)) {
      key <- sub("=.*$", "", key_value)
      value <- sub("^[^=]*=", "", key_value)
      index <- index + 1L
    } else {
      key <- key_value
      if (index < length(argv) && !grepl("^--", argv[[index + 1L]])) {
        value <- argv[[index + 1L]]
        index <- index + 2L
      } else {
        value <- "TRUE"
        index <- index + 1L
      }
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    if (!key %in% names(out)) stop("Unknown argument: --", key, call. = FALSE)
    out[[key]] <- parse_scalar(value, out[[key]])
  }
  out
}

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript Code/in-vivo/figure7/exploratory/treated_vs_vehicle_interval_de.R [options]",
    "",
    "Options:",
    "  --cell_metadata PATH       Exact 2,881-cell CellCycle metadata table",
    "  --seurat_rds PATH          Reviewed Seurat object containing raw RNA counts",
    "  --config PATH              Figure 7 YAML containing the reviewed interval",
    "  --output_root PATH         Isolated exploratory output directory",
    "  --cpm_threshold NUMBER     Retain genes with CPM above this value [1]",
    "  --minimum_mice INTEGER     Required mice above the CPM threshold [4]",
    "  --fdr_threshold NUMBER     Plot significance threshold [0.05]",
    "  --lfc_threshold NUMBER     Plot absolute log2-FC threshold [0.5]",
    "  --label_genes INTEGER      Maximum volcano labels [12]",
    "  --magnitude_bootstrap_replicates INTEGER",
    "                             Whole-mouse wild-bootstrap draws [4999]",
    "  --magnitude_bootstrap_seed INTEGER",
    "                             Seed for the global magnitude comparison [20260818]",
    "  --global_magnitude_only BOOL",
    "                             Skip gene-set loading and write only the magnitude audit [FALSE]",
    "  --verify_input_hashes BOOL Verify inputs against the reviewed config [TRUE]",
    "  --overwrite BOOL           Replace this exact output directory [FALSE]",
    sep = "\n"
  ), "\n")
}

is_absolute_path <- function(path) {
  grepl("^/", path) || grepl("^[A-Za-z]:[/\\\\]", path)
}

resolve_path <- function(path, repo_root, must_work = FALSE) {
  resolved <- if (is_absolute_path(path)) path else file.path(repo_root, path)
  normalizePath(resolved, mustWork = must_work)
}

required_packages <- function() {
  c(
    "Seurat", "SeuratObject", "Matrix", "yaml", "digest", "edgeR",
    "limma", "fgsea", "BiocParallel", "msigdbr", "readr", "ggplot2",
    "ggrepel"
  )
}

global_magnitude_required_packages <- function() {
  c(
    "Seurat", "SeuratObject", "Matrix", "yaml", "digest", "edgeR",
    "limma", "readr"
  )
}

check_packages <- function(packages = required_packages()) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

load_figure7_support <- function(repo_root) {
  support <- new.env(parent = globalenv())
  module_dir <- file.path(repo_root, "Code", "in-vivo", "figure7")
  sys.source(
    file.path(module_dir, "generate_pseudotime_state_pathways_support.R"),
    envir = support
  )
  sys.source(
    file.path(module_dir, "src", "feature_species_policy.R"),
    envir = support
  )
  sys.source(
    file.path(module_dir, "src", "common_io.R"),
    envir = support
  )
  support
}

assert_unique_value <- function(x, label) {
  observed <- unique(x[!is.na(x)])
  if (length(observed) != 1L) {
    stop(label, " must have exactly one nonmissing value per mouse", call. = FALSE)
  }
  observed[[1L]]
}

select_interval_cells <- function(meta, interval) {
  required <- c(
    "cell_id", "sample_id", "initial_ploidy",
    "gemcitabine_dose_mg_per_kg", "pseudotime"
  )
  missing <- setdiff(required, names(meta))
  if (length(missing)) {
    stop("Cell metadata missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(meta$cell_id)) stop("Cell IDs must be unique", call. = FALSE)
  pseudotime <- suppressWarnings(as.numeric(meta$pseudotime))
  if (any(!is.finite(pseudotime))) stop("Pseudotime must be finite", call. = FALSE)
  left <- if (isTRUE(interval$include_start)) {
    pseudotime >= interval$start
  } else {
    pseudotime > interval$start
  }
  right <- if (isTRUE(interval$include_end)) {
    pseudotime <= interval$end
  } else {
    pseudotime < interval$end
  }
  selected <- meta[left & right, , drop = FALSE]
  selected$pseudotime <- pseudotime[left & right]
  if (!nrow(selected)) stop("No cells fall inside the configured interval", call. = FALSE)
  selected
}

summarize_interval_mice <- function(meta) {
  split_meta <- split(meta, as.character(meta$sample_id))
  rows <- lapply(names(split_meta), function(sample_id) {
    x <- split_meta[[sample_id]]
    dose_mg <- as.numeric(assert_unique_value(
      x$gemcitabine_dose_mg_per_kg,
      paste0("Dose for ", sample_id)
    ))
    initial_ploidy <- as.character(assert_unique_value(
      x$initial_ploidy,
      paste0("Injected origin for ", sample_id)
    ))
    data.frame(
      sample_id = sample_id,
      initial_ploidy = initial_ploidy,
      dose_mg = dose_mg,
      treatment = ifelse(dose_mg == 0, "vehicle", "treated"),
      n_cells = nrow(x),
      mean_pseudotime = mean(x$pseudotime),
      median_pseudotime = stats::median(x$pseudotime),
      minimum_pseudotime = min(x$pseudotime),
      maximum_pseudotime = max(x$pseudotime),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$dose_mg, out$initial_ploidy, out$sample_id), , drop = FALSE]
  rownames(out) <- out$sample_id
  out
}

validate_interval_contract <- function(all_meta, interval_meta, mouse_meta) {
  expected_dose_by_origin <- matrix(
    c(4L, 2L, 2L, 4L, 2L, 2L),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("2N", "4N"), c("0", "30", "120"))
  )
  observed <- table(
    factor(mouse_meta$initial_ploidy, levels = rownames(expected_dose_by_origin)),
    factor(as.character(mouse_meta$dose_mg), levels = colnames(expected_dose_by_origin))
  )
  failures <- character()
  if (nrow(all_meta) != 2881L) failures <- c(failures, "full CellCycle universe is not 2,881 cells")
  if (nrow(interval_meta) != 417L) failures <- c(failures, "interval universe is not 417 cells")
  if (nrow(mouse_meta) != 16L) failures <- c(failures, "interval does not contain all 16 mice")
  if (!all(as.integer(observed) == as.integer(expected_dose_by_origin))) {
    failures <- c(failures, "origin-by-dose mouse counts differ from the reviewed 4/2/2 per-origin design")
  }
  if (!setequal(unique(mouse_meta$initial_ploidy), c("2N", "4N"))) {
    failures <- c(failures, "both injected-origin groups are not represented")
  }
  if (!setequal(unique(mouse_meta$dose_mg), c(0, 30, 120))) {
    failures <- c(failures, "vehicle, 30 mg/kg, and 120 mg/kg are not all represented")
  }
  if (length(failures)) {
    stop("Exploratory interval contract failed: ", paste(failures, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

construct_mouse_pseudobulk <- function(meta, counts) {
  if (is.null(colnames(counts)) || is.null(rownames(counts))) {
    stop("Count matrix requires feature and cell names", call. = FALSE)
  }
  if (anyDuplicated(meta$cell_id) || anyDuplicated(colnames(counts))) {
    stop("Cell IDs must be unique before pseudobulk aggregation", call. = FALSE)
  }
  missing_cells <- setdiff(meta$cell_id, colnames(counts))
  if (length(missing_cells)) {
    stop(length(missing_cells), " interval cell(s) are absent from the count matrix", call. = FALSE)
  }
  counts <- counts[, meta$cell_id, drop = FALSE]
  mouse_meta <- summarize_interval_mice(meta)
  sample_factor <- factor(meta$sample_id, levels = rownames(mouse_meta))
  aggregation <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
  colnames(aggregation) <- levels(sample_factor)
  pseudobulk <- counts %*% aggregation
  colnames(pseudobulk) <- levels(sample_factor)
  rownames(pseudobulk) <- rownames(counts)
  mouse_meta$library_size <- as.numeric(Matrix::colSums(pseudobulk)[mouse_meta$sample_id])
  if (any(!is.finite(mouse_meta$library_size)) || any(mouse_meta$library_size <= 0)) {
    stop("Every interval mouse must have a positive pseudobulk library size", call. = FALSE)
  }
  list(counts = pseudobulk, metadata = mouse_meta)
}

prepare_interval_design <- function(
  mouse_meta,
  adjust_mean_pseudotime = TRUE,
  adjust_initial_ploidy = TRUE
) {
  data <- as.data.frame(mouse_meta, stringsAsFactors = FALSE)
  data$dose_group <- factor(
    as.character(data$dose_mg),
    levels = c("0", "30", "120"),
    labels = c("vehicle", "dose_30", "dose_120")
  )
  data$initial_ploidy_factor <- factor(data$initial_ploidy, levels = c("2N", "4N"))
  pseudotime_sd <- stats::sd(data$mean_pseudotime)
  if (!is.finite(pseudotime_sd) || pseudotime_sd <= 0) {
    stop("Mouse mean pseudotime has no finite variation", call. = FALSE)
  }
  data$mean_pseudotime_z <- (
    data$mean_pseudotime - mean(data$mean_pseudotime)
  ) / pseudotime_sd
  if (anyNA(data$dose_group) || anyNA(data$initial_ploidy_factor)) {
    stop("Design contains an unsupported dose or injected-origin value", call. = FALSE)
  }
  design_terms <- "dose_group"
  if (isTRUE(adjust_initial_ploidy)) {
    design_terms <- c(design_terms, "initial_ploidy_factor")
  }
  if (isTRUE(adjust_mean_pseudotime)) {
    design_terms <- c(design_terms, "mean_pseudotime_z")
  }
  formula <- stats::reformulate(design_terms, intercept = FALSE)
  design <- stats::model.matrix(formula, data = data)
  rownames(design) <- data$sample_id
  rank <- qr(design)$rank
  if (rank != ncol(design)) {
    stop("Interval differential-expression design is rank deficient", call. = FALSE)
  }
  list(
    data = data,
    design = design,
    formula = formula,
    rank = rank,
    adjust_mean_pseudotime = isTRUE(adjust_mean_pseudotime),
    adjust_initial_ploidy = isTRUE(adjust_initial_ploidy)
  )
}

treatment_contrasts <- function(design) {
  required <- c("dose_groupvehicle", "dose_groupdose_30", "dose_groupdose_120")
  missing <- setdiff(required, colnames(design))
  if (length(missing)) {
    stop("Design is missing dose coefficients: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  make_contrast <- function(weights) {
    out <- setNames(rep(0, ncol(design)), colnames(design))
    out[names(weights)] <- weights
    out
  }
  list(
    treated_equal_dose_minus_vehicle = make_contrast(c(
      dose_groupvehicle = -1,
      dose_groupdose_30 = 0.5,
      dose_groupdose_120 = 0.5
    )),
    dose_30_minus_vehicle = make_contrast(c(
      dose_groupvehicle = -1,
      dose_groupdose_30 = 1
    )),
    dose_120_minus_vehicle = make_contrast(c(
      dose_groupvehicle = -1,
      dose_groupdose_120 = 1
    ))
  )
}

prepare_treatment_origin_interaction_design <- function(
  mouse_meta,
  adjust_mean_pseudotime = TRUE
) {
  data <- as.data.frame(mouse_meta, stringsAsFactors = FALSE)
  data$dose_group <- factor(
    as.character(data$dose_mg),
    levels = c("0", "30", "120"),
    labels = c("vehicle", "dose_30", "dose_120")
  )
  data$initial_ploidy_factor <- factor(
    data$initial_ploidy,
    levels = c("2N", "4N")
  )
  pseudotime_sd <- stats::sd(data$mean_pseudotime)
  if (!is.finite(pseudotime_sd) || pseudotime_sd <= 0) {
    stop("Mouse mean pseudotime has no finite variation", call. = FALSE)
  }
  data$mean_pseudotime_z <- (
    data$mean_pseudotime - mean(data$mean_pseudotime)
  ) / pseudotime_sd
  if (anyNA(data$dose_group) || anyNA(data$initial_ploidy_factor)) {
    stop("Interaction design contains an unsupported dose or injected-origin value", call. = FALSE)
  }
  design_terms <- "initial_ploidy_factor:dose_group"
  if (isTRUE(adjust_mean_pseudotime)) {
    design_terms <- c(design_terms, "mean_pseudotime_z")
  }
  formula <- stats::reformulate(design_terms, intercept = FALSE)
  design <- stats::model.matrix(formula, data = data)
  rownames(design) <- data$sample_id
  rank <- qr(design)$rank
  if (rank != ncol(design)) {
    stop("Treatment-by-origin interaction design is rank deficient", call. = FALSE)
  }
  list(
    data = data,
    design = design,
    formula = formula,
    rank = rank,
    adjust_mean_pseudotime = isTRUE(adjust_mean_pseudotime),
    adjust_initial_ploidy = TRUE,
    treatment_by_origin_interaction = TRUE
  )
}

treatment_origin_interaction_contrasts <- function(design) {
  coefficient <- function(origin, dose) {
    paste0("initial_ploidy_factor", origin, ":dose_group", dose)
  }
  required <- unlist(lapply(c("2N", "4N"), function(origin) {
    coefficient(origin, c("vehicle", "dose_30", "dose_120"))
  }), use.names = FALSE)
  missing <- setdiff(required, colnames(design))
  if (length(missing)) {
    stop(
      "Interaction design is missing origin-by-dose coefficients: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  make_contrast <- function(weights) {
    out <- setNames(rep(0, ncol(design)), colnames(design))
    out[names(weights)] <- weights
    out
  }
  effect_weights <- function(origin) {
    setNames(
      c(-1, 0.5, 0.5),
      coefficient(origin, c("vehicle", "dose_30", "dose_120"))
    )
  }
  effect_2n <- make_contrast(effect_weights("2N"))
  effect_4n <- make_contrast(effect_weights("4N"))
  list(
    treated_equal_dose_minus_vehicle_2N = effect_2n,
    treated_equal_dose_minus_vehicle_4N = effect_4n,
    treatment_by_origin_4N_minus_2N = effect_4n - effect_2n
  )
}

fit_interval_voom <- function(
  pseudobulk_counts,
  design_bundle,
  cpm_threshold = 1,
  minimum_mice = 4L
) {
  design <- design_bundle$design
  if (!identical(colnames(pseudobulk_counts), rownames(design))) {
    stop("Pseudobulk columns and design rows are not identically ordered", call. = FALSE)
  }
  if (minimum_mice < 2L || minimum_mice > ncol(pseudobulk_counts)) {
    stop("minimum_mice is outside the valid sample range", call. = FALSE)
  }
  unfiltered <- edgeR::DGEList(counts = pseudobulk_counts)
  cpm <- edgeR::cpm(unfiltered)
  keep <- rowSums(cpm > cpm_threshold) >= minimum_mice
  if (sum(keep) < 10L) stop("Fewer than ten genes pass the expression filter", call. = FALSE)
  y <- edgeR::DGEList(counts = pseudobulk_counts[keep, , drop = FALSE])
  y <- edgeR::calcNormFactors(y, method = "TMM")
  voom <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(voom, design)
  list(
    fit = fit,
    voom = voom,
    dge = y,
    retained = keep,
    n_input_genes = nrow(pseudobulk_counts),
    n_retained_genes = sum(keep),
    cpm_threshold = cpm_threshold,
    minimum_mice = minimum_mice
  )
}

extract_contrast_table <- function(
  model,
  contrast,
  contrast_id,
  clean_symbols,
  positive_direction = "higher_in_treated",
  negative_direction = "lower_in_treated"
) {
  coefficients <- colnames(model$fit$coefficients)
  if (!identical(names(contrast), coefficients)) {
    stop("Contrast coefficients do not align with the fitted design", call. = FALSE)
  }
  contrast_matrix <- matrix(
    contrast,
    ncol = 1L,
    dimnames = list(coefficients, contrast_id)
  )
  fit <- limma::contrasts.fit(model$fit, contrasts = contrast_matrix)
  fit <- limma::eBayes(fit, robust = TRUE)
  table <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "none")
  table$feature_id <- rownames(table)
  table$gene_symbol <- clean_symbols(table$feature_id)
  out <- data.frame(
    feature_id = table$feature_id,
    gene_symbol = table$gene_symbol,
    contrast_id = contrast_id,
    log2_fold_change = table$logFC,
    average_log2_expression = table$AveExpr,
    moderated_t = table$t,
    p_value = table$P.Value,
    fdr = table$adj.P.Val,
    log_odds_differential_expression = table$B,
    direction = ifelse(
      table$logFC > 0,
      positive_direction,
      ifelse(table$logFC < 0, negative_direction, "no_change")
    ),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$p_value, -abs(out$moderated_t), out$feature_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

select_volcano_labels <- function(table, maximum_labels, fdr_threshold, lfc_threshold) {
  candidates <- table[
    is.finite(table$fdr) & table$fdr <= fdr_threshold &
      is.finite(table$log2_fold_change) & abs(table$log2_fold_change) >= lfc_threshold &
      !is.na(table$gene_symbol) & nzchar(table$gene_symbol),
    ,
    drop = FALSE
  ]
  if (!nrow(candidates) || maximum_labels <= 0L) return(candidates[FALSE, , drop = FALSE])
  candidates <- candidates[!duplicated(candidates$gene_symbol), , drop = FALSE]
  up <- candidates[candidates$log2_fold_change > 0, , drop = FALSE]
  down <- candidates[candidates$log2_fold_change < 0, , drop = FALSE]
  per_direction <- ceiling(maximum_labels / 2)
  labels <- rbind(head(up, per_direction), head(down, per_direction))
  labels <- labels[order(labels$fdr, -abs(labels$log2_fold_change)), , drop = FALSE]
  head(labels, maximum_labels)
}

mouse_expression_audit <- function(model, feature_ids, de_table, mouse_meta) {
  if (!length(feature_ids)) {
    return(data.frame(
      feature_id = character(),
      gene_symbol = character(),
      sample_id = character(),
      normalized_log2_cpm = numeric(),
      voom_precision_weight = numeric(),
      raw_pseudobulk_count = numeric(),
      initial_ploidy = character(),
      dose_mg = numeric(),
      treatment = character(),
      n_interval_cells = integer(),
      mean_interval_pseudotime = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  feature_ids <- intersect(feature_ids, rownames(model$voom$E))
  samples <- colnames(model$voom$E)
  index <- expand.grid(
    feature_id = feature_ids,
    sample_id = samples,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  row_index <- match(index$feature_id, rownames(model$voom$E))
  column_index <- match(index$sample_id, samples)
  meta_index <- match(index$sample_id, mouse_meta$sample_id)
  symbol_index <- match(index$feature_id, de_table$feature_id)
  data.frame(
    feature_id = index$feature_id,
    gene_symbol = de_table$gene_symbol[symbol_index],
    sample_id = index$sample_id,
    normalized_log2_cpm = model$voom$E[cbind(row_index, column_index)],
    voom_precision_weight = model$voom$weights[cbind(row_index, column_index)],
    raw_pseudobulk_count = model$dge$counts[cbind(row_index, column_index)],
    initial_ploidy = mouse_meta$initial_ploidy[meta_index],
    dose_mg = mouse_meta$dose_mg[meta_index],
    treatment = mouse_meta$treatment[meta_index],
    n_interval_cells = mouse_meta$n_cells[meta_index],
    mean_interval_pseudotime = mouse_meta$mean_pseudotime[meta_index],
    stringsAsFactors = FALSE
  )
}

dose_concordance_table <- function(primary, dose_30, dose_120, feature_ids) {
  columns <- c("feature_id", "gene_symbol", "log2_fold_change", "p_value", "fdr")
  primary_table <- primary[match(feature_ids, primary$feature_id), columns, drop = FALSE]
  dose_30_table <- dose_30[match(feature_ids, dose_30$feature_id), columns, drop = FALSE]
  dose_120_table <- dose_120[match(feature_ids, dose_120$feature_id), columns, drop = FALSE]
  data.frame(
    feature_id = primary_table$feature_id,
    gene_symbol = primary_table$gene_symbol,
    pooled_log2_fold_change = primary_table$log2_fold_change,
    pooled_p_value = primary_table$p_value,
    pooled_fdr = primary_table$fdr,
    dose_30_log2_fold_change = dose_30_table$log2_fold_change,
    dose_30_p_value = dose_30_table$p_value,
    dose_30_fdr = dose_30_table$fdr,
    dose_120_log2_fold_change = dose_120_table$log2_fold_change,
    dose_120_p_value = dose_120_table$p_value,
    dose_120_fdr = dose_120_table$fdr,
    concordant_direction_across_doses = sign(dose_30_table$log2_fold_change) ==
      sign(dose_120_table$log2_fold_change),
    stringsAsFactors = FALSE
  )
}

prepare_pathway_ranking <- function(de_table, support) {
  required <- c(
    "feature_id", "gene_symbol", "contrast_id", "moderated_t", "p_value",
    "fdr"
  )
  missing <- setdiff(required, names(de_table))
  if (length(missing)) {
    stop(
      "Differential-expression table is missing pathway-ranking column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  contrast <- data.frame(
    gene = as.character(de_table$feature_id),
    gene_symbol = as.character(de_table$gene_symbol),
    contrast_id = as.character(de_table$contrast_id),
    t_statistic = as.numeric(de_table$moderated_t),
    p_value = as.numeric(de_table$p_value),
    fdr = as.numeric(de_table$fdr),
    stringsAsFactors = FALSE
  )
  resolution <- support$resolve_gene_symbols(contrast)
  stats <- support$ranked_stats(contrast, resolution)
  if (!length(stats) || is.null(names(stats)) || anyNA(stats) ||
      anyNA(names(stats)) || any(!nzchar(names(stats))) ||
      anyDuplicated(names(stats))) {
    stop("Pathway ranking must contain unique, finite gene symbols", call. = FALSE)
  }
  ranking <- data.frame(
    gene_symbol = names(stats),
    moderated_t = as.numeric(stats),
    rank = seq_along(stats),
    stringsAsFactors = FALSE
  )
  list(stats = stats, ranking = ranking, symbol_resolution = resolution)
}

exploratory_pathway_max_per_direction <- function() 3L

interaction_pathway_max_per_direction <- function() 10L

exploratory_pathway_selection_rule <- function() {
  paste(
    "BH-adjusted P <= 0.05;",
    "up to top three per sign and collection"
  )
}

select_pathway_results <- function(
  gsea,
  support,
  config,
  maximum_per_direction = exploratory_pathway_max_per_direction()
) {
  maximum_per_direction <- suppressWarnings(as.integer(maximum_per_direction))
  if (length(maximum_per_direction) != 1L || is.na(maximum_per_direction) ||
      maximum_per_direction < 1L || maximum_per_direction > 4L) {
    stop("maximum_per_direction must be an integer from one to four", call. = FALSE)
  }
  selection_source <- as.data.frame(gsea, stringsAsFactors = FALSE)
  selection_source$collection_id <- as.character(selection_source$collection)
  selection_source$pathway_id <- as.character(selection_source$pathway)
  selected <- support$figure7_select_generated_pathways(
    selection_source,
    config
  )
  selected <- selected[
    selected$selected_rank_within_direction <= maximum_per_direction,
    ,
    drop = FALSE
  ]
  collections <- as.character(unlist(config$state_pathways$collections))
  selected$collection_display_order <- match(
    as.character(selected$collection),
    collections
  )
  selected$pathway_display_order <- seq_len(nrow(selected))
  selected
}

pathway_collection_summary <- function(gsea, selected, support, config) {
  collections <- as.character(unlist(config$state_pathways$collections))
  threshold <- support$figure7_generated_pathway_fdr_threshold()
  rows <- lapply(collections, function(collection) {
    local <- gsea[gsea$collection == collection, , drop = FALSE]
    local_selected <- selected[selected$collection == collection, , drop = FALSE]
    significant <- is.finite(local$padj) & local$padj <= threshold
    data.frame(
      collection_id = collection,
      collection_label = unique(as.character(local$collection_label))[[1L]],
      pathways_tested = nrow(local),
      fdr_significant = sum(significant),
      fdr_significant_positive_nes = sum(significant & local$NES > 0),
      fdr_significant_negative_nes = sum(significant & local$NES < 0),
      pathways_selected_for_panel = nrow(local_selected),
      selected_positive_nes = sum(local_selected$NES > 0),
      selected_negative_nes = sum(local_selected$NES < 0),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_pathway_enrichment <- function(
  primary,
  support,
  config,
  seed = 1L,
  gene_sets = NULL,
  include_gene_set_membership = TRUE,
  contrast_id = "treated_equal_dose_minus_vehicle"
) {
  expected_msigdbr <- as.character(config$gene_sets$package_version)
  observed_msigdbr <- as.character(utils::packageVersion("msigdbr"))
  if (!identical(observed_msigdbr, expected_msigdbr)) {
    stop(
      "Pathway analysis requires msigdbr ", expected_msigdbr,
      "; observed ", observed_msigdbr,
      call. = FALSE
    )
  }
  collections <- as.character(unlist(config$state_pathways$collections))
  species <- as.character(config$gene_sets$species)
  if (is.null(gene_sets)) {
    gene_sets <- support$fetch_gene_sets(
      paste(collections, collapse = ","),
      species = species
    )
  }
  if (!identical(names(gene_sets), collections)) {
    stop("Reused gene sets do not match the configured collection order", call. = FALSE)
  }
  observed_releases <- unique(vapply(
    gene_sets,
    `[[`,
    character(1L),
    "database_release"
  ))
  expected_release <- as.character(config$gene_sets$database_release)
  if (length(observed_releases) != 1L ||
      !identical(observed_releases, expected_release)) {
    stop(
      "Pathway analysis requires MSigDB ", expected_release,
      "; observed ", paste(observed_releases, collapse = ","),
      call. = FALSE
    )
  }

  ranking <- prepare_pathway_ranking(primary, support)
  gsea <- support$run_all_gsea(
    ranking$stats,
    gene_sets,
    contrast_id,
    15L,
    500L,
    as.integer(config$state_pathways$gsea_nperm_simple),
    as.integer(seed),
    as.integer(config$state_pathways$gsea_nperm_simple_max),
    as.integer(config$state_pathways$gsea_nperm_simple_multiplier)
  )
  if (is.null(gsea) || !nrow(gsea)) {
    stop("No eligible pathways were returned by GSEA", call. = FALSE)
  }
  selected <- select_pathway_results(gsea, support, config)
  leading_edge <- support$leading_edge_table(gsea)
  selected_keys <- paste(selected$collection, selected$pathway, sep = "\r")
  leading_edge <- leading_edge[
    paste(leading_edge$collection, leading_edge$pathway, sep = "\r") %in%
      selected_keys,
    ,
    drop = FALSE
  ]
  membership <- if (isTRUE(include_gene_set_membership)) {
    support$gene_set_membership_table(gene_sets)
  } else {
    NULL
  }
  contract <- data.frame(
    key = c(
      "provider", "msigdbr_package_version", "database_release", "species",
      "collections", "ranking_statistic", "gsea_min_size", "gsea_max_size",
      "gsea_nperm_simple", "gsea_nperm_simple_max",
      "gsea_nperm_simple_multiplier", "pathway_selection_rule",
      "ranked_unique_gene_symbols", "gene_set_membership_sha256"
    ),
    value = c(
      "msigdbr", observed_msigdbr, observed_releases, species,
      paste(collections, collapse = ","), "moderated_t", "15", "500",
      as.character(config$state_pathways$gsea_nperm_simple),
      as.character(config$state_pathways$gsea_nperm_simple_max),
      as.character(config$state_pathways$gsea_nperm_simple_multiplier),
      exploratory_pathway_selection_rule(),
      as.character(length(ranking$stats)), NA_character_
    ),
    stringsAsFactors = FALSE
  )
  list(
    complete = gsea,
    selected = selected,
    leading_edge = leading_edge,
    membership = membership,
    gene_sets = if (isTRUE(include_gene_set_membership)) gene_sets else NULL,
    contract = contract,
    ranking = ranking$ranking,
    symbol_resolution = ranking$symbol_resolution,
    collection_summary = pathway_collection_summary(
      gsea,
      selected,
      support,
      config
    )
  )
}

validate_origin_stratum <- function(mouse_meta, origin) {
  origin <- match.arg(origin, c("2N", "4N"))
  expected_dose_counts <- c("0" = 4L, "30" = 2L, "120" = 2L)
  observed_dose_counts <- table(factor(
    as.character(mouse_meta$dose_mg),
    levels = names(expected_dose_counts)
  ))
  failures <- character()
  if (nrow(mouse_meta) != 8L) {
    failures <- c(failures, "stratum does not contain eight mice")
  }
  observed_origins <- unique(as.character(mouse_meta$initial_ploidy))
  if (anyNA(observed_origins) || !identical(observed_origins, origin)) {
    failures <- c(failures, "stratum contains a different injected origin")
  }
  if (!identical(as.integer(observed_dose_counts), unname(expected_dose_counts))) {
    failures <- c(failures, "dose counts are not 4 vehicle, 2 at 30, and 2 at 120 mg/kg")
  }
  if (any(!is.finite(mouse_meta$n_cells)) || any(mouse_meta$n_cells < 1L)) {
    failures <- c(failures, "one or more stratum mice has no interval cells")
  }
  if (length(failures)) {
    stop(
      origin,
      "-origin interval contract failed: ",
      paste(failures, collapse = "; "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

fit_origin_stratified_analysis <- function(
  pseudobulk,
  origin,
  support,
  config,
  gene_sets,
  cpm_threshold = 1,
  seed = 1L
) {
  origin <- match.arg(origin, c("2N", "4N"))
  mouse_meta <- pseudobulk$metadata[
    pseudobulk$metadata$initial_ploidy == origin,
    ,
    drop = FALSE
  ]
  validate_origin_stratum(mouse_meta, origin)
  counts <- pseudobulk$counts[, mouse_meta$sample_id, drop = FALSE]
  design_bundle <- prepare_interval_design(
    mouse_meta,
    adjust_mean_pseudotime = TRUE,
    adjust_initial_ploidy = FALSE
  )
  contrasts <- treatment_contrasts(design_bundle$design)
  model <- fit_interval_voom(
    counts,
    design_bundle,
    cpm_threshold = cpm_threshold,
    minimum_mice = 2L
  )
  de_tables <- lapply(names(contrasts), function(contrast_id) {
    table <- extract_contrast_table(
      model,
      contrasts[[contrast_id]],
      contrast_id,
      support$clean_gene_symbols
    )
    table$origin_stratum <- origin
    table
  })
  names(de_tables) <- names(contrasts)
  primary <- de_tables$treated_equal_dose_minus_vehicle
  pathways <- run_pathway_enrichment(
    primary,
    support,
    config,
    seed = seed,
    gene_sets = gene_sets,
    include_gene_set_membership = FALSE
  )
  list(
    origin = origin,
    mouse_metadata = mouse_meta,
    design = design_bundle,
    contrasts = contrasts,
    model = model,
    de_tables = de_tables,
    primary = primary,
    pathways = pathways
  )
}

fit_global_response_magnitude_model <- function(
  pseudobulk,
  cpm_threshold = 1,
  minimum_mice = 4L
) {
  design_bundle <- prepare_treatment_origin_interaction_design(
    pseudobulk$metadata,
    adjust_mean_pseudotime = TRUE
  )
  contrasts <- treatment_origin_interaction_contrasts(design_bundle$design)
  model <- fit_interval_voom(
    pseudobulk$counts,
    design_bundle,
    cpm_threshold = cpm_threshold,
    minimum_mice = minimum_mice
  )
  list(
    mouse_metadata = pseudobulk$metadata,
    design = design_bundle,
    contrasts = contrasts,
    model = model
  )
}

fit_treatment_origin_interaction_analysis <- function(
  pseudobulk,
  support,
  config,
  gene_sets,
  cpm_threshold = 1,
  minimum_mice = 4L,
  seed = 1L
) {
  base_model <- fit_global_response_magnitude_model(
    pseudobulk,
    cpm_threshold = cpm_threshold,
    minimum_mice = minimum_mice
  )
  design_bundle <- base_model$design
  contrasts <- base_model$contrasts
  model <- base_model$model
  de_tables <- lapply(names(contrasts), function(contrast_id) {
    is_interaction <- identical(contrast_id, "treatment_by_origin_4N_minus_2N")
    extract_contrast_table(
      model,
      contrasts[[contrast_id]],
      contrast_id,
      support$clean_gene_symbols,
      positive_direction = if (is_interaction) {
        "treatment_effect_more_positive_in_4N"
      } else {
        "higher_in_treated"
      },
      negative_direction = if (is_interaction) {
        "treatment_effect_more_positive_in_2N"
      } else {
        "lower_in_treated"
      }
    )
  })
  names(de_tables) <- names(contrasts)
  primary_id <- "treatment_by_origin_4N_minus_2N"
  primary <- de_tables[[primary_id]]
  pathway_contrast_ids <- c(
    "treated_equal_dose_minus_vehicle_2N",
    "treated_equal_dose_minus_vehicle_4N",
    primary_id
  )
  pathway_results <- lapply(pathway_contrast_ids, function(contrast_id) {
    run_pathway_enrichment(
      de_tables[[contrast_id]],
      support,
      config,
      seed = seed,
      gene_sets = gene_sets,
      include_gene_set_membership = FALSE,
      contrast_id = contrast_id
    )
  })
  names(pathway_results) <- pathway_contrast_ids
  list(
    contrast_id = primary_id,
    mouse_metadata = pseudobulk$metadata,
    design = design_bundle,
    contrasts = contrasts,
    model = model,
    de_tables = de_tables,
    primary = primary,
    pathways = pathway_results[[primary_id]],
    origin_effect_pathways = list(
      `2N` = pathway_results[["treated_equal_dose_minus_vehicle_2N"]],
      `4N` = pathway_results[["treated_equal_dose_minus_vehicle_4N"]]
    )
  )
}

global_response_contrast_ids <- function() {
  c(
    `2N` = "treated_equal_dose_minus_vehicle_2N",
    `4N` = "treated_equal_dose_minus_vehicle_4N"
  )
}

validate_magnitude_bootstrap_settings <- function(replicates, seed) {
  replicates <- suppressWarnings(as.integer(replicates))
  seed <- suppressWarnings(as.integer(seed))
  if (length(replicates) != 1L || is.na(replicates) || replicates < 199L) {
    stop("magnitude_bootstrap_replicates must be an integer of at least 199", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed)) {
    stop("magnitude_bootstrap_seed must be a finite integer", call. = FALSE)
  }
  list(replicates = replicates, seed = seed)
}

prepare_global_response_magnitude_components <- function(interaction_analysis) {
  contrast_ids <- global_response_contrast_ids()
  missing_contrasts <- setdiff(unname(contrast_ids), names(interaction_analysis$contrasts))
  if (length(missing_contrasts)) {
    stop(
      "Global response comparison is missing contrast(s): ",
      paste(missing_contrasts, collapse = ", "),
      call. = FALSE
    )
  }

  model <- interaction_analysis$model
  expression <- as.matrix(model$voom$E)
  design <- interaction_analysis$design$design
  coefficients <- model$fit$coefficients
  if (!identical(colnames(expression), rownames(design))) {
    stop("Voom expression columns and interaction-design rows are not aligned", call. = FALSE)
  }
  if (!identical(rownames(expression), rownames(coefficients))) {
    stop("Voom expression and fitted coefficients use different gene orders", call. = FALSE)
  }
  if (!identical(colnames(coefficients), colnames(design))) {
    stop("Fitted coefficients and interaction-design columns are not aligned", call. = FALSE)
  }

  voom_weights <- model$voom$weights
  if (is.null(voom_weights)) {
    voom_weights <- matrix(1, nrow(expression), ncol(expression))
  } else {
    voom_weights <- as.matrix(voom_weights)
    if (nrow(voom_weights) == 1L && nrow(expression) > 1L) {
      voom_weights <- voom_weights[rep(1L, nrow(expression)), , drop = FALSE]
    }
    if (ncol(voom_weights) == 1L && ncol(expression) > 1L) {
      voom_weights <- voom_weights[, rep(1L, ncol(expression)), drop = FALSE]
    }
  }
  if (!identical(dim(voom_weights), dim(expression)) ||
      any(!is.finite(voom_weights)) || any(voom_weights <= 0)) {
    stop("Global response comparison requires positive gene-by-mouse voom weights", call. = FALSE)
  }

  contrast_matrix <- do.call(cbind, lapply(unname(contrast_ids), function(contrast_id) {
    contrast <- interaction_analysis$contrasts[[contrast_id]]
    if (!identical(names(contrast), colnames(design))) {
      stop("Global response contrast coefficients are not aligned", call. = FALSE)
    }
    contrast
  }))
  colnames(contrast_matrix) <- names(contrast_ids)
  contrast_fit <- limma::contrasts.fit(model$fit, contrasts = contrast_matrix)
  contrast_fit <- limma::eBayes(contrast_fit, robust = TRUE)
  effect <- contrast_fit$coefficients[, names(contrast_ids), drop = FALSE]
  moderated_variance <- contrast_fit$stdev.unscaled[, names(contrast_ids), drop = FALSE]^2 *
    contrast_fit$s2.post
  if (any(!is.finite(effect)) || any(!is.finite(moderated_variance)) ||
      any(moderated_variance < 0)) {
    stop("Global response effects or their sampling variances are not finite", call. = FALSE)
  }

  fitted <- coefficients %*% t(design)
  residual <- expression - fitted
  n_genes <- nrow(expression)
  n_mice <- ncol(expression)
  influence <- lapply(names(contrast_ids), function(origin) {
    matrix(NA_real_, nrow = n_genes, ncol = n_mice)
  })
  names(influence) <- names(contrast_ids)
  leverage <- matrix(NA_real_, nrow = n_genes, ncol = n_mice)

  for (gene_index in seq_len(n_genes)) {
    gene_weights <- voom_weights[gene_index, ]
    weighted_design <- design * gene_weights
    information <- crossprod(design, weighted_design)
    information_inverse <- tryCatch(
      chol2inv(chol(information)),
      error = function(error) solve(information)
    )
    design_information_inverse <- design %*% information_inverse
    leverage[gene_index, ] <- gene_weights * rowSums(
      design_information_inverse * design
    )
    for (origin in names(contrast_ids)) {
      influence[[origin]][gene_index, ] <- gene_weights * as.vector(
        design_information_inverse %*% contrast_matrix[, origin]
      )
    }
  }
  if (any(!is.finite(leverage)) || any(leverage < -1e-8) || any(leverage >= 1)) {
    stop("Invalid weighted leverage in the global response comparison", call. = FALSE)
  }

  reconstructed_effect <- vapply(names(contrast_ids), function(origin) {
    rowSums(influence[[origin]] * expression)
  }, numeric(n_genes))
  colnames(reconstructed_effect) <- names(contrast_ids)
  if (!isTRUE(all.equal(
    unname(reconstructed_effect),
    unname(effect),
    tolerance = 1e-7
  ))) {
    stop("Weighted contrast reconstruction did not reproduce the fitted effects", call. = FALSE)
  }

  hc2_residual <- residual / sqrt(pmax(1 - leverage, 1e-8))
  perturbation <- lapply(names(contrast_ids), function(origin) {
    influence[[origin]] * hc2_residual
  })
  names(perturbation) <- names(contrast_ids)
  wild_variance <- vapply(perturbation, function(x) {
    rowSums(x^2)
  }, numeric(n_genes))
  colnames(wild_variance) <- names(contrast_ids)

  origin_summary <- do.call(rbind, lapply(names(contrast_ids), function(origin) {
    naive_energy <- mean(effect[, origin]^2)
    mean_sampling_variance <- mean(moderated_variance[, origin])
    corrected_energy <- naive_energy - mean_sampling_variance
    data.frame(
      origin = origin,
      retained_common_genes = n_genes,
      naive_mean_squared_log2_fold_change = naive_energy,
      mean_moderated_sampling_variance = mean_sampling_variance,
      noise_corrected_mean_squared_log2_fold_change = corrected_energy,
      naive_rms_log2_fold_change = sqrt(naive_energy),
      noise_corrected_rms_log2_fold_change = if (corrected_energy >= 0) {
        sqrt(corrected_energy)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(origin_summary) <- NULL

  effect_cosine <- sum(effect[, "2N"] * effect[, "4N"]) /
    sqrt(sum(effect[, "2N"]^2) * sum(effect[, "4N"]^2))
  list(
    effect = effect,
    moderated_variance = moderated_variance,
    perturbation = perturbation,
    wild_variance = wild_variance,
    origin_summary = origin_summary,
    effect_cosine = effect_cosine,
    n_genes = n_genes,
    n_mice = n_mice
  )
}

bootstrap_global_response_magnitude <- function(
  components,
  replicates = 4999L,
  seed = 20260818L,
  chunk_size = 200L
) {
  settings <- validate_magnitude_bootstrap_settings(replicates, seed)
  replicates <- settings$replicates
  seed <- settings$seed
  chunk_size <- suppressWarnings(as.integer(chunk_size))
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be a positive integer", call. = FALSE)
  }

  effect_2n <- components$effect[, "2N"]
  effect_4n <- components$effect[, "4N"]
  observed_energy <- setNames(
    components$origin_summary$noise_corrected_mean_squared_log2_fold_change,
    components$origin_summary$origin
  )
  observed_difference <- observed_energy[["2N"]] - observed_energy[["4N"]]
  bootstrap_parameter <- mean(effect_2n^2 - effect_4n^2)
  wild_bias_difference <- mean(
    components$wild_variance[, "2N"] - components$wild_variance[, "4N"]
  )

  set.seed(seed)
  bootstrap_estimate <- numeric(replicates)
  starts <- seq.int(1L, replicates, by = chunk_size)
  for (start in starts) {
    stop_index <- min(replicates, start + chunk_size - 1L)
    draw_count <- stop_index - start + 1L
    multipliers <- matrix(
      sample(c(-1, 1), components$n_mice * draw_count, replace = TRUE),
      nrow = components$n_mice,
      ncol = draw_count
    )
    perturbation_2n <- components$perturbation[["2N"]] %*% multipliers
    perturbation_4n <- components$perturbation[["4N"]] %*% multipliers
    bootstrap_estimate[start:stop_index] <- colMeans(
      (effect_2n + perturbation_2n)^2 -
        (effect_4n + perturbation_4n)^2
    ) - wild_bias_difference
  }

  bootstrap_error <- bootstrap_estimate - bootstrap_parameter
  error_quantiles <- stats::quantile(
    bootstrap_error,
    probs = c(0.025, 0.975),
    names = FALSE,
    type = 8
  )
  confidence_interval <- c(
    observed_difference - error_quantiles[[2L]],
    observed_difference - error_quantiles[[1L]]
  )
  upper_tail <- (sum(bootstrap_error >= observed_difference) + 1) /
    (replicates + 1)
  lower_tail <- (sum(bootstrap_error <= observed_difference) + 1) /
    (replicates + 1)

  comparison <- data.frame(
    estimand = paste(
      "noise-corrected transcriptome-wide mean squared log2 fold change:",
      "2N minus 4N"
    ),
    estimate_2N_minus_4N = observed_difference,
    confidence_level = 0.95,
    confidence_interval_lower = confidence_interval[[1L]],
    confidence_interval_upper = confidence_interval[[2L]],
    bootstrap_standard_error = stats::sd(bootstrap_error),
    one_sided_p_2N_greater = upper_tail,
    two_sided_p = min(1, 2 * min(upper_tail, lower_tail)),
    bootstrap_replicates = replicates,
    bootstrap_seed = seed,
    response_vector_cosine_2N_vs_4N = components$effect_cosine,
    common_gene_universe = components$n_genes,
    inference_method = paste(
      "HC2 whole-mouse Rademacher wild bootstrap on fixed-voom residuals;",
      "the same multiplier is applied to every gene within a mouse"
    ),
    interpretation_if_positive = paste(
      "larger total transcriptional displacement in 2N; this does not imply",
      "that the same pathways change in the same direction"
    ),
    stringsAsFactors = FALSE
  )
  bootstrap <- data.frame(
    bootstrap_replicate = seq_len(replicates),
    corrected_energy_difference_2N_minus_4N = bootstrap_estimate,
    centered_bootstrap_error = bootstrap_error,
    stringsAsFactors = FALSE
  )
  list(
    origin_summary = components$origin_summary,
    comparison = comparison,
    bootstrap = bootstrap
  )
}

analyze_global_response_magnitude <- function(
  interaction_analysis,
  replicates = 4999L,
  seed = 20260818L
) {
  components <- prepare_global_response_magnitude_components(interaction_analysis)
  bootstrap_global_response_magnitude(
    components,
    replicates = replicates,
    seed = seed
  )
}

write_global_response_magnitude_outputs <- function(result, table_dir) {
  write_tsv(
    result$origin_summary,
    file.path(table_dir, "global_treatment_response_magnitude_by_origin.tsv")
  )
  write_tsv(
    result$comparison,
    file.path(table_dir, "global_treatment_response_magnitude_comparison.tsv")
  )
  write_tsv(
    result$bootstrap,
    file.path(table_dir, "global_treatment_response_magnitude_wild_bootstrap.tsv")
  )
  invisible(result)
}

run_global_magnitude_only <- function(
  args,
  repo_root,
  support,
  config,
  interval,
  pseudobulk,
  species_audit,
  match_audit,
  input_hashes
) {
  message("Fitting the joint treatment-by-origin model for the magnitude audit")
  interaction_model <- fit_global_response_magnitude_model(
    pseudobulk,
    cpm_threshold = args$cpm_threshold,
    minimum_mice = args$minimum_mice
  )
  magnitude <- analyze_global_response_magnitude(
    interaction_model,
    replicates = args$magnitude_bootstrap_replicates,
    seed = args$magnitude_bootstrap_seed
  )

  output_root <- prepare_output_root(args$output_root, overwrite = args$overwrite)
  table_dir <- file.path(output_root, "tables")
  metadata_dir <- file.path(output_root, "metadata")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  write_global_response_magnitude_outputs(magnitude, table_dir)
  write_tsv(
    pseudobulk$metadata,
    file.path(table_dir, "interval_mouse_coverage.tsv")
  )
  design_output <- cbind(
    sample_id = rownames(interaction_model$design$design),
    as.data.frame(interaction_model$design$design, check.names = FALSE)
  )
  write_tsv(
    design_output,
    file.path(table_dir, "design_matrix_treatment_by_origin_interaction.tsv")
  )
  write_tsv(
    model_contrast_table(interaction_model$contrasts),
    file.path(table_dir, "contrast_definitions_treatment_by_origin_interaction.tsv")
  )
  write_tsv(species_audit, file.path(metadata_dir, "feature_species_audit.tsv"))
  write_tsv(match_audit, file.path(metadata_dir, "cell_expression_match_audit.tsv"))
  write_tsv(data.frame(
    interval_id = interval$name,
    interval_start = interval$start,
    interval_end = interval$end,
    include_start = interval$include_start,
    include_end = interval$include_end,
    interval_source = interval$source,
    interval_selection_caveat = paste(
      "The interval was localized using treated-versus-vehicle density differences;",
      "the global magnitude comparison is exploratory and conditional on that selection."
    ),
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "interval_definition.tsv"))
  write_tsv(data.frame(
    cpm_threshold = interaction_model$model$cpm_threshold,
    minimum_mice_above_threshold = interaction_model$model$minimum_mice,
    rule = sprintf(
      "CPM > %s in at least %d independent mouse pseudobulks",
      interaction_model$model$cpm_threshold,
      interaction_model$model$minimum_mice
    ),
    n_input_features = interaction_model$model$n_input_genes,
    n_retained_features = interaction_model$model$n_retained_genes,
    shared_between_origins = TRUE,
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "expression_filter_audit.tsv"))

  script <- exploratory_script_path()
  manifest_paths <- c(
    cell_metadata = args$cell_metadata,
    seurat_rds = args$seurat_rds,
    config = args$config,
    exploratory_script = script,
    feature_species_policy = file.path(
      repo_root,
      "Code/in-vivo/figure7/src/feature_species_policy.R"
    ),
    shared_figure7_contract = file.path(
      repo_root,
      "Code/in-vivo/figure7/src/common_io.R"
    )
  )
  manifest_hashes <- setNames(rep(NA_character_, length(manifest_paths)), names(manifest_paths))
  manifest_hashes[names(input_hashes)] <- input_hashes
  hashes_to_compute <- names(manifest_hashes)[is.na(manifest_hashes)]
  manifest_hashes[hashes_to_compute] <- vapply(
    manifest_paths[hashes_to_compute],
    sha256_file,
    character(1L)
  )
  manifest_locators <- vapply(
    manifest_paths,
    portable_locator,
    character(1L),
    repo_root = repo_root
  )
  manifest_locators[["seurat_rds"]] <- file.path(
    as.character(config$raw_data$default_root),
    as.character(config$raw_data$seurat_filename)
  )
  write_tsv(data.frame(
    input_id = names(manifest_paths),
    locator = unname(manifest_locators),
    sha256 = unname(manifest_hashes),
    byte_size = as.numeric(file.info(manifest_paths)$size),
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "input_manifest.tsv"))

  packages <- global_magnitude_required_packages()
  write_tsv(data.frame(
    package = packages,
    version = vapply(packages, function(package) {
      as.character(utils::packageVersion(package))
    }, character(1L)),
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "package_versions.tsv"))
  write_tsv(data.frame(
    key = c(
      "git_revision", "git_branch", "git_worktree_dirty",
      "analysis_scope", "feature_species_policy", "normalization",
      "observation_model", "global_response_estimand", "inference"
    ),
    value = c(
      git_value(repo_root, c("rev-parse", "HEAD")),
      git_value(repo_root, c("branch", "--show-current")),
      as.character(nzchar(git_value(
        repo_root,
        c("status", "--porcelain", "--untracked-files=no")
      ) %||% "")),
      "global magnitude audit only; no pathway database or GSEA was loaded",
      support$figure7_human_feature_policy_description(),
      "edgeR TMM and limma voom",
      "joint 16-mouse origin-by-dose model adjusted for mean interval pseudotime",
      "noise-corrected mean squared log2 fold change; 2N minus 4N",
      paste0(
        "HC2 whole-mouse Rademacher wild bootstrap; ",
        args$magnitude_bootstrap_replicates,
        " replicates; seed ",
        args$magnitude_bootstrap_seed
      )
    ),
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "provenance.tsv"))
  session_info <- sub(
    "[[:space:]]+$",
    "",
    capture.output(utils::sessionInfo())
  )
  writeLines(session_info, file.path(metadata_dir, "sessionInfo.txt"))

  message("Completed global response-magnitude audit: ", output_root)
  invisible(list(
    interaction_model = interaction_model,
    global_response_magnitude = magnitude,
    mouse_metadata = pseudobulk$metadata,
    output_root = output_root
  ))
}

prepare_origin_interaction_pathways <- function(
  gsea_2n,
  gsea_4n,
  interaction_gsea,
  origin_effect_gsea_2n = gsea_2n,
  origin_effect_gsea_4n = gsea_4n,
  fdr_threshold = 0.05,
  maximum_per_direction = interaction_pathway_max_per_direction()
) {
  maximum_per_direction <- suppressWarnings(as.integer(maximum_per_direction))
  if (length(maximum_per_direction) != 1L || is.na(maximum_per_direction) ||
      maximum_per_direction < 1L) {
    stop("maximum_per_direction must be a positive integer", call. = FALSE)
  }
  if (length(fdr_threshold) != 1L || !is.finite(fdr_threshold) ||
      fdr_threshold <= 0 || fdr_threshold >= 1) {
    stop("fdr_threshold must be strictly between zero and one", call. = FALSE)
  }
  required <- c(
    "collection", "collection_label", "pathway", "pathway_label", "NES", "padj"
  )
  named_tables <- list(
    stratified_2N = gsea_2n,
    stratified_4N = gsea_4n,
    interaction = interaction_gsea,
    joint_model_2N_effect = origin_effect_gsea_2n,
    joint_model_4N_effect = origin_effect_gsea_4n
  )
  for (table_id in names(named_tables)) {
    entry <- named_tables[[table_id]]
    missing <- setdiff(required, names(entry))
    if (length(missing)) {
      stop(
        table_id,
        " GSEA is missing column(s): ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    keys <- paste(entry$collection, entry$pathway, sep = "\r")
    if (anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
      stop(table_id, " GSEA requires unique pathway keys", call. = FALSE)
    }
  }

  comparison_columns <- c(
    "collection", "pathway", "pathway_label", "NES", "padj"
  )
  paired <- merge(
    gsea_2n[, comparison_columns, drop = FALSE],
    gsea_4n[, comparison_columns, drop = FALSE],
    by = c("collection", "pathway"),
    suffixes = c("_2N", "_4N")
  )
  paired$pathway_label <- ifelse(
    nzchar(as.character(paired$pathway_label_2N)),
    as.character(paired$pathway_label_2N),
    as.character(paired$pathway_label_4N)
  )
  shared <- paired[
    is.finite(paired$padj_2N) & paired$padj_2N <= fdr_threshold &
      is.finite(paired$padj_4N) & paired$padj_4N <= fdr_threshold,
    c(
      "collection", "pathway", "pathway_label",
      "NES_2N", "padj_2N", "NES_4N", "padj_4N"
    ),
    drop = FALSE
  ]
  if (!nrow(shared)) {
    stop("No FDR-significant pathway is shared between the origin models", call. = FALSE)
  }
  shared <- shared[
    order((shared$NES_2N + shared$NES_4N) / 2, shared$collection, shared$pathway),
    ,
    drop = FALSE
  ]

  effect_paired <- merge(
    origin_effect_gsea_2n[, comparison_columns, drop = FALSE],
    origin_effect_gsea_4n[, comparison_columns, drop = FALSE],
    by = c("collection", "pathway"),
    suffixes = c("_2N", "_4N")
  )

  interaction <- as.data.frame(interaction_gsea, stringsAsFactors = FALSE)
  interaction <- interaction[
    is.finite(interaction$padj) & is.finite(interaction$NES) & interaction$NES != 0,
    ,
    drop = FALSE
  ]
  interaction$more_positive_origin <- ifelse(interaction$NES < 0, "2N", "4N")
  interaction$interaction_fdr_significant <- interaction$padj <= fdr_threshold
  selected_rows <- lapply(c("2N", "4N"), function(origin) {
    local <- interaction[interaction$more_positive_origin == origin, , drop = FALSE]
    local <- if (identical(origin, "2N")) {
      local[order(local$NES, local$padj, local$collection, local$pathway), , drop = FALSE]
    } else {
      local[order(-local$NES, local$padj, local$collection, local$pathway), , drop = FALSE]
    }
    local <- head(local, maximum_per_direction)
    local$selected_rank_within_direction <- seq_len(nrow(local))
    local$panel_id <- if (identical(origin, "2N")) {
      "negative_interaction"
    } else {
      "positive_interaction"
    }
    local$selection_scope <- paste0(
      "Top ",
      if (identical(origin, "2N")) "negative" else "positive",
      " formal treatment-by-origin interaction NES across collections"
    )
    paired_index <- match(
      paste(local$collection, local$pathway, sep = "\r"),
      paste(effect_paired$collection, effect_paired$pathway, sep = "\r")
    )
    if (anyNA(paired_index)) {
      stop("A selected interaction pathway is absent from a joint-model origin effect", call. = FALSE)
    }
    local$NES_2N <- effect_paired$NES_2N[paired_index]
    local$padj_2N <- effect_paired$padj_2N[paired_index]
    local$NES_4N <- effect_paired$NES_4N[paired_index]
    local$padj_4N <- effect_paired$padj_4N[paired_index]
    local$origin_effect_nes_order_concordant <- if (identical(origin, "2N")) {
      local$NES_2N > local$NES_4N
    } else {
      local$NES_4N > local$NES_2N
    }
    local
  })
  names(selected_rows) <- c("2N", "4N")
  if (any(vapply(selected_rows, nrow, integer(1L)) == 0L)) {
    stop("Formal interaction GSEA requires at least one pathway in each NES direction", call. = FALSE)
  }
  selected <- do.call(rbind, selected_rows)
  rownames(selected) <- NULL
  list(
    shared = shared,
    selected_interactions = selected,
    selected_negative = selected_rows[["2N"]],
    selected_positive = selected_rows[["4N"]],
    selected_2N = selected_rows[["2N"]],
    selected_4N = selected_rows[["4N"]],
    interaction_complete = interaction_gsea,
    fdr_threshold = fdr_threshold,
    maximum_per_direction = as.integer(maximum_per_direction),
    n_interaction_fdr_significant = sum(
      is.finite(interaction_gsea$padj) & interaction_gsea$padj <= fdr_threshold
    )
  )
}

plot_candidate_pathways <- function(
  selected,
  path_pdf,
  path_png,
  mouse_meta,
  interval,
  selection_rule,
  title = "Pathway enrichment in the selected pseudotime interval",
  adjustment_caption =
    "Adjusted for injected origin and mean within-interval pseudotime."
) {
  plot_data <- selected
  plot_data <- plot_data[
    order(plot_data$collection_display_order, plot_data$NES),
    ,
    drop = FALSE
  ]
  pathway_key <- paste(plot_data$collection, plot_data$pathway, sep = "\r")
  plot_data$pathway_plot_key <- factor(pathway_key, levels = pathway_key)
  pathway_labels <- setNames(
    vapply(
      as.character(plot_data$pathway_label),
      function(label) paste(strwrap(label, width = 44L), collapse = "\n"),
      character(1L)
    ),
    pathway_key
  )
  collection_labels <- c(
    "H" = "Hallmark",
    "C2:CP:REACTOME" = "Reactome",
    "C5:GO:BP" = "GO BP"
  )
  plot_data$collection_label <- factor(
    as.character(plot_data$collection),
    levels = names(collection_labels),
    labels = unname(collection_labels)
  )
  plot_data$enrichment_direction <- factor(
    ifelse(plot_data$NES > 0, "Enriched in treated", "Enriched in vehicle"),
    levels = c("Enriched in vehicle", "Enriched in treated")
  )
  plot_data$minus_log10_fdr <- -log10(
    pmax(as.numeric(plot_data$padj), .Machine$double.xmin)
  )
  n_vehicle <- sum(mouse_meta$treatment == "vehicle")
  n_treated <- sum(mouse_meta$treatment == "treated")
  n_cells <- sum(mouse_meta$n_cells)
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(y = pathway_plot_key)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "#6B7280",
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = NES, yend = pathway_plot_key),
      color = "#9CA3AF",
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = NES,
        fill = enrichment_direction,
        size = minus_log10_fdr
      ),
      shape = 21,
      color = "white",
      stroke = 0.35
    ) +
    ggplot2::facet_grid(
      collection_label ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_y_discrete(labels = pathway_labels) +
    ggplot2::scale_fill_manual(values = c(
      "Enriched in vehicle" = "#2B6CB0",
      "Enriched in treated" = "#C2413B"
    )) +
    ggplot2::scale_size_continuous(range = c(2.4, 5.5)) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        paste0(
          "Equal-dose treated - vehicle contrast; %.3f-%.3f\n",
          "%d cells across %d mice (%d treated, %d vehicle)"
        ),
        interval$start,
        interval$end,
        n_cells,
        nrow(mouse_meta),
        n_treated,
        n_vehicle
      ),
      x = "Normalized enrichment score (treated - vehicle)",
      y = NULL,
      fill = NULL,
      size = expression(-log[10]~"FDR"),
      caption = paste0(
        "Genes ranked by moderated t statistic; ", selection_rule, ".\n",
        adjustment_caption, "\n",
        "Exploratory: the interval was localized using treated-versus-vehicle ",
        "density differences."
      )
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "#374151"),
      plot.caption = ggplot2::element_text(
        size = 7.5,
        hjust = 0,
        color = "#4B5563"
      ),
      legend.position = "top",
      legend.justification = "left",
      strip.text.y = ggplot2::element_text(angle = 270),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  dir.create(dirname(path_pdf), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path_pdf, plot, width = 9, height = 8, units = "in")
  ggplot2::ggsave(
    path_png,
    plot,
    width = 9,
    height = 8,
    units = "in",
    dpi = 320
  )
  invisible(plot)
}

plot_origin_interaction_pathways <- function(
  comparison,
  path_pdf,
  path_png,
  mouse_metadata,
  interval
) {
  collection_labels <- c(
    "H" = "Hallmark",
    "C2:CP:REACTOME" = "Reactome",
    "C5:GO:BP" = "GO BP"
  )
  plot_data <- comparison$selected_interactions
  if (any(!plot_data$collection %in% names(collection_labels))) {
    stop("Origin-interaction plot contains an unknown collection", call. = FALSE)
  }
  panel_labels <- c(
    negative_interaction = "Negative interaction NES\n(more positive in 2N)",
    positive_interaction = "Positive interaction NES\n(more positive in 4N)"
  )
  plot_data$panel_label <- factor(
    unname(panel_labels[plot_data$panel_id]),
    levels = unname(panel_labels)
  )
  plot_data$collection_label <- factor(
    as.character(plot_data$collection),
    levels = names(collection_labels),
    labels = unname(collection_labels)
  )
  plot_data$minus_log10_fdr <- -log10(
    pmax(as.numeric(plot_data$padj), .Machine$double.xmin)
  )
  plot_data$pathway_plot_key <- paste(
    plot_data$panel_id,
    plot_data$collection,
    plot_data$pathway,
    sep = "\r"
  )

  ordered_pathway_keys <- unlist(lapply(names(panel_labels), function(panel_id) {
    local <- plot_data[plot_data$panel_id == panel_id, , drop = FALSE]
    local <- local[
      order(-local$selected_rank_within_direction),
      ,
      drop = FALSE
    ]
    local$pathway_plot_key
  }), use.names = FALSE)
  plot_data$pathway_plot_key <- factor(
    plot_data$pathway_plot_key,
    levels = ordered_pathway_keys
  )
  pathway_label_lookup <- setNames(
    as.character(plot_data$pathway_label),
    as.character(plot_data$pathway_plot_key)
  )
  pathway_labels <- setNames(
    vapply(
      pathway_label_lookup[ordered_pathway_keys],
      function(label) paste(strwrap(label, width = 43L), collapse = "\n"),
      character(1L)
    ),
    ordered_pathway_keys
  )

  cells_2n <- sum(mouse_metadata[["2N"]]$n_cells)
  cells_4n <- sum(mouse_metadata[["4N"]]$n_cells)
  symmetric_limit <- max(abs(plot_data$NES), na.rm = TRUE) * 1.08
  n_significant <- sum(plot_data$interaction_fdr_significant)

  plot <- ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "#6B7280",
      linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(
        x = 0,
        xend = NES,
        y = pathway_plot_key,
        yend = pathway_plot_key
      ),
      color = "#A7AFBC",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(
        x = NES,
        y = pathway_plot_key,
        color = collection_label,
        size = minus_log10_fdr
      )
    ) +
    ggplot2::facet_grid(
      panel_label ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-symmetric_limit, symmetric_limit),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_discrete(labels = pathway_labels) +
    ggplot2::scale_color_manual(values = c(
      "Hallmark" = "#0072B2",
      "Reactome" = "#D55E00",
      "GO BP" = "#009E73"
    )) +
    ggplot2::scale_size_continuous(range = c(2.5, 5.6)) +
    ggplot2::labs(
      title = "Origin-dependent treatment responses",
      subtitle = sprintf(
        paste0(
          "Interaction: [treated - vehicle]4N - [treated - vehicle]2N; %.3f-%.3f pseudotime\n",
          "2N: %d cells across 8 mice; 4N: %d cells across 8 mice"
        ),
        interval$start,
        interval$end,
        cells_2n,
        cells_4n
      ),
      x = "Treatment-by-origin interaction NES",
      y = NULL,
      color = "Collection",
      size = expression(-log[10]~"FDR"),
      caption = paste0(
        "Top 10 negative and top 10 positive formal interaction NES across Hallmark, Reactome, and GO BP; ",
        "ranked by NES within sign.\n",
        "Negative NES denotes a more positive treatment response in 2N; positive NES denotes a more positive response in 4N ",
        "(more positive can also mean less depleted).\n",
        "Point size is based on the interaction FDR; ",
        sprintf(
          "%d/%d displayed pathways pass FDR <= 0.05 (%d significant interaction pathways total).",
          n_significant,
          nrow(plot_data),
          comparison$n_interaction_fdr_significant
        )
      )
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      size = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "#374151"),
      plot.caption = ggplot2::element_text(
        size = 6.9,
        hjust = 0,
        color = "#4B5563"
      ),
      legend.position = "top",
      legend.justification = "left",
      strip.text.y = ggplot2::element_text(angle = 270, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  dir.create(dirname(path_pdf), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path_pdf, plot, width = 10, height = 10.5, units = "in")
  ggplot2::ggsave(
    path_png,
    plot,
    width = 10,
    height = 10.5,
    units = "in",
    dpi = 320
  )
  invisible(plot)
}

plot_candidate_volcano <- function(
  table,
  path_pdf,
  path_png,
  mouse_meta,
  interval,
  fdr_threshold,
  lfc_threshold,
  maximum_labels
) {
  plot_data <- table
  plot_data$minus_log10_fdr <- -log10(pmax(plot_data$fdr, .Machine$double.xmin))
  plot_data$significance <- "Not significant"
  plot_data$significance[
    is.finite(plot_data$fdr) & plot_data$fdr <= fdr_threshold &
      is.finite(plot_data$log2_fold_change) & plot_data$log2_fold_change >= lfc_threshold
  ] <- "Higher in treated"
  plot_data$significance[
    is.finite(plot_data$fdr) & plot_data$fdr <= fdr_threshold &
      is.finite(plot_data$log2_fold_change) & plot_data$log2_fold_change <= -lfc_threshold
  ] <- "Lower in treated"
  plot_data$significance <- factor(
    plot_data$significance,
    levels = c("Not significant", "Lower in treated", "Higher in treated")
  )
  labels <- select_volcano_labels(
    plot_data,
    maximum_labels,
    fdr_threshold,
    lfc_threshold
  )
  n_vehicle <- sum(mouse_meta$treatment == "vehicle")
  n_treated <- sum(mouse_meta$treatment == "treated")
  n_cells <- sum(mouse_meta$n_cells)
  n_up <- sum(plot_data$significance == "Higher in treated", na.rm = TRUE)
  n_down <- sum(plot_data$significance == "Lower in treated", na.rm = TRUE)
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(log2_fold_change, minus_log10_fdr, color = significance)
  ) +
    ggplot2::geom_point(size = 1.25, alpha = 0.72) +
    ggplot2::geom_vline(
      xintercept = c(-lfc_threshold, lfc_threshold),
      color = "#6B7280",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_threshold),
      color = "#6B7280",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(label = gene_symbol),
      size = 3.0,
      min.segment.length = 0,
      box.padding = 0.35,
      point.padding = 0.18,
      max.overlaps = Inf,
      seed = 173L,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = c(
      "Not significant" = "#B8BDC5",
      "Lower in treated" = "#2B6CB0",
      "Higher in treated" = "#C2413B"
    )) +
    ggplot2::labs(
      title = "Differential expression within the pseudotime interval of interest",
      subtitle = sprintf(
        "30 and 120 mg/kg equally weighted minus vehicle; %.3f-%.3f; %d cells across %d mice (%d treated, %d vehicle)",
        interval$start,
        interval$end,
        n_cells,
        nrow(mouse_meta),
        n_treated,
        n_vehicle
      ),
      x = expression(log[2]~fold~change~("gemcitabine" - "vehicle")),
      y = expression(-log[10]~"FDR"),
      color = NULL,
      caption = sprintf(
        paste0(
          "Adjusted for injected origin and mean within-interval pseudotime; ",
          "FDR <= %.2g and |log2 FC| >= %.2g: %d higher, %d lower\n",
          "Exploratory: the interval was localized using treated-versus-vehicle density differences"
        ),
        fdr_threshold,
        lfc_threshold,
        n_up,
        n_down
      )
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "#374151"),
      plot.caption = ggplot2::element_text(size = 7.5, hjust = 0, color = "#4B5563"),
      legend.position = "top",
      legend.justification = "left",
      legend.key.width = grid::unit(0.8, "lines")
    )
  dir.create(dirname(path_pdf), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path_pdf, plot, width = 7.2, height = 5.3, units = "in")
  ggplot2::ggsave(path_png, plot, width = 7.2, height = 5.3, units = "in", dpi = 320)
  invisible(plot)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(as.data.frame(x, stringsAsFactors = FALSE), path, na = "")
  invisible(path)
}

portable_locator <- function(path, repo_root) {
  normalized <- normalizePath(path, mustWork = TRUE)
  root <- normalizePath(repo_root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (startsWith(normalized, prefix)) substring(normalized, nchar(prefix) + 1L) else paste0("external:", basename(normalized))
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

verify_input_hash <- function(path, expected, label) {
  observed <- sha256_file(path)
  if (!identical(tolower(observed), tolower(as.character(expected)))) {
    stop(
      label,
      " SHA-256 differs from the reviewed config: expected ",
      expected,
      "; observed ",
      observed,
      call. = FALSE
    )
  }
  observed
}

git_value <- function(repo_root, arguments) {
  value <- tryCatch(
    system2("git", c("-C", repo_root, arguments), stdout = TRUE, stderr = FALSE),
    error = function(error) character()
  )
  if (!length(value)) NA_character_ else paste(value, collapse = "\n")
}

prepare_output_root <- function(path, overwrite = FALSE) {
  if (dir.exists(path)) {
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(existing) && !overwrite) {
      stop("Output directory is not empty; use --overwrite=TRUE: ", path, call. = FALSE)
    }
    if (length(existing) && overwrite) unlink(path, recursive = TRUE, force = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

model_contrast_table <- function(contrasts) {
  rows <- lapply(names(contrasts), function(contrast_id) {
    weights <- contrasts[[contrast_id]]
    data.frame(
      contrast_id = contrast_id,
      coefficient = names(weights),
      weight = as.numeric(weights),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

write_origin_stratified_outputs <- function(
  result,
  table_dir,
  figure_dir,
  interval,
  fdr_threshold,
  lfc_threshold
) {
  origin <- result$origin
  for (contrast_id in names(result$de_tables)) {
    write_tsv(
      result$de_tables[[contrast_id]],
      file.path(table_dir, paste0(contrast_id, "_", origin, "_de.tsv"))
    )
  }
  prefix <- paste0("treated_equal_dose_minus_vehicle_", origin)
  write_tsv(
    result$pathways$ranking,
    file.path(table_dir, paste0(prefix, "_gsea_ranking.tsv"))
  )
  write_tsv(
    result$pathways$symbol_resolution,
    file.path(table_dir, paste0(prefix, "_gene_symbol_resolution.tsv"))
  )
  write_tsv(
    result$pathways$complete,
    file.path(table_dir, paste0(prefix, "_gsea_complete.tsv"))
  )
  write_tsv(
    result$pathways$selected,
    file.path(table_dir, paste0(prefix, "_gsea_selected.tsv"))
  )
  write_tsv(
    result$pathways$leading_edge,
    file.path(table_dir, paste0(prefix, "_selected_leading_edge_genes.tsv"))
  )
  write_tsv(
    result$pathways$collection_summary,
    file.path(table_dir, paste0(prefix, "_gsea_summary.tsv"))
  )
  write_tsv(
    result$mouse_metadata,
    file.path(table_dir, paste0("interval_mouse_coverage_", origin, ".tsv"))
  )
  design_output <- cbind(
    sample_id = rownames(result$design$design),
    as.data.frame(result$design$design, check.names = FALSE)
  )
  write_tsv(
    design_output,
    file.path(table_dir, paste0("design_matrix_", origin, ".tsv"))
  )
  write_tsv(
    model_contrast_table(result$contrasts),
    file.path(table_dir, paste0("contrast_definitions_", origin, ".tsv"))
  )
  write_tsv(data.frame(
    origin_stratum = origin,
    cpm_threshold = result$model$cpm_threshold,
    minimum_mice_above_threshold = result$model$minimum_mice,
    rule = sprintf(
      "CPM > %s in at least %d independent mouse pseudobulks",
      result$model$cpm_threshold,
      result$model$minimum_mice
    ),
    rationale = "smallest within-origin dose group contains two mice",
    n_input_features = result$model$n_input_genes,
    n_retained_features = result$model$n_retained_genes,
    stringsAsFactors = FALSE
  ), file.path(table_dir, paste0("expression_filter_audit_", origin, ".tsv")))

  significant <- is.finite(result$primary$fdr) &
    result$primary$fdr <= fdr_threshold &
    is.finite(result$primary$log2_fold_change) &
    abs(result$primary$log2_fold_change) >= lfc_threshold
  write_tsv(data.frame(
    origin_stratum = origin,
    interval_cells = sum(result$mouse_metadata$n_cells),
    interval_mice = nrow(result$mouse_metadata),
    vehicle_mice = sum(result$mouse_metadata$dose_mg == 0),
    treated_mice = sum(result$mouse_metadata$dose_mg > 0),
    dose_30_mice = sum(result$mouse_metadata$dose_mg == 30),
    dose_120_mice = sum(result$mouse_metadata$dose_mg == 120),
    model_formula = paste(deparse(result$design$formula), collapse = " "),
    design_rank = result$design$rank,
    design_columns = ncol(result$design$design),
    injected_origin_adjustment = FALSE,
    mean_within_interval_pseudotime_adjustment = TRUE,
    retained_expressed_features = result$model$n_retained_genes,
    significant_features = sum(significant),
    higher_in_treated = sum(significant & result$primary$log2_fold_change > 0),
    lower_in_treated = sum(significant & result$primary$log2_fold_change < 0),
    gsea_fdr_significant = sum(
      result$pathways$complete$padj <= fdr_threshold,
      na.rm = TRUE
    ),
    displayed_pathways = nrow(result$pathways$selected),
    stringsAsFactors = FALSE
  ), file.path(table_dir, paste0("analysis_summary_", origin, ".tsv")))

  plot_candidate_pathways(
    result$pathways$selected,
    file.path(
      figure_dir,
      paste0(
        "panel_7I_candidate_treated_vs_vehicle_interval_pathways_",
        origin,
        ".pdf"
      )
    ),
    file.path(
      figure_dir,
      paste0(
        "panel_7I_candidate_treated_vs_vehicle_interval_pathways_",
        origin,
        ".png"
      )
    ),
    result$mouse_metadata,
    interval,
    exploratory_pathway_selection_rule(),
    title = paste("Pathway enrichment in", origin, "-origin tumors"),
    adjustment_caption = paste0(
      "Restricted to ", origin,
      "-origin tumors; adjusted for mean within-interval pseudotime."
    )
  )
  invisible(result)
}

write_treatment_origin_interaction_outputs <- function(
  interaction_analysis,
  stratified_analyses,
  table_dir,
  figure_dir,
  interval,
  fdr_threshold
) {
  if (!identical(names(stratified_analyses), c("2N", "4N"))) {
    stop("Formal interaction display requires named 2N and 4N analyses", call. = FALSE)
  }
  for (contrast_id in names(interaction_analysis$de_tables)) {
    write_tsv(
      interaction_analysis$de_tables[[contrast_id]],
      file.path(table_dir, paste0(contrast_id, "_interaction_model_de.tsv"))
    )
  }
  primary_id <- interaction_analysis$contrast_id
  prefix <- primary_id
  write_tsv(
    interaction_analysis$pathways$ranking,
    file.path(table_dir, paste0(prefix, "_gsea_ranking.tsv"))
  )
  write_tsv(
    interaction_analysis$pathways$symbol_resolution,
    file.path(table_dir, paste0(prefix, "_gene_symbol_resolution.tsv"))
  )
  write_tsv(
    interaction_analysis$pathways$complete,
    file.path(table_dir, paste0(prefix, "_gsea_complete.tsv"))
  )
  write_tsv(
    interaction_analysis$pathways$collection_summary,
    file.path(table_dir, paste0(prefix, "_gsea_summary.tsv"))
  )
  for (origin in c("2N", "4N")) {
    origin_prefix <- paste0(
      "treated_equal_dose_minus_vehicle_",
      origin,
      "_joint_interaction_model"
    )
    origin_pathways <- interaction_analysis$origin_effect_pathways[[origin]]
    write_tsv(
      origin_pathways$ranking,
      file.path(table_dir, paste0(origin_prefix, "_gsea_ranking.tsv"))
    )
    write_tsv(
      origin_pathways$complete,
      file.path(table_dir, paste0(origin_prefix, "_gsea_complete.tsv"))
    )
    write_tsv(
      origin_pathways$collection_summary,
      file.path(table_dir, paste0(origin_prefix, "_gsea_summary.tsv"))
    )
  }
  design_output <- cbind(
    sample_id = rownames(interaction_analysis$design$design),
    as.data.frame(interaction_analysis$design$design, check.names = FALSE)
  )
  write_tsv(
    design_output,
    file.path(table_dir, "design_matrix_treatment_by_origin_interaction.tsv")
  )
  write_tsv(
    model_contrast_table(interaction_analysis$contrasts),
    file.path(table_dir, "contrast_definitions_treatment_by_origin_interaction.tsv")
  )

  comparison <- prepare_origin_interaction_pathways(
    stratified_analyses[["2N"]]$pathways$complete,
    stratified_analyses[["4N"]]$pathways$complete,
    interaction_analysis$pathways$complete,
    origin_effect_gsea_2n = interaction_analysis$origin_effect_pathways[["2N"]]$complete,
    origin_effect_gsea_4n = interaction_analysis$origin_effect_pathways[["4N"]]$complete,
    fdr_threshold = fdr_threshold,
    maximum_per_direction = interaction_pathway_max_per_direction()
  )
  write_tsv(
    comparison$shared,
    file.path(
      table_dir,
      "treated_equal_dose_minus_vehicle_shared_significant_pathways_2N_4N.tsv"
    )
  )
  selected_columns <- c(
    "panel_id", "more_positive_origin", "collection", "collection_label", "pathway",
    "pathway_label", "NES", "padj", "NES_2N", "padj_2N", "NES_4N", "padj_4N",
    "origin_effect_nes_order_concordant", "interaction_fdr_significant",
    "selected_rank_within_direction", "selection_scope"
  )
  write_tsv(
    comparison$selected_interactions[, selected_columns, drop = FALSE],
    file.path(table_dir, paste0(prefix, "_selected_for_panel.tsv"))
  )
  write_tsv(data.frame(
    contrast_id = primary_id,
    contrast_definition = paste0(
      "[0.5*(30 mg/kg) + 0.5*(120 mg/kg) - vehicle]_4N - ",
      "[0.5*(30 mg/kg) + 0.5*(120 mg/kg) - vehicle]_2N"
    ),
    positive_gene_level_effect = "treated-minus-vehicle effect more positive in 4N-origin tumors",
    negative_gene_level_effect = "treated-minus-vehicle effect more positive in 2N-origin tumors",
    model_formula = paste(deparse(interaction_analysis$design$formula), collapse = " "),
    n_mouse_pseudobulks = nrow(interaction_analysis$design$design),
    design_rank = interaction_analysis$design$rank,
    design_columns = ncol(interaction_analysis$design$design),
    retained_expressed_features = interaction_analysis$model$n_retained_genes,
    gsea_pathways_fdr_significant = comparison$n_interaction_fdr_significant,
    selected_2N_more_positive = nrow(comparison$selected_2N),
    selected_2N_more_positive_fdr_significant = sum(
      comparison$selected_2N$interaction_fdr_significant
    ),
    selected_4N_more_positive = nrow(comparison$selected_4N),
    selected_4N_more_positive_fdr_significant = sum(
      comparison$selected_4N$interaction_fdr_significant
    ),
    stringsAsFactors = FALSE
  ), file.path(table_dir, "treatment_by_origin_interaction_summary.tsv"))

  plot_origin_interaction_pathways(
    comparison,
    file.path(
      figure_dir,
      "panel_7I_candidate_treated_vs_vehicle_interval_pathways_formal_origin_interaction.pdf"
    ),
    file.path(
      figure_dir,
      "panel_7I_candidate_treated_vs_vehicle_interval_pathways_formal_origin_interaction.png"
    ),
    lapply(stratified_analyses, `[[`, "mouse_metadata"),
    interval
  )
  invisible(comparison)
}

run_interval_de <- function(args, repo_root) {
  check_packages(if (isTRUE(args$global_magnitude_only)) {
    global_magnitude_required_packages()
  } else {
    required_packages()
  })
  support <- load_figure7_support(repo_root)
  config <- support$read_config(args$config)
  interval <- config$interval_list$primary_accumulated_state
  if (!isTRUE(all.equal(interval$start, 0.296, tolerance = 0)) ||
      !isTRUE(all.equal(interval$end, 0.486, tolerance = 0)) ||
      !isTRUE(interval$include_start) || !isTRUE(interval$include_end)) {
    stop("The exploratory analysis requires the reviewed inclusive interval [0.296, 0.486]", call. = FALSE)
  }

  message("Verifying reviewed inputs")
  input_hashes <- c(
    cell_metadata = if (args$verify_input_hashes) {
      verify_input_hash(args$cell_metadata, config$inputs$cellcycle_sha256, "CellCycle metadata")
    } else {
      sha256_file(args$cell_metadata)
    },
    seurat_rds = if (args$verify_input_hashes) {
      verify_input_hash(args$seurat_rds, config$raw_data$seurat_rds_sha256, "Seurat RDS")
    } else {
      sha256_file(args$seurat_rds)
    },
    config = sha256_file(args$config)
  )

  message("Reading exact CellCycle metadata")
  all_meta <- support$read_cell_metadata(args$cell_metadata)
  support$check_pseudotime(all_meta)
  metadata_audit <- support$audit_metadata(all_meta)
  if (nrow(metadata_audit$duplicate_cells) || nrow(metadata_audit$inconsistent_cells)) {
    stop("CellCycle metadata contains duplicate or inconsistent cell records", call. = FALSE)
  }

  message("Loading RNA counts and applying the fail-closed GRCh38-only boundary")
  loaded <- support$load_counts(
    args$seurat_rds,
    assay = args$assay,
    counts_layer = args$counts_layer,
    human_prefix = config$feature_species$human_prefix,
    mouse_prefix = config$feature_species$mouse_prefix
  )
  species_audit <- loaded$species_audit
  matched <- support$match_cells(all_meta, loaded$counts, min_match_rate = 1)
  match_audit <- matched$audit
  rm(loaded)
  invisible(gc())

  interval_meta <- select_interval_cells(matched$meta, interval)
  interval_counts <- matched$counts[, interval_meta$cell_id, drop = FALSE]
  rm(matched)
  invisible(gc())
  pseudobulk <- construct_mouse_pseudobulk(interval_meta, interval_counts)
  validate_interval_contract(all_meta, interval_meta, pseudobulk$metadata)
  if (isTRUE(args$global_magnitude_only)) {
    return(run_global_magnitude_only(
      args = args,
      repo_root = repo_root,
      support = support,
      config = config,
      interval = interval,
      pseudobulk = pseudobulk,
      species_audit = species_audit,
      match_audit = match_audit,
      input_hashes = input_hashes
    ))
  }
  design_bundle <- prepare_interval_design(pseudobulk$metadata, adjust_mean_pseudotime = TRUE)
  contrasts <- treatment_contrasts(design_bundle$design)

  message("Fitting the mouse-level TMM/voom model")
  model <- fit_interval_voom(
    pseudobulk$counts,
    design_bundle,
    cpm_threshold = args$cpm_threshold,
    minimum_mice = args$minimum_mice
  )
  de_tables <- lapply(names(contrasts), function(contrast_id) {
    extract_contrast_table(
      model,
      contrasts[[contrast_id]],
      contrast_id,
      support$clean_gene_symbols
    )
  })
  names(de_tables) <- names(contrasts)
  primary_id <- "treated_equal_dose_minus_vehicle"
  primary <- de_tables[[primary_id]]
  sensitivity_design <- prepare_interval_design(
    pseudobulk$metadata,
    adjust_mean_pseudotime = FALSE
  )
  sensitivity_contrasts <- treatment_contrasts(sensitivity_design$design)
  sensitivity_model <- fit_interval_voom(
    pseudobulk$counts,
    sensitivity_design,
    cpm_threshold = args$cpm_threshold,
    minimum_mice = args$minimum_mice
  )
  if (!identical(rownames(model$dge$counts), rownames(sensitivity_model$dge$counts))) {
    stop("Primary and sensitivity models retained different gene universes", call. = FALSE)
  }
  sensitivity_primary <- extract_contrast_table(
    sensitivity_model,
    sensitivity_contrasts[[primary_id]],
    paste0(primary_id, "_without_mean_pseudotime_adjustment"),
    support$clean_gene_symbols
  )

  message("Running the reused human Hallmark/Reactome/GO-BP GSEA")
  pathway_enrichment <- run_pathway_enrichment(
    primary,
    support,
    config,
    seed = 1L
  )

  message("Fitting the origin-stratified mouse-level models and GSEA")
  origins <- c("2N", "4N")
  stratified_analyses <- setNames(lapply(origins, function(origin) {
    fit_origin_stratified_analysis(
      pseudobulk,
      origin,
      support,
      config,
      gene_sets = pathway_enrichment$gene_sets,
      cpm_threshold = args$cpm_threshold,
      seed = 1L
    )
  }), origins)

  message("Fitting the formal treatment-by-origin interaction model and GSEA")
  interaction_analysis <- fit_treatment_origin_interaction_analysis(
    pseudobulk,
    support,
    config,
    gene_sets = pathway_enrichment$gene_sets,
    cpm_threshold = args$cpm_threshold,
    minimum_mice = args$minimum_mice,
    seed = 1L
  )

  message("Comparing global treatment-response magnitude between origins")
  global_response_magnitude <- analyze_global_response_magnitude(
    interaction_analysis,
    replicates = args$magnitude_bootstrap_replicates,
    seed = args$magnitude_bootstrap_seed
  )

  output_root <- prepare_output_root(args$output_root, overwrite = args$overwrite)
  table_dir <- file.path(output_root, "tables")
  figure_dir <- file.path(output_root, "figures")
  metadata_dir <- file.path(output_root, "metadata")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  for (contrast_id in names(de_tables)) {
    write_tsv(
      de_tables[[contrast_id]],
      file.path(table_dir, paste0(contrast_id, "_de.tsv"))
    )
  }
  write_tsv(
    sensitivity_primary,
    file.path(
      table_dir,
      "treated_equal_dose_minus_vehicle_sensitivity_without_mean_pseudotime_adjustment_de.tsv"
    )
  )
  write_tsv(
    pathway_enrichment$ranking,
    file.path(table_dir, "treated_equal_dose_minus_vehicle_gsea_ranking.tsv")
  )
  write_tsv(
    pathway_enrichment$symbol_resolution,
    file.path(table_dir, "treated_equal_dose_minus_vehicle_gene_symbol_resolution.tsv")
  )
  write_tsv(
    pathway_enrichment$complete,
    file.path(table_dir, "treated_equal_dose_minus_vehicle_gsea_complete.tsv")
  )
  write_tsv(
    pathway_enrichment$selected,
    file.path(table_dir, "treated_equal_dose_minus_vehicle_gsea_selected.tsv")
  )
  write_tsv(
    pathway_enrichment$leading_edge,
    file.path(table_dir, "treated_equal_dose_minus_vehicle_selected_leading_edge_genes.tsv")
  )
  write_tsv(
    pathway_enrichment$collection_summary,
    file.path(table_dir, "treated_equal_dose_minus_vehicle_gsea_summary.tsv")
  )
  write_tsv(pseudobulk$metadata, file.path(table_dir, "interval_mouse_coverage.tsv"))
  design_output <- cbind(
    sample_id = rownames(design_bundle$design),
    as.data.frame(design_bundle$design, check.names = FALSE)
  )
  write_tsv(design_output, file.path(table_dir, "design_matrix.tsv"))
  write_tsv(model_contrast_table(contrasts), file.path(table_dir, "contrast_definitions.tsv"))

  significant <- is.finite(primary$fdr) & primary$fdr <= args$fdr_threshold &
    is.finite(primary$log2_fold_change) & abs(primary$log2_fold_change) >= args$lfc_threshold
  significant_features <- primary$feature_id[significant]
  sensitivity_significant <- is.finite(sensitivity_primary$fdr) &
    sensitivity_primary$fdr <= args$fdr_threshold &
    is.finite(sensitivity_primary$log2_fold_change) &
    abs(sensitivity_primary$log2_fold_change) >= args$lfc_threshold
  sensitivity_significant_features <- sensitivity_primary$feature_id[sensitivity_significant]
  sensitivity_effect_for_primary <- sensitivity_primary$log2_fold_change[
    match(significant_features, sensitivity_primary$feature_id)
  ]
  write_tsv(
    dose_concordance_table(
      primary,
      de_tables$dose_30_minus_vehicle,
      de_tables$dose_120_minus_vehicle,
      significant_features
    ),
    file.path(table_dir, "primary_significant_gene_dose_concordance.tsv")
  )
  write_tsv(
    mouse_expression_audit(
      model,
      significant_features,
      primary,
      pseudobulk$metadata
    ),
    file.path(table_dir, "primary_significant_gene_mouse_expression.tsv")
  )
  summary <- data.frame(
    metric = c(
      "full_cellcycle_cells", "interval_cells", "interval_mice",
      "vehicle_mice", "treated_mice", "dose_30_mice", "dose_120_mice",
      "input_human_features", "retained_expressed_features",
      "significant_features", "higher_in_treated", "lower_in_treated",
      "sensitivity_significant_features_without_mean_pseudotime_adjustment",
      "significant_features_shared_with_sensitivity",
      "primary_significant_features_direction_concordant_in_sensitivity"
    ),
    value = c(
      nrow(all_meta), nrow(interval_meta), nrow(pseudobulk$metadata),
      sum(pseudobulk$metadata$dose_mg == 0),
      sum(pseudobulk$metadata$dose_mg > 0),
      sum(pseudobulk$metadata$dose_mg == 30),
      sum(pseudobulk$metadata$dose_mg == 120),
      model$n_input_genes,
      model$n_retained_genes,
      sum(significant),
      sum(significant & primary$log2_fold_change > 0),
      sum(significant & primary$log2_fold_change < 0),
      sum(sensitivity_significant),
      length(intersect(significant_features, sensitivity_significant_features)),
      sum(
        sign(primary$log2_fold_change[significant]) ==
          sign(sensitivity_effect_for_primary),
        na.rm = TRUE
      )
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(summary, file.path(table_dir, "analysis_summary.tsv"))

  plot_candidate_volcano(
    primary,
    file.path(figure_dir, "panel_7I_candidate_treated_vs_vehicle_interval_de.pdf"),
    file.path(figure_dir, "panel_7I_candidate_treated_vs_vehicle_interval_de.png"),
    pseudobulk$metadata,
    interval,
    args$fdr_threshold,
    args$lfc_threshold,
    args$label_genes
  )
  plot_candidate_pathways(
    pathway_enrichment$selected,
    file.path(
      figure_dir,
      "panel_7I_candidate_treated_vs_vehicle_interval_pathways.pdf"
    ),
    file.path(
      figure_dir,
      "panel_7I_candidate_treated_vs_vehicle_interval_pathways.png"
    ),
    pseudobulk$metadata,
    interval,
    exploratory_pathway_selection_rule()
  )
  for (origin in names(stratified_analyses)) {
    write_origin_stratified_outputs(
      stratified_analyses[[origin]],
      table_dir,
      figure_dir,
      interval,
      args$fdr_threshold,
      args$lfc_threshold
    )
  }
  origin_comparison <- write_treatment_origin_interaction_outputs(
    interaction_analysis,
    stratified_analyses,
    table_dir,
    figure_dir,
    interval,
    args$fdr_threshold
  )
  write_global_response_magnitude_outputs(global_response_magnitude, table_dir)

  write_tsv(species_audit, file.path(metadata_dir, "feature_species_audit.tsv"))
  write_tsv(match_audit, file.path(metadata_dir, "cell_expression_match_audit.tsv"))
  membership_path <- file.path(metadata_dir, "gene_set_membership.tsv")
  write_tsv(pathway_enrichment$membership, membership_path)
  pathway_enrichment$contract$value[
    pathway_enrichment$contract$key == "gene_set_membership_sha256"
  ] <- sha256_file(membership_path)
  write_tsv(
    pathway_enrichment$contract,
    file.path(metadata_dir, "gene_set_contract.tsv")
  )
  design_audit_row <- function(model_id, bundle, primary_model, adjust_pseudotime) {
    data.frame(
      model_id = model_id,
      primary_model = primary_model,
      formula = paste(deparse(bundle$formula), collapse = " "),
      n_observations = nrow(bundle$design),
      n_coefficients = ncol(bundle$design),
      design_rank = bundle$rank,
      design_condition_number = kappa(bundle$design),
      full_rank = bundle$rank == ncol(bundle$design),
      experimental_unit = "mouse",
      pseudobulk_unit = "one library per mouse within exact interval",
      primary_contrast = primary_id,
      primary_contrast_definition = "0.5*(30 mg/kg) + 0.5*(120 mg/kg) - vehicle",
      injected_origin_adjustment = TRUE,
      mean_within_interval_pseudotime_adjustment = adjust_pseudotime,
      treated_mean_pseudotime_correlation = stats::cor(
        as.numeric(pseudobulk$metadata$treatment == "treated"),
        pseudobulk$metadata$mean_pseudotime
      ),
      duplicate_correlation = FALSE,
      rationale_duplicate_correlation = "one independent pseudobulk observation per mouse",
      stringsAsFactors = FALSE
    )
  }
  write_tsv(rbind(
    design_audit_row(
      "mouse_pseudobulk_treated_equal_dose_minus_vehicle_v1",
      design_bundle,
      TRUE,
      TRUE
    ),
    design_audit_row(
      "sensitivity_without_mean_pseudotime_adjustment",
      sensitivity_design,
      FALSE,
      FALSE
    )
  ), file.path(metadata_dir, "design_audit.tsv"))
  write_tsv(data.frame(
    interval_id = interval$name,
    interval_start = interval$start,
    interval_end = interval$end,
    include_start = interval$include_start,
    include_end = interval$include_end,
    interval_source = interval$source,
    interval_selection_caveat = paste(
      "The interval was localized using treated-versus-vehicle density differences;",
      "the within-interval differential-expression analysis is exploratory and conditional on that selection."
    ),
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "interval_definition.tsv"))
  write_tsv(data.frame(
    cpm_threshold = model$cpm_threshold,
    minimum_mice_above_threshold = model$minimum_mice,
    rule = sprintf("CPM > %s in at least %d independent mouse pseudobulks", model$cpm_threshold, model$minimum_mice),
    n_input_features = model$n_input_genes,
    n_retained_features = model$n_retained_genes,
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "expression_filter_audit.tsv"))

  script <- exploratory_script_path()
  manifest_paths <- c(
    cell_metadata = args$cell_metadata,
    seurat_rds = args$seurat_rds,
    config = args$config,
    exploratory_script = script,
    state_pathway_support = file.path(
      repo_root,
      "Code/in-vivo/figure7/generate_pseudotime_state_pathways_support.R"
    ),
    feature_species_policy = file.path(
      repo_root,
      "Code/in-vivo/figure7/src/feature_species_policy.R"
    ),
    shared_figure7_contract = file.path(
      repo_root,
      "Code/in-vivo/figure7/src/common_io.R"
    )
  )
  manifest_hashes <- setNames(rep(NA_character_, length(manifest_paths)), names(manifest_paths))
  manifest_hashes[names(input_hashes)] <- input_hashes
  hashes_to_compute <- names(manifest_hashes)[is.na(manifest_hashes)]
  manifest_hashes[hashes_to_compute] <- vapply(
    manifest_paths[hashes_to_compute],
    sha256_file,
    character(1L)
  )
  manifest <- data.frame(
    input_id = names(manifest_paths),
    locator = vapply(manifest_paths, portable_locator, character(1L), repo_root = repo_root),
    sha256 = unname(manifest_hashes),
    byte_size = as.numeric(file.info(manifest_paths)$size),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(metadata_dir, "input_manifest.tsv"))
  write_tsv(data.frame(
    package = required_packages(),
    version = vapply(required_packages(), function(package) {
      as.character(utils::packageVersion(package))
    }, character(1L)),
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "package_versions.tsv"))
  write_tsv(data.frame(
    key = c(
      "git_revision", "git_branch", "git_worktree_dirty",
      "feature_species_policy", "normalization", "observation_model",
      "multiple_testing", "contrast_sign_convention", "pathway_ranking",
      "pathway_collections", "pathway_selection",
      "treatment_by_origin_interaction", "global_response_magnitude"
    ),
    value = c(
      git_value(repo_root, c("rev-parse", "HEAD")),
      git_value(repo_root, c("branch", "--show-current")),
      as.character(nzchar(git_value(
        repo_root,
        c("status", "--porcelain", "--untracked-files=no")
      ) %||% "")),
      support$figure7_human_feature_policy_description(),
      "edgeR TMM",
      "limma voom with robust empirical Bayes moderation",
      "Benjamini-Hochberg across retained human features within each contrast",
      "positive log2 fold change means higher expression in gemcitabine-treated tumors",
      "moderated t statistic from each stated gene-level contrast",
      paste(as.character(unlist(config$state_pathways$collections)), collapse = ","),
      paste0(
        exploratory_pathway_selection_rule(),
        " for pooled/stratified panels; formal comparison uses the 10 most ",
        "negative and 10 most positive interaction NES across collections"
      ),
      paste0(
        "formal contrast: [equal-dose treated - vehicle]_4N - ",
        "[equal-dose treated - vehicle]_2N; positive values are more positive in 4N"
      ),
      paste0(
        "noise-corrected mean squared log2 fold change on the joint-model gene universe; ",
        "2N minus 4N with HC2 whole-mouse Rademacher wild-bootstrap uncertainty"
      )
    ),
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "provenance.tsv"))
  session_info <- sub(
    "[[:space:]]+$",
    "",
    capture.output(utils::sessionInfo())
  )
  writeLines(session_info, file.path(metadata_dir, "sessionInfo.txt"))

  message("Completed exploratory interval gene and pathway analysis: ", output_root)
  invisible(list(
    primary = primary,
    pathways = pathway_enrichment,
    stratified = stratified_analyses,
    origin_comparison = origin_comparison,
    interaction = interaction_analysis,
    global_response_magnitude = global_response_magnitude,
    summaries = summary,
    mouse_metadata = pseudobulk$metadata,
    output_root = output_root
  ))
}

main <- function() {
  repo_root <- exploratory_repo_root()
  defaults <- list(
    cell_metadata = file.path(
      repo_root,
      "Data/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    ),
    seurat_rds = file.path(
      repo_root,
      "Results/in-vivo/figure7/raw/zenodo_21463392/integrated_sct_cca_seurat_final_reclustered.rds"
    ),
    config = file.path(repo_root, "Code/in-vivo/figure7/figure7_config.yaml"),
    output_root = file.path(
      repo_root,
      "Results/in-vivo/figure7/exploratory/treated_vs_vehicle_interval_de"
    ),
    assay = "RNA",
    counts_layer = "counts",
    cpm_threshold = 1,
    minimum_mice = 4L,
    fdr_threshold = 0.05,
    lfc_threshold = 0.5,
    label_genes = 12L,
    magnitude_bootstrap_replicates = 4999L,
    magnitude_bootstrap_seed = 20260818L,
    global_magnitude_only = FALSE,
    verify_input_hashes = TRUE,
    overwrite = FALSE
  )
  argv <- commandArgs(trailingOnly = TRUE)
  if (any(argv %in% c("--help", "-h"))) {
    usage()
    return(invisible(NULL))
  }
  args <- parse_cli_args(argv, defaults)
  for (field in c("cell_metadata", "seurat_rds", "config")) {
    args[[field]] <- resolve_path(args[[field]], repo_root, must_work = TRUE)
  }
  args$output_root <- resolve_path(args$output_root, repo_root, must_work = FALSE)
  run_interval_de(args, repo_root)
}

if (identical(environment(), globalenv())) main()
