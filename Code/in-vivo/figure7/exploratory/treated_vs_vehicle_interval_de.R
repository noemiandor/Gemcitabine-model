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

prepare_interval_design <- function(mouse_meta, adjust_mean_pseudotime = TRUE) {
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
  formula <- if (isTRUE(adjust_mean_pseudotime)) {
    ~ 0 + dose_group + initial_ploidy_factor + mean_pseudotime_z
  } else {
    ~ 0 + dose_group + initial_ploidy_factor
  }
  design <- stats::model.matrix(formula, data = data)
  rownames(design) <- data$sample_id
  rank <- qr(design)$rank
  if (rank != ncol(design)) {
    stop("Interval differential-expression design is rank deficient", call. = FALSE)
  }
  list(data = data, design = design, formula = formula, rank = rank)
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

extract_contrast_table <- function(model, contrast, contrast_id, clean_symbols) {
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
      "higher_in_treated",
      ifelse(table$logFC < 0, "lower_in_treated", "no_change")
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

run_pathway_enrichment <- function(primary, support, config, seed = 1L) {
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
  gene_sets <- support$fetch_gene_sets(
    paste(collections, collapse = ","),
    species = species
  )
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
    "treated_equal_dose_minus_vehicle",
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
  membership <- support$gene_set_membership_table(gene_sets)
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

plot_candidate_pathways <- function(
  selected,
  path_pdf,
  path_png,
  mouse_meta,
  interval,
  selection_rule
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
      title = "Pathway enrichment in the selected pseudotime interval",
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
        "Adjusted for injected origin and mean within-interval pseudotime.\n",
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

run_interval_de <- function(args, repo_root) {
  check_packages()
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
      "pathway_collections", "pathway_selection"
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
      "moderated t statistic from the primary treated-minus-vehicle contrast",
      paste(as.character(unlist(config$state_pathways$collections)), collapse = ","),
      exploratory_pathway_selection_rule()
    ),
    stringsAsFactors = FALSE
  ), file.path(metadata_dir, "provenance.tsv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(metadata_dir, "sessionInfo.txt"))

  message("Completed exploratory interval gene and pathway analysis: ", output_root)
  invisible(list(
    primary = primary,
    pathways = pathway_enrichment,
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
