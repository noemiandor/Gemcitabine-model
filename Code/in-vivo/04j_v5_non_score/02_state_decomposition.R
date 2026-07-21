ptpv5_group_summary <- function(Y, available, group_mask) {
  group_mask <- as.matrix(group_mask)
  available <- as.numeric(available)
  W <- sweep(group_mask, 1L, available, "*")
  n <- colSums(W)
  Y0 <- Y
  Y0[!is.finite(Y0)] <- 0
  mean <- sweep(Y0 %*% W, 2L, pmax(n, 1), "/")
  mean[, n < 1L] <- NA_real_
  list(mean = mean, n = n)
}

ptpv5_mouse_interaction_matrix <- function(Y, mouse_meta, A, min_group = 1L) {
  p2 <- mouse_meta$initial_ploidy == "2N"
  p4 <- mouse_meta$initial_ploidy == "4N"
  available <- is.finite(colSums(Y))
  one <- function(mask) {
    treated <- ptpv5_group_summary(Y, available & mask, A)
    control <- ptpv5_group_summary(Y, available & mask, !A)
    delta <- treated$mean - control$mean
    valid <- treated$n >= min_group & control$n >= min_group
    delta[, !valid] <- NA_real_
    list(
      delta = delta,
      treated_mean = treated$mean,
      control_mean = control$mean,
      treated_n = treated$n,
      control_n = control$n,
      valid = valid
    )
  }
  z2 <- one(p2)
  z4 <- one(p4)
  list(
    interaction = z4$delta - z2$delta,
    ploidy2 = z2,
    ploidy4 = z4
  )
}

ptpv5_standardize_rows_by_randomization <- function(X) {
  center <- matrixStats::rowMedians(X, na.rm = TRUE)
  scale <- matrixStats::rowSds(X, na.rm = TRUE)
  fallback <- stats::median(scale[is.finite(scale) & scale > 0], na.rm = TRUE)
  if (!is.finite(fallback) || fallback <= 0) fallback <- 1
  scale[!is.finite(scale) | scale <= 0] <- fallback
  Z <- sweep(X, 1L, center, "-")
  sweep(Z, 1L, scale, "/")
}

ptpv5_vector_stat_matrix <- function(Z) {
  Z <- as.matrix(Z)
  finite <- is.finite(Z)
  denom <- colSums(finite)
  denom[denom < 1L] <- NA_real_
  Z0 <- Z
  Z0[!finite] <- 0
  pos <- colSums(pmax(Z0, 0)) / denom
  neg <- colSums(pmax(-Z0, 0)) / denom
  maxabs <- matrixStats::colMaxs(abs(Z), na.rm = TRUE)
  maxabs[!is.finite(maxabs)] <- NA_real_
  rbind(
    true_maxmean = pmax(pos, neg),
    mean_square = colSums(Z0^2) / denom,
    max_absolute = maxabs
  )
}

ptpv5_load_state_expression <- function(context, needed_symbols) {
  cfg <- context$cfg
  dirs <- context$dirs
  cache_path <- file.path(dirs$checkpoints, "state_cell_expression_cache.rds")
  cache_key <- ptpv5_object_sha(list(
    schema = "state_cell_expression_cache_v3_family_gene_tmm_counts_and_cpm",
    seurat = ptpv5_file_sha(cfg$inputs$seurat_rds),
    cell_metadata = ptpv5_file_sha(cfg$inputs$cell_metadata),
    genes = sort(needed_symbols)
  ))
  if (file.exists(cache_path)) {
    x <- readRDS(cache_path)
    if (identical(x$cache_key, cache_key)) return(x$data)
  }
  cells <- read.csv(
    cfg$inputs$cell_metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  cells$cell_id <- as.character(cells$cell_id)
  cells$sample_id <- as.character(cells$sample_id)
  cells$cluster <- as.character(cells$cluster)
  cells$initial_ploidy <- as.character(cells$initial_ploidy)
  cells$treatment <- ifelse(
    suppressWarnings(as.numeric(cells$gemcitabine_dose_mg_per_kg)) > 0,
    "gemcitabine",
    "control"
  )
  cells <- cells[
    is.finite(cells$pseudotime) &
      cells$pseudotime >= 0 &
      cells$pseudotime <= 1,
    ,
    drop = FALSE
  ]
  if (!nrow(cells)) {
    stop("No finite target cells remain for state analyses.", call. = FALSE)
  }
  ptpv5_message(
    "Loading Seurat object for state-specific expression (one-time step).",
    log_file = context$log_file
  )
  obj <- readRDS(cfg$inputs$seurat_rds)
  md <- obj[[]]
  md$cell_id <- rownames(md)
  phase_col <- intersect(c("Phase", "phase", "cell_cycle_phase"), names(md))
  if (!length(phase_col)) {
    stop("The frozen Seurat object has no Phase metadata column.", call. = FALSE)
  }
  cells$Phase <- as.character(md[[phase_col[[1L]]]][
    match(cells$cell_id, md$cell_id)
  ])
  if (anyNA(cells$Phase)) {
    stop("Some frozen primary cells lack Seurat Phase annotations.", call. = FALSE)
  }
  matched_cells <- intersect(cells$cell_id, colnames(obj))
  if (length(matched_cells) != nrow(cells)) {
    missing <- setdiff(cells$cell_id, matched_cells)
    stop(
      "Primary cells missing from Seurat object: ",
      paste(head(missing, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  counts <- try(
    SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"),
    silent = TRUE
  )
  if (inherits(counts, "try-error")) {
    counts <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
  }
  counts <- counts[, cells$cell_id, drop = FALSE]
  clean_features <- ptp_clean_gene_symbols(rownames(counts))
  feature_first <- !duplicated(clean_features)
  feature_lookup <- setNames(which(feature_first), clean_features[feature_first])
  gene_row <- unname(feature_lookup[needed_symbols])
  keep_gene <- is.finite(gene_row)
  if (!any(keep_gene)) {
    stop("No frozen response-family genes map to Seurat RNA counts.", call. = FALSE)
  }
  needed_symbols <- needed_symbols[keep_gene]
  gene_row <- gene_row[keep_gene]
  lib <- Matrix::colSums(counts)
  selected <- counts[gene_row, , drop = FALSE]
  rownames(selected) <- needed_symbols
  dge <- edgeR::DGEList(counts = selected, lib.size = lib)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  effective_lib <- dge$samples$lib.size * dge$samples$norm.factors
  cpm <- sweep(as.matrix(selected), 2L, effective_lib / 1e6, "/")
  selected_counts <- selected
  rm(obj, counts, selected, dge)
  invisible(gc())
  data <- list(
    expression = cpm,
    counts = selected_counts,
    cells = cells,
    gene_symbols = needed_symbols,
    normalization = data.frame(
      cell_id = cells$cell_id,
      library_size = lib,
      effective_library_size = effective_lib,
      stringsAsFactors = FALSE
    )
  )
  saveRDS(
    list(cache_key = cache_key, data = data),
    cache_path,
    compress = "xz"
  )
  data
}

ptpv5_build_state_partitions <- function(cells, cfg) {
  cutpoint <- as.numeric(cfg$state_decomposition$primary_subbin_break)
  cells$pseudotime_subbin <- ifelse(
    cells$pseudotime < cutpoint,
    "pt030_0395",
    "pt0395_049"
  )
  cells$phase_state <- ifelse(
    cells$Phase %in% c("S", "G2M", "G1"),
    cells$Phase,
    "other"
  )
  cells$cluster_state <- as.character(cells$cluster)
  cells$cluster_phase_state <- paste(
    cells$cluster_state, cells$phase_state, sep = "__"
  )
  list(
    cells = cells,
    definitions = list(
      pseudotime_subbin = sort(unique(cells$pseudotime_subbin)),
      cluster = intersect(
        as.character(cfg$state_decomposition$cluster_levels),
        sort(unique(cells$cluster_state))
      ),
      phase = intersect(
        c(
          as.character(cfg$state_decomposition$phase_inference_levels),
          as.character(cfg$state_decomposition$phase_descriptive_levels)
        ),
        sort(unique(cells$phase_state))
      ),
      cluster_phase = sort(unique(cells$cluster_phase_state))
    ),
    columns = c(
      pseudotime_subbin = "pseudotime_subbin",
      cluster = "cluster_state",
      phase = "phase_state",
      cluster_phase = "cluster_phase_state"
    )
  )
}

ptpv5_state_mouse_tables <- function(
  expression,
  cells,
  state_col,
  state_levels,
  sample_ids,
  min_cells
) {
  n_genes <- nrow(expression)
  counts <- matrix(
    0L,
    nrow = length(state_levels),
    ncol = length(sample_ids),
    dimnames = list(state_levels, sample_ids)
  )
  means <- setNames(vector("list", length(state_levels)), state_levels)
  for (state in state_levels) {
    M <- matrix(
      NA_real_,
      nrow = n_genes,
      ncol = length(sample_ids),
      dimnames = list(rownames(expression), sample_ids)
    )
    for (sample in sample_ids) {
      idx <- which(
        cells$sample_id == sample &
          as.character(cells[[state_col]]) == state
      )
      counts[state, sample] <- length(idx)
      if (length(idx) >= min_cells) {
        M[, sample] <- rowMeans(expression[, idx, drop = FALSE])
      }
    }
    means[[state]] <- M
  }
  total <- matrix(
    NA_real_,
    nrow = n_genes,
    ncol = length(sample_ids),
    dimnames = list(rownames(expression), sample_ids)
  )
  for (sample in sample_ids) {
    idx <- which(cells$sample_id == sample)
    if (length(idx) >= min_cells) {
      total[, sample] <- rowMeans(expression[, idx, drop = FALSE])
    }
  }
  list(counts = counts, means = means, total = total)
}

ptpv5_support_gate_state_tables <- function(
  tables,
  mouse_meta,
  min_cells,
  min_mice
) {
  group <- paste(mouse_meta$initial_ploidy, mouse_meta$treatment, sep = "__")
  groups <- unique(group)
  state_pass <- vapply(seq_len(nrow(tables$counts)), function(i) {
    all(vapply(groups, function(g) {
      sum(tables$counts[i, group == g] >= min_cells) >= min_mice
    }, logical(1L)))
  }, logical(1L))
  counts <- tables$counts[state_pass, , drop = FALSE]
  means <- tables$means[state_pass]
  total <- matrix(
    NA_real_,
    nrow = nrow(tables$total),
    ncol = ncol(tables$total),
    dimnames = dimnames(tables$total)
  )
  if (nrow(counts)) {
    for (j in seq_len(ncol(counts))) {
      numerator <- rep(0, nrow(total))
      denominator <- 0
      for (s in seq_len(nrow(counts))) {
        y <- means[[s]][, j]
        weight <- counts[s, j]
        if (weight >= min_cells && all(is.finite(y))) {
          numerator <- numerator + weight * y
          denominator <- denominator + weight
        }
      }
      if (denominator > 0) total[, j] <- numerator / denominator
    }
  }
  list(
    tables = list(counts = counts, means = means, total = total),
    audit = data.frame(
      state = rownames(tables$counts),
      expression_decomposition_support_pass = state_pass,
      minimum_cells_per_mouse_state = min_cells,
      minimum_mice_per_observed_group = min_mice,
      stringsAsFactors = FALSE
    )
  )
}

ptpv5_state_da <- function(counts, mouse_meta, A, pseudocount, observed_index) {
  proportions <- sweep(counts, 2L, pmax(colSums(counts), 1), "/")
  log_prop <- log(counts + pseudocount)
  clr <- sweep(log_prop, 2L, colMeans(log_prop), "-")
  interaction <- ptpv5_mouse_interaction_matrix(
    clr, mouse_meta, A, min_group = 2L
  )$interaction
  Z <- ptpv5_standardize_rows_by_randomization(interaction)
  global <- colSums(Z^2, na.rm = TRUE)
  state_p <- vapply(seq_len(nrow(Z)), function(i) {
    ptpv5_finite_p(abs(Z[i, observed_index]), abs(Z[i, ]))
  }, numeric(1L))
  list(
    clr = clr,
    proportions = proportions,
    interaction = interaction,
    global_statistic = global,
    global_p = ptpv5_finite_p(
      global[[observed_index]], global
    ),
    state_results = data.frame(
      state = rownames(counts),
      observed_clr_interaction = interaction[, observed_index],
      exact_p = state_p,
      stringsAsFactors = FALSE
    )
  )
}

ptpv5_decompose_family <- function(
  gene_symbols,
  expression,
  state_tables,
  mouse_meta,
  A,
  observed_index,
  min_mice
) {
  idx <- match(gene_symbols, rownames(expression))
  idx <- idx[is.finite(idx)]
  if (length(idx) < 2L) return(NULL)
  total <- state_tables$total[idx, , drop = FALSE]
  marginal <- ptpv5_mouse_interaction_matrix(
    total, mouse_meta, A, min_group = min_mice
  )$interaction
  B <- ncol(A)
  within_p <- list(
    `2N` = matrix(0, nrow = length(idx), ncol = B),
    `4N` = matrix(0, nrow = length(idx), ncol = B)
  )
  composition_p <- within_p
  pooled_fraction <- rowSums(state_tables$counts) /
    sum(state_tables$counts)
  pooled_fraction[!is.finite(pooled_fraction)] <- 0
  for (s in seq_along(state_tables$means)) {
    Y <- state_tables$means[[s]][idx, , drop = FALSE]
    z <- ptpv5_mouse_interaction_matrix(
      Y, mouse_meta, A, min_group = min_mice
    )
    pi <- state_tables$counts[s, ] /
      pmax(colSums(state_tables$counts), 1)
    pi_mat <- matrix(
      pi,
      nrow = 1L,
      dimnames = list("state_fraction", colnames(state_tables$counts))
    )
    piz <- ptpv5_mouse_interaction_matrix(
      pi_mat, mouse_meta, A, min_group = min_mice
    )
    for (ploidy in c("2N", "4N")) {
      zz <- if (ploidy == "2N") z$ploidy2 else z$ploidy4
      pp <- if (ploidy == "2N") piz$ploidy2 else piz$ploidy4
      valid <- zz$valid & pp$valid
      within_piece <- pooled_fraction[[s]] * zz$delta
      midpoint <- 0.5 * (zz$treated_mean + zz$control_mean)
      composition_piece <- sweep(
        midpoint,
        2L,
        as.numeric(pp$delta),
        "*"
      )
      within_piece[, !valid] <- NA_real_
      composition_piece[, !valid] <- NA_real_
      within_p[[ploidy]] <- within_p[[ploidy]] + within_piece
      composition_p[[ploidy]] <- composition_p[[ploidy]] +
        composition_piece
    }
  }
  within <- within_p[["4N"]] - within_p[["2N"]]
  composition <- composition_p[["4N"]] - composition_p[["2N"]]
  residual <- marginal - within - composition
  components <- list(
    marginal_total = marginal,
    within_state = within,
    composition = composition,
    decomposition_residual = residual
  )
  distributions <- list()
  observed_rows <- list()
  for (component_id in names(components)) {
    Z <- ptpv5_standardize_rows_by_randomization(components[[component_id]])
    stat <- t(ptpv5_vector_stat_matrix(Z))
    estimable <- all(is.finite(stat[observed_index, ])) &&
      sum(rowSums(is.finite(stat)) == ncol(stat)) >= 2L
    if (estimable) {
      nested <- ptpv5_nested_component_test(stat, observed_index)
      distributions[[component_id]] <- list(
        statistic = stat,
        adaptive_assignment_p = nested$adaptive_assignment_p
      )
      observed_p <- nested$observed_adaptive_p
    } else {
      distributions[[component_id]] <- list(
        statistic = stat,
        adaptive_assignment_p = rep(NA_real_, nrow(stat))
      )
      observed_p <- NA_real_
    }
    observed_rows[[component_id]] <- data.frame(
      decomposition_component = component_id,
      estimability_status = if (estimable) "estimable" else "not_estimable",
      n_genes = nrow(Z),
      observed_true_maxmean = stat[observed_index, "true_maxmean"],
      observed_mean_square = stat[observed_index, "mean_square"],
      observed_max_absolute = stat[observed_index, "max_absolute"],
      adaptive_exact_p = observed_p,
      stringsAsFactors = FALSE
    )
  }
  list(
    observed = ptpv5_bind_rows(observed_rows),
    distributions = distributions
  )
}

ptpv5_run_state_decomposition <- function(context) {
  cfg <- context$cfg
  args <- context$args
  dirs <- context$dirs
  bundle <- context$bundle
  checkpoint_id <- "method02_state_decomposition"
  cached <- ptpv5_load_checkpoint(
    dirs, checkpoint_id, cfg, context$config_sha, context$input_sha
  )
  if (!is.null(cached)) return(cached)
  sets <- ptpv5_frozen_sets(bundle, args)
  family_genes <- sort(unique(unlist(lapply(
    sets$families$union_genes,
    ptpv5_split_genes
  ))))
  if (isTRUE(args$smoke)) {
    family_genes <- head(family_genes, as.integer(args$smoke_genes))
  }
  state_data <- ptpv5_load_state_expression(context, family_genes)
  ptp_write_csv(
    state_data$normalization,
    file.path(dirs$state, "cell_level_tmm_normalization.csv")
  )
  interval <- as.numeric(cfg$state_decomposition$primary_interval)
  primary_keep <- state_data$cells$pseudotime >= interval[[1L]] &
    state_data$cells$pseudotime < interval[[2L]]
  primary_cells <- state_data$cells[primary_keep, , drop = FALSE]
  partition <- ptpv5_build_state_partitions(primary_cells, cfg)
  cells <- partition$cells
  expression <- state_data$expression[
    ,
    match(cells$cell_id, state_data$cells$cell_id),
    drop = FALSE
  ]
  mouse_meta <- bundle$mouse$metadata[
    match(bundle$sample_ids, bundle$mouse$metadata$sample_id),
    ,
    drop = FALSE
  ]
  A <- ptpv5_assignment_matrix(bundle$assignments, bundle$sample_ids)
  min_cells <- as.integer(cfg$state_decomposition$minimum_cells_per_mouse_state)
  min_mice <- as.integer(cfg$state_decomposition$minimum_mice_per_group)
  da_rows <- list()
  da_state_rows <- list()
  support_rows <- list()
  decomposition_rows <- list()
  distributions <- list()
  family_meta <- sets$families
  for (partition_id in names(partition$definitions)) {
    levels <- partition$definitions[[partition_id]]
    if (length(levels) < 2L) next
    tables <- ptpv5_state_mouse_tables(
      expression,
      cells,
      partition$columns[[partition_id]],
      levels,
      bundle$sample_ids,
      min_cells
    )
    da <- ptpv5_state_da(
      tables$counts,
      mouse_meta,
      A,
      as.numeric(cfg$state_decomposition$clr_count_replacement),
      bundle$observed_index
    )
    da_rows[[partition_id]] <- data.frame(
      partition = partition_id,
      n_states = nrow(tables$counts),
      observed_global_statistic =
        da$global_statistic[[bundle$observed_index]],
      exact_global_p = da$global_p,
      role = if (partition_id == "cluster_phase") {
        "support_gated_sensitivity"
      } else {
        "predefined_inferential"
      },
      stringsAsFactors = FALSE
    )
    da$state_results$partition <- partition_id
    da$state_results$inferential_state <- !(
      partition_id == "phase" &
        da$state_results$state %in%
          as.character(cfg$state_decomposition$phase_descriptive_levels)
    )
    da_state_rows[[partition_id]] <- da$state_results
    support <- as.data.frame(as.table(tables$counts), stringsAsFactors = FALSE)
    names(support) <- c("state", "sample_id", "n_cells")
    support$partition <- partition_id
    support$support_pass <- support$n_cells >= min_cells
    support_rows[[partition_id]] <- support
    distributions[[paste0(partition_id, "__DA")]] <- list(
      statistic = da$global_statistic,
      assignment_ids = bundle$assignments$assignment_id
    )
    gated <- ptpv5_support_gate_state_tables(
      tables, mouse_meta, min_cells, min_mice
    )
    gated$audit$partition <- partition_id
    support_rows[[paste0(partition_id, "__gate")]] <- transform(
      gated$audit,
      sample_id = NA_character_,
      n_cells = NA_integer_,
      support_pass = expression_decomposition_support_pass
    )
    if (nrow(gated$tables$counts) < 2L) {
      ptpv5_message(
        "State decomposition expression branch not estimable for partition ",
        partition_id, ": fewer than two support-gated states.",
        log_file = context$log_file
      )
      next
    }
    for (i in seq_len(nrow(family_meta))) {
      genes <- ptpv5_split_genes(family_meta$union_genes[[i]])
      fit <- ptpv5_decompose_family(
        genes,
        expression,
        gated$tables,
        mouse_meta,
        A,
        bundle$observed_index,
        min_mice
      )
      if (is.null(fit)) next
      z <- fit$observed
      z$partition <- partition_id
      z$response_family <- family_meta$response_family[[i]]
      z$response_family_label <- family_meta$response_family_label[[i]]
      z$tiers <- family_meta$tiers[[i]]
      decomposition_rows[[paste(partition_id, i)]] <- z
      distributions[[paste0(
        partition_id, "__", family_meta$response_family[[i]]
      )]] <- fit$distributions
    }
    ptpv5_message(
      "State decomposition partition completed: ", partition_id,
      log_file = context$log_file
    )
  }
  da_results <- ptpv5_bind_rows(da_rows)
  da_results$fdr_across_partitions <- stats::p.adjust(
    da_results$exact_global_p, "BH"
  )
  da_state <- ptpv5_bind_rows(da_state_rows)
  da_state$fdr_descriptive <- ave(
    da_state$exact_p,
    da_state$partition,
    FUN = function(x) stats::p.adjust(x, "BH")
  )
  decomposition <- ptpv5_bind_rows(decomposition_rows)
  decomposition <- ptpv5_add_fdr(
    decomposition, "adaptive_exact_p", "decomposition_fdr"
  )
  support <- ptpv5_bind_rows(support_rows)
  ptp_write_csv(
    da_results,
    file.path(dirs$state, "state_differential_abundance_global.csv")
  )
  ptp_write_csv(
    da_state,
    file.path(dirs$state, "state_differential_abundance_states.csv")
  )
  ptp_write_csv(
    decomposition,
    file.path(dirs$state, "state_expression_decomposition_results.csv")
  )
  ptp_write_csv(
    support,
    file.path(dirs$state, "mouse_state_support.csv")
  )
  saveRDS(
    distributions,
    file.path(dirs$state, "state_decomposition_exact_distributions.rds"),
    compress = "xz"
  )
  status <- ptpv5_write_method_status(
    dirs,
    "state_decomposition",
    "complete",
    paste0(
      "Differential abundance and expression decomposition completed for ",
      nrow(da_results), " partitions and ", nrow(decomposition),
      " family-component tests."
    ),
    nrow(da_results) + nrow(decomposition)
  )
  out <- list(
    differential_abundance = da_results,
    state_results = da_state,
    decomposition = decomposition,
    support = support,
    distributions = distributions,
    status = status
  )
  ptpv5_save_checkpoint(
    dirs, checkpoint_id, out, cfg, context$config_sha, context$input_sha
  )
  out
}
