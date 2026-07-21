ptpv5_assign_merged_bin <- function(pseudotime, bins) {
  out <- rep(NA_character_, length(pseudotime))
  for (i in seq_len(nrow(bins))) {
    include <- pseudotime >= bins$start[[i]] &
      if (isTRUE(bins$include_end[[i]])) {
        pseudotime <= bins$end[[i]]
      } else {
        pseudotime < bins$end[[i]]
      }
    out[include] <- bins$region_id[[i]]
  }
  out
}

ptpv5_trajectory_pseudobulk <- function(state_data, bins, sample_ids, minimum_cells) {
  cells <- state_data$cells
  cells$region_id <- ptpv5_assign_merged_bin(cells$pseudotime, bins)
  keys <- expand.grid(
    sample_id = sample_ids,
    region_id = bins$region_id,
    stringsAsFactors = FALSE
  )
  columns <- list()
  meta_rows <- list()
  for (i in seq_len(nrow(keys))) {
    idx <- which(
      cells$sample_id == keys$sample_id[[i]] &
        cells$region_id == keys$region_id[[i]]
    )
    if (length(idx) < minimum_cells) next
    columns[[length(columns) + 1L]] <- Matrix::rowSums(
      state_data$counts[, idx, drop = FALSE]
    )
    meta_rows[[length(meta_rows) + 1L]] <- data.frame(
      sample_id = keys$sample_id[[i]],
      region_id = keys$region_id[[i]],
      n_cells = length(idx),
      stringsAsFactors = FALSE
    )
  }
  if (!length(columns)) {
    stop("No trajectory pseudobulk columns passed the cell threshold.", call. = FALSE)
  }
  counts <- do.call(cbind, columns)
  rownames(counts) <- rownames(state_data$counts)
  meta <- ptpv5_bind_rows(meta_rows)
  colnames(counts) <- paste(meta$sample_id, meta$region_id, sep = "__")
  dge <- edgeR::DGEList(counts = counts)
  keep <- rowSums(counts) > 0
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 0.5)
  logcpm[!keep, ] <- NA_real_
  list(counts = counts, logcpm = logcpm, metadata = meta)
}

ptpv5_trajectory_interactions <- function(
  logcpm,
  metadata,
  bins,
  mouse_meta,
  A,
  min_mice
) {
  output <- list()
  for (region_id in bins$region_id) {
    Y <- matrix(
      NA_real_,
      nrow = nrow(logcpm),
      ncol = nrow(mouse_meta),
      dimnames = list(rownames(logcpm), mouse_meta$sample_id)
    )
    local <- which(metadata$region_id == region_id)
    if (length(local)) {
      sample_idx <- match(metadata$sample_id[local], mouse_meta$sample_id)
      Y[, sample_idx] <- logcpm[, local, drop = FALSE]
    }
    output[[region_id]] <- ptpv5_mouse_interaction_matrix(
      Y, mouse_meta, A, min_group = min_mice
    )$interaction
  }
  output
}

ptpv5_compile_higher_criticism_cpp <- function() {
  if (exists("ptpv5_higher_criticism_cpp", mode = "function", inherits = TRUE)) {
    return(invisible(TRUE))
  }
  Rcpp::cppFunction(
    plugins = "cpp11",
    code = '
      Rcpp::NumericVector ptpv5_higher_criticism_cpp(
        Rcpp::NumericMatrix pvalue
      ) {
        const int G = pvalue.nrow();
        const int B = pvalue.ncol();
        Rcpp::NumericVector out(B, NA_REAL);
        for (int b = 0; b < B; ++b) {
          std::vector<double> p;
          p.reserve(G);
          for (int g = 0; g < G; ++g) {
            double value = pvalue(g, b);
            if (R_finite(value) && value > 0.0 && value < 1.0) {
              p.push_back(value);
            }
          }
          const int m = p.size();
          if (m < 2) continue;
          std::sort(p.begin(), p.end());
          double best = 0.0;
          const int upper = std::max(1, static_cast<int>(std::floor(0.5 * m)));
          for (int i = 0; i < upper; ++i) {
            const double pv = p[i];
            const double empirical = static_cast<double>(i + 1) / m;
            const double denom = std::sqrt(pv * (1.0 - pv));
            if (denom > 0.0) {
              const double hc = std::sqrt(static_cast<double>(m)) *
                (empirical - pv) / denom;
              if (R_finite(hc) && hc > best) best = hc;
            }
          }
          out[b] = best;
        }
        return out;
      }
    ',
    env = .GlobalEnv
  )
  invisible(TRUE)
}

ptpv5_gene_spline_statistics <- function(interactions, bins, df_candidates) {
  M <- length(interactions)
  G <- nrow(interactions[[1L]])
  B <- ncol(interactions[[1L]])
  Z <- array(NA_real_, dim = c(G, B, M))
  raw_mean <- matrix(0, nrow = G, ncol = B)
  for (m in seq_len(M)) {
    Z[, , m] <- ptpv5_standardize_rows_by_randomization(
      interactions[[m]]
    )
    raw_mean <- raw_mean + interactions[[m]] / M
  }
  q_list <- list()
  basis_rows <- list()
  x <- bins$midpoint
  for (df in df_candidates) {
    if (df > M) next
    basis <- splines::ns(x, df = df, intercept = TRUE)
    P <- basis %*% solve(crossprod(basis)) %*% t(basis)
    Q <- matrix(0, nrow = G, ncol = B)
    for (a in seq_len(M)) {
      for (b in seq_len(M)) {
        Q <- Q + P[a, b] * Z[, , a] * Z[, , b]
      }
    }
    q_list[[paste0("df", df)]] <- Q / df
    basis_rows[[paste0("df", df)]] <- data.frame(
      spline_df = df,
      region_id = bins$region_id,
      midpoint = bins$midpoint,
      basis,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (!length(q_list)) {
    stop("No candidate spline degree is estimable.", call. = FALSE)
  }
  selected <- q_list[[1L]]
  selected_df <- matrix(
    names(q_list)[[1L]],
    nrow = G,
    ncol = B
  )
  if (length(q_list) > 1L) {
    for (i in 2:length(q_list)) {
      replace <- q_list[[i]] > selected
      selected[replace] <- q_list[[i]][replace]
      selected_df[replace] <- names(q_list)[[i]]
    }
  }
  signed <- sign(raw_mean) * sqrt(pmax(selected, 0))
  list(
    statistic = selected,
    signed_statistic = signed,
    selected_df = selected_df,
    candidate_statistics = q_list,
    basis = ptpv5_bind_rows(basis_rows)
  )
}

ptpv5_run_gene_wise_trajectory <- function(context) {
  cfg <- context$cfg
  args <- context$args
  dirs <- context$dirs
  bundle <- context$bundle
  checkpoint_id <- "method05_gene_wise_trajectory"
  cached <- ptpv5_load_checkpoint(
    dirs, checkpoint_id, cfg, context$config_sha, context$input_sha
  )
  if (!is.null(cached)) return(cached)
  sets <- ptpv5_frozen_sets(bundle, args)
  family_meta <- sets$families
  family_genes <- sort(unique(unlist(lapply(
    family_meta$union_genes,
    ptpv5_split_genes
  ))))
  if (isTRUE(args$smoke)) {
    family_genes <- head(family_genes, as.integer(args$smoke_genes))
  }
  state_data <- ptpv5_load_state_expression(context, family_genes)
  trajectory_v4 <- readRDS(cfg$inputs$v4_trajectory_checkpoint)$data
  bins <- trajectory_v4$merge$bins
  support <- trajectory_v4$support
  support$support_pass <- tolower(as.character(support$support_pass)) %in%
    c("true", "t", "1")
  supported <- vapply(bins$region_id, function(id) {
    x <- support[support$region_id == id, , drop = FALSE]
    nrow(x) == 4L && all(x$support_pass) &&
      all(x$n_mice_qc >= as.integer(
        cfg$gene_wise_trajectory$support_min_mice_per_group
      ))
  }, logical(1L))
  supported_bins <- bins[supported, , drop = FALSE]
  minimum_supported <- as.integer(
    cfg$gene_wise_trajectory$minimum_supported_bins
  )
  primary_interval <- as.numeric(
    cfg$gene_wise_trajectory$primary_interval
  )
  primary_bins <- supported_bins[
    supported_bins$midpoint >= primary_interval[[1L]] &
      supported_bins$midpoint <= primary_interval[[2L]],
    ,
    drop = FALSE
  ]
  estimability <- data.frame(
    branch = c("primary_0.30_0.49", "whole_trajectory_0.00_1.00"),
    n_supported_bins = c(nrow(primary_bins), nrow(supported_bins)),
    minimum_required_bins = minimum_supported,
    status = c(
      if (nrow(primary_bins) >= minimum_supported) "estimable" else "not_estimable",
      if (nrow(supported_bins) >= minimum_supported) "estimable" else "not_estimable"
    ),
    reason = c(
      if (nrow(primary_bins) >= minimum_supported) {
        "minimum support satisfied"
      } else {
        "fewer than three four-group-supported bins; no extrapolation"
      },
      if (nrow(supported_bins) >= minimum_supported) {
        "minimum support satisfied"
      } else {
        "fewer than three four-group-supported bins"
      }
    ),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(
    estimability,
    file.path(dirs$trajectory, "trajectory_estimability.csv")
  )
  if (nrow(supported_bins) < minimum_supported) {
    status <- ptpv5_write_method_status(
      dirs,
      "gene_wise_trajectory",
      "not_estimable",
      estimability$reason[[2L]],
      0L
    )
    out <- list(
      estimability = estimability,
      gene_results = data.frame(),
      family_results = data.frame(),
      status = status
    )
    ptpv5_save_checkpoint(
      dirs, checkpoint_id, out, cfg, context$config_sha, context$input_sha
    )
    return(out)
  }
  minimum_cells <- min(
    support$minimum_cells[is.finite(support$minimum_cells)],
    na.rm = TRUE
  )
  if (!is.finite(minimum_cells)) minimum_cells <- 5L
  pb <- ptpv5_trajectory_pseudobulk(
    state_data,
    supported_bins,
    bundle$sample_ids,
    as.integer(minimum_cells)
  )
  rm(state_data)
  invisible(gc())
  mouse_meta <- bundle$mouse$metadata[
    match(bundle$sample_ids, bundle$mouse$metadata$sample_id),
    ,
    drop = FALSE
  ]
  A <- ptpv5_assignment_matrix(bundle$assignments, bundle$sample_ids)
  interactions <- ptpv5_trajectory_interactions(
    pb$logcpm,
    pb$metadata,
    supported_bins,
    mouse_meta,
    A,
    as.integer(cfg$gene_wise_trajectory$support_min_mice_per_group)
  )
  fit <- ptpv5_gene_spline_statistics(
    interactions,
    supported_bins,
    as.integer(cfg$gene_wise_trajectory$spline_df_candidates)
  )
  obs <- bundle$observed_index
  gene_stat <- fit$statistic
  signed_stat <- fit$signed_statistic
  gene_rank <- matrixStats::rowRanks(
    -gene_stat,
    ties.method = "max",
    na.last = "keep",
    preserveShape = TRUE
  )
  gene_finite_n <- rowSums(is.finite(gene_stat))
  gene_assignment_p <- sweep(
    1 + gene_rank,
    1L,
    1 + gene_finite_n,
    "/"
  )
  gene_results <- data.frame(
    gene_symbol = rownames(pb$logcpm),
    n_supported_bins = nrow(supported_bins),
    observed_global_statistic = gene_stat[, obs],
    observed_signed_statistic = signed_stat[, obs],
    observed_selected_spline_df = fit$selected_df[, obs],
    exact_whole_mouse_permutation_p = gene_assignment_p[, obs],
    stringsAsFactors = FALSE
  )
  gene_results$fdr_frozen_family_gene_universe <- stats::p.adjust(
    gene_results$exact_whole_mouse_permutation_p, "BH"
  )
  ptpv5_compile_higher_criticism_cpp()
  family_rows <- list()
  family_distributions <- list()
  signed_Z <- ptpv5_standardize_rows_by_randomization(signed_stat)
  for (i in seq_len(nrow(family_meta))) {
    genes <- ptpv5_split_genes(family_meta$union_genes[[i]])
    idx <- match(genes, rownames(pb$logcpm))
    idx <- idx[is.finite(idx)]
    if (length(idx) < 2L) next
    Z <- signed_Z[idx, , drop = FALSE]
    stat <- rbind(
      true_maxmean = ptpv5_vector_stat_matrix(Z)["true_maxmean", ],
      mean_square = ptpv5_vector_stat_matrix(Z)["mean_square", ],
      higher_criticism = ptpv5_higher_criticism_cpp(
        gene_assignment_p[idx, , drop = FALSE]
      )
    )
    nested <- ptpv5_nested_component_test(t(stat), obs)
    family_rows[[family_meta$response_family[[i]]]] <- data.frame(
      response_family = family_meta$response_family[[i]],
      response_family_label = family_meta$response_family_label[[i]],
      tiers = family_meta$tiers[[i]],
      n_genes = length(idx),
      observed_true_maxmean = stat["true_maxmean", obs],
      observed_mean_square = stat["mean_square", obs],
      observed_higher_criticism = stat["higher_criticism", obs],
      adaptive_exact_p = nested$observed_adaptive_p,
      stringsAsFactors = FALSE
    )
    family_distributions[[family_meta$response_family[[i]]]] <- list(
      statistic = t(stat),
      adaptive_assignment_p = nested$adaptive_assignment_p
    )
  }
  family_results <- ptpv5_bind_rows(family_rows)
  family_results <- ptpv5_add_fdr(
    family_results, "adaptive_exact_p", "trajectory_fdr"
  )
  ptp_write_csv(
    pb$metadata,
    file.path(dirs$trajectory, "trajectory_pseudobulk_support.csv")
  )
  ptp_write_csv(
    fit$basis,
    file.path(dirs$trajectory, "gene_wise_spline_basis_audit.csv")
  )
  ptp_write_csv(
    gene_results,
    file.path(dirs$trajectory, "gene_wise_trajectory_results.csv")
  )
  ptp_write_csv(
    family_results,
    file.path(dirs$trajectory, "family_gene_wise_trajectory_results.csv")
  )
  saveRDS(
    list(
      family = family_distributions,
      assignment_ids = bundle$assignments$assignment_id,
      observed_index = obs
    ),
    file.path(dirs$trajectory, "trajectory_family_exact_distributions.rds"),
    compress = "xz"
  )
  status <- ptpv5_write_method_status(
    dirs,
    "gene_wise_trajectory",
    "complete_with_primary_not_estimable",
    paste0(
      "Whole-trajectory gene-wise spline tests completed for ",
      nrow(gene_results), " frozen family genes; the primary interval branch ",
      "is explicitly not estimable because it has ", nrow(primary_bins),
      " supported bin(s)."
    ),
    nrow(gene_results) + nrow(family_results)
  )
  out <- list(
    estimability = estimability,
    gene_results = gene_results,
    family_results = family_results,
    family_distributions = family_distributions,
    support = pb$metadata,
    basis = fit$basis,
    status = status
  )
  ptpv5_save_checkpoint(
    dirs, checkpoint_id, out, cfg, context$config_sha, context$input_sha
  )
  out
}
