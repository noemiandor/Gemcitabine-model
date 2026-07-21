ptpv5_compile_gsea_cpp <- function() {
  if (exists("ptpv5_gsea_matrix_cpp", mode = "function", inherits = TRUE)) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp is required for exact ranked enrichment.", call. = FALSE)
  }
  Rcpp::cppFunction(
    depends = "Rcpp",
    plugins = "cpp11",
    code = '
      Rcpp::List ptpv5_gsea_matrix_cpp(
        Rcpp::NumericMatrix statistic,
        Rcpp::List gene_sets
      ) {
        const int G = statistic.nrow();
        const int B = statistic.ncol();
        const int S = gene_sets.size();
        Rcpp::NumericMatrix weighted(S, B);
        Rcpp::NumericMatrix unweighted(S, B);
        std::vector< std::vector<int> > sets(S);
        for (int s = 0; s < S; ++s) {
          Rcpp::IntegerVector idx = gene_sets[s];
          sets[s].reserve(idx.size());
          for (int k = 0; k < idx.size(); ++k) {
            int z = idx[k] - 1;
            if (z >= 0 && z < G) sets[s].push_back(z);
          }
        }
        std::vector<int> order(G), rank(G);
        for (int b = 0; b < B; ++b) {
          for (int g = 0; g < G; ++g) order[g] = g;
          std::sort(
            order.begin(),
            order.end(),
            [&](int a, int c) {
              double va = statistic(a, b);
              double vc = statistic(c, b);
              if (!R_finite(va)) va = -INFINITY;
              if (!R_finite(vc)) vc = -INFINITY;
              if (va == vc) return a < c;
              return va > vc;
            }
          );
          for (int r = 0; r < G; ++r) rank[order[r]] = r;
          for (int s = 0; s < S; ++s) {
            const int K = sets[s].size();
            if (K < 1 || K >= G) {
              weighted(s, b) = NA_REAL;
              unweighted(s, b) = NA_REAL;
              continue;
            }
            std::vector< std::pair<int, double> > hits;
            hits.reserve(K);
            double weight_sum = 0.0;
            for (int k = 0; k < K; ++k) {
              const int g = sets[s][k];
              double w = std::fabs(statistic(g, b));
              if (!R_finite(w)) w = 0.0;
              hits.push_back(std::make_pair(rank[g], w));
              weight_sum += w;
            }
            std::sort(hits.begin(), hits.end());
            if (!(weight_sum > 0.0) || !R_finite(weight_sum)) {
              weight_sum = static_cast<double>(K);
              for (int k = 0; k < K; ++k) hits[k].second = 1.0;
            }
            const double miss = 1.0 / static_cast<double>(G - K);
            double run_w = 0.0, run_u = 0.0;
            double max_w = -INFINITY, min_w = INFINITY;
            double max_u = -INFINITY, min_u = INFINITY;
            int previous_rank = -1;
            for (int k = 0; k < K; ++k) {
              const int gap = hits[k].first - previous_rank - 1;
              run_w -= gap * miss;
              run_u -= gap * miss;
              if (run_w > max_w) max_w = run_w;
              if (run_w < min_w) min_w = run_w;
              if (run_u > max_u) max_u = run_u;
              if (run_u < min_u) min_u = run_u;
              run_w += hits[k].second / weight_sum;
              run_u += 1.0 / static_cast<double>(K);
              if (run_w > max_w) max_w = run_w;
              if (run_w < min_w) min_w = run_w;
              if (run_u > max_u) max_u = run_u;
              if (run_u < min_u) min_u = run_u;
              previous_rank = hits[k].first;
            }
            const int tail = G - previous_rank - 1;
            run_w -= tail * miss;
            run_u -= tail * miss;
            if (run_w > max_w) max_w = run_w;
            if (run_w < min_w) min_w = run_w;
            if (run_u > max_u) max_u = run_u;
            if (run_u < min_u) min_u = run_u;
            weighted(s, b) =
              std::fabs(max_w) >= std::fabs(min_w) ? max_w : min_w;
            unweighted(s, b) =
              std::fabs(max_u) >= std::fabs(min_u) ? max_u : min_u;
          }
        }
        return Rcpp::List::create(
          Rcpp::Named("weighted") = weighted,
          Rcpp::Named("unweighted") = unweighted
        );
      }
    ',
    env = .GlobalEnv
  )
  invisible(TRUE)
}

ptpv5_gsea_leading_edge <- function(statistic, set_indices, gene_symbols) {
  statistic <- as.numeric(statistic)
  ord <- order(statistic, decreasing = TRUE, na.last = TRUE)
  rank <- integer(length(statistic))
  rank[ord] <- seq_along(ord)
  set_indices <- set_indices[is.finite(set_indices)]
  set_indices <- set_indices[set_indices >= 1L & set_indices <= length(statistic)]
  hit <- sort(rank[set_indices])
  K <- length(hit)
  G <- length(statistic)
  if (!K || K >= G) {
    return(list(direction = NA_character_, genes = character(), rank = NA_integer_))
  }
  weights <- abs(statistic[ord[hit]])
  if (!sum(weights, na.rm = TRUE) > 0) weights[] <- 1
  weights <- weights / sum(weights)
  miss <- 1 / (G - K)
  running <- numeric(G)
  is_hit <- logical(G)
  is_hit[hit] <- TRUE
  w_by_rank <- numeric(G)
  w_by_rank[hit] <- weights
  running <- cumsum(ifelse(is_hit, w_by_rank, -miss))
  peak_pos <- which.max(abs(running))
  direction <- if (running[[peak_pos]] >= 0) "positive" else "negative"
  edge_rank <- if (direction == "positive") {
    hit[hit <= peak_pos]
  } else {
    hit[hit >= peak_pos]
  }
  edge_idx <- ord[edge_rank]
  list(
    direction = direction,
    genes = gene_symbols[edge_idx],
    rank = peak_pos
  )
}

ptpv5_run_exact_ranked_enrichment <- function(context) {
  cfg <- context$cfg
  args <- context$args
  dirs <- context$dirs
  bundle <- context$bundle
  checkpoint_id <- "method01_exact_ranked_enrichment"
  cached <- ptpv5_load_checkpoint(
    dirs, checkpoint_id, cfg, context$config_sha, context$input_sha
  )
  if (!is.null(cached)) return(cached)
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("rhdf5 is required for V5 method 1.", call. = FALSE)
  }
  ptpv5_compile_gsea_cpp()
  sets <- ptpv5_frozen_sets(bundle, args)
  all_sets <- c(sets$membership_sets, sets$family_sets)
  set_type <- c(
    rep("unique_membership", length(sets$membership_sets)),
    rep("response_family_union", length(sets$family_sets))
  )
  set_id <- names(all_sets)
  minimum <- as.integer(cfg$ranked_enrichment$minimum_set_genes)
  keep <- lengths(all_sets) >= minimum
  all_sets <- all_sets[keep]
  set_type <- set_type[keep]
  set_id <- set_id[keep]
  B <- nrow(bundle$assignments)
  weighted <- matrix(
    NA_real_, nrow = length(all_sets), ncol = B,
    dimnames = list(set_id, bundle$assignments$assignment_id)
  )
  unweighted <- weighted
  chunk_size <- as.integer(cfg$ranked_enrichment$assignment_chunk_size)
  h5_columns <- bundle$assignment_indices_hdf5
  chunks <- split(seq_len(B), ceiling(seq_len(B) / chunk_size))
  for (ii in seq_along(chunks)) {
    local_cols <- chunks[[ii]]
    h5_cols <- h5_columns[local_cols]
    stat <- rhdf5::h5read(
      cfg$inputs$v4_hdf5,
      "t_statistic",
      index = list(seq_along(bundle$gene_ids), h5_cols)
    )
    stat <- as.matrix(stat)
    if (length(local_cols) == 1L) stat <- matrix(stat, ncol = 1L)
    score <- ptpv5_gsea_matrix_cpp(
      stat,
      unname(lapply(all_sets, as.integer))
    )
    weighted[, local_cols] <- score$weighted
    unweighted[, local_cols] <- score$unweighted
    ptpv5_message(
      "Ranked enrichment assignment chunk ", ii, "/", length(chunks),
      log_file = context$log_file
    )
  }
  obs <- bundle$observed_index
  abs_weighted <- abs(weighted)
  abs_unweighted <- abs(unweighted)
  maxT_weighted <- apply(abs_weighted, 2L, max, na.rm = TRUE)
  maxT_unweighted <- apply(abs_unweighted, 2L, max, na.rm = TRUE)
  rows <- lapply(seq_along(all_sets), function(i) {
    meta <- if (set_type[[i]] == "unique_membership") {
      x <- sets$memberships[
        match(set_id[[i]], sets$memberships$membership_id),
        ,
        drop = FALSE
      ]
      data.frame(
        tiers = x$tiers,
        response_family = x$response_families,
        aliases = x$aliases,
        stringsAsFactors = FALSE
      )
    } else {
      x <- sets$families[
        match(set_id[[i]], sets$families$response_family),
        ,
        drop = FALSE
      ]
      data.frame(
        tiers = x$tiers,
        response_family = x$response_family,
        aliases = x$response_family_label,
        stringsAsFactors = FALSE
      )
    }
    data.frame(
      set_type = set_type[[i]],
      set_id = set_id[[i]],
      n_genes = length(all_sets[[i]]),
      weighted_es = weighted[i, obs],
      weighted_exact_p = ptpv5_finite_p(
        abs_weighted[i, obs], abs_weighted[i, ]
      ),
      weighted_westfall_young_maxT_p = ptpv5_finite_p(
        abs_weighted[i, obs], maxT_weighted
      ),
      unweighted_es = unweighted[i, obs],
      unweighted_exact_p = ptpv5_finite_p(
        abs_unweighted[i, obs], abs_unweighted[i, ]
      ),
      unweighted_westfall_young_maxT_p = ptpv5_finite_p(
        abs_unweighted[i, obs], maxT_unweighted
      ),
      meta,
      stringsAsFactors = FALSE
    )
  })
  results <- ptpv5_bind_rows(rows)
  results <- ptpv5_add_fdr(results, "weighted_exact_p", "weighted_fdr")
  results <- ptpv5_add_fdr(results, "unweighted_exact_p", "unweighted_fdr")
  obs_stat <- rhdf5::h5read(
    cfg$inputs$v4_hdf5,
    "t_statistic",
    index = list(
      seq_along(bundle$gene_ids),
      bundle$assignment_indices_hdf5[[obs]]
    )
  )
  leading <- lapply(seq_along(all_sets), function(i) {
    le <- ptpv5_gsea_leading_edge(
      obs_stat, all_sets[[i]], bundle$gene_symbols
    )
    data.frame(
      set_type = set_type[[i]],
      set_id = set_id[[i]],
      direction = le$direction,
      peak_rank = le$rank,
      n_leading_edge_genes = length(le$genes),
      leading_edge_genes = paste(le$genes, collapse = ";"),
      inferential_role = "descriptive_after_set_level_inference_only",
      stringsAsFactors = FALSE
    )
  })
  leading <- ptpv5_bind_rows(leading)
  distribution <- list(
    weighted_es = weighted,
    unweighted_es = unweighted,
    assignment_ids = bundle$assignments$assignment_id,
    observed_index = obs
  )
  ptp_write_csv(
    results,
    file.path(dirs$ranked, "exact_ranked_enrichment_results.csv")
  )
  ptp_write_csv(
    leading,
    file.path(dirs$ranked, "observed_leading_edges_descriptive.csv")
  )
  saveRDS(
    distribution,
    file.path(dirs$ranked, "exact_ranked_enrichment_distributions.rds"),
    compress = "xz"
  )
  status <- ptpv5_write_method_status(
    dirs,
    "exact_ranked_enrichment",
    "complete",
    paste0(
      "Weighted and unweighted whole-mouse ranked enrichment completed for ",
      nrow(results), " frozen gene sets."
    ),
    nrow(results)
  )
  out <- list(
    results = results,
    leading_edge = leading,
    distribution = distribution,
    status = status
  )
  ptpv5_save_checkpoint(
    dirs, checkpoint_id, out, cfg, context$config_sha, context$input_sha
  )
  out
}
