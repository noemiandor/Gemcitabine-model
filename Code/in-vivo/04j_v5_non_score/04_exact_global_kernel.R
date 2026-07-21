ptpv5_assignment_contrasts <- function(mouse_meta, A) {
  p2 <- mouse_meta$initial_ploidy == "2N"
  p4 <- mouse_meta$initial_ploidy == "4N"
  C <- matrix(0, nrow = nrow(A), ncol = ncol(A), dimnames = dimnames(A))
  C[p4, ] <- ifelse(A[p4, , drop = FALSE], 1 / 4, -1 / 4)
  C[p2, ] <- ifelse(A[p2, , drop = FALSE], -1 / 4, 1 / 4)
  C
}

ptpv5_kernel_distribution <- function(Y, mouse_meta, C, ridge_lambda) {
  if (nrow(Y) < 1L) return(rep(NA_real_, ncol(C)))
  Y <- as.matrix(Y)
  center <- rowMeans(Y, na.rm = TRUE)
  scale <- matrixStats::rowSds(Y, na.rm = TRUE)
  keep <- is.finite(scale) & scale > 0
  Y <- Y[keep, , drop = FALSE]
  if (!nrow(Y)) return(rep(NA_real_, ncol(C)))
  Y <- sweep(Y, 1L, center[keep], "-")
  Y <- sweep(Y, 1L, scale[keep], "/")
  X0 <- cbind(
    intercept = 1,
    ploidy_4N = as.numeric(mouse_meta$initial_ploidy == "4N")
  )
  H0 <- X0 %*% solve(crossprod(X0)) %*% t(X0)
  Y <- Y %*% (diag(ncol(Y)) - H0)
  K <- crossprod(Y) / nrow(Y)
  eigen_fit <- eigen(K, symmetric = TRUE)
  positive <- eigen_fit$values[eigen_fit$values > 1e-10]
  ridge <- as.numeric(ridge_lambda) *
    stats::median(positive, na.rm = TRUE)
  if (!is.finite(ridge) || ridge <= 0) ridge <- as.numeric(ridge_lambda)
  adjusted <- pmax(eigen_fit$values, 0)^2 /
    (pmax(eigen_fit$values, 0) + ridge)
  K_ridge <- eigen_fit$vectors %*%
    (adjusted * t(eigen_fit$vectors))
  colSums(C * (K_ridge %*% C))
}

ptpv5_kernel_set_test <- function(
  set_indices,
  Y,
  mouse_meta,
  C,
  ridge_lambda,
  observed_index
) {
  set_indices <- set_indices[
    is.finite(set_indices) &
      set_indices >= 1L &
      set_indices <= nrow(Y)
  ]
  stat <- ptpv5_kernel_distribution(
    Y[set_indices, , drop = FALSE],
    mouse_meta,
    C,
    ridge_lambda
  )
  list(
    statistic = stat,
    observed = stat[[observed_index]],
    exact_p = ptpv5_finite_p(stat[[observed_index]], stat)
  )
}

ptpv5_run_exact_global_kernel <- function(context) {
  cfg <- context$cfg
  args <- context$args
  dirs <- context$dirs
  bundle <- context$bundle
  checkpoint_id <- "method04_exact_global_kernel"
  cached <- ptpv5_load_checkpoint(
    dirs, checkpoint_id, cfg, context$config_sha, context$input_sha
  )
  if (!is.null(cached)) return(cached)
  target_families <- as.character(
    cfg$mechanistic_kernel$response_families
  )
  sets <- ptpv5_frozen_sets(bundle, args)
  family_meta <- bundle$families[
    match(target_families, bundle$families$response_family),
    ,
    drop = FALSE
  ]
  family_sets <- lapply(family_meta$union_genes, function(x) {
    ptpv5_set_indices(ptpv5_split_genes(x), bundle)
  })
  names(family_sets) <- family_meta$response_family
  Y <- bundle$mouse$expression
  mouse_meta <- bundle$mouse$metadata[
    match(colnames(Y), bundle$mouse$metadata$sample_id),
    ,
    drop = FALSE
  ]
  A <- ptpv5_assignment_matrix(bundle$assignments, bundle$sample_ids)
  C <- ptpv5_assignment_contrasts(mouse_meta, A)
  ridge <- as.numeric(cfg$mechanistic_kernel$ridge_lambda)
  obs <- bundle$observed_index
  family_rows <- list()
  distributions <- list()
  fragility_rows <- list()
  for (i in seq_along(family_sets)) {
    family_id <- names(family_sets)[[i]]
    idx <- family_sets[[i]]
    fit <- ptpv5_kernel_set_test(
      idx, Y, mouse_meta, C, ridge, obs
    )
    family_rows[[family_id]] <- data.frame(
      response_family = family_id,
      response_family_label =
        family_meta$response_family_label[[i]],
      tiers = family_meta$tiers[[i]],
      n_genes = length(idx),
      observed_kernel_statistic = fit$observed,
      exact_kernel_p = fit$exact_p,
      stringsAsFactors = FALSE
    )
    distributions[[family_id]] <- fit$statistic
    gene_contribution <- (Y[idx, , drop = FALSE] %*% C[, obs])^2
    drop_local <- which.max(gene_contribution)
    drop_idx <- idx[[drop_local]]
    reduced <- ptpv5_kernel_set_test(
      setdiff(idx, drop_idx), Y, mouse_meta, C, ridge, obs
    )
    fragility_rows[[family_id]] <- data.frame(
      response_family = family_id,
      dropped_gene = bundle$gene_symbols[[drop_idx]],
      dropped_gene_contribution = gene_contribution[[drop_local]],
      full_exact_p = fit$exact_p,
      leave_top_gene_out_exact_p = reduced$exact_p,
      full_observed_statistic = fit$observed,
      leave_top_gene_out_observed_statistic = reduced$observed,
      stringsAsFactors = FALSE
    )
  }
  family_results <- ptpv5_bind_rows(family_rows)
  family_results$holm_five_families <- stats::p.adjust(
    family_results$exact_kernel_p, "holm"
  )
  family_results$family_gate_pass <- (
    family_results$holm_five_families <=
      as.numeric(cfg$mechanistic_kernel$family_threshold)
  )

  alias <- bundle$alias
  membership_meta <- bundle$memberships
  membership_sets <- lapply(membership_meta$genes, function(x) {
    ptpv5_set_indices(ptpv5_split_genes(x), bundle)
  })
  names(membership_sets) <- membership_meta$membership_id
  child_rows <- list()
  family_minp_rows <- list()
  child_distributions <- list()
  for (family_id in target_families) {
    mids <- unique(alias$membership_id[alias$response_family == family_id])
    mids <- intersect(mids, names(membership_sets))
    if (!length(mids)) next
    child_stat <- matrix(
      NA_real_,
      nrow = nrow(bundle$assignments),
      ncol = length(mids),
      dimnames = list(bundle$assignments$assignment_id, mids)
    )
    for (j in seq_along(mids)) {
      fit <- ptpv5_kernel_set_test(
        membership_sets[[mids[[j]]]],
        Y,
        mouse_meta,
        C,
        ridge,
        obs
      )
      child_stat[, j] <- fit$statistic
    }
    child_p_assignment <- ptpv5_empirical_rank_p_matrix(
      child_stat, larger = TRUE
    )
    observed_child_p <- child_p_assignment[obs, ]
    minp <- apply(child_p_assignment, 1L, min, na.rm = TRUE)
    observed_minp <- minp[[obs]]
    family_minp_p <- (
      1 + sum(minp <= observed_minp, na.rm = TRUE)
    ) / (1 + sum(is.finite(minp)))
    min_child_p_by_assignment <- apply(
      child_p_assignment, 1L, min, na.rm = TRUE
    )
    for (j in seq_along(mids)) {
      wy <- (
        1 + sum(
          min_child_p_by_assignment <= observed_child_p[[j]],
          na.rm = TRUE
        )
      ) / (1 + sum(is.finite(min_child_p_by_assignment)))
      m <- membership_meta[
        match(mids[[j]], membership_meta$membership_id),
        ,
        drop = FALSE
      ]
      child_rows[[paste(family_id, mids[[j]])]] <- data.frame(
        response_family = family_id,
        membership_id = mids[[j]],
        aliases = m$aliases,
        n_genes = m$n_genes,
        observed_kernel_statistic = child_stat[obs, j],
        raw_exact_p = observed_child_p[[j]],
        within_family_westfall_young_p = wy,
        stringsAsFactors = FALSE
      )
    }
    family_minp_rows[[family_id]] <- data.frame(
      response_family = family_id,
      observed_child_min_p = observed_minp,
      nested_exact_child_minP_p = family_minp_p,
      n_child_memberships = length(mids),
      stringsAsFactors = FALSE
    )
    child_distributions[[family_id]] <- list(
      statistic = child_stat,
      assignment_p = child_p_assignment,
      min_p = minp
    )
  }
  children <- ptpv5_bind_rows(child_rows)
  family_minp <- ptpv5_bind_rows(family_minp_rows)
  children <- merge(
    children,
    family_results[, c(
      "response_family", "family_gate_pass", "holm_five_families"
    )],
    by = "response_family",
    all.x = TRUE,
    sort = FALSE
  )
  children$child_gate_pass <- children$family_gate_pass &
    children$within_family_westfall_young_p <=
      as.numeric(cfg$mechanistic_kernel$child_threshold)
  family_results <- merge(
    family_results,
    family_minp,
    by = "response_family",
    all.x = TRUE,
    sort = FALSE
  )

  genes <- bundle$gene_checkpoint$observed
  genes$gene_symbol <- ptp_clean_gene_symbols(genes$gene_symbol)
  gene_rows <- list()
  for (family_id in target_families) {
    family_genes <- ptpv5_split_genes(
      family_meta$union_genes[
        match(family_id, family_meta$response_family)
      ]
    )
    x <- genes[genes$gene_symbol %in% family_genes, , drop = FALSE]
    x$holm_within_family <- stats::p.adjust(x$exact_empirical_p, "holm")
    family_gate <- family_results$family_gate_pass[
      match(family_id, family_results$response_family)
    ]
    child_gate <- any(
      children$response_family == family_id &
        children$child_gate_pass,
      na.rm = TRUE
    )
    x$response_family <- family_id
    x$family_gate_pass <- family_gate
    x$any_child_gate_pass <- child_gate
    x$closed_testing_gene_claim <- family_gate & child_gate &
      x$holm_within_family <=
        as.numeric(cfg$mechanistic_kernel$gene_threshold)
    gene_rows[[family_id]] <- x
  }
  gene_results <- ptpv5_bind_rows(gene_rows)

  directional_id <- as.character(
    cfg$mechanistic_kernel$directional_family
  )
  directional <- data.frame()
  if (directional_id %in% target_families) {
    program_id <- bundle$memberships$directional_program_id[
      grepl(
        directional_id,
        bundle$memberships$response_families,
        fixed = TRUE
      ) &
        bundle$memberships$directional_eligible %in% TRUE
    ][1L]
    w <- bundle$weights[
      bundle$weights$program_id == program_id &
        tolower(as.character(bundle$weights$observed)) %in%
          c("true", "t", "1"),
      ,
      drop = FALSE
    ]
    w$gene <- ptp_clean_gene_symbols(w$gene)
    idx <- unname(bundle$symbol_to_row[w$gene])
    keep <- is.finite(idx) & is.finite(w$signed_weight)
    idx <- idx[keep]
    ww <- w$signed_weight[keep]
    if (length(idx)) {
      if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("rhdf5 is required for the directional kernel child.", call. = FALSE)
      }
      stat <- rhdf5::h5read(
        cfg$inputs$v4_hdf5,
        "t_statistic",
        index = list(idx, bundle$assignment_indices_hdf5)
      )
      stat <- as.matrix(stat)
      direction_stat <- as.numeric(crossprod(ww / sum(abs(ww)), stat))
      directional <- data.frame(
        response_family = directional_id,
        program_id = program_id,
        n_genes = length(idx),
        observed_signed_gene_vector_contrast = direction_stat[[obs]],
        exact_two_sided_p = ptpv5_finite_p(
          abs(direction_stat[[obs]]), abs(direction_stat)
        ),
        inferential_unit =
          "predefined_signed_gene_vector_contrast_not_mouse_program_score",
        stringsAsFactors = FALSE
      )
      distributions[[paste0(directional_id, "__directional")]] <-
        direction_stat
    }
  }

  fragility <- ptpv5_bind_rows(fragility_rows)
  ptp_write_csv(
    family_results,
    file.path(dirs$kernel, "mechanistic_family_kernel_results.csv")
  )
  ptp_write_csv(
    children,
    file.path(dirs$kernel, "mechanistic_membership_kernel_results.csv")
  )
  ptp_write_csv(
    gene_results,
    file.path(dirs$kernel, "mechanistic_gene_closed_testing.csv")
  )
  ptp_write_csv(
    directional,
    file.path(dirs$kernel, "gemcitabine_directional_gene_vector_test.csv")
  )
  ptp_write_csv(
    fragility,
    file.path(dirs$kernel, "single_gene_fragility.csv")
  )
  saveRDS(
    list(family = distributions, children = child_distributions),
    file.path(dirs$kernel, "mechanistic_kernel_exact_distributions.rds"),
    compress = "xz"
  )
  status <- ptpv5_write_method_status(
    dirs,
    "exact_global_kernel",
    "complete",
    paste0(
      "Exact ridge-kernel tests and hierarchical closed testing completed ",
      "for ", nrow(family_results), " mechanistic families."
    ),
    nrow(family_results)
  )
  out <- list(
    family_results = family_results,
    children = children,
    gene_results = gene_results,
    directional = directional,
    fragility = fragility,
    distributions = list(
      family = distributions,
      children = child_distributions
    ),
    status = status
  )
  ptpv5_save_checkpoint(
    dirs, checkpoint_id, out, cfg, context$config_sha, context$input_sha
  )
  out
}
