ptpv5_rank_p_vector <- function(statistic, larger = TRUE) {
  x <- as.numeric(statistic)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  values <- if (larger) -x[ok] else x[ok]
  ranks <- rank(values, ties.method = "max")
  out[ok] <- (1 + ranks) / (1 + sum(ok))
  out
}

ptpv5_binomial_interval <- function(success, total, level = 0.95) {
  if (!is.finite(total) || total < 1L) return(c(NA_real_, NA_real_))
  stats::binom.test(success, total, conf.level = level)$conf.int
}

ptpv5_collect_exact_null_p <- function(results) {
  rows <- list()
  m1 <- results$exact_ranked_enrichment$distribution$weighted_es
  for (i in seq_len(nrow(m1))) {
    rows[[length(rows) + 1L]] <- data.frame(
      method_id = "exact_ranked_enrichment",
      test_id = rownames(m1)[[i]],
      assignment_index = seq_len(ncol(m1)),
      exact_p = ptpv5_rank_p_vector(abs(m1[i, ])),
      stringsAsFactors = FALSE
    )
  }
  m2 <- results$state_decomposition$distributions
  for (id in names(m2)) {
    x <- m2[[id]]
    if (grepl("__DA$", id)) {
      p <- ptpv5_rank_p_vector(x$statistic)
      rows[[length(rows) + 1L]] <- data.frame(
        method_id = "state_decomposition",
        test_id = id,
        assignment_index = seq_along(p),
        exact_p = p,
        stringsAsFactors = FALSE
      )
    } else {
      for (component in names(x)) {
        p <- x[[component]]$adaptive_assignment_p
        rows[[length(rows) + 1L]] <- data.frame(
          method_id = "state_decomposition",
          test_id = paste(id, component, sep = "__"),
          assignment_index = seq_along(p),
          exact_p = p,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  m4 <- results$exact_global_kernel$distributions$family
  target <- results$exact_global_kernel$family_results$response_family
  for (id in intersect(target, names(m4))) {
    p <- ptpv5_rank_p_vector(m4[[id]])
    rows[[length(rows) + 1L]] <- data.frame(
      method_id = "exact_global_kernel",
      test_id = id,
      assignment_index = seq_along(p),
      exact_p = p,
      stringsAsFactors = FALSE
    )
  }
  m5 <- results$gene_wise_trajectory$family_distributions
  for (id in names(m5)) {
    p <- m5[[id]]$adaptive_assignment_p
    rows[[length(rows) + 1L]] <- data.frame(
      method_id = "gene_wise_trajectory",
      test_id = id,
      assignment_index = seq_along(p),
      exact_p = p,
      stringsAsFactors = FALSE
    )
  }
  ptpv5_bind_rows(rows)
}

ptpv5_bayes_sbc <- function(context, bayes_result) {
  cfg <- context$cfg
  set.seed(as.integer(cfg$reproducibility$simulation_seed) + 303L)
  R <- as.integer(cfg$bayesian_gene_vector$sbc_replicates)
  if (isTRUE(context$args$smoke)) R <- min(R, 6L)
  se_pool <- bayes_result$estimates$standard_error
  se_pool <- se_pool[is.finite(se_pool) & se_pool > 0]
  n <- min(80L, max(20L, length(se_pool)))
  prior_scale <- as.numeric(
    cfg$bayesian_gene_vector$family_mean_prior_scale
  )
  rows <- vector("list", R)
  for (r in seq_len(R)) {
    mu <- stats::rt(1L, df = as.numeric(
      cfg$bayesian_gene_vector$family_mean_prior_df
    )) * prior_scale
    se <- sample(se_pool, n, replace = TRUE)
    deviation <- stats::rt(
      n,
      df = as.numeric(cfg$bayesian_gene_vector$gene_deviation_df)
    ) * as.numeric(cfg$bayesian_gene_vector$global_scale_prior)
    y <- stats::rnorm(n, mu + deviation, se)
    fit <- ptpv5_bayes_analytic_fallback(
      y,
      se,
      prior_scale,
      as.numeric(cfg$bayesian_gene_vector$rope_half_width_residual_sd)
    )
    mean <- fit$family$posterior_mean
    low <- fit$family$ci_low_95
    high <- fit$family$ci_high_95
    sd <- (high - low) / (2 * 1.96)
    draws <- stats::rnorm(1000L, mean, sd)
    rows[[r]] <- data.frame(
      replicate = r,
      true_family_mean = mu,
      posterior_mean = mean,
      covered_95 = mu >= low & mu <= high,
      posterior_rank_fraction = mean(draws < mu),
      calibration_scope =
        "family_location_submodel_representative_SBC",
      stringsAsFactors = FALSE
    )
  }
  ptpv5_bind_rows(rows)
}

ptpv5_power_benchmark <- function(context, results) {
  cfg <- context$cfg
  set.seed(as.integer(cfg$reproducibility$simulation_seed) + 707L)
  R <- as.integer(cfg$simulation$power_replicates)
  if (isTRUE(context$args$smoke)) R <- min(R, 40L)
  shifts <- cfg$simulation$empirical_null_shift_benchmark
  routes <- data.frame(
    method_id = c(
      rep("exact_ranked_enrichment", 3),
      rep("state_decomposition", 2),
      rep("bayesian_gene_vector", 3),
      rep("exact_global_kernel", 3),
      rep("gene_wise_trajectory", 2)
    ),
    scenario = c(
      "coherent_family_shift", "sparse_family_shift",
      "mixed_direction_shift", "composition_only_shift",
      "within_state_only_shift", "coherent_family_shift",
      "sparse_family_shift", "mixed_direction_shift",
      "coherent_family_shift", "sparse_family_shift",
      "mixed_direction_shift", "local_trajectory_peak",
      "early_late_reversal"
    ),
    shift_sd = c(
      shifts$coherent_family_shift_sd,
      shifts$sparse_family_shift_sd,
      shifts$mixed_direction_shift_sd,
      shifts$composition_only_shift_sd,
      shifts$within_state_only_shift_sd,
      shifts$coherent_family_shift_sd,
      shifts$sparse_family_shift_sd,
      shifts$mixed_direction_shift_sd,
      shifts$coherent_family_shift_sd,
      shifts$sparse_family_shift_sd,
      shifts$mixed_direction_shift_sd,
      shifts$local_trajectory_peak_sd,
      shifts$early_late_reversal_sd
    ),
    stringsAsFactors = FALSE
  )
  null_sources <- list(
    exact_ranked_enrichment = abs(
      as.numeric(results$exact_ranked_enrichment$distribution$weighted_es)
    ),
    state_decomposition = unlist(lapply(
      results$state_decomposition$distributions,
      function(x) {
        if (!is.null(x$statistic)) return(as.numeric(x$statistic))
        unlist(lapply(x, function(y) {
          if (!is.null(y$statistic)) as.numeric(y$statistic) else numeric()
        }))
      }
    )),
    bayesian_gene_vector = stats::rnorm(4900L),
    exact_global_kernel = unlist(
      results$exact_global_kernel$distributions$family[
        results$exact_global_kernel$family_results$response_family
      ]
    ),
    gene_wise_trajectory = unlist(lapply(
      results$gene_wise_trajectory$family_distributions,
      function(x) as.numeric(x$statistic)
    ))
  )
  rows <- lapply(seq_len(nrow(routes)), function(i) {
    id <- routes$method_id[[i]]
    null <- as.numeric(null_sources[[id]])
    null <- null[is.finite(null)]
    if (length(null) < 20L) {
      return(data.frame(
        routes[i, ],
        replicates = R,
        power_at_0_05 = NA_real_,
        benchmark_role = shifts$role,
        stringsAsFactors = FALSE
      ))
    }
    z <- as.numeric(scale(null))
    z <- z[is.finite(z)]
    threshold <- stats::quantile(z, 0.95, names = FALSE)
    simulated <- sample(z, R, replace = TRUE) + routes$shift_sd[[i]]
    data.frame(
      routes[i, ],
      replicates = R,
      power_at_0_05 = mean(simulated >= threshold),
      benchmark_role = shifts$role,
      stringsAsFactors = FALSE
    )
  })
  ptpv5_bind_rows(rows)
}

ptpv5_build_family_evidence <- function(context, results) {
  families <- context$bundle$families
  rows <- list()
  m1 <- results$exact_ranked_enrichment$results
  m1 <- m1[m1$set_type == "response_family_union", , drop = FALSE]
  rows[[1L]] <- data.frame(
    response_family = m1$response_family,
    method_id = "exact_ranked_enrichment",
    evidence_value = m1$weighted_exact_p,
    multiplicity_value = m1$weighted_fdr_all,
    evidence_type = "exact_p_and_BH",
    pass_predefined_threshold = m1$weighted_fdr_all <= 0.10,
    stringsAsFactors = FALSE
  )
  m2 <- results$state_decomposition$decomposition
  m2 <- m2[
    m2$partition != "cluster_phase" &
      m2$decomposition_component %in%
        c("within_state", "composition"),
    ,
    drop = FALSE
  ]
  if (nrow(m2)) {
    split_m2 <- split(m2, m2$response_family)
    rows[[2L]] <- ptpv5_bind_rows(lapply(split_m2, function(x) {
      j <- which.min(x$adaptive_exact_p)
      data.frame(
        response_family = x$response_family[[j]],
        method_id = "state_decomposition",
        evidence_value = x$adaptive_exact_p[[j]],
        multiplicity_value = x$decomposition_fdr_all[[j]],
        evidence_type = paste(
          "minimum_predefined_partition_component",
          x$partition[[j]], x$decomposition_component[[j]], sep = ":"
        ),
        pass_predefined_threshold =
          x$decomposition_fdr_all[[j]] <= 0.10,
        stringsAsFactors = FALSE
      )
    }))
  }
  m3 <- results$bayesian_gene_vector$family_results
  rows[[3L]] <- data.frame(
    response_family = m3$response_family,
    method_id = "bayesian_gene_vector",
    evidence_value = m3$probability_outside_rope,
    multiplicity_value = NA_real_,
    evidence_type = "posterior_probability_outside_ROPE",
    pass_predefined_threshold = m3$probability_outside_rope >= 0.90,
    stringsAsFactors = FALSE
  )
  m4 <- results$exact_global_kernel$family_results
  rows[[4L]] <- data.frame(
    response_family = m4$response_family,
    method_id = "exact_global_kernel",
    evidence_value = m4$exact_kernel_p,
    multiplicity_value = m4$holm_five_families,
    evidence_type = "exact_p_and_Holm_five_mechanistic_families",
    pass_predefined_threshold = m4$family_gate_pass,
    stringsAsFactors = FALSE
  )
  m5 <- results$gene_wise_trajectory$family_results
  if (nrow(m5)) {
    rows[[5L]] <- data.frame(
      response_family = m5$response_family,
      method_id = "gene_wise_trajectory",
      evidence_value = m5$adaptive_exact_p,
      multiplicity_value = m5$trajectory_fdr_all,
      evidence_type = "whole_trajectory_exact_p_and_BH",
      pass_predefined_threshold = m5$trajectory_fdr_all <= 0.10,
      stringsAsFactors = FALSE
    )
  }
  long <- ptpv5_bind_rows(rows)
  long <- merge(
    families[, c(
      "response_family", "response_family_label", "tiers", "n_union_genes"
    )],
    long,
    by = "response_family",
    all.x = TRUE,
    sort = FALSE
  )
  summary <- ptpv5_bind_rows(lapply(
    split(long, long$response_family),
    function(x) {
      passes <- sum(x$pass_predefined_threshold %in% TRUE, na.rm = TRUE)
      data.frame(
        response_family = x$response_family[[1L]],
        response_family_label = x$response_family_label[[1L]],
        tiers = x$tiers[[1L]],
        n_methods_available = sum(!is.na(x$method_id)),
        n_methods_passing = passes,
        synthesis_class = if (passes >= 2L) {
          "concordant_non_score_support"
        } else if (passes == 1L) {
          "isolated_method_support_requires_caution"
        } else {
          "no_predefined_non_score_threshold_pass"
        },
        stringsAsFactors = FALSE
      )
    }
  ))
  list(long = long, summary = summary)
}

ptpv5_v4_comparison <- function(context, evidence) {
  path <- file.path(
    context$cfg$analysis$v4_result_root,
    "primary_initial_ploidy",
    "response_family_hierarchical",
    "response_family_hierarchical_results.csv"
  )
  if (!file.exists(path)) return(data.frame())
  v4 <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- intersect(
    c(
      "response_family", "response_family_label", "family_omnibus_p",
      "fdr_all_families", "fdr_focused_11", "fdr_tier3_3",
      "fdr_controls"
    ),
    names(v4)
  )
  v4 <- v4[, keep, drop = FALSE]
  merge(
    v4,
    evidence$summary,
    by = intersect(c("response_family", "response_family_label"), names(v4)),
    all = TRUE,
    sort = FALSE
  )
}

ptpv5_run_simulation_and_synthesis <- function(context, results) {
  cfg <- context$cfg
  dirs <- context$dirs
  checkpoint_id <- "method06_simulation_and_synthesis"
  cached <- ptpv5_load_checkpoint(
    dirs, checkpoint_id, cfg, context$config_sha, context$input_sha
  )
  if (!is.null(cached)) return(cached)
  null_p <- ptpv5_collect_exact_null_p(results)
  set.seed(as.integer(cfg$reproducibility$simulation_seed))
  R <- as.integer(cfg$simulation$null_replicates)
  if (isTRUE(context$args$smoke)) R <- min(R, 100L)
  sampled <- ptpv5_bind_rows(lapply(
    split(null_p, interaction(null_p$method_id, null_p$test_id, drop = TRUE)),
    function(x) {
      idx <- sample(seq_len(nrow(x)), R, replace = TRUE)
      x[idx, , drop = FALSE]
    }
  ))
  alpha <- as.numeric(cfg$simulation$alpha)
  calibration <- ptpv5_bind_rows(lapply(
    split(sampled, sampled$method_id),
    function(x) {
      reject <- x$exact_p <= alpha
      interval <- ptpv5_binomial_interval(sum(reject, na.rm = TRUE), sum(is.finite(reject)))
      data.frame(
        method_id = x$method_id[[1L]],
        null_draws = sum(is.finite(reject)),
        empirical_type1_at_0_05 = mean(reject, na.rm = TRUE),
        ci_low_95 = interval[[1L]],
        ci_high_95 = interval[[2L]],
        calibration_pass = mean(reject, na.rm = TRUE) <=
          as.numeric(cfg$simulation$type1_reference_interval[[2L]]),
        calibration_basis =
          "pseudo-observed assignments sampled from the exact randomization distribution",
        stringsAsFactors = FALSE
      )
    }
  ))
  sbc <- ptpv5_bayes_sbc(context, results$bayesian_gene_vector)
  coverage <- mean(sbc$covered_95)
  bayes_cal <- data.frame(
    method_id = "bayesian_gene_vector",
    null_draws = nrow(sbc),
    empirical_type1_at_0_05 = NA_real_,
    ci_low_95 = NA_real_,
    ci_high_95 = NA_real_,
    calibration_pass = coverage >=
      as.numeric(cfg$simulation$coverage_reference_interval[[1L]]) &&
      coverage <=
        as.numeric(cfg$simulation$coverage_reference_interval[[2L]]),
    calibration_basis = paste0(
      "representative family-location SBC; 95% coverage=", signif(coverage, 4)
    ),
    stringsAsFactors = FALSE
  )
  calibration <- ptpv5_bind_rows(list(calibration, bayes_cal))
  power <- ptpv5_power_benchmark(context, results)
  evidence <- ptpv5_build_family_evidence(context, results)
  comparison <- ptpv5_v4_comparison(context, evidence)
  method_summary <- ptpv5_bind_rows(lapply(names(results), function(id) {
    status <- results[[id]]$status
    data.frame(
      method_id = id,
      status = status$status[[1L]],
      note = status$note[[1L]],
      n_tests = status$n_tests[[1L]],
      stringsAsFactors = FALSE
    )
  }))
  ptp_write_csv(
    calibration,
    file.path(dirs$simulation, "method_calibration_summary.csv")
  )
  ptp_write_csv(
    sbc,
    file.path(dirs$simulation, "bayesian_location_sbc.csv")
  )
  ptp_write_csv(
    power,
    file.path(dirs$simulation, "generic_power_stress_test.csv")
  )
  ptp_write_csv(
    evidence$long,
    file.path(dirs$synthesis, "family_evidence_long.csv")
  )
  ptp_write_csv(
    evidence$summary,
    file.path(dirs$synthesis, "family_evidence_synthesis.csv")
  )
  ptp_write_csv(
    comparison,
    file.path(dirs$synthesis, "v4_v5_family_conclusion_comparison.csv")
  )
  ptp_write_csv(
    method_summary,
    file.path(dirs$synthesis, "five_method_status_summary.csv")
  )
  status <- ptpv5_write_method_status(
    dirs,
    "simulation_and_synthesis",
    "complete",
    paste0(
      "Calibration, generic power stress tests, and family-level synthesis ",
      "completed for five non-score methods."
    ),
    nrow(calibration) + nrow(evidence$long)
  )
  out <- list(
    calibration = calibration,
    sbc = sbc,
    power = power,
    evidence = evidence,
    v4_comparison = comparison,
    method_summary = method_summary,
    status = status
  )
  ptpv5_save_checkpoint(
    dirs, checkpoint_id, out, cfg, context$config_sha, context$input_sha
  )
  out
}
