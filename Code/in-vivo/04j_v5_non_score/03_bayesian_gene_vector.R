ptpv5_gene_ols <- function(Y, mouse_meta, factors = NULL) {
  ploidy <- as.numeric(mouse_meta$initial_ploidy == "4N")
  treated <- as.numeric(mouse_meta$treatment != "control")
  X <- cbind(
    intercept = 1,
    ploidy_4N = ploidy,
    gemcitabine = treated,
    interaction = ploidy * treated
  )
  if (!is.null(factors) && ncol(as.matrix(factors)) > 0L) {
    X <- cbind(X, as.matrix(factors))
  }
  qrX <- qr(X)
  rank <- qrX$rank
  if (rank < ncol(X)) {
    stop("Bayesian gene-vector design matrix is rank deficient.", call. = FALSE)
  }
  inv <- solve(crossprod(X))
  beta <- Y %*% X %*% inv
  fitted <- beta %*% t(X)
  residual <- Y - fitted
  df <- ncol(Y) - ncol(X)
  sigma2 <- rowSums(residual^2) / df
  i <- match("interaction", colnames(X))
  data.frame(
    gene = rownames(Y),
    estimate = beta[, i],
    standard_error = sqrt(pmax(sigma2 * inv[i, i], 0)),
    residual_sd = sqrt(pmax(sigma2, 0)),
    df_residual = df,
    stringsAsFactors = FALSE
  )
}

ptpv5_residual_factors <- function(Y, mouse_meta, cfg) {
  ploidy <- as.numeric(mouse_meta$initial_ploidy == "4N")
  reduced <- cbind(intercept = 1, ploidy_4N = ploidy)
  H <- reduced %*% solve(crossprod(reduced)) %*% t(reduced)
  residual <- Y %*% (diag(ncol(Y)) - H)
  pc <- stats::prcomp(t(residual), center = TRUE, scale. = FALSE)
  variance <- pc$sdev^2
  fraction <- variance / sum(variance)
  cumulative <- cumsum(fraction)
  target <- as.numeric(cfg$bayesian_gene_vector$residual_factor_variance_target)
  k <- which(cumulative >= target)[1L]
  if (!is.finite(k)) k <- 1L
  k <- min(
    as.integer(cfg$bayesian_gene_vector$maximum_residual_factors),
    k,
    ncol(pc$x)
  )
  full_design <- cbind(
    intercept = 1,
    ploidy_4N = ploidy,
    gemcitabine = as.numeric(mouse_meta$treatment != "control"),
    interaction = ploidy * as.numeric(mouse_meta$treatment != "control")
  )
  F <- pc$x[, seq_len(k), drop = FALSE]
  F <- qr.resid(qr(full_design), F)
  F <- scale(F, center = TRUE, scale = TRUE)
  F[!is.finite(F)] <- 0
  colnames(F) <- paste0("residual_factor_", seq_len(ncol(F)))
  list(
    scores = F,
    diagnostics = data.frame(
      factor = seq_along(variance),
      variance_fraction = fraction,
      cumulative_variance_fraction = cumulative,
      retained = seq_along(variance) <= k,
      stringsAsFactors = FALSE
    ),
    n_factors = k
  )
}

ptpv5_bayes_stan_code <- function() {
  "
  data {
    int<lower=2> N;
    vector[N] y;
    vector<lower=1e-6>[N] se;
    real<lower=0> mean_scale;
    real<lower=0> global_scale;
    real<lower=0> slab_scale;
    real<lower=1> slab_df;
    real<lower=1> mean_df;
    real<lower=1> deviation_df;
  }
  parameters {
    real mu;
    real<lower=0> tau;
    vector<lower=0>[N] lambda;
    vector[N] z;
    real<lower=0> caux;
  }
  transformed parameters {
    real<lower=0> c2 = square(slab_scale) * caux;
    vector[N] lambda_tilde =
      sqrt(c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)));
    vector[N] theta = mu + tau * lambda_tilde .* z;
  }
  model {
    mu ~ student_t(mean_df, 0, mean_scale);
    tau ~ cauchy(0, global_scale);
    lambda ~ cauchy(0, 1);
    z ~ student_t(deviation_df, 0, 1);
    caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
    y ~ normal(theta, se);
  }
  "
}

ptpv5_bayes_analytic_fallback <- function(y, se, prior_scale, rope) {
  w <- 1 / pmax(se^2, 1e-8)
  likelihood_var <- 1 / sum(w)
  likelihood_mean <- sum(w * y) / sum(w)
  posterior_var <- 1 / (1 / likelihood_var + 1 / prior_scale^2)
  posterior_mean <- posterior_var * likelihood_mean / likelihood_var
  posterior_sd <- sqrt(posterior_var)
  shrink <- pmax(0, 1 - se^2 / pmax(stats::var(y), se^2))
  theta <- posterior_mean + shrink * (y - posterior_mean)
  list(
    engine = "analytic_normal_normal_fallback",
    family = data.frame(
      posterior_mean = posterior_mean,
      posterior_median = posterior_mean,
      ci_low_95 = stats::qnorm(0.025, posterior_mean, posterior_sd),
      ci_high_95 = stats::qnorm(0.975, posterior_mean, posterior_sd),
      probability_positive = stats::pnorm(
        0, posterior_mean, posterior_sd, lower.tail = FALSE
      ),
      probability_negative = stats::pnorm(
        0, posterior_mean, posterior_sd, lower.tail = TRUE
      ),
      probability_outside_rope =
        stats::pnorm(-rope, posterior_mean, posterior_sd) +
        stats::pnorm(rope, posterior_mean, posterior_sd, lower.tail = FALSE),
      probability_above_rope = stats::pnorm(
        rope, posterior_mean, posterior_sd, lower.tail = FALSE
      ),
      probability_below_rope = stats::pnorm(
        -rope, posterior_mean, posterior_sd
      ),
      stringsAsFactors = FALSE
    ),
    genes = data.frame(
      posterior_mean = theta,
      posterior_median = theta,
      ci_low_95 = theta - 1.96 * se,
      ci_high_95 = theta + 1.96 * se,
      stringsAsFactors = FALSE
    ),
    diagnostics = data.frame(
      parameter = "mu",
      rhat = NA_real_,
      ess_bulk = NA_real_,
      ess_tail = NA_real_,
      divergences = NA_integer_,
      stringsAsFactors = FALSE
    )
  )
}

ptpv5_fit_bayes_family <- function(
  model,
  y,
  se,
  cfg,
  args,
  seed
) {
  rope <- as.numeric(cfg$bayesian_gene_vector$rope_half_width_residual_sd)
  stan_data <- list(
    N = length(y),
    y = as.numeric(y),
    se = pmax(as.numeric(se), 1e-4),
    mean_scale = as.numeric(cfg$bayesian_gene_vector$family_mean_prior_scale),
    global_scale = as.numeric(cfg$bayesian_gene_vector$global_scale_prior),
    slab_scale = as.numeric(cfg$bayesian_gene_vector$slab_scale),
    slab_df = as.numeric(cfg$bayesian_gene_vector$slab_df),
    mean_df = as.numeric(cfg$bayesian_gene_vector$family_mean_prior_df),
    deviation_df = as.numeric(cfg$bayesian_gene_vector$gene_deviation_df)
  )
  iterations <- as.integer(cfg$bayesian_gene_vector$iterations)
  warmup <- as.integer(cfg$bayesian_gene_vector$warmup)
  chains <- as.integer(cfg$bayesian_gene_vector$chains)
  if (isTRUE(args$smoke)) {
    iterations <- min(iterations, 400L)
    warmup <- min(warmup, 200L)
    chains <- min(chains, 2L)
  }
  fit <- try(
    rstan::sampling(
      model,
      data = stan_data,
      chains = chains,
      iter = iterations,
      warmup = warmup,
      thin = as.integer(cfg$bayesian_gene_vector$thin),
      seed = seed,
      cores = min(chains, as.integer(args$workers)),
      control = list(
        adapt_delta = as.numeric(cfg$bayesian_gene_vector$adapt_delta),
        max_treedepth = as.integer(cfg$bayesian_gene_vector$max_treedepth)
      ),
      refresh = 0
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    fallback <- ptpv5_bayes_analytic_fallback(
      y,
      se,
      as.numeric(cfg$bayesian_gene_vector$family_mean_prior_scale),
      rope
    )
    fallback$error <- as.character(fit)
    return(fallback)
  }
  draws <- rstan::extract(fit, pars = c("mu", "theta"), permuted = TRUE)
  mu <- as.numeric(draws$mu)
  theta <- as.matrix(draws$theta)
  family <- data.frame(
    posterior_mean = mean(mu),
    posterior_median = stats::median(mu),
    ci_low_95 = stats::quantile(mu, 0.025, names = FALSE),
    ci_high_95 = stats::quantile(mu, 0.975, names = FALSE),
    probability_positive = mean(mu > 0),
    probability_negative = mean(mu < 0),
    probability_outside_rope = mean(abs(mu) > rope),
    probability_above_rope = mean(mu > rope),
    probability_below_rope = mean(mu < -rope),
    stringsAsFactors = FALSE
  )
  genes <- data.frame(
    posterior_mean = colMeans(theta),
    posterior_median = apply(theta, 2L, stats::median),
    ci_low_95 = apply(theta, 2L, stats::quantile, probs = 0.025),
    ci_high_95 = apply(theta, 2L, stats::quantile, probs = 0.975),
    stringsAsFactors = FALSE
  )
  summary <- rstan::summary(fit, pars = c("mu", "tau"))$summary
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(
    sampler,
    function(x) sum(x[, "divergent__"]),
    numeric(1L)
  ))
  diagnostics <- data.frame(
    parameter = rownames(summary),
    rhat = summary[, "Rhat"],
    ess_bulk = summary[, "n_eff"],
    ess_tail = NA_real_,
    divergences = divergences,
    stringsAsFactors = FALSE
  )
  list(
    engine = "rstan_regularized_horseshoe",
    family = family,
    genes = genes,
    diagnostics = diagnostics,
    error = NA_character_
  )
}

ptpv5_prior_sensitivity <- function(y, se, scales) {
  w <- 1 / pmax(se^2, 1e-8)
  likelihood_var <- 1 / sum(w)
  likelihood_mean <- sum(w * y) / sum(w)
  ptpv5_bind_rows(lapply(scales, function(scale) {
    posterior_var <- 1 / (1 / likelihood_var + 1 / scale^2)
    posterior_mean <- posterior_var * likelihood_mean / likelihood_var
    data.frame(
      prior_scale = scale,
      posterior_mean_location_approximation = posterior_mean,
      posterior_sd_location_approximation = sqrt(posterior_var),
      approximation_role =
        "deterministic_location_submodel_prior_sensitivity",
      stringsAsFactors = FALSE
    )
  }))
}

ptpv5_family_location <- function(estimates, family_genes) {
  x <- estimates[estimates$gene_symbol %in% family_genes, , drop = FALSE]
  x <- x[is.finite(x$estimate) & is.finite(x$standard_error) &
           x$standard_error > 0, , drop = FALSE]
  if (nrow(x) < 2L) return(c(estimate = NA_real_, se = NA_real_, n = nrow(x)))
  w <- 1 / x$standard_error^2
  c(
    estimate = sum(w * x$estimate) / sum(w),
    se = sqrt(1 / sum(w)),
    n = nrow(x)
  )
}

ptpv5_family_effective_rank <- function(Y) {
  Y <- as.matrix(Y)
  Y <- sweep(Y, 1L, rowMeans(Y, na.rm = TRUE), "-")
  scale <- matrixStats::rowSds(Y, na.rm = TRUE)
  keep <- is.finite(scale) & scale > 0
  Y <- Y[keep, , drop = FALSE]
  if (nrow(Y) < 2L) return(1)
  Y <- sweep(Y, 1L, scale[keep], "/")
  singular <- svd(Y, nu = 0L, nv = 0L)$d^2
  singular <- singular[is.finite(singular) & singular > 1e-10]
  if (!length(singular)) return(1)
  max(1, (sum(singular)^2) / sum(singular^2))
}

ptpv5_run_bayesian_gene_vector <- function(context) {
  cfg <- context$cfg
  args <- context$args
  dirs <- context$dirs
  bundle <- context$bundle
  checkpoint_id <- "method03_bayesian_gene_vector"
  cached <- ptpv5_load_checkpoint(
    dirs, checkpoint_id, cfg, context$config_sha, context$input_sha
  )
  if (!is.null(cached)) return(cached)
  Y <- bundle$mouse$expression
  gene_sd <- matrixStats::rowSds(Y, na.rm = TRUE)
  fallback_sd <- stats::median(gene_sd[is.finite(gene_sd) & gene_sd > 0])
  gene_sd[!is.finite(gene_sd) | gene_sd <= 0] <- fallback_sd
  Y <- sweep(Y, 1L, gene_sd, "/")
  mouse_meta <- bundle$mouse$metadata[
    match(colnames(Y), bundle$mouse$metadata$sample_id),
    ,
    drop = FALSE
  ]
  factors <- ptpv5_residual_factors(Y, mouse_meta, cfg)
  estimates <- ptpv5_gene_ols(Y, mouse_meta, factors$scores)
  estimates$gene_symbol <- bundle$gene_symbols[
    match(estimates$gene, bundle$gene_ids)
  ]
  ptp_write_csv(
    estimates,
    file.path(dirs$bayes, "factor_adjusted_gene_interaction_estimates.csv")
  )
  ptp_write_csv(
    factors$diagnostics,
    file.path(dirs$bayes, "residual_factor_diagnostics.csv")
  )
  sets <- ptpv5_frozen_sets(bundle, args)
  family_meta <- sets$families
  model <- NULL
  compile_error <- NA_character_
  if (requireNamespace("rstan", quietly = TRUE)) {
    options(mc.cores = min(
      as.integer(args$workers),
      as.integer(cfg$bayesian_gene_vector$chains)
    ))
    model_try <- try(
      rstan::stan_model(
        model_code = ptpv5_bayes_stan_code(),
        model_name = "ptpv5_regularized_horseshoe_gene_vector",
        auto_write = TRUE
      ),
      silent = TRUE
    )
    if (!inherits(model_try, "try-error")) {
      model <- model_try
    } else {
      compile_error <- as.character(model_try)
    }
  } else {
    compile_error <- "rstan package unavailable"
  }
  family_rows <- list()
  gene_rows <- list()
  diagnostic_rows <- list()
  sensitivity_rows <- list()
  fit_dir <- ptp_ensure_dir(file.path(dirs$bayes, "family_summary_checkpoints"))
  for (i in seq_len(nrow(family_meta))) {
    family_id <- family_meta$response_family[[i]]
    family_genes <- ptpv5_split_genes(family_meta$union_genes[[i]])
    x <- estimates[
      estimates$gene_symbol %in% family_genes &
        is.finite(estimates$estimate) &
        is.finite(estimates$standard_error) &
        estimates$standard_error > 0,
      ,
      drop = FALSE
    ]
    if (isTRUE(args$smoke) && nrow(x) > as.integer(args$smoke_genes)) {
      x <- head(x, as.integer(args$smoke_genes))
    }
    if (nrow(x) < 2L) next
    family_row_index <- match(x$gene, rownames(Y))
    effective_rank <- ptpv5_family_effective_rank(
      Y[family_row_index, , drop = FALSE]
    )
    correlation_inflation <- sqrt(nrow(x) / effective_rank)
    x$model_standard_error <- x$standard_error * correlation_inflation
    family_path <- file.path(fit_dir, paste0(family_id, ".rds"))
    old <- if (file.exists(family_path)) try(readRDS(family_path), silent = TRUE) else NULL
    valid_old <- is.list(old) &&
      identical(old$config_sha, context$config_sha) &&
      identical(old$input_sha, context$input_sha) &&
      identical(old$genes, x$gene)
    if (valid_old) {
      fit <- old$fit
    } else {
      fit <- if (is.null(model)) {
        z <- ptpv5_bayes_analytic_fallback(
          x$estimate,
          x$standard_error,
          as.numeric(cfg$bayesian_gene_vector$family_mean_prior_scale),
          as.numeric(cfg$bayesian_gene_vector$rope_half_width_residual_sd)
        )
        z$error <- compile_error
        z
      } else {
        ptpv5_fit_bayes_family(
          model,
          x$estimate,
          x$model_standard_error,
          cfg,
          args,
          as.integer(cfg$reproducibility$bayesian_seed) + i
        )
      }
      saveRDS(
        list(
          config_sha = context$config_sha,
          input_sha = context$input_sha,
          genes = x$gene,
          fit = fit
        ),
        family_path,
        compress = "xz"
      )
    }
    family_row <- cbind(
      data.frame(
        response_family = family_id,
        response_family_label = family_meta$response_family_label[[i]],
        tiers = family_meta$tiers[[i]],
        n_genes = nrow(x),
        family_expression_effective_rank = effective_rank,
        gene_likelihood_se_inflation = correlation_inflation,
        engine = fit$engine,
        model_error = fit$error %||% NA_character_,
        stringsAsFactors = FALSE
      ),
      fit$family
    )
    family_rows[[family_id]] <- family_row
    genes <- cbind(
      data.frame(
        response_family = family_id,
        gene = x$gene,
        gene_symbol = x$gene_symbol,
        observed_estimate = x$estimate,
        observed_standard_error = x$standard_error,
        model_standard_error = x$model_standard_error,
        stringsAsFactors = FALSE
      ),
      fit$genes
    )
    gene_rows[[family_id]] <- genes
    diag <- fit$diagnostics
    diag$response_family <- family_id
    diagnostic_rows[[family_id]] <- diag
    sens <- ptpv5_prior_sensitivity(
      x$estimate,
      x$model_standard_error,
      as.numeric(cfg$bayesian_gene_vector$prior_sensitivity_scales)
    )
    sens$response_family <- family_id
    sensitivity_rows[[family_id]] <- sens
    ptpv5_message(
      "Bayesian family completed: ", family_id, " [", fit$engine, "]",
      log_file = context$log_file
    )
  }
  family_results <- ptpv5_bind_rows(family_rows)
  gene_results <- ptpv5_bind_rows(gene_rows)
  diagnostics <- ptpv5_bind_rows(diagnostic_rows)
  sensitivity <- ptpv5_bind_rows(sensitivity_rows)

  loo_rows <- list()
  for (j in seq_len(ncol(Y))) {
    keep <- setdiff(seq_len(ncol(Y)), j)
    local_factors <- factors$scores[keep, , drop = FALSE]
    local <- ptpv5_gene_ols(
      Y[, keep, drop = FALSE],
      mouse_meta[keep, , drop = FALSE],
      local_factors
    )
    local$gene_symbol <- bundle$gene_symbols[match(local$gene, bundle$gene_ids)]
    for (i in seq_len(nrow(family_meta))) {
      loc <- ptpv5_family_location(
        local,
        ptpv5_split_genes(family_meta$union_genes[[i]])
      )
      loo_rows[[length(loo_rows) + 1L]] <- data.frame(
        omitted_sample = colnames(Y)[[j]],
        response_family = family_meta$response_family[[i]],
        location_estimate = loc[["estimate"]],
        location_se = loc[["se"]],
        n_genes = loc[["n"]],
        diagnostic_role =
          "leave_one_mouse_out_inverse_variance_location_diagnostic",
        stringsAsFactors = FALSE
      )
    }
  }
  loo <- ptpv5_bind_rows(loo_rows)
  ptp_write_csv(
    family_results,
    file.path(dirs$bayes, "bayesian_family_posterior_summary.csv")
  )
  ptp_write_csv(
    gene_results,
    file.path(dirs$bayes, "bayesian_gene_posterior_summary.csv")
  )
  ptp_write_csv(
    diagnostics,
    file.path(dirs$bayes, "bayesian_sampler_diagnostics.csv")
  )
  ptp_write_csv(
    sensitivity,
    file.path(dirs$bayes, "bayesian_prior_sensitivity.csv")
  )
  ptp_write_csv(
    loo,
    file.path(dirs$bayes, "bayesian_leave_one_mouse_out_diagnostic.csv")
  )
  all_rstan <- nrow(family_results) > 0L &&
    all(family_results$engine == "rstan_regularized_horseshoe")
  status <- ptpv5_write_method_status(
    dirs,
    "bayesian_gene_vector",
    if (all_rstan) "complete" else "complete_with_declared_fallback",
    if (all_rstan) {
      paste0(
        "Regularized-horseshoe posterior summaries completed for ",
        nrow(family_results), " response families."
      )
    } else {
      paste0(
        "Bayesian family analysis completed for ", nrow(family_results),
        " families; at least one family used the explicitly recorded ",
        "normal-normal fallback because Stan compilation or sampling failed."
      )
    },
    nrow(family_results)
  )
  out <- list(
    family_results = family_results,
    gene_results = gene_results,
    diagnostics = diagnostics,
    sensitivity = sensitivity,
    loo = loo,
    estimates = estimates,
    factors = factors,
    status = status
  )
  ptpv5_save_checkpoint(
    dirs, checkpoint_id, out, cfg, context$config_sha, context$input_sha
  )
  out
}
