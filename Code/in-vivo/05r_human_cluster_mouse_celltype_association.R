#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 280)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

cluster_order <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  if (all(is.finite(numeric_x))) x[order(numeric_x)] else sort(x)
}

p_value_stars <- function(p) {
  ifelse(
    !is.finite(p), "NA",
    ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
  )
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  out <- lapply(seq_along(x), function(i) {
    rest <- all_permutations(x[-i])
    rbind(x[[i]], rest)
  })
  do.call(cbind, out)
}

build_within_dose_permutations <- function(sample_metadata_one_ploidy) {
  z <- sample_metadata_one_ploidy[order(sample_metadata_one_ploidy$sample_order), , drop = FALSE]
  blocks <- split(seq_len(nrow(z)), as.character(z$Dose))
  block_permutations <- lapply(blocks, all_permutations)
  choice_grid <- expand.grid(lapply(block_permutations, function(x) seq_len(ncol(x))), KEEP.OUT.ATTRS = FALSE)
  permutations <- matrix(rep(seq_len(nrow(z)), nrow(choice_grid)), nrow = nrow(z))
  for (j in seq_len(nrow(choice_grid))) {
    for (b in seq_along(blocks)) {
      permutations[blocks[[b]], j] <- block_permutations[[b]][, choice_grid[j, b]]
    }
  }
  rownames(permutations) <- z$sample
  permutations
}

add_clr <- function(counts, category_column, clr_column) {
  samples <- unique(counts$sample[order(counts$sample_order)])
  categories <- unique(as.character(counts[[category_column]])[order(counts$category_order)])
  count_matrix <- matrix(0, nrow = length(samples), ncol = length(categories), dimnames = list(samples, categories))
  count_matrix[cbind(match(counts$sample, samples), match(as.character(counts[[category_column]]), categories))] <- counts$n_cells
  adjusted <- count_matrix + 0.5
  composition <- adjusted / rowSums(adjusted)
  log_composition <- log(composition)
  clr <- log_composition - rowMeans(log_composition)
  counts[[clr_column]] <- clr[cbind(match(counts$sample, samples), match(as.character(counts[[category_column]]), categories))]
  counts
}

dose_adjusted_spearman <- function(x, y, dose) {
  if (length(unique(x)) < 3L || length(unique(y)) < 3L) return(NA_real_)
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  x_residual <- stats::residuals(stats::lm(rx ~ factor(dose)))
  y_residual <- stats::residuals(stats::lm(ry ~ factor(dose)))
  stats::cor(x_residual, y_residual, method = "pearson")
}

permutation_partial_spearman <- function(x, y, dose, permutations) {
  observed <- dose_adjusted_spearman(x, y, dose)
  if (!is.finite(observed)) return(list(rho = NA_real_, p_value = NA_real_))
  permuted <- vapply(seq_len(ncol(permutations)), function(j) {
    dose_adjusted_spearman(x, y[permutations[, j]], dose)
  }, numeric(1))
  list(
    rho = observed,
    p_value = mean(abs(permuted) >= abs(observed) - sqrt(.Machine$double.eps), na.rm = TRUE)
  )
}

run_association_family <- function(human_counts, mouse_counts, mouse_category_column, outcome_column, outcome_family, resolution, sample_metadata) {
  ploidy_levels <- c("2N", "4N")
  human_levels <- cluster_order(human_counts$human_cluster)
  mouse_levels <- unique(as.character(mouse_counts[[mouse_category_column]])[order(mouse_counts$category_order)])
  results <- list()
  k <- 0L
  for (ploidy_i in ploidy_levels) {
    metadata_i <- sample_metadata[as.character(sample_metadata$Ploidy) == ploidy_i, , drop = FALSE]
    metadata_i <- metadata_i[order(metadata_i$sample_order), , drop = FALSE]
    permutations <- build_within_dose_permutations(metadata_i)
    if (nrow(metadata_i) != 8L || ncol(permutations) != 96L) {
      stop("Each Ploidy-specific association requires eight samples and 96 within-Dose permutations.", call. = FALSE)
    }
    for (human_i in human_levels) {
      h <- human_counts[as.character(human_counts$human_cluster) == human_i, , drop = FALSE]
      h <- h[match(metadata_i$sample, h$sample), , drop = FALSE]
      for (mouse_i in mouse_levels) {
        m <- mouse_counts[as.character(mouse_counts[[mouse_category_column]]) == mouse_i, , drop = FALSE]
        m <- m[match(metadata_i$sample, m$sample), , drop = FALSE]
        if (anyNA(h$sample) || anyNA(m$sample)) stop("Association table is missing a required sample.", call. = FALSE)
        x <- as.numeric(h$human_clr)
        y <- as.numeric(m[[outcome_column]])
        human_nonzero <- sum(h$n_cells > 0L)
        mouse_nonzero <- sum(m$n_cells > 0L)
        testable <- human_nonzero >= 4L && mouse_nonzero >= 4L && length(unique(x)) >= 3L && length(unique(y)) >= 3L
        test <- if (testable) permutation_partial_spearman(x, y, metadata_i$Dose, permutations) else list(rho = NA_real_, p_value = NA_real_)
        k <- k + 1L
        results[[k]] <- data.frame(
          resolution = resolution,
          outcome_family = outcome_family,
          Ploidy = ploidy_i,
          human_cluster = human_i,
          mouse_category = mouse_i,
          n_samples = nrow(metadata_i),
          human_nonzero_samples = human_nonzero,
          mouse_nonzero_samples = mouse_nonzero,
          unadjusted_spearman_rho = if (testable) stats::cor(x, y, method = "spearman") else NA_real_,
          dose_adjusted_partial_spearman_rho = test$rho,
          n_exact_within_dose_permutations = ncol(permutations),
          p_value = test$p_value,
          test_status = if (testable) "tested" else "insufficient_prevalence_or_variation",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, results)
  out$p_value_bh <- NA_real_
  for (ploidy_i in ploidy_levels) {
    idx <- which(out$Ploidy == ploidy_i & out$test_status == "tested")
    out$p_value_bh[idx] <- stats::p.adjust(out$p_value[idx], method = "BH")
  }
  out$significance_fdr <- p_value_stars(out$p_value_bh)
  out$human_cluster_order <- match(out$human_cluster, human_levels)
  out$mouse_category_order <- match(out$mouse_category, mouse_levels)
  out[order(match(out$Ploidy, ploidy_levels), out$mouse_category_order, out$human_cluster_order), , drop = FALSE]
}

run_interaction_family <- function(human_counts, mouse_counts, mouse_category_column, outcome_column, outcome_family, resolution, associations, sample_metadata) {
  human_levels <- cluster_order(human_counts$human_cluster)
  mouse_levels <- unique(as.character(mouse_counts[[mouse_category_column]])[order(mouse_counts$category_order)])
  rows <- list()
  k <- 0L
  for (human_i in human_levels) {
    h <- human_counts[as.character(human_counts$human_cluster) == human_i, , drop = FALSE]
    h <- h[match(sample_metadata$sample, h$sample), , drop = FALSE]
    for (mouse_i in mouse_levels) {
      m <- mouse_counts[as.character(mouse_counts[[mouse_category_column]]) == mouse_i, , drop = FALSE]
      m <- m[match(sample_metadata$sample, m$sample), , drop = FALSE]
      analysis_data <- data.frame(
        x = h$human_clr,
        y = m[[outcome_column]],
        Ploidy = factor(as.character(sample_metadata$Ploidy), levels = c("2N", "4N")),
        Dose = factor(as.character(sample_metadata$Dose)),
        stringsAsFactors = FALSE
      )
      testable_2n <- sum(h$n_cells[analysis_data$Ploidy == "2N"] > 0L) >= 4L && sum(m$n_cells[analysis_data$Ploidy == "2N"] > 0L) >= 4L
      testable_4n <- sum(h$n_cells[analysis_data$Ploidy == "4N"] > 0L) >= 4L && sum(m$n_cells[analysis_data$Ploidy == "4N"] > 0L) >= 4L
      testable <- testable_2n && testable_4n && length(unique(analysis_data$x)) >= 3L && length(unique(analysis_data$y)) >= 3L
      beta <- standard_error <- p_value <- NA_real_
      if (testable) {
        fit <- stats::lm(y ~ x * Ploidy + Dose, data = analysis_data)
        coefficient_table <- summary(fit)$coefficients
        coefficient_name <- "x:Ploidy4N"
        if (coefficient_name %in% rownames(coefficient_table)) {
          beta <- unname(coefficient_table[coefficient_name, "Estimate"])
          standard_error <- unname(coefficient_table[coefficient_name, "Std. Error"])
          p_value <- unname(coefficient_table[coefficient_name, "Pr(>|t|)"])
        }
      }
      rho_2n <- associations$dose_adjusted_partial_spearman_rho[
        associations$Ploidy == "2N" & associations$human_cluster == human_i & associations$mouse_category == mouse_i
      ]
      rho_4n <- associations$dose_adjusted_partial_spearman_rho[
        associations$Ploidy == "4N" & associations$human_cluster == human_i & associations$mouse_category == mouse_i
      ]
      k <- k + 1L
      rows[[k]] <- data.frame(
        resolution = resolution,
        outcome_family = outcome_family,
        human_cluster = human_i,
        mouse_category = mouse_i,
        n_samples = nrow(sample_metadata),
        partial_spearman_rho_2N = if (length(rho_2n) == 1L) rho_2n else NA_real_,
        partial_spearman_rho_4N = if (length(rho_4n) == 1L) rho_4n else NA_real_,
        difference_rho_4N_minus_2N = if (length(rho_2n) == 1L && length(rho_4n) == 1L) rho_4n - rho_2n else NA_real_,
        interaction_beta_4N_minus_2N = beta,
        interaction_standard_error = standard_error,
        interaction_p_value = p_value,
        interaction_test = "linear model on transformed values: mouse ~ human * Ploidy + Dose; exploratory because n=16",
        test_status = if (testable && is.finite(p_value)) "tested_exploratory" else "insufficient_prevalence_variation_or_model_rank",
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  tested <- out$test_status == "tested_exploratory"
  out$interaction_p_value_bh <- NA_real_
  out$interaction_p_value_bh[tested] <- stats::p.adjust(out$interaction_p_value[tested], method = "BH")
  out$interaction_significance_fdr <- p_value_stars(out$interaction_p_value_bh)
  out$human_cluster_order <- match(out$human_cluster, human_levels)
  out$mouse_category_order <- match(out$mouse_category, mouse_levels)
  out[order(out$mouse_category_order, out$human_cluster_order), , drop = FALSE]
}

write_plot_pair <- function(plot, stem, out_dir, width, height) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
}

make_heatmap <- function(x, title, subtitle, facet = TRUE) {
  x$human_cluster <- factor(x$human_cluster, levels = cluster_order(x$human_cluster))
  mouse_levels <- unique(x$mouse_category[order(x$mouse_category_order)])
  x$mouse_category <- factor(x$mouse_category, levels = rev(mouse_levels))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = human_cluster, y = mouse_category, fill = dose_adjusted_partial_spearman_rho)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.22) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(significance_fdr == "ns", "", significance_fdr)), size = 2.6, fontface = "bold") +
    ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "#F4F4F2", high = "#E67E22", midpoint = 0, limits = c(-1, 1), na.value = "#D9DDE3") +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "Human tumor cluster", y = "Mouse host category", fill = "Dose-adjusted\npartial Spearman rho",
      caption = "Stars denote BH-adjusted permutation p-values; blank cells without stars may be non-significant or not testable."
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "right"
    )
  if (facet) p <- p + ggplot2::facet_grid(~Ploidy)
  p
}

script_dir <- resolve_script_dir()
script_path <- file.path(script_dir, "05r_human_cluster_mouse_celltype_association.R")
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "ggplot2", "digest", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05p_latest_human_mouse_sample_matrix", "latest_human_mouse_sample_matrix.rds")
input_completion <- file.path(results_root, "05p_latest_human_mouse_sample_matrix", "completion.txt")
out_dir <- file.path(results_root, "05r_human_cluster_mouse_celltype_association")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds) || !file.exists(input_completion)) stop("05p current sample matrix is incomplete.", call. = FALSE)

dat <- readRDS(input_rds)
sample_metadata <- dat$sample_metadata[order(dat$sample_metadata$sample_order), , drop = FALSE]
human_counts <- add_clr(dat$human_cluster_counts, "human_cluster", "human_clr")
primary_counts <- add_clr(dat$mouse_cell_type_counts, "mouse_cell_type_latest", "mouse_clr")
fine_counts <- add_clr(dat$mouse_fine_category_counts, "mouse_cell_type_state_latest", "mouse_clr")
primary_counts$log2_mouse_type_per_human_tumor <- log2((primary_counts$n_cells + 0.5) / (primary_counts$sample_total_human_tumor_cells + 0.5))
fine_counts$log2_mouse_type_per_human_tumor <- log2((fine_counts$n_cells + 0.5) / (fine_counts$sample_total_human_tumor_cells + 0.5))

association_primary_composition <- run_association_family(
  human_counts, primary_counts, "mouse_cell_type_latest", "mouse_clr",
  "within_mouse_composition_clr", "primary_cell_type", sample_metadata
)
association_primary_recruitment <- run_association_family(
  human_counts, primary_counts, "mouse_cell_type_latest", "log2_mouse_type_per_human_tumor",
  "mouse_type_per_human_tumor_log2", "primary_cell_type", sample_metadata
)
association_fine_composition <- run_association_family(
  human_counts, fine_counts, "mouse_cell_type_state_latest", "mouse_clr",
  "within_mouse_composition_clr", "fine_cell_type_state", sample_metadata
)
association_fine_recruitment <- run_association_family(
  human_counts, fine_counts, "mouse_cell_type_state_latest", "log2_mouse_type_per_human_tumor",
  "mouse_type_per_human_tumor_log2", "fine_cell_type_state", sample_metadata
)
associations <- rbind(
  association_primary_composition, association_primary_recruitment,
  association_fine_composition, association_fine_recruitment
)

interaction_primary_composition <- run_interaction_family(
  human_counts, primary_counts, "mouse_cell_type_latest", "mouse_clr",
  "within_mouse_composition_clr", "primary_cell_type", association_primary_composition, sample_metadata
)
interaction_primary_recruitment <- run_interaction_family(
  human_counts, primary_counts, "mouse_cell_type_latest", "log2_mouse_type_per_human_tumor",
  "mouse_type_per_human_tumor_log2", "primary_cell_type", association_primary_recruitment, sample_metadata
)
interaction_fine_composition <- run_interaction_family(
  human_counts, fine_counts, "mouse_cell_type_state_latest", "mouse_clr",
  "within_mouse_composition_clr", "fine_cell_type_state", association_fine_composition, sample_metadata
)
interaction_fine_recruitment <- run_interaction_family(
  human_counts, fine_counts, "mouse_cell_type_state_latest", "log2_mouse_type_per_human_tumor",
  "mouse_type_per_human_tumor_log2", "fine_cell_type_state", association_fine_recruitment, sample_metadata
)
interactions <- rbind(
  interaction_primary_composition, interaction_primary_recruitment,
  interaction_fine_composition, interaction_fine_recruitment
)

heat_primary_composition <- make_heatmap(
  association_primary_composition,
  "Human tumor-cluster versus mouse cell-type composition associations",
  "Separate 2N and 4N panels; eight biological samples per panel; association is adjusted for Dose"
)
heat_primary_recruitment <- make_heatmap(
  association_primary_recruitment,
  "Human tumor-cluster versus mouse cell-type recruitment-proxy associations",
  "Outcome is log2(mouse type cells / all retained human tumor cells); separate 2N and 4N panels"
)
heat_fine_composition <- make_heatmap(
  association_fine_composition,
  "Human tumor-cluster versus fine mouse cell-type/state associations",
  "Exploratory fine-resolution analysis; separate 2N and 4N panels; association is adjusted for Dose"
)

difference_plot_data <- interaction_primary_composition
difference_plot_data$human_cluster <- factor(difference_plot_data$human_cluster, levels = cluster_order(difference_plot_data$human_cluster))
mouse_levels <- unique(difference_plot_data$mouse_category[order(difference_plot_data$mouse_category_order)])
difference_plot_data$mouse_category <- factor(difference_plot_data$mouse_category, levels = rev(mouse_levels))
heat_difference <- ggplot2::ggplot(difference_plot_data, ggplot2::aes(x = human_cluster, y = mouse_category, fill = difference_rho_4N_minus_2N)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.22) +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(interaction_significance_fdr == "ns", "", interaction_significance_fdr)), size = 2.6, fontface = "bold") +
  ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "#F4F4F2", high = "#E67E22", midpoint = 0, limits = c(-2, 2), na.value = "#D9DDE3", oob = scales::squish) +
  ggplot2::labs(
    title = "Difference in human-cluster-mouse-cell-type association between 4N and 2N",
    subtitle = "Tiles show partial Spearman rho(4N) - rho(2N); interaction tests are exploratory because n=16",
    x = "Human tumor cluster", y = "Mouse host cell type", fill = "Delta rho\n4N - 2N",
    caption = "Stars denote BH-adjusted Human-cluster x Ploidy interaction p-values from Dose-adjusted models."
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

top_pairs <- do.call(rbind, lapply(c("2N", "4N"), function(ploidy_i) {
  z <- association_primary_composition[association_primary_composition$Ploidy == ploidy_i & association_primary_composition$test_status == "tested", , drop = FALSE]
  z <- z[order(is.na(z$p_value_bh), z$p_value_bh, -abs(z$dose_adjusted_partial_spearman_rho)), , drop = FALSE]
  utils::head(z, 6L)
}))
top_pair_samples <- list()
k <- 0L
for (i in seq_len(nrow(top_pairs))) {
  pair <- top_pairs[i, , drop = FALSE]
  h <- human_counts[human_counts$human_cluster == pair$human_cluster & as.character(human_counts$Ploidy) == pair$Ploidy, , drop = FALSE]
  m <- primary_counts[primary_counts$mouse_cell_type_latest == pair$mouse_category & as.character(primary_counts$Ploidy) == pair$Ploidy, , drop = FALSE]
  z <- merge(
    h[c("sample", "Ploidy", "Dose", "human_cluster", "n_cells", "fraction_within_human_tumor", "human_clr")],
    m[c("sample", "mouse_cell_type_latest", "n_cells", "fraction_within_mouse", "mouse_clr")],
    by = "sample", suffixes = c("_human", "_mouse"), all = FALSE, sort = FALSE
  )
  z$pair_label <- paste0(
    pair$Ploidy, ": H", pair$human_cluster, " vs ", pair$mouse_category,
    "\nrho=", format(round(pair$dose_adjusted_partial_spearman_rho, 2), nsmall = 2),
    "; FDR=", format(round(pair$p_value_bh, 3), nsmall = 3)
  )
  z$rho <- pair$dose_adjusted_partial_spearman_rho
  z$fdr <- pair$p_value_bh
  k <- k + 1L
  top_pair_samples[[k]] <- z
}
top_pair_samples <- do.call(rbind, top_pair_samples)
top_pair_samples$Ploidy <- factor(as.character(top_pair_samples$Ploidy), levels = c("2N", "4N"))
top_pair_samples$Dose <- factor(as.character(top_pair_samples$Dose), levels = as.character(unlist(cfg$immune_composition$dose_order, use.names = FALSE)))
ploidy_colors <- c(`2N` = "#3B6FB6", `4N` = "#E67E22")
p_top_scatter <- ggplot2::ggplot(
  top_pair_samples,
  ggplot2::aes(x = fraction_within_human_tumor, y = fraction_within_mouse, color = Ploidy, shape = Dose)
) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      x = fraction_within_human_tumor, y = fraction_within_mouse,
      color = Ploidy, group = 1
    ),
    inherit.aes = FALSE, method = "lm", se = FALSE, linewidth = 0.6,
    linetype = "dashed", show.legend = FALSE
  ) +
  ggplot2::geom_point(size = 2.4, alpha = 0.9) +
  ggplot2::geom_text(ggplot2::aes(label = sample), size = 2.2, vjust = -0.75, check_overlap = TRUE, show.legend = FALSE) +
  ggplot2::facet_wrap(~pair_label, scales = "free", ncol = 4) +
  ggplot2::scale_color_manual(values = ploidy_colors, drop = FALSE) +
  ggplot2::scale_shape_manual(values = c(`0` = 21, `30` = 22, `120` = 24), drop = FALSE) +
  ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  ggplot2::labs(
    title = "Top human tumor-cluster-mouse cell-type associations",
    subtitle = "Top six tested associations per Ploidy by permutation FDR and effect size; each point is one biological sample",
    x = "Human cluster fraction within retained human tumor cells",
    y = "Mouse cell-type fraction within retained mouse cells",
    color = "Ploidy", shape = "Dose (mg/kg)",
    caption = "Dashed lines are descriptive linear trends; inference uses Dose-adjusted partial Spearman correlations and within-Dose exact permutations."
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"), legend.position = "top", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

utils::write.csv(human_counts, file.path(out_dir, "sample_human_cluster_composition_clr.csv"), row.names = FALSE)
utils::write.csv(primary_counts, file.path(out_dir, "sample_mouse_cell_type_composition_and_recruitment_proxy.csv"), row.names = FALSE)
utils::write.csv(fine_counts, file.path(out_dir, "sample_mouse_fine_category_composition_and_recruitment_proxy.csv"), row.names = FALSE)
utils::write.csv(associations, file.path(out_dir, "human_cluster_mouse_category_associations.csv"), row.names = FALSE, na = "NA")
utils::write.csv(interactions, file.path(out_dir, "human_cluster_mouse_category_ploidy_interactions.csv"), row.names = FALSE, na = "NA")
utils::write.csv(top_pairs, file.path(out_dir, "top_primary_composition_association_pairs.csv"), row.names = FALSE, na = "NA")
utils::write.csv(top_pair_samples, file.path(out_dir, "top_primary_composition_association_sample_values.csv"), row.names = FALSE, na = "NA")

human_width <- max(12, length(cluster_order(human_counts$human_cluster)) * 0.55 + 7)
write_plot_pair(heat_primary_composition, "heatmap_human_cluster_mouse_cell_type_composition_association", out_dir, human_width, 8)
write_plot_pair(heat_primary_recruitment, "heatmap_human_cluster_mouse_cell_type_recruitment_proxy_association", out_dir, human_width, 8)
write_plot_pair(heat_fine_composition, "heatmap_human_cluster_mouse_fine_category_association", out_dir, human_width, 13)
write_plot_pair(heat_difference, "heatmap_primary_association_difference_4N_minus_2N", out_dir, human_width, 7)
write_plot_pair(p_top_scatter, "top_human_cluster_mouse_cell_type_association_scatterplots", out_dir, 18, 12)

contract <- list(
  question = "Are human tumor-cluster compositions associated with mouse host cell-type compositions or recruitment proxies, and do associations differ between 2N and 4N tumors?",
  biological_unit = "sample",
  human_predictor = "0.5-cell-pseudocount CLR of each human tumor cluster within all retained human tumor cells in the same sample",
  primary_mouse_outcome = "0.5-cell-pseudocount CLR of each mouse category within all retained mouse cells in the same sample",
  sensitivity_mouse_outcome = "log2((mouse category cells + 0.5) / (all retained human tumor cells + 0.5))",
  separate_ploidy_analysis = "Dose-adjusted partial Spearman correlation in eight samples per Ploidy",
  exact_significance_test = "all 96 permutations of outcome values within Dose blocks for each Ploidy",
  association_multiple_testing = "Benjamini-Hochberg across all human-cluster-by-mouse-category pairs separately within each Ploidy, resolution, and outcome family",
  ploidy_difference = "exploratory Human-cluster CLR by Ploidy interaction in a Dose-adjusted linear model; n=16",
  testability_rule = "at least four nonzero samples for both categories within a Ploidy and at least three distinct transformed values",
  joint_expression_integration_or_clustering = FALSE,
  causal_interpretation = FALSE,
  required_language = "association consistent with differential recruitment; not proof of recruitment",
  captured_fraction_caveat = "scRNA-seq captured-cell proportions are affected by dissociation and recovery"
)
yaml::write_yaml(contract, file.path(out_dir, "analysis_contract.yaml"))

input_files <- c(input_rds, input_completion, normalizePath(config_path, mustWork = TRUE), script_path)
input_manifest <- data.frame(
  role = c("05p_current_sample_matrix", "05p_completion", "analysis_config", "analysis_script"),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05r.txt"))
writeLines(c(
  "05r human tumor-cluster versus mouse host-category association analysis completed.",
  paste0("Human tumor clusters: ", length(cluster_order(human_counts$human_cluster)), "."),
  paste0("Primary mouse cell types: ", length(unique(primary_counts$mouse_cell_type_latest)), "; fine categories: ", length(unique(fine_counts$mouse_cell_type_state_latest)), "."),
  "Separate 2N and 4N associations use eight biological samples and 96 exact within-Dose permutations per Ploidy.",
  paste0("Primary-composition associations with BH-FDR<0.05: ", sum(association_primary_composition$p_value_bh < 0.05, na.rm = TRUE), "."),
  paste0("Primary recruitment-proxy associations with BH-FDR<0.05: ", sum(association_primary_recruitment$p_value_bh < 0.05, na.rm = TRUE), "."),
  paste0("Exploratory primary-composition Ploidy interactions with BH-FDR<0.05: ", sum(interaction_primary_composition$interaction_p_value_bh < 0.05, na.rm = TRUE), "."),
  "Human and mouse matrices were not jointly integrated or reclustered.",
  "Results are associations consistent with differential recruitment hypotheses, not causal proof of recruitment."
), file.path(out_dir, "completion.txt"))

output_files <- setdiff(list.files(out_dir, full.names = TRUE), file.path(out_dir, "output_manifest_sha256.csv"))
output_manifest <- data.frame(
  file = basename(output_files), size_bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, file.path(out_dir, "output_manifest_sha256.csv"), row.names = FALSE)
message("[05r] Completed. Output: ", out_dir)
