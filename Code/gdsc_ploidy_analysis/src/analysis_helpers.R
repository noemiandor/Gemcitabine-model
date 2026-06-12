run_enrichment_or_stop <- function(values,
                                   annotations,
                                   cancer,
                                   direction,
                                   permute_n,
                                   pvalue_cutoff,
                                   enrichment_fn = enrichment) {
  tryCatch(
    enrichment_fn(
      values,
      annotations,
      permute.n = permute_n,
      normalize = FALSE,
      pvalue.cutoff = pvalue_cutoff
    )$pvalue,
    error = function(e) {
      stop(
        sprintf(
          "Enrichment failed for cancer=%s direction=%s drugs=%d groups=%d: %s",
          cancer,
          direction,
          length(values),
          length(unique(annotations$group)),
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
}

plot_barplot_or_stop <- function(plot_values,
                                 colors,
                                 cancer,
                                 xlab,
                                 plotting_fn = barplot) {
  tryCatch(
    plotting_fn(
      plot_values,
      col = colors,
      main = cancer,
      horiz = TRUE,
      las = 2,
      cex.lab = 0.7,
      cex.names = 0.35,
      xlab = xlab
    ),
    error = function(e) {
      stop(
        sprintf(
          "Plotting failed for cancer=%s plotted_drugs=%d: %s",
          cancer,
          length(plot_values),
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
}

normalize_cell_line_name <- function(x) {
  x <- toupper(as.character(x))
  gsub("[^A-Z0-9]", "", x)
}

audit_normalized_key_collisions <- function(dt, raw_col, key_col, source_name) {
  missing_cols <- setdiff(c(raw_col, key_col), colnames(dt))
  if (length(missing_cols) > 0) {
    stop("Collision audit input is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  tmp <- unique(as.data.frame(dt)[, c(raw_col, key_col), drop = FALSE])
  names(tmp) <- c("raw_name", "normalized_key")
  tmp <- tmp[!is.na(tmp$normalized_key) & nzchar(tmp$normalized_key), , drop = FALSE]

  out <- data.table::as.data.table(tmp)
  out <- out[, .(
    n_raw_names = data.table::uniqueN(raw_name),
    raw_names = paste(sort(unique(raw_name)), collapse = " | ")
  ), by = normalized_key]
  out <- out[n_raw_names > 1]
  out[, source := source_name]
  out[]
}

validate_ploidy_key_collisions <- function(appCL,
                                           key_col = "CELL_LINE_KEY",
                                           raw_col = "Cell iname",
                                           ploidy_col = "ploidy",
                                           tolerance = 0.05,
                                           tables_dir = NULL) {
  missing_cols <- setdiff(c(key_col, raw_col, ploidy_col), colnames(appCL))
  if (length(missing_cols) > 0) {
    stop("Ploidy collision validation is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  appCL_dt <- data.table::as.data.table(appCL)
  summary_dt <- appCL_dt[, .(
    n_rows = .N,
    n_raw_names = data.table::uniqueN(get(raw_col)),
    ploidy_min = min(get(ploidy_col), na.rm = TRUE),
    ploidy_max = max(get(ploidy_col), na.rm = TRUE),
    ploidy_range = max(get(ploidy_col), na.rm = TRUE) - min(get(ploidy_col), na.rm = TRUE),
    raw_names = paste(sort(unique(get(raw_col))), collapse = " | ")
  ), by = key_col]
  data.table::setnames(summary_dt, key_col, "normalized_key")

  collapsed <- summary_dt[n_raw_names > 1 & ploidy_range <= tolerance]
  failures <- summary_dt[n_raw_names > 1 & ploidy_range > tolerance]

  if (!is.null(tables_dir)) {
    write_tsv(as.data.frame(summary_dt), file.path(tables_dir, "cell_line_key_ploidy_summary.tsv"))
    write_tsv(as.data.frame(collapsed), file.path(tables_dir, "cell_line_key_ploidy_near_duplicate_collapses.tsv"))
    if (nrow(failures) > 0) {
      write_tsv(as.data.frame(failures), file.path(tables_dir, "cell_line_key_collisions_ploidy_FAIL.tsv"))
    }
  }

  if (nrow(failures) > 0) {
    stop(
      "Normalized ploidy cell-line keys are not unique and have discordant ploidy values. Review cell_line_key_collisions_ploidy_FAIL.tsv.",
      call. = FALSE
    )
  }

  invisible(list(summary = summary_dt, collapsed = collapsed, failures = failures))
}

write_cell_line_matching_delta <- function(dr, appCL, path) {
  raw_names <- if ("CELL_LINE_NAME_RAW" %in% colnames(dr)) dr$CELL_LINE_NAME_RAW else dr$CELL_LINE_NAME
  gdsc <- unique(data.frame(
    gdsc_cell_line_name = dr$CELL_LINE_NAME,
    gdsc_cell_line_name_raw = raw_names,
    normalized_key = dr$CELL_LINE_KEY,
    stringsAsFactors = FALSE
  ))

  ploidy <- unique(data.frame(
    ploidy_cell_line_name = appCL$`Cell iname`,
    normalized_key = appCL$CELL_LINE_KEY,
    stringsAsFactors = FALSE
  ))

  ploidy_by_key <- data.table::as.data.table(ploidy)
  ploidy_by_key <- ploidy_by_key[, .(
    normalized_match_count = .N,
    normalized_match_raw_names = paste(sort(unique(ploidy_cell_line_name)), collapse = " | ")
  ), by = normalized_key]

  delta <- data.table::as.data.table(gdsc)
  delta[, canonical_raw_match := gdsc_cell_line_name %in% appCL$`Cell iname`]
  delta <- merge(delta, ploidy_by_key, by = "normalized_key", all.x = TRUE, sort = FALSE)
  delta[, normalized_key_match := !is.na(normalized_match_count)]
  delta[is.na(normalized_match_count), normalized_match_count := 0L]
  delta[is.na(normalized_match_raw_names), normalized_match_raw_names := ""]
  delta[, match_status := data.table::fifelse(
    canonical_raw_match & normalized_key_match,
    "matched_by_both",
    data.table::fifelse(
      canonical_raw_match & !normalized_key_match,
      "raw_only",
      data.table::fifelse(!canonical_raw_match & normalized_key_match, "normalized_only", "unmatched")
    )
  )]
  delta <- delta[order(match_status, gdsc_cell_line_name)]
  write_tsv(as.data.frame(delta), path)
  invisible(delta)
}

compute_drug_ploidy_correlation <- function(response, ploidy, min_n = 10) {
  ok <- is.finite(response) & is.finite(ploidy)
  n <- sum(ok)

  empty <- function(method) {
    data.frame(
      n = n,
      pearson_r = NA_real_,
      pearson_p = NA_real_,
      pearson_ci_low = NA_real_,
      pearson_ci_high = NA_real_,
      spearman_rho = NA_real_,
      spearman_p = NA_real_,
      spearman_ci_low = NA_real_,
      spearman_ci_high = NA_real_,
      spearman_ci_method = method,
      stringsAsFactors = FALSE
    )
  }

  if (n < min_n) {
    return(empty("not_estimated_n_below_minimum"))
  }

  response_ok <- response[ok]
  ploidy_ok <- ploidy[ok]
  if (stats::sd(response_ok) == 0 || stats::sd(ploidy_ok) == 0) {
    return(empty("not_estimated_no_variation"))
  }

  pearson <- tryCatch(
    suppressWarnings(stats::cor.test(response_ok, ploidy_ok, method = "pearson")),
    error = function(e) NULL
  )
  spearman <- tryCatch(
    suppressWarnings(stats::cor.test(response_ok, ploidy_ok, method = "spearman", exact = FALSE)),
    error = function(e) NULL
  )

  pearson_ci <- c(NA_real_, NA_real_)
  pearson_r <- pearson_p <- NA_real_
  if (!is.null(pearson)) {
    pearson_r <- unname(pearson$estimate)
    pearson_p <- pearson$p.value
    pearson_ci <- pearson$conf.int
  }

  spearman_rho <- spearman_p <- NA_real_
  spearman_ci <- c(NA_real_, NA_real_)
  spearman_ci_method <- "not_estimated"
  if (!is.null(spearman)) {
    spearman_rho <- unname(spearman$estimate)
    spearman_p <- spearman$p.value
    if (is.finite(spearman_rho) && abs(spearman_rho) < 1 && n > 3) {
      z <- atanh(spearman_rho)
      se <- 1 / sqrt(n - 3)
      spearman_ci <- tanh(c(z - 1.96 * se, z + 1.96 * se))
      spearman_ci_method <- "approximate_fisher_transform"
    }
  }

  data.frame(
    n = n,
    pearson_r = pearson_r,
    pearson_p = pearson_p,
    pearson_ci_low = pearson_ci[1],
    pearson_ci_high = pearson_ci[2],
    spearman_rho = spearman_rho,
    spearman_p = spearman_p,
    spearman_ci_low = spearman_ci[1],
    spearman_ci_high = spearman_ci[2],
    spearman_ci_method = spearman_ci_method,
    stringsAsFactors = FALSE
  )
}

response_direction_note <- function(metric) {
  switch(
    metric,
    Z_SCORE = "Positive r means higher ploidy is associated with higher GDSC Z-score; sensitivity direction requires metric-specific interpretation.",
    LN_IC50 = "Positive r means higher ploidy is associated with higher LN_IC50; interpret sensitivity direction separately from the correlation sign.",
    AUC = "Positive r means higher ploidy is associated with higher AUC; interpret sensitivity direction separately from the correlation sign.",
    "Positive r means higher ploidy is associated with a higher response metric value; interpret sensitivity direction separately."
  )
}

ploidy_map_for_matching <- function(appCL, matching_mode) {
  matching_mode <- match.arg(matching_mode, c("raw", "normalized"))
  if (matching_mode == "raw") {
    out <- unique(data.frame(
      match_key = appCL$`Cell iname`,
      ploidy = appCL$ploidy,
      stringsAsFactors = FALSE
    ))
  } else {
    appCL_dt <- data.table::as.data.table(appCL)
    out_dt <- appCL_dt[, .(
      ploidy = mean(ploidy, na.rm = TRUE),
      ploidy_rows = .N,
      ploidy_raw_names = paste(sort(unique(`Cell iname`)), collapse = " | ")
    ), by = CELL_LINE_KEY]
    data.table::setnames(out_dt, "CELL_LINE_KEY", "match_key")
    out <- as.data.frame(out_dt)
  }
  out[!is.na(out$match_key) & nzchar(out$match_key), , drop = FALSE]
}

build_drug_ploidy_correlation_table <- function(dr,
                                                appCL,
                                                metrics,
                                                matching_mode = c("raw", "normalized"),
                                                min_n = 10) {
  matching_mode <- match.arg(matching_mode)
  metrics <- intersect(metrics, colnames(dr))
  cancer_types <- c("allcancers", sort(unique(dr$TCGA_DESC)))
  ploidy_map <- ploidy_map_for_matching(appCL, matching_mode)
  match_col <- if (matching_mode == "raw") "CELL_LINE_NAME" else "CELL_LINE_KEY"
  rows <- list()
  idx <- 0L

  for (cancer_type in cancer_types) {
    dr_sub <- dr
    if (cancer_type != "allcancers") {
      dr_sub <- dr[dr$TCGA_DESC == cancer_type, , drop = FALSE]
    }
    if (nrow(dr_sub) == 0) {
      next
    }

    for (drug in sort(unique(dr_sub$DRUG_NAME))) {
      dr_drug <- dr_sub[dr_sub$DRUG_NAME == drug, , drop = FALSE]
      merged <- merge(
        dr_drug,
        ploidy_map,
        by.x = match_col,
        by.y = "match_key",
        all = FALSE,
        sort = FALSE
      )
      for (metric in metrics) {
        stats <- compute_drug_ploidy_correlation(merged[[metric]], merged$ploidy, min_n = min_n)
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          cancer_type = cancer_type,
          drug = drug,
          metric = metric,
          stats,
          response_direction_note = response_direction_note(metric),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }
  res <- data.table::as.data.table(do.call(rbind, rows))
  res[, pearson_fdr := stats::p.adjust(pearson_p, method = "BH"), by = .(cancer_type, metric)]
  res[, spearman_fdr := stats::p.adjust(spearman_p, method = "BH"), by = .(cancer_type, metric)]
  res[, rank_by_pearson_desc := data.table::frank(-pearson_r, ties.method = "min", na.last = "keep"), by = .(cancer_type, metric)]
  res[, rank_by_pearson_asc := data.table::frank(pearson_r, ties.method = "min", na.last = "keep"), by = .(cancer_type, metric)]
  res[, rank_by_spearman_desc := data.table::frank(-spearman_rho, ties.method = "min", na.last = "keep"), by = .(cancer_type, metric)]
  res[, rank_by_spearman_asc := data.table::frank(spearman_rho, ties.method = "min", na.last = "keep"), by = .(cancer_type, metric)]

  ordered_cols <- c(
    "cancer_type",
    "drug",
    "metric",
    "n",
    "pearson_r",
    "pearson_p",
    "pearson_ci_low",
    "pearson_ci_high",
    "pearson_fdr",
    "spearman_rho",
    "spearman_p",
    "spearman_ci_low",
    "spearman_ci_high",
    "spearman_ci_method",
    "spearman_fdr",
    "rank_by_pearson_desc",
    "rank_by_pearson_asc",
    "rank_by_spearman_desc",
    "rank_by_spearman_asc",
    "response_direction_note"
  )
  as.data.frame(res[, ..ordered_cols])
}

safe_sheet_name <- function(x) {
  x <- gsub("[\\[\\]\\:\\*\\?\\/\\\\]", "_", as.character(x))
  substr(x, 1, 31)
}

write_correlations_xlsx <- function(cor_dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  wb <- openxlsx::createWorkbook()
  if (nrow(cor_dt) == 0) {
    openxlsx::addWorksheet(wb, "no_results")
    openxlsx::writeData(wb, "no_results", data.frame(message = "No correlation results"))
  } else {
    for (cancer_type in sort(unique(cor_dt$cancer_type))) {
      sheet <- safe_sheet_name(cancer_type)
      openxlsx::addWorksheet(wb, sheet)
      sheet_dt <- cor_dt[cor_dt$cancer_type == cancer_type, , drop = FALSE]
      openxlsx::writeData(wb, sheet, sheet_dt)
      openxlsx::freezePane(wb, sheet, firstRow = TRUE)
      openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(sheet_dt)), widths = "auto")
    }
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

write_correlation_delta <- function(raw_dt, normalized_dt, path) {
  raw <- raw_dt
  normalized <- normalized_dt
  suffix_cols <- c(
    "n",
    "pearson_r",
    "pearson_p",
    "pearson_fdr",
    "spearman_rho",
    "spearman_p",
    "spearman_fdr"
  )
  raw <- raw[, c("cancer_type", "drug", "metric", suffix_cols), drop = FALSE]
  normalized <- normalized[, c("cancer_type", "drug", "metric", suffix_cols), drop = FALSE]
  names(raw)[match(suffix_cols, names(raw))] <- paste0(suffix_cols, "_raw")
  names(normalized)[match(suffix_cols, names(normalized))] <- paste0(suffix_cols, "_normalized")

  delta <- merge(raw, normalized, by = c("cancer_type", "drug", "metric"), all = TRUE, sort = FALSE)
  delta$n_delta <- delta$n_normalized - delta$n_raw
  delta$pearson_r_delta <- delta$pearson_r_normalized - delta$pearson_r_raw
  delta$spearman_rho_delta <- delta$spearman_rho_normalized - delta$spearman_rho_raw
  delta$delta_status <- ifelse(
    is.na(delta$n_raw),
    "normalized_only",
    ifelse(
      is.na(delta$n_normalized),
      "raw_only",
      ifelse(
        delta$n_delta == 0 &
          (is.na(delta$pearson_r_delta) | abs(delta$pearson_r_delta) < .Machine$double.eps^0.5) &
          (is.na(delta$spearman_rho_delta) | abs(delta$spearman_rho_delta) < .Machine$double.eps^0.5),
        "unchanged",
        "changed"
      )
    )
  )
  write_tsv(delta, path)
  invisible(delta)
}

normalize_drug_key <- function(x) {
  x <- toupper(as.character(x))
  gsub("[^A-Z0-9]", "", x)
}

write_gemcitabine_rank_summary <- function(cor_dt,
                                           path,
                                           aliases = c("GEMCITABINE", "GEMZAR")) {
  aliases <- normalize_drug_key(aliases)
  out <- cor_dt[normalize_drug_key(cor_dt$drug) %in% aliases, , drop = FALSE]
  write_tsv(out, path)
  invisible(out)
}

build_tissue_adjusted_ploidy_models <- function(dr,
                                                appCL,
                                                metrics,
                                                min_n = 20L,
                                                min_tissues = 3L,
                                                min_rows_per_tissue = 2L) {
  metrics <- intersect(metrics, colnames(dr))
  ploidy_map <- unique(data.frame(
    CELL_LINE_NAME = appCL$`Cell iname`,
    ploidy = appCL$ploidy,
    stringsAsFactors = FALSE
  ))
  merged_all <- merge(dr, ploidy_map, by = "CELL_LINE_NAME", all = FALSE, sort = FALSE)
  rows <- list()
  idx <- 0L

  empty_row <- function(drug,
                        metric,
                        data,
                        status,
                        skipped_reason = "",
                        rank_deficient = NA,
                        beta = NA_real_,
                        se = NA_real_,
                        tval = NA_real_,
                        pval = NA_real_,
                        notes = "") {
    data.frame(
      drug = drug,
      metric = metric,
      n = nrow(data),
      n_cancer_types = length(unique(data$TCGA_DESC)),
      n_tissues = length(unique(data$TCGA_DESC)),
      beta_ploidy = beta,
      se_ploidy = se,
      t_ploidy = tval,
      p_ploidy = pval,
      fdr_ploidy = NA_real_,
      model_formula = "response ~ ploidy + TCGA_DESC",
      model_status = status,
      rank_deficient = rank_deficient,
      skipped_reason = skipped_reason,
      notes = notes,
      stringsAsFactors = FALSE
    )
  }

  for (drug in sort(unique(merged_all$DRUG_NAME))) {
    drug_dt <- merged_all[merged_all$DRUG_NAME == drug, , drop = FALSE]
    for (metric in metrics) {
      fit_data <- data.frame(
        response = drug_dt[[metric]],
        ploidy = drug_dt$ploidy,
        TCGA_DESC = drug_dt$TCGA_DESC,
        stringsAsFactors = FALSE
      )
      fit_data <- fit_data[
        is.finite(fit_data$response) &
          is.finite(fit_data$ploidy) &
          !is.na(fit_data$TCGA_DESC) &
          nzchar(fit_data$TCGA_DESC),
        ,
        drop = FALSE
      ]
      fit_data$TCGA_DESC <- droplevels(factor(fit_data$TCGA_DESC))

      idx <- idx + 1L
      if (nrow(fit_data) < min_n) {
        rows[[idx]] <- empty_row(drug, metric, fit_data, "skipped_insufficient_n", "n below tissue_model_min_n")
        next
      }
      tissue_counts <- table(fit_data$TCGA_DESC)
      if (length(tissue_counts) < min_tissues) {
        rows[[idx]] <- empty_row(drug, metric, fit_data, "skipped_insufficient_tissue_count", "fewer tissues than tissue_model_min_tissues")
        next
      }
      if (any(tissue_counts < min_rows_per_tissue)) {
        rows[[idx]] <- empty_row(drug, metric, fit_data, "skipped_sparse_tissue_levels", "one or more tissues below tissue_model_min_rows_per_tissue")
        next
      }

      fit <- tryCatch(
        stats::lm(response ~ ploidy + TCGA_DESC, data = fit_data),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        rows[[idx]] <- empty_row(drug, metric, fit_data, "fit_error", conditionMessage(fit))
        next
      }

      coef_table <- summary(fit)$coefficients
      rank_deficient <- fit$rank < length(stats::coef(fit))
      if (!"ploidy" %in% rownames(coef_table) || is.na(coef_table["ploidy", "Estimate"])) {
        rows[[idx]] <- empty_row(
          drug,
          metric,
          fit_data,
          "rank_deficient_no_ploidy_coef",
          "ploidy coefficient absent or NA",
          rank_deficient = TRUE
        )
        next
      }

      status <- if (rank_deficient) "rank_deficient_ploidy_estimable" else "ok"
      rows[[idx]] <- empty_row(
        drug,
        metric,
        fit_data,
        status,
        rank_deficient = rank_deficient,
        beta = coef_table["ploidy", "Estimate"],
        se = coef_table["ploidy", "Std. Error"],
        tval = coef_table["ploidy", "t value"],
        pval = coef_table["ploidy", "Pr(>|t|)"],
        notes = if (rank_deficient) "Model is rank-deficient, but ploidy coefficient is estimable." else ""
      )
    }
  }

  out <- do.call(rbind, rows)
  out_dt <- data.table::as.data.table(out)
  out_dt[, fdr_ploidy := stats::p.adjust(p_ploidy, method = "BH"), by = metric]
  as.data.frame(out_dt)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

validate_required_groups <- function(groups, required_groups = c("SIGNALING", "CYTOTOXIC")) {
  missing_groups <- setdiff(required_groups, unique(groups))
  if (length(missing_groups) > 0) {
    stop(
      "Missing required annotation groups after filtering: ",
      paste(missing_groups, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

resolve_duplicate_drug_cell_lines <- function(drug_table,
                                             metric,
                                             strategy = "lowest_rmse",
                                             qc_dir = NULL) {
  required_cols <- c("DRUG_NAME", "CELL_LINE_NAME", metric)
  missing_cols <- setdiff(required_cols, colnames(drug_table))
  if (length(missing_cols) > 0) {
    stop("Drug table is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!strategy %in% c("lowest_rmse", "first_input_order")) {
    stop("Unsupported duplicate resolution strategy: ", strategy, call. = FALSE)
  }

  key_cols <- c("DRUG_NAME", "CELL_LINE_NAME")
  duplicate_mask <- duplicated(drug_table[, key_cols]) | duplicated(drug_table[, key_cols], fromLast = TRUE)
  duplicate_rows <- drug_table[duplicate_mask, , drop = FALSE]

  if (!is.null(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(duplicate_rows, file.path(qc_dir, "duplicate_drug_cell_line_records.tsv"))
  }

  if (strategy == "first_input_order") {
    resolved <- drug_table[!duplicated(drug_table[, key_cols]), , drop = FALSE]
  } else {
    order_cols <- list(drug_table$DRUG_NAME, drug_table$CELL_LINE_NAME)
    if ("RMSE" %in% colnames(drug_table)) {
      order_cols <- c(order_cols, list(is.na(drug_table$RMSE), drug_table$RMSE))
    }
    for (col in c("NLME_RESULT_ID", "NLME_CURVE_ID", "DRUG_ID", "COSMIC_ID")) {
      if (col %in% colnames(drug_table)) {
        order_cols <- c(order_cols, list(drug_table[[col]]))
      }
    }
    ordered <- drug_table[do.call(order, order_cols), , drop = FALSE]
    resolved <- ordered[!duplicated(ordered[, key_cols]), , drop = FALSE]
  }

  remaining_duplicates <- duplicated(resolved[, key_cols]) | duplicated(resolved[, key_cols], fromLast = TRUE)
  if (any(remaining_duplicates)) {
    stop("Duplicate resolution failed to produce one row per DRUG_NAME and CELL_LINE_NAME.", call. = FALSE)
  }

  summary <- data.frame(
    metric = metric,
    strategy = strategy,
    input_rows = nrow(drug_table),
    duplicate_rows = nrow(duplicate_rows),
    duplicate_keys = if (nrow(duplicate_rows) == 0) 0 else nrow(unique(duplicate_rows[, key_cols])),
    resolved_rows = nrow(resolved),
    removed_rows = nrow(drug_table) - nrow(resolved),
    stringsAsFactors = FALSE
  )
  if (!is.null(qc_dir)) {
    write_tsv(summary, file.path(qc_dir, "duplicate_resolution_summary.tsv"))
  }

  attr(resolved, "duplicate_resolution_summary") <- summary
  resolved
}
