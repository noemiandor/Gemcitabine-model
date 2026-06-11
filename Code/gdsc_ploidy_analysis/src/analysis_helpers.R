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
