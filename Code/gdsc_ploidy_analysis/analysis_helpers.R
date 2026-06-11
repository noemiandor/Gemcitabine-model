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
