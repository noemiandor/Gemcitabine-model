options(stringsAsFactors = FALSE)

arg_value <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(prefix, args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[1], fixed = TRUE)
}

split_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    return(character())
  }
  trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

default_paths <- function(base_dir) {
  list(
    app_cl = file.path(base_dir, "data", "raw", "Cell_app_export.txt"),
    expression_columns = file.path(base_dir, "data", "derived", "ccle_expression_columns.tsv"),
    ploidy = file.path(base_dir, "data", "manual", "breast_ccle_ploidy.tsv"),
    drug_aliases = file.path(base_dir, "data", "raw", "DrugAliases.txt"),
    breast_drug_sensitivity = file.path(base_dir, "data", "raw", "breast_cancer_drug_sensitivity"),
    grbrowser = file.path(base_dir, "data", "raw", "grbrowser")
  )
}

default_output_dir <- function(base_dir) {
  repo_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = TRUE)
  file.path(
    repo_root,
    "Results",
    "public_data",
    "ccle_ploidy_analysis",
    "runs",
    paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_ccle")
  )
}

load_common_inputs <- function(base_dir, require_primary_adherent = FALSE) {
  paths <- default_paths(base_dir)
  required <- unlist(paths[c("app_cl", "expression_columns", "ploidy", "drug_aliases")])
  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop("Missing required input files: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  appCL <- read_tsv(paths$app_cl)
  ploidy <- read_ploidy_table(paths$ploidy)
  expression_columns <- read_expression_columns(paths$expression_columns)
  ploidy <- prepare_ploidy(ploidy, appCL, expression_columns, require_primary_adherent)

  list(paths = paths, appCL = appCL, ploidy = ploidy, expression_columns = expression_columns)
}

matrix_from_named_vectors <- function(named_vectors, rows) {
  mat <- matrix(
    NA_real_,
    nrow = length(rows),
    ncol = length(named_vectors),
    dimnames = list(rows, names(named_vectors))
  )
  for (nm in names(named_vectors)) {
    values <- as.numeric(named_vectors[[nm]])
    names(values) <- names(named_vectors[[nm]])
    overlap <- intersect(rows, names(values))
    mat[overlap, nm] <- values[overlap]
  }
  mat
}

run_legacy_z_score <- function(base_dir, out_dir, correlation_threshold, drug_filter) {
  inputs <- load_common_inputs(base_dir, require_primary_adherent = TRUE)
  paths <- inputs$paths
  ploidy <- inputs$ploidy
  files <- list.files(paths$breast_drug_sensitivity, full.names = TRUE)
  if (length(files) == 0) {
    stop("No legacy breast drug-sensitivity files found in ", paths$breast_drug_sensitivity, call. = FALSE)
  }

  clin <- list()
  for (path in files) {
    sName <- strsplit(tools::file_path_sans_ext(basename(path)), "_")[[1]][1]
    if (!paste0(sName, "_BREAST") %in% names(ploidy)) {
      next
    }
    tmp <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE)
    if (!all(c("Drug Name", "Z Score") %in% colnames(tmp))) {
      next
    }
    z_score <- suppressWarnings(as.numeric(tmp$`Z Score`))
    means <- tapply(z_score, as.character(tmp$`Drug Name`), mean, na.rm = TRUE)
    means <- means[sort(names(means))]
    clin[[sName]] <- means
  }
  if (length(clin) == 0) {
    stop("No legacy breast drug-sensitivity records overlapped with the filtered ploidy table.", call. = FALSE)
  }

  fr <- as.data.frame(table(unlist(lapply(clin, names))), stringsAsFactors = FALSE)
  colnames(fr) <- c("drug", "freq")
  fr <- fr[order(fr$freq, decreasing = TRUE, fr$drug), ]
  coverage_cutoff <- 0.8 * max(fr$freq)
  doi <- as.character(fr$drug[fr$freq >= coverage_cutoff])
  clin_mat <- matrix_from_named_vectors(clin, doi)
  colnames(clin_mat) <- paste0(colnames(clin_mat), "_BREAST")

  correlations <- do.call(rbind, lapply(rownames(clin_mat), function(drug) {
    cbind(drug = drug, compute_pearson(clin_mat[drug, ], ploidy[colnames(clin_mat)], min_n = 3))
  }))
  correlations$estimate <- as.numeric(correlations$estimate)
  correlations$p.value <- as.numeric(correlations$p.value)
  correlations$n <- as.integer(correlations$n)
  correlations$pathway_name <- match_pathways_from_aliases(correlations$drug, paths$drug_aliases)
  correlations$pathway_name <- finalize_pathway_labels(correlations$pathway_name)
  correlations <- correlations[order(correlations$estimate), , drop = FALSE]
  plotted <- correlations[is.finite(correlations$estimate) & abs(correlations$estimate) >= correlation_threshold, , drop = FALSE]

  if (length(drug_filter) > 0) {
    keep <- canonical_drug_name(plotted$drug) %in% canonical_drug_name(drug_filter)
    plotted <- plotted[keep, , drop = FALSE]
  }
  if (nrow(plotted) == 0) {
    stop("No drugs passed the legacy Z_SCORE plotting filters.", call. = FALSE)
  }

  plotted$display_drug <- plotted$drug
  plotted$plot_label <- gsub("Inhibitor", "I.", plotted$display_drug)
  plotted$direction <- ifelse(plotted$estimate > 0, "Low ploidy is sensitive", "High ploidy is sensitive")

  list(
    metric = "Z_SCORE",
    metric_label = "BreastCancerDrugSensitivity Z Score",
    metric_source = "legacy_breast_drug_sensitivity_z_score",
    metric_note = "Historical CCLE breast drug-sensitivity values from the BreastCancerDrugSensitivity Z Score column.",
    response_column = "Z Score",
    lower_metric_more_sensitive = TRUE,
    sensitivity_direction_assumption = "Lower BreastCancerDrugSensitivity Z Score values are documented as more sensitive.",
    ploidy = ploidy,
    response_matrix = clin_mat,
    drug_coverage = fr,
    coverage_cutoff = coverage_cutoff,
    correlations = correlations,
    plotted = plotted,
    input_files = c(paths$app_cl, paths$expression_columns, paths$ploidy, paths$drug_aliases, files)
  )
}

read_grbrowser_file <- function(path, metric) {
  x <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "")
  colnames(x) <- gsub("Small_Molecule", "Perturbagen", colnames(x), fixed = TRUE)
  metric_idx <- match(toupper(metric), toupper(colnames(x)))
  if (!all(c("Cell_Line", "Perturbagen") %in% colnames(x)) || is.na(metric_idx)) {
    return(NULL)
  }
  out <- x[, c("Cell_Line", "Perturbagen", colnames(x)[metric_idx])]
  colnames(out)[3] <- "response"
  out$Cell_Line <- gsub("-", "", out$Cell_Line)
  out$response <- suppressWarnings(as.numeric(out$response))
  out <- out[is.finite(out$response), , drop = FALSE]
  out
}

run_grbrowser_metric <- function(base_dir, out_dir, metric, correlation_threshold, drug_filter) {
  inputs <- load_common_inputs(base_dir, require_primary_adherent = FALSE)
  paths <- inputs$paths
  ploidy <- inputs$ploidy
  files <- list.files(paths$grbrowser, full.names = TRUE, pattern = "[.]tsv$")
  if (length(files) == 0) {
    stop("No grbrowser TSV files found in ", paths$grbrowser, call. = FALSE)
  }

  clin <- list()
  used_files <- character()
  for (path in files) {
    tmp <- read_grbrowser_file(path, metric)
    if (is.null(tmp) || nrow(tmp) == 0) {
      next
    }
    used_files <- c(used_files, path)
    sNames <- intersect(unique(tmp$Cell_Line), gsub("_BREAST", "", names(ploidy), fixed = TRUE))
    for (sName in sNames) {
      ii <- which(tmp$Cell_Line == sName)
      if (sName %in% names(clin)) {
        clin[[sName]] <- rbind(clin[[sName]], tmp[ii, c("Perturbagen", "response")])
      } else {
        clin[[sName]] <- tmp[ii, c("Perturbagen", "response")]
      }
    }
  }
  if (length(clin) == 0) {
    stop("No grbrowser records for metric ", metric, " overlapped with the filtered ploidy table.", call. = FALSE)
  }

  for (sName in names(clin)) {
    tmp <- clin[[sName]]
    clin[[sName]] <- tapply(tmp$response, as.character(tmp$Perturbagen), median, na.rm = TRUE)
  }

  fr <- as.data.frame(table(unlist(lapply(clin, names))), stringsAsFactors = FALSE)
  colnames(fr) <- c("drug", "freq")
  fr <- fr[order(fr$freq, decreasing = TRUE, fr$drug), ]
  coverage_cutoff <- 0.8 * max(fr$freq)
  doi <- as.character(fr$drug[fr$freq >= coverage_cutoff])
  clin_mat <- matrix_from_named_vectors(clin, doi)
  colnames(clin_mat) <- paste0(colnames(clin_mat), "_BREAST")

  correlations <- do.call(rbind, lapply(rownames(clin_mat), function(drug) {
    cbind(drug = drug, compute_pearson(clin_mat[drug, ], ploidy[colnames(clin_mat)], min_n = 3))
  }))
  correlations$estimate <- as.numeric(correlations$estimate)
  correlations$p.value <- as.numeric(correlations$p.value)
  correlations$n <- as.integer(correlations$n)
  correlations$drugName <- correlations$drug

  manual <- manual_grbrowser_annotations(correlations$drug)
  correlations$pathway_name <- manual$pathway_name
  correlations$drugCategory <- manual$drugCategory
  alias_pathways <- match_pathways_from_aliases(correlations$drug, paths$drug_aliases)
  correlations$pathway_name[!is.na(alias_pathways)] <- alias_pathways[!is.na(alias_pathways)]
  correlations$pathway_name <- finalize_pathway_labels(correlations$pathway_name, correlations$drugCategory)
  correlations <- correlations[order(abs(correlations$estimate)), , drop = FALSE]

  display_names <- caps(correlations$drug)
  keep_unique <- !duplicated(display_names, fromLast = TRUE)
  correlations <- correlations[keep_unique, , drop = FALSE]
  correlations$display_drug <- display_names[keep_unique]
  correlations <- correlations[order(correlations$estimate), , drop = FALSE]

  plotted <- correlations[is.finite(correlations$estimate) & abs(correlations$estimate) >= correlation_threshold, , drop = FALSE]
  if (length(drug_filter) > 0) {
    keep <- canonical_drug_name(plotted$display_drug) %in% canonical_drug_name(drug_filter) |
      canonical_drug_name(plotted$drug) %in% canonical_drug_name(drug_filter)
    plotted <- plotted[keep, , drop = FALSE]
  }
  if (nrow(plotted) == 0) {
    stop("No drugs passed the grbrowser plotting filters.", call. = FALSE)
  }

  plotted$plot_label <- gsub("Inhibitor", "I.", plotted$display_drug)
  lower_metric_more_sensitive <- toupper(metric) %in% c("IC50", "EC50", "GR50", "GEC50")
  plotted$direction <- if (lower_metric_more_sensitive) {
    ifelse(plotted$estimate > 0, "Low ploidy is sensitive", "High ploidy is sensitive")
  } else {
    ifelse(plotted$estimate > 0, "High ploidy is sensitive", "Low ploidy is sensitive")
  }

  list(
    metric = metric,
    metric_label = metric,
    metric_source = "grbrowser",
    metric_note = "MEP-LINCS grbrowser metric.",
    response_column = metric,
    lower_metric_more_sensitive = lower_metric_more_sensitive,
    sensitivity_direction_assumption = if (lower_metric_more_sensitive) {
      paste0("Lower ", metric, " values are treated as more sensitive.")
    } else {
      paste0("Higher ", metric, " values are treated as more sensitive.")
    },
    ploidy = ploidy,
    response_matrix = clin_mat,
    drug_coverage = fr,
    coverage_cutoff = coverage_cutoff,
    correlations = correlations,
    plotted = plotted,
    input_files = c(paths$app_cl, paths$expression_columns, paths$ploidy, paths$drug_aliases, used_files)
  )
}

write_outputs <- function(result, base_dir, out_dir, correlation_threshold, drug_filter) {
  slug <- metric_slug(result$metric)
  metadata_dir <- file.path(out_dir, "metadata")
  tables_dir <- file.path(out_dir, "tables")
  qc_dir <- file.path(out_dir, "qc")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  write_run_metadata(
    data.frame(
      key = c(
        "metric",
        "metric_label",
        "metric_source",
        "metric_note",
        "response_column",
        "lower_metric_more_sensitive",
        "sensitivity_direction_assumption",
        "correlation_threshold",
        "drug_coverage_threshold_fraction",
        "drug_coverage_cutoff",
        "drug_filter",
        "module_dir"
      ),
      value = c(
        result$metric,
        result$metric_label,
        result$metric_source,
        result$metric_note,
        result$response_column,
        result$lower_metric_more_sensitive,
        result$sensitivity_direction_assumption,
        correlation_threshold,
        0.8,
        result$coverage_cutoff,
        paste(drug_filter, collapse = ","),
        normalizePath(base_dir, mustWork = FALSE)
      ),
      stringsAsFactors = FALSE
    ),
    metadata_dir
  )
  write_session_metadata(file.path(metadata_dir, "session_info.txt"))
  write_input_manifest(result$input_files, file.path(metadata_dir, "input_manifest.tsv"))

  write_tsv(result$drug_coverage, file.path(tables_dir, paste0("drug_coverage_", slug, ".tsv")))
  write_tsv(result$correlations, file.path(tables_dir, paste0("drug_ploidy_correlations_all_", slug, ".tsv")))
  write_tsv(result$plotted, file.path(tables_dir, paste0("drug_ploidy_correlations_plotted_", slug, ".tsv")))
  write_tsv(result$plotted, file.path(out_dir, paste0("ccle_drug_ploidy_correlations_", slug, ".tsv")))

  plot_file <- file.path(out_dir, paste0("ccle_drug_ploidy_correlations_", slug, ".pdf"))
  plot_correlation_barplot(
    result$plotted,
    plot_file,
    metric_label = result$metric_label,
    lower_metric_more_sensitive = result$lower_metric_more_sensitive
  )
  save(result, file = file.path(out_dir, paste0("drug_ploidy_correlations_", slug, ".RData")))

  summary <- data.frame(
    metric = result$metric,
    metric_label = result$metric_label,
    metric_source = result$metric_source,
    response_column = result$response_column,
    n_cell_lines = ncol(result$response_matrix),
    n_drugs_after_coverage_filter = nrow(result$response_matrix),
    n_plotted_drugs = nrow(result$plotted),
    min_plotted_estimate = min(result$plotted$estimate, na.rm = TRUE),
    max_plotted_estimate = max(result$plotted$estimate, na.rm = TRUE),
    plot_file = basename(plot_file),
    stringsAsFactors = FALSE
  )
  write_tsv(summary, file.path(out_dir, "result_summary.tsv"))

  invisible(list(result = result, summary = summary, plot_file = plot_file))
}

run_ccle_ploidy_analysis <- function(base_dir,
                                     out_dir = default_output_dir(base_dir),
                                     metric = "Z_SCORE",
                                     metric_source = NULL,
                                     correlation_threshold = 0.2,
                                     drug_filter = character()) {
  source(file.path(base_dir, "src", "dependencies.R"))
  source(file.path(base_dir, "src", "analysis_helpers.R"))
  setup_local_lib(base_dir)
  require_packages(character())

  metric_key <- toupper(metric)
  legacy_metric_keys <- c("Z_SCORE", "ZSCORE", "LEGACY_Z_SCORE")
  metric_source <- metric_source %||% if (metric_key %in% legacy_metric_keys) "legacy" else "grbrowser"
  metric_source_key <- tolower(metric_source)

  if (!metric_source_key %in% c("legacy", "grbrowser")) {
    stop("metric_source must be either 'legacy' or 'grbrowser'.", call. = FALSE)
  }

  if (metric_source_key == "legacy") {
    if (metric_key == "IC50") {
      warning(
        "--metric=IC50 --metric-source=legacy is deprecated because the legacy files contain Z Score values. ",
        "Running the legacy analysis as Z_SCORE.",
        call. = FALSE
      )
    } else if (!metric_key %in% legacy_metric_keys) {
      stop("The legacy metric source only supports Z_SCORE.", call. = FALSE)
    }
    result <- run_legacy_z_score(base_dir, out_dir, correlation_threshold, drug_filter)
  } else {
    result <- run_grbrowser_metric(base_dir, out_dir, metric, correlation_threshold, drug_filter)
  }
  write_outputs(result, base_dir, out_dir, correlation_threshold, drug_filter)
}

ccle_ploidy_analysis_main <- function(base_dir) {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- normalizePath(arg_value(args, "output-dir", default_output_dir(base_dir)), mustWork = FALSE)
  metric <- arg_value(args, "metric", "Z_SCORE")
  metric_source <- arg_value(args, "metric-source", arg_value(args, "ic50-source", NULL))
  correlation_threshold <- as.numeric(arg_value(args, "correlation-threshold", "0.2"))
  drug_filter <- split_arg(arg_value(args, "drugs", ""))
  run_ccle_ploidy_analysis(
    base_dir = base_dir,
    out_dir = out_dir,
    metric = metric,
    metric_source = metric_source,
    correlation_threshold = correlation_threshold,
    drug_filter = drug_filter
  )
}
