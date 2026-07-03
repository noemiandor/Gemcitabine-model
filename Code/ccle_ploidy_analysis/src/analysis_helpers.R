`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

read_tsv <- function(path, ...) {
  read.table(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    ...
  )
}

canonical_drug_name <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(x))
}

caps <- function(x) {
  vapply(strsplit(x, " "), function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
  }, character(1))
}

metric_slug <- function(metric) {
  tolower(gsub("[^A-Za-z0-9]+", "_", metric))
}

read_ploidy_table <- function(path) {
  x <- read_tsv(path)
  missing_cols <- setdiff(c("ccle_name", "ploidy"), colnames(x))
  if (length(missing_cols) > 0) {
    stop("Ploidy table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  ploidy <- as.numeric(x$ploidy)
  names(ploidy) <- x$ccle_name
  sort(ploidy)
}

read_expression_columns <- function(path) {
  x <- read_tsv(path)
  if (!"ccle_name" %in% colnames(x)) {
    stop("Expression column table must contain 'ccle_name'.", call. = FALSE)
  }
  x$ccle_name
}

prepare_ploidy <- function(ploidy,
                           appCL,
                           expression_columns,
                           require_primary_adherent = FALSE) {
  if (!"CCLE name" %in% colnames(appCL)) {
    stop("Cell_app_export table is missing 'CCLE name'.", call. = FALSE)
  }
  breast_rows <- grep("BREAST", appCL$`CCLE name`)
  expression_overlap <- intersect(expression_columns, appCL[breast_rows, ]$`CCLE name`)
  ploidy <- sort(ploidy[names(ploidy) %in% expression_overlap])

  if (require_primary_adherent) {
    missing_cols <- setdiff(c("Donor tumor phase", "Growth pattern"), colnames(appCL))
    if (length(missing_cols) > 0) {
      stop("Cell_app_export table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    appCL_idx <- match(names(ploidy), appCL$`CCLE name`)
    keep <- appCL[appCL_idx, "Donor tumor phase"] == "primary" &
      appCL[appCL_idx, "Growth pattern"] != "suspension"
    keep[is.na(keep)] <- FALSE
    ploidy <- sort(ploidy[keep])
  }

  ploidy
}

compute_pearson <- function(response, ploidy, min_n = 3) {
  ok <- is.finite(response) & is.finite(ploidy)
  n <- sum(ok)
  empty <- data.frame(estimate = NA_real_, p.value = NA_real_, n = n)
  if (n < min_n) {
    return(empty)
  }
  response <- response[ok]
  ploidy <- ploidy[ok]
  if (stats::sd(response) == 0 || stats::sd(ploidy) == 0) {
    return(empty)
  }
  test <- stats::cor.test(ploidy, response, method = "pearson")
  data.frame(estimate = unname(test$estimate), p.value = test$p.value, n = n)
}

match_pathways_from_aliases <- function(drugs, alias_file) {
  out <- rep(NA_character_, length(drugs))
  names(out) <- drugs
  if (!file.exists(alias_file)) {
    return(out)
  }

  anno <- read_tsv(alias_file)
  required <- c("drug_name", "synonyms", "pathway_name")
  missing_cols <- setdiff(required, colnames(anno))
  if (length(missing_cols) > 0) {
    stop("DrugAliases file is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  drug_names <- anno$drug_name
  synonyms <- anno$synonyms
  drug_names[is.na(drug_names)] <- ""
  synonyms[is.na(synonyms)] <- ""
  synonyms_no_dash <- gsub("-", "", synonyms)
  for (drug in drugs) {
    ix <- unique(c(
      grep(drug, drug_names, fixed = TRUE),
      grep(drug, synonyms, fixed = TRUE),
      grep(drug, synonyms_no_dash, fixed = TRUE)
    ))
    if (length(ix) > 0) {
      out[drug] <- anno$pathway_name[ix[1]]
    }
  }
  out
}

manual_grbrowser_annotations <- function(drugs) {
  pathway <- rep(NA_character_, length(drugs))
  category <- rep(NA_character_, length(drugs))
  names(pathway) <- names(category) <- drugs

  manual <- list(
    Dichloroacetate = c("pyruvate dehydrogenase kinase", "Signaling"),
    Oxamflatin = c("HDAC inhibitor", "Cytotoxic"),
    AG1024 = c("MAPK/ERK2 signaling", "Signaling"),
    `TPCA-1` = c("IKB kinase inhibitor", "Signaling"),
    IKK16 = c("IKK-16 kinase inhibitor", "Signaling"),
    XRP44X = c("Ras/Erk inhibitor", "Signaling"),
    `TCS 2312` = c("Chk1 inhibitor", "Cytotoxic"),
    `Sigma A6730` = c("Akt1/2 kinase inhibitor", "Signaling"),
    `SB-3CT` = c("MMP2 inhibitor", "Signaling"),
    PD184352 = c("MEK1/2 inhibitor", "Signaling"),
    `PD 98059` = c("MEK inhibitor", "Signaling"),
    `Olomoucine II` = c("CDK/cyclin inhibitor", "Cytotoxic"),
    `L-779450` = c("Raf kinase inhibitor", "Signaling"),
    `JNK-IN-5A` = c("JNK2/JNK3 inhibitor", "Signaling"),
    `Glycyl-H-1152` = c("Rho-kinase inhibitor", "Signaling"),
    `AS-252424` = c("PI3K Inhibitor", "Signaling"),
    `Trichostatin A` = c("HDAC Inhibitor", "Cytotoxic"),
    `CGC-11047` = c("DNA replication", NA_character_),
    `PS-1145` = c("IKB kinase", "Signaling"),
    `Nutlin 3a` = c("MDM2 inhibitor, Apoptosis", "Cytotoxic")
  )

  for (drug in intersect(names(manual), drugs)) {
    pathway[drug] <- manual[[drug]][1]
    category[drug] <- manual[[drug]][2]
  }

  pathway_only <- c(
    lapatinib = "HER1/EGFR/ERBB1 Inhibitor",
    erlotinib = "EGFR Inhibitor",
    pd0325901 = "MEK Inhibitor",
    crizotinib = "ALK/HGFR Inhibitor"
  )
  for (drug in intersect(names(pathway_only), drugs)) {
    pathway[drug] <- pathway_only[[drug]]
  }

  data.frame(drug = drugs, pathway_name = pathway, drugCategory = category, stringsAsFactors = FALSE)
}

finalize_pathway_labels <- function(pathway_name, drug_category = NULL) {
  if (!is.null(drug_category)) {
    pathway_name[is.na(pathway_name)] <- drug_category[is.na(pathway_name)]
  }
  pathway_name[is.na(pathway_name)] <- "Other"
  pathway_name[pathway_name == "Microtubul"] <- "Microtubuli"
  pathway_name[tolower(pathway_name) == "mitotic"] <- "Mitosis"
  pathway_name
}

plot_correlation_barplot <- function(results,
                                     plot_file,
                                     metric_label,
                                     lower_metric_more_sensitive) {
  dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)
  x_axis_label <- if (identical(metric_label, "BreastCancerDrugSensitivity Z Score")) {
    "Pearson r(ploidy, drug-resistance Z-score)"
  } else {
    paste0("Pearson (ploidy, ", metric_label, ")")
  }
  colors <- if (lower_metric_more_sensitive) {
    c("purple", "orange")[1 + (results$estimate > 0)]
  } else {
    c("orange", "purple")[1 + (results$estimate > 0)]
  }

  grDevices::pdf(plot_file, width = 9, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mai = c(1, 1, 0.04, 0.9))
  b <- graphics::barplot(
    results$estimate,
    horiz = TRUE,
    names = results$plot_label,
    las = 2,
    xlab = x_axis_label,
    col = colors,
    border = "white",
    cex.names = 0.7,
    cex.lab = 0.8,
    cex.axis = 0.7
  )
  graphics::axis(
    side = 4,
    labels = gsub(" and ", "/", results$pathway_name),
    las = 2,
    at = b[, 1],
    cex.axis = 0.7
  )
  graphics::plot(1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  graphics::legend(
    "topright",
    c("Low ploidy is sensitive", "High ploidy is sensitive"),
    fill = c("orange", "purple"),
    bty = "n",
    cex = 0.69
  )
  invisible(plot_file)
}
