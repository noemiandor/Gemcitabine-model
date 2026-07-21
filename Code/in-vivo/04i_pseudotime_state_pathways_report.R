#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

default_results_root <- "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04i_pseudotime_state_pathways"
default_hpc_results_root <- "/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels/Results/04i_pseudotime_state_pathways"
default_plan_md <- "/Users/4482173/Documents/Slack/pseudotime_accumulation_state_pathway_analysis_plan.md"
default_repo_root <- "/Users/4482173/Documents/GitHub/Gemcitabine-model"
default_hpc_repo_root <- "/share/lab_crd/lab_crd/taoli/Gemcitabine-model"
default_project_root <- "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels"
default_hpc_project_root <- "/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels"

parse_args <- function(args) {
  out <- list(
    results_root = default_results_root,
    hpc_results_root = default_hpc_results_root,
    output_dir = file.path(default_results_root, "report"),
    plan_md = default_plan_md,
    repo_root = default_repo_root,
    hpc_repo_root = default_hpc_repo_root,
    project_root = default_project_root,
    hpc_project_root = default_hpc_project_root,
    figure_dpi = "160",
    embed_figures = TRUE,
    overwrite = TRUE
  )
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- kv[[1L]]
    value <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
    key <- gsub("-", "_", key)
    if (!key %in% names(out)) stop("Unknown argument: --", key, call. = FALSE)
    if (key %in% c("overwrite", "embed_figures")) {
      out[[key]] <- toupper(value) %in% c("TRUE", "T", "1", "YES", "Y")
    } else {
      out[[key]] <- value
    }
  }
  out
}

normalize_existing <- function(path) {
  if (file.exists(path)) normalizePath(path, winslash = "/", mustWork = TRUE) else path
}

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

replace_prefix <- function(x, from, to) {
  if (!nzchar(from) || !nzchar(to)) return(x)
  from_norm <- normalize_existing(from)
  x <- gsub(from_norm, to, x, fixed = TRUE)
  x <- gsub(from, to, x, fixed = TRUE)
  x
}

sanitize_text <- function(x, cfg) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- replace_prefix(x, cfg$results_root, cfg$hpc_results_root)
  x <- replace_prefix(x, cfg$project_root, cfg$hpc_project_root)
  x <- replace_prefix(x, cfg$repo_root, cfg$hpc_repo_root)
  x <- gsub("/Users/4482173/Documents/Slack", "/share/lab_crd/lab_crd/taoli/Project/Slack", x, fixed = TRUE)
  x
}

display_path <- function(path, cfg) {
  sanitize_text(path, cfg)
}

safe_read_csv <- function(path, strings_as_factors = FALSE) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(
    utils::read.csv(path, stringsAsFactors = strings_as_factors, check.names = FALSE),
    error = function(e) data.frame()
  )
}

safe_read_lines <- function(path) {
  if (!file.exists(path)) return(character())
  tryCatch(readLines(path, warn = FALSE), error = function(e) character())
}

format_value <- function(x) {
  if (length(x) != 1L) x <- x[[1L]]
  if (is.na(x)) return("")
  if (is.numeric(x)) {
    if (!is.finite(x)) return(as.character(x))
    ax <- abs(x)
    if (ax > 0 && ax < 0.001) return(formatC(x, format = "e", digits = 2L))
    if (ax >= 1000) return(formatC(x, format = "f", digits = 0L, big.mark = ","))
    return(formatC(x, format = "f", digits = 3L))
  }
  as.character(x)
}

format_table <- function(df, cfg, max_rows = 20L, caption = NULL, class = "report-table") {
  if (is.null(df) || nrow(df) == 0L || ncol(df) == 0L) {
    return("<p class=\"missing\">No table data available.</p>")
  }
  total_rows <- nrow(df)
  if (total_rows > max_rows) df <- utils::head(df, max_rows)
  df[] <- lapply(df, function(col) {
    vapply(col, function(x) format_value(x), character(1L))
  })
  df[] <- lapply(df, sanitize_text, cfg = cfg)
  header <- paste0("<tr>", paste0("<th>", html_escape(names(df)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1L, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  cap <- if (!is.null(caption)) paste0("<caption>", html_escape(caption), "</caption>") else ""
  note <- if (total_rows > max_rows) {
    paste0("<p class=\"table-note\">Showing first ", max_rows, " of ", total_rows, " rows.</p>")
  } else {
    ""
  }
  paste0(
    "<div class=\"table-wrap\"><table class=\"", class, "\">", cap,
    "<thead>", header, "</thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table></div>", note
  )
}

slugify <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
  x <- gsub("(^-+|-+$)", "", x)
  if (!nzchar(x)) "item" else x
}

section_builder <- function(cfg) {
  env <- new.env(parent = emptyenv())
  env$content <- character()
  env$nav <- data.frame(level = integer(), number = character(), title = character(), id = character(), stringsAsFactors = FALSE)
  env$counters <- integer(6L)
  env$figure_counter <- 0L
  env$used_ids <- character()
  env$section_open <- FALSE
  env
}

unique_id <- function(builder, title, prefix = "sec") {
  base <- paste0(prefix, "-", slugify(title))
  id <- base
  i <- 2L
  while (id %in% builder$used_ids) {
    id <- paste0(base, "-", i)
    i <- i + 1L
  }
  builder$used_ids <- c(builder$used_ids, id)
  id
}

add_html <- function(builder, ...) {
  builder$content <- c(builder$content, paste0(...))
  invisible(NULL)
}

add_heading <- function(builder, level, title, card = FALSE) {
  if (level < 1L || level > 6L) stop("Invalid heading level.", call. = FALSE)
  builder$counters[level] <- builder$counters[level] + 1L
  if (level < 6L) builder$counters[(level + 1L):6L] <- 0L
  number <- paste(builder$counters[seq_len(level)][builder$counters[seq_len(level)] > 0L], collapse = ".")
  id <- unique_id(builder, paste(number, title))
  builder$nav <- rbind(
    builder$nav,
    data.frame(level = level, number = number, title = title, id = id, stringsAsFactors = FALSE)
  )
  cls <- if (card) "report-card" else "report-section"
  if (level == 1L) {
    if (isTRUE(builder$section_open)) add_html(builder, "</section>")
    add_html(builder, "<section class=\"", cls, "\" id=\"", id, "-section\">")
    builder$section_open <- TRUE
  }
  add_html(builder, "<h", level, " id=\"", id, "\"><span class=\"secno\">", number, "</span> ", html_escape(title), "</h", level, ">")
  invisible(id)
}

add_paragraph <- function(builder, text, class = NULL) {
  cls <- if (!is.null(class)) paste0(" class=\"", class, "\"") else ""
  add_html(builder, "<p", cls, ">", text, "</p>")
}

add_bullets <- function(builder, items) {
  items <- items[nzchar(items)]
  if (length(items) == 0L) return(invisible(NULL))
  add_html(builder, "<ul>", paste0("<li>", items, "</li>", collapse = "\n"), "</ul>")
}

base64_encode_file <- function(path) {
  size <- file.info(path)$size
  if (!is.finite(size) || size <= 0L) return("")
  bytes <- as.integer(readBin(path, what = "raw", n = size))
  n <- length(bytes)
  if (n == 0L) return("")
  pad <- (3L - (n %% 3L)) %% 3L
  if (pad > 0L) bytes <- c(bytes, rep(0L, pad))
  mat <- matrix(bytes, ncol = 3L, byrow = TRUE)
  alphabet <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "", fixed = TRUE)[[1L]]
  idx1 <- bitwShiftR(mat[, 1L], 2L)
  idx2 <- bitwOr(bitwShiftL(bitwAnd(mat[, 1L], 3L), 4L), bitwShiftR(mat[, 2L], 4L))
  idx3 <- bitwOr(bitwShiftL(bitwAnd(mat[, 2L], 15L), 2L), bitwShiftR(mat[, 3L], 6L))
  idx4 <- bitwAnd(mat[, 3L], 63L)
  out <- character(nrow(mat) * 4L)
  out[seq(1L, length(out), by = 4L)] <- alphabet[idx1 + 1L]
  out[seq(2L, length(out), by = 4L)] <- alphabet[idx2 + 1L]
  out[seq(3L, length(out), by = 4L)] <- alphabet[idx3 + 1L]
  out[seq(4L, length(out), by = 4L)] <- alphabet[idx4 + 1L]
  if (pad > 0L) out[(length(out) - pad + 1L):length(out)] <- "="
  paste(out, collapse = "")
}

image_mime_type <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("jpg", "jpeg")) return("image/jpeg")
  if (identical(ext, "gif")) return("image/gif")
  if (identical(ext, "svg")) return("image/svg+xml")
  if (identical(ext, "webp")) return("image/webp")
  "image/png"
}

image_data_uri <- function(path) {
  paste0("data:", image_mime_type(path), ";base64,", base64_encode_file(path))
}

copy_figure_asset <- function(src, cfg, label) {
  if (!file.exists(src)) return(NULL)
  if (isTRUE(cfg$embed_figures)) {
    dest_dir <- file.path(tempdir(), "04i_pseudotime_state_pathways_report_figures")
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    ext <- tools::file_ext(src)
    if (tolower(ext) == "pdf") {
      png <- convert_pdf_to_png(src, dest_dir, label, as.integer(cfg$figure_dpi %||% "160"))
      if (!is.null(png) && file.exists(png)) return(image_data_uri(png))
      return(NULL)
    }
    if (tolower(ext) %in% c("png", "jpg", "jpeg", "gif", "svg", "webp")) {
      return(image_data_uri(src))
    }
    return(NULL)
  }
  dest_dir <- file.path(cfg$output_dir, "assets", "figures")
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  ext <- tools::file_ext(src)
  if (tolower(ext) == "pdf") {
    png <- convert_pdf_to_png(src, dest_dir, label, as.integer(cfg$figure_dpi %||% "160"))
    if (!is.null(png)) return(file.path("assets", "figures", basename(png)))
  }
  name <- paste0(slugify(label), ".", ext)
  dest <- file.path(dest_dir, name)
  if (file.exists(dest) && !isTRUE(cfg$overwrite)) return(file.path("assets", "figures", name))
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!ok) return(NULL)
  file.path("assets", "figures", name)
}

convert_pdf_to_png <- function(src, dest_dir, label, dpi = 160L) {
  dpi <- if (is.finite(dpi) && dpi > 0L) dpi else 160L
  base <- file.path(dest_dir, slugify(label))
  dest <- paste0(base, ".png")
  if (file.exists(dest)) return(dest)

  magick <- Sys.which("magick")
  if (nzchar(magick)) {
    cmd <- paste(
      shQuote(magick),
      "-density", shQuote(as.character(dpi)),
      shQuote(paste0(src, "[0]")),
      "-quality 95",
      shQuote(dest)
    )
    status <- tryCatch(
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE),
      warning = function(w) 1L,
      error = function(e) 1L
    )
    if (identical(status, 0L) && file.exists(dest)) return(dest)
  }

  pdftoppm <- Sys.which("pdftoppm")
  if (nzchar(pdftoppm)) {
    cmd <- paste(
      shQuote(pdftoppm),
      "-png -r", shQuote(as.character(dpi)),
      "-singlefile",
      shQuote(src),
      shQuote(base)
    )
    status <- tryCatch(
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE),
      warning = function(w) 1L,
      error = function(e) 1L
    )
    if (identical(status, 0L) && file.exists(dest)) return(dest)
  }

  NULL
}

add_figure <- function(builder, cfg, src, title, legend, interpretation, data_paths = character(), anchor_label = NULL) {
  builder$figure_counter <- builder$figure_counter + 1L
  fig_no <- builder$figure_counter
  fig_id <- unique_id(builder, paste("figure", fig_no, title), prefix = "fig")
  rel <- copy_figure_asset(src, cfg, paste0("figure-", fig_no, "-", anchor_label %||% title))
  source_text <- display_path(src, cfg)
  data_text <- if (length(data_paths) > 0L) {
    paste0(
      "<ul class=\"path-list\">",
      paste0("<li><code>", html_escape(display_path(data_paths, cfg)), "</code></li>", collapse = "\n"),
      "</ul>"
    )
  } else {
    "<p class=\"muted\">No plot-data table was registered for this figure.</p>"
  }
  media <- if (is.null(rel)) {
    paste0("<div class=\"missing\">Figure file not found: <code>", html_escape(source_text), "</code></div>")
  } else if (grepl("^data:image/", rel) || tolower(tools::file_ext(rel)) %in% c("png", "jpg", "jpeg", "gif", "svg", "webp")) {
    paste0("<img src=\"", html_escape(rel), "\" alt=\"", html_escape(title), "\"/>")
  } else if (tolower(tools::file_ext(rel)) == "pdf") {
    paste0(
      "<object class=\"pdf-object\" data=\"", html_escape(rel), "\" type=\"application/pdf\">",
      "<p>PDF preview is not available in this browser. <a href=\"", html_escape(rel), "\">Open the PDF asset</a>.</p>",
      "</object>"
    )
  } else {
    paste0("<a href=\"", html_escape(rel), "\">Open figure asset</a>")
  }
  add_html(
    builder,
    "<figure class=\"report-figure\" id=\"", fig_id, "\">",
    "<div class=\"figure-media\">", media, "</div>",
    "<figcaption>",
    "<p class=\"figure-title\"><strong>Figure ", fig_no, ". ", html_escape(title), ".</strong></p>",
    "<p><strong>Legend.</strong> ", legend, "</p>",
    "<p><strong>Interpretation.</strong> ", interpretation, "</p>",
    "<p><strong>Source figure.</strong> <code>", html_escape(source_text), "</code></p>",
    "<div><strong>Plot data / statistics.</strong>", data_text, "</div>",
    "</figcaption></figure>"
  )
  invisible(fig_no)
}

model_order <- c(
  "primary_initial_ploidy",
  "etp_continuous_adjusted",
  "ETP_fixed_threshold_2_25",
  "ETP_boundary_stress_threshold_2_375",
  "ETP_reference_balanced_threshold_2_24"
)

model_display_name <- function(model_id, specs) {
  if (!is.null(specs) && nrow(specs) > 0L && model_id %in% specs$model_id) {
    label <- specs$model_label[match(model_id, specs$model_id)]
    if (nzchar(label)) return(label)
  }
  gsub("_", " ", model_id)
}

select_top_pathways <- function(gsea, n_each = 5L) {
  if (nrow(gsea) == 0L) return(data.frame())
  needed <- c("collection_label", "pathway_label", "NES", "padj")
  if (!all(needed %in% names(gsea))) return(data.frame())
  rows <- lapply(split(gsea, gsea$collection_label), function(df) {
    pos <- df[df$NES > 0 & is.finite(df$padj), , drop = FALSE]
    neg <- df[df$NES < 0 & is.finite(df$padj), , drop = FALSE]
    pos <- utils::head(pos[order(pos$padj, -pos$NES), , drop = FALSE], n_each)
    neg <- utils::head(neg[order(neg$padj, neg$NES), , drop = FALSE], n_each)
    if (nrow(pos) > 0L) pos$direction <- "positive"
    if (nrow(neg) > 0L) neg$direction <- "negative"
    rbind(pos, neg)
  })
  out <- do.call(rbind, rows)
  keep <- intersect(c("collection_label", "direction", "pathway_label", "NES", "padj", "size"), names(out))
  out <- out[, keep, drop = FALSE]
  out[order(out$collection_label, out$direction, out$padj, -abs(out$NES)), , drop = FALSE]
}

significant_counts <- function(gsea) {
  if (nrow(gsea) == 0L || !all(c("collection_label", "NES", "padj") %in% names(gsea))) return(data.frame())
  gsea$direction <- ifelse(gsea$NES >= 0, "positive", "negative")
  gsea$significant <- is.finite(gsea$padj) & gsea$padj < 0.05
  agg <- stats::aggregate(significant ~ collection_label + direction, gsea, sum)
  names(agg)[names(agg) == "significant"] <- "FDR_lt_0_05_pathways"
  total <- stats::aggregate(pathway_label ~ collection_label, gsea, length)
  names(total)[names(total) == "pathway_label"] <- "tested_pathways"
  merge(agg, total, by = "collection_label", all.x = TRUE, sort = FALSE)
}

top_gene_table <- function(genes, n_each = 8L) {
  if (nrow(genes) == 0L || !"t_statistic" %in% names(genes)) return(data.frame())
  pos <- genes[genes$t_statistic > 0, , drop = FALSE]
  neg <- genes[genes$t_statistic < 0, , drop = FALSE]
  pos <- utils::head(pos[order(pos$fdr, -pos$t_statistic), , drop = FALSE], n_each)
  neg <- utils::head(neg[order(neg$fdr, neg$t_statistic), , drop = FALSE], n_each)
  if (nrow(pos) > 0L) pos$direction <- "positive"
  if (nrow(neg) > 0L) neg$direction <- "negative"
  out <- rbind(pos, neg)
  keep <- intersect(c("direction", "gene_symbol", "fitted_log_expr_difference", "t_statistic", "p_value", "fdr"), names(out))
  out[, keep, drop = FALSE]
}

pathway_summary_sentence <- function(gsea) {
  if (nrow(gsea) == 0L) return("No primary-state GSEA table was available for this model.")
  pos <- gsea[gsea$NES > 0 & is.finite(gsea$padj), , drop = FALSE]
  neg <- gsea[gsea$NES < 0 & is.finite(gsea$padj), , drop = FALSE]
  pos <- utils::head(pos[order(pos$padj, -pos$NES), , drop = FALSE], 3L)
  neg <- utils::head(neg[order(neg$padj, neg$NES), , drop = FALSE], 3L)
  fmt <- function(df) {
    if (nrow(df) == 0L) return("none")
    paste0(
      html_escape(df$pathway_label), " (NES ", vapply(df$NES, format_value, character(1L)),
      ", FDR ", vapply(df$padj, format_value, character(1L)), ")",
      collapse = "; "
    )
  }
  paste0(
    "The strongest positive primary-state pathways are ", fmt(pos),
    ". The strongest negative primary-state pathways are ", fmt(neg),
    ". Positive NES means higher fitted expression in the frozen primary state than in the adjacent flanks; negative NES means lower fitted expression in the frozen primary state."
  )
}

registered_model_figures <- function(model_dir) {
  fig_dir <- file.path(model_dir, "05_figures")
  list(
    list(
      file = file.path(fig_dir, "pseudotime_density_with_frozen_state.pdf"),
      title = "Pseudotime density with frozen state intervals",
      legend = "Colored curves summarize equal-mouse pseudotime densities by treatment group. The black curve is treated minus control density difference. Shaded regions mark the frozen primary, sensitivity, and negative-control intervals.",
      interpretation = "This figure documents the interval-selection evidence. It is used once to freeze the accumulated state; the downstream expression model does not re-test treatment versus control expression.",
      data = file.path(fig_dir, c(
        "pseudotime_density_with_frozen_state_mean_density_plot_data.csv",
        "pseudotime_density_with_frozen_state_difference_plot_data.csv",
        "pseudotime_density_with_frozen_state_sample_density_plot_data.csv"
      ))
    ),
    list(
      file = file.path(fig_dir, "sample_region_coverage.pdf"),
      title = "Sample coverage in frozen pseudotime regions",
      legend = "Bars show retained cells or bins for each mouse across the frozen primary, flank, sensitivity, and negative-control regions.",
      interpretation = "This figure checks whether the primary contrast is estimable without relying on a single sample, treatment group, or ploidy group.",
      data = file.path(fig_dir, "sample_region_coverage_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_hallmark_gsea.pdf"),
      title = "Primary-state Hallmark GSEA",
      legend = "Bars show Hallmark normalized enrichment scores for the primary adjacent-state contrast. Fill encodes FDR.",
      interpretation = "Hallmark pathways give a compact view of broad biological programs that characterize, or are depleted from, the accumulated state.",
      data = file.path(fig_dir, "primary_state_hallmark_gsea_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_reactome_gsea.pdf"),
      title = "Primary-state Reactome GSEA",
      legend = "Bars show Reactome normalized enrichment scores for the primary adjacent-state contrast. Fill encodes FDR.",
      interpretation = "Reactome pathways provide mechanistic pathway-level annotation of the accumulated state.",
      data = file.path(fig_dir, "primary_state_reactome_gsea_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_go_bp_gsea.pdf"),
      title = "Primary-state GO Biological Process GSEA",
      legend = "Bars show GO Biological Process normalized enrichment scores for the primary adjacent-state contrast. Fill encodes FDR.",
      interpretation = "GO Biological Process terms provide a broad ontology-level view of state-associated biological processes.",
      data = file.path(fig_dir, "primary_state_go_bp_gsea_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_pathway_activity_heatmap.pdf"),
      title = "Primary-state pathway activity heatmap",
      legend = "Columns are pseudotime grid points, rows are top positive and negative pathways, and color is the mean standardized fitted expression of each pathway leading-edge gene set. Dashed vertical lines mark pseudotime 0.30 and 0.49.",
      interpretation = "This heatmap shows where each enriched pathway is active along the trajectory. Primary-state-positive pathways should peak inside the frozen interval; primary-state-negative pathways are expected to peak outside it or be low within it.",
      data = file.path(fig_dir, "primary_state_pathway_activity_heatmap_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_top_pathway_curves.pdf"),
      title = "Top pathway fitted activity curves",
      legend = "Each curve is the mean standardized fitted activity of leading-edge genes for a top pathway. The gray band marks the frozen primary interval.",
      interpretation = "The curves show whether pathway activity is localized to the accumulated state or reflects a broader trajectory trend.",
      data = file.path(fig_dir, "primary_state_top_pathway_curves_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_leading_edge_gene_curves.pdf"),
      title = "Leading-edge gene fitted curves",
      legend = "Each curve is a standardized fitted expression curve for a leading-edge gene from top positive primary-state pathways. The gray band marks the frozen primary interval.",
      interpretation = "This figure checks whether the positive pathway signal is supported by coherent gene-level curves rather than a single outlier gene.",
      data = file.path(fig_dir, "primary_state_leading_edge_gene_curves_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "pathway_robustness_summary.pdf"),
      title = "Pathway robustness across sensitivity analyses",
      legend = "Points show NES values for pre-specified sensitivity intervals and modeling variants, grouped by pathway.",
      interpretation = "Robust state annotations should preserve the sign and approximate rank of the primary-state signal across sensitivity checks.",
      data = file.path(fig_dir, "pathway_robustness_summary_plot_data.csv")
    )
  )
}

registered_comparison_figures <- function(root) {
  fig_dir <- file.path(root, "model_comparison", "figures")
  list(
    list(
      file = file.path(fig_dir, "gsea_NES_correlation_primary_vs_ETP.pdf"),
      title = "GSEA NES correlation across covariate models",
      legend = "Points compare pathway NES values between the reference initial-ploidy model and each ETP-adjusted or ETP-threshold model.",
      interpretation = "High correlations indicate that adding ETP as a continuous covariate or threshold group does not materially change the primary-state pathway signature.",
      data = file.path(root, "model_comparison", "gsea_NES_correlations.csv")
    ),
    list(
      file = file.path(fig_dir, "top_pathway_rank_shift_primary_vs_ETP.pdf"),
      title = "Top pathway rank shifts across covariate models",
      legend = "Points show rank movement for pathways when comparing ETP models against the initial-ploidy reference model.",
      interpretation = "Small rank shifts among top pathways support model robustness; large isolated shifts identify model-sensitive annotations that should be interpreted cautiously.",
      data = file.path(root, "model_comparison", "pathway_rank_shift.csv")
    )
  )
}

write_css_js <- function(cfg) {
  css <- paste0(
    "body{margin:0;font-family:-apple-system,BlinkMacSystemFont,\"Segoe UI\",sans-serif;background:#f4f7fa;color:#1b2a38;}",
    ".report-shell{display:flex;gap:28px;max-width:1680px;margin:0 auto;padding:24px;}",
    ".report-sidebar{position:sticky;top:24px;align-self:flex-start;width:330px;max-height:calc(100vh - 48px);border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,.08);overflow:auto;}",
    ".report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff;}",
    ".report-kicker{font-size:11px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;opacity:.78;}",
    ".report-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15;}",
    ".report-subtitle{margin-top:4px;font-size:12px;opacity:.85;}",
    ".report-nav{padding:10px 8px 12px 8px;}.report-nav-list{margin:0;padding:0;list-style:none;}.report-nav-item{margin:3px 0;}",
    ".nav-row{display:flex;align-items:center;gap:4px;border-radius:8px;}.nav-toggle{border:0;background:transparent;color:#536577;cursor:pointer;width:18px;padding:0;font-size:11px;}.nav-spacer{display:inline-block;width:18px;}",
    ".report-nav-link{display:block;flex:1;padding:8px 8px;border-radius:8px;text-decoration:none;color:#17324c;font-size:13px;font-weight:600;line-height:1.35;}",
    ".report-nav-link:hover,.report-nav-link.active{background:rgba(47,110,164,.12);}.report-nav-2{padding-left:6px;}.report-nav-3{padding-left:20px;font-size:12px;font-weight:500;}.report-nav-4{padding-left:34px;font-size:12px;font-weight:500;color:#536577;}.report-nav-5{padding-left:48px;font-size:11px;font-weight:500;color:#6a7a89;}",
    ".children.collapsed{display:none;}.report-main{flex:1;min-width:0;max-width:1180px;}.report-card,.report-section{margin-bottom:24px;padding:20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.05);}",
    ".report-card h1{margin:0 0 8px 0;font-size:28px;line-height:1.15;}.report-section h1{margin-top:0;font-size:25px;line-height:1.18;}h1,h2,h3,h4,h5,h6{scroll-margin-top:24px;color:#17324c;}h2{font-size:20px;margin-top:22px;}h3{font-size:17px;margin-top:20px;}h4{font-size:15px;margin-top:18px;color:#284662;}.secno{color:#68798a;font-weight:700;margin-right:4px;}",
    ".muted{color:#617184;}.missing{padding:10px 12px;border:1px solid #f0c36a;background:#fff8e7;border-radius:8px;color:#6d4b00;}.callout{padding:12px 14px;border-left:4px solid #2f6ea4;background:#eef6fc;border-radius:8px;}",
    ".table-wrap{overflow:auto;margin:12px 0;}table.report-table{border-collapse:collapse;width:100%;font-size:12px;}table.report-table caption{text-align:left;font-weight:700;margin-bottom:6px;color:#284662;}table.report-table th,table.report-table td{border:1px solid #dfe6ee;padding:6px 8px;vertical-align:top;}table.report-table th{background:#edf3f8;color:#17324c;position:sticky;top:0;}table.report-table tr:nth-child(even){background:#fafcff;}.table-note{font-size:12px;color:#617184;margin-top:-4px;}",
    ".report-figure{margin:20px 0;padding:14px;border:1px solid #dfe6ee;border-radius:12px;background:#fbfdff;}.figure-media{border:1px solid #ccd7e2;border-radius:8px;background:#fff;overflow:hidden;}.figure-media img{max-width:100%;display:block;margin:0 auto;}.pdf-object{width:100%;height:760px;display:block;}figcaption{font-size:13px;color:#2b3f52;line-height:1.45;margin-top:10px;}.figure-title{font-size:14px;color:#17324c;}.path-list{margin:6px 0 0 0;padding-left:20px;}code{font-family:\"SFMono-Regular\",Consolas,monospace;font-size:12px;background:#eef2f6;padding:1px 4px;border-radius:4px;}",
    ".figure-grid{display:grid;gap:14px;margin:20px 0;align-items:start;}.figure-grid-3{grid-template-columns:repeat(3,minmax(0,1fr));}.figure-grid .report-figure{margin:0;height:100%;}.figure-grid .figure-media img{width:100%;}.figure-grid figcaption{font-size:12px;}.figure-grid .figure-title{font-size:13px;}@media(max-width:1200px){.figure-grid-3{grid-template-columns:1fr;}}",
    ".pill{display:inline-block;padding:2px 7px;border-radius:999px;background:#e8f0f7;color:#284662;font-size:12px;font-weight:700;margin-right:4px;}.footer{font-size:12px;color:#617184;margin:28px 0 8px 0;}"
  )
  js <- paste0(
    "document.addEventListener('DOMContentLoaded',function(){",
    "function childList(li){for(var i=0;i<li.children.length;i++){var el=li.children[i];if(el.tagName==='UL'&&el.classList.contains('children'))return el;}return null;}",
    "function directToggle(li){var row=li.firstElementChild;if(!row)return null;for(var i=0;i<row.children.length;i++){var el=row.children[i];if(el.classList&&el.classList.contains('nav-toggle'))return el;}return null;}",
    "function setCollapsed(li,collapsed){var kids=childList(li);var btn=directToggle(li);if(!kids)return;if(collapsed){kids.classList.add('collapsed');if(btn)btn.textContent='+';}else{kids.classList.remove('collapsed');if(btn)btn.textContent='-';}}",
    "function topItem(li){var cur=li,last=li;while(cur){if(cur.parentElement&&cur.parentElement.classList&&cur.parentElement.classList.contains('report-nav-list'))return cur;last=cur;cur=cur.parentElement?cur.parentElement.closest('li.report-nav-item'):null;}return last;}",
    "function expandBranch(link){var li=link.closest('li.report-nav-item');if(!li)return;var top=topItem(li);document.querySelectorAll('.report-nav-list>li.report-nav-item').forEach(function(item){if(item!==top)setCollapsed(item,true);});var cur=li;while(cur){setCollapsed(cur,false);cur=cur.parentElement?cur.parentElement.closest('li.report-nav-item'):null;}}",
    "function activate(link){document.querySelectorAll('.report-nav-link').forEach(function(a){a.classList.remove('active');});if(link){link.classList.add('active');expandBranch(link);}}",
    "document.querySelectorAll('li.report-nav-item').forEach(function(li){if(childList(li))setCollapsed(li,true);});",
    "document.querySelectorAll('.nav-toggle').forEach(function(btn){btn.addEventListener('click',function(e){e.preventDefault();e.stopPropagation();var li=btn.closest('li.report-nav-item');var kids=li?childList(li):null;if(kids)setCollapsed(li,!kids.classList.contains('collapsed'));});});",
    "var links=[].slice.call(document.querySelectorAll('.report-nav-link'));",
    "links.forEach(function(a){a.addEventListener('click',function(e){var id=a.getAttribute('href').slice(1);var target=document.getElementById(id);if(target){e.preventDefault();activate(a);target.scrollIntoView({behavior:'smooth',block:'start'});if(history.pushState)history.pushState(null,'','#'+id);else location.hash=id;}});});",
    "var first=links[0];var hash=window.location.hash?window.location.hash.slice(1):'';var initial=hash?document.querySelector('.report-nav-link[href=\"#'+hash+'\"]'):first;if(initial)activate(initial);",
    "var ids=links.map(function(a){return document.getElementById(a.getAttribute('href').slice(1));}).filter(Boolean);",
    "if('IntersectionObserver' in window){var obs=new IntersectionObserver(function(entries){entries.forEach(function(e){if(e.isIntersecting){var a=document.querySelector('.report-nav-link[href=\"#'+e.target.id+'\"]');if(a)activate(a);}});},{rootMargin:'-18% 0px -72% 0px',threshold:0});ids.forEach(function(id){obs.observe(id);});}",
    "});"
  )
  list(css = css, js = js)
}

build_nav <- function(nav) {
  if (nrow(nav) == 0L) return("")
  cursor <- new.env(parent = emptyenv())
  cursor$i <- 1L
  render_level <- function(level) {
    cls <- if (identical(level, 1L)) "report-nav-list" else "children"
    out <- c(paste0("<ul class=\"", cls, "\">"))
    while (cursor$i <= nrow(nav)) {
      current_level <- nav$level[[cursor$i]]
      if (current_level < level) break
      if (current_level > level) {
        out <- c(out, render_level(current_level))
        next
      }
      row <- cursor$i
      cursor$i <- cursor$i + 1L
      has_child <- cursor$i <= nrow(nav) && nav$level[[cursor$i]] > level
      toggle <- if (has_child) "<button class=\"nav-toggle\" type=\"button\" aria-label=\"Toggle section\">-</button>" else "<span class=\"nav-spacer\"></span>"
      link <- paste0(
        "<a class=\"report-nav-link report-nav-", level, "\" href=\"#", html_escape(nav$id[[row]]), "\">",
        html_escape(nav$number[[row]]), " ", html_escape(nav$title[[row]]), "</a>"
      )
      out <- c(out, paste0("<li class=\"report-nav-item report-nav-level-", level, "\"><div class=\"nav-row\">", toggle, link, "</div>"))
      if (has_child) out <- c(out, render_level(nav$level[[cursor$i]]))
      out <- c(out, "</li>")
    }
    c(out, "</ul>")
  }
  paste(render_level(nav$level[[1L]]), collapse = "\n")
}

add_model_section <- function(builder, cfg, model_id, specs) {
  model_dir <- file.path(cfg$results_root, model_id)
  label <- model_display_name(model_id, specs)
  add_heading(builder, 1L, paste0(model_id, ": ", label))
  if (!dir.exists(model_dir)) {
    add_paragraph(builder, paste0("The model output directory is missing: <code>", html_escape(display_path(model_dir, cfg)), "</code>"), "missing")
    return(invisible(NULL))
  }
  model_spec <- if (nrow(specs) > 0L && model_id %in% specs$model_id) specs[match(model_id, specs$model_id), , drop = FALSE] else data.frame()
  model_params <- safe_read_csv(file.path(model_dir, "00_manifest", "model_parameters.csv"))
  design <- safe_read_csv(file.path(model_dir, "01_qc", "model_design_rank_audit.csv"))
  genes <- safe_read_csv(file.path(model_dir, "03_gene_models", "gene_primary_adjacent_state_contrast.csv"))
  gsea <- safe_read_csv(file.path(model_dir, "04_gsea", "all_collections_primary_adjacent_state_gsea.csv"))
  robust <- safe_read_csv(file.path(model_dir, "06_sensitivity", "pathway_robustness_summary.csv"))
  loo <- safe_read_csv(file.path(model_dir, "06_sensitivity", "leave_one_mouse_out_pathway_stability.csv"))

  add_heading(builder, 2L, "Model definition and design audit")
  add_paragraph(
    builder,
    paste0(
      "This section reports the primary adjacent-state contrast for <span class=\"pill\">", html_escape(model_id),
      "</span>. The model treats dose and the selected ploidy/ETP variable as additive nuisance terms; it does not include treatment-by-pseudotime, dose-by-pseudotime, or ploidy/ETP-by-pseudotime interactions."
    )
  )
  if (nrow(model_spec) > 0L) add_html(builder, format_table(model_spec, cfg, max_rows = 5L, caption = "Model specification."))
  if (nrow(model_params) > 0L) add_html(builder, format_table(model_params, cfg, max_rows = 20L, caption = "Model-level run parameters."))
  if (nrow(design) > 0L) add_html(builder, format_table(design, cfg, max_rows = 5L, caption = "Design-rank audit."))

  add_heading(builder, 2L, "Primary-state gene and pathway summary")
  add_paragraph(builder, pathway_summary_sentence(gsea), "callout")
  sig <- significant_counts(gsea)
  if (nrow(sig) > 0L) add_html(builder, format_table(sig, cfg, max_rows = 20L, caption = "FDR < 0.05 pathway counts by collection and direction."))
  top_pw <- select_top_pathways(gsea, n_each = 5L)
  if (nrow(top_pw) > 0L) add_html(builder, format_table(top_pw, cfg, max_rows = 40L, caption = "Top primary-state pathways by collection and direction."))
  top_genes <- top_gene_table(genes, n_each = 8L)
  if (nrow(top_genes) > 0L) add_html(builder, format_table(top_genes, cfg, max_rows = 20L, caption = "Top primary-state genes by t statistic and FDR."))

  add_heading(builder, 2L, "Key figures")
  figs <- registered_model_figures(model_dir)
  gsea_files <- c(
    "primary_state_hallmark_gsea.pdf",
    "primary_state_reactome_gsea.pdf",
    "primary_state_go_bp_gsea.pdf"
  )
  gsea_idx <- which(basename(vapply(figs, `[[`, character(1L), "file")) %in% gsea_files)
  for (i in seq_along(figs)) {
    fig <- figs[[i]]
    if (i %in% gsea_idx && identical(i, gsea_idx[[1L]])) {
      add_html(builder, "<div class=\"figure-grid figure-grid-3\">")
      for (j in gsea_idx) {
        local_fig <- figs[[j]]
        add_figure(
          builder, cfg,
          src = local_fig$file,
          title = paste0(model_id, ": ", local_fig$title),
          legend = local_fig$legend,
          interpretation = local_fig$interpretation,
          data_paths = local_fig$data,
          anchor_label = paste0(model_id, "-", basename(local_fig$file))
        )
      }
      add_html(builder, "</div>")
      next
    }
    if (i %in% gsea_idx) next
    add_figure(
      builder, cfg,
      src = fig$file,
      title = paste0(model_id, ": ", fig$title),
      legend = fig$legend,
      interpretation = fig$interpretation,
      data_paths = fig$data,
      anchor_label = paste0(model_id, "-", basename(fig$file))
    )
  }

  add_heading(builder, 2L, "Sensitivity and leave-one-mouse-out checks")
  if (nrow(robust) > 0L) {
    keep <- intersect(c("collection_label", "pathway_label", "analysis_id", "direction", "NES", "padj", "primary_NES", "primary_padj", "sign_concordant_with_primary"), names(robust))
    add_html(builder, format_table(utils::head(robust[, keep, drop = FALSE], 30L), cfg, max_rows = 30L, caption = "Sensitivity GSEA summary, first rows."))
  } else {
    add_paragraph(builder, "No sensitivity table was available for this model.", "missing")
  }
  if (nrow(loo) > 0L) {
    keep <- intersect(c("collection_label", "pathway_label", "primary_NES", "primary_padj", "n_leave_one_out_runs", "fraction_sign_concordant", "min_NES", "max_NES"), names(loo))
    add_html(builder, format_table(utils::head(loo[, keep, drop = FALSE], 30L), cfg, max_rows = 30L, caption = "Leave-one-mouse-out pathway stability, first rows."))
  } else {
    add_paragraph(builder, "No leave-one-mouse-out stability table was available for this model.", "missing")
  }
  invisible(NULL)
}

summarize_comparison <- function(gene_cor, nes_cor, overlap) {
  parts <- character()
  if (nrow(gene_cor) > 0L && "model_id" %in% names(gene_cor)) {
    x <- gene_cor[gene_cor$model_id != "primary_initial_ploidy", , drop = FALSE]
    if (nrow(x) > 0L) {
      parts <- c(parts, paste0("Gene-level primary-contrast t statistics remain highly concordant with the reference model: minimum Pearson r = ", format_value(min(x$pearson_r, na.rm = TRUE)), "."))
    }
  }
  if (nrow(nes_cor) > 0L && "model_id" %in% names(nes_cor)) {
    x <- nes_cor[nes_cor$model_id != "primary_initial_ploidy", , drop = FALSE]
    if (nrow(x) > 0L) {
      by_coll <- stats::aggregate(pearson_r ~ collection_label, x, min, na.rm = TRUE)
      txt <- paste0(by_coll$collection_label, " min r = ", vapply(by_coll$pearson_r, format_value, character(1L)), collapse = "; ")
      parts <- c(parts, paste0("Pathway NES concordance is also high across collections: ", txt, "."))
    }
  }
  if (nrow(overlap) > 0L && "model_id" %in% names(overlap)) {
    x <- overlap[overlap$model_id != "primary_initial_ploidy", , drop = FALSE]
    if (nrow(x) > 0L) {
      by_coll <- stats::aggregate(jaccard ~ collection_label, x, min, na.rm = TRUE)
      txt <- paste0(by_coll$collection_label, " min Jaccard = ", vapply(by_coll$jaccard, format_value, character(1L)), collapse = "; ")
      parts <- c(parts, paste0("The FDR-significant pathway set is stable under ETP adjustment/grouping: ", txt, "."))
    }
  }
  if (length(parts) == 0L) "Model-comparison summary tables were not available."
  else paste(parts, collapse = " ")
}

build_report <- function(cfg) {
  cfg$results_root <- normalize_existing(cfg$results_root)
  cfg$output_dir <- normalizePath(cfg$output_dir, winslash = "/", mustWork = FALSE)
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(cfg$overwrite)) {
    unlink(file.path(cfg$output_dir, "assets"), recursive = TRUE, force = TRUE)
  }
  page_assets <- write_css_js(cfg)

  b <- section_builder(cfg)
  params <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "analysis_parameters.csv"))
  input_checksums <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "input_checksums.csv"))
  package_versions <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "package_versions.csv"))
  intervals <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "frozen_interval_definition.csv"))
  specs <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "selected_model_specifications.csv"))
  design_all <- safe_read_csv(file.path(cfg$results_root, "model_comparison", "model_design_rank_audit_all_models.csv"))
  expression_audit <- safe_read_csv(file.path(cfg$results_root, "01_qc", "expression_source_audit.csv"))
  join_audit <- safe_read_csv(file.path(cfg$results_root, "01_qc", "cell_id_join_audit.csv"))
  coverage <- safe_read_csv(file.path(cfg$results_root, "01_qc", "primary_coverage_check.csv"))
  etp_cross <- safe_read_csv(file.path(cfg$results_root, "01_qc", "etp_group_cross_tab.csv"))
  plan_lines <- safe_read_lines(cfg$plan_md)

  if (nrow(specs) == 0L) {
    specs <- data.frame(model_id = model_order, model_label = model_order, stringsAsFactors = FALSE)
  }
  ordered_models <- model_order[model_order %in% specs$model_id]
  extras <- setdiff(specs$model_id, ordered_models)
  ordered_models <- c(ordered_models, extras)

  add_html(
    b,
    "<section class=\"report-card\" id=\"report-title\">",
    "<h1>Pseudotime accumulation-state pathway analysis report</h1>",
    "<p class=\"muted\">Generated ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    ". Reported result root: <code>", html_escape(display_path(cfg$results_root, cfg)), "</code>.</p>",
    "<p>This report summarizes the 04i workflow outputs. It annotates the pre-specified CellCycle pseudotime accumulation state using a sample-aware pseudobulk expression model and GSEA. The report is generated from local result files, but all displayed absolute paths are mapped to the HPC path namespace.</p>",
    "</section>"
  )

  add_heading(b, 1L, "Executive summary")
  add_bullets(
    b,
    c(
      "The primary estimand is fitted expression inside the frozen primary accumulated-state interval minus the equally weighted left and right neighboring intervals.",
      "The pathway results should be interpreted as state-associated expression programs, not as direct treatment-induced up- or down-regulation at fixed pseudotime.",
      "Across ETP continuous and ETP threshold models, the primary-state pathway signature is expected to be interpreted through concordance with the primary initial-ploidy-adjusted reference model."
    )
  )
  add_heading(b, 2L, "Analysis questions from the frozen plan")
  plan_map <- data.frame(
    Plan_question_or_constraint = c(
      "Freeze the pseudotime interval before inspecting expression",
      "Avoid a new treated-versus-control expression contrast",
      "Use sample-aware pseudobulk rather than treating cells as independent replicates",
      "Fit a common pseudotime effect with additive nuisance covariates",
      "Annotate the accumulated state with Hallmark, Reactome, and GO BP GSEA",
      "Assess sensitivity to interval and model choices",
      "Evaluate whether ETP adjustment or ETP threshold grouping changes the conclusion"
    ),
    Report_evidence = c(
      "Frozen interval table and density/coverage figures",
      "Primary adjacent-state contrast description and GSEA interpretation",
      "Input audit, pseudobulk metadata, and coverage checks",
      "Model specification and design-rank audit for each model",
      "Primary-state GSEA tables and pathway figures for each collection",
      "Robustness summary and leave-one-mouse-out stability tables",
      "Model-comparison correlations, significant-overlap tables, and rank-shift figures"
    ),
    Interpretation_rule = c(
      "Do not move endpoints after viewing pathway results",
      "Positive NES means state-characterizing, not treatment-upregulated",
      "Mouse/sample is the biological unit; cell counts only support pseudobulk bins",
      "Dose, initial ploidy, and ETP terms are nuisance covariates unless explicitly compared",
      "Prioritize pathways with strong NES, low FDR, coherent activity curves, and stable leading-edge support",
      "Robust findings preserve sign and approximate rank across sensitivity checks",
      "If correlations and FDR-set overlap are high, the biological conclusion is robust to ETP handling"
    ),
    stringsAsFactors = FALSE
  )
  add_html(b, format_table(plan_map, cfg, max_rows = 20L, caption = "Plan-to-report mapping."))
  if (length(plan_lines) > 0L) {
    headings <- grep("^#{1,3} ", plan_lines, value = TRUE)
    headings <- gsub("^#{1,3} +", "", headings)
    if (length(headings) > 0L) {
      add_paragraph(b, paste0("Plan headings detected: ", html_escape(paste(headings, collapse = " | ")), "."))
    }
  }

  add_heading(b, 1L, "Run parameters and input audit")
  add_heading(b, 2L, "Run parameters")
  if (nrow(params) > 0L) add_html(b, format_table(params, cfg, max_rows = 100L, caption = "Workflow command-line parameters."))
  add_heading(b, 2L, "Frozen intervals")
  if (nrow(intervals) > 0L) add_html(b, format_table(intervals, cfg, max_rows = 20L, caption = "Frozen pseudotime intervals."))
  add_heading(b, 2L, "Input and expression audit")
  if (nrow(expression_audit) > 0L) add_html(b, format_table(expression_audit, cfg, max_rows = 5L, caption = "Expression source audit."))
  if (nrow(join_audit) > 0L) add_html(b, format_table(join_audit, cfg, max_rows = 20L, caption = "Cell identity join audit."))
  if (nrow(coverage) > 0L) add_html(b, format_table(coverage, cfg, max_rows = 10L, caption = "Primary coverage check."))
  if (nrow(input_checksums) > 0L) add_html(b, format_table(input_checksums, cfg, max_rows = 30L, caption = "Input checksums."))
  add_heading(b, 2L, "Model set and design-rank audit")
  if (nrow(specs) > 0L) add_html(b, format_table(specs, cfg, max_rows = 20L, caption = "Selected model specifications."))
  if (nrow(design_all) > 0L) add_html(b, format_table(design_all, cfg, max_rows = 20L, caption = "Design-rank audit across all models."))
  if (nrow(etp_cross) > 0L) add_html(b, format_table(etp_cross, cfg, max_rows = 40L, caption = "ETP threshold group cross-tabulation."))
  if (nrow(package_versions) > 0L) {
    add_heading(b, 2L, "Package versions")
    add_html(b, format_table(package_versions, cfg, max_rows = 80L, caption = "Recorded package versions."))
  }

  for (model_id in ordered_models) {
    add_model_section(b, cfg, model_id, specs)
  }

  add_heading(b, 1L, "Model comparison")
  gene_cor <- safe_read_csv(file.path(cfg$results_root, "model_comparison", "gene_t_stat_correlations.csv"))
  nes_cor <- safe_read_csv(file.path(cfg$results_root, "model_comparison", "gsea_NES_correlations.csv"))
  overlap <- safe_read_csv(file.path(cfg$results_root, "model_comparison", "gsea_significant_overlap.csv"))
  rank_shift <- safe_read_csv(file.path(cfg$results_root, "model_comparison", "pathway_rank_shift.csv"))
  add_paragraph(b, summarize_comparison(gene_cor, nes_cor, overlap), "callout")
  if (nrow(gene_cor) > 0L) add_html(b, format_table(gene_cor, cfg, max_rows = 20L, caption = "Gene-level t-statistic correlations versus the reference model."))
  if (nrow(nes_cor) > 0L) add_html(b, format_table(nes_cor, cfg, max_rows = 30L, caption = "GSEA NES correlations versus the reference model."))
  if (nrow(overlap) > 0L) add_html(b, format_table(overlap, cfg, max_rows = 30L, caption = "FDR-significant pathway overlap versus the reference model."))
  if (nrow(rank_shift) > 0L) {
    keep <- intersect(c("collection_label", "pathway_label", "model_id", "reference_rank", "model_rank", "rank_shift", "reference_NES", "model_NES", "reference_padj", "model_padj"), names(rank_shift))
    add_html(b, format_table(utils::head(rank_shift[, keep, drop = FALSE], 40L), cfg, max_rows = 40L, caption = "Top pathway rank shifts, first rows."))
  }
  for (fig in registered_comparison_figures(cfg$results_root)) {
    add_figure(
      b, cfg,
      src = fig$file,
      title = fig$title,
      legend = fig$legend,
      interpretation = fig$interpretation,
      data_paths = fig$data,
      anchor_label = basename(fig$file)
    )
  }

  add_heading(b, 1L, "Interpretation limits")
  add_bullets(
    b,
    c(
      "The primary pathway signature describes the frozen accumulated pseudotime state relative to adjacent states; it is not a treatment differential-expression test.",
      "A positive NES means that the pathway characterizes the accumulated state. It does not by itself imply that gemcitabine directly activates the pathway.",
      "A negative NES means that the pathway is lower in the accumulated state than in neighboring states. It does not by itself imply direct treatment repression.",
      "ETP-adjusted and ETP-threshold models should be used to assess robustness of the state annotation, not to redefine the frozen pseudotime interval.",
      "Pathways with weak FDR support or unstable rank shifts across sensitivity analyses should be treated as descriptive rather than primary conclusions."
    )
  )

  add_heading(b, 1L, "Appendix: output file registry")
  all_files <- list.files(cfg$results_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  all_files <- all_files[!grepl("/report(/|$)", all_files)]
  registry <- data.frame(
    file = display_path(all_files, cfg),
    type = tools::file_ext(all_files),
    stringsAsFactors = FALSE
  )
  add_html(b, format_table(registry, cfg, max_rows = 250L, caption = "Result files used or available to the report."))

  add_html(b, "<p class=\"footer\">End of report.</p>")
  if (isTRUE(b$section_open)) add_html(b, "</section>")

  nav <- build_nav(b$nav)
  html <- paste0(
    "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"/>",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>",
    "<title>04i Pseudotime State Pathway Report</title>",
    "<style>", page_assets$css, "</style>",
    "</head><body><div class=\"report-shell\">",
    "<aside class=\"report-sidebar\"><div class=\"report-sidebar-header\">",
    "<div class=\"report-kicker\">04i report</div>",
    "<div class=\"report-title\">Pseudotime accumulation-state pathways</div>",
    "<div class=\"report-subtitle\">Numbered sections with collapsible navigation</div>",
    "</div><nav class=\"report-nav\">", nav, "</nav></aside>",
    "<main class=\"report-main\">", paste(b$content, collapse = "\n"), "</main>",
    "</div><script>", page_assets$js, "</script></body></html>"
  )
  html <- sanitize_text(html, cfg)
  out <- file.path(cfg$output_dir, "04i_pseudotime_state_pathways_report.html")
  writeLines(html, out, useBytes = TRUE)
  out
}

workflow_display_name <- function(workflow_id) {
  switch(
    workflow_id,
    binning = "Binning workflow",
    non_binning = "Non-binning workflow",
    workflow_id
  )
}

workflow_description <- function(workflow_id) {
  switch(
    workflow_id,
    binning = "This workflow builds sample-by-pseudotime-bin pseudobulk profiles and fits a shared smooth pseudotime model. The minimum retained sample-bin threshold is applied before region coverage and model fitting.",
    non_binning = "This workflow avoids pseudotime binning by aggregating cells directly within the frozen left-neighbor, primary accumulated-state, and right-neighbor intervals. The contrast is primary minus the average of the two neighboring intervals.",
    ""
  )
}

available_figure_specs <- function(specs) {
  specs[vapply(specs, function(x) file.exists(x$file), logical(1L))]
}

workflow_model_figures <- function(model_dir, workflow_id) {
  fig_dir <- file.path(model_dir, "05_figures")
  base <- list(
    list(
      file = file.path(fig_dir, "sample_region_coverage.pdf"),
      title = "Sample coverage in frozen pseudotime regions",
      legend = "Bars show retained cells or retained sample bins for each mouse across the frozen primary and neighboring regions.",
      interpretation = "This figure checks whether the primary contrast is supported by multiple mice, treatment groups, and ploidy or ETP strata.",
      data = file.path(fig_dir, "sample_region_coverage_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_hallmark_gsea.pdf"),
      title = "Primary-state Hallmark GSEA",
      legend = "Bars show Hallmark normalized enrichment scores for the primary adjacent-state contrast. Fill encodes FDR.",
      interpretation = "Hallmark pathways summarize broad biological programs that characterize or are depleted from the frozen accumulated state.",
      data = file.path(fig_dir, "primary_state_hallmark_gsea_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_reactome_gsea.pdf"),
      title = "Primary-state Reactome GSEA",
      legend = "Bars show Reactome normalized enrichment scores for the primary adjacent-state contrast. Fill encodes FDR.",
      interpretation = "Reactome pathways provide mechanistic pathway-level annotation of the accumulated state.",
      data = file.path(fig_dir, "primary_state_reactome_gsea_plot_data.csv")
    ),
    list(
      file = file.path(fig_dir, "primary_state_go_bp_gsea.pdf"),
      title = "Primary-state GO Biological Process GSEA",
      legend = "Bars show GO Biological Process normalized enrichment scores for the primary adjacent-state contrast. Fill encodes FDR.",
      interpretation = "GO Biological Process terms provide ontology-level annotation of state-associated biological processes.",
      data = file.path(fig_dir, "primary_state_go_bp_gsea_plot_data.csv")
    )
  )
  if (identical(workflow_id, "binning")) {
    base <- c(
      list(
        list(
          file = file.path(fig_dir, "pseudotime_density_with_frozen_state.pdf"),
          title = "Pseudotime density with frozen state intervals",
          legend = "Equal-mouse pseudotime densities are shown by treatment group. Shaded regions mark the frozen primary, sensitivity, and negative-control intervals.",
          interpretation = "This is the evidence used to freeze the accumulated-state interval before expression modeling. It should not be re-tuned based on downstream pathway results.",
          data = file.path(fig_dir, c(
            "pseudotime_density_with_frozen_state_mean_density_plot_data.csv",
            "pseudotime_density_with_frozen_state_difference_plot_data.csv",
            "pseudotime_density_with_frozen_state_sample_density_plot_data.csv"
          ))
        )
      ),
      base,
      list(
        list(
          file = file.path(fig_dir, "primary_state_pathway_activity_heatmap.pdf"),
          title = "Primary-state pathway activity heatmap",
          legend = "Columns are pseudotime grid points, rows are top positive and negative pathways, and color is mean standardized fitted activity of leading-edge genes.",
          interpretation = "This heatmap shows whether enriched pathways are localized inside the frozen primary interval or reflect broader trajectory trends.",
          data = file.path(fig_dir, "primary_state_pathway_activity_heatmap_plot_data.csv")
        ),
        list(
          file = file.path(fig_dir, "primary_state_top_pathway_curves.pdf"),
          title = "Top pathway fitted activity curves",
          legend = "Each curve is mean standardized fitted activity for leading-edge genes of a top pathway. The gray band marks the frozen primary interval.",
          interpretation = "Localized peaks inside the frozen interval support state-specific pathway annotation.",
          data = file.path(fig_dir, "primary_state_top_pathway_curves_plot_data.csv")
        ),
        list(
          file = file.path(fig_dir, "primary_state_leading_edge_gene_curves.pdf"),
          title = "Leading-edge gene fitted curves",
          legend = "Each curve is a standardized fitted expression curve for a leading-edge gene from top positive primary-state pathways.",
          interpretation = "This checks whether pathway-level enrichment is supported by coherent gene-level behavior rather than a single outlier gene.",
          data = file.path(fig_dir, "primary_state_leading_edge_gene_curves_plot_data.csv")
        ),
        list(
          file = file.path(fig_dir, "pathway_robustness_summary.pdf"),
          title = "Pathway robustness across sensitivity analyses",
          legend = "Points show NES values across pre-specified sensitivity intervals, modeling variants, and leave-one-mouse-out checks where available.",
          interpretation = "Robust annotations should preserve sign and approximate rank across sensitivity checks.",
          data = file.path(fig_dir, "pathway_robustness_summary_plot_data.csv")
        )
      )
    )
  }
  available_figure_specs(base)
}

add_existing_figure <- function(builder, cfg, fig, title_prefix = "", anchor_prefix = "") {
  if (is.null(fig$file) || !file.exists(fig$file)) return(invisible(FALSE))
  add_figure(
    builder, cfg,
    src = fig$file,
    title = paste0(title_prefix, fig$title),
    legend = fig$legend,
    interpretation = fig$interpretation,
    data_paths = fig$data,
    anchor_label = paste0(anchor_prefix, "-", basename(fig$file))
  )
  invisible(TRUE)
}

add_workflow_model_section <- function(builder, cfg, workflow_id, workflow_root, model_id, specs) {
  model_dir <- file.path(workflow_root, model_id)
  label <- model_display_name(model_id, specs)
  add_heading(builder, 2L, paste0(model_id, ": ", label))
  if (!dir.exists(model_dir)) {
    add_paragraph(builder, paste0("The model output directory is missing: <code>", html_escape(display_path(model_dir, cfg)), "</code>"), "missing")
    return(invisible(NULL))
  }

  model_spec <- if (nrow(specs) > 0L && model_id %in% specs$model_id) specs[match(model_id, specs$model_id), , drop = FALSE] else data.frame()
  model_params <- safe_read_csv(file.path(model_dir, "00_manifest", "model_parameters.csv"))
  design <- safe_read_csv(file.path(model_dir, "01_qc", "model_design_rank_audit.csv"))
  genes <- safe_read_csv(file.path(model_dir, "03_gene_models", "gene_primary_adjacent_state_contrast.csv"))
  gsea <- safe_read_csv(file.path(model_dir, "04_gsea", "all_collections_primary_adjacent_state_gsea.csv"))
  robust <- safe_read_csv(file.path(model_dir, "06_sensitivity", "pathway_robustness_summary.csv"))
  loo <- safe_read_csv(file.path(model_dir, "06_sensitivity", "leave_one_mouse_out_pathway_stability.csv"))

  add_heading(builder, 3L, "Model definition and design audit")
  add_paragraph(
    builder,
    paste0(
      "Model <span class=\"pill\">", html_escape(model_id), "</span> was run inside the ",
      html_escape(workflow_display_name(workflow_id)), ". The primary estimand is the frozen accumulated state relative to its adjacent regions; dose and ploidy or ETP terms are nuisance covariates."
    )
  )
  if (nrow(model_spec) > 0L) add_html(builder, format_table(model_spec, cfg, max_rows = 5L, caption = "Model specification."))
  if (nrow(model_params) > 0L) add_html(builder, format_table(model_params, cfg, max_rows = 20L, caption = "Model-level run parameters."))
  if (nrow(design) > 0L) add_html(builder, format_table(design, cfg, max_rows = 5L, caption = "Design-rank audit."))

  add_heading(builder, 3L, "Primary-state gene and pathway summary")
  add_paragraph(builder, pathway_summary_sentence(gsea), "callout")
  sig <- significant_counts(gsea)
  if (nrow(sig) > 0L) add_html(builder, format_table(sig, cfg, max_rows = 20L, caption = "FDR < 0.05 pathway counts by collection and direction."))
  top_pw <- select_top_pathways(gsea, n_each = 5L)
  if (nrow(top_pw) > 0L) add_html(builder, format_table(top_pw, cfg, max_rows = 40L, caption = "Top primary-state pathways by collection and direction."))
  top_genes <- top_gene_table(genes, n_each = 8L)
  if (nrow(top_genes) > 0L) add_html(builder, format_table(top_genes, cfg, max_rows = 20L, caption = "Top primary-state genes by t statistic and FDR."))

  add_heading(builder, 3L, "Key figures")
  figs <- workflow_model_figures(model_dir, workflow_id)
  gsea_files <- c("primary_state_hallmark_gsea.pdf", "primary_state_reactome_gsea.pdf", "primary_state_go_bp_gsea.pdf")
  gsea_idx <- which(basename(vapply(figs, `[[`, character(1L), "file")) %in% gsea_files)
  for (i in seq_along(figs)) {
    fig <- figs[[i]]
    if (length(gsea_idx) > 0L && i %in% gsea_idx && identical(i, gsea_idx[[1L]])) {
      add_html(builder, "<div class=\"figure-grid figure-grid-3\">")
      for (j in gsea_idx) {
        add_existing_figure(
          builder, cfg, figs[[j]],
          title_prefix = paste0(model_id, ": "),
          anchor_prefix = paste0(workflow_id, "-", model_id)
        )
      }
      add_html(builder, "</div>")
      next
    }
    if (i %in% gsea_idx) next
    add_existing_figure(
      builder, cfg, fig,
      title_prefix = paste0(model_id, ": "),
      anchor_prefix = paste0(workflow_id, "-", model_id)
    )
  }

  if (identical(workflow_id, "binning")) {
    add_heading(builder, 3L, "Sensitivity and leave-one-mouse-out checks")
    if (nrow(robust) > 0L) {
      keep <- intersect(c("collection_label", "pathway_label", "analysis_id", "direction", "NES", "padj", "primary_NES", "primary_padj", "sign_concordant_with_primary"), names(robust))
      add_html(builder, format_table(utils::head(robust[, keep, drop = FALSE], 30L), cfg, max_rows = 30L, caption = "Sensitivity GSEA summary, first rows."))
    } else {
      add_paragraph(builder, "No sensitivity table was available for this model.", "missing")
    }
    if (nrow(loo) > 0L) {
      keep <- intersect(c("collection_label", "pathway_label", "primary_NES", "primary_padj", "n_leave_one_out_runs", "fraction_sign_concordant", "min_NES", "max_NES"), names(loo))
      add_html(builder, format_table(utils::head(loo[, keep, drop = FALSE], 30L), cfg, max_rows = 30L, caption = "Leave-one-mouse-out pathway stability, first rows."))
    } else {
      add_paragraph(builder, "No leave-one-mouse-out stability table was available for this model.", "missing")
    }
  }
  invisible(NULL)
}

summarize_comparison_table <- function(gene_cor, nes_cor, overlap, comparison_label) {
  parts <- character()
  if (nrow(gene_cor) > 0L && "pearson_r" %in% names(gene_cor)) {
    x <- gene_cor
    if ("model_id" %in% names(x)) x <- x[x$model_id != "primary_initial_ploidy", , drop = FALSE]
    if (nrow(x) > 0L) parts <- c(parts, paste0(comparison_label, " gene-level t-statistic concordance: minimum Pearson r = ", format_value(min(x$pearson_r, na.rm = TRUE)), "."))
  }
  if (nrow(nes_cor) > 0L && all(c("collection_label", "pearson_r") %in% names(nes_cor))) {
    x <- nes_cor
    if ("model_id" %in% names(x)) x <- x[x$model_id != "primary_initial_ploidy", , drop = FALSE]
    if (nrow(x) == 0L) x <- nes_cor
    by_coll <- stats::aggregate(pearson_r ~ collection_label, x, min, na.rm = TRUE)
    txt <- paste0(by_coll$collection_label, " min r = ", vapply(by_coll$pearson_r, format_value, character(1L)), collapse = "; ")
    parts <- c(parts, paste0(comparison_label, " pathway NES concordance: ", txt, "."))
  }
  if (nrow(overlap) > 0L && all(c("collection_label", "jaccard") %in% names(overlap))) {
    x <- overlap
    if ("model_id" %in% names(x)) x <- x[x$model_id != "primary_initial_ploidy", , drop = FALSE]
    if (nrow(x) == 0L) x <- overlap
    by_coll <- stats::aggregate(jaccard ~ collection_label, x, min, na.rm = TRUE)
    txt <- paste0(by_coll$collection_label, " min Jaccard = ", vapply(by_coll$jaccard, format_value, character(1L)), collapse = "; ")
    parts <- c(parts, paste0(comparison_label, " FDR-significant pathway overlap: ", txt, "."))
  }
  if (length(parts) == 0L) "Comparison summary tables were not available." else paste(parts, collapse = " ")
}

add_internal_comparison_section <- function(builder, cfg, workflow_id, workflow_root) {
  add_heading(builder, 2L, "Internal model comparison")
  comp_root <- file.path(workflow_root, "model_comparison")
  gene_cor <- safe_read_csv(file.path(comp_root, "gene_t_stat_correlations.csv"))
  nes_cor <- safe_read_csv(file.path(comp_root, "gsea_NES_correlations.csv"))
  overlap <- safe_read_csv(file.path(comp_root, "gsea_significant_overlap.csv"))
  rank_shift <- safe_read_csv(file.path(comp_root, "pathway_rank_shift.csv"))
  add_paragraph(builder, summarize_comparison_table(gene_cor, nes_cor, overlap, workflow_display_name(workflow_id)), "callout")
  if (nrow(gene_cor) > 0L) add_html(builder, format_table(gene_cor, cfg, max_rows = 20L, caption = "Gene-level t-statistic correlations versus the primary initial-ploidy model."))
  if (nrow(nes_cor) > 0L) add_html(builder, format_table(nes_cor, cfg, max_rows = 30L, caption = "GSEA NES correlations versus the primary initial-ploidy model."))
  if (nrow(overlap) > 0L) add_html(builder, format_table(overlap, cfg, max_rows = 30L, caption = "FDR-significant pathway overlap versus the primary initial-ploidy model."))
  if (nrow(rank_shift) > 0L) {
    keep <- intersect(c("collection_label", "pathway_label", "model_id", "reference_rank", "model_rank", "rank_shift", "reference_NES", "model_NES", "reference_padj", "model_padj"), names(rank_shift))
    add_html(builder, format_table(utils::head(rank_shift[, keep, drop = FALSE], 40L), cfg, max_rows = 40L, caption = "Top pathway rank shifts, first rows."))
  }
  for (fig in registered_comparison_figures(workflow_root)) {
    add_existing_figure(
      builder, cfg, fig,
      title_prefix = paste0(workflow_display_name(workflow_id), ": "),
      anchor_prefix = paste0(workflow_id, "-comparison")
    )
  }
}

add_workflow_section <- function(builder, cfg, workflow_id, specs, ordered_models) {
  workflow_root <- file.path(cfg$results_root, workflow_id)
  add_heading(builder, 1L, workflow_display_name(workflow_id))
  add_paragraph(builder, workflow_description(workflow_id), "callout")
  if (!dir.exists(workflow_root)) {
    add_paragraph(builder, paste0("Workflow directory is missing: <code>", html_escape(display_path(workflow_root, cfg)), "</code>"), "missing")
    return(invisible(NULL))
  }

  add_heading(builder, 2L, "Workflow coverage and task log")
  coverage_check <- safe_read_csv(file.path(workflow_root, "01_qc", "primary_coverage_check.csv"))
  coverage_file <- if (identical(workflow_id, "non_binning")) "interval_region_coverage.csv" else "sample_region_coverage.csv"
  coverage <- safe_read_csv(file.path(workflow_root, "01_qc", coverage_file))
  task_log <- safe_read_csv(file.path(workflow_root, "00_manifest", "parallel_task_log.csv"))
  design_all <- safe_read_csv(file.path(workflow_root, "model_comparison", "model_design_rank_audit_all_models.csv"))
  if (nrow(coverage_check) > 0L) add_html(builder, format_table(coverage_check, cfg, max_rows = 10L, caption = "Primary coverage check."))
  if (nrow(coverage) > 0L) add_html(builder, format_table(coverage, cfg, max_rows = 40L, caption = "Region coverage, first rows."))
  if (nrow(task_log) > 0L) add_html(builder, format_table(task_log, cfg, max_rows = 20L, caption = "Parallel model task log."))
  if (nrow(design_all) > 0L) add_html(builder, format_table(design_all, cfg, max_rows = 20L, caption = "Design-rank audit across workflow models."))

  for (model_id in ordered_models) {
    add_workflow_model_section(builder, cfg, workflow_id, workflow_root, model_id, specs)
  }
  add_internal_comparison_section(builder, cfg, workflow_id, workflow_root)
  invisible(NULL)
}

registered_cross_workflow_figures <- function(root) {
  fig_dir <- file.path(root, "cross_workflow_comparison", "figures")
  list(
    list(
      file = file.path(fig_dir, "gsea_NES_correlation_binning_vs_non_binning.pdf"),
      title = "Primary-state GSEA NES: binning versus non-binning",
      legend = "Each point is a pathway NES value for the same model and collection, comparing the binned and non-binned workflows.",
      interpretation = "High concordance supports treating the no-binning workflow as a complementary sensitivity analysis rather than a competing definition of the accumulated state.",
      data = file.path(root, "cross_workflow_comparison", "gsea_NES_correlations_binning_vs_non_binning.csv")
    ),
    list(
      file = file.path(fig_dir, "top_pathway_rank_shift_binning_vs_non_binning.pdf"),
      title = "Largest pathway rank shifts: binning versus non-binning",
      legend = "Points show pathway rank movement when the same model is fit with interval-level pseudobulk instead of pseudotime bins.",
      interpretation = "Large rank shifts identify annotations that are sensitive to the binning decision; stable top pathways should show small or directionally coherent shifts.",
      data = file.path(root, "cross_workflow_comparison", "pathway_rank_shift_binning_vs_non_binning.csv")
    )
  )
}

add_cross_workflow_section <- function(builder, cfg) {
  add_heading(builder, 1L, "Cross-workflow comparison")
  comp_root <- file.path(cfg$results_root, "cross_workflow_comparison")
  gene_cor <- safe_read_csv(file.path(comp_root, "gene_t_stat_correlations_binning_vs_non_binning.csv"))
  nes_cor <- safe_read_csv(file.path(comp_root, "gsea_NES_correlations_binning_vs_non_binning.csv"))
  overlap <- safe_read_csv(file.path(comp_root, "gsea_significant_overlap_binning_vs_non_binning.csv"))
  rank_shift <- safe_read_csv(file.path(comp_root, "pathway_rank_shift_binning_vs_non_binning.csv"))
  add_paragraph(builder, summarize_comparison_table(gene_cor, nes_cor, overlap, "Binning versus non-binning"), "callout")
  if (nrow(gene_cor) > 0L) add_html(builder, format_table(gene_cor, cfg, max_rows = 20L, caption = "Gene-level t-statistic correlations: binning versus non-binning."))
  if (nrow(nes_cor) > 0L) add_html(builder, format_table(nes_cor, cfg, max_rows = 40L, caption = "GSEA NES correlations: binning versus non-binning."))
  if (nrow(overlap) > 0L) add_html(builder, format_table(overlap, cfg, max_rows = 40L, caption = "FDR-significant pathway overlap: binning versus non-binning."))
  if (nrow(rank_shift) > 0L) {
    keep <- intersect(c("model_id", "collection_label", "pathway_label", "binning_rank", "non_binning_rank", "rank_shift", "NES_binning", "NES_non_binning", "NES_delta", "padj_binning", "padj_non_binning"), names(rank_shift))
    add_html(builder, format_table(utils::head(rank_shift[, keep, drop = FALSE], 60L), cfg, max_rows = 60L, caption = "Largest pathway rank shifts, first rows."))
  }
  for (fig in registered_cross_workflow_figures(cfg$results_root)) {
    add_existing_figure(
      builder, cfg, fig,
      title_prefix = "",
      anchor_prefix = "cross-workflow"
    )
  }
}

minimum_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

format_metric <- function(x, digits = 3L) {
  if (!is.finite(x)) return("not available")
  formatC(x, format = "f", digits = digits)
}

conclusion_workflow_metrics <- function(cfg, workflow_id) {
  comp_root <- file.path(cfg$results_root, workflow_id, "model_comparison")
  coverage <- safe_read_csv(file.path(cfg$results_root, workflow_id, "01_qc", "primary_coverage_check.csv"))
  gene_cor <- safe_read_csv(file.path(comp_root, "gene_t_stat_correlations.csv"))
  nes_cor <- safe_read_csv(file.path(comp_root, "gsea_NES_correlations.csv"))
  overlap <- safe_read_csv(file.path(comp_root, "gsea_significant_overlap.csv"))
  if (nrow(gene_cor) > 0L && "model_id" %in% names(gene_cor)) {
    gene_cor <- gene_cor[gene_cor$model_id != "primary_initial_ploidy", , drop = FALSE]
  }
  if (nrow(nes_cor) > 0L && "model_id" %in% names(nes_cor)) {
    nes_cor <- nes_cor[nes_cor$model_id != "primary_initial_ploidy", , drop = FALSE]
  }
  if (nrow(overlap) > 0L && "model_id" %in% names(overlap)) {
    overlap <- overlap[overlap$model_id != "primary_initial_ploidy", , drop = FALSE]
  }
  data.frame(
    workflow = workflow_display_name(workflow_id),
    contributing_mice = if (nrow(coverage) > 0L && "n_contributing_mice" %in% names(coverage)) coverage$n_contributing_mice[[1L]] else NA_integer_,
    min_gene_t_pearson_across_ETP_models = minimum_or_na(gene_cor$pearson_r),
    min_pathway_NES_pearson_across_ETP_models = minimum_or_na(nes_cor$pearson_r),
    min_FDR_pathway_Jaccard_across_ETP_models = minimum_or_na(overlap$jaccard),
    stringsAsFactors = FALSE
  )
}

conclusion_cross_metrics <- function(cfg) {
  comp_root <- file.path(cfg$results_root, "cross_workflow_comparison")
  gene_cor <- safe_read_csv(file.path(comp_root, "gene_t_stat_correlations_binning_vs_non_binning.csv"))
  nes_cor <- safe_read_csv(file.path(comp_root, "gsea_NES_correlations_binning_vs_non_binning.csv"))
  overlap <- safe_read_csv(file.path(comp_root, "gsea_significant_overlap_binning_vs_non_binning.csv"))
  data.frame(
    comparison = "Binning versus non-binning",
    min_gene_t_pearson = minimum_or_na(gene_cor$pearson_r),
    min_pathway_NES_pearson = minimum_or_na(nes_cor$pearson_r),
    min_FDR_pathway_Jaccard = minimum_or_na(overlap$jaccard),
    stringsAsFactors = FALSE
  )
}

conclusion_direction_metrics <- function(cfg) {
  binning <- safe_read_csv(file.path(cfg$results_root, "binning", "primary_initial_ploidy", "04_gsea", "all_collections_primary_adjacent_state_gsea.csv"))
  non_binning <- safe_read_csv(file.path(cfg$results_root, "non_binning", "primary_initial_ploidy", "04_gsea", "all_collections_primary_adjacent_state_gsea.csv"))
  needed <- c("collection_label", "pathway", "pathway_label", "NES", "padj")
  if (!all(needed %in% names(binning)) || !all(needed %in% names(non_binning))) return(data.frame())
  b <- binning[, needed, drop = FALSE]
  n <- non_binning[, c("collection_label", "pathway", "NES", "padj"), drop = FALSE]
  names(b)[names(b) == "NES"] <- "NES_binning"
  names(b)[names(b) == "padj"] <- "padj_binning"
  names(n)[names(n) == "NES"] <- "NES_non_binning"
  names(n)[names(n) == "padj"] <- "padj_non_binning"
  merged <- merge(b, n, by = c("collection_label", "pathway"), all = FALSE, sort = FALSE)
  merged$binning_significant <- is.finite(merged$padj_binning) & merged$padj_binning < 0.05
  merged$non_binning_significant <- is.finite(merged$padj_non_binning) & merged$padj_non_binning < 0.05
  merged$both_significant <- merged$binning_significant & merged$non_binning_significant
  merged$same_direction <- sign(merged$NES_binning) == sign(merged$NES_non_binning)
  rows <- lapply(split(merged, merged$collection_label), function(df) {
    data.frame(
      collection_label = df$collection_label[[1L]],
      common_tested_pathways = nrow(df),
      binning_FDR_lt_0_05 = sum(df$binning_significant, na.rm = TRUE),
      non_binning_FDR_lt_0_05 = sum(df$non_binning_significant, na.rm = TRUE),
      both_FDR_lt_0_05 = sum(df$both_significant, na.rm = TRUE),
      both_significant_same_direction = sum(df$both_significant & df$same_direction, na.rm = TRUE),
      both_significant_opposite_direction = sum(df$both_significant & !df$same_direction, na.rm = TRUE),
      NES_pearson = minimum_or_na(stats::cor(df$NES_binning, df$NES_non_binning, use = "pairwise.complete.obs")),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(match(out$collection_label, c("hallmark", "reactome", "go_bp"))), , drop = FALSE]
}

add_conclusion_section <- function(builder, cfg) {
  workflow_metrics <- rbind(
    conclusion_workflow_metrics(cfg, "binning"),
    conclusion_workflow_metrics(cfg, "non_binning")
  )
  cross_metrics <- conclusion_cross_metrics(cfg)
  direction_metrics <- conclusion_direction_metrics(cfg)
  add_heading(builder, 1L, "Conclusion")
  add_paragraph(
    builder,
    "The accumulated CellCycle pseudotime state is most consistently characterized by depletion of proliferation, MYC/E2F/G2M, DNA replication and repair, mitochondrial oxidative phosphorylation, mitochondrial gene-expression, and translation-related programs. This is the high-confidence biological conclusion because it is preserved across the initial-ploidy model, continuous ETP adjustment, ETP threshold models, and both binned and non-binned pseudobulk workflows.",
    "callout"
  )
  add_bullets(
    builder,
    c(
      paste0(
        "ETP handling does not drive the result. Within the binning workflow, the minimum gene-level Pearson correlation across ETP models is ",
        format_metric(workflow_metrics$min_gene_t_pearson_across_ETP_models[workflow_metrics$workflow == "Binning workflow"]),
        " and the minimum pathway-NES Pearson correlation is ",
        format_metric(workflow_metrics$min_pathway_NES_pearson_across_ETP_models[workflow_metrics$workflow == "Binning workflow"]),
        ". Within the non-binning workflow, the corresponding minima are ",
        format_metric(workflow_metrics$min_gene_t_pearson_across_ETP_models[workflow_metrics$workflow == "Non-binning workflow"]),
        " and ",
        format_metric(workflow_metrics$min_pathway_NES_pearson_across_ETP_models[workflow_metrics$workflow == "Non-binning workflow"]),
        "."
      ),
      paste0(
        "The main difference is the workflow, not the covariate choice. Binning is more sensitive and calls a broader set of significant pathways, while non-binning is more conservative and avoids choosing a pseudotime bin size. Cross-workflow pathway-NES concordance is moderate rather than near-identical, with minimum Pearson correlation ",
        format_metric(cross_metrics$min_pathway_NES_pearson[[1L]]),
        "."
      ),
      "The shared significant signal is directionally stable. Among pathways significant in both primary workflows, no collection shows opposite-direction enrichment for the shared significant pathways.",
      "Workflow-specific positive pathways, including hypoxia/GPCR/sensory-like signals, should be treated as secondary or exploratory. The primary interpretation should emphasize the shared negative cell-cycle, replication-repair, mitochondrial, and translation programs."
    )
  )
  add_html(builder, format_table(workflow_metrics, cfg, max_rows = 10L, caption = "Internal robustness of ETP handling within each workflow."))
  add_html(builder, format_table(cross_metrics, cfg, max_rows = 5L, caption = "Cross-workflow concordance summary."))
  if (nrow(direction_metrics) > 0L) {
    add_html(builder, format_table(direction_metrics, cfg, max_rows = 10L, caption = "Direction check for pathways tested in both primary workflows."))
  }
}

build_report <- function(cfg) {
  cfg$results_root <- normalize_existing(cfg$results_root)
  cfg$output_dir <- normalizePath(cfg$output_dir, winslash = "/", mustWork = FALSE)
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(cfg$overwrite)) {
    unlink(file.path(cfg$output_dir, "assets"), recursive = TRUE, force = TRUE)
  }
  page_assets <- write_css_js(cfg)

  b <- section_builder(cfg)
  params <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "analysis_parameters.csv"))
  input_checksums <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "input_checksums.csv"))
  package_versions <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "package_versions.csv"))
  intervals <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "frozen_interval_definition.csv"))
  specs <- safe_read_csv(file.path(cfg$results_root, "00_manifest", "selected_model_specifications.csv"))
  expression_audit <- safe_read_csv(file.path(cfg$results_root, "shared_qc", "expression_source_audit.csv"))
  join_audit <- safe_read_csv(file.path(cfg$results_root, "shared_qc", "cell_id_join_audit.csv"))
  etp_cross <- safe_read_csv(file.path(cfg$results_root, "shared_qc", "etp_group_cross_tab.csv"))
  etp_covariates <- safe_read_csv(file.path(cfg$results_root, "shared_qc", "etp_covariate_sample_table.csv"))
  plan_lines <- safe_read_lines(cfg$plan_md)

  if (nrow(specs) == 0L) {
    specs <- data.frame(model_id = model_order, model_label = model_order, stringsAsFactors = FALSE)
  }
  ordered_models <- model_order[model_order %in% specs$model_id]
  ordered_models <- c(ordered_models, setdiff(specs$model_id, ordered_models))

  add_html(
    b,
    "<section class=\"report-card\" id=\"report-title\">",
    "<h1>Pseudotime accumulation-state pathway analysis report</h1>",
    "<p class=\"muted\">Generated ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    ". Reported result root: <code>", html_escape(display_path(cfg$results_root, cfg)), "</code>.</p>",
    "<p>This self-contained HTML report summarizes the 04i binning and non-binning pathway workflows. It is generated from local result files, but all displayed absolute paths are mapped to the HPC namespace.</p>",
    "</section>"
  )

  add_conclusion_section(b, cfg)

  add_heading(b, 1L, "Executive summary")
  add_bullets(
    b,
    c(
      "The primary biological question is which expression programs characterize the frozen CellCycle pseudotime accumulated state relative to adjacent trajectory regions.",
      "The binning workflow is the main smooth pseudotime analysis with min_cells_per_sample_bin set to 5.",
      "The non-binning workflow is a complementary sensitivity analysis that aggregates cells directly inside the frozen intervals and avoids choosing a pseudotime bin size.",
      "Five matched covariate specifications are run in both workflows: initial ploidy, continuous ETP, and three thresholded ETP group models.",
      "Cross-workflow comparison asks whether the biological pathway annotation is stable to the binning decision."
    )
  )

  add_heading(b, 2L, "Analysis questions from the frozen plan")
  plan_map <- data.frame(
    Plan_question_or_constraint = c(
      "Freeze the pseudotime interval before inspecting expression",
      "Avoid a new treated-versus-control expression contrast",
      "Use sample-aware pseudobulk rather than treating cells as independent replicates",
      "Fit nuisance-adjusted models with initial ploidy and ETP alternatives",
      "Annotate the accumulated state with Hallmark, Reactome, and GO BP GSEA",
      "Assess sensitivity to interval, model, and binning choices"
    ),
    Report_evidence = c(
      "Frozen interval table and density figure in the binning workflow",
      "Primary adjacent-state contrast description and pathway interpretation",
      "Pseudobulk metadata, coverage checks, and model task logs",
      "Five model sections in each workflow plus internal model comparison",
      "Three collection-specific GSEA figures shown in one row for each model",
      "Binning internal comparison, non-binning internal comparison, and cross-workflow comparison"
    ),
    Interpretation_rule = c(
      "Do not move endpoints after viewing pathway results",
      "Positive NES means state-characterizing, not treatment-upregulated",
      "Mouse/sample is the biological unit; cells support pseudobulk construction",
      "ETP models evaluate robustness rather than redefining the frozen state",
      "Prioritize pathways with strong NES, low FDR, coherent fitted curves, and stable rank",
      "Treat binning and non-binning concordance as evidence that conclusions are not bin-size artifacts"
    ),
    stringsAsFactors = FALSE
  )
  add_html(b, format_table(plan_map, cfg, max_rows = 20L, caption = "Plan-to-report mapping."))
  if (length(plan_lines) > 0L) {
    headings <- grep("^#{1,3} ", plan_lines, value = TRUE)
    headings <- gsub("^#{1,3} +", "", headings)
    if (length(headings) > 0L) add_paragraph(b, paste0("Plan headings detected: ", html_escape(paste(headings, collapse = " | ")), "."))
  }

  add_heading(b, 1L, "Run parameters and shared input audit")
  add_heading(b, 2L, "Run parameters")
  if (nrow(params) > 0L) add_html(b, format_table(params, cfg, max_rows = 100L, caption = "Workflow command-line parameters."))
  add_heading(b, 2L, "Frozen intervals")
  if (nrow(intervals) > 0L) add_html(b, format_table(intervals, cfg, max_rows = 20L, caption = "Frozen pseudotime intervals."))
  add_heading(b, 2L, "Shared input and expression audit")
  if (nrow(expression_audit) > 0L) add_html(b, format_table(expression_audit, cfg, max_rows = 5L, caption = "Expression source audit."))
  if (nrow(join_audit) > 0L) add_html(b, format_table(join_audit, cfg, max_rows = 20L, caption = "Cell identity join audit."))
  if (nrow(input_checksums) > 0L) add_html(b, format_table(input_checksums, cfg, max_rows = 30L, caption = "Input checksums."))
  add_heading(b, 2L, "Model set and ETP audit")
  if (nrow(specs) > 0L) add_html(b, format_table(specs, cfg, max_rows = 20L, caption = "Selected model specifications."))
  if (nrow(etp_covariates) > 0L) add_html(b, format_table(etp_covariates, cfg, max_rows = 30L, caption = "ETP covariate sample table, first rows."))
  if (nrow(etp_cross) > 0L) add_html(b, format_table(etp_cross, cfg, max_rows = 40L, caption = "ETP threshold group cross-tabulation."))
  if (nrow(package_versions) > 0L) {
    add_heading(b, 2L, "Package versions")
    add_html(b, format_table(package_versions, cfg, max_rows = 80L, caption = "Recorded package versions."))
  }

  add_workflow_section(b, cfg, "binning", specs, ordered_models)
  add_workflow_section(b, cfg, "non_binning", specs, ordered_models)
  add_cross_workflow_section(b, cfg)

  add_heading(b, 1L, "Interpretation limits")
  add_bullets(
    b,
    c(
      "The primary pathway signature describes the frozen accumulated pseudotime state relative to adjacent states; it is not a direct treatment differential-expression test.",
      "A positive NES means the pathway characterizes the accumulated state. It does not by itself imply gemcitabine directly activates that pathway.",
      "A negative NES means the pathway is lower in the accumulated state than in neighboring states. It does not by itself imply direct treatment repression.",
      "ETP-adjusted and ETP-threshold models evaluate robustness of the state annotation; they should not be used to redefine the frozen pseudotime interval.",
      "The non-binning workflow is best interpreted as complementary sensitivity evidence for the binning workflow because it changes the pseudobulk unit and model shape."
    )
  )

  add_heading(b, 1L, "Appendix: output file registry")
  all_files <- list.files(cfg$results_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  all_files <- all_files[!grepl("/report(/|$)", all_files)]
  registry <- data.frame(file = display_path(all_files, cfg), type = tools::file_ext(all_files), stringsAsFactors = FALSE)
  add_html(b, format_table(registry, cfg, max_rows = 300L, caption = "Result files used or available to the report."))

  add_html(b, "<p class=\"footer\">End of report.</p>")
  if (isTRUE(b$section_open)) add_html(b, "</section>")

  nav <- build_nav(b$nav)
  html <- paste0(
    "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"/>",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>",
    "<title>04i Pseudotime State Pathway Report</title>",
    "<style>", page_assets$css, "</style>",
    "</head><body><div class=\"report-shell\">",
    "<aside class=\"report-sidebar\"><div class=\"report-sidebar-header\">",
    "<div class=\"report-kicker\">04i report</div>",
    "<div class=\"report-title\">Pseudotime accumulation-state pathways</div>",
    "<div class=\"report-subtitle\">Binning, non-binning, and cross-workflow comparison</div>",
    "</div><nav class=\"report-nav\">", nav, "</nav></aside>",
    "<main class=\"report-main\">", paste(b$content, collapse = "\n"), "</main>",
    "</div><script>", page_assets$js, "</script></body></html>"
  )
  html <- sanitize_text(html, cfg)
  out <- file.path(cfg$output_dir, "04i_pseudotime_state_pathways_report.html")
  writeLines(html, out, useBytes = TRUE)
  out
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!dir.exists(cfg$results_root)) stop("Results root does not exist: ", cfg$results_root, call. = FALSE)
  if (dir.exists(cfg$output_dir) && !isTRUE(cfg$overwrite)) {
    stop("Output directory already exists and --overwrite=FALSE: ", cfg$output_dir, call. = FALSE)
  }
  report <- build_report(cfg)
  message("Wrote report: ", report)
}

if (identical(environment(), globalenv())) {
  main()
}
