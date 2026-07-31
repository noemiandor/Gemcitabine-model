#!/usr/bin/env Rscript

test_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(test_file_arg)) stop("Cannot resolve test path", call. = FALSE)
test_path <- normalizePath(
  sub("^--file=", "", test_file_arg[[1L]]),
  mustWork = TRUE
)
repo_root <- normalizePath(
  file.path(dirname(test_path), "..", "..", "..", ".."),
  mustWork = TRUE
)

for (package in c("patchwork", "pheatmap")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("SI7 heatmap test requires package: ", package, call. = FALSE)
  }
}

# Evaluate the shared production helpers without running the top-level
# renderer.
helper_path <- file.path(
  repo_root,
  "Code", "in-vivo", "SI_figures", "shared_context_panels.R"
)
helper_environment <- new.env(parent = globalenv())
sys.source(helper_path, envir = helper_environment)

read_si7_matrix <- function(filename) {
  data <- utils::read.delim(
    file.path(repo_root, "Data", "in-vivo", "SIfigures", filename),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  matrix_data <- as.matrix(data[, -1L, drop = FALSE])
  storage.mode(matrix_data) <- "double"
  rownames(matrix_data) <- data[[1L]]
  matrix_data
}

si7_panels <- list(
  A = list(
    matrix = read_si7_matrix(
      "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv"
    ),
    diverging = FALSE
  ),
  B = list(
    matrix = read_si7_matrix(
      "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv"
    ),
    diverging = TRUE
  )
)

for (panel_name in names(si7_panels)) {
  panel <- si7_panels[[panel_name]]
  heatmap <- helper_environment$shared_context_build_heatmap(
    panel$matrix,
    paste("SI7", panel_name),
    panel$diverging,
    panel_name
  )
  repeated_heatmap <- helper_environment$shared_context_build_heatmap(
    panel$matrix,
    paste("SI7", panel_name),
    panel$diverging,
    panel_name
  )
  stopifnot(
    inherits(heatmap$tree_row, "hclust"),
    inherits(heatmap$tree_col, "hclust"),
    length(heatmap$tree_row$order) == nrow(panel$matrix),
    length(heatmap$tree_col$order) == ncol(panel$matrix),
    setequal(heatmap$tree_row$labels, rownames(panel$matrix)),
    setequal(heatmap$tree_col$labels, colnames(panel$matrix)),
    identical(heatmap$tree_row$order, repeated_heatmap$tree_row$order),
    identical(heatmap$tree_col$order, repeated_heatmap$tree_col$order)
  )

  layout_names <- heatmap$gtable$layout$name
  stopifnot(
    sum(layout_names == "row_tree") == 1L,
    sum(layout_names == "col_tree") == 1L
  )
  column_tree <- heatmap$gtable$grobs[[which(layout_names == "col_tree")]]
  if (!inherits(column_tree, "polyline")) {
    stop("SI7", panel_name, " column dendrogram is not visible", call. = FALSE)
  }

  plot <- helper_environment$shared_context_heatmap_plot(
    panel$matrix,
    paste("SI7", panel_name),
    panel$diverging,
    panel_name
  )
  heatmap_gtable <- attr(plot, "grobs")$full
  if (is.null(heatmap_gtable)) {
    stop("SI7", panel_name, " did not retain its heatmap gtable", call. = FALSE)
  }
  dendrograms <- intersect(
    heatmap_gtable$layout$name,
    c("row_tree", "col_tree")
  )
  if (!setequal(dendrograms, c("row_tree", "col_tree"))) {
    stop(
      "SI7", panel_name,
      " must render both row and column dendrograms; observed: ",
      paste(dendrograms, collapse = ", "),
      call. = FALSE
    )
  }
}

si7_built <- helper_environment$shared_context_build_si7_heatmap_panels(
  ora_matrix = si7_panels$A$matrix,
  gsea_matrix = si7_panels$B$matrix
)
stopifnot(
  identical(names(si7_built), "plots"),
  identical(names(si7_built$plots), c("ora", "gsea")),
  all(vapply(si7_built$plots, inherits, logical(1L), what = "ggplot"))
)

# Empty display tags are supported for reuse in a compositor that supplies its
# own uniform A-K labels. They must not leave an embedded title prefix.
untagged <- helper_environment$shared_context_build_heatmap(
  si7_panels$B$matrix,
  "Cluster Hallmark GSEA NES",
  TRUE,
  ""
)
main_grob <- untagged$gtable$grobs[[
  which(untagged$gtable$layout$name == "main")
]]
stopifnot(
  identical(main_grob$label, "Cluster Hallmark GSEA NES"),
  inherits(untagged$tree_row, "hclust"),
  inherits(untagged$tree_col, "hclust")
)

generator_text <- paste(readLines(file.path(
  repo_root,
  "Code", "in-vivo", "SI_figures", "generate_supplementary_figures.R"
), warn = FALSE), collapse = "\n")
stopifnot(
  grepl("sys.source(shared_context_helper_path", generator_text, fixed = TRUE),
  grepl("shared_context_build_si7_heatmap_panels(", generator_text, fixed = TRUE),
  !grepl("build_heatmap <- function", generator_text, fixed = TRUE),
  !grepl("heatmap_plot <- function", generator_text, fixed = TRUE)
)

message("SI7 row/column heatmap clustering tests passed.")
