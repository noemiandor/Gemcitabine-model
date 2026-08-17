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

for (package in c("ggplot2", "yaml")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Shared context-panel test requires package: ", package, call. = FALSE)
  }
}

helper_environment <- new.env(parent = globalenv())
sys.source(
  file.path(
    repo_root,
    "Code", "in-vivo", "SI_figures", "normalized_composition.R"
  ),
  envir = helper_environment
)
sys.source(
  file.path(
    repo_root,
    "Code", "in-vivo", "SI_figures", "shared_context_panels.R"
  ),
  envir = helper_environment
)

config <- yaml::read_yaml(file.path(
  repo_root,
  "Code", "in-vivo", "figure7", "figure7_config.yaml"
))$si_figures
cluster_levels <- as.character(unlist(config$cluster_order))
ploidy_levels <- as.character(unlist(config$ploidy_levels))
dose_levels <- as.character(unlist(config$dose_levels))
plot_seed <- as.integer(config$plot_shuffle_seed)

cells <- utils::read.csv(
  file.path(
    repo_root,
    "Data", "in-vivo", "SIfigures", "si_figures_cell_metadata.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cluster_key <- utils::read.delim(
  file.path(
    repo_root,
    "Data", "in-vivo", "SIfigures", "si_figures_cluster_key.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
data <- data.frame(
  cell = cells$cell_id,
  UMAP_1 = as.numeric(cells$UMAP_1),
  UMAP_2 = as.numeric(cells$UMAP_2),
  mouse = cells$sample_id,
  cluster = factor(cells$cluster_id, levels = cluster_levels),
  initial_ploidy = factor(cells$initial_ploidy, levels = ploidy_levels),
  dose = factor(cells$dose, levels = dose_levels),
  context = factor(cells$context, levels = c("Tumor", "CellLine")),
  included = tolower(as.character(cells$included_in_si_figures)) %in%
    c("true", "t", "1"),
  stringsAsFactors = FALSE
)
n_tumor <- sum(data$included)
point_size <- helper_environment$shared_context_publication_umap_point_size(
  nrow(data)
)
tumor_point_size <-
  helper_environment$shared_context_publication_umap_point_size(n_tumor)
faceted_point_size <-
  helper_environment$shared_context_publication_umap_point_size(
    max(table(data$mouse[data$included])),
    faceted = TRUE
  )
stopifnot(
  point_size >= 0.50,
  tumor_point_size >= 0.62,
  faceted_point_size >= 0.85,
  faceted_point_size > tumor_point_size
)

midpoint_test_data <- data.frame(
  UMAP_1 = c(0, 1, 2),
  UMAP_2 = c(0, 1, 0),
  score = c(-1, 0, 1)
)
midpoint_test_plot <- helper_environment$shared_context_make_umap_continuous(
  data = midpoint_test_data,
  field = "score",
  title = "Midpoint test",
  legend_title = "Score",
  tag = "",
  point_size = 0.70,
  plot_seed = plot_seed,
  limits = c(-1, 1),
  diverging = TRUE,
  mid_color = "#A6A6A6",
  shuffle_key = "D"
)
stopifnot(
  identical(
    toupper(midpoint_test_plot$scales$get_scales("colour")$map(0)),
    "#A6A6A6"
  )
)
cluster_colors <- stats::setNames(
  as.character(cluster_key$color),
  as.character(cluster_key$cluster_id)
)
ploidy_colors <- c("2N" = "#4C78A8", "4N" = "#E45756")
context_colors <- c("Tumor" = "#4C78A8", "CellLine" = "#F2CF5B")

build_panels <- function(tags) {
  helper_environment$shared_context_build_si4_panels(
    data = data,
    cluster_levels = cluster_levels,
    cluster_colors = cluster_colors,
    ploidy_colors = ploidy_colors,
    context_colors = context_colors,
    point_size = point_size,
    plot_seed = plot_seed,
    tags = tags,
    composition_builder = helper_environment$make_normalized_composition_plot
  )
}

supplement_tags <- c(
  cluster = "A",
  initial_ploidy = "B",
  context = "C",
  composition = "E"
)
supplement <- build_panels(supplement_tags)
repeated <- build_panels(supplement_tags)
untagged <- build_panels(stats::setNames(rep("", 4L), names(supplement_tags)))

stopifnot(
  identical(names(supplement), c("plots", "composition_result")),
  identical(
    names(supplement$plots),
    c("cluster", "initial_ploidy", "context", "composition")
  ),
  all(vapply(supplement$plots, inherits, logical(1L), what = "ggplot")),
  identical(
    vapply(
      supplement$plots,
      function(plot) plot$labels$tag,
      character(1L)
    ),
    supplement_tags
  )
)

# Repeated builds and retagged copies retain the exact source-panel draw order.
for (panel_name in c("cluster", "initial_ploidy", "context")) {
  expected_order <- supplement$plots[[panel_name]]$data$cell
  stopifnot(
    identical(expected_order, repeated$plots[[panel_name]]$data$cell),
    identical(expected_order, untagged$plots[[panel_name]]$data$cell),
    is.null(untagged$plots[[panel_name]]$labels$tag)
  )
}
stopifnot(is.null(untagged$plots$composition$labels$tag))
stopifnot(
  sum(vapply(
    supplement$plots$cluster$layers,
    function(layer) identical(class(layer$geom)[[1L]], "GeomLabel"),
    logical(1L)
  )) == 1L,
  !any(vapply(
    supplement$plots$cluster$layers,
    function(layer) inherits(layer$geom, "GeomLabelRepel"),
    logical(1L)
  ))
)

# SI4E remains an equal-sample, initial-ploidy-stratified exact permutation
# analysis. The reviewed cache has 16 Tumor and 2 CellLine samples and 81 exact
# assignments; only positive context enrichments receive stars.
composition <- supplement$composition_result
enriched <- composition$tests[composition$tests$enriched, , drop = FALSE]
stopifnot(
  nrow(composition$tests) == 18L,
  all(composition$tests$exact_permutations == 81L),
  all(composition$tests$strata == "initial_ploidy"),
  all(composition$tests$significance[composition$tests$difference <= 0] == ""),
  setequal(
    paste(enriched$group_value, enriched$cluster, sep = "/"),
    c(
      "Tumor/6", "Tumor/10", "Tumor/13",
      "CellLine/0", "CellLine/8", "CellLine/14"
    )
  ),
  all(abs(enriched$q_value - 0.037037037037037) < 1e-14),
  identical(
    sort(unique(as.character(composition$sample_data$unit_value))),
    sort(unique(as.character(data$mouse)))
  )
)
group_sums <- stats::aggregate(
  mean_proportion ~ group_value,
  data = composition$plot_data,
  FUN = sum
)
stopifnot(all(abs(group_sums$mean_proportion - 1) < 1e-12))

comparison_columns <- setdiff(names(composition$plot_data), "panel_tag")
stopifnot(isTRUE(all.equal(
  composition$plot_data[, comparison_columns, drop = FALSE],
  untagged$composition_result$plot_data[, comparison_columns, drop = FALSE],
  tolerance = 1e-12,
  check.attributes = FALSE
)))

generator_text <- paste(readLines(
  file.path(
    repo_root,
    "Code", "in-vivo", "SI_figures", "generate_supplementary_figures.R"
  ),
  warn = FALSE
), collapse = "\n")
stopifnot(
  grepl("sys.source(shared_context_helper_path", generator_text, fixed = TRUE),
  grepl("shared_context_build_si4_panels(", generator_text, fixed = TRUE),
  grepl('"shared_context_panels_helper"', generator_text, fixed = TRUE),
  grepl('"shared_context_panels_helper_sha256"', generator_text, fixed = TRUE),
  !grepl("make_umap_discrete <- function", generator_text, fixed = TRUE),
  !grepl("figure_theme <- function", generator_text, fixed = TRUE)
)

message("Shared SI4 context-panel tests passed.")
