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

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Normalized-composition test requires ggplot2", call. = FALSE)
}
sys.source(
  file.path(
    repo_root,
    "Code", "in-vivo", "SI_figures", "normalized_composition.R"
  ),
  envir = environment()
)

expect_error <- function(expression, pattern) {
  observed <- tryCatch(
    {
      force(expression)
      ""
    },
    error = function(error) conditionMessage(error)
  )
  if (!nzchar(observed) || !grepl(pattern, observed, fixed = TRUE)) {
    stop(
      "Expected error containing `", pattern, "`; observed: ",
      if (nzchar(observed)) observed else "<no error>",
      call. = FALSE
    )
  }
  invisible(observed)
}

cells_for_sample <- function(unit, group, cluster_1, cluster_2) {
  data.frame(
    unit = rep(unit, cluster_1 + cluster_2),
    group = rep(group, cluster_1 + cluster_2),
    cluster = c(rep("C1", cluster_1), rep("C2", cluster_2)),
    stringsAsFactors = FALSE
  )
}

test_theme <- function(base_size) ggplot2::theme_classic(base_size = base_size)
make_test_plot <- function(
  data,
  group_levels = c("A", "B"),
  strata_cols = character()
) {
  make_normalized_composition_plot(
    data = data,
    unit_col = "unit",
    cluster_col = "cluster",
    group_col = "group",
    cluster_levels = c("C1", "C2"),
    group_levels = group_levels,
    fill_colors = c(A = "#3366AA", B = "#CC6633"),
    bar_axis = "cluster",
    strata_cols = strata_cols,
    title = "Synthetic composition",
    subtitle = "Independent-sample test",
    x_title = "Cluster",
    y_title = "Mean within-sample proportion",
    legend_title = "Group",
    tag = "T",
    theme_function = test_theme,
    x_text_angle = 0
  )
}

# A deeply sequenced sample and a shallow sample contribute exactly one unit of
# weight apiece.  The equal-sample estimate is 0.5, not the pooled-cell value
# 91 / 110, and multiplying one sample's cell count cannot change the result.
depth_data <- rbind(
  cells_for_sample("A1", "A", 90L, 10L),
  cells_for_sample("A2", "A", 1L, 9L),
  cells_for_sample("B1", "B", 10L, 90L),
  cells_for_sample("B2", "B", 9L, 1L)
)
depth_result <- make_test_plot(depth_data)
a_c1 <- depth_result$plot_data[
  as.character(depth_result$plot_data$group_value) == "A" &
    as.character(depth_result$plot_data$cluster) == "C1",
  "mean_proportion"
]
stopifnot(length(a_c1) == 1L, abs(a_c1 - 0.5) < 1e-12)
stopifnot(abs(a_c1 - 91 / 110) > 0.25)

a1 <- depth_data[depth_data$unit == "A1", , drop = FALSE]
depth_duplicated <- rbind(depth_data, a1[rep(seq_len(nrow(a1)), 7L), ])
duplicated_result <- make_test_plot(depth_duplicated)
comparison_columns <- c("group_value", "cluster", "mean_proportion")
stopifnot(isTRUE(all.equal(
  transform(
    depth_result$plot_data[, comparison_columns],
    group_value = as.character(group_value),
    cluster = as.character(cluster)
  ),
  transform(
    duplicated_result$plot_data[, comparison_columns],
    group_value = as.character(group_value),
    cluster = as.character(cluster)
  ),
  tolerance = 1e-12,
  check.attributes = FALSE
)))

group_sums <- stats::aggregate(
  mean_proportion ~ group_value,
  data = depth_result$plot_data,
  FUN = sum
)
stopifnot(all(abs(group_sums$mean_proportion - 1) < 1e-12))
stopifnot(inherits(depth_result$plot, "ggplot"))
invisible(ggplot2::ggplot_build(depth_result$plot))

# Four independent samples per group provide enough exact label permutations
# for the two strong positive group/cluster enrichments to survive panel-wide
# BH adjustment.  Stars are emitted only for positive enrichments.
separated_data <- do.call(rbind, c(
  lapply(paste0("A", 1:4), cells_for_sample, group = "A", cluster_1 = 9L, cluster_2 = 1L),
  lapply(paste0("B", 1:4), cells_for_sample, group = "B", cluster_1 = 1L, cluster_2 = 9L)
))
separated_result <- make_test_plot(separated_data)
enriched <- separated_result$tests[separated_result$tests$enriched, ]
stopifnot(
  nrow(enriched) == 2L,
  setequal(
    paste(enriched$group_value, enriched$cluster, sep = "/"),
    c("A/C1", "B/C2")
  ),
  all(enriched$q_value <= 0.05),
  all(nzchar(enriched$significance)),
  all(
    separated_result$tests$significance[
      separated_result$tests$difference <= 0
    ] == ""
  )
)
separated_build <- ggplot2::ggplot_build(separated_result$plot)
bar_geometry <- separated_build$data[[1L]]
star_geometry <- separated_build$data[[2L]]
star_geometry <- star_geometry[nzchar(star_geometry$label), , drop = FALSE]
stopifnot(nrow(star_geometry) == 2L)
for (row_index in seq_len(nrow(star_geometry))) {
  matching_bar <- bar_geometry[
    bar_geometry$group == star_geometry$group[[row_index]] &
      abs(bar_geometry$x - star_geometry$x[[row_index]]) < 1e-12,
    ,
    drop = FALSE
  ]
  stopifnot(
    nrow(matching_bar) == 1L,
    abs(matching_bar$y - star_geometry$y[[row_index]]) < 1e-12
  )
}

stacked_result <- make_normalized_composition_plot(
  data = separated_data,
  unit_col = "unit",
  cluster_col = "cluster",
  group_col = "group",
  cluster_levels = c("C1", "C2"),
  group_levels = c("A", "B"),
  fill_colors = c(C1 = "#3366AA", C2 = "#CC6633"),
  bar_axis = "group",
  title = "Synthetic stacked composition",
  subtitle = "Stacked-star geometry",
  x_title = "Group",
  y_title = "Mean within-sample proportion",
  legend_title = "Cluster",
  tag = "S",
  theme_function = test_theme,
  x_text_angle = 0
)
stacked_build <- ggplot2::ggplot_build(stacked_result$plot)
stacked_bars <- stacked_build$data[[1L]]
stacked_stars <- stacked_build$data[[2L]]
stacked_stars <- stacked_stars[nzchar(stacked_stars$label), , drop = FALSE]
stopifnot(nrow(stacked_stars) == 2L)
for (row_index in seq_len(nrow(stacked_stars))) {
  matching_segment <- stacked_bars[
    stacked_bars$group == stacked_stars$group[[row_index]] &
      abs(stacked_bars$x - stacked_stars$x[[row_index]]) < 1e-12,
    ,
    drop = FALSE
  ]
  stopifnot(
    nrow(matching_segment) == 1L,
    abs(
      stacked_stars$y[[row_index]] -
        (matching_segment$ymin + matching_segment$ymax) / 2
    ) < 1e-12
  )
}

separated_data$stratum <- factor(
  ifelse(sub("^[AB]", "", separated_data$unit) %in% c("1", "2"), "S1", "S2"),
  levels = c("S1", "S2")
)
stratified_result <- make_test_plot(
  separated_data,
  strata_cols = "stratum"
)
stopifnot(
  all(stratified_result$tests$exact_permutations == 36L),
  all(stratified_result$tests$strata == "stratum")
)

# A display group nested inside a stratum (the SI5G design) must compare only
# with exchangeable samples in that stratum. Samples in another ploidy stratum
# must not reverse the reported effect sign.
nested_spec <- data.frame(
  unit = paste0("N", 1:8),
  group = rep(c("P1|A", "P1|B", "P2|A", "P2|B"), each = 2L),
  stratum = rep(c("P1", "P1", "P2", "P2"), each = 2L),
  c1 = c(6L, 6L, 1L, 1L, 10L, 10L, 10L, 10L),
  c2 = c(4L, 4L, 9L, 9L, 0L, 0L, 0L, 0L),
  stringsAsFactors = FALSE
)
nested_data <- do.call(rbind, lapply(seq_len(nrow(nested_spec)), function(index) {
  output <- cells_for_sample(
    nested_spec$unit[[index]],
    nested_spec$group[[index]],
    nested_spec$c1[[index]],
    nested_spec$c2[[index]]
  )
  output$stratum <- nested_spec$stratum[[index]]
  output
}))
nested_result <- make_normalized_composition_plot(
  data = nested_data,
  unit_col = "unit",
  cluster_col = "cluster",
  group_col = "group",
  cluster_levels = c("C1", "C2"),
  group_levels = c("P1|A", "P1|B", "P2|A", "P2|B"),
  fill_colors = c(
    `P1|A` = "#3366AA", `P1|B` = "#6699CC",
    `P2|A` = "#CC6633", `P2|B` = "#DD9966"
  ),
  bar_axis = "cluster",
  strata_cols = "stratum",
  title = "Nested synthetic composition",
  subtitle = "Within-stratum comparator",
  x_title = "Cluster",
  y_title = "Mean within-sample proportion",
  legend_title = "Group",
  tag = "N",
  theme_function = test_theme,
  x_text_angle = 0
)
nested_contrast <- nested_result$tests[
  nested_result$tests$group_value == "P1|A" &
    nested_result$tests$cluster == "C1",
  ,
  drop = FALSE
]
global_rest_difference <- 0.6 - mean(c(0.1, 0.1, rep(1, 4L)))
stopifnot(
  nrow(nested_contrast) == 1L,
  nested_contrast$n_group_samples == 2L,
  nested_contrast$n_rest_samples == 2L,
  nested_contrast$exact_permutations == 6L,
  abs(nested_contrast$difference - 0.5) < 1e-12,
  global_rest_difference < 0
)

# A sample cannot be split over experimental groups: that would invalidate the
# independent-sample permutation unit.
invalid <- rbind(
  cells_for_sample("shared", "A", 2L, 1L),
  cells_for_sample("shared", "B", 1L, 2L)
)
expect_error(
  make_test_plot(invalid),
  "Each independent sample must map to one display/test group"
)

# SI5F deliberately uses one displayed group per mouse.  The shared helper must
# accept group_col == unit_col while explicitly disabling impossible
# single-mouse inference rather than treating cells as replicates.
mouse_data <- do.call(rbind, lapply(seq_len(4L), function(index) {
  cells_for_sample(
    paste0("M", index),
    paste0("M", index),
    2L + index,
    7L - index
  )
}))
dose_levels <- c("0mg/kg", "30mg/kg", "120mg/kg")
mouse_metadata <- data.frame(
  unit = paste0("M", 1:4),
  ploidy = factor(c("2N", "2N", "2N", "4N"), levels = c("2N", "4N")),
  dose = factor(
    c("0mg/kg", "30mg/kg", "120mg/kg", "0mg/kg"),
    levels = dose_levels
  ),
  stringsAsFactors = FALSE
)
mouse_data <- merge(mouse_data, mouse_metadata, by = "unit", sort = FALSE)
mouse_result <- make_normalized_composition_plot(
  data = mouse_data,
  unit_col = "unit",
  cluster_col = "cluster",
  group_col = "unit",
  cluster_levels = c("C1", "C2"),
  group_levels = paste0("M", 1:4),
  fill_colors = c(C1 = "#3366AA", C2 = "#CC6633"),
  bar_axis = "group",
  x_col = "unit",
  x_levels = paste0("M", 1:4),
  facet_cols = c("ploidy", "dose"),
  facet_formula = stats::as.formula(". ~ ploidy + dose"),
  facet_type = "grid",
  facet_scales = "free_x",
  facet_space = "free_x",
  title = "Per-mouse composition",
  subtitle = "One mouse per displayed group",
  x_title = "Mouse",
  y_title = "Within-mouse proportion",
  legend_title = "Cluster",
  tag = "F",
  theme_function = test_theme,
  x_text_angle = 0,
  test_mode = "descriptive",
  descriptive_reason = "one biological sample per displayed mouse"
)
stopifnot(
  inherits(mouse_result$plot, "ggplot"),
  all(mouse_result$tests$n_group_samples == 1L),
  all(!mouse_result$tests$testable),
  all(is.na(mouse_result$tests$p_value)),
  all(is.na(mouse_result$tests$q_value)),
  all(mouse_result$tests$significance == ""),
  identical(levels(mouse_result$plot_data$dose), dose_levels)
)
mouse_build <- ggplot2::ggplot_build(mouse_result$plot)
observed_dose_order <- unique(as.character(mouse_build$layout$layout$dose))
stopifnot(identical(observed_dose_order, dose_levels))

# The production renderer must route SI4E through the shared context-panel
# builder, while the other five normalized panels continue to call the generic
# normalized-composition implementation directly.
generator <- readLines(
  file.path(
    repo_root,
    "Code", "in-vivo", "SI_figures", "generate_supplementary_figures.R"
  ),
  warn = FALSE
)
generator_text <- paste(generator, collapse = "\n")
shared_helper <- readLines(
  file.path(
    repo_root,
    "Code", "in-vivo", "SI_figures", "shared_context_panels.R"
  ),
  warn = FALSE
)
shared_helper_text <- paste(shared_helper, collapse = "\n")
generator_call_count <- lengths(regmatches(
  generator_text,
  gregexpr("make_normalized_composition_plot\\(", generator_text)
))
shared_call_count <- lengths(regmatches(
  shared_helper_text,
  gregexpr("composition_builder\\(", shared_helper_text)
))
stopifnot(
  generator_call_count == 5L,
  shared_call_count == 1L,
  grepl(
    "composition_builder = make_normalized_composition_plot",
    shared_helper_text,
    fixed = TRUE
  ),
  grepl(
    "shared_context_build_si4_panels(",
    generator_text,
    fixed = TRUE
  ),
  !grepl("s4e_result <- make_normalized_composition_plot(", generator_text, fixed = TRUE)
)
for (result_name in c(
  "s4g_result", "s5f_result", "s5g_result", "s5h_result", "s5i_result"
)) {
  pattern <- paste0(
    result_name,
    "[[:space:]]*<-[[:space:]]*make_normalized_composition_plot\\("
  )
  stopifnot(grepl(pattern, generator_text))
}

extract_production_call <- function(result_name) {
  start <- grep(paste0("^", result_name, " <-"), generator)
  stopifnot(length(start) == 1L)
  relative_end <- which(grepl("^\\)$", generator[start:length(generator)]))[[1L]]
  paste(generator[start:(start + relative_end - 1L)], collapse = "\n")
}
expected_call_contracts <- list(
  s4g_result = 'strata_cols = "context"',
  s5f_result = c(
    'test_mode = "descriptive"',
    'descriptive_reason = "one biological sample per displayed mouse"'
  ),
  s5g_result = 'strata_cols = "initial_ploidy"',
  s5h_result = 'strata_cols = "initial_ploidy"',
  s5i_result = 'strata_cols = "dose"'
)
for (result_name in names(expected_call_contracts)) {
  call_text <- extract_production_call(result_name)
  stopifnot(all(vapply(
    expected_call_contracts[[result_name]],
    grepl,
    logical(1),
    x = call_text,
    fixed = TRUE
  )))
}

shared_composition_start <- grep(
  "^  composition_result <- composition_builder\\($",
  shared_helper
)
stopifnot(length(shared_composition_start) == 1L)
shared_relative_end <- which(grepl(
  "^  \\)$",
  shared_helper[shared_composition_start:length(shared_helper)]
))[[1L]]
shared_composition_call <- paste(
  shared_helper[
    shared_composition_start:
      (shared_composition_start + shared_relative_end - 1L)
  ],
  collapse = "\n"
)
stopifnot(
  grepl('strata_cols = "initial_ploidy"', shared_composition_call, fixed = TRUE),
  grepl('unit_col = "mouse"', shared_composition_call, fixed = TRUE),
  grepl('group_col = "context"', shared_composition_call, fixed = TRUE),
  grepl('bar_axis = "cluster"', shared_composition_call, fixed = TRUE)
)

# Production metadata must retain the intended biological sample counts at all
# six call boundaries.
cells <- utils::read.csv(
  file.path(
    repo_root,
    "Data", "in-vivo", "SIfigures", "si_figures_cell_metadata.csv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
production <- data.frame(
  sample = cells$sample_id,
  cluster = cells$cluster_id,
  context = cells$context,
  initial_ploidy = cells$initial_ploidy,
  dose = cells$dose,
  included = tolower(as.character(cells$included_in_si_figures)) %in%
    c("true", "t", "1"),
  stringsAsFactors = FALSE
)
production$ploidy_dose <- paste(
  production$initial_ploidy,
  production$dose,
  sep = " | "
)
tumor_production <- production[production$included, , drop = FALSE]
sample_counts <- function(data, group_col) {
  metadata <- unique(data[, c("sample", group_col), drop = FALSE])
  stats::setNames(
    as.integer(table(metadata[[group_col]])),
    names(table(metadata[[group_col]]))
  )
}
stopifnot(
  identical(sample_counts(production, "context")[c("Tumor", "CellLine")], c(Tumor = 16L, CellLine = 2L)),
  identical(sample_counts(production, "initial_ploidy")[c("2N", "4N")], c(`2N` = 9L, `4N` = 9L)),
  all(sample_counts(tumor_production, "sample") == 1L),
  identical(
    sample_counts(tumor_production, "ploidy_dose")[c(
      "2N | 0mg/kg", "2N | 30mg/kg", "2N | 120mg/kg",
      "4N | 0mg/kg", "4N | 30mg/kg", "4N | 120mg/kg"
    )],
    c(
      `2N | 0mg/kg` = 4L, `2N | 30mg/kg` = 2L,
      `2N | 120mg/kg` = 2L, `4N | 0mg/kg` = 4L,
      `4N | 30mg/kg` = 2L, `4N | 120mg/kg` = 2L
    )
  ),
  identical(
    sample_counts(tumor_production, "dose")[c("0mg/kg", "30mg/kg", "120mg/kg")],
    c(`0mg/kg` = 8L, `30mg/kg` = 4L, `120mg/kg` = 4L)
  ),
  identical(
    sample_counts(tumor_production, "initial_ploidy")[c("2N", "4N")],
    c(`2N` = 8L, `4N` = 8L)
  )
)

cat("Normalized composition tests passed.\n")
