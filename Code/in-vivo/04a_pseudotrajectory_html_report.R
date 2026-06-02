#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
results_root <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04a_psudo_rajectory"
}

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

summary_dir <- file.path(results_root, "02_summary")
existing_plot_dir <- file.path(results_root, "03_plots")
report_dir <- file.path(results_root, "07_highlight_report")
figure_dir <- file.path(report_dir, "figures")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

must_exist <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  path
}

read_csv_silent <- function(path) {
  readr::read_csv(must_exist(path), show_col_types = FALSE, progress = FALSE)
}

dose_cells <- read_csv_silent(file.path(summary_dir, "pseudotime_cells_all_doses.csv"))
dose_karyotype_cells <- read_csv_silent(file.path(summary_dir, "pseudotime_cells_all_dose_karyotype.csv"))
combined_cells <- read_csv_silent(file.path(summary_dir, "pseudotime_cells_all_cellline_0mg_tumor_combined.csv"))
dose_edges <- read_csv_silent(file.path(summary_dir, "principal_graph_edges_all_doses.csv"))
dose_karyotype_edges <- read_csv_silent(file.path(summary_dir, "principal_graph_edges_all_dose_karyotype.csv"))
dose_counts <- read_csv_silent(file.path(summary_dir, "dose_cell_counts.csv"))
dose_karyotype_counts <- read_csv_silent(file.path(summary_dir, "dose_karyotype_cell_counts.csv"))
origin_counts <- read_csv_silent(file.path(summary_dir, "cellline_vs_0mg_tumor_cell_counts.csv"))

ann_col <- "cluster_final_annotation_primary"
if (!(ann_col %in% colnames(dose_cells)) || !(ann_col %in% colnames(combined_cells))) {
  stop("Missing annotation column: ", ann_col, call. = FALSE)
}

dose_order <- c("0mg/kg", "30mg/kg", "120mg/kg")
karyotype_order <- c("2N", "4N")
origin_order <- c("CellLine", "0mg_Tumor")

set_ordered_factors <- function(df) {
  if ("Dose" %in% colnames(df)) {
    dose_levels <- c(dose_order, sort(setdiff(unique(as.character(df$Dose)), dose_order)))
    df$Dose <- factor(as.character(df$Dose), levels = dose_levels)
  }
  if ("Karyotype" %in% colnames(df)) {
    karyotype_levels <- c(karyotype_order, sort(setdiff(unique(as.character(df$Karyotype)), karyotype_order)))
    df$Karyotype <- factor(as.character(df$Karyotype), levels = karyotype_levels)
  }
  if ("Dose_Karyotype" %in% colnames(df)) {
    dk_levels <- as.vector(t(outer(dose_order, karyotype_order, paste, sep = "_")))
    dk_levels <- c(dk_levels, sort(setdiff(unique(as.character(df$Dose_Karyotype)), dk_levels)))
    df$Dose_Karyotype <- factor(as.character(df$Dose_Karyotype), levels = dk_levels)
  }
  if ("OriginDoseGroup" %in% colnames(df)) {
    origin_levels <- c(origin_order, sort(setdiff(unique(as.character(df$OriginDoseGroup)), origin_order)))
    df$OriginDoseGroup <- factor(as.character(df$OriginDoseGroup), levels = origin_levels)
  }
  df
}

dose_cells <- set_ordered_factors(dose_cells)
dose_karyotype_cells <- set_ordered_factors(dose_karyotype_cells)
combined_cells <- set_ordered_factors(combined_cells)
dose_edges <- set_ordered_factors(dose_edges)
dose_karyotype_edges <- set_ordered_factors(dose_karyotype_edges)

theme_report <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(color = "grey35"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 10)
    )
}

theme_umap <- function(base_size = 11) {
  theme_report(base_size) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}

save_both <- function(plot_obj, name, width, height, dpi = 300) {
  png_path <- file.path(figure_dir, paste0(name, ".png"))
  pdf_path <- file.path(figure_dir, paste0(name, ".pdf"))
  ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi)
  ggsave(pdf_path, plot_obj, width = width, height = height)
  png_path
}

pal_annotations <- function(values) {
  values <- values[!is.na(values)]
  stats::setNames(grDevices::hcl.colors(length(values), palette = "Dark 3"), values)
}

format_pct <- function(x, digits = 1) paste0(formatC(100 * x, format = "f", digits = digits), "%")
format_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

summarise_pseudotime <- function(df, group_cols) {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(pseudotime, na.rm = TRUE),
      median = median(pseudotime, na.rm = TRUE),
      q25 = as.numeric(stats::quantile(pseudotime, 0.25, na.rm = TRUE)),
      q75 = as.numeric(stats::quantile(pseudotime, 0.75, na.rm = TRUE)),
      frac_ge_0_25 = mean(pseudotime >= 0.25, na.rm = TRUE),
      frac_ge_0_5 = mean(pseudotime >= 0.5, na.rm = TRUE),
      frac_ge_0_75 = mean(pseudotime >= 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

annotation_composition <- function(df, group_cols) {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, ann_col)))) %>%
    dplyr::summarise(n = dplyr::n(), median_pseudotime = median(pseudotime, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()
}

dose_pt <- summarise_pseudotime(dose_cells, "Dose")
within_dose_karyotype_pt <- summarise_pseudotime(dose_cells, c("Dose", "Karyotype"))
combined_pt <- summarise_pseudotime(combined_cells, "OriginDoseGroup")

dose_ann <- annotation_composition(dose_cells, "Dose")
within_dose_karyotype_ann <- annotation_composition(dose_cells, c("Dose", "Karyotype"))
combined_ann <- annotation_composition(combined_cells, "OriginDoseGroup")

ann_levels <- within_dose_karyotype_ann %>%
  dplyr::group_by(.data[[ann_col]]) %>%
  dplyr::summarise(n = sum(n), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::pull(.data[[ann_col]])
ann_palette <- pal_annotations(ann_levels)

# Missing figure 1: same Dose trajectory, colored by 2N/4N.
p_within_dose_karyotype_umap <- ggplot(dose_cells, aes(x = UMAP_1, y = UMAP_2, color = Karyotype)) +
  geom_point(size = 0.55, alpha = 0.82, stroke = 0) +
  geom_segment(
    data = dose_edges,
    aes(x = UMAP_1_from, y = UMAP_2_from, xend = UMAP_1_to, yend = UMAP_2_to),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.07, "inches"), type = "closed"),
    linewidth = 0.35,
    color = "black"
  ) +
  facet_wrap(~Dose, nrow = 1) +
  coord_equal() +
  scale_color_manual(values = c("2N" = "#2A7F9E", "4N" = "#D95F02"), na.value = "grey70") +
  labs(
    title = "2N / 4N distribution within each Dose trajectory",
    subtitle = "Each panel uses the Dose-level Monocle3 trajectory, then colors cells by Karyotype",
    color = "Karyotype"
  ) +
  theme_umap()
fig_within_dose_karyotype_umap <- save_both(p_within_dose_karyotype_umap, "within_dose_karyotype_umap", 12, 5.6)

# Corrected Dose + Karyotype grid: first row 2N, second row 4N; columns 0, 30, 120.
p_dose_karyotype_grid <- ggplot(dose_karyotype_cells, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.55, alpha = 0.86, stroke = 0) +
  geom_segment(
    data = dose_karyotype_edges,
    aes(x = UMAP_1_from, y = UMAP_2_from, xend = UMAP_1_to, yend = UMAP_2_to),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.07, "inches"), type = "closed"),
    linewidth = 0.35,
    color = "black"
  ) +
  facet_grid(Karyotype ~ Dose) +
  coord_equal() +
  scale_color_gradientn(colors = c("#1B2A41", "#216869", "#F2C14E", "#C43B3B"), na.value = "grey85") +
  labs(
    title = "Monocle3 pseudotime by Dose + Karyotype",
    subtitle = "Rows are 2N then 4N; columns are 0mg/kg, 30mg/kg, 120mg/kg",
    color = "Pseudotime"
  ) +
  theme_umap()
fig_dose_karyotype_grid <- save_both(p_dose_karyotype_grid, "dose_karyotype_umap_grid", 12, 8)

p_within_dose_karyotype_dist <- ggplot(dose_cells, aes(x = Karyotype, y = pseudotime, fill = Karyotype)) +
  geom_violin(scale = "width", alpha = 0.72, linewidth = 0.25, na.rm = TRUE) +
  geom_boxplot(width = 0.16, outlier.size = 0.1, alpha = 0.85, na.rm = TRUE) +
  facet_wrap(~Dose, nrow = 1) +
  scale_fill_manual(values = c("2N" = "#2A7F9E", "4N" = "#D95F02"), na.value = "grey70") +
  labs(
    title = "Within-Dose pseudotime distribution by Karyotype",
    subtitle = "This is the direct 2N vs 4N comparison on the same Dose-level trajectory",
    x = NULL,
    y = "Pseudotime"
  ) +
  theme_report() +
  theme(legend.position = "none")
fig_within_dose_karyotype_dist <- save_both(p_within_dose_karyotype_dist, "within_dose_karyotype_pseudotime_distribution", 9.5, 5.2)

late_thresholds <- c(">=0.25" = 0.25, ">=0.50" = 0.5, ">=0.75" = 0.75)
late_fraction <- dplyr::bind_rows(lapply(names(late_thresholds), function(label) {
  threshold <- late_thresholds[[label]]
  dose_cells %>%
    dplyr::group_by(Dose, Karyotype) %>%
    dplyr::summarise(frac = mean(pseudotime >= threshold, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(threshold = label)
}))
late_fraction$threshold <- factor(late_fraction$threshold, levels = names(late_thresholds))
p_late_fraction <- ggplot(late_fraction, aes(x = Dose, y = frac, fill = Karyotype)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  facet_wrap(~threshold, nrow = 1) +
  scale_y_continuous(labels = function(x) format_pct(x, 0)) +
  scale_fill_manual(values = c("2N" = "#2A7F9E", "4N" = "#D95F02"), na.value = "grey70") +
  labs(
    title = "Fraction of late-pseudotime cells within each Dose",
    subtitle = "Fractions are shown separately for 2N and 4N cells within each Dose-level trajectory",
    x = NULL,
    y = "Cell fraction",
    fill = "Karyotype"
  ) +
  theme_report()
fig_late_fraction <- save_both(p_late_fraction, "within_dose_karyotype_late_fraction", 10, 4.8)

p_dose_ann <- ggplot(dose_ann, aes(x = Dose, y = frac, fill = .data[[ann_col]])) +
  geom_col(width = 0.72) +
  scale_y_continuous(labels = function(x) format_pct(x, 0)) +
  scale_fill_manual(values = ann_palette, name = "Annotation") +
  labs(
    title = "Functional-state composition by Dose",
    subtitle = "Cell fractions are grouped by the primary cluster annotation",
    x = NULL,
    y = "Cell fraction"
  ) +
  theme_report()
fig_dose_ann <- save_both(p_dose_ann, "dose_annotation_composition", 9.5, 5.8)

p_within_dose_karyotype_ann <- ggplot(within_dose_karyotype_ann, aes(x = Karyotype, y = frac, fill = .data[[ann_col]])) +
  geom_col(width = 0.72) +
  facet_wrap(~Dose, nrow = 1) +
  scale_y_continuous(labels = function(x) format_pct(x, 0)) +
  scale_fill_manual(values = ann_palette, name = "Annotation") +
  labs(
    title = "Functional-state composition by Karyotype within Dose",
    subtitle = "Cell fractions are shown for 2N and 4N cells within each Dose",
    x = NULL,
    y = "Cell fraction"
  ) +
  theme_report()
fig_within_dose_karyotype_ann <- save_both(p_within_dose_karyotype_ann, "within_dose_karyotype_annotation_composition", 12, 5.8)

delta_rows <- list()
for (dose in levels(dose_cells$Dose)) {
  sub <- within_dose_karyotype_ann %>% dplyr::filter(Dose == dose)
  for (ann in ann_levels) {
    frac_2n <- sub %>% dplyr::filter(Karyotype == "2N", .data[[ann_col]] == ann) %>% dplyr::pull(frac)
    frac_4n <- sub %>% dplyr::filter(Karyotype == "4N", .data[[ann_col]] == ann) %>% dplyr::pull(frac)
    delta_rows[[length(delta_rows) + 1]] <- data.frame(
      Dose = dose,
      annotation = ann,
      delta_4N_minus_2N = (ifelse(length(frac_4n) == 0, 0, frac_4n[1]) - ifelse(length(frac_2n) == 0, 0, frac_2n[1])),
      stringsAsFactors = FALSE
    )
  }
}
karyotype_delta <- dplyr::bind_rows(delta_rows)
delta_ann_order <- karyotype_delta %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarise(max_abs_delta = max(abs(delta_4N_minus_2N)), .groups = "drop") %>%
  dplyr::arrange(max_abs_delta) %>%
  dplyr::pull(annotation)
karyotype_delta$annotation <- factor(karyotype_delta$annotation, levels = delta_ann_order)
karyotype_delta$Dose <- factor(karyotype_delta$Dose, levels = dose_order)
p_karyotype_delta <- ggplot(karyotype_delta, aes(x = Dose, y = annotation, fill = delta_4N_minus_2N)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(ifelse(delta_4N_minus_2N >= 0, "+", ""), format_pct(delta_4N_minus_2N, 1))), size = 3.2) +
  scale_fill_gradient2(low = "#2A7F9E", mid = "white", high = "#D95F02", midpoint = 0, labels = function(x) format_pct(x, 0), name = "4N - 2N") +
  labs(
    title = "Annotation fraction difference by Karyotype within each Dose",
    subtitle = "Positive values indicate a higher fraction in 4N; negative values indicate a higher fraction in 2N",
    x = NULL,
    y = NULL
  ) +
  theme_report()
fig_karyotype_delta <- save_both(p_karyotype_delta, "within_dose_karyotype_annotation_delta_heatmap", 9.5, 5.8)

p_combined_ann <- ggplot(combined_ann, aes(x = OriginDoseGroup, y = frac, fill = .data[[ann_col]])) +
  geom_col(width = 0.68) +
  scale_y_continuous(labels = function(x) format_pct(x, 0)) +
  scale_fill_manual(values = ann_palette, name = "Annotation") +
  labs(
    title = "Functional-state composition: CellLine vs 0mg Tumor",
    subtitle = "Cell fractions are grouped by OriginDoseGroup and primary cluster annotation",
    x = NULL,
    y = "Cell fraction"
  ) +
  theme_report()
fig_combined_ann <- save_both(p_combined_ann, "cellline_0mg_tumor_annotation_composition", 8.5, 5.8)

combined_delta <- data.frame(annotation = ann_levels, stringsAsFactors = FALSE)
combined_delta$cellline_frac <- vapply(combined_delta$annotation, function(ann) {
  val <- combined_ann %>% dplyr::filter(OriginDoseGroup == "CellLine", .data[[ann_col]] == ann) %>% dplyr::pull(frac)
  ifelse(length(val) == 0, 0, val[1])
}, numeric(1))
combined_delta$tumor_frac <- vapply(combined_delta$annotation, function(ann) {
  val <- combined_ann %>% dplyr::filter(OriginDoseGroup == "0mg_Tumor", .data[[ann_col]] == ann) %>% dplyr::pull(frac)
  ifelse(length(val) == 0, 0, val[1])
}, numeric(1))
combined_delta$delta_tumor_minus_cellline <- combined_delta$tumor_frac - combined_delta$cellline_frac
combined_delta$annotation <- factor(combined_delta$annotation, levels = combined_delta$annotation[order(combined_delta$delta_tumor_minus_cellline)])
p_combined_delta <- ggplot(combined_delta, aes(x = annotation, y = delta_tumor_minus_cellline, fill = delta_tumor_minus_cellline > 0)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_y_continuous(labels = function(x) format_pct(x, 0)) +
  scale_fill_manual(values = c("TRUE" = "#0FA3B1", "FALSE" = "#E76F51"), guide = "none") +
  labs(
    title = "Annotation fraction difference: 0mg Tumor minus CellLine",
    subtitle = "Positive values indicate higher fractions in 0mg Tumor; negative values indicate higher fractions in CellLine",
    x = NULL,
    y = "Fraction difference"
  ) +
  theme_report()
fig_combined_delta <- save_both(p_combined_delta, "cellline_0mg_tumor_annotation_delta", 8.5, 5.4)

combined_late_fraction <- dplyr::bind_rows(lapply(names(late_thresholds), function(label) {
  threshold <- late_thresholds[[label]]
  combined_cells %>%
    dplyr::group_by(OriginDoseGroup) %>%
    dplyr::summarise(frac = mean(pseudotime >= threshold, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(threshold = label)
}))
combined_late_fraction$threshold <- factor(combined_late_fraction$threshold, levels = names(late_thresholds))
p_combined_late_fraction <- ggplot(combined_late_fraction, aes(x = threshold, y = frac, fill = OriginDoseGroup)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  scale_y_continuous(labels = function(x) format_pct(x, 0)) +
  scale_fill_manual(values = c("CellLine" = "#E76F51", "0mg_Tumor" = "#0FA3B1"), name = "Group") +
  labs(
    title = "Late-pseudotime fraction in the combined CellLine / 0mg Tumor trajectory",
    subtitle = "Fractions are calculated separately for CellLine and 0mg Tumor at each threshold",
    x = "Pseudotime threshold",
    y = "Cell fraction"
  ) +
  theme_report()
fig_combined_late_fraction <- save_both(p_combined_late_fraction, "cellline_0mg_tumor_late_fraction", 8, 4.8)

late_cells <- combined_cells %>%
  dplyr::filter(pseudotime >= 0.5)
late_ann <- annotation_composition(late_cells, "OriginDoseGroup")
p_late_ann <- ggplot(late_ann, aes(x = OriginDoseGroup, y = frac, fill = .data[[ann_col]])) +
  geom_col(width = 0.68) +
  scale_y_continuous(labels = function(x) format_pct(x, 0)) +
  scale_fill_manual(values = ann_palette, name = "Annotation") +
  labs(
    title = "Composition of late cells in the combined trajectory",
    subtitle = "Late cells are defined as pseudotime >= 0.50",
    x = NULL,
    y = "Fraction within late cells"
  ) +
  theme_report()
fig_late_ann <- save_both(p_late_ann, "cellline_0mg_tumor_late_annotation_composition", 8.5, 5.8)

copy_existing <- function(rel_path) {
  file.path("..", rel_path)
}

plot_ref <- function(path, caption, alt = caption) {
  sprintf(
    '<figure><img src="%s" alt="%s"><figcaption>%s</figcaption></figure>',
    path,
    html_escape(alt),
    html_escape(caption)
  )
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

rel_report_figure <- function(path) file.path("figures", basename(path))

count_for_group <- function(df, group_col, count_col, group_label) {
  group_values <- as.character(df[[group_col]])
  hit <- !is.na(group_values) & group_values == group_label
  values <- df[[count_col]][hit]
  if (length(values) == 0 || is.na(values[1])) return("NA")
  as.character(values[1])
}

html_table <- function(df, digits = 3) {
  out <- df
  for (col in colnames(out)) {
    if (is.numeric(out[[col]])) {
      if (grepl("^frac|fraction|delta", col)) {
        out[[col]] <- format_pct(out[[col]], 1)
      } else if (col %in% c("mean", "median", "q25", "q75")) {
        out[[col]] <- format_num(out[[col]], digits)
      } else {
        out[[col]] <- as.character(out[[col]])
      }
    } else {
      out[[col]] <- as.character(out[[col]])
    }
  }
  header <- paste0("<tr>", paste0("<th>", html_escape(colnames(out)), "</th>", collapse = ""), "</tr>")
  rows <- apply(out, 1, function(row) paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>"))
  paste0('<table class="data-table">', header, paste(rows, collapse = "\n"), "</table>")
}

dose_key <- dose_pt %>%
  dplyr::select(Dose, n, mean, median, q25, q75, frac_ge_0_5, frac_ge_0_75)
within_karyotype_key <- within_dose_karyotype_pt %>%
  dplyr::select(Dose, Karyotype, n, mean, median, q25, q75, frac_ge_0_5, frac_ge_0_75)
combined_key <- combined_pt %>%
  dplyr::select(OriginDoseGroup, n, mean, median, q25, q75, frac_ge_0_25, frac_ge_0_5, frac_ge_0_75)

cellline_count <- count_for_group(origin_counts, "OriginDoseGroup", "n_cells", "CellLine")
tumor_0mg_count <- count_for_group(origin_counts, "OriginDoseGroup", "n_cells", "0mg_Tumor")

summary_numbers <- list(
  dose_0_median = dose_key %>% dplyr::filter(Dose == "0mg/kg") %>% dplyr::pull(median),
  dose_30_median = dose_key %>% dplyr::filter(Dose == "30mg/kg") %>% dplyr::pull(median),
  dose_120_median = dose_key %>% dplyr::filter(Dose == "120mg/kg") %>% dplyr::pull(median),
  tumor_hypoxia = combined_delta %>% dplyr::filter(annotation == "Hypoxia") %>% dplyr::pull(tumor_frac),
  cellline_uv = combined_delta %>% dplyr::filter(annotation == "Uv Response Dn") %>% dplyr::pull(cellline_frac)
)

css <- '
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 0; color: #1f2933; background: #f7f8fa; }
header { background: #111827; color: white; padding: 34px 44px; }
header h1 { margin: 0 0 8px 0; font-size: 34px; letter-spacing: 0; }
header p { margin: 0; color: #d1d5db; font-size: 15px; }
main { max-width: 1280px; margin: 0 auto; padding: 28px 36px 56px; }
section { background: white; border: 1px solid #e5e7eb; border-radius: 8px; padding: 24px; margin: 0 0 22px; box-shadow: 0 1px 2px rgba(0,0,0,0.03); }
h2 { margin: 0 0 14px; font-size: 24px; }
h3 { margin: 22px 0 10px; font-size: 18px; }
p { line-height: 1.55; }
ul { line-height: 1.55; }
.cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin: 18px 0; }
.card { border: 1px solid #e5e7eb; border-radius: 8px; padding: 14px; background: #fafafa; }
.metric { font-size: 24px; font-weight: 750; color: #111827; }
.label { font-size: 12px; color: #667085; margin-top: 4px; }
.figure-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px; align-items: start; }
figure { margin: 16px 0; }
figure img { width: 100%; height: auto; border: 1px solid #e5e7eb; border-radius: 6px; background: white; }
figcaption { margin-top: 7px; color: #4b5563; font-size: 13px; line-height: 1.45; }
.data-table { border-collapse: collapse; width: 100%; margin: 12px 0 18px; font-size: 13px; }
.data-table th { text-align: left; background: #eef2f7; border: 1px solid #d8dee9; padding: 7px 8px; }
.data-table td { border: 1px solid #e5e7eb; padding: 7px 8px; }
.note { border-left: 4px solid #f59e0b; background: #fffbeb; padding: 12px 14px; margin: 14px 0; color: #713f12; }
.source { font-size: 12px; color: #6b7280; }
@media (max-width: 900px) { .cards, .figure-grid { grid-template-columns: 1fr; } main { padding: 18px; } header { padding: 26px 24px; } }
'

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>04a Pseudotrajectory Results Report</title>
<style>', css, '</style>
</head>
<body>
<header>
<h1>04a Monocle3 pseudotime results report</h1>
<p>Result directory: ', html_escape(results_root), '</p>
</header>
<main>
<section>
<h2>Result Summary</h2>
<p>This report summarizes the Monocle3 pseudotime outputs for Dose, Dose + Karyotype, direct 2N/4N contrasts within each Dose-level trajectory, and the combined CellLine vs 0mg Tumor trajectory. Additional figures generated for this report are stored in the <code>figures</code> subdirectory.</p>
<div class="cards">
<div class="card"><div class="metric">', format_num(summary_numbers$dose_120_median, 3), '</div><div class="label">120mg/kg median pseudotime</div></div>
<div class="card"><div class="metric">', format_num(summary_numbers$dose_30_median, 3), '</div><div class="label">30mg/kg median pseudotime</div></div>
<div class="card"><div class="metric">', format_pct(summary_numbers$tumor_hypoxia, 1), '</div><div class="label">0mg Tumor Hypoxia fraction</div></div>
<div class="card"><div class="metric">', format_pct(summary_numbers$cellline_uv, 1), '</div><div class="label">CellLine UV Response Dn fraction</div></div>
</div>
<ul>
<li><b>Dose:</b> The median pseudotime is 0.103 for 0mg/kg, 0.185 for 30mg/kg, and 0.332 for 120mg/kg.</li>
<li><b>Karyotype:</b> Within each Dose-level trajectory, the fraction of cells above the late-pseudotime thresholds is tabulated separately for 2N and 4N cells.</li>
<li><b>30mg/kg_2N:</b> In the independently learned Dose + Karyotype trajectory set, this group has a higher median pseudotime than 30mg/kg_4N.</li>
<li><b>CellLine vs Tumor:</b> In the combined trajectory, CellLine and 0mg Tumor occupy overlapping UMAP regions and have different primary annotation fractions.</li>
</ul>
<div class="note">Interpretation note: pseudotime is independently scaled to 0-1 for each separately learned trajectory. Across independent trajectories, the values describe the trajectory learned for that group. Direct numeric comparisons use cells from the same trajectory, such as 2N vs 4N within a Dose-level trajectory or CellLine vs 0mg Tumor within the combined trajectory.</div>
</section>

<section>
<h2>1. Dose-level trajectory</h2>
<p>This section shows the Dose-level trajectories, pseudotime distributions, and primary annotation fractions for 0mg/kg, 30mg/kg, and 120mg/kg.</p>
', html_table(dose_key), '
<div class="figure-grid">
', plot_ref("../03_plots/umap_monocle3_pseudotime_by_Dose.png", "Dose-level pseudotime UMAP."), '
', plot_ref("../03_plots/pseudotime_distribution_by_Dose.png", "Dose-level pseudotime distribution."), '
', plot_ref(rel_report_figure(fig_dose_ann), "Primary annotation fractions by Dose."), '
</div>
</section>

<section>
<h2>2. Dose + Karyotype results</h2>
<p>This section includes the independently learned Dose + Karyotype trajectories and the direct 2N/4N comparison within each Dose-level trajectory.</p>
', html_table(within_karyotype_key), '
<div class="figure-grid">
', plot_ref(rel_report_figure(fig_dose_karyotype_grid), "Dose + Karyotype pseudotime UMAP. Rows are 2N and 4N; columns are 0mg/kg, 30mg/kg, and 120mg/kg."), '
', plot_ref(rel_report_figure(fig_within_dose_karyotype_umap), "Dose-level trajectories colored by 2N/4N."), '
', plot_ref(rel_report_figure(fig_within_dose_karyotype_dist), "Within-Dose pseudotime distribution by Karyotype."), '
', plot_ref(rel_report_figure(fig_late_fraction), "Late-pseudotime fractions by Dose and Karyotype."), '
', plot_ref(rel_report_figure(fig_within_dose_karyotype_ann), "Primary annotation fractions by Karyotype within Dose."), '
', plot_ref(rel_report_figure(fig_karyotype_delta), "Primary annotation fraction difference, 4N minus 2N, within each Dose."), '
</div>
</section>

<section>
<h2>3. CellLine vs 0mg Tumor: combined trajectory</h2>
<p>This section shows the combined trajectory for CellLine and 0mg Tumor cells, together with pseudotime distributions and primary annotation fractions.</p>
', html_table(combined_key), '
<div class="figure-grid">
', plot_ref("../03_plots/cellline_0mg_tumor_combined_by_OriginDoseGroup.png", "Combined trajectory colored by OriginDoseGroup."), '
', plot_ref("../03_plots/cellline_0mg_tumor_combined_by_OriginDoseGroup_pseudotime_distribution.png", "Combined pseudotime distribution by OriginDoseGroup."), '
', plot_ref(rel_report_figure(fig_combined_ann), "Primary annotation fractions for CellLine and 0mg Tumor."), '
', plot_ref(rel_report_figure(fig_combined_delta), "Primary annotation fraction difference, 0mg Tumor minus CellLine."), '
', plot_ref(rel_report_figure(fig_combined_late_fraction), "Late-pseudotime fractions in the combined trajectory."), '
', plot_ref(rel_report_figure(fig_late_ann), "Primary annotation fractions among cells with pseudotime >= 0.50."), '
</div>
</section>

<section>
<h2>4. Data and reproducibility notes</h2>
<ul>
<li>CellLine cell count: ', html_escape(cellline_count), '; 0mg Tumor cell count: ', html_escape(tumor_0mg_count), '.</li>
<li>The Dose analysis does not include CellLine cells because their Dose value is NA. The CellLine vs Tumor comparison uses the separate OriginDoseGroup design.</li>
<li>The root cluster is currently cluster 0 for all trajectory runs. If a different biological starting state is preferred, pseudotime should be rerun with an explicit root.</li>
<li>A balanced downsampling check is recommended for the combined CellLine/Tumor analysis to ensure that the larger CellLine population is not dominating the trajectory structure.</li>
</ul>
<p class="source">Generated from CSV and image outputs under <code>', html_escape(results_root), '</code>.</p>
</section>
</main>
</body>
</html>')

if (length(html) != 1) {
  stop("Internal HTML assembly error: expected one document, got ", length(html), ".", call. = FALSE)
}

html_path <- file.path(report_dir, "index.html")
writeLines(html, con = html_path, useBytes = TRUE)
message("Wrote HTML report: ", html_path)
message("Wrote report figures: ", figure_dir)
