# plot_plate_counts.R
# Run from within dest_root (e.g., setwd("GemDelayKillTerm"))
setwd( "~/Repositories/Gemcitabine-model/Data/GemDelayKillTerm")
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(stringr)
  library(readxl)
  library(ggplot2)
})

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Load counts-by-well-time
# ----------------------------
counts_path <- file.path("processed", "counts_by_well_time_wellAggregated.parquet")
stopifnot(file.exists(counts_path))

dt <- as.data.table(read_parquet(counts_path))
setnames(dt, old = c("Classifier.Phenotype", "N"), new = c("phenotype", "count"), skip_absent = TRUE)

# IMPORTANT:
# Prefer plate_row/plate_col already present (from your prep script).
# Only parse from 'well' if they are missing.
if (!all(c("plate_row", "plate_col") %in% names(dt))) {
  # Parse canonical well A2 from either "A2" or "A2_1"
  m <- str_match(dt$well, "^([A-H])([0-9]{1,2})(?:_([0-9]+))?$")
  dt[, plate_row := m[,2]]
  dt[, plate_col := suppressWarnings(as.integer(m[,3]))]
  dt[, field     := suppressWarnings(as.integer(m[,4]))]  # optional
}

# Time axis
if (!"time_hours" %in% names(dt)) stop("Expected column 'time_hours' in counts table.")
dt[, time_days := time_hours / 24]
dt = dt[dt$time_days<=5,]

# Normalize phenotype labels / ordering
dt[, phenotype := factor(
  phenotype,
  levels = c("Alive", "Dead", "Transitional", "Oversegmentation"),
  ordered = TRUE
)]

# Keep only plate wells with valid row/col
dt <- dt[!is.na(plate_row) & !is.na(plate_col) & plate_col >= 1 & plate_col <= 12]

# ----------------------------
# 2) Load plate map (first xlsx)
# ----------------------------
pm_files <- list.files(file.path("raw", "platemap"), pattern = "\\.xlsx$", full.names = TRUE)
if (length(pm_files) == 0) {
  warning("No plate map xlsx found under raw/platemap/. Plot will be unlabeled.")
  plate_long <- data.table(plate_row = character(), plate_col = integer(),
                           ploidy = character(), gem = character(), condition_label = character())
  plate_source <- NA_character_
} else {
  plate_source <- basename(pm_files[1])
  pm <- read_excel(pm_files[1], sheet = 1, .name_repair = "unique")
  pm <- as.data.table(pm)
  
  # ---- Flexible plate map parsing ----
  cols_lower <- tolower(names(pm))
  
  # Case A: already long (has well / row / col)
  if (any(cols_lower %in% c("well", "wellid")) || (any(cols_lower == "row") && any(cols_lower %in% c("col","column")))) {
    
    if (!any(cols_lower %in% c("well","wellid"))) {
      row_col_name <- names(pm)[match("row", cols_lower)]
      col_col_name <- names(pm)[match(ifelse(any(cols_lower == "col"), "col", "column"), cols_lower)]
      pm[, well := paste0(get(row_col_name), get(col_col_name))]
    } else {
      well_col_name <- names(pm)[match(ifelse(any(cols_lower == "well"), "well", "wellid"), cols_lower)]
      setnames(pm, well_col_name, "well")
    }
    
    # Parse well into row/col (ignore any field suffix if present)
    m <- str_match(toupper(as.character(pm$well)), "^([A-H])([0-9]{1,2})(?:_([0-9]+))?$")
    pm[, plate_row := m[,2]]
    pm[, plate_col := suppressWarnings(as.integer(m[,3]))]
    
    ploidy_col <- names(pm)[which(tolower(names(pm)) %in% c("ploidy","ploidy_state"))][1]
    gem_col    <- names(pm)[which(str_detect(tolower(names(pm)), "gem"))][1]
    cond_col   <- names(pm)[which(tolower(names(pm)) %in% c("condition","treatment","label"))][1]
    
    pm[, ploidy := if (!is.na(ploidy_col)) as.character(get(ploidy_col)) else NA_character_]
    pm[, gem_raw := if (!is.na(gem_col)) as.character(get(gem_col)) else if (!is.na(cond_col)) as.character(get(cond_col)) else NA_character_]
    
    pm[is.na(ploidy) | ploidy == "", ploidy := str_extract(gem_raw, "(?i)\\b[24]N\\b")]
    pm[, ploidy := toupper(ploidy)]
    
    parse_gem_label <- function(x) {
      if (is.na(x)) return(NA_character_)
      m <- str_match(x, "(?i)\\b([0-9]+\\.?[0-9]*)\\s*(nM|uM|µM|pM)\\b")
      if (all(is.na(m))) return(NA_character_)
      val <- suppressWarnings(as.numeric(m[,2]))
      unit <- tolower(m[,3])
      unit <- gsub("µm", "uM", unit, fixed = TRUE)
      unit <- gsub("um", "uM", unit, fixed = TRUE)
      unit <- gsub("nm", "nM", unit, fixed = TRUE)
      if (is.na(val) | is.na(unit)) return(NA_character_)
      paste0(val, " ", unit)
    }
    pm[, gem := vapply(gem_raw, parse_gem_label, character(1))]
    
    plate_long <- unique(pm[, .(plate_row, plate_col, ploidy, gem)])
    
  } else {
    # Case B: matrix-like (rows A-H, columns 1-12)
    row_col <- names(pm)[tolower(names(pm)) == "row"][1]
    if (!is.na(row_col)) {
      setnames(pm, row_col, "Row")
    } else {
      setnames(pm, names(pm)[1], "Row")
    }
    
    pm[, Row := toupper(as.character(Row))]
    pm <- pm[Row %in% LETTERS[1:8]]
    
    colnames_wells <- setdiff(names(pm), "Row")
    plate_long <- melt(pm,
                       id.vars = "Row",
                       measure.vars = colnames_wells,
                       variable.name = "Col",
                       value.name = "condition_raw")
    plate_long[, plate_row := Row]
    plate_long[, plate_col := suppressWarnings(as.integer(str_extract(as.character(Col), "[0-9]+")))]
    
    plate_long[, ploidy := toupper(str_extract(as.character(condition_raw), "(?i)\\b[24]N\\b"))]
    plate_long[, gem := {
      m <- str_match(as.character(condition_raw), "(?i)\\b([0-9]+\\.?[0-9]*)\\s*(nM|uM|µM|pM)\\b")
      val <- suppressWarnings(as.numeric(m[,2]))
      unit <- m[,3]
      unit <- gsub("µM", "uM", unit, fixed = TRUE)
      ifelse(is.na(val) | is.na(unit), NA_character_, paste0(val, " ", unit))
    }]
    
    plate_long <- plate_long[, .(plate_row, plate_col, ploidy, gem)]
  }
  
  plate_long <- plate_long[!is.na(plate_row) & !is.na(plate_col)]
  plate_long[, condition_label := paste0(
    ifelse(!is.na(ploidy) & ploidy != "", ploidy, ""),
    ifelse(!is.na(gem) & gem != "", paste0(" | Gem ", gem), "")
  )]
  plate_long[condition_label == " | Gem " | condition_label == "", condition_label := NA_character_]
}

# ----------------------------
# 3) Join conditions onto counts
# ----------------------------
dt <- merge(dt, plate_long[, .(plate_row, plate_col, condition_label, ploidy, gem)],
            by = c("plate_row","plate_col"), all.x = TRUE)

# Panel annotation: one row per well (use max count for y-position)
labels <- dt[, .(
  time_days0 = min(time_days, na.rm = TRUE),
  y0 = max(count, na.rm = TRUE),
  condition_label = unique(na.omit(condition_label))[1]
), by = .(plate_row, plate_col)]

# ----------------------------
# 3b) Build title/subtitle with condition details
# ----------------------------
ploidy_levels <- sort(unique(na.omit(plate_long$ploidy)))
gem_levels <- sort(unique(na.omit(plate_long$gem)))
n_conditions <- plate_long[!is.na(condition_label), uniqueN(condition_label)]

subtitle <- paste(
  c(
    if (!is.na(plate_source)) paste0("Plate map: ", plate_source) else NULL,
    if (length(ploidy_levels) > 0) paste0("Ploidy: ", paste(ploidy_levels, collapse = ", ")) else NULL,
    if (length(gem_levels) > 0) paste0("Gem doses: ", paste(gem_levels, collapse = ", ")) else NULL,
    if (!is.na(n_conditions) && n_conditions > 0) paste0("Unique conditions: ", n_conditions) else NULL
  ),
  collapse = " • "
)
subtitle <- str_wrap(subtitle, width = 140)

# ----------------------------
# 4) Plot as plate layout
# ----------------------------
dt[, plate_row := factor(plate_row, levels = LETTERS[1:8])]
dt[, plate_col := factor(plate_col, levels = 1:12)]

dt_plot <- dt[phenotype %in% c("Alive","Dead","Transitional")]

p <- ggplot(dt_plot, aes(x = time_days, y = count, color = phenotype, group = phenotype)) +
  geom_line(linewidth = 0.35, na.rm = TRUE) +
  facet_grid(plate_row ~ plate_col) +
  geom_text(
    data = labels,
    aes(x = time_days0, y = y0, label = condition_label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 2.0, lineheight = 0.9
  ) +
  scale_x_continuous(name = "Time (days)") +
  scale_y_continuous(name = "Cell count") + #scale_y_log10() + 
  labs(
    color = NULL,
    title = "IncuCyte cell counts by well (Alive / Dead / Transitional)",
    subtitle = subtitle
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 8)
  )

out_pdf <- file.path("plots", "plate_counts_live_dead_transitional.pdf")
ggsave(out_pdf, p, width = 16, height = 10, units = "in", dpi = 300)

message("Saved: ", normalizePath(out_pdf))
