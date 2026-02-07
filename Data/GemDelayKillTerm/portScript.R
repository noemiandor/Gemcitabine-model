# prep_gem_data.R
# هدف: copy the two raw datasets (tracks + object tables) into a clean analysis structure,
# and write tidy "processed" tables for modeling.

suppressPackageStartupMessages({
  library(fs)
  library(data.table)
  library(stringr)
  library(arrow)
})

# -------------------------
# CONFIG: set these paths
# -------------------------

# 1) Track repo root (where your newer scripts assume things live)
src_track_repo <- path_expand("~/Repositories/Gemcitabine-model")

# 2) HALO/Incucyte classifier output directory (where ObjectData*.csv lives)
#    This should be the folder you used as setwd(...) in the older script, e.g. ".../Results_Classifier_2/"
src_object_dir  <- "/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/K01_20250314_HALOCellClassification"

# 3) Destination root for the NEW analysis structure
dest_root <- "~/Repositories/Gemcitabine-model/Data/GemDelayKillTerm"

# -------------------------
# Helpers
# -------------------------

dir_create(path(dest_root, "raw", "tracks"), recurse = TRUE)
dir_create(path(dest_root, "raw", "object_tables"), recurse = TRUE)
dir_create(path(dest_root, "raw", "platemap"), recurse = TRUE)
dir_create(path(dest_root, "processed"), recurse = TRUE)
dir_create(path(dest_root, "manifest"), recurse = TRUE)

copy_with_rel <- function(files, from_root, to_root) {
  if (length(files) == 0) {
    return(data.table(src = character(), dest = character()))
  }
  rel  <- as.character(path_rel(files, start = from_root))
  dest <- as.character(path(to_root, rel))
  dir_create(path_dir(dest), recurse = TRUE)
  file_copy(files, dest, overwrite = TRUE)
  data.table(src = as.character(files), dest = dest)
}


# Robust well/time parser for HALO Image.Location
# Your older code assumes:
#   well is token #6 after splitting by "_"
#   time is token #8 after splitting by "_", formatted like "DD-HH" (day-hour)
to_hours_from_DD_HH <- function(x) {
  # expects "DD-HH" (e.g., "02-05") -> 2*24 + 5
  # returns NA if parse fails
  x <- str_replace(x, "\\.tif$", "")
  if (!str_detect(x, "^[0-9]{2}[-_][0-9]{2}$")) return(NA_real_)
  day  <- as.numeric(str_sub(x, 1, 2))
  hour <- as.numeric(str_sub(x, 4, 5))
  day * 24 + hour
}

parse_well_time_from_image_location <- function(img_loc) {
  base <- basename(img_loc)
  
  # Parse well + optional field: A2_1 or A2
  m <- str_match(img_loc, "([A-H])([0-9]{1,2})(?:_([0-9]+))?")
  plate_row <- m[,2]
  plate_col <- suppressWarnings(as.integer(m[,3]))
  field     <- suppressWarnings(as.integer(m[,4]))  # NA if absent
  
  well <- ifelse(!is.na(plate_row) & !is.na(plate_col),
                 paste0(plate_row, plate_col),
                 NA_character_)
  
  # Parse time like 02d02h00m from basename
  tm <- str_match(base, "([0-9]+)d([0-9]+)h([0-9]+)m")
  days <- suppressWarnings(as.numeric(tm[,2]))
  hrs  <- suppressWarnings(as.numeric(tm[,3]))
  mins <- suppressWarnings(as.numeric(tm[,4]))
  
  time_hours <- days*24 + hrs + mins/60
  time_hours[is.na(days) | is.na(hrs) | is.na(mins)] <- NA_real_
  
  data.table(
    well_raw   = m[,1],     # the matched full string, e.g. "A2_1" or "A2"
    well       = well,      # canonical "A2"
    plate_row  = plate_row, # "A"
    plate_col  = plate_col, # 2
    field      = field,     # 1 (or NA)
    time_hours = time_hours,
    days       = days,
    hours      = hrs,
    minutes    = mins,
    basename   = base
  )
}




infer_track_type <- function(p) {
  # uses path segments to label track type
  if (str_detect(p, "pre-division"))  return("pre-division")
  if (str_detect(p, "post-division")) return("post-division")
  if (str_detect(p, "inter-division"))return("inter-division")
  return(NA_character_)
}

infer_well_from_filename <- function(fname) {
  # tries to pull e.g. "E2" or "E_2" from filename
  m <- str_match(fname, "([A-H]_?[0-9]{1,2})")
  m[,2]
}

# -------------------------
# 1) Collect and copy TRACK files
# -------------------------

track_data_root <- path(src_track_repo, "Data/in-vitro")

track_dirs <- dir_ls(
  track_data_root,
  recurse = TRUE,
  type = "directory",
  regexp = "CellCycleClassification"
)


# 2) Now list txt files only inside those dirs (much smaller search space)
track_files <- unlist(lapply(track_dirs, function(d) {
  dir_ls(d, type = "file", glob = "*.txt")
}), use.names = FALSE)

length(track_files)

manifest_tracks <- copy_with_rel(
  files     = track_files,
  from_root = track_data_root,
  to_root   = path(dest_root, "raw", "tracks")
)

# also copy platemap(s) if present
platemaps <- dir_ls(track_data_root, recurse = TRUE, type = "file", glob = "*.xlsx")
manifest_platemap <- copy_with_rel(
  files     = platemaps,
  from_root = track_data_root,
  to_root   = path(dest_root, "raw", "platemap")
)

# -------------------------
# 2) Collect and copy OBJECT TABLE CSV files (Alive/Dead/Transitional)
# -------------------------

object_files <- dir_ls(src_object_dir, type = "file", glob = "*.csv")

# Your older script used pattern "ata.csv" (to match "*Data.csv")
# We'll emulate that filter to avoid unrelated CSVs.
object_files <- object_files[str_detect(path_file(object_files), "ata\\.csv$|Data\\.csv$")]

manifest_objects <- copy_with_rel(
  files     = object_files,
  from_root = src_object_dir,
  to_root   = path(dest_root, "raw", "object_tables")
)

# -------------------------
# 3) Write a manifest
# -------------------------

manifest <- rbindlist(list(
  cbind(manifest_tracks, kind = "tracks_txt"),
  cbind(manifest_objects, kind = "object_csv"),
  cbind(manifest_platemap, kind = "platemap_xlsx")
), fill = TRUE)

fwrite(manifest, path(dest_root, "manifest", "file_manifest.csv"))

# -------------------------
# 4) OPTIONAL: build tidy PROCESSED tables
# -------------------------

# 4a) Object tables -> one long table
# This is what you’ll likely fit first (10-day coverage)
if (nrow(manifest_objects) > 0) {
  obj_paths <- manifest_objects$dest
  
  obj_dt <- rbindlist(lapply(obj_paths, function(fp) {
    dt <- fread(fp, select = c("Image Location","Classifier Phenotype",
                           "XMin","XMax","YMin","YMax",
                           "Region Area  μm  ","Region Perimeter  μm  "))

    dt[, source_file := path_file(fp)]
    # parse well/time from Image.Location if present
    if ("Image Location" %in% names(dt)) {
      wt <- parse_well_time_from_image_location(dt[["Image Location"]])
      dt[, well_raw   := wt$well_raw]
      dt[, well       := wt$well]
      dt[, plate_row  := wt$plate_row]
      dt[, plate_col  := wt$plate_col]
      dt[, field      := wt$field]
      dt[, time_hours := wt$time_hours]
    } else {
      dt[, `:=`(well_raw=NA_character_, well=NA_character_,
                plate_row=NA_character_, plate_col=NA_integer_,
                field=NA_integer_, time_hours=NA_real_)]
    }
    
    dt
  }), fill = TRUE)
  
  # Keep the columns most relevant for the new delay-aware kill term
  keep_cols <- intersect(names(obj_dt), c(
    "source_file","Image Location",
    "well_raw","well","plate_row","plate_col","field","time_hours",
    "Classifier Phenotype",
    "XMin","XMax","YMin","YMax",
    "Region Area  μm  ","Region Perimeter  μm  "
  ))
  
  
  obj_slim <- obj_dt[, ..keep_cols]
  
  write_parquet(obj_slim, path(dest_root, "processed", "object_cells_long.parquet"))
  
  # 4b) counts by well/time/phenotype (what your PD fitting will probably consume)
  colnames(obj_slim) = gsub(" ",".",colnames(obj_slim))
  if ("Classifier.Phenotype" %in% names(obj_slim)) {
    counts <- obj_slim[
      !is.na(well) & !is.na(time_hours) & !is.na(plate_col),
      .N,
      by = .(well, plate_row, plate_col, field, time_hours, Classifier.Phenotype)
    ]
    
    setorder(counts, well, field, time_hours, Classifier.Phenotype)
    write_parquet(counts, path(dest_root, "processed", "counts_by_well_time.parquet"))
    
    counts_well <- counts[, .(N = sum(N)), by=.(well, plate_row, plate_col, time_hours, Classifier.Phenotype)]
    write_parquet(counts_well, path(dest_root, "processed", "counts_by_well_time_wellAggregated.parquet"))
    
  }
}

# 4c) Tracks -> one long table (short window; use mainly for validation / single-cell delay)
if (nrow(manifest_tracks) > 0) {
  tr_paths <- manifest_tracks$dest
  
  tracks_long <- rbindlist(lapply(tr_paths, function(fp) {
    dt <- fread(fp, showProgress = FALSE)
    dt[, source_file := path_file(fp)]
    dt[, source_path := fp]
    dt[, track_type := infer_track_type(fp)]
    dt[, well := infer_well_from_filename(path_file(fp))]
    # If there is a time index 't' but no real time, you can convert later
    dt
  }), fill = TRUE)
  
  write_parquet(tracks_long, path(dest_root, "processed", "tracks_long.parquet"))
}

message("Done. New analysis folder created at: ", normalizePath(dest_root))
