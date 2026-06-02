# 
# resect(id="SUM159-2N-0-O", from="SUM159-2N-0-O", weight_mg=15.60, size_cubicmm=1729.65, tx = "2024-10-18 06:56:26")
# 
# 
# inject(mouseID="SUM159-2N-0-O",from="SUM-159_NLS_2N_A6M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# inject(mouseID="SUM159-2N-0-R",from="SUM-159_NLS_2N_A6M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# inject(mouseID="SUM159-2N-0-L",from="SUM-159_NLS_2N_A6M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# inject(mouseID="SUM159-2N-0-RR",from="SUM-159_NLS_2N_A6M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# inject(mouseID="SUM159-2N-0-RL",from="SUM-159_NLS_2N_A6M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# 
# 
# inject(mouseID="SUM159-4N-0-O",from="SUM-159_NLS_4N_A4M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# inject(mouseID="SUM159-4N-0-R",from="SUM-159_NLS_4N_A4M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# inject(mouseID="SUM159-4N-0-L",from="SUM-159_NLS_4N_A4M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# inject(mouseID="SUM159-4N-0-RR",from="SUM-159_NLS_4N_A4M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
# inject(mouseID="SUM159-4N-0-RL",from="SUM-159_NLS_4N_A4M_harvest", cellCount=1E6, tx = "2024-08-22 04:56:26", strain=133, injection_type=23)
resolve_in_vivo_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_match <- grep("--file=", cmd_args, value = TRUE)
  candidate_files <- character(0)
  if (length(file_match) > 0) {
    candidate_files <- c(candidate_files, sub("--file=", "", file_match[1]))
  }

  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) x$ofile else NA_character_
    },
    character(1)
  )
  candidate_files <- c(candidate_files, frame_files[!is.na(frame_files)])

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(active_path)) candidate_files <- c(candidate_files, active_path)
  }

  candidate_files <- candidate_files[!is.na(candidate_files) & nzchar(candidate_files)]
  candidate_dirs <- unique(dirname(normalizePath(candidate_files, mustWork = FALSE)))
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  cwd_parts <- strsplit(cwd, .Platform$file.sep, fixed = TRUE)[[1]]
  parent_dirs <- vapply(
    seq_along(cwd_parts),
    function(i) {
      paste(c(cwd_parts[seq_len(length(cwd_parts) - i + 1)]), collapse = .Platform$file.sep)
    },
    character(1)
  )
  parent_dirs <- parent_dirs[nzchar(parent_dirs)]
  parent_dirs <- if (grepl("^/", cwd)) paste0("/", sub("^/+", "", parent_dirs)) else parent_dirs
  candidate_dirs <- unique(c(
    candidate_dirs,
    cwd,
    file.path(cwd, "Code", "in-vivo"),
    parent_dirs,
    file.path(parent_dirs, "Code", "in-vivo")
  ))

  utils_paths <- file.path(candidate_dirs, "Utils.R")
  hit <- candidate_dirs[file.exists(utils_paths)]
  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))
  stop("Cannot locate Code/in-vivo/Utils.R from script path or working directory: ", cwd, call. = FALSE)
}

script_dir <- resolve_in_vivo_script_dir()
setwd(script_dir)

mydb = cloneid::connect2DB()
q <- "SELECT id, event,media, passaged_from_id1,cellLine, correctedCount,passage, date,lastModified,owner from Passaging"
x <- dbGetQuery(mydb,q)
rownames(x) <- x$id
x=x[order(x$date,decreasing=T),]
x=x[!is.na(x$media),]
x=x[x$media==133,]

x[x$event=="harvest","passaged_from_id1"]
alive = setdiff(x$id[x$event=="seeding"], x[x$event=="harvest","passaged_from_id1"])
saced=intersect(x$id[x$event=="seeding"], x[x$event=="harvest","passaged_from_id1"])

## compute Kaplan Meier curves
library(survminer)
library(survival)
library(xlsx)
ii=match(x$passaged_from_id1, x$id)
dt=cbind(x[,c('id','passaged_from_id1','date')],x[ii,c('id','date')])
colnames(dt)[3]="date2"
dt=dt[!is.na(dt$date),]
dt$date2=as.POSIXct(dt$date2, format = "%Y-%m-%d %H:%M:%S")
dt$date=as.POSIXct(dt$date, format = "%Y-%m-%d %H:%M:%S")
dt$deltaDays=dt$date2-dt$date
write.xlsx(dt, "~/Downloads/dt_Gem.xlsx")

# Fit a Kaplan-Meier survival model
my_data <- data.frame(
  time = dt$deltaDays,  # Time differences
  event = rep(1,nrow(dt)),      # Event occurrence (1 = occurred, 0 = censored)
  group= as.numeric(grepl('4N', dt$id))
)
# Create a survival object
surv_object <- Surv(my_data$time, my_data$event)
# Fit a Kaplan-Meier survival model, stratified by group
km_fit <- survfit(surv_object ~ group, data = my_data)
# Plot the Kaplan-Meier curves
ggsurvplot(
  km_fit,
  conf.int = TRUE,           # Show confidence intervals
  risk.table = TRUE,         # Include a risk table
  xlab = "Time", 
  ylab = "Survival Probability",
  title = "Kaplan-Meier Survival Curve by Group",
  legend.title = "Group",
  legend.labs = c("Group A", "Group B") # Rename legend labels
)


# ## stratify by ploidy at time of resection
# seqStats=read.xlsx("../../Data/in-vivo/scRNAseq_CICPT-4990-Andor-10x-BAtch01-Manifest-11062024.xlsx", sheetName = '4R')
# rownames(seqStats) = seqStats$Sample
# dt=read.xlsx("../../Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx", sheetIndex = 1)
# dt$Date.of.first.treatment=as.POSIXct(dt$Date.of.first.treatment, format = "%Y-%m-%d %H:%M:%S")
# dt$deltaDaysSinceTreatment=as.numeric(dt$date_harvest-dt$Date.of.first.treatment)
# dt$dose=as.numeric(sapply(strsplit(dt$seed,'-'),"[[",3))
# dt$group=sapply(strsplit(dt$seed,'-'),"[[",2)
# col=fliplr(rainbow(5)[1:4])
# names(col)=as.character(sort(unique(dt$dose)))
# plot(seqStats[dt$Sequencing.IDs,]$Fraction.GEMs.with..1.Cell,jitter(dt$deltaDaysSinceTreatment),col=col[as.character(dt$dose)], pch=20+as.numeric(dt$group=="4N"),cex=2, xlab='GEMs with >1 cell', ylab="Days to death")
# legend("bottomright", names(col), fill=col)
# legend("topright", c("2N group", "4N group"), pch=c(20,21))
# cor.test(seqStats[dt$Sequencing.IDs,]$Fraction.GEMs.with..1.Cell,dt$Day_27)
