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
setwd("~/Repositories/Gemcitabine-model/Code/in-vivo")

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
  legend.labs = c("2N", "4N") # Rename legend labels
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
