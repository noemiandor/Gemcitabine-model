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
dt=read.xlsx("../../Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx", sheetIndex = 1)
dt <- rbind(setNames(data.frame(matrix(NA, nrow = 2, ncol = ncol(dt))), names(dt)),dt)
dt$harvest[1:2] = c("SUM-159_NLS_2N_A7M_K_harvest","SUM-159_NLS_4N_A5M_K_harvest")
dt$Sequencing.IDs[1:2] = c("2N-Cell-Culture","4N-Cell-Culture")
dt=dt[!is.na(dt$Sequencing.IDs),];
rownames(dt)=dt$Sequencing.IDs
unambiguous_matches <- list(
  # A5 group maps to 0
  "SUM159-4N-A5-0"  = "SUM159-4N-0-0",
  "SUM159-4N-A5-R"  = "SUM159-4N-0-R",
  "SUM159-4N-A5-L"  = "SUM159-4N-0-L",
  "SUM159-4N-A5-RL" = "SUM159-4N-0-RL",
  "SUM159-4N-A5-RR" = "SUM159-4N-0-RR",
  
  # A6 group maps to 30
  "SUM159-4N-A6-0"  = "SUM159-4N-30-0",
  "SUM159-4N-A6-R"  = "SUM159-4N-30-R",
  "SUM159-4N-A6-L"  = "SUM159-4N-30-L",
  "SUM159-4N-A6-RL" = "SUM159-4N-30-RL",
  "SUM159-4N-A6-RR" = "SUM159-4N-30-RR",
  
  # A8 group maps to 120
  "SUM159-4N-A8-0"  = "SUM159-4N-120-0",
  "SUM159-4N-A8-R"  = "SUM159-4N-120-R",
  "SUM159-4N-A8-L"  = "SUM159-4N-120-L",
  "SUM159-4N-A8-RL" = "SUM159-4N-120-RL",
  "SUM159-4N-A8-RR" = "SUM159-4N-120-RR",
  
  # A7 group maps to 60
  "SUM159-4N-A7-0" =  "SUM159-4N-60-0",
  "SUM159-4N-A7-R" = "SUM159-4N-60-R",
  "SUM159-4N-A7-L" =  "SUM159-4N-60-L",
  "SUM159-4N-A7-RL" =  "SUM159-4N-60-RL",
  "SUM159-4N-A7-RR" = "SUM159-4N-60-RR"
)

mydb = cloneid::connect2DB()
q <- "SELECT id, event,media, passaged_from_id1,cellLine, correctedCount,passage, date,lastModified,owner from Passaging"
x <- dbGetQuery(mydb,q)
rownames(x) <- x$id
x=x[order(x$date,decreasing=T),]
x=x[!is.na(x$media),]
x=x[x$media==133 | x$media==134,]
x=x[grep("Gem",x$id, invert = T),]; ## exclude treatment start events

##Enter treatment start event as follow up information
seeds=x$id[x$event=="seeding"]
seeds2N=grep("2N",seeds,value=T)
seeds4N=grep("4N",seeds,value=T)
for(id in seeds2N){
  event_id = paste0(id,'_GemStart');
  diagnosis_id = id;
  treatment = 134;
  tx = "2024-09-17 00:00:00"
  x=.seed_or_harvest(event = "harvest", id=event_id, from=diagnosis_id, cellCount = 1E9, tx = tx,  media = treatment, preprocessing=F, param=NULL, inject=1)
}
tmp=read.table("~/Downloads/4N_start.txt")
tmp$V1=paste0("SUM159-4N-",gsub(":","",tmp$V1))
tmp$ids=unlist(unambiguous_matches[tmp$V1])
rownames(tmp)=tmp$ids
tmp$V2=paste(tmp$V2,tmp$V3)
for(id in seeds4N){
  event_id = paste0(id,'_GemStart');
  diagnosis_id = id;
  treatment = 134;
  tx = tmp[id,"V2"]
  x=.seed_or_harvest(event = "harvest", id=event_id, from=diagnosis_id, cellCount = 1E9, tx = tx,  media = treatment, preprocessing=F, param=NULL, inject=1)
}


## gather harvests only
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
  conf.int.style = "step",  
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
