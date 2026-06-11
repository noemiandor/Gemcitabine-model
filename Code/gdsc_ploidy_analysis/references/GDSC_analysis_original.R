
options(java.parameters = "-Xmx7g")
library("EnrichIntersect")
library(xlsx)
dr=read.table("/Users/4470246/Projects/Resources/data/Databases/GDSC/GDSC2_fitted_dose_response_24Jul22.txt",sep="\t", header = T)
dr$CELL_LINE_NAME=toupper(gsub("-","",dr$CELL_LINE_NAME))

appCL=read.table("/Users/4470246/Projects/Resources/data/ploidyAcrossCellLines/ploidyAcrossCellLines_V1.txt",sep="\t",check.names = F, header = T)
appCL=appCL[!is.na(appCL$ploidy),];
appCL=appCL[!duplicated(appCL$`Cell iname`),]
rownames(appCL)=appCL$`Cell iname`

## correlation between ploidy and IC50 for each drug:
R=list()
WHAT="Z_SCORE"; #"LN_IC50"; #AUC
pdf("~/Downloads/drugsVsPloidyCorr.pdf",width = 15,height = 7)
par(mfrow=c(3,7))
for(CAN in c("allcancers",unique(dr$TCGA_DESC))){
  dr_=dr
  if(CAN!="allcancers"){
    dr_=dr[dr$TCGA_DESC==CAN,]
  }
  ii=intersect(dr_$CELL_LINE_NAME,rownames(appCL))
  if(length(ii)<10){
    next
  }
  r=list()
  for(drug in unique(dr_$DRUG_NAME)){
    dr__=dr_[dr_$DRUG_NAME==drug,]
    dr__=dr__[!duplicated(dr__$CELL_LINE_NAME),]
    rownames(dr__)=dr__$CELL_LINE_NAME
    if(sum(!is.na(dr__[ii,WHAT]))<10){
      next
    }
    r[[drug]]=cor(dr__[ii,WHAT],appCL[ii,"ploidy"],use="pairwise.complete.obs")
    if(abs(r[[drug]])>0.6){
      plot(dr__[ii,WHAT],appCL[ii,"ploidy"],main=paste(CAN,drug))
    }
  }
  R[[CAN]]=sort(unlist(r))
}
dev.off()
save(file="~/Downloads/drugsVsPloidyCorr.RObj","R")

## annotate drug category
coxIn=as.data.frame(matrix(data = unique(unlist(sapply(R,names)))))
colnames(coxIn)="drug"
library(matlab)
source("/Users/4470246/Projects/code/RCode/scripts/annotateFromDrugBank.R")
source("/Users/4470246/Projects/code/RCode/scripts/annotateFromPubchem.R")
# coxIn=annotateFromDrugBank(coxIn)
coxIn=annotateFromPubchem(coxIn)
save(file="~/Downloads/coxIn.RObj","coxIn")

## select drugs
load(file="~/Downloads/drugsVsPloidyCorr.RObj")
load(file="~/Downloads/coxIn.RObj")
colnames(coxIn)=c("drug","group")
rownames(coxIn)=toupper(coxIn$drug)
rownames(custom.set)=custom.set$drug
ii=intersect(rownames(custom.set),rownames(coxIn[is.na(coxIn$group),]))
coxIn[ii,"group"]=custom.set[ii,"group"]
coxIn_other=coxIn[is.na(coxIn$group),]
coxIn_other$group="NOTCLASSIFIED"
coxIn=coxIn[!is.na(coxIn$group),]
tmp=sapply(coxIn$group, function (x) strsplit(x,"; ")[[1]])
coxIn$group=sapply(tmp, function(x) x[length(x)])
coxIn$group=gsub("Cytotoxic medicines","Cytotoxic",gsub(";","",gsub(",","",coxIn$group)))
coxIn$group[grep("Alkylating",coxIn$group)]="Alkylating"
coxIn$group[grep("Topoisomerase",coxIn$group)]="Cytotoxic"
coxIn$group[grep("Tubulin",coxIn$group)]="Cytotoxic"
coxIn$group[grep("Antimitotic",coxIn$group)]="Cytotoxic"
coxIn$group[grep("Antineoplastic Agents",coxIn$group)]="Antineoplastic Agents"
coxIn$group=toupper(coxIn$group)
coxIn$group=gsub("(ANTI-)INFLAMMATORY","IMMUNOSUPPRESSIVE AGENTS",coxIn$group, fixed = T)
coxIn$group[coxIn$group %in% c("PARP INHIBITORS","SIGNAL TRANSDUCTION INHIBITORS","JAK INHIBITORS","ENZYME INHIBITORS","TARGETED THERAPIES")]="SIGNALING"
coxIn=coxIn[nchar(coxIn$group)>0,]
fr=plyr::count(coxIn$group)
coxIn=coxIn[coxIn$group %in% fr$x[fr$freq>1],]
print(fr)
# R_=sapply(R, function(x) x[abs(x)>0.2])
# R_=sapply(R_, function(x) x[c(1:min(10,sum(x< -0.2)), length(x):max(length(x)-10,sum(x>0.2)))])
coxIn$drug=toupper(coxIn$drug)
rownames(coxIn)=coxIn$drug
for(CAN in names(R)){
  names(R[[CAN]])=toupper(names(R[[CAN]]))
}
R_=sapply(R, function(x) x[names(x) %in% coxIn$drug])
# R_=sapply(R_, function(x) x[abs(x)>0.1])


## Enrichment analysis
x <- sapply(names(R_), function(x) matrix(R_[[x]],dimnames = list(names(R_[[x]]),x)))
# x <- sapply(names(R), function(x) matrix(sign(R[[x]])*as.numeric(abs(R[[x]])>0.1),dimnames = list(names(R[[x]]),x)))
lowpIsSens<-highpIsSens<- list()
for(can in names(x)){
  # lowpIsSens[[can]] <- try(enrichment(x[[can]], rbind(coxIn,coxIn_other), permute.n = 200,normalize=F,pvalue.cutoff = 0.075)$pvalue)
  # highpIsSens[[can]] <- try(enrichment(-x[[can]], rbind(coxIn,coxIn_other),permute.n = 200,normalize=F,pvalue.cutoff = 0.001)$pvalue)
  lowpIsSens[[can]] <- try(enrichment(x[[can]], coxIn, permute.n = 200,normalize=F,pvalue.cutoff = 0.05)$pvalue)
  highpIsSens[[can]] <- try(enrichment(-x[[can]], coxIn, permute.n = 200,normalize=F,pvalue.cutoff = 0.1)$pvalue)
}
groups=unique(coxIn$group)
lowpIsSens=sapply(lowpIsSens, function(x) as.data.frame(x)[groups,])
highpIsSens=sapply(highpIsSens, function(x) as.data.frame(x)[groups,])
rownames(lowpIsSens)<-rownames(highpIsSens)<- groups
lowpIsSens=lowpIsSens[,order(lowpIsSens["SIGNALING",])]
lowpIsSens=lowpIsSens[,order(lowpIsSens["CYTOTOXIC",])]
highpIsSens=highpIsSens[,order(highpIsSens["CYTOTOXIC",])]
highpIsSens=highpIsSens[,order(highpIsSens["SIGNALING",])]
write.xlsx(t(lowpIsSens), file="~/Downloads/drugsVsPloidyCorr.xlsx",sheetName = "lowpIsSens") 
write.xlsx(t(highpIsSens), file="~/Downloads/drugsVsPloidyCorr.xlsx",sheetName = "highpIsSens",append = T) 

## display
tmp=sort(unique(coxIn$group))
col=rainbow(length(tmp)*1.3)[1:length(tmp)]
names(col)=tmp
pdf("~/Downloads/ploidyVsDrugSensitivity.pdf", width = 3, height = 6)
for(x in colnames(lowpIsSens)){
  try(barplot(R_[[x]], col=col[coxIn[names(R_[[x]]),"group"]], main=x,horiz = T,las=2,cex.lab=0.7,cex.names  = 0.35, xlab="Pearson r between ploidy and drug sensitivity (IC50)"))
  # legend("topleft", names(col),fill=col,cex=0.3)
}
dev.off()


