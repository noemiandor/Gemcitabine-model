library(ggplot2)
library(GSVA)
library(RColorBrewer)
library(ggplot2)
library(matlab)
source("~/Projects/code/RCode/scripts/plotCorr.R")
source("~/Projects/code/RCode/scripts/getAllPathways.R")
OUTDIR = "~/Projects/PMO/InferringMultiSamplePhylo/results/MEP-LINCS"
LOCCOLS=c("chr","startpos","endpos")
HMCOLS=fliplr(brewer.pal(11,"RdBu")); 
appCL=read.table("/Users/4470246/Projects/Resources/data/Databases/Cmap_LINCS/L1000/Cell_app_export.txt",sep="\t",check.names = F, stringsAsFactors = F, header = T, quote = "", comment.char = "")

ii = grep("BREAST", appCL$`CCLE name`)

## Expression
ex=read.table("~/Projects/Resources/data/Databases/CCLE/CCLE.rpkm.2016-06-17f.gct", skip = 2, header = T, check.names = F, stringsAsFactors = F, sep="\t")
ex=ex[!duplicated(ex$Description),  ]
rownames(ex)=ex$Description
colnames(ex)=gsub("-Tumor$","",gsub("^fh_","",colnames(ex)))
ex=ex[,-(1:2)]
coi = intersect(colnames(ex), appCL[ii,]$`CCLE name`)
## Quantify pathway activity
gs=getAllPathways(include_genes=T, loadPresaved = T);     
gs=gs[sapply(gs, length)>=5]
pq <- gsva(as.matrix(ex[,coi]), gs, kcdf="Poisson", mx.diff=T, verbose=FALSE, parallel.sz=2, min.sz=10)

## Ploidy
ploidy = list()
ploidy[["AU565_BREAST"]] = "3.614"           
ploidy[["BT20_BREAST"]] = "2.445"
ploidy[["BT474_BREAST"]] = "2.693"            
ploidy[["BT483_BREAST"]] = "3.84"  
ploidy[["BT549_BREAST"]] = "2.966"       
ploidy[["CAL120_BREAST"]] = "3.198"      
ploidy[["CAL148_BREAST"]] = "2.066"     
ploidy[["CAL51_BREAST"]] = "1.988"      
ploidy[["CAL851_BREAST"]] = "2.675"      
ploidy[["CAMA1_BREAST"]] = "1.933"      
ploidy[["DU4475_BREAST"]] = "3.744"      
ploidy[["EFM19_BREAST"]] = "2.838"        
ploidy[["EFM192A_BREAST"]] = "3.149"     
ploidy[["EVSAT_BREAST"]] = "2.839"       
ploidy[["HCC1143_BREAST"]] = "3.362"     
ploidy[["HCC1187_BREAST"]] = "2.644"    
ploidy[["HCC1395_BREAST"]] = "2.692"    
ploidy[["HCC1419_BREAST"]] = "3.555"     
ploidy[["HCC1428_BREAST"]] = "3.515"      
ploidy[["HCC1500_BREAST"]] = "1.661"     
ploidy[["HCC1569_BREAST"]] = "2.948"    
ploidy[["HCC1599_BREAST"]] = "2.995"     
ploidy[["HCC1806_BREAST"]] = "2.215"     
ploidy[["HCC1937_BREAST"]] = "4.238"     
ploidy[["HCC1954_BREAST"]] = "4.204"      
ploidy[["HCC202_BREAST"]] = "2.919"       
ploidy[["HCC2157_BREAST"]] = "2.69"    
ploidy[["HCC2218_BREAST"]] = "3.925"    
ploidy[["HCC38_BREAST"]] = "1.823"      
ploidy[["HCC70_BREAST"]] = "4.243"      
ploidy[["HDQP1_BREAST"]] = "4.139"    
ploidy[["HS578T_BREAST"]] = "2.603"    
ploidy[["JIMT1_BREAST"]] = "2.511"
ploidy[["MDAMB157_BREAST"]] = "2.812"
ploidy[["MCF7_BREAST"]] = "2.974"
ploidy[["MDAMB175VII_BREAST"]] = "3.627"
ploidy[["MDAMB231_BREAST"]] = "2.777"   
ploidy[["MDAMB361_BREAST"]] = "2.655"   
ploidy[["MDAMB415_BREAST"]] = "2.875"   
ploidy[["MDAMB436_BREAST"]] = "2.973"   
ploidy[["MDAMB453_BREAST"]] = "4.199"   
ploidy[["MDAMB468_BREAST"]] = "2.842"   
ploidy[["T47D_BREAST"]] = "2.777"
ploidy[["UACC812_BREAST"]] = "2.912"    
ploidy[["UACC893_BREAST"]] = "2.947"    
ploidy[["ZR7530_BREAST"]] = "3.741" 
ploidy = ploidy[names(ploidy) %in% coi]
ploidy = ploidy[appCL[match(names(ploidy), appCL$`CCLE name`),"Donor tumor phase"]=="primary" ]
ploidy = ploidy[appCL[match(names(ploidy), appCL$`CCLE name`),"Growth pattern"]!="suspension" ]
p=as.numeric(unlist(ploidy))
names(p) = names(ploidy)
ploidy = sort(p)
lowp = names(ploidy)[ploidy<=3.4]
hip = names(ploidy)[ploidy>3.4]
plotCorr(as.numeric(appCL$`Cellosaurus Doubling time`[match(names(ploidy), appCL$`CCLE name`)]),ploidy, xlab="Doubling Time", ylab="Ploidy", log = "")

## Correlation of ploidy to pathway activity
c2p = apply(pq[, names(ploidy)],1, function(x) unlist(cor.test(as.numeric(x), ploidy)[c("estimate","p.value")]))
c2p = as.data.frame(t(c2p))
# c2p$p.value.adj = p.adjust(c2p$p.value, method = "fdr" )
c2p = c2p[sort(abs(c2p$estimate.cor), decreasing = T, index.return=T)$ix,]
for(x in c("Hyaluronan metabolism", "Metabolism of vitamins and cofactors", "RHO GTPases Activate NADPH Oxidases")){
  tmp = as.data.frame(cbind(pq[x, names(ploidy)], ploidy))
  tmp$COI = "blue"
  tmp[names(ploidy)=="HCC1954_BREAST","COI"] = 'red'
  colnames(tmp)[1:2] = c("pathway","ploidy")
  p <- ggplot(tmp,aes(pathway, ploidy)) +
    geom_point(color=tmp$COI) + labs(x=x) +
    geom_smooth(method='lm', formula= y~x) 
  ggsave(paste0(OUTDIR, filesep,x,".pdf"), plot = p, width = 3, height = 3)
}
hist(as.numeric(unlist(ploidy)), 10, col="blue", border = "white", xlab="ploidy", ylab="Number of cell lines", cex.lab=1.4)
write.table(c2p[c2p$p.value<=0.05,], paste0(OUTDIR, filesep,"Ploidy_PathwayActivity_Corr_CCLE.txt"), sep="\t", quote = F)

# heatmap: pathways involved in nutrient sensing (X) and motility (Y)
poi = names(which(sapply(gs, function(x) !isempty(grep("GCK",x)))))
poi = c(poi,names(which(sapply(gs, function(x) !isempty(grep("SLC2A2",x))))))
poi = c(poi,names(which(sapply(gs, function(x) !isempty(grep("AMPK",x))))))
poi = c(poi,grep("tegrin", rownames(pq), value=T))
poi = c(poi, grep("igrat", rownames(pq), value=T))
poi = c(poi, grep("otil", rownames(pq), value=T))
poi = c(poi, grep("cose", rownames(pq), value=T))
poi = c(poi, grep("RHO", rownames(pq), value=T))
# poi = c(poi, grep("port", rownames(pq), value=T))
poi = c(poi, grep("TOR", rownames(pq), value=T))
poi = unique(poi)
gplots::heatmap.2(t(pq[poi,]),margins=c(16,10),trace ="none",symm = F, cexCol = 0.73, col=HMCOLS);
cc=matrix(NA, length(poi), length(poi))
rownames(cc) <- colnames(cc) <- poi
for(p1 in poi){
  for(p2 in poi){
    cc[p1,p2] = cor(pq[p1,names(ploidy)] / pq[p2,names(ploidy)], ploidy)    
  }
}
gplots::heatmap.2(cc,margins=c(16,10),trace ="none",symm  = F, cexCol = 0.53, col=HMCOLS);
poi = which(abs(cc)>0.6, arr.ind = T)
poi = cbind(rownames(cc)[poi[,"row"]], rownames(cc)[poi[,"col"]] )
B = "SEMA3A-Plexin repulsion signaling by inhibiting Integrin adhesion";
A="SLC-mediated transmembrane transport"
plot(pq[B, names(ploidy)]/ pq[A, names(ploidy)], ploidy, pch=20, xlab=paste(B ,"/",A))
cor.test(pq[B, names(ploidy)]/ pq[A, names(ploidy)], ploidy)

# Pairwise correlations
poi=rownames(pq)
cc = cor(t(pq[poi,lowp]))
cc_ = cor(t(pq[poi,hip]))
poi = which(abs(cc) - abs(cc_) > 0.9, arr.ind = T)
poi = cbind(rownames(cc)[poi[,"row"]], rownames(cc)[poi[,"col"]] )
cc[cc==1]=NA
gplots::heatmap.2(cc[poi,poi],margins=c(16,10),trace ="none",symm  = F, cexCol = 0.73, col=HMCOLS);

# Interactions
ucols=rainbow(max(ploidy)*10)
A = "RHO GTPases Activate Formins"; #"Integrin cell surface interactions";
B = "Facilitative Na+-independent glucose transporters"; #"Transport of vitamins, nucleosides, and related molecules";
plot(pq[A , names(ploidy)],  pq[B ,names(ploidy)], pch=20, col=ucols[round(ploidy*10)], xlab=A, ylab=B)
phytools::add.color.bar(0.4,ucols[sort(round(ploidy*10))], subtitle = "", title="Ploidy", lim=quantile(ploidy,c(0,1)))
text(pq[A , hip],  pq[B ,hip], labels = sapply(strsplit(hip,"_"),"[[",1), cex = 0.5)
# 



### Compare ECM specific inferred model paramters to pathway and receptor expression in HCC1954 cell line
ecm_param = read.table(paste0(OUTDIR, filesep,"BestFits_per_ECM.txt"), sep="\t", check.names = F, stringsAsFactors = F, header = T);
grep("ITGA", rownames(ecm_param))
## Split ECMs with same best fits
cols = colnames(ecm_param)
tmp = sapply(ecm_param$ECMs, strsplit,", ")
if(max(sapply(tmp, length))>1){
  ecm_param = sapply(1:nrow(ecm_param), function(i) matrix(rep(ecm_param[i,,drop=F],each=length(tmp[[i]])),nrow=length(tmp[[i]]))  )
  ecm_param = as.data.frame(do.call(rbind, ecm_param))
  colnames(ecm_param) = cols
  ecm_param$ECMs = unlist(tmp)
}
rownames(ecm_param) = ecm_param$ECMs
## Done splitting
## Fix names
ii = grep("ITGA", rownames(ecm_param))
rownames(ecm_param)[ii]=sapply(strsplit(rownames(ecm_param)[ii],"B"),"[[",1)
## For each ECM, quantify cummulative expression of interacting genes in HCC1954
ii = match(rownames(ecm_param), rownames(ex))
ecm_param$ECM_RNAseq_Expression = ex[ii,"HCC1954_BREAST"]
ecm_param = ecm_param[!is.na(ecm_param$ECM_RNAseq_Expression),]
coi = c("eta","xi_u","k_u","a","Xi")
## For each ECM, quantify cummulative expression of pathways including it
gs = gs[rownames(pq)]
tmp = sapply(rownames(ecm_param), function(x) names(gs)[sapply(gs, function(y) any(x==y))])
ecm_param$ECMpathway_RNAseq = sapply(tmp, function(pathways) sum(pq[pathways,"HCC1954_BREAST"], na.rm=T))
ecm_param$ECMpathway_RNAseq[!is.finite(ecm_param$ECMpathway_RNAseq)] = NA
ecm_param$pathways = sapply(tmp[rownames(ecm_param)], paste, collapse=", ")
write.table(ecm_param, paste0(OUTDIR, filesep,"BestFits_per_ECM+pathways.txt"), sep="\t", quote=F, row.names = F);
## Visualize and get statistics
for(rna in c("ECM_RNAseq_Expression","ECMpathway_RNAseq")){
  cc = sapply(coi, function(ii) cor.test(as.numeric(ecm_param[,ii]),ecm_param[, rna], na.rm=T)[c("estimate","p.value")])
  la=sapply(coi, function(ii) plot(as.numeric(ecm_param[,ii]),as.numeric(ecm_param[, rna]), pch=20, xlab = ii, ylab=rna))
  print(rna)
  print(cc)
}
p <- ggplot(ecm_param,aes(xi_u, ECMpathway_RNAseq, label= rownames(ecm_param))) +  geom_point() + labs(x=expression(xi), y="RNA-seq derived pathway activity involving ECM") +  geom_smooth(method='lm', formula= y~x)  + geom_text(size = 2, hjust=0.5, vjust=-0.5) +
  ggtitle(paste("r =", round(unlist(cc["estimate","xi_u"]),3),"; p =", round(unlist(cc["p.value","xi_u"]),5))) 
ggsave(paste0(OUTDIR, filesep,"ECMspecificCalibration_vs_RNAseq.pdf"), plot = p, width = 4.05, height = 3.95)
## Transcriptome subtypes
st = read.table("~/Projects/Resources/data/Databases/MEP_LINCS/JWGray_BCCL_classifications_v5.txt", sep="\t", check.names = F, stringsAsFactors = F, header=T)
rownames(st) = paste0(st$`Cell Line`,"_BREAST")
missing = setdiff(names(ploidy),rownames(st))
missing = data.frame("Cell Line" = missing, "Classification-3class" = c("Luminal","Basal"), "Classification-4class" = rep(NA, length(missing)), check.names = F)
rownames(missing) = missing$`Cell Line`
st=rbind(st, missing)
coi = intersect(names(ploidy),rownames(st))
data = cbind(as.data.frame(ploidy[coi]), st[coi,]$`Classification-3class`)
colnames(data) = c("ploidy", "Subtype")
stat_box_data <- function(y, upper_limit = max(data$ploidy) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('n =', length(y))
    )
  )
}
p <- ggplot(data, aes(factor(Subtype), ploidy, fill=T)) + 
  geom_violin() + geom_jitter(height = 0, width = 0.1) + xlab("") + stat_summary(
  fun.data = stat_box_data, 
  geom = "text", 
  hjust = 0.5,
  vjust = 0.9
) 
ggsave(paste0(OUTDIR, filesep,"BreastCancerSubtype_Ploidy.pdf"), plot = p, width = 3.65, height = 3)



## Drug sensitivity: GDSC (Genomics of Drug Sensitivity in Cancer) database
f = list.files("~/Projects/Resources/data/Databases/CCLE/BreastCancerDrugSensitivity/", full.names = T)
#Drug sensitivities for each CL
clin  = list()
for(x in f){
  sName = strsplit(fileparts(x)$name,"_")[[1]][1]
  tmp = read.table(x, sep="\t", header=T)
  tmp = grpstats(as.matrix(tmp$Z.Score), as.character(tmp$Drug.Name), "mean")$mean
  ix = sort(rownames(tmp), index.return=T)$ix
  clin[[sName]] = tmp[ix,, drop=F]
}
## Count drugs:
fr = plyr::count(unlist(sapply(clin,rownames)))
fr = fr[sort(fr$freq, index.return=T, decreasing = T)$ix,]
doi = as.character(fr$x[fr$freq>=0.8*max(fr$freq)])
#Rows = drugs; Columns = CLs; entries = z-score
clin = sapply(clin, function(x) as.data.frame(x)[doi,])
rownames(clin) = doi
colnames(clin) = paste0(colnames(clin), "_BREAST")
cc = sapply(rownames(clin), function(x) cor.test(ploidy[colnames(clin)],clin[x,], na.rm=T)[c("estimate", "p.value")])
## Annotate
cc = as.data.frame(t(cc))
rownames(cc) = gsub(".cor","",rownames(cc))
anno = read.table("~/Projects/Resources/data/Databases/CCLE/DrugAliases.txt", sep="\t", check.names = F, stringsAsFactors = F, header=T)
ix = sapply(rownames(cc), function(x) unique(c(grep(x, anno$drug_name, fixed = T), grep(x,anno$synonyms, fixed = T), grep(x,gsub("-","", anno$synonyms, fixed = T) )))[1] )
cc$pathway_name = anno$pathway_name[unlist(ix)]
cc$estimate = as.numeric(cc$estimate)
cc=cc[sort(cc$estimate, index.return=T)$ix,]
cc = cc[abs(cc$estimate)>=0.2,]
## Visualize
par(mfrow=c(1,2), mai=c(1,1,0.04,0.9))
b = barplot(cc$estimate, horiz=T, names=gsub("Inhibitor","I.",rownames(cc)), las=2, xlab="Pearson (ploidy, IC50)", col=c("purple","orange")[1+(cc$estimate>0)], border="white", cex.names = 0.7, cex.lab=0.8, cex.axis=0.7)
axis(side = 4, labels=gsub(" and ","/",cc$pathway_name), las=2, at = b[,1], cex.axis=0.7)
plot(1,xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
legend("topright",c("Low ploidy is sensitive","High ploidy is sensitive"), fill=c("orange","purple"),bty="n", cex=0.69)


## Predict drug sensitivity from interaction between drug category and ploidy
library(survival)
anno$category = NA
anno$category[anno$pathway_name %in% c("DNA replication", "p53 pathway", "Genome integrity", "Mitosis","Cell cycle")] = "Cytotoxic"
anno$category[grep("signaling",anno$pathway_name)] = "Signaling"
coxIn = as.data.frame(matrix(NA,nrow(clin)*ncol(clin), 6));
colnames(coxIn) = c("drugName", "drugCategory", "ploidy", "DoublingTime", "IC50", "Subtype")
rownames(coxIn)=sapply(colnames(clin), function(x) paste(rownames(clin),x))
for(sName in colnames(clin)){
  for(drug in rownames(clin)){
    ij = paste(drug,sName)
    coxIn[ij,"drugName"] = drug
    coxIn[ij,"ploidy"] = ploidy[sName]
    coxIn[ij,"IC50"] = clin[drug,sName] 
    coxIn[ij,"Subtype"] = st[sName,"Classification-3class"] 
    coxIn[ij,"DoublingTime"] = appCL$`Cellosaurus Doubling time`[match(sName, appCL$`CCLE name`)]
    ii = unlist(sapply(c("drug_name","synonyms","pathway_name","targets"), function(x) grep(drug, anno[,x])))
    if(!isempty(ii)){
      coxIn[ij,"drugCategory"] = unique(anno[ii,"category"])
    }
  }
}
coxIn = coxIn[!is.na(coxIn$drugCategory),]
coxIn$drugCategory = as.factor(coxIn$drugCategory)
## Visualize and get statistics
for(X in c("ploidy", "DoublingTime")){
  coxIn$X = coxIn[,X]
  m = summary(lm(IC50 ~ Subtype + X * drugCategory, data = coxIn[coxIn$Subtype!="LA",])); ##Basal
  print(m)
  write.table(m$coefficients, file=paste0(OUTDIR, filesep, "X=",X ,"_DrugCategory_IC50.txt"), sep="\t", quote = F)
  write.table(paste("adj.r.squared:",m$adj.r.squared), file=paste0(OUTDIR, filesep, "X=",X ,"_DrugCategory_IC50.txt"), sep="\t", quote = F, append = T, row.names = F, col.names = F)
  
  
  plot.interact <- qplot(x = X, y = IC50, color = drugCategory, data = coxIn[coxIn$Subtype!="LA",]) + #Basal
    stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, size=3) + xlab(X)  +
    theme(axis.text.x = element_text(face="plain", size = 16),axis.text.y = element_text(face="plain", size = 16),axis.title.x = element_text(face="plain", size = 15.5),axis.title.y = element_text(face="plain", size = 15.5))
  ggsave(paste0(OUTDIR, filesep, "X=",X ,"_DrugCategory_IC50.pdf"), plot = plot.interact, width = 5.75, height = 6.65)
}

print(plyr::count(coxIn$drugCategory[!duplicated(coxIn$drugName)]))




