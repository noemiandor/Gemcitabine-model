library(xlsx)
library(gplots)
library(data.table)
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
source("~/Repositories/Gemcitabine-model/Code/in-vivo/Utils.R")
#################################
## Expected chromosome lengths ##
x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
chrarms=x[ , .(length = sum(chromEnd - chromStart)),by = .(chrom, arm = substring(name, 1, 1)) ]
chrwhole=grpstats(as.matrix(chrarms$length),chrarms$chrom, "sum")$sum
rownames(chrwhole)=gsub("chrY","chr24",gsub("chrX","chr23",rownames(chrwhole)))

# weighted Manhattan distance for an *entire* matrix
chrWeightedDist <- function(mat) {
  # vector of chromosome weights
  w <- chrwhole[paste0("chr", 1:22), 1]
  mat.w <- sweep(mat, 2, w, `*`)         # weight every column
  dist(mat.w, method = "manhattan") / sum(w)
}

# correlation‑based distance (no chromosome weighting)
chrCorrDist <- function(mat) {
  ## Distance = 1 – Pearson correlation
  ##  (values lie in [0,2]; divide by 2 if you prefer [0,1])
  as.dist(1 - cor(t(mat), use = "pairwise.complete.obs", method = "pearson"))
}

# correlation‑based distance **with chromosome‑length weights**
chrWeightedCorrDist <- function(mat) {
  # chromosome weights (lengths in 'chrwhole')
  w <- chrwhole[paste0("chr", 1:22), 1]
  w <- w / sum(w)                       # normalise so Σw = 1  (optional but convenient)
  
  # ----- internal: weighted Pearson correlation for two vectors -----
  wcor <- function(x, y) {
    mx <- sum(w * x)
    my <- sum(w * y)
    cov_xy <- sum(w * (x - mx) * (y - my))
    sx <- sqrt(sum(w * (x - mx)^2))
    sy <- sqrt(sum(w * (y - my)^2))
    if (sx == 0 || sy == 0) return(0)   # guard against zero variance
    cov_xy / (sx * sy)
  }
  
  # ----- pairwise distance matrix -----
  n <- nrow(mat)
  D <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      r <- wcor(mat[i, ], mat[j, ])
      D[i, j] <- D[j, i] <- 1 - r       # distance = 1 – correlation
    }
  }
  as.dist(D)
}



setwd("~/Repositories/Gemcitabine-model/Code/in-vivo")
source("Utils.R")
dt=read.xlsx("../../Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx", sheetIndex = 1)
dt <- rbind(setNames(data.frame(matrix(NA, nrow = 2, ncol = ncol(dt))), names(dt)),dt)
dt$harvest[1:2] = c("SUM-159_NLS_2N_A7M_K_harvest","SUM-159_NLS_4N_A5M_K_harvest")
dt$Sequencing.IDs[1:2] = c("2N-Cell-Culture","4N-Cell-Culture")
dt=dt[!is.na(dt$Sequencing.IDs),];
rownames(dt)=dt$Sequencing.IDs
## Path on workstation is:
# setwd("/mnt/ix1/Shared_Folders/lab_crd/HighPloidy_CostBenefits/data/BreastCancerOrthotopicModels/SUM-159")
setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerOrthotopicModels/SUM-159")

#############################
#### scRNA-seq results ######
f=list.files("A03_Numbat")#, full.names = T)#, pattern = "tsv", recursive = T)
f=grep("15",f,invert = T, value = T)
f=grep("Cell-Culture",f,invert = T, value = T)
f=grep("_Numbat",f,invert = T, value = T)
cn=NumbatPostProcess(DATASETID="A03_Numbat/", mpoi=f, path2karyo="./")

numbatRun="C4N_chr13"
origin=unique(cn[[numbatRun]]$cells)
col =RColorBrewer::brewer.pal(length(origin),"Paired")
names(col) = origin
hm=heatmap.2(cn[[numbatRun]]$cn, Colv = NULL, trace = "n", RowSideColors = col[cn[[numbatRun]]$cells])
legend("topright",names(col),fill = col,cex=0.5)
print(dt[origin,1:8])

#############################
#### Karyotyping results ####
f = list.files("B02_Karyotyping/", pattern = ".csv", recursive = T, full.names = T)
kn=sapply(f, read.csv, simplify = F)
kn=do.call(rbind,kn)
kn$Images.Name=paste0(kn$Dataset.Name,"_",kn$Images.Name)
cells=unique(kn$Images.Name)
chr= 1:22
# cells= sapply(cells, function(x) kn[kn$Images.Name==x,], simplify = F)
# fusions <- arm_level_karyo <-list()
# coi=c('Chromosome.Upper.Arm.area','Chromosome.Lower.Arm.area')
samples= sapply( unique(kn$Dataset.Name), function(x) kn[kn$Dataset.Name==x,], simplify = F)

for(s in names(samples)){
  lineage=samples[[s]]
  cells=sort(unique(lineage$Images.Name))
  ploidy <- numchrcopies <-rep(NA,length(cells))
  karyo=matrix(0, length(cells), length(chr))
  names(ploidy) <- rownames(karyo) <- cells
  colnames(karyo)=as.character(chr)
  for(id in cells){
    cell= lineage[lineage$Images.Name==id,]
    autosomes = cell[cell$Group.ID<=22,]
    fr=plyr::count(autosomes$Group.ID)
    karyo[id,as.character(fr$x)]=fr$freq
    # ploidy[id]=sum(cell$Chromosome.Length)
    numchrcopies[id]=nrow(cell)
    # ploidy[id]=sum(cell$Chromosome.Area)
    ploidy[id]=sum(chrwhole[paste0("chr",cell$Group.ID),])/sum(chrwhole)
  }
  heatmap.2(karyo,trace='n',main = s)
  # rownames(karyo)= formatNames(rownames(karyo))
  # names(ploidy)= formatNames(names(ploidy))
  samples[[s]]=list(karyo=karyo,ploidy=ploidy, numchrcopies=numchrcopies)
}
ploidy =sapply(samples, function(x) x$ploidy, simplify = F)
whole_chr_karyo =do.call(rbind,sapply(samples, function(x) x$karyo, simplify = F))
write.table(whole_chr_karyo, file="~/Downloads/whole_chr_karyo.txt", sep="\t", quote=F)
par(mai=c(0.5,2,0.5,0.5)); boxplot(ploidy[-c(1:2)],las=2, horizontal = T)

###############################
#### Karyotyping arm level ####
f = list.files("B02_Karyotyping/", pattern = "ArmLevel.xlsx", recursive = T, full.names = T)
kan=sapply(f, function(x) read.xlsx(x,sheetIndex =1, check.names=F), simplify = F)
whole_chrarm_karyo=do.call(rbind,kan)
# whole_chrarm_karyo = read.xlsx("B02_Karyotyping/SUM159-4N-parental/SUM159_4N_Karyotyping_ArmLevel.xlsx",sheetIndex =1, check.names=F);
whole_chrarm_karyo=whole_chrarm_karyo[apply(!is.na(whole_chrarm_karyo),1,all),]
ii = which(!colnames(whole_chrarm_karyo) %in% c("Cell","MARKER"))
heatmap.2(as.matrix(whole_chrarm_karyo[,ii]),trace='n')

#######################################
#### Merge scRNA-seq & Karyotyping ####
numbat2keep=c("C4N_chr11","C2N_chr23")
cn_B=cn; ## make a copy
cn = cn_B[numbat2keep]
# cn=cn[-2]
# subset="none"
subset="4N"
whole_chr_karyo_ = whole_chr_karyo;
whole_chrarm_karyo_ = whole_chrarm_karyo
if(subset!="none"){
  ii=grep(subset,names(cn));
  cn_ = sapply(cn[ii], function(x) list(cn=x$cn[grep(subset,x$cells),]), simplify = F)
  cn_[[1]]$cells = grep(subset,cn[[ii]]$cells, value=T)
  cn_[[1]]$anno = cn[[ii]]$anno
  whole_chr_karyo_=whole_chr_karyo[grep(subset,rownames(whole_chr_karyo)),]
  whole_chrarm_karyo_ = whole_chrarm_karyo[grep(subset, rownames(whole_chrarm_karyo)),]
}
## Align segments
la=alignCNmatrices(whole_chrarm_karyo_, cn_[[1]])
whole_chr_karyo_ = la$cn_karyo[whole_chrarm_karyo_$MARKER<10,]
cn_[[1]]$cn = la$cn_scRNAseq

## select cells of interest
whole_chr_scRNA=do.call(rbind,sapply(cn_, function(x) x$cn[grep("Cell-Culture",x$cells),], simplify = F))
if(subset=="4N"){
  ## Numbat output needs to be adjusted for tetraploidy
  whole_chr_scRNA = whole_chr_scRNA * 2
}

## reorder chr (arms)
ii=sort(colnames(whole_chr_karyo_))
whole_chr_karyo_ = whole_chr_karyo_[,ii]
whole_chr_scRNA = whole_chr_scRNA[,ii]


mar=c(20,5)
dfun = chrCorrDist
# dfun = function(x) dist(x, method="manhattan")
pdf(paste0("~/Downloads/",subset,".pdf"))
## scRNAseq
cells=sapply(cn_, function(x) grep("Cell-Culture",x$cells, value=T), simplify = F)
cells=unlist(cells)
origin=unique(cells)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = origin
hm_s=heatmap.2(whole_chr_scRNA, Colv = NULL, trace = "n", RowSideColors =col[cells], hclustfun=function(x) hclust(x, method="ward.D2"),distfun=dfun, margins = mar) 
legend("topright",names(col),fill=col, cex=0.65)
cl_s=cutree(as.hclust(hm_s$rowDendrogram), k=8)
fr=plyr::count(cl_s);
fr$freq=round(fr$freq/nrow(whole_chr_scRNA),3)
origin=unique(cl_s)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = as.character(origin)
heatmap.2(whole_chr_scRNA, Colv = NULL, trace = "n", RowSideColors =col[as.character(cl_s)], hclustfun=function(x) hclust(x, method="ward.D2"),distfun=dfun, margins = mar) 
cl_s = grpstats(whole_chr_scRNA,cl_s,"mean")$mean
rownames(cl_s) = paste0(rownames(cl_s),"_", fr[rownames(cl_s),"freq"])

## Karyotyping
cells=sapply(strsplit(rownames(whole_chr_karyo_),"_"),"[[",1)
origin=unique(cells)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = origin
hm_k=heatmap.2(whole_chr_karyo_, Colv = NULL, trace = "n", RowSideColors =col[cells],distfun=chrWeightedDist, margins = mar) 
legend("topright",names(col),fill=col, cex=0.65)
cl_k=cutree(as.hclust(hm_k$rowDendrogram), k=8)
heatmap.2(whole_chr_karyo_, Colv = NULL, trace = "n", RowSideColors =col[cells], colRow =  cl_k,distfun=chrWeightedDist, margins = mar) 
legend("topright",names(col),fill=col, cex=0.65)
cl_k = grpstats(whole_chr_karyo_,cl_k,"mean")$mean


cn_comb = rbind(cl_k,cl_s)
# cn_comb = cn_comb[,-5]
rownames(cn_comb) = paste0(rownames(cn_comb),"_",c(rep("karyo",nrow(cl_k)),rep("sc",nrow(cl_s))))
hm=heatmap.2(cn_comb, Colv = NULL, trace = "n", hclustfun=function(x) hclust(x, method="ward.D2"),distfun=dfun, margins = mar) 
cl=cutree(as.hclust(hm$rowDendrogram), k=8)
origin=unique(cl)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = as.character(origin)
heatmap.2(cn_comb, Colv = NULL, trace = "n", RowSideColors =col[as.character(cl)], hclustfun=function(x) hclust(x, method="ward.D2"),distfun=dfun, main=subset, margins = mar) 

dev.off()

##chatGPT: define reusable functions
## conclusions: refine chr arm level calls for 2N
## rerun numbat for 4N.
## annotate cell representation or % on combined heatmap