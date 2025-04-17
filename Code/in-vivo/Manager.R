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

setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerOrthotopicModels/SUM-159")



#############################
#### scRNA-seq results ######
f=list.files("A03_Numbat")#, full.names = T)#, pattern = "tsv", recursive = T)
f=grep("15",f,invert = T, value = T)
f=grep("_Numbat",f,invert = T, value = T)
cn=NumbatPostProcess(DATASETID="A03_Numbat/", mpoi=f, path2karyo="./")

batch="C2N_chr2"
origin=unique(cn[[batch]]$cells)
col =RColorBrewer::brewer.pal(length(origin),"Paired")
names(col) = origin
hm=heatmap.2(cn[[batch]]$cn, Colv = NULL, trace = "n", RowSideColors = col[cn[[batch]]$cells])
legend("topright",names(col),fill = col,cex=0.5)

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
    ploidy[id]=sum(cell$Chromosome.Area)
  }
  heatmap.2(karyo,trace='n',main = s)
  # rownames(karyo)= formatNames(rownames(karyo))
  # names(ploidy)= formatNames(names(ploidy))
  samples[[s]]=list(karyo=karyo,ploidy=ploidy, numchrcopies=numchrcopies)
}
whole_chr_karyo =do.call(rbind,sapply(samples, function(x) x$karyo, simplify = F))
write.table(whole_chr_karyo, file="~/Downloads/whole_chr_karyo.txt", sep="\t", quote=F)


#######################################
#### Merge scRNA-seq & Karyotyping ####
# cn=cn[-2]
# subset="none"
subset="2N"
cn_=cn;
whole_chr_karyo_ = whole_chr_karyo;
if(subset!="none"){
  cn_ = sapply(cn, function(x) list(cn=x$cn[grep(subset,x$cells),]), simplify = F)
  cn_[[1]]$cells = grep(subset,cn[[1]]$cells, value=T)
  cn_[[2]]$cells = grep(subset,cn[[2]]$cells, value=T)
  whole_chr_karyo_=whole_chr_karyo[grep(subset,rownames(whole_chr_karyo_)),]
}

whole_chr_scRNA=do.call(rbind,sapply(cn_, function(x) x$cn[grep("Cell-Culture",x$cells),], simplify = F))
if(subset=="4N"){
  ## Numbat output needs to be adjusted for tetraploidy
  whole_chr_scRNA = whole_chr_scRNA * 2
}
# pdf(paste0("~/Downloads/",subset,".pdf"))

cells=sapply(cn_, function(x) grep("Cell-Culture",x$cells, value=T), simplify = F)
cells=unlist(cells)
origin=unique(cells)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = origin
hm_s=heatmap.2(whole_chr_scRNA, Colv = NULL, trace = "n", RowSideColors =col[cells], hclustfun=function(x) hclust(x, method="ward.D2"),distfun=function(x) dist(x, method="manhattan"))
legend("topright",names(col),fill=col, cex=0.65)
cl_s=cutree(as.hclust(hm_s$rowDendrogram), k=8)
fr=plyr::count(cl_s);
fr$freq=round(fr$freq/nrow(whole_chr_scRNA),3)
origin=unique(cl_s)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = as.character(origin)
heatmap.2(whole_chr_scRNA, Colv = NULL, trace = "n", RowSideColors =col[as.character(cl_s)], hclustfun=function(x) hclust(x, method="ward.D2"),distfun=function(x) dist(x, method="manhattan"))
cl_s = grpstats(whole_chr_scRNA,cl_s,"mean")$mean
rownames(cl_s) = paste0(rownames(cl_s),"_", fr[rownames(cl_s),"freq"])

cells=sapply(strsplit(rownames(whole_chr_karyo_),"_"),"[[",1)
origin=unique(cells)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = origin
hm_k=heatmap.2(whole_chr_karyo_, Colv = NULL, trace = "n", RowSideColors =col[cells], hclustfun=function(x) hclust(x, method="ward.D2"),distfun=function(x) dist(x, method="manhattan"))
legend("topright",names(col),fill=col, cex=0.65)
cl_k=cutree(as.hclust(hm_k$rowDendrogram), k=8)
heatmap.2(whole_chr_karyo_, Colv = NULL, trace = "n", RowSideColors =col[cells], colRow =  cl_k, hclustfun=function(x) hclust(x, method="ward.D2"),distfun=function(x) dist(x, method="manhattan"))
legend("topright",names(col),fill=col, cex=0.65)
cl_k = grpstats(whole_chr_karyo_,cl_k,"mean")$mean

cn_comb = rbind(cl_k,cl_s)
# cn_comb = cn_comb[,-5]
rownames(cn_comb) = paste0(rownames(cn_comb),"_",c(rep("karyo",nrow(cl_k)),rep("sc",nrow(cl_s))))
hm=heatmap.2(cn_comb, Colv = NULL, trace = "n", hclustfun=function(x) hclust(x, method="ward.D2"),distfun=function(x) dist(x, method="manhattan"))
cl=cutree(as.hclust(hm$rowDendrogram), k=8)
origin=unique(cl)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = as.character(origin)
heatmap.2(cn_comb, Colv = NULL, trace = "n", RowSideColors =col[as.character(cl)], hclustfun=function(x) hclust(x, method="ward.D2"),distfun=function(x) dist(x, method="manhattan"), main=subset)

dev.off()

##chatGPT: define reusable functions
## conclusions: refine chr arm level calls for 2N
## rerun numbat for 4N.
## annotate cell representation or % on combined heatmap