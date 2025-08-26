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

# ===================== clustering-by-origin utilities ======================

# tiny entropy helper (natural log)
.entropy <- function(p) { p <- p[p > 0]; if (!length(p)) return(0); -sum(p * log(p)) }

# score a clustering against sample labels
.score_cut <- function(labels, clusters) {
  S <- factor(labels)
  C <- factor(clusters)
  tab <- table(S, C)
  P   <- prop.table(tab)                 # joint p(s,c)
  ps  <- rowSums(P)                      # p(s)
  pc  <- colSums(P)                      # p(c)
  
  # Mutual information
  MI <- sum(ifelse(P > 0, P * log(P / (ps %o% pc)), 0))
  HS <- .entropy(ps)
  HC <- .entropy(pc)
  
  # Normalized MI (symmetric): balances intra- and inter- variability
  NMI <- if (HS > 0 && HC > 0) MI / sqrt(HS * HC) else 0
  
  # Weighted within-cluster entropy of sample origin
  w_intra_entropy <- sum(pc * apply(tab, 2, function(col) .entropy(col / sum(col))))
  
  data.frame(N = length(S), K = nlevels(C), MI = MI, HS = HS, HC = HC,
             NMI = NMI, IntraEntropy = w_intra_entropy, check.names = FALSE)
}

# choose best cut across a range of K by maximizing NMI (tie-break: min IntraEntropy)
choose_best_cut <- function(hc_or_dend, sample_labels, k_range = 2:20) {
  hc <- if (inherits(hc_or_dend, "dendrogram")) stats::as.hclust(hc_or_dend) else hc_or_dend
  
  res <- lapply(k_range, function(k) {
    cl <- cutree(hc, k = k)
    cbind(K = k, .score_cut(sample_labels, cl))
  })
  scores <- do.call(rbind, res)
  
  # best by NMI; tie-break with lower IntraEntropy, then fewer clusters
  ord <- with(scores, order(-NMI, IntraEntropy, K))
  best_row <- scores[ord[1], ]
  best_k   <- best_row$K
  best_cl  <- cutree(hc, k = best_k)
  
  list(k = best_k, clusters = best_cl, scores = scores[order(scores$K), ])
}

# ------------------------------------------------------------
# Barplot B: distribution of each sample across clusters
#   - one bar per sample
#   - stacks = clusters (proportion within each sample)
#   - bars ordered by diversity (Shannon entropy of cluster composition)
# ------------------------------------------------------------
plot_sample_distribution <- function(clusters, sample_labels, cluster_colors = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  stopifnot(length(clusters) == length(sample_labels))
  
  df <- data.frame(sample = factor(sample_labels), cluster = factor(clusters))
  prop <- as.data.frame(prop.table(table(df$sample, df$cluster), 1))
  colnames(prop) <- c("sample", "cluster", "prop")
  
  # compute entropy per sample
  ent <- tapply(prop$prop, prop$sample, .entropy)
  prop$sample <- factor(prop$sample, levels = names(sort(ent, decreasing = TRUE)))
  
  # optional: order clusters by global size
  cl_sizes <- sort(table(df$cluster), decreasing = TRUE)
  prop$cluster <- factor(prop$cluster, levels = names(cl_sizes))
  
  p <- ggplot2::ggplot(prop, ggplot2::aes(x = sample, y = prop, fill = cluster)) +
    ggplot2::geom_bar(stat = "identity", width = 0.85) +
    ggplot2::ylab("Proportion of cells in sample") +
    ggplot2::xlab("Sample (ordered by diversity)") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  if (!is.null(cluster_colors)) {
    lev <- levels(prop$cluster)
    p <- p + ggplot2::scale_fill_manual(values = cluster_colors[lev], drop = FALSE)
  }
  p
}


# ------------------------------------------------------------
# Barplot: distribution of each sample across clusters
#   - one bar per sample
#   - stacks are clusters (proportion within each sample)
# ------------------------------------------------------------
plot_sample_distribution <- function(clusters, sample_labels) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2")
  }
  stopifnot(length(clusters) == length(sample_labels))

  df <- data.frame(
    sample  = factor(sample_labels),
    cluster = factor(clusters)
  )

  # proportions of clusters within each sample (rows sum to 1)
  prop <- as.data.frame(prop.table(table(df$sample, df$cluster), 1))
  colnames(prop) <- c("sample", "cluster", "prop")

  # order samples by their dominant cluster proportion (just for a nicer look)
  dom <- aggregate(prop ~ sample, prop[ave(prop$prop, prop$sample, FUN = max) == prop$prop, ], max)
  prop$sample <- factor(prop$sample, levels = dom$sample[order(-dom$prop)])

  # optional: order clusters by global size
  cl_sizes <- sort(table(df$cluster), decreasing = TRUE)
  prop$cluster <- factor(prop$cluster, levels = names(cl_sizes))

  ggplot2::ggplot(prop, ggplot2::aes(x = sample, y = prop, fill = cluster)) +
    ggplot2::geom_bar(stat = "identity", width = 0.85) +
    ggplot2::ylab("Proportion of cells in sample") +
    ggplot2::xlab("Sample") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "right")
}

# stacked bar plot: one bar per cluster, fill by sample origin proportion
plot_cluster_composition <- function(clusters, sample_labels, sample_colors = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  stopifnot(length(clusters) == length(sample_labels))
  
  df <- data.frame(cluster = factor(clusters), sample = factor(sample_labels))
  prop <- as.data.frame(prop.table(table(df$cluster, df$sample), 1))
  colnames(prop) <- c("cluster", "sample", "prop")
  
  # compute entropy per cluster
  ent <- tapply(prop$prop, prop$cluster, .entropy)
  prop$cluster <- factor(prop$cluster, levels = names(sort(ent, decreasing = TRUE)))
  
  p <- ggplot2::ggplot(prop, ggplot2::aes(x = cluster, y = prop, fill = sample)) +
    ggplot2::geom_bar(stat = "identity", width = 0.85) +
    ggplot2::ylab("Proportion of cells in cluster") +
    ggplot2::xlab("Cluster (ordered by diversity)") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::theme_classic(base_size = 12)
  
  if (!is.null(sample_colors)) {
    lev <- levels(prop$sample)
    p <- p + ggplot2::scale_fill_manual(values = sample_colors[lev], drop = FALSE)
  }
  p
}

plot_sample_distribution <- function(clusters, sample_labels, cluster_colors = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  df <- data.frame(sample = factor(sample_labels), cluster = factor(clusters))
  prop <- as.data.frame(prop.table(table(df$sample, df$cluster), 1))
  colnames(prop) <- c("sample", "cluster", "prop")
  
  gg <- ggplot2::ggplot(prop, ggplot2::aes(x = sample, y = prop, fill = cluster)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::ylab("Proportion of cells in sample") +
    ggplot2::xlab("Sample") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  if (!is.null(cluster_colors)) {
    gg <- gg + ggplot2::scale_fill_manual(values = cluster_colors)
  }
  gg
}

# weighted Manhattan distance for an *entire* matrix
chrWeightedDist <- function(mat, w=NULL) {
  # vector of chromosome weights
  if(is.null(w)){
    w <- chrwhole[paste0("chr", 1:22), 1]
  }
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
#### Karyotyping results ####
f = list.files("B02_Karyotyping/", pattern = ".csv", recursive = T, full.names = T)
f=grep("Labcorp",f,invert = T,value = T)
f=grep("parental",f,invert = T,value = T)
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
par(mai=c(0.5,2,0.5,0.5)); boxplot(ploidy,las=2, horizontal = T)
N2_whole_chr_karyo = whole_chr_karyo[grep("N2",rownames(whole_chr_karyo)),]
N4_whole_chr_karyo = whole_chr_karyo[grep("N4",rownames(whole_chr_karyo)),]

# ###############################
# #### Karyotyping arm level ####
# f = list.files("B02_Karyotyping/", pattern = "ArmLevel.xlsx", recursive = T, full.names = T)
# kan=sapply(f, function(x) read.xlsx(x,sheetIndex =1, check.names=F), simplify = F)
# whole_chrarm_karyo=do.call(rbind,kan)
# # whole_chrarm_karyo = read.xlsx("B02_Karyotyping/SUM159-4N-parental/SUM159_4N_Karyotyping_ArmLevel.xlsx",sheetIndex =1, check.names=F);
# whole_chrarm_karyo=whole_chrarm_karyo[apply(!is.na(whole_chrarm_karyo),1,all),]
# ii = which(!colnames(whole_chrarm_karyo) %in% c("Cell","MARKER"))
# heatmap.2(as.matrix(whole_chrarm_karyo[,ii]),trace='n')


#############################
#### scRNA-seq results ######
f=list.files("A03_Numbat")#, full.names = T)#, pattern = "tsv", recursive = T)
f=grep("15",f,invert = T, value = T)
f=grep("Cell-Culture",f,invert = T, value = T)
f=grep("_Numbat",f,invert = T, value = T)
f_2N=grep("C2N",f,invert = F, value = T)
cn=NumbatPostProcess(DATASETID="A03_Numbat/", mpoi=f_2N, path2karyo="./", gBandedKaryo2align=N2_whole_chr_karyo, scRNAseqCells2align="2N-Cell-Culture")
cn_B=cn; ## make a copy
f_4N=grep("C4N_chr18",f,invert = F, value = T)
f_4N=grep("C4N_chr11",f_4N,invert = T, value = T)
cn=NumbatPostProcess(DATASETID="A03_Numbat/", mpoi=f_4N, path2karyo="./", mergeRun2C4N.tumor_chr6 = T, gBandedKaryo2align=N4_whole_chr_karyo, scRNAseqCells2align="4N-Cell-Culture", lambda_conflict=0)
# cn_=NumbatPostProcess(DATASETID="A03_Numbat/", mpoi="C4N.tumor_chr6", path2karyo="./")
cn_B[names(cn)] = cn


numbatRun="C4N_chr18"
# numbatRun="C4N.tumor_chr6"
# numbatRun="C2N_chr2"
origin=unique(cn[[numbatRun]]$cells)
col =RColorBrewer::brewer.pal(length(origin),"Paired")
names(col) = origin
hm=heatmap.2(cn[[numbatRun]]$cn, Colv = NULL, trace = "n", RowSideColors = col[cn[[numbatRun]]$cells], hclustfun = function(x) hclust(x, method = "ward.D"), margins = c(15,5), col=(rainbow(8))[1:6])
legend("topright",names(col),fill = col,cex=0.75)
print(dt[origin,1:8])
## ploidy
tmp=colnames(cn[[numbatRun]]$cn)
loci=parseLOCUS(sapply(strsplit(tmp,"_"),"[[",1))
segweight = loci[,"seglength"]/sum(loci[,"seglength"])
cn_wholechr = grpstats(t(cn[[numbatRun]]$cn), loci[,"chr"],"mean")$mean
ploidy = apply(sweep(cn[[numbatRun]]$cn,2,segweight,"*"),1,sum, na.rm=T)
ploidy = apply(cn_wholechr,2,sum, na.rm=T)
vioplot::vioplot(ploidy~cn[[numbatRun]]$cells,las=2,horizontal=F, xlab="")
vioplot::vioplot(ploidy~as.numeric(cn[[numbatRun]]$cells=="4N-Cell-Culture"),las=2,horizontal=F, xlab="")

## Cluster and enrichment analysis:
# D <- chrWeightedDist(cn[[numbatRun]]$cn, loci[,"seglength"])
D <- dist(cn[[numbatRun]]$cn, method = "manhattan")
dend <- hclust(D, method = "ward.D")
sample_origin <- cn[[numbatRun]]$cells[dend$order]
best <- choose_best_cut(dend, sample_labels = sample_origin,
                        k_range = 7:min(7, nrow(cn[[numbatRun]]$cn) - 1))
head(best$scores[order(-best$scores$NMI), ], 5)  # top scoring cuts
# -----------------------------------------------------------------
# Create a consistent cluster color palette
# -----------------------------------------------------------------
uniq_clusters <- sort(unique(best$clusters))
cluster_colors <- setNames(rainbow(length(uniq_clusters)), uniq_clusters)
hm <- gplots::heatmap.2(
  cn[[numbatRun]]$cn,
  Rowv = as.dendrogram(dend),    # or just Rowv = dend
  Colv = NA,                     # don't cluster columns
  dendrogram = "row",            # draw only the row dendrogram
  reorderfun = function(d, w) d, # keep *exact* input dendrogram order
  trace = "n",
  RowSideColors = col[cn[[numbatRun]]$cells],
  margins = c(15, 5),
  col = (rainbow(8))[1:6],
  colRow = cluster_colors[best$clusters[names(cn[[numbatRun]]$cells)]]
)
legend("topright",names(col),fill = col,cex=0.75)
p1 <- plot_cluster_composition(best$clusters, sample_origin[names(best$clusters)], col)
p2 <- plot_sample_distribution(best$clusters, sample_origin[names(best$clusters)], cluster_colors)
print(p1)
print(p2)
# Fisher's test:
enr <- fisher_cluster_enrichment(best$clusters, sample_origin[names(best$clusters)])
enr$padj_table
plot_enrichment_heatmap(enr, cluster_rows = FALSE, cluster_cols = TRUE)

# ## LIAYSON ##
# library(Seurat)
# library(liayson)
# library(matlab);
# library(RColorBrewer)
# HOST="https://may2021.archive.ensembl.org"
#   
# anno=read.table("A03_Numbat/C4N_chr6/C4N_chr6.cell4numbat.anno.txt", header = T);
# expression_data <- Read10X(data.dir = "A02_cellRanger/4N-Cell-Culture-Count-HM/outs/filtered_feature_bc_matrix")
# cells <- anno$cell[anno$sample=="4N-Cell-Culture"]
# object1 <- CreateSeuratObject(counts = expression_data, project = "MyProject", min.cells = 3, min.features = 200)
# object1 = object1[,cells]
# 
# anno=read.table("A03_Numbat/C2N_chr2/C2N_chr2.cell4numbat.anno.txt", header = T);
# expression_data <- Read10X(data.dir = "A02_cellRanger/2N-Cell-Culture-Count-HM/outs/filtered_feature_bc_matrix")
# cells <- anno$cell[anno$sample=="2N-Cell-Culture"]
# object2 <- CreateSeuratObject(counts = expression_data, project = "MyProject", min.cells = 3, min.features = 200)
# object2 = object2[,cells]
# 
# # Merge the two objects: We add cell IDs to make sure every cell name is unique
# combined_object <- merge(x = object1, y = object2, add.cell.ids = c("4N", "2N"), project = "CombinedAnalysis")
# combined_object[["RNA"]] <- JoinLayers(combined_object[["RNA"]])
# print(combined_object)
# 
# ## scRNAseq derived population- average copy number assigned to each segment
# # ii=names(cn[[numbatRun]]$cells) %in% c(colnames(object1), colnames(object2))
# # segments$CN_Estimate=apply(cn[[numbatRun]]$cn[ii,],2,mean) * 1.25
# 
# ## Karyo derived population- average copy number assigned to each segment
# n2_frac_karyo = length(grep("N2",rownames(whole_chr_karyo)))/nrow(whole_chr_karyo);
# n2_frac_seq = ncol(object2)/(ncol(object2) + ncol(object1))
# apply(whole_chr_karyo,2,mean)
# ## Align segments karyo vs scRNAseq
# la=alignCNmatrices(whole_chr_karyo, cn[[numbatRun]], arm_level_karyo = F)
# segments=as.data.frame(parseLOCUS(colnames(la$cn_scRNAseq)))
# rownames(segments)=paste0(segments$chr,":",segments$startpos,"-",segments$endpos)
# segments$CN_Estimate=apply(la$cn_karyo,2,mean, na.rm=T)
# plot(segments$CN_Estimate)
# lines(segments$CN_Estimate)
# 
# ## Now run Liayson
# epg = as.matrix(combined_object@assays$RNA$counts)
# rownames(epg) = gsub("GRCh38-","",rownames(epg))
# epg=epg[grep("GRCm39", rownames(epg), invert = T),]; ## exclude mouse genes
# eps = aggregateSegmentExpression(epg,as.matrix(segments),host=HOST,mingps = 20,GRCh=38)$eps
# gpc=apply(epg>0,2,sum); 
# names(gpc)=colnames(epg)
# cps=segmentExpression2CopyNumber(eps,gpc,cn=as.matrix(segments)[rownames(eps),"CN_Estimate"],seed = 0.75, nCores = 2, stdOUT="~/Downloads/log.liayson")
# cps=cps[!apply(is.na(cps),1,all),]
# ## Plot heatmap
# origin=sapply(strsplit(colnames(cps),"_"),"[[",1)
# col =RColorBrewer::brewer.pal(2,"Paired")[1:2]
# names(col) = unique(origin)
# hm=heatmap.2((cps),trace = "n", ColSideColors = col[origin], hclustfun = function(x) hclust(x, method = "ward.D2"))
# legend("topright",names(col),fill = col,cex=1.25)
# 
# ## Plot
# ploidy=apply(cps,2,sum)
# vioplot::vioplot(ploidy~origin,log="")



#######################################
#### Merge scRNA-seq & Karyotyping ####
numbat2keep=c("C4N_chr11","C2N_chr23")
cn = cn[numbat2keep]
# cn=cn[-2]
# subset="none"
subset="2N"
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
dfun = chrWeightedCorrDist
# dfun = chrCorrDist
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
seg2exclude=c("13:113164695-114314503_q", "13:19633659-48007418_q")
cn_comb = cn_comb[,!colnames(cn_comb) %in% seg2exclude]
rownames(cn_comb) = paste0(rownames(cn_comb),"_",c(rep("karyo",nrow(cl_k)),rep("sc",nrow(cl_s))))
hm=heatmap.2(round(cn_comb), Colv = NULL, trace = "n", hclustfun=function(x) hclust(x, method="ward.D2"),distfun=dfun, margins = mar) 
cl=cutree(as.hclust(hm$rowDendrogram), k=8)
origin=unique(cl)
col =rainbow(length(origin)*1.1)[1:length(origin)]
names(col) = as.character(origin)
heatmap.2(round(cn_comb), Colv = NULL, trace = "n", RowSideColors =col[as.character(cl)], hclustfun=function(x) hclust(x, method="ward.D2"),distfun=dfun, main=subset, margins = mar) 

dev.off()

##chatGPT: define reusable functions
## conclusions: refine chr arm level calls for 2N
## rerun numbat for 4N.
## annotate cell representation or % on combined heatmap