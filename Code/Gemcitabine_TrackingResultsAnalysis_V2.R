library(mclust)
library(mnormt)
library(matlab)
library( celltrackR )
library( xlsx )
setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")

assignWGDstatus <- function(track, distr){
  track = track[order(track$t),]
  track$WGD=0
  track_orig=track;
  i=1
  WGD=0;
  while(i<=nrow(track)){
    # resetting cell cycle clock to zero (as if cell had just divided):
    track = track[i:nrow(track),,drop=F]
    track$AreFC=track$Object_Area_0/track$Object_Area_0[1]
    track_orig[rownames(track),"WGD"]=WGD
    
    # ## univariate:
    # pdiv=pnorm(track$AreFC,mean = distr$estimate["mean"], sd = distr$estimate["sd"])
    ## multivariate:
    dat=t(sapply(timepoints2include:nrow(track), function(x) track$AreFC[(x-(timepoints2include-1)):x] ))
    pdiv=pmnorm(dat, mean = as.numeric(distr$parameters$mean), varcov = distr$parameters$variance$Sigma)
      
    if(all(pdiv<0.5)){
      break
    }
    # Probability of dividing in the subsequent timestep exceeds 0.5:
    i=which(pdiv>=0.5)[1]+1
    WGD=WGD+1;
  }
  plot(track_orig$t, track_orig$Object_Area_0, col=track_orig$WGD+1)
  return(track_orig)
}


timepoints2include=3 ; ## for multivariate gaussian fit
plateMap=read.xlsx("Gemcitabine_PlateMap_20240111.xlsx", sheetIndex = 1)

f=list.files("J01_20240111_CellTracking_Ilastik/F_row/F6_1_inter-division/", full.names = T)
# f=list.files("J01_20240111_CellTracking_Ilastik/A_row/A_row_inter-division/", full.names = T)
daughterParentCells=sapply(f, function(x) read.table(x), simplify = F)
dt_hours=sapply(daughterParentCells, function(x) quantile(x$t,c(0,1)))
dt_hours= 2*(dt_hours[2,]-dt_hours[1,])
hist(dt_hours)

##  inter-division tracks which lasted at least 18 hours (doubling time for SUM-159 is 22 hours).
daughterParentCells = daughterParentCells[dt_hours>=18];
pearson_R=sapply(daughterParentCells, function(x) cor(x$t,x$Object_Area_0))
hist(pearson_R)

## a strong correlation between time and cell area (Pearson r>=0.1), suggesting these are indeed cells that progress through the cell cycle.
daughterParentCells = daughterParentCells[pearson_R>=0.3];

## low dead cell representation:
x=sapply(daughterParentCells, function(x) sum(x$Classifier.Phenotype=="Dead")/nrow(x))
daughterParentCells=daughterParentCells[x<0.1]

# univariate:
sizeFoldChange=sapply(daughterParentCells, function(x) x$Object_Area_0[which.max(x$t)]/x$Object_Area_0[which.min(x$t)])
hist((sizeFoldChange))
hist(log(sizeFoldChange))
d=fitdist(sizeFoldChange,"norm")
# multivariate:
sizeFoldChange=sapply(daughterParentCells, function(x) x$Object_Area_0[nrow(x):(nrow(x)-(timepoints2include-1))]/x$Object_Area_0[which.min(x$t)])
sizeFoldChange=sizeFoldChange[,apply(is.finite(sizeFoldChange),2,all)]
d=mvn("XXX",t(sizeFoldChange)); 

# @TODO: use post-division tracks (not just inter-division tracks) and read them in here
WGD=sapply(daughterParentCells, function(x) try(assignWGDstatus(x, d)), simplify = F)
WGD=WGD[sapply(WGD, class)!="try-error"]

## WGD distribution per timepoint
cells=list()
wgdMax=max(sapply(WGD, function(x) max(x$WGD)))
tp=1:40
for(t in tp){
  x=sapply(WGD, function(x) x[x$t==t,], simplify = F)
  x=do.call(rbind,x)
  x=plyr::count(c(x$WGD,0:wgdMax))
  rownames(x) = x$x
  x=x[as.character(0:wgdMax),-1,drop=F]
  cells[[t]] = x
}
cells=do.call(cbind,cells)
colnames(cells)=as.character(tp)

## plot: @TODO save plot
barplot(as.matrix(cells), col=rainbow(nrow(cells)), xlab="timepoint", ylab="cell count");
cells_=sweep(cells,MARGIN = 2, STATS=apply(cells,2,sum), FUN = "/")
barplot(as.matrix(cells_), col=rainbow(nrow(cells)), xlab="timepoint", ylab="cell fraction");
legend("bottomleft", as.character(0:wgdMax), fill=rainbow(nrow(cells)),title="WGD",bty = "")