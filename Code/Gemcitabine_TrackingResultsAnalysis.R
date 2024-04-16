library(matlab)
library( celltrackR )
library( fitdistrplus )

# setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/B02_20230614_CellTracking_Ilastik")
# setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/C02_20230726_CellTracking_Ilastik")
# setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/D02_20230804_CellTracking_Ilastik/Cell_Tracking_Results(Dead_Cell_Exclusion)_disaprear_cost-100")
# setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/E02_20230817_CellTracking_Ilastik")
# setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/F02_20230825_CellTracking_Ilastik")
# setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/G02_20230831_CellTracking_Ilastik")
setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/H02_20230831_CellTracking_Ilastik")
f=list.files(pattern=".csv")
# i=17;  #2N
i=grep("E2_1_",f);  #4N
dat=read.csv(f[i]);
# well=gsub("CSV-Table.h5.csv","",f[i])
well=gsub("CSV-Table.tiff.csv","",f[i])
ii=match(c("frame","Bounding_Box_Maximum_0","Bounding_Box_Maximum_1"),colnames(dat))
colnames(dat)[ii]=c("t","x","y")
la=sapply(unique(dat$lineageId), function(x) dat[dat$lineageId==x,c("t","x","y","Object_Area_0")], simplify = F)
names(la)=as.character(unique(dat$lineageId))
la=la[!names(la) %in% c("-1")]
plot(as.tracks(la))

## consider only cells that were not yet there at timepoint 0
daughterCells=la[sapply(la, function(x) min(x$t)>0)]
daughterParentCells = daughterCells[sapply(daughterCells, function(x) any(duplicated(x$t)) ) ]
daughterParentCells = sapply(daughterParentCells, function(x) x[!duplicated(x$t),,drop=F], simplify = F)
sizeFoldChange=sapply(daughterParentCells, function(x) x$Object_Area_0[which.max(x$t)]/x$Object_Area_0[which.min(x$t)])
hist(sizeFoldChange[sizeFoldChange<quantile(sizeFoldChange,0.75)],10)
hist(log(sizeFoldChange))
d=fitdist(sizeFoldChange,"norm")


## Integrate with cell classification results
seg=read.csv("../B01_20230407_Incucyte_Images_woGFP_Analysis_QI_Core/Results_Classifier_2/E_objectdata.csv")
seg=seg[grep(well,seg$Image.Location),]
seg$time=sapply(strsplit(seg$Image.Location,"_"), function(x) x[length(x)])
seg$time=gsub(".tif","",seg$time)
## replace time with frame:
map=0:(length(unique(seg$time))-1)
names(map)=sort(unique(seg$time))
seg$t=map[seg$time]
seg=seg[order(seg$t),]
boxplot(seg$t~seg$time, horizontal = T,las=2,ylab="",cex.axis=0.54)
## x, y comparison
ii=match(c("XMax","YMax"),colnames(seg))
colnames(seg)[ii]=c("x","y")
coi=intersect(colnames(seg),colnames(dat))
ii=sapply(map, function(t) list(seg=which(seg$t==t), dat=which(dat$t==t)), simplify = F)
ii=ii[sapply(ii, function(x) min(length(x[[1]]),length(x[[2]])))>0]
d=sapply(ii, function(t) flexclust::dist2(seg[t$seg,coi,drop=F],dat[t$dat,coi,drop=F]))
## assigned each tracked entry the closest matching segmentation entry
seg_matched=sapply(names(d), function(x) seg[ii[[x]]$seg,][apply(d[[x]],2,which.min),],simplify = F)
merged=sapply(names(d), function(x) cbind(dat[ii[[x]]$dat,], seg_matched[[x]][,c("time","Classifier.Phenotype")]), simplify = F)
merged_B=merged
merged=do.call(rbind, merged)
merged=merged[order(merged$t),];

## center of object:
merged$y=(merged$Bounding_Box_Minimum_1+merged$y)/2
merged$x=(merged$Bounding_Box_Minimum_0+merged$x)/2
merged$y=-merged$y
merged$y=merged$y-min(merged$y)


## Alive = 1; transitional = 0; Dead = -1
merged$Class=1
merged$Class[merged$Classifier.Phenotype=="Transitional"]=0
merged$Class[merged$Classifier.Phenotype=="Dead"]=-1


# Is number of tracked cells same as number of segmented cells:
sapply(ii, function(x) c(length(x$dat),length(x$seg)))


## Visualize merged dataset to double-check:
la=sapply(unique(merged$lineageId), function(x) merged[merged$lineageId==x,c("t","x","y","Object_Area_0","time","Classifier.Phenotype","Class")], simplify = F)
names(la)=as.character(unique(merged$lineageId))
la=la[!names(la) %in% c("-1")]
la_=la
# la_=sapply(la, function(x) x[x$t<(24*10)/4,],simplify = F)
# la_=la_[sapply(la_,nrow)>=10 & sapply(la_, function(x) min(x$t))<100]
te=sapply(la_, function(x) cor.test(x$t,x$Class)$estimate)
hist(te, xlab="Pearson correlation(time, cell state)",main="Correlation between time and assigned class")
mtext("Dead cells should not come back to live (corr should be negative for dying cells)")
mtext("corr should be near-zero for cells that never die",line = 3)
# la_=la_[te<0.1];; ## dead cells
la_=la_[is.na(te) | te<0.05]; ## cells that are alive at first
la_=la_[!sapply(la_,is.null)];
print(paste(length(la_)/length(te),"cells included (which do not come back to live from dead)!"))


## look at delta size over time
# la_=sapply(la_, function(x) x[x$t<(24*5)/4,],simplify = F)
tmp=do.call(rbind,la_)
plot(tmp$t,tmp$Object_Area_0,col="white")
sapply(names(la_), function(x) points(la_[[x]]$t,la_[[x]]$Object_Area_0,col=as.numeric(x),pch=20))

# la_=sapply(la_, function(x) x[x$t<(24*5)/4,],simplify = F)
te=sapply(la_, function(x) cor.test(x$t,x$Object_Area_0)$estimate)
hist(te)


## combine two tracks
cids=as.character(c(87,1005))
toi=rbind(la_[[cids[1]]], la_[[cids[2]]])
trackID="combined"


## visualize one track at a time on raw images
fi=list.files("../B01_20230407_Incucyte_Images_woGFP_Analysis_QI_Core",pattern=well,full.names = T)
# trackIDs=names(la_)[order(sapply(la_,nrow),decreasing = T)]
trackIDs=names(la_)[order(te,decreasing = T)]
for(trackID in trackIDs){
  toi=la_[[trackID]]
  pdf(paste0("~/Downloads/testTrack_",trackID,".pdf"))
  
  plot(toi$t,toi$Object_Area_0, col=toi$Class+2,pch=20); 
  legend("topleft",unique(toi$Classifier.Phenotype), fill=unique(toi$Class+2))
  
  fidx=grep(toi$time[1],fi,value=F)
  center=NULL
  # for(i in 1:nrow(toi)){
  while(fidx<=length(fi)){  
    thistime=strsplit(fileparts(fi[fidx])$name,"_")[[1]][4]
    i=which(toi$time==thistime)
    if(is.null(center)){
      center=c(toi$x[i], toi$y[i])
    }
    img=bioimagetools::readTIF(fi[fidx],as.is = T)
    img <- EBImage::resize(img, dim(img)[1]/1)
    plot(raster::as.raster(img[,,,1]), xlim=c(center[1]-200,center[1]+200), ylim=c(center[2]-50,center[2]+50))
    mtext(fileparts(fi[fidx])$name)
    ## Just for testing
    # tmp=merged[merged$t==toi$t[i],]
    # points(tmp$x,tmp$y,pch=3,col="green")
    if(!isempty(i)){
      points(toi$x[i], toi$y[i],pch=3,cex=3)
    }
    fidx=fidx+1;
  }
  dev.off()
}


## plot overlays along with tracking ID:
pdf(paste0("~/Downloads/trackIDs.pdf"))
fidx=1
while(fidx<=length(fi)){  
  img=bioimagetools::readTIF(fi[fidx],as.is = T)
  img <- EBImage::resize(img, dim(img)[1]/1)
  plot(raster::as.raster(img[,,,1]), xlim=c(211,700), ylim=c(150,650))
  mtext(fileparts(fi[fidx])$name)
  
  thistime=strsplit(fileparts(fi[fidx])$name,"_")[[1]][4]
  tmp=sapply(names(la_), function(x) la_[[x]][la[[x]]$time==thistime,,drop=F], simplify = F)
  tmp=sapply(names(tmp), function(x) cbind(tmp[[x]], rep(x, nrow(tmp[[x]]))), simplify = F)
  tmp=do.call(rbind, tmp)
  tmp=as.data.frame(tmp)
  tmp=tmp[which(!sapply(tmp$x,isempty)),]
  colnames(tmp)[ncol(tmp)]="lineageID"
  text(unlist(tmp$x), unlist(tmp$y),labels = unlist(tmp$lineageID), cex=0.6)
  fidx=fidx+1;
}
dev.off()


#Plot all points all times
fi <- list.files("../B01_20230407_Incucyte_Images_woGFP_Analysis_QI_Core",pattern=well,full.names = T)
for(i in unique(merged1$t)){
  p <- merged1 %>% filter(t == i, trackId != -1)
  pdf(paste0("~/Documents/", "tracking_", i, ".pdf"))
  img=bioimagetools::readTIF(fi[i+1],as.is = T)
  img <- EBImage::resize(img, dim(img)[1]/1)
  plot(raster::as.raster(img[,,,1]))
  points(p$x, p$y,pch=(p$Class+2),cex=0.5)
  text(p$x+20, p$y, labels=p$trackId, cex=0.25)
  legend("topleft",unique(merged$Classifier.Phenotype), pch=unique(merged$Class+2), cex = 0.5)
  dev.off()
}



