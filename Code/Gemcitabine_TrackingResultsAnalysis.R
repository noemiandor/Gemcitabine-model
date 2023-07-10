
library( celltrackR )

setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/B02_20230614_CellTracking_Ilastik")
f=list.files(pattern=".csv")
# i=17;  #2N
i=209;  #4N
dat=read.csv(f[i]);
well=gsub("CSV-Table.h5.csv","",f[i])
ii=match(c("frame","Bounding_Box_Maximum_0","Bounding_Box_Maximum_1"),colnames(dat))
colnames(dat)[ii]=c("t","x","y")
la=sapply(unique(dat$trackId), function(x) dat[dat$trackId==x,c("t","x","y","Object_Area_0")], simplify = F)
names(la)=as.character(unique(dat$trackId))
la=la[!names(la) %in% c("-1")]
plot(as.tracks(la))


## Integrate with cell classification results
seg=read.csv("../B01_20230407_Incucyte_Images_woGFP_Analysis_QI_Core/Results_Classifier_1/E_objectdata.csv")
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
## assigned each tarcked entry the closest matchig segmentation entry
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
la=sapply(unique(merged$trackId), function(x) merged[merged$trackId==x,c("t","x","y","Object_Area_0","time","Classifier.Phenotype","Class")], simplify = F)
names(la)=as.character(unique(merged$trackId))
la=la[!names(la) %in% c("-1")]
la_=la[sapply(la,nrow)>=50 & sapply(la, function(x) min(x$t))<10]
la_=sapply(la_, function(x) x[x$t<(24*5)/4,],simplify = F)
te=sapply(la_, function(x) cor.test(x$t,x$Class)$estimate)
hist(te, xlab="Pearson correlation(time, cell state)",main="Correlation between time and assigned class")
mtext("Dead cells should not come back to live (corr should be negative for dying cells)")
mtext("corr should be near-zero for cells that never dye",line = 3)
# la_=la_[te<0.1];; ## dead cells
la_=la_[te>-0.1]; ## live cells

## look at delta size over time
# la_=sapply(la_, function(x) x[x$t<(24*5)/4,],simplify = F)
tmp=do.call(rbind,la_)
plot(tmp$t,tmp$Object_Area_0,col="white")
sapply(names(la_), function(x) points(la_[[x]]$t,la_[[x]]$Object_Area_0,col=as.numeric(x),pch=20))

# la_=sapply(la_, function(x) x[x$t<(24*5)/4,],simplify = F)
te=sapply(la_, function(x) cor.test(x$t,x$Object_Area_0)$estimate)
hist(te)


## visualize one track at a time on raw images
fi=list.files("../B01_20230407_Incucyte_Images_woGFP_Analysis_QI_Core",pattern=well,full.names = T)
trackID="114"
for(trackID in names(la_)){
  toi=la_[[trackID]]
  pdf(paste0("~/Downloads/testTrack_",trackID,".pdf"))
  for(i in 1:nrow(toi)){
    t = toi$time[i]
    img=bioimagetools::readTIF(grep(t,fi,value=T),as.is = T)
    img <- EBImage::resize(img, dim(img)[1]/1)
    plot(raster::as.raster(img[,,,1]))
    mtext(paste(well,t))
    ## Just for testing
    # tmp=merged[merged$t==toi$t[i],]
    # points(tmp$x,tmp$y,pch=3,col="green")
    points(toi$x[i], toi$y[i],pch=3,cex=1)
  }
  dev.off()
}
