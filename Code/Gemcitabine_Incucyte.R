library(xlsx)
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
setwd("/Users/4470246/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/20230407_Incucyte_Images_woGFP_Analysis_QI_Core/Results_Classifier_2/")
# setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023/ObjectData")
OUTDIR="../../../K01_SkippedMitosisClassification_042523"
dir.create(OUTDIR)

toHours<-function(x){
  hours=as.numeric(substr(x,1,2))*24 + as.numeric(substr(x,4,5))
  return(hours)
}
foi=c("Region.Area..μm..","Region.Perimeter..μm.")
maxMitosesSkipped=2
COORD=c("XMin", "XMax", "YMin", "YMax")
calculateNeighborDistance<-function(x, radius = 40){
  CLOSESTLIVECELL=paste0("ClosestLiveCell_",radius)
  CLOSESTNONLIVECELL=paste0("ClosestNonLiveCell_",radius)
  x[,CLOSESTLIVECELL] <- x[,CLOSESTNONLIVECELL] <- 0
  ii=x$Classifier.Phenotype=="Alive"
  if(!any(ii)){
    return(x)
  }
  d_alive=flexclust::dist2(x[ii,COORD,drop=F],x[ii,COORD,drop=F])
  d_alive[d_alive==0]=Inf
  localCells=apply(d_alive, 1, function(y) which(ii)[which(y<radius)], simplify=F)
  names(localCells)=as.character(which(ii))
  localCells=localCells[sapply(localCells, length)>0]
  x[as.numeric(names(localCells)),CLOSESTLIVECELL]=sapply(localCells,function(i) sum(x$Region.Area..μm..[i]))

  if(all(ii)){
    return(x)
  }
  # print(head(x[!ii,COORD,drop=F]))
  d_notalive=flexclust::dist2(x[ii,COORD,drop=F],x[!ii,COORD,drop=F])
  d_notalive[d_notalive==0]=Inf
  localCells=apply(d_notalive, 1, function(y) which(!ii)[which(y<radius)], simplify=F)
  names(localCells)=as.character(which(ii))
  localCells=localCells[sapply(localCells, length)>0]
  x[as.numeric(names(localCells)),CLOSESTNONLIVECELL]=sapply(localCells,function(i) sum(x$Region.Area..μm..[i]))
  # x[ii,CLOSESTNONLIVECELL] = apply(d_notalive, 1,function(x) sum(x<radius, na.rm=T))
  return(x)
}

classifyCells2SkippedMitoCount<-function(t0_2N_4N,maxMitosesSkipped, freq=T){    
  ## @TODO: cont here! <-- rewrite classifier to look at cells within same regions of same frame! consider cell tracking -- ask Mahmoud!
  load(file="~/Downloads/binnedClassifier.RObj")
  t0_2N_4N$bin=round(t0_2N_4N$ClosestLiveCell/1000)
  t0_2N_4N$mitoSkipped=NA
  allcounts=list()
  for(i in which(!is.na(binnedClassifier$N4_N2_ratio))){
    bin =binnedClassifier$Bin[i]
    ii=which(t0_2N_4N$bin==bin)
    baseline=binnedClassifier$N2[i]
    
    classes=0:maxMitosesSkipped*(binnedClassifier$N4_N2_ratio[i]-1)
    # classes=0:maxMitosesSkipped*(1.2-1)
    
    foldChange=t0_2N_4N$Region.Area..μm..[ii]/baseline  -1
    d=flexclust::dist2(classes,foldChange)
    t0_2N_4N$mitoSkipped[ii]=apply(d,2,which.min)-1
    counts=sapply(unique(t0_2N_4N$well), function(well) t0_2N_4N$mitoSkipped[intersect(ii,which(t0_2N_4N$well==well))], simplify = F)
    counts=sapply(counts, function(x) plyr::count(c(0:maxMitosesSkipped,x))$freq, simplify = F)
    if(freq){
      counts=sapply(counts, function(x) x/sum(x), simplify = F)
    }
    allcounts[[as.character(bin)]]=do.call(c,counts)
  }
  allcounts=allcounts[order(as.numeric(names(allcounts)))]
  allcounts=as.matrix(do.call(cbind,allcounts))
  return(allcounts)
}


##Read plate map
plateMap=read.xlsx("../../Gemcitabine_PlateMap_20230327.xlsx", sheetIndex = 1,check.names=F)
rownames(plateMap)=plateMap$Row
plateMap=plateMap[,-1]
ploidy=plateMap[,"Ploidy",drop=F]
plateMap=plateMap[,-1]

##Read data
f=list.files(pattern = "ata.csv")
all=list()
for(f_ in f){
  la=read.csv(f_)
  la$well =sapply(strsplit(la$Image.Location,"_"),"[[",6)
  la$time =sapply(strsplit(la$Image.Location,"_"),"[[",8)
  la$time=gsub(".tif","",la$time)
  la$time=toHours(la$time)
  la_orig=la;
  
  timewell=unlist(sapply(unique(la_orig$time), function(x) paste(x,unique(la_orig$well)), simplify=F))
  all[timewell]=sapply(timewell, function(x)  la[la$time==as.numeric(strsplit(x," ")[[1]][1]) & la$well==strsplit(x," ")[[1]][2],], simplify=F)
}
save(file="~/Downloads/gemcitabine_class1.RObj","all")
# save(file="~/Downloads/gemcitabine_class2.RObj","all")

## Choose between classifiers based on frequency of suspiciously small cells
##load("~/Downloads/gemcitabine_class1.RObj")
alt<-tooSmall<-list()
for (obj in list.files("~/Downloads/",pattern="gemcitabine_class",full.names = T)){
  load(obj)
  # # ##keep only center of image
  # all=sapply(all, function(x) x[x$XMin>200 & x$XMax<1200 & x$YMin>200 & x$YMax<800,], simplify = F)
  
  ## exclude timepoints after 5 days
  keep=which(sapply(all, function(x) x$time[1]<=24*5))
  all=all[keep]
  
  for(well in names(all)){
    x=all[[well]]
    x$radius=2*sqrt(x$Region.Area..μm../pi)
    all[[well]]=x
  }
  alt[[obj]]=all
  tooSmall[[obj]]=sapply(all, function(x) sum(x$Classifier.Phenotype=="Alive" & x$radius<10)/sum(x$Classifier.Phenotype=="Alive"))
}
tooSmall[[2]]=tooSmall[[2]][names(tooSmall[[1]])]
multiNucI=names(which(tooSmall[[1]]>tooSmall[[2]]))
plot(tooSmall[[1]],tooSmall[[2]],xlab=names(tooSmall)[1],ylab=names(tooSmall)[2])
## Generally use classifier 1 except when it predicts a lot of small cells
la=alt$`/Users/4470246/Downloads//gemcitabine_class1.RObj`
la[multiNucI]=alt$`/Users/4470246/Downloads//gemcitabine_class2.RObj`[multiNucI]



## Relabel suspiciously small cells
for(well in names(la)){
  x=la[[well]]
  # col=c("red","purple","cyan")
  # names(col)=unique(x$Classifier.Phenotype)
  # plot(x$XMax,-x$YMin, col=col[x$Classifier.Phenotype])
  
  x$Classifier.Phenotype[x$Classifier.Phenotype== "Alive" & x$radius<15]="Oversegmentation"
  la[[well]]=x
}
plyr::count(la$`40 H1`$Classifier.Phenotype)
sort(la$`40 H1`$radius[la$`40 H3`$Classifier.Phenotype=="Alive"])
sapply(la, function(x) plyr::count(x$Classifier.Phenotype), simplify = F)




## separate live, label time and well
la=la[sapply(la,nrow)>0]
## one replicate only: @TODO for testing purposes only <-- remove once pipeline works
# ii=which(wellCol %in% 1:6) # & wellRow %in% c("H","A") 
# la=la[ii]
time=sapply(la, function(x) x[1,"time"])
wells = sapply(la, function(x) x[1,"well"])
wellRow=sapply(wells, substr,1,1)
wellCol=as.numeric(gsub("_","",sapply(wells,function(x) substr(x,2,nchar(x)))))
la=sapply(la, calculateNeighborDistance, radius=100,simplify=F)
la=sapply(la, calculateNeighborDistance, radius=60,simplify=F)
la=sapply(la, calculateNeighborDistance, radius=40,simplify=F)
la=sapply(la, calculateNeighborDistance, radius=20,simplify=F)
live=sapply(la, function(x) x[x[,"Classifier.Phenotype"] %in% c("Alive"),,drop=F], simplify = F)
dead=sapply(la, function(x) x[x[,"Classifier.Phenotype"] %in% c("Dead","Transitional"),,drop=F], simplify = F)
# ##keep only cells detected as live with high confidence:
# live=sapply(live, function(x) x[x$Alive>0.999,], simplify = F)
names(wellCol) = names(wells)=names(time)=names(live)

# ## exclude cells from dense regions
# live=sapply(live, function(x) x[x$ClosestLiveCell<10,,drop=F], simplify = F)

live=live[sapply(live,nrow)>0]
dead=dead[sapply(dead,nrow)>0]

## add total area occupied to live cell count:
for(x in names(live)){
  live[[x]]$occupiedArea= rep(sum(live[[x]][,"Region.Area..μm.."]), nrow(live[[x]]));
  live[[x]]$time= rep(time[x], nrow(live[[x]]));
  live[[x]]$well= rep(wells[x], nrow(live[[x]]));
  live[[x]]$well= rep(wells[x], nrow(live[[x]]));
  live[[x]]$ploidy= rep(ploidy[wellRow[x],], nrow(live[[x]]));
}
save(file="~/Downloads/gemcitabine_CellsAfterFilter_Class1.RObj", list =c("live","dead","la","plateMap","time","wells","wellRow","wellCol","ploidy") )


# load(file="~/Downloads/gemcitabine_CellsAfterFilter_Class1.RObj")
# load(file="~/Downloads/gemcitabine_CellsAfterFilter_Class2.RObj")
## exclude first timepoints when cells haven't fully attached?
live=live[time[names(live)]>4]
# live=live[time[names(live)]>60]
live=live[time[names(live)]<5*24]


## compare 2N and 4N cells at first timepoint: area
# woi=grep("1$", grep("SUM159",colnames(plateMap)), value = T,invert = F)
woi=grep("*", grep("SUM159",colnames(plateMap)), value = T)
woi_2N=unlist(sapply(rownames(ploidy)[ploidy$Ploidy=="2N"], function(x) paste0(x,woi),simplify = F))
woi_2N=grep("A", woi_2N, value=T)
# woi=grep("1$", grep("SUM159",colnames(plateMap)), value = T,invert = F)
woi=grep("1$", grep("SUM159",colnames(plateMap)), value = T,invert = T)
woi_4N=unlist(sapply(rownames(ploidy)[ploidy$Ploidy=="4N"], function(x) paste0(x,woi),simplify = F))
woi_4N=grep("H", woi_4N, value=T)
## include first 16 hours (before drug effect kicks in) for 4N cells only
ii=which((wells %in% woi_2N & time[names(wells)]<=16) | (wells %in% woi_4N & time[names(wells)]<=16) )
# ii=which((wells %in% woi_2N ) | (wells %in% woi_4N ) )
woi=intersect(names(wells[ii]), names(live))
## Plot by time: each well separately
t0_2N_4N=sapply(woi, function(x) live[[x]][,c("Region.Area..μm..","ClosestLiveCell_60","time","well")], simplify = F)
t0_2N_4N=sapply(unique(wells[woi]), function(well) t0_2N_4N[grep(paste0(well,"$"), names(t0_2N_4N))], simplify = F)
t0_2N_4N=sapply(t0_2N_4N, function(x) do.call(rbind,x), simplify = F)
t0_2N_4N=do.call(rbind,t0_2N_4N)
t0_2N_4N$bin=round(t0_2N_4N$ClosestLiveCell_60/1000)
boxplot(t0_2N_4N$Region.Area..μm..~t0_2N_4N$well+t0_2N_4N$bin,col=c("blue","red"),log="x",horizontal = T, las=2,cex=0.5)
legend("topright",c("2N","4N"),fill=c("blue","red"))
# ## @TODO: ggplot all three (n2, n4, ratio)
# ii=sapply(unique(t0_2N_4N$well), function(well)  sapply(unique(t0_2N_4N$bin), function(bin) which(t0_2N_4N$bin==bin & t0_2N_4N$well==well)))
# n2=sapply(ii[,"A1"], function(i) median(t0_2N_4N$Region.Area..μm..[i]))
# n4=sapply(ii[,"H1"], function(i) median(t0_2N_4N$Region.Area..μm..[i]))
# bin=sapply(ii[,"H1"], function(i) median(t0_2N_4N$bin[i]))
# plot(bin,n4,pch=20,col="red")
# points(bin,n4/n2)
# points(bin,n2,pch=20,col="blue")
# ## @TODO: save n2,n4 per each bin
# binnedClassifier=as.data.frame(cbind(bin,n2, n4,n4/n2))
# colnames(binnedClassifier)=c("Bin","N2","N4","N4_N2_ratio")
# save(file="~/Downloads/binnedClassifier.RObj","binnedClassifier")
# 
# 
# ## classification
# allcounts=classifyCells2SkippedMitoCount(t0_2N_4N, maxMitosesSkipped)
# par(mfrow=c(2,1))
# barplot(allcounts[1:(maxMitosesSkipped+1),])
# barplot(allcounts[(maxMitosesSkipped+2):(maxMitosesSkipped*2+2),])




## predict skipped mitosis from area occupied + nucleus area
t0_2N_4N=sapply(woi, function(x) live[[x]][,c("Region.Area..μm..","Region.Perimeter..μm.","occupiedArea","time",grep("Closest",colnames(live[[x]]),value=T),"well")], simplify = F)
t0_2N_4N=do.call(rbind,t0_2N_4N)

t0_2N_4N$mitosesSkipped=t0_2N_4N$well %in% woi_4N
# te=sapply(foi, function(x) t.test(t0_2N_4N[,x]~t0_2N_4N$mitosesSkipped)$estimate, simplify = T)
## Artificially generate additional ploidy groups (more rounds of ployploidization)
# ratio=te[2,]/te[1,]
# for(level in 1:2){
#   newData=t0_2N_4N[t0_2N_4N$mitosesSkipped==level,]
#   newData[,foi]=sweep(newData[,foi],MARGIN = 2,STATS = ratio,FUN = "*")
#   newData$mitosesSkipped=newData$mitosesSkipped+1
#   t0_2N_4N=rbind(t0_2N_4N,newData)
# }

## @TODO: continue here <-- open up all relevant timepoints for A2 (2N) and H2 (4N) and look at how classifier works for cells that are "alone" (sole survivors in local neighborhood)
## @TODO: continue here <-- problem is there are two separate reasons why local cell density is low, each with different implications for area~ploidy relationship: 
## i) early timepoints, cells didn't have time yet to grow, in which case larger inter-cell distance correctly counteracts large nucleus area
## ii) treatment killing most cells, in which case larger inter-cell distance incorrectly overshadows large nucleus area 
## @TODO: continue here <-- solution: need to include 1st round of treatment-induced polyploidization into training, for low Gemcitabine concentrations 
# <-- look at relation between time after Gemcitabine and area among 4N cells
# <-- find regions of similar density among treatment-naive 2N and 4N cells --> see Notability: Emails: todo next gemcitabine
t0_2N_4N__=t0_2N_4N
print(plyr::count(t0_2N_4N$mitosesSkipped))
# t0_2N_4N__=t0_2N_4N__[t0_2N_4N__$time<60,]
# t0_2N_4N__$Region.Area..μm..=0.9*t0_2N_4N__$Region.Area..μm..; ## Cells under treatment will have very low division rates, and thus will be smaller
m=lm(mitosesSkipped ~ Region.Area..μm..*ClosestLiveCell_100*ClosestNonLiveCell_40*ClosestNonLiveCell_100*ClosestNonLiveCell_60*time, data = t0_2N_4N__)
## Including occupiedArea works much better when treatment is off, but fails when treatment is on:
# m=lm(mitosesSkipped ~ Region.Area..μm..*occupiedArea*ClosestLiveCell_100*ClosestNonLiveCell_40*ClosestNonLiveCell_100*ClosestNonLiveCell_60*time, data = t0_2N_4N__)
print(summary(m))
out=predict(m,data=t0_2N_4N)
# plot(out,col=t0_2N_4N$mitosesSkipped+1)
# plot(out,t0_2N_4N$mitosesSkipped,col=t0_2N_4N$mitosesSkipped+1)
# boxplot(t0_2N_4N$Region.Area..μm..~t0_2N_4N$mitosesSkipped,log="y")
# boxplot(out~t0_2N_4N$mitosesSkipped,log="")

## use model on all data
for (x in names(live)){
  live[[x]]$mitosesSkipped=predict(m,newdata=live[[x]])
  ## correct classification <-> @TODO: should not be necessary
  live[[x]][live[[x]]$ploidy[1]=="4N" & live[[x]]$mitosesSkipped<1,"mitosesSkipped"]=1 
  ii=which(live[[x]]$ploidy=="2N")
  live[[x]][ii,"mitosesSkipped"]=live[[x]][ii,"mitosesSkipped"]-1 
  ## set min and max
  live[[x]][live[[x]]$mitosesSkipped<0,"mitosesSkipped"]=0
  ii=which(live[[x]]$mitosesSkipped>min(live[[x]]$mitosesSkipped))
  live[[x]]$mitosesSkipped[ii]=log(live[[x]]$mitosesSkipped[ii]+1)
  live[[x]]$mitosesSkipped[live[[x]]$mitosesSkipped>maxMitosesSkipped]=maxMitosesSkipped
}



## violin plots in each well across time
live_=sapply(unique(wells), function(x) live[grep(paste0(x,"$"),names(live))], simplify = F)
names(live_)=gsub("_","",names(live_))
pdf("~/Downloads/Gemcitabine.pdf")
par(mfrow=c(1,2))
sapply(names(live_), function(x) try(boxplot(sapply(live_[[x]], function(y) y$Region.Area..μm..), main=x, las=2,horizontal=T)))
dev.off()

## histogram/scatterplots across timepoints per well
dead_=sapply(unique(wells), function(x) dead[grep(paste0(x,"$"),names(dead))], simplify = F)
names(dead_)=gsub("_","",names(dead_))
area=sapply(names(live_), function(x) sapply(live_[[x]], function(y) y[,c(foi,"mitosesSkipped")], simplify = F), simplify = F)
pdf("~/Downloads/Gemcitabine_hist_perTp.pdf", width = 12,height = 5)
par(mfrow=c(2,4))
# ii=order(gsub("11","911",gsub("10","910",gsub("12","912",names(area)))))
ii=order(as.numeric(sapply(names(area), function(x) substr(x,2,nchar(x)))))
for(well in names(area)[ii]){
  ## Dead cells:
  cnt=sapply(dead_[[well]], nrow)
  cnt=cnt[order(time[names(cnt)])]
  wellRow_=wellRow[names(cnt)[1]]
  wellCol_=wellCol[names(cnt)[1]]
  plot(time[names(cnt)],cnt,ylab="dead cells", xlab="hour",pch=20,main=paste(well,ploidy[wellRow_,1],"Gem:",plateMap[wellRow_,wellCol_]))
  
  ## Live cells:
  area_=area[[well]]
  # sapply(names(area_), function(x) hist(area_[[x]][,"Region.Area..μm.."],breaks = 50, xlim=c(0,4000),main=paste(well,time[[x]])))
  cellcomp=sapply(names(area_), function(x) plyr::count(c(0:maxMitosesSkipped,round(area_[[x]]$mitosesSkipped) ) )$freq, simplify = T)
  cellcomp=cellcomp[,order(time[colnames(cellcomp)])]
  rownames(cellcomp)=paste0(0:(nrow(cellcomp)-1),"_MitosesSkipped")
  ## Save data to be read by matlab
  livePlusDead=rbind(cnt[colnames(cellcomp)],cellcomp)
  rownames(livePlusDead)[1]="Dead"
  ##@TODO: don't exclude first 16 hours if we manage to get rid of blue artifact background for early timepoints
  livePlusDead=livePlusDead[,time[colnames(livePlusDead)]>28]
  write.table(livePlusDead, file=paste0(OUTDIR,filesep,well,".txt"),sep="\t",quote = F)
  ## Normalize to show frequency instead of absolute cell counts
  # cellcomp=sweep(cellcomp,MARGIN = 2,FUN = "/",STATS = apply(cellcomp,2,sum))
  ## Multiple nymber of skipped mitosis by two for 4N cells: @TODO -- should not be necessary
  # if(ploidy[wellRow,1]=="4N"){
  #   cellcomp[2:nrow(cellcomp),]=cellcomp[1:(nrow(cellcomp)-1),]
  #   cellcomp[1,]=0
  # }
  ##Plot live cells
  barplot(cellcomp,beside = F, xlab="Number of cells",horiz = T,las=2,names=paste(wells[colnames(cellcomp)], time[colnames(cellcomp)]),col = rainbow(12),main=paste(well,ploidy[wellRow_,1],"Gem:",plateMap[wellRow_,wellCol_]))
}
dev.off()
plot(1)
legend("topleft",paste("Skipped mitoses:",0:maxMitosesSkipped),fill=rainbow(12)[1:(1+maxMitosesSkipped)])


