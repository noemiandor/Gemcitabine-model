library(e1071)
library(dyno)
library(matlab)
setwd("/Users/4470246/Downloads/Post-division")

classifyCellCyclePhase <- function(x, y, main="", svmfit=NULL){
  dat = data.frame(as.matrix(x), y = as.factor(y))
  if(is.null(svmfit)){
    svmfit = svm(y ~ ., data = dat, kernel = "radial", cost = 10, scale = F)
    # print(svmfit)
    try(plot(svmfit, dat,area_nucleus.p~area_mito.p, main=main),silent = T)
  }
  out=predict(svmfit, dat)
  confMat=caret::confusionMatrix(dat$y, out)$byClass
  # print(confMat)
  return(list(svmfit=svmfit, out=out, confusionMatrix=confMat))
}

asDataset<-function(imgStats, imgStats_raw, FoF=NULL, coi){
  ii=1:nrow(imgStats)
  if(!is.null(FoF)){
    ii=which(imgStats$FoF==FoF)
  }
  tmp=apply(imgStats[ii,coi],2,as.numeric)
  tmp = tmp + min(tmp[tmp>0])*0.1
  tmp2=apply(imgStats_raw[ii,coi],2,as.numeric)
  rownames(tmp)<- rownames(tmp2) <- rownames(imgStats)[ii]
  dataset <- wrap_expression(
    expression = log2(tmp),
    counts = tmp2
  )
  return(dataset)
}



f=list.files(pattern = "E")
# f=sample(f,900)
dm=sapply(f, function(x) read.table(x), simplify = F)
dm = dm[sapply(dm, function(x) any(x$Classifier.Phenotype!="Dead"))]
for(x in names(dm)){
  dm[[x]] = dm[[x]][dm[[x]]$Classifier.Phenotype!="Dead",,drop=F] 
  dm[[x]]$cellCycle = "G2/M"
  dm[[x]]$cellCycle[1:round(0.5*nrow(dm[[x]]))] = "G1/S"
  dm[[x]]$lifetime_frac=dm[[x]]$lifetime/max(dm[[x]]$lifetime)
  rownames(dm[[x]])=paste0(x,".",rownames(dm[[x]]))
}
dm = dm[sapply(dm,nrow)>=10]

imgStats_raw=do.call(rbind,dm)
rownames(imgStats_raw) = do.call(c,sapply(dm, rownames))
coi=setdiff(colnames(imgStats_raw),c("lifetime_frac","lifetime","cellCycle","Classifier.Phenotype","Class","t","labelimageId", "trackId", "lineageId", "parentTrackId", "mergerLabelId", "x","y", "parent", "dividing","well_info","time")); 
coi=grep("Center", coi,value=T, invert = T)
coi = grep("Bounding", coi,value=T, invert = T)
coi=coi[apply(imgStats_raw[,coi], 2, function(x) !all(x==0 | is.na(x)))]

imgStats = imgStats_raw
tmp= abs(as.numeric(as.matrix(imgStats_raw[,coi])))
imgStats[,coi] = 0.5*min(tmp[tmp>0]) + sweep(imgStats_raw[,coi], 2, STATS = apply(imgStats_raw[,coi],2,min, na.rm=T),FUN = "-")
imgStats[,coi]=sweep(imgStats[,coi], 2, STATS = apply(imgStats[,coi],2,median, na.rm=T),FUN = "/")
imgStats_=asDataset(imgStats, imgStats_raw, coi=coi)


## Feature selection
svmFeatures=coi
NREP=5;
acc=-Inf
e2cells = grep("E2",rownames(imgStats_raw), value = T);
while(!isempty(svmFeatures)){
  svmFeatures_=combn(svmFeatures,length(svmFeatures)-1)
  ccPred=rep(0,ncol(svmFeatures_))
  print(paste0("training SVM with ", length(svmFeatures)," features"))
  for(i in 1:ncol(svmFeatures_)){
    for(rep in 1:NREP){
      trainCells=sample(e2cells, round(length(e2cells)*0.5) )
      testCells=setdiff(e2cells, trainCells)
      x=imgStats_$expression[trainCells,svmFeatures_[,i]]
      y=imgStats_raw[trainCells,"cellCycle"]
      
      # train SVM to classify cells to cell cycle phase based on allen Model features:
      svm=classifyCellCyclePhase(x, y, main="", svmfit=NULL)
      x=imgStats_$expression[testCells,svmFeatures_[,i]]
      y=imgStats_raw[testCells,"cellCycle"]
      ccPred[i]=ccPred[i]+classifyCellCyclePhase(x, y, main="", svmfit=svm$svmfit)$confusionMatrix["Balanced Accuracy"]
    }
    ccPred[i]=ccPred[i]/NREP
  }
  if(acc - max(ccPred) < 0.02 ){
    acc=max(ccPred)
    print(paste("excluding",setdiff(svmFeatures,svmFeatures_[,which.max(ccPred)]),"; accuracy=",acc))
    svmFeatures = svmFeatures_[,which.max(ccPred)]
  }else{
    break
  }
}
# # ii_B=order(ccPred)
# # ii=order(ccPred)
# plot(ii_B,ii)
# cor.test(ii_B,ii)
# svmFeatures=c(svmFeatures,"Skewness_of_Intensity_1","Terminal_1_1","Maximum_intensity_0","Object_Area_0","Variance_of_Intensity_0","Mean_Defect_Displacement_0","Minimum_intensity_1")

##Training:
cells = sapply(c("E2","E4","E6"), function(x) grep(x,rownames(imgStats_$expression), value = T), simplify = F )
ii=sample(cells$E2, min(10000, length(cells$E2)) )
svm=classifyCellCyclePhase(imgStats_$expression[ii,svmFeatures], imgStats_raw[ii,"cellCycle"], main="", svmfit=NULL)

## Application of trained SVM
output=sapply(cells, function(ii) classifyCellCyclePhase(imgStats_$expression[ii,svmFeatures], imgStats_raw[ii,"cellCycle"], main="", svmfit=svm$svmfit), simplify = F)
out=do.call(c,sapply(output, function(x) x$out) )
names(out) = unlist(sapply(output, function(x) names(x$out)) )
imgStats$cellCycleSVM = as.character(out[rownames(imgStats)])
for(x in names(dm)){
  dm[[x]]$cellCycleSVM = as.character(out[rownames(dm[[x]])])
}


## Pseudotime inference
imgStats__=asDataset(imgStats, imgStats_raw, coi=svmFeatures)
guidelines <- guidelines(
  imgStats__,
  answers = answer_questions(
    imgStats__,
    multiple_disconnected = FALSE,
    expect_topology = TRUE,
    expected_topology = "cycle"
  )
)
model_ <- infer_trajectory(imgStats__, ti_angle(dimred = "pca"))
plot_dimred(model_,color_cells = "pseudotime")


## Compare pseudotime to real time within each track
imgStats$pseudotime =model_$pseudotime[rownames(imgStats)]
plot(imgStats$lifetime_frac,imgStats$pseudotime)
cor.test(imgStats$lifetime_frac,imgStats$pseudotime)
te=sapply(dm[sapply(dm,nrow)>=10], function(x) try(cor.test(x$lifetime, model_$pseudotime[rownames(x)]) ), simplify = F)
te=te[sapply(te, function(x) class(x)!="try-error")]
te=as.data.frame(t(sapply(te, function(x) c(x$estimate,x$p.value))))
colnames(te)=c("r","P")
te$P.adjust=p.adjust(te$P,method = "fdr")
te=te[order(te$P),]
hist(te$r[te$P.adjust<0.1],30)


## Only keep tracks where real time is well correlated to pseudotime
goodTracks = rownames(te)[te$r> 0.4]
badTracks = rownames(te)[te$r<= 0.05]
if(mean(te$r[te$P.adjust<0.1])<0){ ## pseudotime directionality is unknown
  goodTracks = rownames(te)[te$r< -0.4]
  badTracks = rownames(te)[te$r>= 0.05]
}
write.table(goodTracks, "~/Downloads/goodTracks.txt", row.names=F, col.names=F, quote=F)
write.table(badTracks, "~/Downloads/badTracks.txt", row.names=F, col.names=F, quote=F)

## Plot cell cycle composition per timepoint
par(mfrow=c(3,2))
for (well in names(cells)){
  # ii=unlist(sapply(goodTracks, function(x) grep(x, cells[[well]], value=T)))
  ii = cells[[well]]
  imgStats__=imgStats[ii,]
  imgStats__$hour=imgStats__$t*2
  tmp=sapply(unique(imgStats__$hour), function(x) c(unique(imgStats__$cellCycle),imgStats__[imgStats__$hour==x,"cellCycleSVM"]))
  tmp=sapply(tmp, function(x) plyr::count(x)$freq, simplify = T)
  colnames(tmp) = unique(imgStats__$hour)
  tmp=tmp[,order(as.numeric(colnames(tmp)))]
  barplot(tmp, col=rainbow(nrow(tmp)), xlab="timepoint", ylab="cell count", main=well);
  legend("topleft", unique(imgStats__$cellCycle), fill=rainbow(nrow(tmp)),title="Cell Cycle SVM",bty = "")
  tmp = sweep(tmp,MARGIN=2,STATS=apply(tmp,2,sum), FUN="/")
  barplot(tmp, col=rainbow(nrow(tmp)), xlab="timepoint", ylab="cell count", main=well);
}

