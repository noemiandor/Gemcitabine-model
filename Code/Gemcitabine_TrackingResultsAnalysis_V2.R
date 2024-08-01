library(mclust)
library(mnormt)
library(matlab)
library( celltrackR )
library( xlsx )
library(ggplot2)


maindir="/Users/4477116/Documents/projects/polyploidization/Gemcitabine_model/"
#maindir="~/Repositories/Gemcitabine-model/"

setwd(paste0(maindir,"Data/E_row/E_row_inter-division/"))
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")


OUTD=paste0(maindir,"Data/matlab")

assignWGDstatus <- function(track, distr){
  track = track[order(track$t),]
  track$WGD=0
  track_orig=track;
  i=1
  WGD=0;
  while(i<=nrow(track)){
    # resetting cell cycle clock to zero (as if cell had just divided):
    track = track[i:nrow(track),,drop=F]
    track$AreFC=track$Size_in_pixels_0/track$Size_in_pixels_0[1]
    track_orig[rownames(track),"WGD"]=WGD
    
    # ## univariate:
    # pdiv=pnorm(track$AreFC,mean = distr$estimate["mean"], sd = distr$estimate["sd"])
    ## multivariate:
    dat=t(sapply(timepoints2include:nrow(track), function(x) track$AreFC[(x-(timepoints2include-1)):x] ))
    # print(dim(dat))
    if(timepoints2include==1){
      dat=t(dat)
      pdiv=pnorm(dat, mean = as.numeric(distr$parameters$mean), sd = sqrt(distr$parameters$variance$sigmasq))
    }else{
      pdiv=pmnorm(dat, mean = as.numeric(distr$parameters$mean), varcov = distr$parameters$variance$Sigma)
    }
    
    if(all(is.na(pdiv)) | all(pdiv<0.5)){
      break
    }
    # Probability of dividing in the subsequent timestep exceeds 0.5:
    i=which(pdiv>=0.5)[1]+1
    WGD=WGD+1;
  }
  plot(track_orig$t, track_orig$Size_in_pixels_0, col=track_orig$WGD+1)
  return(track_orig)
}

plateMap=read.xlsx(paste0(maindir,"Data/Gemcitabine_PlateMap_20240111.xlsx"), sheetIndex = 1)
tp=1:40; ## timepoints of interest to be recorded for matlab model fit
timepoints2include=1; ## for multivariate gaussian fit
train_test<-list(test=c("E6"), train=c("E2"))
#train<-apply(expand.grid(Rows, Cols$train), 1, function(x) paste(x, collapse = ""))
#test<-apply(expand.grid(Rows, Cols$test), 1, function(x) paste(x, collapse = ""))


for(i in 1:length(train_test$train)){
  trainRow=substr(train_test$train[i],1,1)
  testRow=substr(train_test$test[i],1,1)
  f=list(train=list.files(paste0(maindir,"Data/",trainRow,"_row/",trainRow,"_row_inter-division_new/"),pattern=paste0("*",train_test$train[i],".*","\\.txt"), full.names =TRUE))
  f$test=list.files(paste0(maindir,"Data/",testRow,"_row/",testRow,"_row_post-division/"),pattern=paste0("*",train_test$test[i],".*","\\.txt"), full.names = TRUE)
  
  daughterParentCells <- deadCells <- list()
  for(what in names(train_test)){
    print(what)
    tmp=sapply(f[[what]], function(x) read.table(x), simplify = F)
    names(tmp) = sapply(names(tmp), function(x) fileparts(x)$name)
    daughterParentCells[[what]] = tmp
    
    dt_hours=sapply(daughterParentCells[[what]], function(x) quantile(x$t,c(0,1)))
    dt_hours= 2*(dt_hours[2,]-dt_hours[1,])
    hist(dt_hours)
    
    ## Dead cell count
    dead=list()
    for(t in tp){
      x=sapply(daughterParentCells[[what]], function(x) x[x$t==t & x$Classifier.Phenotype=="Dead",], simplify = F)
      x=do.call(rbind,x)
      dead[[as.character(t)]]=nrow(x)
    }
    deadCells[[what]]= dead
    
    print(paste("daughterParent cells total:",length(daughterParentCells[[what]])))
    ##  inter-division tracks which lasted at least 18 hours (doubling time for SUM-159 is 22 hours).
    daughterParentCells[[what]] = daughterParentCells[[what]][dt_hours>=18];
    print(paste("daughterParent cells surviving at least 18 hours:",length(daughterParentCells[[what]])))
    
    ## low dead cell representation:
    x=sapply(daughterParentCells[[what]], function(x) sum(x$Classifier.Phenotype=="Dead")/nrow(x))
    daughterParentCells[[what]]=daughterParentCells[[what]][x<0.1]
    print(paste("daughterParent cells with low fraction dead cells:",length(daughterParentCells[[what]])))
    
    ## Strong correlation between time and cell area (Pearson r>=0.1), suggesting these are indeed cells that progress through the cell cycle.
    pearson_R=sapply(daughterParentCells[[what]], function(x) cor(x$t,x$Size_in_pixels_0))
    hist(pearson_R)
    daughterParentCells[[what]] = daughterParentCells[[what]][pearson_R>=0.3];
    print(paste("daughterParent cells with strong correlation between size and time:",length(daughterParentCells[[what]])))
  }
  
  # multivariate:
  sizeFoldChange=sapply(daughterParentCells$train, function(x) x$Size_in_pixels_0[nrow(x):(nrow(x)-(timepoints2include-1))]/x$Size_in_pixels_0[which.min(x$t)])
  print(paste("timepoints2include",timepoints2include))
  ## univariate:
  if(timepoints2include==1){
    sizeFoldChange=matrix(sizeFoldChange,nrow=1, ncol=length(sizeFoldChange))
  }
  print(paste("sizeFoldChange","nrow", nrow(sizeFoldChange), "ncol", ncol( sizeFoldChange), "length", length( sizeFoldChange), "class",class( sizeFoldChange)))
  sizeFoldChange<-na.omit(sizeFoldChange) #remove unkonwn vaues
  
  if(timepoints2include>1){
    d=mvn("XXX",t(sizeFoldChange)); 
  }else{
    d=mvn("X",t(sizeFoldChange)); 
  }
  
  ## Now apply trained model on test set:
  well=train_test$test[i]
  WGD=sapply(daughterParentCells$test, function(x) try(assignWGDstatus(x, d)), simplify = F)
  WGD=WGD[sapply(WGD, class)!="try-error"]
  
  ####################################################
  ###  Extracting Cells with certain WGD events  #####
  ####################################################
  track_ids <- names(WGD)
  
  #Cells with zero WGD events at ALL time points
  track_0wgd=track_ids[sapply(WGD,function(x) all(x$WGD==0))]
  writeLines(track_0wgd, paste0(maindir,"Data/",train_test$test[i],"_0_wgd_track_IDs.txt"))
  
  wgd_0_df<-do.call(rbind.data.frame, WGD[sapply(WGD,function(x) all(x$WGD==0))])
  write.table(  wgd_0_df, paste0(maindir,"Data/",train_test$test[i],"_","0_wgd_df.txt"),sep = "\t",quote=F,row.names=T,col.names = FALSE)
  
  
  
  # Tracks with  j WGD events at at least one time point
  for(j in c(1,2,3,4,5)){
    temp_ids<-paste0("track_",j,"wgd")
    temp_all_Data<-paste0("wgd",j)
    
    assign(temp_ids,track_ids[sapply(WGD,function(x) any(x$WGD==j))])
    assign(temp_all_Data, do.call(rbind.data.frame, WGD[sapply(WGD,function(x) any(x$WGD==j)),drop=F]))
    
    if(!is.null(as.character(get(temp_ids)))&& length(as.character(get(temp_ids)))>0){
      writeLines(as.character(get(temp_ids)), paste0(maindir,"Data/",train_test$test[i],"_",j,"_wgd_track_IDs.txt"))
      write.table(get(temp_all_Data), paste0(maindir,"Data/",train_test$test[i],"_",j,"_wgd_df.txt"),sep = "\t",quote=F,row.names=T,col.names = T)
    }
  }
  
  
  
  ## WGD distribution per timepoint
  cells=list()
  wgdMax=max(sapply(WGD, function(x) max(x$WGD)))
  for(t in tp){
    x=sapply(WGD, function(x) x[x$t==t,], simplify = F)
    x=do.call(rbind,x)
    x=plyr::count(c(x$WGD,0:wgdMax))
    rownames(x) = x$x
    x=x[as.character(0:wgdMax),-1,drop=F]
    cells[[t]] = x
  }
  cells=do.call(cbind,cells)
  colnames(cells)=paste(as.character(tp*2), well)
  rownames(cells)=paste0(rownames(cells),"_MitosesSkipped")
  cells=rbind(unlist(deadCells$test),cells)
  rownames(cells)[1]="Dead"
  
  ## Save output for Matlab code
  tmp=paste(strsplit(train_test$test[i],"")[[1]],collapse = "_")
  write.table(cells[4:ncol(cells)],paste0(OUTD,filesep,tmp,".txt"),row.names = TRUE,quote = F)
  
  ## Everything below is plotting only: exclude dead cells from plots
  cells=cells[-1,]
  
  ## @TODO: remove first @timepoints2include timepoints from plot since they are meaningless
  png(file=paste0(maindir,"Figs/trained_",trainRow,"2_ID_applied_",well,"_PD_count.png"),
      width=600, height=538)
  barplot(as.matrix(cells)[,4:ncol(cells)], col=rainbow(nrow(cells)), xlab="timepoint", ylab="cell count");
  legend("bottomleft", as.character(0:wgdMax), fill=rainbow(nrow(cells)),title="WGD",bty = "")
  dev.off()
  
  cells_=sweep(cells,MARGIN = 2, STATS=apply(cells,2,sum), FUN = "/")
  png(file=paste0(maindir,"Figs/trained_",trainRow,"2_ID_applied_",well,"_PD_fraction.png"),
      width=600, height=538)
  barplot(as.matrix(cells_)[,4:ncol(cells)], col=rainbow(nrow(cells)), xlab="timepoint", ylab="cell fraction");
  legend("bottomleft", as.character(0:wgdMax), fill=rainbow(nrow(cells)),title="WGD",bty = "")
  dev.off()
}
