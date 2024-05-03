library(mclust)
library(mnormt)
library(matlab)
library( celltrackR )
library( xlsx )
setwd("/Users/john/Documents/IMO/polyploidization/Gemcitabine_model/Data/A_row/A_row_inter-division/")
devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")


OUTD="/Users/john/Documents/IMO/polyploidization/Gemcitabine_model/Data/matlab"

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
    
    if(all(!is.na(pdiv)<0.5)){
      break
    }
    # Probability of dividing in the subsequent timestep exceeds 0.5:
    i=which(!is.na(pdiv)>=0.5)[1]+1
    WGD=WGD+1;
  }
  plot(track_orig$t, track_orig$Object_Area_0, col=track_orig$WGD+1)
  return(track_orig)
}

plateMap=read.xlsx("/Users/john/Documents/IMO/polyploidization/Gemcitabine_model/Data/Gemcitabine_PlateMap_20240111.xlsx", sheetIndex = 1)

train_test<-list(train=c("2"),test=c("3","4","5","6","7","8","9","10","11"))

for(train_row in c("A","F")){
    tp=1:40; ## timepoints of interest to be recorded for matlab model fit
    timepoints2include=3 ; ## for multivariate gaussian fit
  
  
    f=list.files(paste0("/Users/john/Documents/IMO/polyploidization/Gemcitabine_model/Data/",train_row,"_row/",train_row,"_row_inter-division/"),
                 pattern=paste0("*",train_row,"2.*","\\.txt"), full.names = F)
   
   setwd(paste0("/Users/john/Documents/IMO/polyploidization/Gemcitabine_model/Data/",train_row,"_row/",train_row,"_row_inter-division/"))
    daughterParentCells=sapply(f, function(x) read.table(x), simplify = F)
    
  
    
    
    dt_hours=sapply(daughterParentCells, function(x) quantile(x$t,c(0,1)))
    dt_hours= 2*(dt_hours[2,]-dt_hours[1,])
    hist(dt_hours)
    
    ## Dead cell count
    dead=list()
    for(t in tp){
      x=sapply(daughterParentCells, function(x) x[x$t==t & x$Classifier.Phenotype=="Dead",], simplify = F)
      x=do.call(rbind,x)
      dead[[as.character(t)]]=nrow(x)
    }
    
    ##  inter-division tracks which lasted at least 18 hours (doubling time for SUM-159 is 22 hours).
    daughterParentCells = daughterParentCells[dt_hours>=18];
    pearson_R=sapply(daughterParentCells, function(x) cor(x$t,x$Object_Area_0))
    hist(pearson_R)
    
    ## a strong correlation between time and cell area (Pearson r>=0.1), suggesting these are indeed cells that progress through the cell cycle.
    daughterParentCells = daughterParentCells[pearson_R>=0.3];
    
    ## low dead cell representation:
    x=sapply(daughterParentCells, function(x) sum(x$Classifier.Phenotype=="Dead")/nrow(x))
    daughterParentCells=daughterParentCells[x<0.1]
    
    ## univariate:
    # sizeFoldChange=sapply(daughterParentCells, function(x) x$Object_Area_0[which.max(x$t)]/x$Object_Area_0[which.min(x$t)])
    # hist((sizeFoldChange))
    # hist(log(sizeFoldChange))
    # d=fitdist(sizeFoldChange,"norm")
    # multivariate:
    sizeFoldChange=sapply(daughterParentCells, function(x) x$Object_Area_0[nrow(x):(nrow(x)-(timepoints2include-1))]/x$Object_Area_0[which.min(x$t)])
    sizeFoldChange<-na.omit(sizeFoldChange) #remove unkonwn vaues
    
    sizeFoldChange<-sizeFoldChange[is.finite(sizeFoldChange)] #remove inf vaues
    #sizeFoldChange=sizeFoldChange[,apply(is.finite(sizeFoldChange),2,all)]
    d=mvn("XXX",t(sizeFoldChange)); 
    
    for(test_row in c("F")){
      for(col in train_test$test){
      well=paste0(test_row,col)
      # @TODO: use post-division tracks (not just inter-division tracks) and read them in here
      setwd(paste0("/Users/john/Documents/IMO/polyploidization/Gemcitabine_model/Data/",test_row,"_row/",test_row,"_row_post-division/"))
      post_divsion=list.files(path = paste0("/Users/john/Documents/IMO/polyploidization/Gemcitabine_model/Data/",test_row,"_row/",test_row,"_row_post-division/"),
                              pattern=paste0("*",well,".*","\\.txt"), full.names = F)
      
          # confirm there is data for this row and column
        daughterParentCells=sapply(post_divsion, function(x)read.table(x),simplify=F)
      
      
      
      
        if(length(daughterParentCells)!=0){
      WGD=sapply(daughterParentCells, function(x) try(assignWGDstatus(x, d)), simplify = F)
      WGD=WGD[sapply(WGD, class)!="try-error"]
      
      ## WGD distribution per timepoint
      cells=list()
      if(length(WGD)!=0){
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
      cells=rbind(unlist(dead),cells)
      rownames(cells)[1]="Dead"
      
      ## Save output for Matlab code
      train_applied=paste0(test_row,"_",well)
      write.table(cells,paste0(OUTD,filesep,train_applied,".txt"),row.names = T,quote = F)
      
      ## plot: @TODO save plot
      barplot(as.matrix(cells), col=rainbow(nrow(cells)), xlab="timepoint", ylab="cell count");
      cells_=sweep(cells,MARGIN = 2, STATS=apply(cells,2,sum), FUN = "/")
      barplot(as.matrix(cells_), col=rainbow(nrow(cells)), xlab="timepoint", ylab="cell fraction");
      legend("bottomleft", as.character(0:wgdMax), fill=rainbow(nrow(cells)),title="WGD",bty = "")
      
      
      
      png(file=paste0("/Users/john/Documents/IMO/polyploidization/inter_results/trained_",train_row,"2_ID_applied_",well,"_PD_count.png"),
          width=600, height=538)
      barplot(as.matrix(cells), col=rainbow(nrow(cells)), xlab="timepoint", ylab="cell count");
      dev.off()
      
      cells_=sweep(cells,MARGIN = 2, STATS=apply(cells,2,sum), FUN = "/")
      png(file=paste0("/Users/john/Documents/IMO/polyploidization/inter_results/trained_",train_row,"2_ID_applied_",well,"_PD_fraction.png"),
          width=600, height=538)
      barplot(as.matrix(cells_), col=rainbow(nrow(cells)), xlab="timepoint", ylab="cell fraction");
      legend("bottomleft", as.character(0:wgdMax), fill=rainbow(nrow(cells)),title="WGD",bty = "")
      dev.off()
        }
       }
      }
   }
}







