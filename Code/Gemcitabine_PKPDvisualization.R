setwd("~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/M00_GemcitabinePKPD_101823")

devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")

## Gemcitabine PKPD
dm=list()
coi=c("Sample" ,"Gemcitabine..ng.mL.", "dFdU...ng.mL.","dFdCTP...ng.mL." , "cellType" , "time" )
for(sheet in c("2N","4N","2N_lowInitialGemcitabine","4N_lowInitialGemcitabine")){
  la=read.xlsx("GemcitabineExposure_PKPD.xlsx",sheetName = sheet)
  tmp=strsplit(la$Sample,"-")
  la$cellType=sapply(tmp, "[[",1)
  la$time=as.numeric(gsub("h","",sapply(tmp, "[[",2)))
  ii=grep("dFdCTP", colnames(la))
  la[,ii]=as.numeric(la[,ii])
  boxplot(la[,ii] ~ la$time)
  dm[[sheet]]=la[,coi]
}
dm=list(nM1000=do.call(rbind,dm[1:2]), nM100=do.call(rbind,dm[3:4]))
for(dm_ in dm){
  boxplot(dm_[,grep("dFdCTP", colnames(dm_))] ~ paste0(dm_$time,"-",dm_$cellType),horizontal = T,ylab="",las=2,col=c("red","blue"))
}


## Save as input to ODE fitting


dmx=list()
for(replicate in 1:3){
  for (dat in c("nM1000","nM100")){
    dm_=dm[[dat]][grep(paste0("-",replicate,"$"),dm[[dat]]$Sample),]
    rownames(dm_)=dm_$Sample
    # ## Mean across replicates:
    # dm_=grpstats(dm[[dat]][,c("time","dFdCTP...ng.mL.")],paste0(dm[[dat]]$time,"-",dm[[dat]]$cellType),"mean")$mean
    dmx[[paste0(dat,"_2N-",replicate)]]=as.data.frame(dm_[grep("2N",rownames(dm_)),])
    dmx[[paste0(dat,"_4N-",replicate)]]=as.data.frame(dm_[grep("4N",rownames(dm_)),])
  }
  for(cellType in c("2N","4N")){
    ii=paste0("nM1000_",cellType,"-",replicate)
    ij=paste0("nM100_",cellType,"-",replicate)
    dmx[[ii]][rownames(dmx[[ij]]),c("time","dFdCTP...ng.mL._low")]=dmx[[ij]][,c("time","dFdCTP...ng.mL.")]
    dmx[[ii]]=dmx[[ii]][order(dmx[[ii]]$time),]
    dmx=dmx[-which(names(dmx)==ij)]
  }
}
sapply(names(dmx), function(x) write.table(dmx[[x]],paste0(x,".txt"),sep="\t", quote = F,row.names = F) )
       
# dmx$nM1000_2N[rownames(dmx$nM100_2N),c("time","dFdCTP...ng.mL._low")]=dmx$nM100_2N[,c("time","dFdCTP...ng.mL.")]
# dmx$nM1000_4N[rownames(dmx$nM100_4N),c("time","dFdCTP...ng.mL._low")]=dmx$nM100_4N[,c("time","dFdCTP...ng.mL.")]
# dm_2N=dmx$nM1000_2N
# dm_4N=dmx$nM1000_4N
# dm_2N=dm_2N[order($time),]
# dm_4N=dm_4N[order(dm_4N$time),]
# write.table(dm_2N,"dm_2N.txt",sep="\t", quote = F,row.names = F)
# write.table(dm_4N,"dm_4N.txt",sep="\t", quote = F,row.names = F)
