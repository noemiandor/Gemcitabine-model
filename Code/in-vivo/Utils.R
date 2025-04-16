
## Extracts cell-by-segment matrix of copy number states from Numbat results
getCNVmatrix<-function(path2jointtsvfile, iteration=NA){
  library(tidyr)
  library(dplyr)
  
  # setwd(path2jointtsvfile)
  # iteration=5 # varies by patients
  # genotype.order<-c(1,2,3,4) # varies by patients
  if(is.na(iteration)){
    iteration=list.files(path2jointtsvfile,pattern="joint_post_")
    split_list <- strsplit(iteration, "[_.]")
    iteration <- sort(as.numeric(sapply(split_list, function(x) x[3])))
    iteration=iteration[length(iteration)]
    # iteration=strsplit(iteration,"_")[[1]]
    # iteration=as.numeric(gsub(".tsv","",iteration[length(iteration)]))
    print(paste("Numbat converged at iteration",iteration,"for patient",path2jointtsvfile))
  }
  
  # read in numbat results
  joint_post<-read.delim(paste0(path2jointtsvfile,"/joint_post_",iteration,".tsv"))
  segs_consensus<-read.delim(paste0(path2jointtsvfile,"/segs_consensus_",iteration,".tsv"))
  clone_post<-read.delim(paste0(path2jointtsvfile,"/clone_post_",iteration,".tsv"))
  exp_post<-read.delim(paste0(path2jointtsvfile,"/exp_post_",iteration,".tsv"))
  
  joint_post = joint_post %>% mutate(CHROM = as.integer(as.character(CHROM)))
  segs_consensus = segs_consensus %>% mutate(CHROM = as.integer(as.character(CHROM)))
  
  # if no multi allelic CNVs
  if (!'n_states' %in% colnames(joint_post)) {
    joint_post = joint_post %>% mutate(
      n_states = ifelse(cnv_state == 'neu', 1, 0),
      cnv_states = cnv_state
    )
  } else {
    # only keep one record per CNV and color by most likely state
    joint_post = joint_post %>%
      group_by(cell, CHROM, seg_end, seg_start) %>%
      mutate(p_cnv = sum(p_cnv)) %>%
      ungroup() %>%
      distinct(cell, CHROM, seg_end, seg_start, .keep_all = TRUE) %>%
      mutate(cnv_states = ifelse(n_states > 1, cnv_state_map, cnv_state))########cnv_state -> cnv_states
  }
  
  
  clone_post$clone_opt<-as.numeric(clone_post$clone_opt)
  # order cells by genotypes
  genotype.order<-sort(unique(clone_post$clone_opt)) # varies by patients
  genotype.ordered.cell<-c()
  for(i in genotype.order){
    genotype.cells<-clone_post %>% filter(clone_opt==i)
    genotype.cells$Final<-genotype.cells[,paste("p",i,sep="_")]
    genotype.cells<-genotype.cells %>% arrange(-Final)
    cat("i is", i, " ",nrow(genotype.cells),"\n")
    genotype.ordered.cell<-rbind(genotype.ordered.cell,genotype.cells)
  }
  
  cell_order=rev(genotype.ordered.cell$cell)
  
  # reformat the joint_post table
  joint_post = joint_post %>%
    mutate(cell = factor(cell, cell_order)) %>%
    mutate(cell_index = as.integer(droplevels(cell))) %>%
    mutate(region=paste0(CHROM,":",seg_start,"-",seg_end))
  
  # get the cnv matrix
  joint_post$cnv_Final<-as.character(joint_post$cnv_states)
  joint_post$cnv_Final[which(joint_post$p_cnv<0.9)]<-"neu"
  
  ## The MLE for fold-change in single cells is available in exp_post$phi_mle.
  # see also: https://github.com/kharchenkolab/numbat/issues/147
  ## log2(phi_mle) in exp_post is approximately the copy changes according to  https://github.com/kharchenkolab/numbat/issues/147
  joint_post<-joint_post %>% left_join(exp_post %>% select(cell,seg,phi_mle),by=c("cell"="cell","seg"="seg")) %>%
    mutate(cnv_integer=round(log2(phi_mle)))
  joint_post$cnv_integer[which(joint_post$p_cnv<0.9)]<-0
  
  cnv_matrix<-joint_post %>% select(cell,region,cnv_Final) %>%
    tidyr::pivot_wider(names_from = region, values_from=cnv_Final) 
  
  cnv_matrix_integer<-joint_post %>% select(cell,region,cnv_integer) %>%
    tidyr::pivot_wider(names_from = region, values_from=cnv_integer)  
  
  return(list(cnv_matrix = cnv_matrix, cnv_matrix_integer=cnv_matrix_integer))
}


# DATASETID="Melanoma_GSE174401"
# DATASETID="HNSC_GSE181919"
## Returns Numbat copy number calls in a format that can be used as input to ALFA-K
NumbatPostProcess <- function(DATASETID="CNV", mpoi=NULL, path2karyo="/Users/4482173/Documents/Project/GBM/data/"){
  useLogFC=F
  
  ## patient cohort
  path2numbat=paste0(path2karyo,DATASETID)
  samples=list.dirs(path2numbat,recursive = F,full.names = F)
  ploidies=rep(2,length(samples))
  names(ploidies) =samples
  if(!is.null(mpoi)){
    ploidies = ploidies[mpoi]
  }
  # ploidies["PatientB"]=2.8
  outputs=list()
  for(patient in names(ploidies)){
    print(patient)
    # ## read METADATA downloaded from GEO
    # library(rjson)
    # x <- fromJSON(file="../data/Melanoma_GSE174401/METADATA/all.json")
    # sample = x$covs$sample
    # names(sample)=x$covs$id
    
    ## read METADATA for sample origin directly from Numbat
    la=read.table(paste0(path2numbat,patient,matlab::filesep,patient,".cell4numbat.anno.txt"),header = T)
    names(la)=gsub("cell","Cell",gsub("Hash","Sample",gsub("anno","Sample", names(la))))
    sample =la$sample
    names(sample)=la$Cell
    #### extra Primary and Recurrent sample name
    la_unique <- la[!duplicated(la$Sample), ]
    desired_order<-c('Primary','Recurrent')
    la_unique$Stage <- factor(la_unique$Stage, levels = desired_order)
    la_unique <- la_unique[order(la_unique$Stage),]
    
    ## Numbat results
    la=getCNVmatrix(paste0(path2numbat,patient))
    if(!useLogFC){
      la=as.data.frame(la$cnv_matrix)
    }else{
      la=as.data.frame(la$cnv_matrix_integer)
    }
    rownames(la)=la$cell
    
    if(ncol(la)==2){
      la1<-as.data.frame(matrix(data = NA,nrow = nrow(la),ncol = 1))
      la1$V1<-la[,2]
      colnames(la1)<-colnames(la)[2]
      rownames(la1)<-rownames(la)
      cn=la1
      la<-la1
    }else{
      la=la[,-1]
      cn=la
    }
    
    
    ## convert states to integer CNs:
    if(!useLogFC){
      #NA/empty values: unchanged copy number (here we include missing values and "loh" encoding in the original tsv files)
      #1: single copy gain (encoded as ‘amp’ in the original tsv files)
      #-1: single copy loss (encoded as ‘del’ in the original tsv files)
      #2: duplicate copies, (coded as 'bamp' in the original tsv files). Be careful when trying to find cellular karyotypes, we should double the number of copies and not just add 2 copies.
      cn[T]=0
      cn[la=="amp"]=1
      cn[la=="bamp"]=2
      cn[la=="del"]=-1
    }
    
    
    
    ## center to ploidy
    cn=cn+ploidies[patient];
    cn[is.na(cn)]=ploidies[patient];
    cn=as.matrix(cn)
    # p<-gplots::heatmap.2(cn,trace='n',cexCol = 0.4,symbreaks = F,symkey=F)
    # pdf(paste0(path2figures,,"/",patient,"ploidy_center.pdf"),width = 5, height= 5)
    # on.exit(dev.off(), add = TRUE)
    # print(p)
    # dev.off()
    
    ## set copy number of chromosome to copy number of largest segment for that chromosome
    segments= sapply(sapply(strsplit(colnames(cn),":"),"[[",2), function(x) strsplit(x[[1]],"-")[[1]],simplify = F)
    segments= as.data.frame(do.call(rbind,sapply(segments, as.numeric,simplify = F)))
    rownames(segments) = colnames(cn)
    colnames(segments) = c("start","end")
    segments$length=1+segments$end-segments$start
    segments$chr = as.numeric(sapply(strsplit(colnames(cn),":"),"[[",1))
    chrsegments=sapply(unique(segments$chr), function(x) segments[segments$chr==x,,drop=FALSE],simplify = F)
    chrsegments=sapply(chrsegments, function(x) x[which.max(x$length),,drop=F],simplify = F)
    chrsegments = do.call(rbind,chrsegments)
    
    cn=cn[,rownames(chrsegments)]
    colnames(cn)=chrsegments$chr
    
    # gplots::heatmap.2(cn,trace='n',symbreaks = F,symkey=F)
    
    
    ## all other chromosomes have copy number equal to ploidy for all cells
    otherchr = setdiff(1:22,colnames(cn))
    cn_ = matrix(ploidies[patient],nrow(cn),length(otherchr))
    colnames(cn_)=otherchr
    cn = cbind(cn,cn_)
    # pdf(paste0(path2figures,matlab::filesep,patient,".pdf"))
    # on.exit(dev.off(), add = TRUE)
    # gplots::heatmap.2(cn,trace='n',symbreaks = F,symkey=F)
    # dev.off()
    # graphics.off()
    
    ## karyotype frequency across timepoints
    cn=round(cn)
    Chr_seq<-as.character(seq(1,22,1))
    cn <- cn[, Chr_seq]
    #cn %>% select(all_of(Chr_seq))
    
    outputs[[patient]] = list(cn=cn, cells=sample)
    
    # ## Read expression data:
    # la=readRDS(paste0(path2numbat,patient,"/",patient,".count.rds"))
    
  }
  return(outputs)
}
