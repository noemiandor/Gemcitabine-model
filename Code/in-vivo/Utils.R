
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
  joint_post<-joint_post %>% left_join(exp_post %>% dplyr::select(cell,seg,phi_mle),by=c("cell"="cell","seg"="seg")) %>%
    mutate(cnv_integer=round(log2(phi_mle)))
  joint_post$cnv_integer[which(joint_post$p_cnv<0.9)]<-0
  
  cnv_matrix<-joint_post %>% dplyr::select(cell,region,cnv_Final) %>%
    tidyr::pivot_wider(names_from = region, values_from=cnv_Final) 
  
  cnv_matrix_integer<-joint_post %>% dplyr::select(cell,region,cnv_integer) %>%
    tidyr::pivot_wider(names_from = region, values_from=cnv_integer)  
  
  return(list(cnv_matrix = cnv_matrix, cnv_matrix_integer=cnv_matrix_integer))
}


# DATASETID="Melanoma_GSE174401"
# DATASETID="HNSC_GSE181919"
## Returns Numbat copy number calls in a format that can be used as input to ALFA-K
NumbatPostProcess <- function(DATASETID="CNV", mpoi=NULL, path2karyo="/Users/4482173/Documents/Project/GBM/data/", gBandedKaryo2align = NULL, scRNAseqCells2align = NULL, lambda_conflict=0.25){
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
    laa=getCNVmatrix(paste0(path2numbat,patient))
    if(!useLogFC || !is.null(gBandedKaryo2align)){
      la=as.data.frame(laa$cnv_matrix)
    }else{
      la=as.data.frame(laa$cnv_matrix_integer)
    }
    rownames(la)=la$cell
    sample=sample[rownames(la)]
    
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
    if(!useLogFC && is.null(gBandedKaryo2align)){
      #NA/empty values: unchanged copy number (here we include missing values and "loh" encoding in the original tsv files)
      #1: single copy gain (encoded as ‘amp’ in the original tsv files)
      #-1: single copy loss (encoded as ‘del’ in the original tsv files)
      #2: duplicate copies, (coded as 'bamp' in the original tsv files). Be careful when trying to find cellular karyotypes, we should double the number of copies and not just add 2 copies.
      cn[T]=0
      cn[la=="amp"]=1
      cn[la=="bamp"]=2
      cn[la=="del"]=-1
      
      
      
      
      ## center to ploidy
      cn=cn+ploidies[patient];
      cn[is.na(cn)]=ploidies[patient];
      cn=as.matrix(cn)
      # p<-gplots::heatmap.2(cn,trace='n',cexCol = 0.4,symbreaks = F,symkey=F)
      # pdf(paste0(path2figures,,"/",patient,"ploidy_center.pdf"),width = 5, height= 5)
      # on.exit(dev.off(), add = TRUE)
      # print(p)
      # dev.off()
    }
    
    anno=parseLOCUS(colnames(cn))
    anno=addChrAnnotation(anno)
    anno=as.data.frame(anno)
    rownames(anno)=paste0(anno$chr,":",anno$startpos,"-",anno$endpos)
    
    # ## set copy number of chromosome to copy number of largest segment for that chromosome
    # segments= sapply(sapply(strsplit(colnames(cn),":"),"[[",2), function(x) strsplit(x[[1]],"-")[[1]],simplify = F)
    # segments= as.data.frame(do.call(rbind,sapply(segments, as.numeric,simplify = F)))
    # rownames(segments) = colnames(cn)
    # colnames(segments) = c("start","end")
    # segments$length=1+segments$end-segments$start
    # segments$chr = as.numeric(sapply(strsplit(colnames(cn),":"),"[[",1))
    # chrsegments=sapply(unique(segments$chr), function(x) segments[segments$chr==x,,drop=FALSE],simplify = F)
    # chrsegments=sapply(chrsegments, function(x) x[which.max(x$length),,drop=F],simplify = F)
    # chrsegments = do.call(rbind,chrsegments)
    # cn=cn[,rownames(chrsegments)]
    # colnames(cn)=chrsegments$chr
    
    ## all other chromosomes have copy number equal to ploidy/2 for all cells (assuming signal is not there because of copy number loss)
    otherchr = setdiff(rownames(anno),colnames(cn))
    if(is.null(gBandedKaryo2align)){
      cn_ = matrix(ploidies[patient]/2,nrow(cn),length(otherchr))
    }else{
      cn_ = matrix("neu",nrow(cn),length(otherchr))
    }
    colnames(cn_) <- otherchr
    cn = cbind(cn,cn_)
    # pdf(paste0(path2figures,matlab::filesep,patient,".pdf"))
    # on.exit(dev.off(), add = TRUE)
    # gplots::heatmap.2(cn,trace='n',symbreaks = F,symkey=F)
    # dev.off()
    # graphics.off()
    
    
    
    ## Parse cell names for matching to CLONEID DB entries
    ii = grep("2N",sample)
    ii = setdiff(ii, grep("Cell-Culture",sample))
    sample[ii] = paste0(sample[ii],"-HM")
    
    
    ## align to karyo if exists
    theout=list(cn=cn, cells=sample, anno=anno)
    if(!is.null(gBandedKaryo2align)){
      jnt=alignCNmatrices(gBandedKaryo2align, theout, arm_level_karyo = F)
      ii2align=names(theout$cells)[theout$cells==scRNAseqCells2align]
      
      ## identify rules of mapping neu, amp, del, loh, etc to absolute copy number bu aligning to copy number distribution from karyotyping: 
      alg=map_scRNA_to_absolute_cn(jnt$cn_scRNAseq[ii2align,], jnt$cn_karyo, karyo_orig =  gBandedKaryo2align, amp_max_increase = 3, lambda_conflict = lambda_conflict)
      # alg=map_scRNA_to_absolute_cn(jnt$cn_scRNAseq[ii2align,], jnt$cn_karyo, karyo_orig =  gBandedKaryo2align, amp_max_increase = 3, lambda_conflict = lambda_conflict)
      tmp=rbind(alg$cn_matrix[sample(ii2align,255),],jnt$cn_karyo)
      # hm=heatmap.2(t(tmp), Rowv = NULL, trace = "n")
      
      # apply rules to all cells:
      alg$fixed <- apply_cn_rules_to_scRNA(jnt$cn_scRNAseq, alg$rules)
      # tmp=rbind(alg$fixed[sample(rownames(alg$fixed),255),],jnt$cn_karyo)
      tmp=rbind(alg$fixed[sample(ii2align,255),],jnt$cn_karyo)
      col=rep("blue",nrow(tmp))
      col[rownames(tmp) %in% rownames(jnt$cn_karyo)] = "red"
      hm=heatmap.2(tmp, Colv = NULL, trace = "n", RowSideColors = col, hclustfun = function(x) hclust(x, method = "ward.D"), margins = c(15,5))
      
      theout$alignment=alg
      
      theout$cn=alg$fixed
    }
    
    outputs[[patient]] = theout
    
    # ## Read expression data:
    # la=readRDS(paste0(path2numbat,patient,"/",patient,".count.rds"))
    
  }
  return(outputs)
}


parseLOCUS <- function(loci) {
  chr = as.numeric(sapply(strsplit(loci, ":"), "[[", 1))
  startend = sapply(strsplit(loci, ":"), "[[", 2)
  startp = as.numeric(sapply(strsplit(startend, "-"), "[[", 1))
  endp = as.numeric(sapply(strsplit(startend, "-"), "[[", 2))
  seglength = 1 + endp - startp
  dm = cbind(chr, startp, endp, seglength) 
  colnames(dm) = c("chr", "startpos", "endpos", "seglength")
  return(dm)  
}           


addChrAnnotation <- function(y_mat){
  ## Expected chromosome lengths ##
  x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", 
             col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
  
  
  # --- Start Processing ---
  
  # A. Get Chromosome Arm Boundaries from x
  x[, arm := substr(name, 1, 1)]
  x[, chr_std := chrom] # Use existing 'chr' column name if already standardized
  # Handle potential numeric chr names in x if necessary (our example uses 'chr' prefix)
  # x[grep("^\\d+$", chrom), chr_std := paste0("chr", chrom)]
  # x[chrom == "X", chr_std := "chrX"]
  # x[chrom == "Y", chr_std := "chrY"]
  
  # Calculate boundaries for p and q arms for all chromosomes in x
  arm_boundaries <- x[arm %in% c('p', 'q'),
                      .(p_start = min(fifelse(arm == 'p', chromStart, Inf)), # Min start for p-bands
                        p_end   = max(fifelse(arm == 'p', chromEnd, -Inf)), # Max end for p-bands
                        q_start = min(fifelse(arm == 'q', chromStart, Inf)), # Min start for q-bands
                        q_end   = max(fifelse(arm == 'q', chromEnd, -Inf))  # Max end for q-bands
                      ), by = .(chr_std)]
  
  # Clean up infinite values if an arm is entirely missing in x
  arm_boundaries[is.infinite(p_start), p_start := NA]
  arm_boundaries[is.infinite(p_end), p_end := NA]
  arm_boundaries[is.infinite(q_start), q_start := NA]
  arm_boundaries[is.infinite(q_end), q_end := NA]
  
  # B. Initial Processing of y (Assign Arms)
  # Prepare x_arms for overlap (only p and q bands)
  x_arms <- x[arm %in% c('p', 'q')]
  setkey(x_arms, chr_std, chromStart, chromEnd)
  
  # Prepare y
  y <- as.data.table(y_mat)
  y[, y_id := .I]
  if (!any(grepl("^chr", y$chr))) {
    y[, chr_std_y := paste0("chr", chr)]
  } else {
    y[, chr_std_y := as.character(chr)] # Ensure character
  }
  y[chr == 23, chr_std_y := "chrX"]
  y[chr == 24, chr_std_y := "chrY"]
  setkey(y, chr_std_y, startpos, endpos)
  
  # Find overlaps
  overlaps <- foverlaps(y, x_arms,
                        by.x = c("chr_std_y", "startpos", "endpos"),
                        by.y = c("chr_std", "chromStart", "chromEnd"),
                        type = "any", mult = "all", nomatch = NULL)
  
  # Initialize y with 'unknown' arm - will be overwritten if overlaps found
  y[, arm := "unknown"]
  
  if (!is.null(overlaps) && nrow(overlaps) > 0) {
    # Calculate overlap length
    overlaps[, overlap_len := pmax(0, pmin(endpos, chromEnd) - pmax(startpos, chromStart))]
    
    # Aggregate overlap length by original segment and arm
    arm_overlaps <- overlaps[, .(total_overlap = sum(overlap_len)), by = .(y_id, arm)] # Use arm from x_arms
    
    # Reshape data
    arm_summary <- dcast(arm_overlaps, y_id ~ arm, value.var = "total_overlap", fill = 0)
    
    # Ensure p and q columns exist
    if (!"p" %in% names(arm_summary)) arm_summary[, p := 0]
    if (!"q" %in% names(arm_summary)) arm_summary[, q := 0]
    
    # Assign Arm Classification
    arm_summary[, arm_assigned := fifelse(p > 0 & q > 0, "whole",
                                          fifelse(p > 0, "p",
                                                  fifelse(q > 0, "q", "unknown")))]
    
    # Merge assigned arm back to y
    # Need to handle potential duplicates in y_id if merge creates them (shouldn't with y_id)
    y <- merge(y, arm_summary[, .(y_id, arm_assigned)], by = "y_id", all.x = TRUE)
    
    # Update the 'arm' column, keeping 'unknown' for those with no p/q overlap found
    y[!is.na(arm_assigned), arm := arm_assigned]
    y[, arm_assigned := NULL] # Remove temporary merge column
    
  } else {
    warning("No overlaps found between y segments and p/q arms in x. Initial arms set to 'unknown'.")
  }
  
  
  # C. Identify and Create Missing Arm Placeholders
  # Find which chromosome arms are represented in the processed y
  # Using original numeric 'chr' column for grouping
  y_arm_presence <- y[, .(has_p = any(arm %in% c('p', 'whole')),
                          has_q = any(arm %in% c('q', 'whole'))),
                      by = .(chr, chr_std_y)] # Keep both chr formats
  
  
  # List to store the placeholders for missing arms
  missing_arms_list <- list()
  
  # Iterate through chromosomes present in y
  for (i in 1:nrow(y_arm_presence)) {
    current_chr_num <- y_arm_presence$chr[i]
    current_chr_std <- y_arm_presence$chr_std_y[i]
    has_p_coverage <- y_arm_presence$has_p[i]
    has_q_coverage <- y_arm_presence$has_q[i]
    
    # Look up boundaries for this chromosome
    boundaries <- arm_boundaries[chr_std == current_chr_std]
    
    if (nrow(boundaries) == 1) { # Ensure boundaries were found
      # Check if p arm is missing
      if (!has_p_coverage && !is.na(boundaries$p_start) && !is.na(boundaries$p_end)) {
        message(paste("Creating placeholder for missing 'p' arm on chromosome", current_chr_num))
        p_placeholder <- data.table(
          chr = current_chr_num,
          startpos = boundaries$p_start,
          endpos = boundaries$p_end,
          seglength = boundaries$p_end - boundaries$p_start,
          arm = "p"
          # y_id = NA, # Or assign a unique negative ID? Let's omit for now.
          # chr_std_y = current_chr_std # Keep consistent? Added below if needed.
        )
        missing_arms_list[[length(missing_arms_list) + 1]] <- p_placeholder
      }
      
      # Check if q arm is missing
      if (!has_q_coverage && !is.na(boundaries$q_start) && !is.na(boundaries$q_end)) {
        message(paste("Creating placeholder for missing 'q' arm on chromosome", current_chr_num))
        q_placeholder <- data.table(
          chr = current_chr_num,
          startpos = boundaries$q_start,
          endpos = boundaries$q_end,
          seglength = boundaries$q_end - boundaries$q_start,
          arm = "q"
        )
        missing_arms_list[[length(missing_arms_list) + 1]] <- q_placeholder
      }
    } else {
      warning(paste("Could not find arm boundary information in 'x' for chromosome:", current_chr_std))
    }
  }
  
  # D. Append Placeholders
  # Remove temporary columns from y before binding
  y[, c("y_id", "chr_std_y") := NULL]
  
  if (length(missing_arms_list) > 0) {
    # Combine all placeholders
    missing_arms_dt <- rbindlist(missing_arms_list)
    
    # Ensure column order consistency (optional, but good practice)
    # setcolorder(missing_arms_dt, names(y)) # This might fail if columns aren't exact matches yet
    
    # Append to the original y data
    # Use fill=TRUE to handle potential missing columns (like y_id if it wasn't removed)
    y_final <- rbindlist(list(y, missing_arms_dt), use.names = TRUE, fill = TRUE)
    message(paste("Appended", nrow(missing_arms_dt), "placeholder arm segments."))
  } else {
    y_final <- y # No placeholders to add
    message("No missing arms found to append.")
  }
  
  # Clean up potential NA values from fill=TRUE if necessary (e.g., y_id column if kept)
  # y_final[is.na(y_id), y_id := -seq_len(.N)] # Example: assign negative IDs
  
  # --- Display Final Result ---
  print(y_final)
  
  # Display Arm Boundaries for reference
  print("Calculated Arm Boundaries:")
  print(arm_boundaries)
  return(y_final)
}


alignCNmatrices <- function(cn_karyo, scRNAseq, arm_level_karyo=T){
  
  a=colnames(cn_karyo)
  b=scRNAseq$anno
  B = scRNAseq$cn
  A = cn_karyo
  
  cat("Original Matrix A dimensions:", dim(A), "\n")
  cat("Original Matrix A colnames:", colnames(A), "\n\n")
  cat("Original Matrix B dimensions:", dim(B), "\n")
  cat("Original Matrix B colnames:", colnames(B), "\n\n")
  
  
  # --- 1. & 2. Define Target Columns and Create Mapping ---
  map_list <- list()
  special_cols <- c("Cell", "MARKER")
  
  # Iterate through each high-resolution segment in b
  for (i in 1:nrow(b)) {
    b_segment_name <- rownames(b)[i]
    chr_num <- b$chr[i]
    arm_type <- b$arm[i]
    
    # Handle potential X chromosome naming
    chr_char <- as.character(chr_num)
    # if (chr_num == 23) chr_char <- "X" # Adapt if needed
    
    original_a_col_p <- paste0(chr_char, "p")
    original_a_col_q <- paste0(chr_char, "q")
    
    if (arm_type == "p") {
      map_list[[length(map_list) + 1]] <- data.table(
        target_col = paste0(b_segment_name, "_p"),
        original_a_col = original_a_col_p,
        original_b_col = b_segment_name
      )
    } else if (arm_type == "q") {
      map_list[[length(map_list) + 1]] <- data.table(
        target_col = paste0(b_segment_name, "_q"),
        original_a_col = original_a_col_q,
        original_b_col = b_segment_name
      )
    } else if (arm_type == "whole") {
      # Check if corresponding p/q arms exist in the low-res annotation 'a'
      p_exists_in_a <- original_a_col_p %in% a
      q_exists_in_a <- original_a_col_q %in% a
      whole_exists_in_a <- chr_char %in% a
      
      # If p arm exists, create a target column for it
      if (p_exists_in_a) {
        target_p_name <- paste0(b_segment_name, "_p")
        map_list[[length(map_list) + 1]] <- data.table(
          target_col = target_p_name,
          original_a_col = original_a_col_p,
          original_b_col = b_segment_name # Maps back to the 'whole' column in B
        )
      }
      # If q arm exists, create a target column for it
      if (q_exists_in_a) {
        target_q_name <- paste0(b_segment_name, "_q")
        map_list[[length(map_list) + 1]] <- data.table(
          target_col = target_q_name,
          original_a_col = original_a_col_q,
          original_b_col = b_segment_name # Maps back to the 'whole' column in B
        )
      }
      # If whole chr exists, create a target column for it
      if (whole_exists_in_a) {
        target_wc_name <- paste0(b_segment_name, "_whole")
        map_list[[length(map_list) + 1]] <- data.table(
          target_col = target_wc_name,
          original_a_col = chr_char,
          original_b_col = b_segment_name # Maps back to the 'whole' column in B
        )
      }
      # If neither exists (unlikely for standard chromosomes), skip? Or log warning.
      if (!p_exists_in_a && !q_exists_in_a && !whole_exists_in_a) {
        warning(paste("Segment", b_segment_name, "is 'whole', but neither",
                      original_a_col_p, "nor", original_a_col_q, "found in 'a'. Skipping."))
      }
    }
  }
  
  # Combine map entries
  mapping_table <- rbindlist(map_list)
  
  # Define the final target column names
  target_colnames <- unique(c(special_cols, mapping_table$target_col))
  
  cat("--- Mapping Created ---\n")
  print(mapping_table)
  cat("\nTarget Columns (", length(target_colnames), "):", target_colnames, "\n\n")
  
  
  # --- 3. Reconstruct Matrix A ---
  cat("--- Reconstructing Matrix A ---\n")
  # Create the new matrix structure with NAs
  A_new <- matrix(NA_real_,
                  nrow = nrow(A),
                  ncol = length(target_colnames),
                  dimnames = list(rownames(A), target_colnames))
  
  # Copy special columns directly
  for (scol in special_cols) {
    if (scol %in% colnames(A) && scol %in% colnames(A_new)) {
      A_new[, scol] <- A[, scol]
      cat("Copied special column:", scol, "\n")
    }
  }
  
  if(!arm_level_karyo){
    mapping_table$original_a_col = gsub("q","",gsub("p","",mapping_table$original_a_col))
  }
  
  # Populate other columns using the mapping table
  # Iterate through unique low-resolution columns needed from A
  a_cols_needed <- unique(mapping_table$original_a_col)
  for (a_col in unique(a_cols_needed)) {
    # Find all target columns that map to this original A column
    target_cols_for_a <- mapping_table[original_a_col == a_col, target_col]
    
    # Check if the column exists in the original A
    if (a_col %in% colnames(A)) {
      cat("Mapping data from A[,'", a_col, "'] to A_new columns: ", paste(target_cols_for_a, collapse=", "), "\n", sep="")
      # Copy data from the original A column to all corresponding new target columns
      # Matrix subsetting handles the duplication automatically if target_cols_for_a has multiple elements
      A_new[, target_cols_for_a] <- A[, a_col]
    } else {
      warning(paste("Original column", a_col, "needed for mapping not found in Matrix A. Target columns will be NA."))
      # Columns in A_new remain NA by default
    }
  }
  cat("Finished Reconstructing Matrix A.\n\n")
  
  
  # --- 4. Reconstruct Matrix B ---
  cat("--- Reconstructing Matrix B ---\n")
  # Create the new matrix structure with NAs
  B_new <- matrix(NA_real_,
                  nrow = nrow(B),
                  ncol = length(target_colnames),
                  dimnames = list(rownames(B), target_colnames))
  
  # Handle special columns (initialize, e.g., with 0 or NA)
  # Here initializing with NA, adjust if 0 or another value is more appropriate
  for (scol in special_cols) {
    if (scol %in% colnames(B_new)) {
      B_new[, scol] <- NA # Or 0.0, depending on expected data type/meaning
      cat("Initialized special column:", scol, "with NA\n")
    }
  }
  
  
  # Populate other columns using the mapping table
  # Iterate through unique high-resolution columns needed from B
  b_cols_needed <- unique(mapping_table$original_b_col)
  for (b_col in b_cols_needed) {
    # Find all target columns that map to this original B column
    # (this handles the 1-to-many mapping for split 'whole' columns)
    target_cols_for_b <- mapping_table[original_b_col == b_col, target_col]
    
    # Check if the column exists in the original B
    if (b_col %in% colnames(B)) {
      cat("Mapping data from B[,'", b_col, "'] to B_new columns: ", paste(target_cols_for_b, collapse=", "), "\n", sep="")
      # Copy data from the original B column to all corresponding new target columns
      # Matrix subsetting handles the duplication automatically
      B_new[, target_cols_for_b] <- B[, b_col]
    } else {
      warning(paste("Original column", b_col, "needed for mapping not found in Matrix B. Target columns will be NA."))
      # Columns in B_new remain NA
    }
  }
  cat("Finished Reconstructing Matrix B.\n\n")
  
  A_new = apply(A_new,2,as.numeric)
  A_new = A_new[,!colnames(A_new) %in% c("Cell", "MARKER" )]
  B_new = B_new[,!colnames(B_new) %in% c("Cell", "MARKER" )]
  rownames(A_new) = rownames(A)
  rownames(B_new) = rownames(B)
  
  # --- Final Check ---
  cat("Final Matrix A_new dimensions:", dim(A_new), "\n")
  # print(head(A_new[, 1:min(10, ncol(A_new))])) # Print head of first few columns
  cat("Final Matrix B_new dimensions:", dim(B_new), "\n")
  # print(head(B_new[, 1:min(10, ncol(B_new))])) # Print head of first few columns
  
  # Verify columns are identical
  all_cols_match <- identical(colnames(A_new), colnames(B_new))
  cat("\nColumn names of A_new and B_new are identical:", all_cols_match, "\n")
  
  # Optional: Assign the new matrices back to A and B
  # A <- A_new
  # B <- B_new
  # rm(A_new, B_new, mapping_table) # Clean up intermediate objects
  return(list(cn_karyo=A_new, cn_scRNAseq=B_new))
}



# ------------------------------------------------------------------------------
# apply_cn_rules_to_scRNA
# ------------------------------------------------------------------------------
# Inputs:
#   scRNA_df_new : data.frame/matrix [cells x segments] with values in
#                  {"neu","loh","amp","del","bamp","bdel"} (case-insensitive; extra states ignored)
#   rules_df     : data.frame with rownames = segment IDs (matching column names in scRNA_df_new),
#                  columns = c("neu","loh","amp","del","bamp","bdel"), integer CNs per segment
#
# Returns:
#   Integer matrix [cells x segments_in_overlap] with absolute CN per cell and segment.
#   Segments without rules are dropped (or kept as NA if keep_all_cols=TRUE).
#
# Options:
#   keep_all_cols : if TRUE, keep all columns from scRNA_df_new; columns without rules become NA
#   verbose       : print a brief summary of how many segments were mapped
# ------------------------------------------------------------------------------
apply_cn_rules_to_scRNA <- function(scRNA_df_new,
                                    rules_df,
                                    keep_all_cols = FALSE,
                                    verbose = TRUE) {
  # coerce inputs
  if (is.matrix(scRNA_df_new)) scRNA_df_new <- as.data.frame(scRNA_df_new, stringsAsFactors = FALSE)
  stopifnot(is.data.frame(scRNA_df_new), is.data.frame(rules_df))
  
  # required rule columns (any missing will be ignored safely)
  rule_states <- c("neu","loh","amp","del","bamp","bdel")
  
  # align segments
  segs_rules <- rownames(rules_df)
  if (is.null(segs_rules)) stop("rules_df must have rownames = segment IDs.")
  segs_new   <- colnames(scRNA_df_new)
  if (is.null(segs_new)) stop("scRNA_df_new must have column names = segment IDs.")
  
  overlap <- intersect(segs_new, segs_rules)
  
  if (length(overlap) == 0L) {
    stop("No overlapping segments between scRNA_df_new columns and rules_df rownames.")
  }
  
  # decide output columns
  out_cols <- if (keep_all_cols) segs_new else overlap
  
  # helper: standardize state labels (lower-case, trim)
  norm_states <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    tolower(x)
  }
  
  # prepare output matrix
  out <- matrix(NA_integer_, nrow = nrow(scRNA_df_new), ncol = length(out_cols),
                dimnames = list(rownames(scRNA_df_new), out_cols))
  
  # map only for segments with rules; others remain NA (if keep_all_cols=TRUE)
  for (seg in overlap) {
    # rule row -> named integer vector
    rr <- as.integer(rules_df[seg, intersect(colnames(rules_df), rule_states), drop = TRUE])
    names(rr) <- intersect(colnames(rules_df), rule_states)
    
    # column of states
    v <- norm_states(scRNA_df_new[[seg]])
    
    # vectorized mapping
    mapped <- rep(NA_integer_, length(v))
    for (st in names(rr)) {
      mapped[v == st] <- rr[[st]]
    }
    # write to output (respecting potential keep_all_cols)
    out[, seg] <- mapped
  }
  
  if (verbose) {
    cat(sprintf("[apply_cn_rules_to_scRNA] mapped %d/%d segments (%.1f%% overlap).\n",
                length(overlap), length(segs_new), 100 * length(overlap) / length(segs_new)))
    if (keep_all_cols && length(setdiff(segs_new, overlap)) > 0) {
      cat(sprintf("  %d segments without rules kept as NA (set keep_all_cols=FALSE to drop them).\n",
                  length(setdiff(segs_new, overlap))))
    }
  }
  out
}


# ------------------------------------------------------------------------------
# map_scRNA_to_absolute_cn  (supports neu, loh, amp, del, bamp, bdel)
# ------------------------------------------------------------------------------
# Inputs:
#   scRNA_df : data.frame/matrix [cells x segments], values in
#              {"neu","loh","amp","del","bamp","bdel"} (strings; some may be absent)
#   karyo_df : data.frame/matrix [metaphases x segments] with integer copy numbers
#
# Tunables:
#   amp_max_increase: max steps above CN0 allowed for "amp" (CN0+1..CN0+K)
#   lambda_conflict : weight for conflict penalty (small; enforces sign rules)
#   allow_loh_minus2_if_supported: allow LOH=CN0-2 only if karyotype has mass there
#   allow_bamp_plus2_if_supported: allow bamp=CN0+2 only if karyotype has mass there
#   allow_bdel_minus2_if_supported: allow bdel=CN0-2 only if karyotype has mass there
#
# Returns list:
#   $CN0         : global neutral CN
#   $rules       : data.frame [segments x {neu,loh,amp,del,bamp,bdel}] integer mapping
#   $cn_matrix   : matrix [cells x segments] with assigned absolute CNs
#   $diagnostics : data.frame with distances/penalties per segment
# ------------------------------------------------------------------------------
map_scRNA_to_absolute_cn <- function(scRNA_df,
                                     karyo_df,
                                     karyo_orig,
                                     amp_max_increase = 3,
                                     lambda_conflict = 0.25,
                                     beta_nonneutral = 0,
                                     allow_loh_minus2_if_supported = TRUE,
                                     allow_bamp_plus2_if_supported = TRUE,
                                     allow_bdel_minus2_if_supported = TRUE,
                                     # existing knobs...
                                     min_karyo_cells = 5,
                                     tau_karyo_neutral = 0.85,
                                     tau_karyo_nonneutral_min = 0.10,
                                     r_scRNA_sign_dominance = 2.0,
                                     tau_karyo_side_dominance = 0.60,
                                     enforce_sign_constraints = TRUE,
                                     consensus_correct_uninformed = TRUE,
                                     # NEW per-state informativity knobs
                                     tau_karyo_side_min = 0.15,
                                     min_scRNA_state_mass = 0.01) {
  
  ## ---- helpers -------------------------------------------------------------
  
  to_df <- function(x) {
    if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
    x
  }
  
  ssum <- function(x) { x <- as.numeric(x); if (!length(x)) return(0); sum(x) }
  
  mode_int <- function(v) {
    v <- suppressWarnings(as.integer(v)); v <- v[!is.na(v)]
    if (!length(v)) return(integer(0))
    tb <- table(v, useNA = "no"); tb_num <- as.numeric(tb)
    mx <- max(tb_num); as.integer(names(tb)[tb_num == mx])
  }
  
  prop_table_int <- function(v) {
    v <- suppressWarnings(as.integer(v)); v <- v[!is.na(v)]
    if (!length(v)) return(setNames(numeric(0), character(0)))
    tb <- table(v, useNA = "no"); counts <- as.numeric(tb)
    p <- counts / ssum(counts); names(p) <- names(tb); p
  }
  
  # Wasserstein on non-neutral mass (remove CN0 bin and renormalize)
  w1_nonneutral <- function(kd_full, cd_full, CN0) {
    idx <- which(names(kd_full) == as.character(CN0))
    if (length(idx)) { kd_nn <- kd_full[-idx]; cd_nn <- cd_full[-idx] } else { kd_nn <- kd_full; cd_nn <- cd_full }
    sk <- sum(kd_nn); sc <- sum(cd_nn)
    if (sk == 0 && sc == 0) return(0)
    if (sk > 0) kd_nn <- kd_nn / sk
    if (sc > 0) cd_nn <- cd_nn / sc
    sum(abs(cumsum(kd_nn - cd_nn)))
  }
  
  # Correct Wasserstein on a discrete integer grid (preserves names)
  wasserstein1d_discrete <- function(p_dist, q_dist) {
    p_levels <- as.integer(names(p_dist))
    q_levels <- as.integer(names(q_dist))
    all_levels <- sort(unique(c(p_levels, q_levels)))
    p_vec <- numeric(length(all_levels)); names(p_vec) <- all_levels
    q_vec <- numeric(length(all_levels)); names(q_vec) <- all_levels
    if (length(p_dist)) p_vec[match(p_levels, all_levels)] <- as.numeric(p_dist)
    if (length(q_dist)) q_vec[match(q_levels, all_levels)] <- as.numeric(q_dist)
    sp <- sum(p_vec); if (sp > 0) p_vec <- p_vec / sp
    sq <- sum(q_vec); if (sq > 0) q_vec <- q_vec / sq
    sum(abs(cumsum(p_vec - q_vec)))
  }
  
  build_candidate_dist <- function(state_props, map_cn) {
    if (!length(state_props)) return(setNames(numeric(0), character(0)))
    valid_states <- intersect(names(state_props), names(map_cn))
    if (!length(valid_states)) return(setNames(numeric(0), character(0)))
    cn_vals <- as.integer(map_cn[valid_states])
    probs   <- as.numeric(state_props[valid_states])
    df <- data.frame(cn = cn_vals, p = probs, stringsAsFactors = FALSE)
    agg <- rowsum(df$p, group = df$cn, reorder = FALSE)
    cn_names <- rownames(agg); agg_vec <- as.numeric(agg[, 1]); names(agg_vec) <- cn_names
    tot <- ssum(agg_vec); if (tot > 0) agg_vec <- agg_vec / tot
    agg_vec
  }
  
  conflict_penalty <- function(state_props, map_cn, CN0) {
    gp <- function(nm) if (nm %in% names(state_props)) as.numeric(state_props[nm]) else 0
    p_amp  <- gp("amp");  p_del  <- gp("del")
    p_loh  <- gp("loh");  p_bamp <- gp("bamp"); p_bdel <- gp("bdel")
    pen <- 0
    if (!is.na(map_cn["amp"])  && as.integer(map_cn["amp"])  <= CN0) pen <- pen + p_amp
    if (!is.na(map_cn["del"])  && as.integer(map_cn["del"])  >= CN0) pen <- pen + p_del
    if (!is.na(map_cn["bamp"]) && as.integer(map_cn["bamp"]) <= CN0) pen <- pen + 0.75 * p_bamp
    if (!is.na(map_cn["bdel"]) && as.integer(map_cn["bdel"]) >= CN0) pen <- pen + 0.75 * p_bdel
    if (!is.na(map_cn["loh"])) {
      loh_cn <- as.integer(map_cn["loh"])
      if (loh_cn > CN0)       pen <- pen + 0.5 * p_loh
      if (loh_cn < (CN0 - 2)) pen <- pen + 0.5 * p_loh
    }
    pen
  }
  
  state_proportions <- function(vec) {
    vec <- as.character(vec); vec <- vec[!is.na(vec)]
    if (!length(vec)) return(setNames(numeric(0), character(0)))
    tb <- table(vec, useNA = "no"); counts <- as.numeric(tb)
    p <- counts / ssum(counts); names(p) <- names(tb); p
  }
  
  ## ---- sanitize inputs & align segments ------------------------------------
  scRNA_df <- to_df(scRNA_df); karyo_df <- to_df(karyo_df)
  segs <- intersect(colnames(scRNA_df), colnames(karyo_df))
  if (length(segs) == 0L) stop("No overlapping segment columns between scRNA_df and karyo_df.")
  scRNA_df <- scRNA_df[, segs, drop = FALSE]
  karyo_df <- karyo_df[, segs, drop = FALSE]
  
  ## ---- Step A: Global neutral CN0 ------------------------------------------
  all_karyo_vals <- suppressWarnings(as.integer(unlist(karyo_orig, use.names = FALSE)))
  all_karyo_vals <- all_karyo_vals[!is.na(all_karyo_vals)]
  if (!length(all_karyo_vals)) stop("karyo_orig contains no numeric entries to estimate CN0.")
  modes <- mode_int(all_karyo_vals)
  if (length(modes) == 1L) {
    CN0 <- modes
  } else {
    global_states <- state_proportions(as.character(unlist(scRNA_df, use.names = FALSE)))
    best <- Inf; best_m <- modes[1]
    for (m in modes) {
      tmp_map <- c(neu = m, loh = m, amp = m + 1, del = max(m - 1, 0), bamp = m + 1, bdel = max(m - 1, 0))
      sc <- conflict_penalty(global_states, tmp_map, m)
      if (sc < best) { best <- sc; best_m <- m }
    }
    CN0 <- best_m
  }
  
  ## ---- Step B/C: Per-segment grid search -----------------------------------
  state_names <- c("neu","loh","amp","del","bamp","bdel")
  rule_mat <- matrix(NA_integer_, nrow = length(segs), ncol = length(state_names),
                     dimnames = list(segs, state_names))
  diag_list <- vector("list", length(segs)); names(diag_list) <- segs
  uninformed_segments <- character(0)
  
  CNmax_global <- max(CN0 + amp_max_increase,
                      suppressWarnings(max(all_karyo_vals, na.rm = TRUE)),
                      na.rm = TRUE)
  for (seg in segs) {
      
      kd <- prop_table_int(karyo_df[[seg]])
      sp <- state_proportions(scRNA_df[[seg]])
      for (s in state_names) if (!(s %in% names(sp))) sp[s] <- 0
      sp <- sp[state_names]
      
      # ---- karyotype/scRNA summaries
      k_counts <- suppressWarnings(as.integer(karyo_df[[seg]])); k_counts <- k_counts[!is.na(k_counts)]
      n_k <- length(k_counts)
      
      kd_neu <- if (as.character(CN0) %in% names(kd)) kd[as.character(CN0)] else 0
      kd_nonneutral <- 1 - kd_neu
      above <- as.integer(names(kd))[as.integer(names(kd)) >  CN0]
      below <- as.integer(names(kd))[as.integer(names(kd)) <  CN0]
      kd_amp_side <- if (length(above)) ssum(kd[as.character(above)]) else 0
      kd_del_side <- if (length(below)) ssum(kd[as.character(below)]) else 0
      
      sc_amp  <- as.numeric(sp["amp"])  + as.numeric(sp["bamp"])
      sc_del  <- as.numeric(sp["del"])  + as.numeric(sp["bdel"])
      sc_loh  <- as.numeric(sp["loh"])
      sc_neu  <- as.numeric(sp["neu"])
      
      # segment-level informativity (as you had)
      sc_amp_dom <- sc_amp >= r_scRNA_sign_dominance * max(sc_del, 1e-12)
      sc_del_dom <- sc_del >= r_scRNA_sign_dominance * max(sc_amp, 1e-12)
      kd_amp_dom <- kd_amp_side >= tau_karyo_side_dominance
      kd_del_dom <- kd_del_side >= tau_karyo_side_dominance
      contradictory_seg <- !((sc_amp_dom && kd_amp_dom) || (sc_del_dom && kd_del_dom) || (sc_neu >= 0.5 && kd_neu >= 0.5))
      neutral_dominated <- (kd_neu >= tau_karyo_neutral) || (kd_nonneutral < tau_karyo_nonneutral_min)
      too_few_karyo     <- (n_k < min_karyo_cells)
      karyo_informable_segment <- !(too_few_karyo || neutral_dominated || contradictory_seg)
      
      # ---------- NEW: per-state informativity flags ---------------------------
      # A state is "informable" only if:
      #  - the karyotype has enough mass on the relevant side (or for LOH, enough del-side to support CN<CN0),
      #  - AND scRNA actually exhibits that state in >= min_scRNA_state_mass of cells (otherwise distance can’t learn it),
      #  - AND the sign isn't contradicted by dominant signal on the opposite side.
      inf_amp  <- (kd_amp_side >= tau_karyo_side_min) && (as.numeric(sp["amp"])  >= min_scRNA_state_mass) && !(sc_del_dom && kd_amp_dom)
      inf_bamp <- (kd_amp_side >= tau_karyo_side_min) && (as.numeric(sp["bamp"]) >= min_scRNA_state_mass) && !(sc_del_dom && kd_amp_dom)
      inf_del  <- (kd_del_side >= tau_karyo_side_min) && (as.numeric(sp["del"])  >= min_scRNA_state_mass) && !(sc_amp_dom && kd_del_dom)
      inf_bdel <- (kd_del_side >= tau_karyo_side_min) && (as.numeric(sp["bdel"]) >= min_scRNA_state_mass) && !(sc_amp_dom && kd_del_dom)
      # For LOH, only consider karyotype-informative if there is substantial del-side signal,
      # otherwise LOH is ambiguous and we use fallback.
      inf_loh  <- (kd_del_side >= tau_karyo_side_min) && (sc_loh >= min_scRNA_state_mass)
      
      # ---------- fallback rule (same semantics you used; use base R 'ceiling') --
      fallback_map <- c(
        neu  = CN0,
        loh  = ceiling((CN0 + 0.1)/2),         # your chosen default
        amp  = CN0 + 1,
        del  = max(floor(CN0/2), 0),           # your chosen default
        bamp = CN0 + 1,
        bdel = max(floor(CN0/2), 0)
      )
      
      # If the whole segment is not informable → assign all fallback and skip search
      if (!karyo_informable_segment) {
        rule_mat[seg, ] <- fallback_map
        diag_list[[seg]] <- data.frame(
          segment = seg,
          neu  = fallback_map["neu"],
          loh  = fallback_map["loh"],
          amp  = fallback_map["amp"],
          del  = fallback_map["del"],
          bamp = fallback_map["bamp"],
          bdel = fallback_map["bdel"],
          distance = NA_real_, penalty = NA_real_,
          n_karyo = n_k, kd_neutral = kd_neu, kd_amp_side = kd_amp_side, kd_del_side = kd_del_side,
          sc_amp = sc_amp, sc_del = sc_del,
          karyo_informed = FALSE,
          inf_loh = FALSE, inf_amp = FALSE, inf_del = FALSE, inf_bamp = FALSE, inf_bdel = FALSE,
          fallback_loh = TRUE, fallback_amp = TRUE, fallback_del = TRUE, fallback_bamp = TRUE, fallback_bdel = TRUE,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # ---------- candidate sets (lock non-informable states to fallback) -------
      amp_candidates   <- if (inf_amp)   (CN0 + 1):(CN0 + max(1, amp_max_increase)) else fallback_map["amp"]
      amp_candidates   <- amp_candidates[amp_candidates > CN0]
      bamp_candidates  <- if (inf_bamp)  amp_candidates else fallback_map["bamp"]
      
      del_candidates   <- if (inf_del)   unique(c(max(CN0 - 1, 0), max(CN0 - 2, 0), max(CN0 - 3, 0))) else fallback_map["del"]
      del_candidates   <- del_candidates[del_candidates >= 0]
      bdel_candidates  <- if (inf_bdel)  del_candidates else fallback_map["bdel"]
      
      loh_candidates   <- if (inf_loh)   unique(c(CN0, del_candidates)) else fallback_map["loh"]
      loh_candidates   <- loh_candidates[loh_candidates >= 0]
      
      # --- grid search (only states marked informable can vary) -----------------
      best_obj <- Inf
      best_map <- c(neu = CN0, loh = NA_integer_, amp = NA_integer_, del = NA_integer_, bamp = NA_integer_, bdel = NA_integer_)
      best_parts <- c(distance = NA_real_, penalty = NA_real_)
      
      for (loh_val in loh_candidates) {
        for (amp_val in amp_candidates) {
          for (del_val in del_candidates) {
            for (bamp_val in bamp_candidates) {
              for (bdel_val in bdel_candidates) {
                
                if (enforce_sign_constraints) {
                  if (!(amp_val > CN0 && del_val < CN0 && bamp_val > CN0 && bdel_val < CN0)) next
                }
                
                map_cn <- c(neu = CN0, loh = loh_val, amp = amp_val, del = del_val,
                            bamp = bamp_val, bdel = bdel_val)
                
                cand_dist <- build_candidate_dist(sp, map_cn)
                
                CNmax <- max(CN0 + amp_max_increase,
                             suppressWarnings(max(as.integer(names(cand_dist)), na.rm = TRUE)),
                             suppressWarnings(max(as.integer(names(kd)), na.rm = TRUE)),
                             na.rm = TRUE)
                support <- as.character(0:CNmax)
                kd_full <- numeric(length(support)); names(kd_full) <- support
                cd_full <- kd_full
                if (length(kd))        kd_full[names(kd)]        <- as.numeric(kd)
                if (length(cand_dist)) cd_full[names(cand_dist)] <- as.numeric(cand_dist)
                
                spk <- ssum(kd_full); if (spk > 0) kd_full <- kd_full / spk
                spc <- ssum(cd_full); if (spc > 0) cd_full <- cd_full / spc
                
                dist_total  <- wasserstein1d_discrete(kd_full, cd_full)
                dist_nonneu <- w1_nonneutral(kd_full, cd_full, CN0)
                
                pen_val  <- conflict_penalty(sp, map_cn, CN0)
                dist_val <- dist_total + beta_nonneutral * dist_nonneu
                obj <- dist_val + lambda_conflict * pen_val
                
                if (obj < best_obj) {
                  best_obj <- obj
                  best_map <- map_cn
                  best_parts <- c(distance = dist_val, penalty = pen_val)
                }
              }
            }
          }
        }
      }
      
      # After search, enforce fallback again for any non-informable state (safety net)
      if (!inf_loh)  best_map["loh"]  <- fallback_map["loh"]
      if (!inf_amp)  best_map["amp"]  <- fallback_map["amp"]
      if (!inf_del)  best_map["del"]  <- fallback_map["del"]
      if (!inf_bamp) best_map["bamp"] <- fallback_map["bamp"]
      if (!inf_bdel) best_map["bdel"] <- fallback_map["bdel"]
      
      rule_mat[seg, ] <- best_map
      diag_list[[seg]] <- data.frame(
        segment = seg,
        neu  = as.integer(best_map["neu"]),
        loh  = as.integer(best_map["loh"]),
        amp  = as.integer(best_map["amp"]),
        del  = as.integer(best_map["del"]),
        bamp = as.integer(best_map["bamp"]),
        bdel = as.integer(best_map["bdel"]),
        distance = as.numeric(best_parts["distance"]),
        penalty  = as.numeric(best_parts["penalty"]),
        n_karyo = n_k, kd_neutral = kd_neu, kd_amp_side = kd_amp_side, kd_del_side = kd_del_side,
        sc_amp = sc_amp, sc_del = sc_del,
        karyo_informed = TRUE,
        # per-state informativity + whether fallback was enforced
        inf_loh = inf_loh, inf_amp = inf_amp, inf_del = inf_del, inf_bamp = inf_bamp, inf_bdel = inf_bdel,
        fallback_loh = !inf_loh, fallback_amp = !inf_amp, fallback_del = !inf_del, fallback_bamp = !inf_bamp, fallback_bdel = !inf_bdel,
        stringsAsFactors = FALSE
      )
    }
  
  diagnostics <- do.call(rbind, diag_list)
  rules_df <- as.data.frame(rule_mat, stringsAsFactors = FALSE)
  for (nm in colnames(rules_df)) rules_df[[nm]] <- as.integer(rules_df[[nm]])
  
  # ## ---- Stage 2: consensus correction for uninformed segments ---------------
  # if (consensus_correct_uninformed) {
  #   informed <- rownames(rules_df)[!is.na(diagnostics$karyo_informed) & diagnostics$karyo_informed]
  #   uninformed <- rownames(rules_df)[!is.na(diagnostics$karyo_informed) & !diagnostics$karyo_informed]
  #   
  #   if (length(uninformed) > 0 && length(informed) > 0) {
  #     informed_rules <- rules_df[informed, , drop = FALSE]
  #     modal_vals <- vapply(informed_rules, function(col) {
  #       tab <- table(col, useNA = "no"); as.integer(names(tab)[which.max(tab)])
  #     }, integer(1))
  #     # overwrite only non-neu states; neu already global (CN0)
  #     for (seg in uninformed) {
  #       rules_df[seg, c("loh","amp","del","bamp","bdel")] <- modal_vals[c("loh","amp","del","bamp","bdel")]
  #     }
  #     if (!"consensus_corrected" %in% names(diagnostics)) diagnostics$consensus_corrected <- FALSE
  #     diagnostics$consensus_corrected[diagnostics$segment %in% uninformed] <- TRUE
  #   } else {
  #     diagnostics$consensus_corrected <- FALSE
  #   }
  # } else {
  #   diagnostics$consensus_corrected <- FALSE
  # }
  
  ## ---- Step D: Apply rules --------------------------------------------------
  map_column <- function(states_vec, rule_row_df) {
    if (is.data.frame(rule_row_df)) {
      rr_vals <- as.integer(as.vector(unlist(rule_row_df[1, , drop = TRUE])))
      rr_names <- colnames(rule_row_df)
      rr <- rr_vals; names(rr) <- rr_names
    } else {
      rr <- as.integer(as.vector(unlist(rule_row_df)))
      names(rr) <- names(rule_row_df)
    }
    states_vec <- as.character(states_vec)
    out <- rep(NA_integer_, length(states_vec))
    if ("neu"  %in% names(rr)) out[states_vec == "neu"]  <- rr[["neu"]]
    if ("loh"  %in% names(rr)) out[states_vec == "loh"]  <- rr[["loh"]]
    if ("amp"  %in% names(rr)) out[states_vec == "amp"]  <- rr[["amp"]]
    if ("del"  %in% names(rr)) out[states_vec == "del"]  <- rr[["del"]]
    if ("bamp" %in% names(rr)) out[states_vec == "bamp"] <- rr[["bamp"]]
    if ("bdel" %in% names(rr)) out[states_vec == "bdel"] <- rr[["bdel"]]
    out
  }
  
  cn_mat <- matrix(NA_integer_, nrow = nrow(scRNA_df), ncol = ncol(scRNA_df),
                   dimnames = list(rownames(scRNA_df), colnames(scRNA_df)))
  for (j in seq_along(segs)) {
    seg <- segs[j]
    cn_mat[, j] <- map_column(scRNA_df[[seg]], rules_df[seg, , drop = FALSE])
  }
  
  list(
    CN0 = CN0,
    rules = rules_df,
    cn_matrix = cn_mat,
    diagnostics = diagnostics
  )
}
