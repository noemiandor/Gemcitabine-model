
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
  #samples=list.dirs(path2numbat,recursive = F,full.names = F)
  #samples<-sapply(strsplit(samples, "_"), `[`, 1)
  subdirs <- list.dirs(path = path2numbat, full.names = FALSE, recursive = FALSE)
  prefixes <- sapply(strsplit(subdirs, "_"), `[`, 1)
  samples <- data.frame(
    folder = subdirs,
    prefix = prefixes,
    stringsAsFactors = FALSE
  )
  
  ploidies=rep(2,length(samples))
  names(ploidies) =samples$prefix
  if(!is.null(mpoi)){
    ploidies = ploidies[mpoi]
  }
  # ploidies["PatientB"]=2.8
  outputs=list()
  for(patient in names(ploidies)){
    print(patient)
    SP<-samples[which(samples$prefix == patient),1]
    # ## read METADATA downloaded from GEO
    # library(rjson)
    # x <- fromJSON(file="../data/Melanoma_GSE174401/METADATA/all.json")
    # sample = x$covs$sample
    # names(sample)=x$covs$id
    
    ## read METADATA for sample origin directly from Numbat
    la=read.table(paste0(path2numbat,SP,matlab::filesep,patient,".cell4numbat.anno.txt"),header = T)
    names(la)=gsub("cell","Cell",gsub("Hash","Sample",gsub("anno","Sample", names(la))))
    sample =la$sample
    names(sample)=la$Cell
    #### extra Primary and Recurrent sample name
    la_unique <- la[!duplicated(la$Sample), ]
    desired_order<-c('Primary','Recurrent')
    la_unique$Stage <- factor(la_unique$Stage, levels = desired_order)
    la_unique <- la_unique[order(la_unique$Stage),]
    
    ## Numbat results
    la=getCNVmatrix(paste0(path2numbat,SP))
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
    cn_ = matrix(ploidies[patient]/2,nrow(cn),length(otherchr))
    colnames(cn_)=otherchr
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
    
    outputs[[SP]] = list(cn=cn, cells=sample, anno=anno)
    
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


alignCNmatrices <- function(cn_karyo, scRNAseq){
  
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
      # If neither exists (unlikely for standard chromosomes), skip? Or log warning.
      if (!p_exists_in_a && !q_exists_in_a) {
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
  
  # Populate other columns using the mapping table
  # Iterate through unique low-resolution columns needed from A
  a_cols_needed <- unique(mapping_table$original_a_col)
  for (a_col in a_cols_needed) {
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


## -----------------------------
## Shared helpers for in-vivo scRNA workflows
## -----------------------------

load_in_vivo_config <- function(config_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read in_vivo_config.yaml.", call. = FALSE)
  }
  if (missing(config_path) || is.null(config_path) || !nzchar(config_path)) {
    stop("A config_path must be provided.", call. = FALSE)
  }

  config_path <- normalizePath(config_path, mustWork = TRUE)
  config <- yaml::read_yaml(config_path)

  if (is.null(config$Results_root) || !nzchar(config$Results_root)) {
    stop("Config must define a non-empty 'Results_root'.", call. = FALSE)
  }

  if (!grepl("^/", config$Results_root)) {
    config$Results_root <- normalizePath(
      file.path(dirname(config_path), config$Results_root),
      mustWork = FALSE
    )
  } else {
    config$Results_root <- normalizePath(config$Results_root, mustWork = FALSE)
  }

  attr(config, "config_path") <- config_path
  config
}

get_results_root <- function(config) {
  root <- config$Results_root
  if (is.null(root) || !nzchar(root)) {
    stop("Config does not contain 'Results_root'.", call. = FALSE)
  }
  normalizePath(root, mustWork = FALSE)
}

results_path <- function(config, ...) {
  file.path(get_results_root(config), ...)
}

.ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

sort_maybe_numeric <- function(x) {
  x <- as.character(x)
  ux <- unique(x)
  nums <- suppressWarnings(as.numeric(ux))
  if (!any(is.na(nums))) ux[order(nums)] else sort(ux)
}

save_plot_pdf_png <- function(plot_obj, file_stub, width = 9, height = 7, dpi = 300) {
  ggplot2::ggsave(paste0(file_stub, ".pdf"), plot_obj, width = width, height = height)
  ggplot2::ggsave(paste0(file_stub, ".png"), plot_obj, width = width, height = height, dpi = dpi)
}

resolve_col_case_insensitive <- function(df, candidates) {
  nms <- names(df)
  for (cand in candidates) {
    idx <- which(tolower(nms) == tolower(cand))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  NA_character_
}

infer_in_vivo_sample_type <- function(sample_values) {
  sample_values <- as.character(sample_values)
  out <- rep(NA_character_, length(sample_values))
  keep <- !is.na(sample_values) & nzchar(sample_values)

  out[keep & grepl("^2N", sample_values)] <- "2N-tumor"
  out[keep & grepl("^4N", sample_values)] <- "4N-tumor"
  out[keep & grepl("^A5", sample_values)] <- "4N-tumor"
  out[keep & grepl("^A6", sample_values)] <- "4N-tumor"
  out[keep & sample_values == "2N-Cell-Culture"] <- "2N-cellline"
  out[keep & sample_values == "4N-Cell-Culture"] <- "4N-cellline"

  out
}

standardize_in_vivo_dose <- function(dose_values) {
  dose_chr <- trimws(as.character(dose_values))
  dose_chr[dose_chr %in% c("0", "0mg", "0 mg/kg", "0mg/kg")] <- "0mg/kg"
  dose_chr[dose_chr %in% c("30", "30mg", "30 mg/kg", "30mg/kg")] <- "30mg/kg"
  dose_chr[dose_chr %in% c("120", "120mg", "120 mg/kg", "120mg/kg")] <- "120mg/kg"
  dose_chr[dose_chr == ""] <- NA_character_
  dose_chr
}

normalize_in_vivo_sample_id <- function(sample_values) {
  sample_chr <- trimws(as.character(sample_values))
  sample_chr[sample_chr %in% c("", "NA", "NaN", "NULL", "None")] <- NA_character_
  sample_chr <- sub("-Count-HM$", "", sample_chr)
  sample_chr <- sub("-Count$", "", sample_chr)
  sample_chr <- sub("-HM$", "", sample_chr)
  sample_chr
}

infer_in_vivo_initial_ploidy <- function(sample_values, harvest_values = NULL) {
  sample_chr <- trimws(as.character(sample_values))
  harvest_chr <- if (is.null(harvest_values)) rep(NA_character_, length(sample_chr)) else trimws(as.character(harvest_values))
  signal <- paste(sample_chr, harvest_chr)
  out <- rep(NA_character_, length(sample_chr))
  out[grepl("(^|-)2N($|-)", signal, ignore.case = TRUE) | grepl("^2N($|-)", sample_chr, ignore.case = TRUE)] <- "2N"
  out[
    grepl("(^|-)4N($|-)", signal, ignore.case = TRUE) |
      grepl("^4N($|-)", sample_chr, ignore.case = TRUE) |
      grepl("^A[0-9]+-4N($|-)", sample_chr, ignore.case = TRUE)
  ] <- "4N"
  out
}

infer_in_vivo_dose_from_harvest <- function(harvest_values) {
  harvest_chr <- trimws(as.character(harvest_values))
  dose_chr <- rep(NA_character_, length(harvest_chr))
  hit <- grepl("^SUM159-(2N|4N)-[0-9]+-", harvest_chr, ignore.case = TRUE)
  dose_chr[hit] <- sub("^SUM159-(2N|4N)-([0-9]+)-.*$", "\\2", harvest_chr[hit], ignore.case = TRUE)
  standardize_in_vivo_dose(dose_chr)
}

calculate_in_vivo_tgi <- function(
  growth_curve_file,
  sheet = NULL,
  control_dose = "0mg/kg",
  baseline_day = "Day_0",
  final_day = NULL,
  endpoint_days = NULL,
  match_control_by = "initial_ploidy"
) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to calculate TGI from xlsx growth curves.", call. = FALSE)
  }
  if (!file.exists(growth_curve_file)) {
    stop("Growth curve file does not exist: ", growth_curve_file, call. = FALSE)
  }

  sheet_use <- if (is.null(sheet)) 1 else sheet
  growth_raw <- readxl::read_excel(growth_curve_file, sheet = sheet_use)
  growth_raw <- as.data.frame(growth_raw, stringsAsFactors = FALSE)

  harvest_col <- resolve_col_case_insensitive(growth_raw, c("harvest"))
  sample_col <- resolve_col_case_insensitive(growth_raw, c("Sequencing IDs", "Sequencing.IDs", "Sequencing ID", "SequencingIDs"))
  if (is.na(harvest_col)) stop("Growth curve file is missing a harvest column.", call. = FALSE)
  if (is.na(sample_col)) stop("Growth curve file is missing a Sequencing IDs column.", call. = FALSE)

  day_cols <- grep("^Day_[0-9]+$", names(growth_raw), value = TRUE)
  if (length(day_cols) == 0) {
    stop("Growth curve file does not contain Day_* tumor volume columns.", call. = FALSE)
  }
  day_nums <- suppressWarnings(as.numeric(sub("^Day_", "", day_cols)))
  day_cols <- day_cols[order(day_nums)]
  if (!(baseline_day %in% day_cols)) {
    stop("Requested baseline_day is not present in growth curve file: ", baseline_day, call. = FALSE)
  }
  if (!is.null(final_day) && !(final_day %in% day_cols)) {
    stop("Requested final_day is not present in growth curve file: ", final_day, call. = FALSE)
  }
  baseline_time <- suppressWarnings(as.numeric(sub("^Day_", "", baseline_day)))
  endpoint_day_cols <- endpoint_days
  if (is.null(endpoint_day_cols)) {
    endpoint_day_cols <- day_cols[day_nums > baseline_time]
  } else {
    endpoint_day_cols <- as.character(endpoint_day_cols)
    endpoint_day_cols <- endpoint_day_cols[endpoint_day_cols %in% day_cols]
    endpoint_day_cols <- endpoint_day_cols[match(endpoint_day_cols, day_cols, nomatch = 0L) > 0L]
  }
  endpoint_day_cols <- endpoint_day_cols[endpoint_day_cols != baseline_day]

  rows <- seq_len(nrow(growth_raw))
  volume_rows <- lapply(rows, function(i) {
    values <- suppressWarnings(as.numeric(unlist(growth_raw[i, day_cols], use.names = FALSE)))
    names(values) <- day_cols
    baseline_value <- values[[baseline_day]]
    final_day_i <- final_day
    if (is.null(final_day_i)) {
      finite_idx <- which(is.finite(values))
      final_day_i <- if (length(finite_idx) == 0) NA_character_ else day_cols[finite_idx[length(finite_idx)]]
    }
    final_value <- if (!is.na(final_day_i)) values[[final_day_i]] else NA_real_
    finite_auc <- is.finite(day_nums) & is.finite(values) & is.finite(baseline_value)
    auc_delta <- NA_real_
    auc_raw <- NA_real_
    auc_n_days <- sum(finite_auc, na.rm = TRUE)
    auc_start_day <- NA_character_
    auc_end_day <- NA_character_
    if (sum(finite_auc, na.rm = TRUE) >= 2) {
      t <- day_nums[finite_auc]
      y <- values[finite_auc]
      ord <- order(t)
      t <- t[ord]
      y <- y[ord]
      auc_raw <- sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
      y_delta <- y - baseline_value
      auc_delta <- sum(diff(t) * (head(y_delta, -1) + tail(y_delta, -1)) / 2)
      auc_start_day <- paste0("Day_", t[1])
      auc_end_day <- paste0("Day_", t[length(t)])
    }
    row <- data.frame(
      growth_curve_row = i,
      sample_id_raw = as.character(growth_raw[[sample_col]][i]),
      sample_id = normalize_in_vivo_sample_id(growth_raw[[sample_col]][i]),
      harvest = as.character(growth_raw[[harvest_col]][i]),
      initial_ploidy = infer_in_vivo_initial_ploidy(growth_raw[[sample_col]][i], growth_raw[[harvest_col]][i]),
      dose = infer_in_vivo_dose_from_harvest(growth_raw[[harvest_col]][i]),
      baseline_day = baseline_day,
      baseline_volume = baseline_value,
      final_day = final_day_i,
      final_volume = final_value,
      tumor_volume_delta = final_value - baseline_value,
      tumor_volume_auc = auc_raw,
      tumor_volume_auc_delta = auc_delta,
      auc_start_day = auc_start_day,
      auc_end_day = auc_end_day,
      auc_n_days = auc_n_days,
      stringsAsFactors = FALSE
    )
    for (day_col in endpoint_day_cols) {
      suffix <- day_col
      day_value <- values[[day_col]]
      row[[paste0("tumor_volume_", suffix)]] <- day_value
      row[[paste0("tumor_volume_delta_", suffix)]] <- day_value - baseline_value
    }
    row
  })
  out <- do.call(rbind, volume_rows)
  out <- out[!is.na(out$sample_id) & nzchar(out$sample_id) & is.finite(out$tumor_volume_delta), , drop = FALSE]
  out$dose <- standardize_in_vivo_dose(out$dose)
  control_dose <- standardize_in_vivo_dose(control_dose)[1]

  out$control_match_group <- if (identical(match_control_by, "initial_ploidy")) out$initial_ploidy else "all"
  out$matched_control_mean_delta <- NA_real_
  out$matched_control_n <- NA_integer_
  out$matched_control_mean_auc_delta <- NA_real_
  out$matched_control_n_auc <- NA_integer_
  for (day_col in endpoint_day_cols) {
    out[[paste0("matched_control_mean_delta_", day_col)]] <- NA_real_
    out[[paste0("matched_control_n_", day_col)]] <- NA_integer_
  }
  groups <- unique(out$control_match_group[!is.na(out$control_match_group)])
  for (group_value in groups) {
    control_idx <- out$control_match_group == group_value & out$dose == control_dose & is.finite(out$tumor_volume_delta)
    mean_delta <- mean(out$tumor_volume_delta[control_idx], na.rm = TRUE)
    n_control <- sum(control_idx, na.rm = TRUE)
    target_idx <- out$control_match_group == group_value
    out$matched_control_mean_delta[target_idx] <- if (is.finite(mean_delta)) mean_delta else NA_real_
    out$matched_control_n[target_idx] <- n_control

    control_auc_idx <- out$control_match_group == group_value & out$dose == control_dose & is.finite(out$tumor_volume_auc_delta)
    mean_auc_delta <- mean(out$tumor_volume_auc_delta[control_auc_idx], na.rm = TRUE)
    n_control_auc <- sum(control_auc_idx, na.rm = TRUE)
    out$matched_control_mean_auc_delta[target_idx] <- if (is.finite(mean_auc_delta)) mean_auc_delta else NA_real_
    out$matched_control_n_auc[target_idx] <- n_control_auc

    for (day_col in endpoint_day_cols) {
      delta_col <- paste0("tumor_volume_delta_", day_col)
      mean_col <- paste0("matched_control_mean_delta_", day_col)
      n_col <- paste0("matched_control_n_", day_col)
      control_day_idx <- out$control_match_group == group_value & out$dose == control_dose & is.finite(out[[delta_col]])
      mean_day_delta <- mean(out[[delta_col]][control_day_idx], na.rm = TRUE)
      n_day_control <- sum(control_day_idx, na.rm = TRUE)
      out[[mean_col]][target_idx] <- if (is.finite(mean_day_delta)) mean_day_delta else NA_real_
      out[[n_col]][target_idx] <- n_day_control
    }
  }
  out$TGI_percent <- 100 * (1 - out$tumor_volume_delta / out$matched_control_mean_delta)
  out$TGI_percent[!is.finite(out$TGI_percent)] <- NA_real_
  out$TGI_percent_auc <- 100 * (1 - out$tumor_volume_auc_delta / out$matched_control_mean_auc_delta)
  out$TGI_percent_auc[!is.finite(out$TGI_percent_auc)] <- NA_real_
  for (day_col in endpoint_day_cols) {
    delta_col <- paste0("tumor_volume_delta_", day_col)
    mean_col <- paste0("matched_control_mean_delta_", day_col)
    tgi_col <- paste0("TGI_percent_", day_col)
    out[[tgi_col]] <- 100 * (1 - out[[delta_col]] / out[[mean_col]])
    out[[tgi_col]][!is.finite(out[[tgi_col]])] <- NA_real_
  }

  out[order(out$initial_ploidy, out$dose, out$sample_id), , drop = FALSE]
}

make_in_vivo_tgi_process_markdown <- function(
  growth_curve_file,
  output_csv_file,
  analysis_label,
  target_cell_description,
  cluster_description,
  code_location,
  sample_tgi
) {
  format_num <- function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    out <- rep(NA_character_, length(x_num))
    finite <- is.finite(x_num)
    out[finite] <- format(signif(x_num[finite], 6), trim = TRUE, scientific = FALSE)
    out[!finite] <- "NA"
    out
  }
  escape_md <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    gsub("[|]", "\\\\|", x)
  }

  endpoint_days <- sub("^TGI_percent_", "", grep("^TGI_percent_Day_", names(sample_tgi), value = TRUE))
  endpoint_days <- endpoint_days[order(suppressWarnings(as.numeric(sub("^Day_", "", endpoint_days))))]
  endpoint_label <- if (length(endpoint_days) == 0) "None" else paste(endpoint_days, collapse = ", ")
  baseline_label <- paste(unique(sample_tgi$baseline_day), collapse = ", ")
  final_day_label <- paste(unique(sample_tgi$final_day), collapse = ", ")

  control_summary <- sample_tgi[sample_tgi$dose == "0mg/kg", , drop = FALSE]
  control_lines <- if (nrow(control_summary) > 0) {
    control_summary <- control_summary[!duplicated(control_summary$initial_ploidy), , drop = FALSE]
    c(
      "| initial ploidy | control dose | control n | final endpoint mean control delta | AUC mean control delta |",
      "|---|---:|---:|---:|---:|",
      apply(control_summary, 1, function(row) {
        paste0(
          "| ", row[["initial_ploidy"]],
          " | 0mg/kg",
          " | ", row[["matched_control_n"]],
          " | ", format_num(row[["matched_control_mean_delta"]]),
          " | ", format_num(row[["matched_control_mean_auc_delta"]]),
          " |"
        )
      })
    )
  } else {
    "No matched control summary was available."
  }

  endpoint_control_rows <- lapply(endpoint_days, function(day_col) {
    mean_col <- paste0("matched_control_mean_delta_", day_col)
    n_col <- paste0("matched_control_n_", day_col)
    if (!all(c(mean_col, n_col) %in% names(sample_tgi))) return(NULL)
    day_summary <- sample_tgi[sample_tgi$dose == "0mg/kg", c("initial_ploidy", mean_col, n_col), drop = FALSE]
    day_summary <- day_summary[!duplicated(day_summary$initial_ploidy), , drop = FALSE]
    data.frame(
      endpoint_day = day_col,
      initial_ploidy = day_summary$initial_ploidy,
      matched_control_n = day_summary[[n_col]],
      matched_control_mean_delta = day_summary[[mean_col]],
      stringsAsFactors = FALSE
    )
  })
  endpoint_control_summary <- do.call(rbind, endpoint_control_rows)
  endpoint_control_lines <- if (!is.null(endpoint_control_summary) && nrow(endpoint_control_summary) > 0) {
    c(
      "| endpoint day | initial ploidy | control n | matched control mean delta |",
      "|---|---:|---:|---:|",
      apply(endpoint_control_summary, 1, function(row) {
        paste0(
          "| ", row[["endpoint_day"]],
          " | ", row[["initial_ploidy"]],
          " | ", row[["matched_control_n"]],
          " | ", format_num(row[["matched_control_mean_delta"]]),
          " |"
        )
      })
    )
  } else {
    "No endpoint-specific control summary was available."
  }

  export_colnames <- character(0)
  if (!is.null(output_csv_file) && file.exists(output_csv_file)) {
    export_colnames <- tryCatch(
      names(utils::read.csv(output_csv_file, nrows = 0, check.names = FALSE)),
      error = function(e) character(0)
    )
  }
  if (length(export_colnames) == 0) {
    export_colnames <- c(
      "cell_id", "sample_id", "cluster", "initial_ploidy", "gemcitabine_dose",
      "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy", "average_ploidy",
      "median_ploidy", "sample_mean_pseudotime", "sample_median_pseudotime",
      "n_target_cells", "n_tumor_cells", "target_cluster_fraction",
      "target_cluster_percent", "growth_curve_sample_id_raw", "growth_curve_harvest",
      "tgi_initial_ploidy", "tgi_dose", "tumor_volume_baseline_day",
      "tumor_volume_baseline", "tumor_volume_final_day", "tumor_volume_final",
      "tumor_volume_delta", "matched_control_mean_delta", "matched_control_n",
      "TGI_percent", "tumor_volume_auc", "tumor_volume_auc_delta", "auc_start_day",
      "auc_end_day", "auc_n_days", "matched_control_mean_auc_delta",
      "matched_control_n_auc", "TGI_percent_auc"
    )
  }
  describe_export_column <- function(column_name) {
    fixed <- list(
      cell_id = c("cell-level scRNA/scVelo", "Cell barcode or cell ID from the scVelo all-cells metrics table."),
      sample_id = c("cell-level join key", "Normalized sample ID used in the exported fig4 data and for joining sample-level TGI values."),
      cluster = c("cell-level scRNA/scVelo", "Tumor cluster assignment for the exported target cell."),
      initial_ploidy = c("sample-level scRNA metadata", "Initial ploidy group for the sample, `2N` or `4N`, from the scRNA metadata."),
      gemcitabine_dose = c("sample-level scRNA metadata", "Gemcitabine dose label from the scRNA metadata."),
      gemcitabine_dose_mg_per_kg = c("sample-level scRNA metadata", "Numeric Gemcitabine dose in mg/kg."),
      pseudotime = c("cell-level scVelo", "scVelo pseudotime value used to draw `fig4_pseudotime_distribution_per_sample.pdf`."),
      cell_ploidy = c("cell-level ploidy", "Cell-level ploidy value assigned before the 04f/04g export."),
      average_ploidy = c("sample-level summary", "Mean `cell_ploidy` across exported target cells from the same sample."),
      median_ploidy = c("sample-level summary", "Median `cell_ploidy` across exported target cells from the same sample."),
      sample_mean_pseudotime = c("sample-level summary", "Mean pseudotime across exported target cells from the same sample."),
      sample_median_pseudotime = c("sample-level summary", "Median pseudotime across exported target cells from the same sample."),
      n_target_cells = c("sample-level summary", "Number of exported target Tumor cells from the same sample."),
      n_tumor_cells = c("sample-level summary", "Number of Tumor cells with finite pseudotime in the same sample before the 04f/04g target-cluster filter."),
      target_cluster_fraction = c("sample-level summary", "Fraction `n_target_cells / n_tumor_cells` for the same sample."),
      target_cluster_percent = c("sample-level summary", "Percent `100 * target_cluster_fraction` for the same sample."),
      growth_curve_sample_id_raw = c("growth curve/sample-level", "Original `Sequencing IDs` value from the tumor growth-curve spreadsheet."),
      growth_curve_harvest = c("growth curve/sample-level", "Original `harvest` label from the tumor growth-curve spreadsheet."),
      tgi_initial_ploidy = c("growth curve/sample-level", "Initial ploidy inferred by the TGI helper from `Sequencing IDs` and `harvest`; used for matched-control grouping."),
      tgi_dose = c("growth curve/sample-level", "Gemcitabine dose inferred by the TGI helper from `harvest`; standardized to `0mg/kg`, `30mg/kg`, or `120mg/kg`."),
      tumor_volume_baseline_day = c("TGI/sample-level", "Baseline day used for endpoint delta and AUC baseline adjustment; currently `Day_0`."),
      tumor_volume_baseline = c("TGI/sample-level", "Tumor volume at `tumor_volume_baseline_day`."),
      tumor_volume_final_day = c("TGI/sample-level", "Final finite day used for the backward-compatible endpoint `TGI_percent`; currently `Day_38` for the exported samples."),
      tumor_volume_final = c("TGI/sample-level", "Tumor volume at `tumor_volume_final_day`."),
      tumor_volume_delta = c("TGI/sample-level", "`tumor_volume_final - tumor_volume_baseline`; numerator for backward-compatible `TGI_percent`."),
      matched_control_mean_delta = c("TGI/sample-level", "Mean `tumor_volume_delta` among `0mg/kg` controls with the same initial ploidy."),
      matched_control_n = c("TGI/sample-level", "Number of same-initial-ploidy `0mg/kg` controls used for `matched_control_mean_delta`."),
      TGI_percent = c("TGI/sample-level", "Backward-compatible final endpoint TGI: `100 * (1 - tumor_volume_delta / matched_control_mean_delta)`."),
      tumor_volume_auc = c("AUC TGI/sample-level", "Raw trapezoidal AUC of tumor volume over all finite `Day_*` measurements for the sample."),
      tumor_volume_auc_delta = c("AUC TGI/sample-level", "Baseline-adjusted trapezoidal AUC of `tumor_volume_Day_X - tumor_volume_baseline` over all finite `Day_*` measurements."),
      auc_start_day = c("AUC TGI/sample-level", "First finite `Day_*` included in the AUC calculation."),
      auc_end_day = c("AUC TGI/sample-level", "Last finite `Day_*` included in the AUC calculation."),
      auc_n_days = c("AUC TGI/sample-level", "Number of finite `Day_*` tumor-volume measurements included in the AUC calculation."),
      matched_control_mean_auc_delta = c("AUC TGI/sample-level", "Mean `tumor_volume_auc_delta` among `0mg/kg` controls with the same initial ploidy."),
      matched_control_n_auc = c("AUC TGI/sample-level", "Number of same-initial-ploidy `0mg/kg` controls used for `matched_control_mean_auc_delta`."),
      TGI_percent_auc = c("AUC TGI/sample-level", "AUC-based TGI: `100 * (1 - tumor_volume_auc_delta / matched_control_mean_auc_delta)`.")
    )
    if (column_name %in% names(fixed)) return(fixed[[column_name]])
    if (grepl("^tumor_volume_Day_[0-9]+$", column_name)) {
      day <- sub("^tumor_volume_", "", column_name)
      return(c("growth curve/sample-level", paste0("Observed tumor volume at `", day, "`; repeated for every exported cell from the same sample.")))
    }
    if (grepl("^tumor_volume_delta_Day_[0-9]+$", column_name)) {
      day <- sub("^tumor_volume_delta_", "", column_name)
      return(c("endpoint TGI/sample-level", paste0("Endpoint tumor-volume delta at `", day, "`: `tumor_volume_", day, " - tumor_volume_baseline`.")))
    }
    if (grepl("^matched_control_mean_delta_Day_[0-9]+$", column_name)) {
      day <- sub("^matched_control_mean_delta_", "", column_name)
      return(c("endpoint TGI/sample-level", paste0("Mean `tumor_volume_delta_", day, "` among `0mg/kg` controls with the same initial ploidy.")))
    }
    if (grepl("^matched_control_n_Day_[0-9]+$", column_name)) {
      day <- sub("^matched_control_n_", "", column_name)
      return(c("endpoint TGI/sample-level", paste0("Number of same-initial-ploidy `0mg/kg` controls used for `matched_control_mean_delta_", day, "`.")))
    }
    if (grepl("^TGI_percent_Day_[0-9]+$", column_name)) {
      day <- sub("^TGI_percent_", "", column_name)
      return(c("endpoint TGI/sample-level", paste0("Endpoint TGI at `", day, "`: `100 * (1 - tumor_volume_delta_", day, " / matched_control_mean_delta_", day, ")`.")))
    }
    c("exported column", "Column exported by the analysis script; no specific TGI dictionary entry was assigned.")
  }
  column_description_lines <- c(
    "| column | level/source | description |",
    "|---|---|---|",
    vapply(export_colnames, function(column_name) {
      desc <- describe_export_column(column_name)
      paste0("| `", escape_md(column_name), "` | ", escape_md(desc[[1]]), " | ", escape_md(desc[[2]]), " |")
    }, character(1))
  )

  c(
    "# TGI Calculation Process",
    "",
    paste0("Analysis: ", analysis_label),
    "",
    "This document describes how tumor growth inhibition (TGI) was calculated for the exported cell-level data used with `fig4_pseudotime_distribution_per_sample.pdf`.",
    "",
    "Output CSV:",
    "",
    paste0("`", output_csv_file, "`"),
    "",
    "## Input Data",
    "",
    "The tumor growth curve data were read from:",
    "",
    paste0("`", growth_curve_file, "`"),
    "",
    "The analysis used `Sheet1`. Tumor volume columns were detected dynamically from all columns matching `Day_*`.",
    "",
    "Endpoint TGI was calculated for every detected follow-up day after `Day_0`:",
    "",
    paste0("`", endpoint_label, "`"),
    "",
    paste0("Baseline day detected/used: `", baseline_label, "`."),
    paste0("Final day used for backward-compatible `TGI_percent`: `", final_day_label, "`."),
    "",
    "## Sample ID Matching",
    "",
    "Growth-curve sample IDs were normalized before joining to the scRNA-seq sample IDs:",
    "",
    "- trailing `-HM`, `-Count`, and `-Count-HM` suffixes were removed.",
    "- initial ploidy was inferred from `Sequencing IDs` and `harvest`.",
    "- Gemcitabine dose was inferred from the `harvest` label and standardized to `0mg/kg`, `30mg/kg`, or `120mg/kg`.",
    "",
    "Rows were kept for TGI calculation only when the normalized sample ID was non-missing and the final endpoint delta was finite.",
    "",
    "Sample-level TGI values were then joined back to the exported cell-level table by normalized `sample_id`; therefore every cell from the same sample carries the same TGI and tumor-volume columns.",
    "",
    "## Endpoint TGI",
    "",
    "For sample `i`, baseline tumor volume is:",
    "",
    "```text",
    "B_i = V_i,Day_0",
    "```",
    "",
    "For every detected follow-up day `Day_X`:",
    "",
    "```text",
    "delta_i,Day_X = V_i,Day_X - B_i",
    "matched_control_mean_delta_g,Day_X = mean(delta_j,Day_X for controls j with initial_ploidy = g and dose = 0mg/kg)",
    "TGI_percent_Day_X = 100 * (1 - delta_i,Day_X / matched_control_mean_delta_g,Day_X)",
    "```",
    "",
    "Controls were matched separately by initial ploidy group `g`, so 2N samples were compared with 2N `0mg/kg` controls and 4N samples were compared with 4N `0mg/kg` controls.",
    "",
    "The existing `TGI_percent` column is retained for backward compatibility and represents the final endpoint TGI using the final finite day for each sample, currently `Day_38` in this dataset.",
    "",
    "Endpoint-specific matched control means used in this run:",
    "",
    endpoint_control_lines,
    "",
    "## AUC-Based TGI",
    "",
    "AUC was calculated with the trapezoidal rule over all finite `Day_*` tumor-volume points for each sample. If a sample had finite tumor volumes at ordered times `t_1, ..., t_k` with corresponding volumes `V_1, ..., V_k`, raw tumor-volume AUC was:",
    "",
    "```text",
    "tumor_volume_auc = sum over m = 1...(k - 1) of (t_(m+1) - t_m) * (V_m + V_(m+1)) / 2",
    "```",
    "",
    "The AUC-based TGI uses baseline-adjusted volume:",
    "",
    "```text",
    "B_i = V_i,Day_0",
    "tumor_volume_auc_delta = sum over m = 1...(k - 1) of (t_(m+1) - t_m) * ((V_m - B_i) + (V_(m+1) - B_i)) / 2",
    "matched_control_mean_auc_delta_g = mean(tumor_volume_auc_delta_j for controls j with initial_ploidy = g and dose = 0mg/kg)",
    "TGI_percent_auc = 100 * (1 - tumor_volume_auc_delta / matched_control_mean_auc_delta)",
    "```",
    "",
    "`auc_start_day`, `auc_end_day`, and `auc_n_days` record which finite timepoints were used for each sample. The matched control AUC mean was calculated from the same initial-ploidy-matched `0mg/kg` control group.",
    "",
    "## Final Endpoint and AUC Control Means",
    "",
    control_lines,
    "",
    "## Exported CSV Structure and Complete Column Dictionary",
    "",
    paste0("The exported CSV is cell-level. Each row is one ", target_cell_description, "."),
    "",
    paste0("Cluster scope: ", cluster_description),
    "",
    "All sample-level fields, including TGI columns, are repeated across all exported cells from the same sample.",
    "",
    column_description_lines,
    "",
    "## Code Location",
    "",
    "Reusable helper:",
    "",
    "`calculate_in_vivo_tgi()` in `/Users/4482173/Documents/GitHub/Gemcitabine-model/Code/in-vivo/Utils.R`",
    "",
    "Export logic:",
    "",
    paste0("`", code_location, "`")
  )
}

sanitize_path_component <- function(x, prefix = NULL) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | x == ""] <- "NA"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == ""] <- "NA"
  if (!is.null(prefix) && nzchar(prefix)) {
    x <- paste0(prefix, x)
  }
  x
}

write_stackfig_outputs <- function(
  obj,
  output_dir,
  group_col_candidates = c("sample_type", "sample.type", "SampleType", "sampleType"),
  fallback_sample_col_candidates = c("sample", "Sample"),
  cluster_col_candidates = c("cluster", "clusters", "Cluster", "seurat_clusters"),
  group_label = "sample_type",
  cluster_label = "cluster",
  group_title = "Cluster Proportion Within Each sample_type",
  cluster_title = "sample_type Proportion Within Each Cluster",
  group_by_cluster_stub = "stack_sample_type_by_cluster",
  cluster_by_group_stub = "stack_cluster_by_sample_type",
  group_fill_colors = NULL,
  cluster_order_by_group_sum = NULL,
  cluster_order_tiebreak_groups = NULL,
  width = 10,
  height = 6,
  dpi = 300
) {
  obj_slots <- tryCatch(methods::slotNames(obj), error = function(e) character(0))
  if (!inherits(obj, "Seurat") && !("meta.data" %in% obj_slots)) {
    stop("`obj` must be a Seurat-like object with a meta.data slot.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for stack figures.", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for stack figures.", call. = FALSE)
  }

  .ensure_dir(output_dir)

  meta <- obj@meta.data
  if (nrow(meta) == 0) {
    stop("No cells in meta.data.", call. = FALSE)
  }

  group_col <- resolve_col_case_insensitive(meta, group_col_candidates)
  if (!is.na(group_col)) {
    group_vals <- as.character(meta[[group_col]])
  } else {
    sample_col <- resolve_col_case_insensitive(meta, fallback_sample_col_candidates)
    if (is.na(sample_col)) {
      stop("Cannot find a grouping column or fallback sample column in meta.data.", call. = FALSE)
    }
    group_vals <- infer_in_vivo_sample_type(meta[[sample_col]])
  }

  cluster_col <- resolve_col_case_insensitive(meta, cluster_col_candidates)
  if (is.na(cluster_col)) {
    cluster_vals <- as.character(Seurat::Idents(obj))
  } else {
    cluster_vals <- as.character(meta[[cluster_col]])
  }

  group_vals[is.na(group_vals) | group_vals == ""] <- "NA"
  cluster_vals[is.na(cluster_vals) | cluster_vals == ""] <- "NA"

  preferred_group_levels <- c("2N-cellline", "4N-cellline", "2N-tumor", "4N-tumor")
  observed_group_levels <- unique(group_vals)
  group_levels <- c(
    intersect(preferred_group_levels, observed_group_levels),
    setdiff(observed_group_levels, preferred_group_levels)
  )
  cluster_levels <- sort_maybe_numeric(cluster_vals)

  plot_df <- data.frame(
    group = factor(group_vals, levels = group_levels),
    cluster = factor(cluster_vals, levels = cluster_levels),
    stringsAsFactors = FALSE
  )

  tab_group_cluster <- as.data.frame(
    table(group = plot_df$group, cluster = plot_df$cluster),
    stringsAsFactors = FALSE
  )
  names(tab_group_cluster)[names(tab_group_cluster) == "Freq"] <- "count"
  tab_group_cluster$proportion <- with(
    tab_group_cluster,
    count / ave(count, group, FUN = sum)
  )
  names(tab_group_cluster)[names(tab_group_cluster) == "group"] <- group_label
  names(tab_group_cluster)[names(tab_group_cluster) == "cluster"] <- cluster_label

  utils::write.csv(
    tab_group_cluster,
    file = file.path(output_dir, paste0(group_label, "_", cluster_label, "_proportion.csv")),
    row.names = FALSE
  )

  p1 <- ggplot2::ggplot(
    tab_group_cluster,
    ggplot2::aes(
      x = !!rlang::sym(group_label),
      y = proportion,
      fill = !!rlang::sym(cluster_label)
    )
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    ggplot2::labs(
      title = group_title,
      x = group_label,
      y = "Proportion",
      fill = cluster_label
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  save_plot_pdf_png(
    p1,
    file.path(output_dir, group_by_cluster_stub),
    width = width,
    height = height,
    dpi = dpi
  )

  tab_cluster_group <- as.data.frame(
    table(cluster = plot_df$cluster, group = plot_df$group),
    stringsAsFactors = FALSE
  )
  names(tab_cluster_group)[names(tab_cluster_group) == "Freq"] <- "count"
  tab_cluster_group$proportion <- with(
    tab_cluster_group,
    count / ave(count, cluster, FUN = sum)
  )
  names(tab_cluster_group)[names(tab_cluster_group) == "cluster"] <- cluster_label
  names(tab_cluster_group)[names(tab_cluster_group) == "group"] <- group_label

  if (!is.null(cluster_order_by_group_sum) && length(cluster_order_by_group_sum) > 0) {
    cluster_stats_df <- data.frame(
      cluster_chr = as.character(tab_cluster_group[[cluster_label]]),
      group_chr = as.character(tab_cluster_group[[group_label]]),
      proportion = tab_cluster_group$proportion,
      stringsAsFactors = FALSE
    )

    split_cluster_stats <- split(cluster_stats_df, cluster_stats_df$cluster_chr)
    cluster_order_stats <- do.call(
      rbind,
      lapply(names(split_cluster_stats), function(cluster_i) {
        sub_df <- split_cluster_stats[[cluster_i]]
        data.frame(
          cluster_chr = cluster_i,
          primary_score = sum(sub_df$proportion[sub_df$group_chr %in% cluster_order_by_group_sum], na.rm = TRUE),
          tiebreak_score = if (!is.null(cluster_order_tiebreak_groups) && length(cluster_order_tiebreak_groups) > 0) {
            sum(sub_df$proportion[sub_df$group_chr %in% cluster_order_tiebreak_groups], na.rm = TRUE)
          } else {
            0
          },
          stringsAsFactors = FALSE
        )
      })
    )

    base_cluster_levels <- cluster_levels
    cluster_order_stats$base_rank <- match(cluster_order_stats$cluster_chr, base_cluster_levels)
    cluster_levels <- cluster_order_stats$cluster_chr[
      order(-cluster_order_stats$primary_score, -cluster_order_stats$tiebreak_score, cluster_order_stats$base_rank)
    ]
  }

  tab_cluster_group[[cluster_label]] <- factor(tab_cluster_group[[cluster_label]], levels = cluster_levels)
  tab_cluster_group[[group_label]] <- factor(tab_cluster_group[[group_label]], levels = group_levels)

  utils::write.csv(
    tab_cluster_group,
    file = file.path(output_dir, paste0(cluster_label, "_", group_label, "_proportion.csv")),
    row.names = FALSE
  )

  p2 <- ggplot2::ggplot(
    tab_cluster_group,
    ggplot2::aes(
      x = !!rlang::sym(cluster_label),
      y = proportion,
      fill = !!rlang::sym(group_label)
    )
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    ggplot2::labs(
      title = cluster_title,
      x = cluster_label,
      y = "Proportion",
      fill = group_label
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(group_fill_colors)) {
    valid_group_colors <- group_fill_colors[intersect(names(group_fill_colors), group_levels)]
    if (length(valid_group_colors) > 0) {
      p2 <- p2 + ggplot2::scale_fill_manual(values = valid_group_colors, drop = FALSE)
    }
  }

  save_plot_pdf_png(
    p2,
    file.path(output_dir, cluster_by_group_stub),
    width = width,
    height = height,
    dpi = dpi
  )

  invisible(
    list(
      group_cluster = tab_group_cluster,
      cluster_group = tab_cluster_group
    )
  )
}

assert_required_column <- function(df, col_name) {
  if (!(col_name %in% colnames(df))) {
    stop("Required metadata column is missing: ", col_name, call. = FALSE)
  }
}

configure_future_for_seurat <- function(max_size_gb = 30, strategy = "sequential", workers = NULL) {
  max_size_bytes <- as.numeric(max_size_gb) * 1024^3
  current_max <- getOption("future.globals.maxSize")
  if (is.null(current_max) || is.na(current_max) || current_max < max_size_bytes) {
    options(future.globals.maxSize = max_size_bytes)
  }

  strategy <- tolower(trimws(as.character(strategy)[1]))
  valid_strategies <- c("auto", "sequential", "multisession", "multicore")
  if (!(strategy %in% valid_strategies)) {
    stop(
      "Invalid future strategy '", strategy,
      "'. Valid values: ", paste(valid_strategies, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(workers)) {
    workers <- suppressWarnings(as.integer(workers)[1])
    if (is.na(workers) || workers < 1L) {
      stop("'workers' must be a positive integer.", call. = FALSE)
    }
  }

  if (!requireNamespace("future", quietly = TRUE)) {
    if (strategy != "sequential") {
      stop("Package 'future' is required for parallel Seurat execution.", call. = FALSE)
    }
    return(invisible(list(strategy = "sequential", workers = 1L, max_size_gb = as.numeric(max_size_gb))))
  }

  resolved_strategy <- strategy
  if (resolved_strategy == "auto") {
    supports_multicore <- FALSE
    if (.Platform$OS.type != "windows") {
      supports_multicore <- tryCatch(
        isTRUE(future::supportsMulticore()),
        error = function(e) FALSE
      )
    }
    in_rstudio <- nzchar(Sys.getenv("RSTUDIO")) ||
      nzchar(Sys.getenv("RSTUDIO_SESSION_PORT")) ||
      identical(Sys.getenv("TERM_PROGRAM"), "RStudio")
    resolved_strategy <- if (supports_multicore && !in_rstudio) "multicore" else "multisession"
  }

  if (resolved_strategy == "sequential") {
    future::plan(future::sequential)
    return(invisible(list(strategy = resolved_strategy, workers = 1L, max_size_gb = as.numeric(max_size_gb))))
  }

  if (is.null(workers)) {
    detected_workers <- suppressWarnings(as.integer(tryCatch(future::availableCores(), error = function(e) 2L))[1])
    if (is.na(detected_workers) || detected_workers < 1L) {
      detected_workers <- 2L
    }
    workers <- max(1L, detected_workers - 1L)
  }

  if (resolved_strategy == "multisession") {
    future::plan(future::multisession, workers = workers)
  } else if (resolved_strategy == "multicore") {
    future::plan(future::multicore, workers = workers)
  }

  invisible(list(strategy = resolved_strategy, workers = workers, max_size_gb = as.numeric(max_size_gb)))
}

maybe_join_layers <- function(obj, assay = "RNA") {
  if (!("JoinLayers" %in% getNamespaceExports("Seurat"))) return(obj)
  if (!(assay %in% names(obj@assays))) return(obj)
  tryCatch(
    Seurat::JoinLayers(obj, assay = assay),
    error = function(e) obj
  )
}

get_assay_data_slot <- function(obj, assay = "RNA", slot_name = "data") {
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, slot = slot_name),
    error = function(e1) {
      tryCatch(
        Seurat::GetAssayData(obj, assay = assay, layer = slot_name),
        error = function(e2) NULL
      )
    }
  )
}

get_assay_matrix <- function(obj, assay = "RNA", slot_name = "counts") {
  get_assay_data_slot(obj, assay = assay, slot_name = slot_name)
}

resolve_lfc_col <- function(df) {
  candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  for (x in candidates) {
    if (x %in% colnames(df)) return(x)
  }
  stop("Cannot find logFC column in DEG result.", call. = FALSE)
}

safe_jaccard <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  a <- a[!is.na(a) & a != ""]
  b <- b[!is.na(b) & b != ""]
  denom <- length(union(a, b))
  if (denom == 0) return(NA_real_)
  length(intersect(a, b)) / denom
}

write_table_csv <- function(df, file_path) {
  readr::write_csv(as.data.frame(df, stringsAsFactors = FALSE), file_path)
}

clean_gene_symbols <- function(genes) {
  g <- as.character(genes)
  g <- trimws(g)
  g <- sub("^GRCh[0-9]+[-_]", "", g, ignore.case = TRUE)
  g <- sub("^GRCm39[-_]", "", g, ignore.case = TRUE)
  g <- sub("^hg38[-_]", "", g, ignore.case = TRUE)
  g <- sub("\\.[0-9]+$", "", g)
  g[g == ""] <- NA_character_
  g
}

prepare_obj_for_markers <- function(obj, assay = "RNA") {
  if (!(assay %in% names(obj@assays))) {
    stop("Assay not found in Seurat object: ", assay, call. = FALSE)
  }
  DefaultAssay(obj) <- assay
  obj <- maybe_join_layers(obj, assay = assay)
  assay_data <- get_assay_data_slot(obj, assay = assay, slot_name = "data")
  if (is.null(assay_data) || nrow(assay_data) == 0 || ncol(assay_data) == 0) {
    message(assay, " data slot is empty. Running NormalizeData.")
    obj <- Seurat::NormalizeData(obj, assay = assay, verbose = FALSE)
  }
  obj
}

run_seurat_cluster_markers <- function(obj, cluster_id, assay, min_pct, logfc_threshold) {
  message("Running DEG for cluster ", cluster_id, " vs rest.")
  de <- Seurat::FindMarkers(
    object = obj,
    ident.1 = cluster_id,
    assay = assay,
    slot = "data",
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    verbose = FALSE
  )
  de$gene <- rownames(de)
  de$gene_symbol <- clean_gene_symbols(de$gene)
  lfc_col <- resolve_lfc_col(de)
  de <- de[order(de$p_val_adj, -abs(de[[lfc_col]])), , drop = FALSE]
  rownames(de) <- NULL
  de$cluster <- as.character(cluster_id)
  front_cols <- c("cluster", "gene", "gene_symbol")
  de[, c(front_cols, setdiff(colnames(de), front_cols)), drop = FALSE]
}

resolve_deg_workers <- function(requested_workers, n_tasks) {
  requested_workers <- suppressWarnings(as.integer(requested_workers))
  if (is.na(requested_workers) || requested_workers < 1L) requested_workers <- 1L
  max(1L, min(requested_workers, n_tasks))
}

run_seurat_cluster_markers_parallel <- function(
  obj,
  cluster_levels,
  assay,
  min_pct,
  logfc_threshold,
  workers
) {
  workers <- resolve_deg_workers(workers, length(cluster_levels))
  if (.Platform$OS.type != "unix" && workers > 1L) {
    warning("Forked DEG parallelization is only available on Unix-like systems. Falling back to 1 worker.")
    workers <- 1L
  }

  message(
    "Running DEG across ",
    length(cluster_levels),
    " clusters with Seurat::FindMarkers using ",
    workers,
    " forked worker(s)."
  )

  run_one <- function(cluster_id) {
    tryCatch(
      list(
        cluster = as.character(cluster_id),
        markers = run_seurat_cluster_markers(
          obj = obj,
          cluster_id = cluster_id,
          assay = assay,
          min_pct = min_pct,
          logfc_threshold = logfc_threshold
        ),
        error = NULL
      ),
      error = function(e) {
        list(
          cluster = as.character(cluster_id),
          markers = NULL,
          error = conditionMessage(e)
        )
      }
    )
  }

  marker_results <- if (workers > 1L) {
    parallel::mclapply(
      X = cluster_levels,
      FUN = run_one,
      mc.cores = workers,
      mc.preschedule = FALSE
    )
  } else {
    lapply(cluster_levels, run_one)
  }

  names(marker_results) <- cluster_levels
  attr(marker_results, "workers") <- workers
  marker_results
}

default_deg_similarity_feature_specs <- function() {
  list(
    list(
      name = "all_ranked",
      type = "all",
      description = "All genes retained with signed avg_log2FC."
    ),
    list(
      name = "lenient",
      type = "filtered",
      padj_max = 0.05,
      abs_log2fc_min = 0.25,
      abs_delta_pct_min = 0.05,
      description = "padj < 0.05, abs(avg_log2FC) >= 0.25, abs(pct.1 - pct.2) >= 0.05."
    ),
    list(
      name = "moderate",
      type = "filtered",
      padj_max = 0.01,
      abs_log2fc_min = 0.50,
      abs_delta_pct_min = 0.10,
      description = "padj < 0.01, abs(avg_log2FC) >= 0.50, abs(pct.1 - pct.2) >= 0.10."
    ),
    list(
      name = "strict",
      type = "filtered",
      padj_max = 0.001,
      abs_log2fc_min = 1.00,
      abs_delta_pct_min = 0.20,
      description = "padj < 0.001, abs(avg_log2FC) >= 1.00, abs(pct.1 - pct.2) >= 0.20."
    ),
    list(
      name = "top50_each",
      type = "top_n",
      padj_max = 0.05,
      abs_log2fc_min = 0.25,
      abs_delta_pct_min = 0.05,
      top_n_up = 50L,
      top_n_down = 50L,
      description = "After lenient filtering, keep top 50 up and top 50 down genes by avg_log2FC."
    ),
    list(
      name = "top30_each",
      type = "top_n",
      padj_max = 0.05,
      abs_log2fc_min = 0.25,
      abs_delta_pct_min = 0.05,
      top_n_up = 30L,
      top_n_down = 30L,
      description = "After lenient filtering, keep top 30 up and top 30 down genes by avg_log2FC."
    ),
    list(
      name = "top20_each",
      type = "top_n",
      padj_max = 0.05,
      abs_log2fc_min = 0.25,
      abs_delta_pct_min = 0.05,
      top_n_up = 20L,
      top_n_down = 20L,
      description = "After lenient filtering, keep top 20 up and top 20 down genes by avg_log2FC."
    ),
    list(
      name = "top10_each",
      type = "top_n",
      padj_max = 0.05,
      abs_log2fc_min = 0.25,
      abs_delta_pct_min = 0.05,
      top_n_up = 10L,
      top_n_down = 10L,
      description = "After lenient filtering, keep top 10 up and top 10 down genes by avg_log2FC."
    )
  )
}

normalize_gene_key <- function(x) {
  x <- clean_gene_symbols(x)
  x <- toupper(x)
  x[x == ""] <- NA_character_
  x
}

resolve_gene_label <- function(df) {
  gene_symbol <- if ("gene_symbol" %in% colnames(df)) as.character(df$gene_symbol) else rep(NA_character_, nrow(df))
  gene <- if ("gene" %in% colnames(df)) as.character(df$gene) else rep(NA_character_, nrow(df))
  out <- ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene)
  out[is.na(out) | out == ""] <- gene[is.na(out) | out == ""]
  out
}

collapse_deg_table <- function(df) {
  df$gene_label <- resolve_gene_label(df)
  df$gene_key <- normalize_gene_key(df$gene_label)
  lfc_col <- resolve_lfc_col(df)
  df$avg_log2FC <- as.numeric(df[[lfc_col]])
  df$delta_pct <- as.numeric(df$pct.1) - as.numeric(df$pct.2)
  df$abs_log2fc <- abs(df$avg_log2FC)
  df$abs_delta_pct <- abs(df$delta_pct)
  df$p_val_adj_num <- as.numeric(df$p_val_adj)

  df <- dplyr::filter(df, !is.na(gene_key), !is.na(avg_log2FC), !is.na(p_val_adj_num))
  df <- dplyr::arrange(df, dplyr::desc(abs_log2fc), p_val_adj_num)
  df <- dplyr::group_by(df, gene_key)
  df <- dplyr::slice(df, 1)
  dplyr::ungroup(df)
}

read_cluster_deg <- function(cluster_dir, deg_subdir = NULL, marker_file = NULL) {
  cluster_name <- basename(cluster_dir)
  cluster_id <- sub("^cluster_", "", cluster_name)
  if (!is.null(marker_file)) {
    deg_path <- file.path(cluster_dir, marker_file)
  } else if (!is.null(deg_subdir) && nzchar(deg_subdir)) {
    deg_path <- file.path(cluster_dir, deg_subdir, paste0(cluster_name, "_vs_rest_markers.csv"))
  } else {
    deg_path <- file.path(cluster_dir, paste0(cluster_name, "_vs_rest_markers.csv"))
  }

  if (!file.exists(deg_path)) {
    stop("Missing DEG file: ", deg_path, call. = FALSE)
  }

  df <- readr::read_csv(deg_path, show_col_types = FALSE)
  required_cols <- c("pct.1", "pct.2", "p_val_adj")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required DEG columns in ", deg_path, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  resolve_lfc_col(df)

  df <- collapse_deg_table(df)
  df$cluster <- cluster_id
  df
}

filter_feature_table <- function(df, spec) {
  if (spec$type == "all") {
    keep <- rep(TRUE, nrow(df))
  } else {
    keep <- df$p_val_adj_num < spec$padj_max &
      df$abs_log2fc >= spec$abs_log2fc_min &
      df$abs_delta_pct >= spec$abs_delta_pct_min
  }

  filtered <- df[keep, , drop = FALSE]

  if (spec$type == "top_n") {
    up_tbl <- dplyr::filter(filtered, avg_log2FC > 0)
    up_tbl <- dplyr::arrange(up_tbl, dplyr::desc(avg_log2FC), p_val_adj_num)
    up_tbl <- dplyr::slice_head(up_tbl, n = spec$top_n_up)
    down_tbl <- dplyr::filter(filtered, avg_log2FC < 0)
    down_tbl <- dplyr::arrange(down_tbl, avg_log2FC, p_val_adj_num)
    down_tbl <- dplyr::slice_head(down_tbl, n = spec$top_n_down)
    filtered <- dplyr::bind_rows(up_tbl, down_tbl)
    filtered <- dplyr::arrange(filtered, dplyr::desc(abs_log2fc), p_val_adj_num)
    filtered <- dplyr::distinct(filtered, gene_key, .keep_all = TRUE)
  }

  filtered
}

make_feature_object <- function(df, spec) {
  feature_tbl <- filter_feature_table(df, spec)
  up_tbl <- dplyr::filter(feature_tbl, avg_log2FC > 0)
  down_tbl <- dplyr::filter(feature_tbl, avg_log2FC < 0)

  signed_vec <- feature_tbl$avg_log2FC
  names(signed_vec) <- feature_tbl$gene_key

  list(
    table = feature_tbl,
    up = unique(up_tbl$gene_key),
    down = unique(down_tbl$gene_key),
    signed_vec = signed_vec
  )
}

safe_directional_jaccard <- function(f1, f2) {
  same_up <- safe_jaccard(f1$up, f2$up)
  same_down <- safe_jaccard(f1$down, f2$down)
  cross_up_down <- safe_jaccard(f1$up, f2$down)
  cross_down_up <- safe_jaccard(f1$down, f2$up)

  if (all(is.na(c(same_up, same_down, cross_up_down, cross_down_up)))) {
    return(NA_real_)
  }

  0.5 * dplyr::coalesce(same_up, 0) +
    0.5 * dplyr::coalesce(same_down, 0) -
    0.5 * dplyr::coalesce(cross_up_down, 0) -
    0.5 * dplyr::coalesce(cross_down_up, 0)
}

make_union_vectors <- function(v1, v2) {
  genes <- union(names(v1), names(v2))
  if (length(genes) == 0) {
    return(list(x = numeric(0), y = numeric(0)))
  }
  x <- stats::setNames(rep(0, length(genes)), genes)
  y <- stats::setNames(rep(0, length(genes)), genes)
  x[names(v1)] <- v1
  y[names(v2)] <- v2
  list(x = unname(x), y = unname(y))
}

safe_cosine <- function(v1, v2) {
  vv <- make_union_vectors(v1, v2)
  x <- vv$x
  y <- vv$y
  if (length(x) == 0) return(NA_real_)
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(x * y) / denom
}

safe_spearman <- function(v1, v2) {
  vv <- make_union_vectors(v1, v2)
  x <- vv$x
  y <- vv$y
  if (length(x) < 3) return(NA_real_)
  if (length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
}

compute_pair_metrics <- function(cluster_ids, feature_map, threshold_name) {
  out <- list()
  idx <- 1L

  for (i in seq_len(length(cluster_ids))) {
    for (j in seq(i, length(cluster_ids))) {
      c1 <- cluster_ids[i]
      c2 <- cluster_ids[j]

      if (c1 == c2) {
        directional_jaccard <- 1
        signed_lfc_cosine <- 1
        signed_lfc_spearman <- 1
      } else {
        f1 <- feature_map[[c1]]
        f2 <- feature_map[[c2]]
        directional_jaccard <- safe_directional_jaccard(f1, f2)
        signed_lfc_cosine <- safe_cosine(f1$signed_vec, f2$signed_vec)
        signed_lfc_spearman <- safe_spearman(f1$signed_vec, f2$signed_vec)
      }

      out[[idx]] <- data.frame(
        threshold = threshold_name,
        cluster_1 = c1,
        cluster_2 = c2,
        directional_jaccard = directional_jaccard,
        signed_lfc_cosine = signed_lfc_cosine,
        signed_lfc_spearman = signed_lfc_spearman,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  dplyr::bind_rows(out)
}

pair_to_matrix <- function(pair_df, cluster_ids, metric_name) {
  mat <- matrix(
    NA_real_,
    nrow = length(cluster_ids),
    ncol = length(cluster_ids),
    dimnames = list(cluster_ids, cluster_ids)
  )
  for (i in seq_len(nrow(pair_df))) {
    c1 <- as.character(pair_df$cluster_1[i])
    c2 <- as.character(pair_df$cluster_2[i])
    value <- as.numeric(pair_df[[metric_name]][i])
    mat[c1, c2] <- value
    mat[c2, c1] <- value
  }
  diag(mat) <- 1
  mat
}

write_matrix_csv <- function(mat, file_path, row_id = "cluster") {
  df <- as.data.frame(mat, check.names = FALSE)
  df <- tibble::rownames_to_column(df, row_id)
  readr::write_csv(df, file_path)
}

plot_similarity_heatmap <- function(mat, title, file_pdf, na_to_zero = FALSE, width = 8, height = 7) {
  mat_plot <- mat
  if (isTRUE(na_to_zero)) {
    mat_plot[is.na(mat_plot)] <- 0
  }

  grDevices::pdf(file_pdf, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  heatmap_args <- list(
    mat = mat_plot,
    color = grDevices::colorRampPalette(c("#1f4e79", "white", "#b03a2e"))(101),
    breaks = seq(-1, 1, length.out = 102),
    border_color = NA,
    main = title
  )
  if (!isTRUE(na_to_zero)) {
    heatmap_args$na_col <- "grey90"
  }
  do.call(pheatmap::pheatmap, heatmap_args)
  invisible(mat_plot)
}

manual_merge_rules_from_map <- function(merge_map) {
  if (is.null(merge_map) || length(merge_map) == 0) {
    stop("`merge_map` must be a non-empty named list.", call. = FALSE)
  }
  merge_names <- names(merge_map)
  if (is.null(merge_names) || any(is.na(merge_names) | merge_names == "")) {
    stop("`merge_map` must be a named list.", call. = FALSE)
  }

  source_clusters <- as.character(unlist(merge_map, use.names = FALSE))
  source_clusters <- source_clusters[!is.na(source_clusters) & source_clusters != ""]
  if (length(source_clusters) != length(unique(source_clusters))) {
    dup <- unique(source_clusters[duplicated(source_clusters)])
    stop("Duplicate source cluster(s) in merge_map: ", paste(dup, collapse = ", "), call. = FALSE)
  }

  target_clusters <- rep(as.character(merge_names), lengths(merge_map))
  stats::setNames(target_clusters, source_clusters)
}

build_manual_merge_labels <- function(cluster_vec, merge_map, cluster_order = NULL) {
  cluster_vec <- as.character(cluster_vec)
  rules <- manual_merge_rules_from_map(merge_map)

  out <- cluster_vec
  idx <- cluster_vec %in% names(rules)
  out[idx] <- unname(rules[cluster_vec[idx]])

  original_order <- if (!is.null(cluster_order)) as.character(cluster_order) else sort_maybe_numeric(cluster_vec)
  merged_level_order <- unique(vapply(
    original_order,
    function(cl) {
      if (cl %in% names(rules)) unname(rules[[cl]]) else cl
    },
    character(1)
  ))
  merged_level_order <- c(merged_level_order, setdiff(unique(out), merged_level_order))

  factor(out, levels = merged_level_order)
}

make_manual_merge_mapping_table <- function(original_clusters, merge_map) {
  original_clusters <- sort_maybe_numeric(original_clusters)
  rules <- manual_merge_rules_from_map(merge_map)
  merged_group <- ifelse(original_clusters %in% names(rules), unname(rules[original_clusters]), original_clusters)

  data.frame(
    original_cluster = original_clusters,
    merged_group = merged_group,
    merge_status = ifelse(
      original_clusters %in% names(rules),
      ifelse(original_clusters == merged_group, "anchor", "absorbed"),
      "retained"
    ),
    stringsAsFactors = FALSE
  )
}

make_internal_pair_df_from_merge_map <- function(merge_map) {
  pair_list <- lapply(names(merge_map), function(merge_name) {
    members <- as.character(merge_map[[merge_name]])
    if (length(members) < 2) return(NULL)
    cmb <- utils::combn(members, 2)
    data.frame(
      merged_group = merge_name,
      group_1 = cmb[1, ],
      group_2 = cmb[2, ],
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(pair_list)
  if (is.null(out) || nrow(out) == 0) {
    out <- data.frame(
      merged_group = character(0),
      group_1 = character(0),
      group_2 = character(0),
      stringsAsFactors = FALSE
    )
  }
  out
}

run_markers_vs_rest_default <- function(
  obj,
  ident_col,
  groups,
  out_dir,
  min_pct = 0.10,
  logfc_threshold = 0.25,
  assay = "RNA",
  top_n = 20L
) {
  .ensure_dir(out_dir)
  if (length(groups) == 0) {
    stop("`groups` must contain at least one cluster/group.", call. = FALSE)
  }
  if (!(ident_col %in% colnames(obj@meta.data))) {
    stop("Metadata column is missing: ", ident_col, call. = FALSE)
  }

  obj <- Seurat::SetIdent(obj, value = obj@meta.data[[ident_col]])
  all_markers <- list()
  top_markers <- list()
  marker_tables <- list()

  for (grp in as.character(groups)) {
    marker_file <- file.path(out_dir, paste0("markers_", grp, "_vs_rest.csv"))
    top_file <- file.path(out_dir, paste0("top_positive_markers_", grp, "_vs_rest.csv"))

    if (file.exists(marker_file)) {
      message("[markers vs rest] Reusing existing DEG file: ", marker_file)
      de <- readr::read_csv(marker_file, show_col_types = FALSE)
      de <- as.data.frame(de, stringsAsFactors = FALSE)
    } else {
      message("[markers vs rest] ", ident_col, " = ", grp)
      de <- Seurat::FindMarkers(
        object = obj,
        ident.1 = grp,
        assay = assay,
        slot = "data",
        min.pct = min_pct,
        logfc.threshold = logfc_threshold,
        verbose = FALSE
      )
      de <- tibble::rownames_to_column(as.data.frame(de, stringsAsFactors = FALSE), "gene")
      de$gene_symbol <- clean_gene_symbols(de$gene)
    }
    if (!("gene" %in% colnames(de))) {
      de$gene <- rownames(de)
    }
    if (!("gene_symbol" %in% colnames(de))) {
      de$gene_symbol <- clean_gene_symbols(de$gene)
    }
    lfc_col <- resolve_lfc_col(de)
    de$group <- grp
    de$comparison <- paste0(grp, "_vs_rest")
    de <- de[, c("group", "comparison", "gene", "gene_symbol", setdiff(colnames(de), c("group", "comparison", "gene", "gene_symbol"))), drop = FALSE]

    if (!file.exists(marker_file)) {
      write_table_csv(de, marker_file)
    }

    if (file.exists(top_file)) {
      top_pos <- readr::read_csv(top_file, show_col_types = FALSE)
      top_pos <- as.data.frame(top_pos, stringsAsFactors = FALSE)
    } else {
      top_pos <- de |>
        dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05, .data[[lfc_col]] > 0) |>
        dplyr::arrange(dplyr::desc(.data[[lfc_col]]), p_val_adj, gene) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::mutate(rank_within_group = dplyr::row_number())

      write_table_csv(top_pos, top_file)
    }

    all_markers[[grp]] <- de
    top_markers[[grp]] <- top_pos
    marker_tables[[grp]] <- list(full = de, top_positive = top_pos, lfc_col = lfc_col)
  }

  list(
    full = dplyr::bind_rows(all_markers),
    top_positive = dplyr::bind_rows(top_markers),
    by_group = marker_tables
  )
}

run_pairwise_markers_default <- function(
  obj,
  ident_col,
  pair_df,
  out_dir,
  min_pct = 0.10,
  logfc_threshold = 0.25,
  assay = "RNA"
) {
  .ensure_dir(out_dir)
  required_cols <- c("group_1", "group_2")
  missing_cols <- setdiff(required_cols, colnames(pair_df))
  if (length(missing_cols) > 0) {
    stop("pair_df is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (nrow(pair_df) == 0) {
    return(invisible(list(summary = data.frame(), by_pair = list())))
  }

  obj <- Seurat::SetIdent(obj, value = obj@meta.data[[ident_col]])
  summary_rows <- list()
  pair_tables <- list()

  for (i in seq_len(nrow(pair_df))) {
    g1 <- as.character(pair_df$group_1[i])
    g2 <- as.character(pair_df$group_2[i])
    pair_name <- paste0(g1, "_vs_", g2)
    marker_file <- file.path(out_dir, paste0("markers_", pair_name, ".csv"))
    if (file.exists(marker_file)) {
      message("[pairwise markers] Reusing existing DEG file: ", marker_file)
      de <- readr::read_csv(marker_file, show_col_types = FALSE)
      de <- as.data.frame(de, stringsAsFactors = FALSE)
    } else {
      message("[pairwise markers] ", ident_col, " : ", pair_name)
      de <- Seurat::FindMarkers(
        object = obj,
        ident.1 = g1,
        ident.2 = g2,
        assay = assay,
        slot = "data",
        min.pct = min_pct,
        logfc.threshold = logfc_threshold,
        verbose = FALSE
      )
      de <- tibble::rownames_to_column(as.data.frame(de, stringsAsFactors = FALSE), "gene")
      de$gene_symbol <- clean_gene_symbols(de$gene)
    }
    if (!("gene" %in% colnames(de))) {
      de$gene <- rownames(de)
    }
    if (!("gene_symbol" %in% colnames(de))) {
      de$gene_symbol <- clean_gene_symbols(de$gene)
    }
    lfc_col <- resolve_lfc_col(de)
    if (!file.exists(marker_file)) {
      write_table_csv(de, marker_file)
    }

    n_sig <- sum(!is.na(de$p_val_adj) & de$p_val_adj < 0.05)
    n_sig_abs025 <- sum(!is.na(de$p_val_adj) & de$p_val_adj < 0.05 & abs(de[[lfc_col]]) >= 0.25)
    top_up_gene <- de |>
      dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05, .data[[lfc_col]] > 0) |>
      dplyr::arrange(dplyr::desc(.data[[lfc_col]]), p_val_adj, gene) |>
      dplyr::slice_head(n = 1) |>
      dplyr::pull(gene)
    top_down_gene <- de |>
      dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05, .data[[lfc_col]] < 0) |>
      dplyr::arrange(.data[[lfc_col]], p_val_adj, gene) |>
      dplyr::slice_head(n = 1) |>
      dplyr::pull(gene)

    summary_rows[[pair_name]] <- data.frame(
      ident_col = ident_col,
      group_1 = g1,
      group_2 = g2,
      n_sig = n_sig,
      n_sig_abs_log2fc_0.25 = n_sig_abs025,
      median_abs_log2fc = if (nrow(de) > 0) median(abs(de[[lfc_col]]), na.rm = TRUE) else NA_real_,
      top_up_gene = ifelse(length(top_up_gene) == 0, NA_character_, top_up_gene[1]),
      top_down_gene = ifelse(length(top_down_gene) == 0, NA_character_, top_down_gene[1]),
      stringsAsFactors = FALSE
    )
    pair_tables[[pair_name]] <- de
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  write_table_csv(summary_df, file.path(out_dir, "pairwise_marker_summary.csv"))
  invisible(list(summary = summary_df, by_pair = pair_tables))
}

extract_top_positive_genes <- function(marker_tbl, n_top = 50L) {
  lfc_col <- resolve_lfc_col(marker_tbl)
  marker_tbl |>
    dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05, .data[[lfc_col]] > 0) |>
    dplyr::arrange(dplyr::desc(.data[[lfc_col]]), p_val_adj, gene) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::pull(gene) |>
    unique()
}

build_marker_overlap_summary <- function(original_markers, merged_markers, merge_map, n_top = 50L) {
  out <- list()
  idx <- 1L

  for (merge_name in names(merge_map)) {
    merged_tbl <- merged_markers$by_group[[merge_name]]$full
    merged_genes <- extract_top_positive_genes(merged_tbl, n_top = n_top)

    source_clusters <- as.character(merge_map[[merge_name]])
    source_gene_sets <- lapply(source_clusters, function(cl) {
      extract_top_positive_genes(original_markers$by_group[[cl]]$full, n_top = n_top)
    })
    names(source_gene_sets) <- source_clusters

    union_genes <- unique(unlist(source_gene_sets, use.names = FALSE))
    intersect_genes <- if (length(source_gene_sets) > 0) Reduce(intersect, source_gene_sets) else character(0)

    compare_sets <- c(source_gene_sets, list(union_of_sources = union_genes, intersect_of_sources = intersect_genes))
    for (ref_name in names(compare_sets)) {
      ref_genes <- compare_sets[[ref_name]]
      out[[idx]] <- data.frame(
        merged_group = merge_name,
        reference_set = ref_name,
        merged_top_n = length(merged_genes),
        reference_top_n = length(ref_genes),
        overlap_n = length(intersect(merged_genes, ref_genes)),
        jaccard = safe_jaccard(merged_genes, ref_genes),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  dplyr::bind_rows(out)
}

make_manual_merge_feature_object <- function(marker_tbl, spec) {
  spec_use <- spec
  if (is.null(spec_use$type)) spec_use$type <- "top_n"
  feature_tbl <- collapse_deg_table(marker_tbl)
  feature_tbl <- filter_feature_table(feature_tbl, spec_use)
  signed_vec <- feature_tbl$avg_log2FC
  names(signed_vec) <- feature_tbl$gene_key

  list(
    table = feature_tbl,
    signed_vec = signed_vec
  )
}

build_manual_merge_cosine_summary <- function(
  original_markers,
  merged_markers,
  original_groups,
  merged_groups,
  merge_map,
  spec
) {
  original_groups <- as.character(original_groups)
  merged_groups <- as.character(merged_groups)
  original_feature_map <- stats::setNames(
    lapply(original_groups, function(grp) make_manual_merge_feature_object(original_markers$by_group[[grp]]$full, spec)),
    original_groups
  )
  merged_feature_map <- stats::setNames(
    lapply(merged_groups, function(grp) make_manual_merge_feature_object(merged_markers$by_group[[grp]]$full, spec)),
    merged_groups
  )

  rows <- list()
  idx <- 1L
  for (merge_name in merged_groups) {
    for (original_name in original_groups) {
      merged_feature <- merged_feature_map[[merge_name]]
      original_feature <- original_feature_map[[original_name]]

      rows[[idx]] <- data.frame(
        merged_group = merge_name,
        original_cluster = original_name,
        is_source_cluster = original_name %in% as.character(merge_map[[merge_name]]),
        merged_n_features = nrow(merged_feature$table),
        original_n_features = nrow(original_feature$table),
        signed_lfc_cosine = safe_cosine(merged_feature$signed_vec, original_feature$signed_vec),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  dplyr::bind_rows(rows)
}

cosine_summary_to_matrix <- function(cosine_df, merged_groups, original_groups) {
  mat <- matrix(
    NA_real_,
    nrow = length(merged_groups),
    ncol = length(original_groups),
    dimnames = list(as.character(merged_groups), as.character(original_groups))
  )

  for (i in seq_len(nrow(cosine_df))) {
    mat[
      as.character(cosine_df$merged_group[i]),
      as.character(cosine_df$original_cluster[i])
    ] <- as.numeric(cosine_df$signed_lfc_cosine[i])
  }

  mat
}

plot_cosine_heatmap <- function(mat, title_text, file_stub) {
  plot_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(plot_df) <- c("merged_group", "original_cluster", "signed_lfc_cosine")
  plot_df$merged_group <- factor(plot_df$merged_group, levels = rev(rownames(mat)))
  plot_df$original_cluster <- factor(plot_df$original_cluster, levels = colnames(mat))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = original_cluster, y = merged_group, fill = signed_lfc_cosine)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(signed_lfc_cosine), "NA", sprintf("%.2f", signed_lfc_cosine))),
      size = 3
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#1f4e79",
      mid = "white",
      high = "#b03a2e",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey90"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = title_text,
      x = "Original cluster",
      y = "Merged group",
      fill = "Cosine"
    )

  save_plot_pdf_png(
    p,
    file_stub,
    width = max(8, 0.8 * ncol(mat) + 3),
    height = max(4, 0.8 * nrow(mat) + 2)
  )
  invisible(p)
}

extract_features_for_dotplot <- function(top_marker_df, group_col = "group", n_per_group = 5L) {
  top_marker_df |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::arrange(rank_within_group, .by_group = TRUE) |>
    dplyr::slice_head(n = n_per_group) |>
    dplyr::ungroup() |>
    dplyr::pull(gene) |>
    unique()
}

plot_dotplot_safe <- function(obj, features, group.by, title_text, file_stub, assay = "RNA") {
  features <- unique(features[!is.na(features) & features != ""])
  if (length(features) == 0) {
    message("Skipping dotplot for ", group.by, " because no features were selected.")
    return(invisible(NULL))
  }
  p <- Seurat::DotPlot(obj, features = features, group.by = group.by, assay = assay) +
    Seurat::RotatedAxis() +
    ggplot2::labs(title = title_text)
  save_plot_pdf_png(p, file_stub, width = max(8, 0.32 * length(features) + 3), height = 6)
  invisible(p)
}

normalize_ora_gene_key <- function(x) {
  normalize_gene_key(x)
}

resolve_deg_gene_label <- function(df) {
  gene_symbol <- if ("gene_symbol" %in% colnames(df)) as.character(df$gene_symbol) else rep(NA_character_, nrow(df))
  gene <- if ("gene" %in% colnames(df)) as.character(df$gene) else rownames(df)
  out <- ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene)
  out[is.na(out) | out == ""] <- gene[is.na(out) | out == ""]
  out
}

prepare_deg_table_for_ora <- function(df, lfc_col) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$gene_label <- resolve_deg_gene_label(df)
  df$gene_key <- normalize_ora_gene_key(df$gene_label)
  df$gene_symbol <- clean_gene_symbols(df$gene_label)
  missing_gene_symbol <- is.na(df$gene_symbol) | df$gene_symbol == ""
  df$gene_symbol[missing_gene_symbol] <- df$gene_key[missing_gene_symbol]
  df$lfc_value <- as.numeric(df[[lfc_col]])
  df$p_val_adj_num <- as.numeric(df$p_val_adj)
  df$delta_pct <- as.numeric(df$pct.1) - as.numeric(df$pct.2)
  df$abs_logfc <- abs(df$lfc_value)
  df$abs_delta_pct <- abs(df$delta_pct)

  df <- df |>
    dplyr::filter(
      !is.na(gene_key),
      !is.na(gene_symbol),
      !is.na(lfc_value),
      !is.na(p_val_adj_num),
      !is.na(delta_pct)
    ) |>
    dplyr::arrange(dplyr::desc(abs_logfc), p_val_adj_num)

  df |>
    dplyr::group_by(gene_key) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
}

select_deg_top_up_for_ora <- function(
  df,
  lfc_col,
  padj_max,
  abs_logfc_min,
  abs_delta_pct_min,
  top_n
) {
  filtered <- prepare_deg_table_for_ora(df, lfc_col) |>
    dplyr::filter(
      p_val_adj_num < padj_max,
      abs_logfc >= abs_logfc_min,
      abs_delta_pct >= abs_delta_pct_min
    )

  filtered |>
    dplyr::filter(lfc_value > 0) |>
    dplyr::arrange(dplyr::desc(lfc_value), p_val_adj_num) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(direction = "up") |>
    dplyr::arrange(dplyr::desc(abs_logfc), p_val_adj_num) |>
    dplyr::distinct(gene_key, .keep_all = TRUE)
}

beautify_hallmark_name <- function(x) {
  x <- as.character(x)
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

get_hallmark_sets <- function(species = "Homo sapiens") {
  df_h <- tryCatch(
    msigdbr::msigdbr(species = species, collection = "H"),
    error = function(e) msigdbr::msigdbr(species = species, category = "H")
  )
  sets <- split(df_h$gene_symbol, df_h$gs_name)
  lapply(sets, unique)
}

run_ora_hypergeom <- function(query_genes, universe_genes, pathways, min_size = 15, max_size = 500, min_overlap = 3) {
  query <- unique(stats::na.omit(as.character(query_genes)))
  universe <- unique(stats::na.omit(as.character(universe_genes)))
  if (length(query) == 0 || length(universe) == 0) return(data.frame())
  query <- intersect(query, universe)
  if (length(query) == 0) return(data.frame())

  out <- lapply(names(pathways), function(pw_name) {
    pw_genes <- unique(intersect(as.character(pathways[[pw_name]]), universe))
    M <- length(pw_genes)
    if (M < min_size || M > max_size) return(NULL)
    overlap <- intersect(query, pw_genes)
    k <- length(overlap)
    if (k < min_overlap) return(NULL)

    U <- length(universe)
    N <- length(query)
    pval <- stats::phyper(q = k - 1, m = M, n = U - M, k = N, lower.tail = FALSE)
    odds <- suppressWarnings((k / max(1, N - k)) / (M / max(1, U - M)))

    data.frame(
      pathway = pw_name,
      hallmark_label = beautify_hallmark_name(pw_name),
      set_size = M,
      query_size = N,
      overlap = k,
      gene_ratio = sprintf("%d/%d", k, N),
      bg_ratio = sprintf("%d/%d", M, U),
      odds_ratio = odds,
      p_value = pval,
      p_adj = NA_real_,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) return(data.frame())
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[order(out$p_adj, out$p_value, -out$overlap), , drop = FALSE]
  rownames(out) <- NULL
  out
}

safe_scale01 <- function(x) {
  x <- as.numeric(x)
  out <- rep(0, length(x))
  keep <- is.finite(x)
  if (!any(keep)) return(out)

  rng <- range(x[keep], na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    out[keep] <- 1
    return(out)
  }

  out[keep] <- (x[keep] - rng[1]) / diff(rng)
  out
}

arrange_ora_for_annotation <- function(ora_df) {
  if (is.null(ora_df) || nrow(ora_df) == 0) return(ora_df)

  df <- as.data.frame(ora_df, stringsAsFactors = FALSE)
  primary <- if ("annotation_score" %in% colnames(df)) as.numeric(df$annotation_score) else rep(NA_real_, nrow(df))
  secondary <- if ("expression_score_raw" %in% colnames(df)) as.numeric(df$expression_score_raw) else rep(NA_real_, nrow(df))
  tertiary <- if ("detection_score_raw" %in% colnames(df)) as.numeric(df$detection_score_raw) else rep(NA_real_, nrow(df))
  quaternary <- if ("overlap_score_raw" %in% colnames(df)) as.numeric(df$overlap_score_raw) else rep(NA_real_, nrow(df))

  primary[!is.finite(primary)] <- -Inf
  secondary[!is.finite(secondary)] <- -Inf
  tertiary[!is.finite(tertiary)] <- -Inf
  quaternary[!is.finite(quaternary)] <- -Inf

  df[order(-primary, -secondary, -tertiary, -quaternary, df$pathway), , drop = FALSE]
}

score_ora_for_annotation <- function(
  ora_df,
  deg_df,
  weight_expression = 1 / 3,
  weight_detection = 1 / 3,
  weight_overlap = 1 / 3,
  fdr_cutoff = 0.05
) {
  if (is.null(ora_df) || nrow(ora_df) == 0) return(data.frame())

  ora_df <- as.data.frame(ora_df, stringsAsFactors = FALSE)
  ora_df <- ora_df[!is.na(ora_df$p_adj) & ora_df$p_adj < fdr_cutoff, , drop = FALSE]
  if (nrow(ora_df) == 0) return(data.frame())
  deg_df <- as.data.frame(deg_df, stringsAsFactors = FALSE)

  if (!("gene_key" %in% colnames(deg_df))) {
    deg_df$gene_key <- normalize_ora_gene_key(deg_df$gene_symbol)
  }
  if (!("lfc_value" %in% colnames(deg_df))) {
    lfc_col <- resolve_lfc_col(deg_df)
    deg_df$lfc_value <- as.numeric(deg_df[[lfc_col]])
  }
  if (!("abs_logfc" %in% colnames(deg_df))) {
    deg_df$abs_logfc <- abs(as.numeric(deg_df$lfc_value))
  }
  if (!("delta_pct" %in% colnames(deg_df))) {
    deg_df$delta_pct <- as.numeric(deg_df$pct.1) - as.numeric(deg_df$pct.2)
  }
  if (!("abs_delta_pct" %in% colnames(deg_df))) {
    deg_df$abs_delta_pct <- abs(as.numeric(deg_df$delta_pct))
  }
  if (!("direction" %in% colnames(deg_df))) {
    deg_df$direction <- ifelse(deg_df$lfc_value > 0, "up", ifelse(deg_df$lfc_value < 0, "down", "flat"))
  }

  overlap_stats <- lapply(seq_len(nrow(ora_df)), function(i) {
    overlap_genes <- ora_df$overlap_genes[i]
    overlap_symbols <- if (is.na(overlap_genes) || overlap_genes == "") {
      character(0)
    } else {
      trimws(unlist(strsplit(overlap_genes, ";", fixed = TRUE)))
    }
    overlap_keys <- normalize_ora_gene_key(overlap_symbols)
    overlap_keys <- unique(overlap_keys[!is.na(overlap_keys)])

    if (length(overlap_keys) == 0) {
      return(data.frame(
        marker_mean_abs_logfc = NA_real_,
        marker_mean_pct1 = NA_real_,
        marker_mean_abs_delta_pct = NA_real_,
        n_overlap_up = 0L,
        n_overlap_down = 0L,
        stringsAsFactors = FALSE
      ))
    }

    overlap_deg <- deg_df[deg_df$gene_key %in% overlap_keys, , drop = FALSE]
    overlap_deg <- overlap_deg[!duplicated(overlap_deg$gene_key), , drop = FALSE]

    data.frame(
      marker_mean_abs_logfc = if (nrow(overlap_deg) > 0) mean(overlap_deg$abs_logfc, na.rm = TRUE) else NA_real_,
      marker_mean_pct1 = if (nrow(overlap_deg) > 0) mean(as.numeric(overlap_deg$pct.1), na.rm = TRUE) else NA_real_,
      marker_mean_abs_delta_pct = if (nrow(overlap_deg) > 0) mean(overlap_deg$abs_delta_pct, na.rm = TRUE) else NA_real_,
      n_overlap_up = sum(overlap_deg$direction == "up", na.rm = TRUE),
      n_overlap_down = sum(overlap_deg$direction == "down", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  overlap_stats_df <- dplyr::bind_rows(overlap_stats)
  ora_df <- cbind(ora_df, overlap_stats_df)

  ora_df$overlap_query_fraction <- ifelse(ora_df$query_size > 0, ora_df$overlap / ora_df$query_size, NA_real_)
  ora_df$overlap_set_fraction <- ifelse(ora_df$set_size > 0, ora_df$overlap / ora_df$set_size, NA_real_)
  overlap_fraction_mat <- cbind(ora_df$overlap_query_fraction, ora_df$overlap_set_fraction)
  ora_df$overlap_score_raw <- rowMeans(overlap_fraction_mat, na.rm = TRUE)
  ora_df$overlap_score_raw[!is.finite(ora_df$overlap_score_raw)] <- NA_real_

  ora_df$expression_score_raw <- ora_df$marker_mean_abs_logfc
  ora_df$detection_score_raw <- ora_df$marker_mean_pct1

  ora_df$expression_score_scaled <- safe_scale01(ora_df$expression_score_raw)
  ora_df$detection_score_scaled <- safe_scale01(ora_df$detection_score_raw)
  ora_df$overlap_score_scaled <- safe_scale01(ora_df$overlap_score_raw)

  total_weight <- weight_expression + weight_detection + weight_overlap
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("Annotation score weights must sum to a positive finite value.", call. = FALSE)
  }

  ora_df$annotation_score <- (
    weight_expression * ora_df$expression_score_scaled +
      weight_detection * ora_df$detection_score_scaled +
      weight_overlap * ora_df$overlap_score_scaled
  ) / total_weight

  ora_df <- arrange_ora_for_annotation(ora_df)
  ora_df$annotation_rank <- seq_len(nrow(ora_df))
  rownames(ora_df) <- NULL
  ora_df
}

plot_ora_top <- function(ora_df, out_pdf, out_png, title, top_n = 15) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    grDevices::pdf(out_pdf, width = 10, height = 5)
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste0(title, "\nNo Hallmark ORA results"))
    grDevices::dev.off()
    grDevices::png(out_png, width = 3000, height = 1500, res = 300)
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste0(title, "\nNo Hallmark ORA results"))
    grDevices::dev.off()
    return(invisible(NULL))
  }

  df <- ora_df
  df <- df[is.finite(df$p_adj) & !is.na(df$p_adj), , drop = FALSE]
  if (nrow(df) == 0) {
    grDevices::pdf(out_pdf, width = 10, height = 5)
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste0(title, "\nNo valid Hallmark ORA rows"))
    grDevices::dev.off()
    grDevices::png(out_png, width = 3000, height = 1500, res = 300)
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste0(title, "\nNo valid Hallmark ORA rows"))
    grDevices::dev.off()
    return(invisible(NULL))
  }

  use_annotation_score <- "annotation_score" %in% colnames(df) && any(is.finite(df$annotation_score))
  if (use_annotation_score) {
    df <- arrange_ora_for_annotation(df)
    df$plot_value <- df$annotation_score
    y_label <- "Integrated annotation score"
  } else {
    df_sig <- df[df$p_adj < 0.05, , drop = FALSE]
    if (nrow(df_sig) > 0) df <- df_sig
    df <- df[order(df$p_adj, df$p_value, -df$overlap, df$pathway), , drop = FALSE]
    df$plot_value <- -log10(pmax(df$p_adj, 1e-300))
    y_label <- "-log10(FDR)"
  }
  df <- head(df, top_n)
  df$hallmark_label <- factor(df$hallmark_label, levels = rev(df$hallmark_label))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = hallmark_label, y = plot_value)) +
    ggplot2::geom_col(fill = "#2c7fb8", width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = y_label) +
    ggplot2::theme_classic(base_size = 11)

  ggplot2::ggsave(out_pdf, p, width = 10, height = max(5, 0.28 * nrow(df)))
  ggplot2::ggsave(out_png, p, width = 10, height = max(5, 0.28 * nrow(df)), dpi = 300)
  invisible(df)
}

summarize_cluster_annotation <- function(cluster_id, n_cells, deg_df, ora_df, top_n = 3L, fdr_cutoff = 0.05) {
  deg_genes_n <- if (is.null(deg_df) || nrow(deg_df) == 0) 0L else nrow(deg_df)
  up_genes_n <- if (is.null(deg_df) || nrow(deg_df) == 0 || !("direction" %in% colnames(deg_df))) {
    0L
  } else {
    sum(deg_df$direction == "up", na.rm = TRUE)
  }
  down_genes_n <- if (is.null(deg_df) || nrow(deg_df) == 0 || !("direction" %in% colnames(deg_df))) {
    0L
  } else {
    sum(deg_df$direction == "down", na.rm = TRUE)
  }

  empty_row <- function(note) {
    data.frame(
      cluster = as.character(cluster_id),
      n_cells = n_cells,
      n_deg_for_ora = deg_genes_n,
      n_up_deg_for_ora = up_genes_n,
      n_down_deg_for_ora = down_genes_n,
      n_significant_hallmarks = 0L,
      annotation_primary = NA_character_,
      annotation_secondary = NA_character_,
      annotation_tertiary = NA_character_,
      annotation_multi = NA_character_,
      top_hallmark_1_annotation_score = NA_real_,
      top_hallmark_2_annotation_score = NA_real_,
      top_hallmark_3_annotation_score = NA_real_,
      top_hallmark_1_fdr = NA_real_,
      top_hallmark_2_fdr = NA_real_,
      top_hallmark_3_fdr = NA_real_,
      note = note,
      stringsAsFactors = FALSE
    )
  }

  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(empty_row("no Hallmark ORA result"))
  }

  ora_sig <- ora_df |>
    dplyr::filter(!is.na(p_adj), p_adj < fdr_cutoff) |>
    dplyr::arrange(p_adj, dplyr::desc(overlap), pathway)

  if (nrow(ora_sig) == 0) {
    return(empty_row(paste0("no Hallmark passes FDR < ", fdr_cutoff, "; no annotation assigned")))
  }

  note <- paste0(
    "annotation based on integrated score combining marker strength, marker detection fraction, and overlap degree after Hallmark FDR < ",
    fdr_cutoff,
    " filtering"
  )

  ora_use <- arrange_ora_for_annotation(ora_df)
  ora_use <- head(ora_use, top_n)

  labels <- ora_use$hallmark_label
  padj_values <- ora_use$p_adj
  annotation_scores <- if ("annotation_score" %in% colnames(ora_use)) ora_use$annotation_score else rep(NA_real_, nrow(ora_use))
  labels <- c(labels, rep(NA_character_, max(0, top_n - length(labels))))
  padj_values <- c(padj_values, rep(NA_real_, max(0, top_n - length(padj_values))))
  annotation_scores <- c(annotation_scores, rep(NA_real_, max(0, top_n - length(annotation_scores))))

  data.frame(
    cluster = as.character(cluster_id),
    n_cells = n_cells,
    n_deg_for_ora = deg_genes_n,
    n_up_deg_for_ora = up_genes_n,
    n_down_deg_for_ora = down_genes_n,
    n_significant_hallmarks = nrow(ora_sig),
    annotation_primary = labels[1],
    annotation_secondary = labels[2],
    annotation_tertiary = labels[3],
    annotation_multi = paste(labels[!is.na(labels)], collapse = "; "),
    top_hallmark_1_annotation_score = annotation_scores[1],
    top_hallmark_2_annotation_score = annotation_scores[2],
    top_hallmark_3_annotation_score = annotation_scores[3],
    top_hallmark_1_fdr = padj_values[1],
    top_hallmark_2_fdr = padj_values[2],
    top_hallmark_3_fdr = padj_values[3],
    note = note,
    stringsAsFactors = FALSE
  )
}
