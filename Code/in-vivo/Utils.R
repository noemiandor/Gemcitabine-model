
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
    
    outputs[[patient]] = list(cn=cn, cells=sample, anno=anno)
    
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