library(tidyverse)

# root="~/Projects/PMO/HighPloidy_DoubleEdgedSword/data/BreastCancerCLs/SUM159/K00_GemcitabineExposure_033023"
# path2repo="~/Repositories/Gemcitabine-model"

root <- "C:/Users/80027908/Desktop/K00_GemcitabineExposure_033023"
path2repo <- "~/Gemcitabine-model"
filesep <- "/"

# Load in tracking results
tracking_files <- paste0(root, filesep, "Tracking_CSVs")
tracking_csvs <- list.files(tracking_files, full.names = T)

# Read in the classifier file
class <- read.csv(paste0(root, filesep, "F6_Object_data.csv"))

# Set the well
well <- "F6"

# Get the well information by splitting the image file paths
class_1 <- class %>% 
  filter(grepl(well, Image.Location)) %>% 
  separate(col = Image.Location, 
           into = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "well_info", NA), 
           sep = "\\\\", 
           remove = FALSE)

# Load in each field
for(i in 1:4) {
  dat <- read.csv(tracking_csvs[grepl(paste0(well, "_", i), tracking_csvs)]) %>% 
    mutate(well_info = paste0(well, "_", i))
  assign(paste0("dat", i), dat)
}

data <- list(dat1 = dat1, dat2 = dat2, dat3 = dat3, dat4 = dat4)

for(frame in names(data)) {
  dat <- get(frame)
  # Rename columns
  ii=match(c("frame","Bounding_Box_Maximum_0","Bounding_Box_Maximum_1"),colnames(dat))
  colnames(dat)[ii]=c("t","x","y")
  # Save data
  assign(frame, dat)
}

# Pull in the data for the well of interest
class_1$time <- sapply(strsplit(class_1$Image.Location,"_"), function(x) x[length(x)])
class_1$time <- gsub(".tif","",class_1$time)

# Change times to frames
map <- 0:(length(unique(class_1$time))-1)
names(map) <- sort(unique(class_1$time))
class_1$t=map[class_1$time]
class_1 <- class_1[order(class_1$t),]

# Check changes
# boxplot(class_1$t~class_1$time, horizontal = T, las=2, ylab="", cex.axis=0.54)

for(frame in names(data)) {
  dat <- get(frame)
  seg <- class_1 %>% 
    filter(well_info == dat$well_info[1]) %>% 
    group_by(t)
  
  # Rename classifier column names to match data
  ii=match(c("XMax","YMax"),colnames(seg))
  colnames(seg)[ii]=c("x","y")
  # Get the intersection of column names -> should only be t, x, and y
  coi=intersect(colnames(seg),colnames(dat))[2:4]
  # Align data by time, matching to the previous "map" data
  ii=sapply(map, function(t) list(seg=which(seg$t==t), dat=which(dat$t==t)), simplify = F)
  ii=ii[sapply(ii, function(x) min(length(x[[1]]),length(x[[2]])))>0]
  d=sapply(ii, function(t) flexclust::dist2(seg[t$seg,coi,drop=F],dat[t$dat,coi,drop=F]))
  seg_matched=sapply(names(d), function(x) seg[ii[[x]]$seg,][apply(d[[x]],2,which.min),],simplify = F)
  merged=sapply(names(d), function(x) cbind(dat[ii[[x]]$dat,], seg_matched[[x]][,c("time","Classifier.Phenotype")]), simplify = F)
  merged=do.call(rbind, merged)
  merged=merged[order(merged$t),];
  # Find the center of the cells and name them x & y
  merged$y=(merged$Bounding_Box_Minimum_1+merged$y)/2
  merged$x=(merged$Bounding_Box_Minimum_0+merged$x)/2
  merged$y=-merged$y
  merged$y=merged$y-min(merged$y)
  # Set classification to a numerical value
  merged$Class=1
  merged$Class[merged$Classifier.Phenotype=="Transitional"]=0
  merged$Class[merged$Classifier.Phenotype=="Dead"]=-1
  # Save dataframe
  assign(paste("merged", frame, sep = "_"), merged)
}

#Create an "is not in" function
'%!in%' <- function(x,y)!('%in%'(x,y))

merged_data <- list(merged_dat1 = merged_dat1, 
                    merged_dat2 = merged_dat2, 
                    merged_dat3 = merged_dat3, 
                    merged_dat4 = merged_dat4)

# Find cells at last time points, then label as parents/dividing.
for(data in names(merged_data)) {
  df <- get(data)
  #Find end point class
  df_1 <- filter(df, trackId != -1)
  cells_list <- unique(df_1$trackId)
  #Recreate parents list
  parents_list <- unique(df$parentTrackId)
  parents_list <- parents_list[!parents_list %in% 0]
  allcells_max <- data.frame()
  for(i in cells_list) {
    p <- filter(df, trackId == i)
    q <- max(p$t)
    r <- filter(p, t == q)
    if(i %in% parents_list) {
      r <- mutate(r, parent = TRUE)
      allcells_max <- rbind(allcells_max, r)
    }
    if(i %!in% parents_list) {
      r <- mutate(r, parent = FALSE)
      allcells_max <- rbind(allcells_max, r)    
    }
  }
  # Visualize class of cells at end points
  allcells_max <- filter(allcells_max, trackId != -1)
  print(ggplot(allcells_max) + 
          geom_bar(aes(Classifier.Phenotype, fill = parent), position = 'dodge') + 
          theme_light() + 
          xlab("class") + 
          ggtitle(paste0("Class of Cells at Last Time Point ", df[1,63])))
  # Get list of all time points labeled with dividing and parent statuses
  # Time consuming step
  allcells <- data.frame()
  for(i in cells_list) {
    p <- filter(df, trackId == i)
    q <- max(p$t)
    p <- cbind(p, lifetime = c(1:nrow(p)))
    r <- filter(p, t == q)
    s <- filter(p, t != q) %>% mutate(dividing = FALSE)
    if(i %in% parents_list) {
      r <- mutate(r, dividing = TRUE)
      r <- rbind(s, r)
      r <- mutate(r, parent = TRUE)
      allcells <- rbind(allcells, r)
    }
    if(i %!in% parents_list) {
      r <- mutate(r, dividing = FALSE)
      r <- rbind(s, r)
      r <- mutate(r, parent = FALSE)
      allcells <- rbind(allcells, r)    
    }
  }
  assign(paste("allcells", data, sep = "_") , allcells)
}

##################### The graphs etc for F6_1 only
library(matlab)
library(ggpubr)

all_dat1_filt <- allcells_merged_dat1 %>%
  filter(t < 41)
cor_dat1 <- data.frame()
for(cell in sort(unique(all_dat1_filt$trackId))) {
  df_1 <- filter(all_dat1_filt, trackId == cell)
  if(max(df_1$lifetime) > 24) {
    df_2 <- df_1 %>%
      mutate(cor = cor(t, Object_Area_0, method = "pearson"))
    cor_dat1 <- rbind(cor_dat1, df_2)
    print(cell)
  }
}

cor_dat1_list <- cor_dat1[order(-cor_dat1$cor),]
cor_dat1_list <- unique(cor_dat1_list$trackId)
targets <- cor_dat1_list[1:10]

cor_dat1_alltime <- allcells_merged_dat1 %>% 
  filter(trackId %in% targets | parentTrackId %in% targets)

cor_dat2 <- data.frame()
for(cell in sort(unique(cor_dat1_alltime$trackId))) {
  df_1 <- filter(cor_dat1_alltime, trackId == cell | parentTrackId == cell)
  if(max(df_1$lifetime) > 24) {
    df_2 <- df_1 %>%
      mutate(cor = cor(t, Object_Area_0, method = "pearson"))
    cor_dat2 <- rbind(cor_dat2, df_2)
    print(cell)
  }
}

well <- "F6_1"
col_values <- c("Dead" = "cyan", "Unstained" = "green", "Alive" = "orangered", "Transitional" = "magenta")

# fi <- list.files(paste0(root, filesep, well),full.names = T)
fi <- list.files(paste0(root, filesep, well, "/HALO Markup/output/overlay_output"),full.names = T)
for(cell in targets) {
  filt <- cor_dat2 %>%
    filter(trackId == cell | parentTrackId == cell)
  for(i in min(filt$t):(max(filt$t)+5)) {
    tiff(paste0("~/corr_", well, "_", cell, "_", i, ".tif"), width = 1468, height = 1100)
    img=bioimagetools::readTIF(fi[i+1],as.is = T)
    img <- EBImage::resize(img, dim(img)[1]/1)
    plot(raster::as.raster(img[,,,1]))
    if(i %in% unique(filt$t)) {
      p <- filt %>% filter(t == i)
      p$Class[p$Classifier.Phenotype=="UnStained"]=2
      points(p$x, p$y+30, pch=(p$Class+2), col="black", cex=1.5, lwd = 3)
      text(p$x+50, p$y+30, labels=p$trackId, cex=1.5)
    }
    legend("topleft",
           legend = c("Dead", "Transitional", "Alive", "Unstained"),
           pch=c(1,2,3,4),
           cex = 1)
    mtext(fileparts(fi[i+1])$name, cex = 2)
    cor_num <- filt %>% 
      filter(trackId == cell)
    cor_num <- cor_num$cor
    mtext(paste0("Pearson correlation of track ", cell, ":  ", cor_num), side = 1, cex = 2)
    dev.off()
  }
}

# fi <- list.files(paste0(root, filesep, well, filesep, "HALO MARKUP", filesep, "output"),full.names = T)
# for(cell in targets) {
#   filt <- cor_dat2 %>%
#     filter(trackId == cell | parentTrackId == cell)
#   for(i in min(filt$t):(max(filt$t)+5)) {
#     tiff(paste0("~/corr_HALO_", well, "_", cell, "_", i, ".tif"), width = 1468, height = 1100)
#     img=bioimagetools::readTIF(fi[i+1],as.is = T)
#     img <- EBImage::resize(img, dim(img)[1]/1)
#     plot(raster::as.raster(img[,,,1]))
#     if(i %in% unique(filt$t)) {
#       p <- filt %>% filter(t == i)
#       p$Class[p$Classifier.Phenotype=="UnStained"]=2
#       points(p$x, p$y+30, pch=(p$Class+2), col="white", cex=1.5, lwd = 3)
#       text(p$x+50, p$y+30, labels=p$trackId, col="white", cex=1.5)
#     }
#     mtext(fileparts(fi[i+1])$name, cex = 2)
#     dev.off()
#   }
# }

for(cell in targets){
  p <- cor_dat2 %>% 
    filter(trackId == cell)
  for(i in unique(p$t)) {
    size <- p %>% 
      filter(t == i)
    ggplot(p) +
      geom_line(aes(t, Object_Area_0)) +
      geom_point(aes(t, Object_Area_0, color = Classifier.Phenotype)) +
      scale_color_manual(values = col_values) +
      geom_point(aes(x = i, y = size[1,43]), color = "royalblue", size = 3) +
      theme_light() +
      ggtitle(paste0("Area over time for track ", cell))
    ggsave(paste0("corr_", cell, "_", i, ".tiff"), units="px", width=1468, height=1100, dpi=300)
  }
  q <- cor_dat2 %>% 
    filter(parentTrackId == cell)
  if(nrow(q) == 2) {
    for(i in 1:6) {
      ggplot(p) +
        geom_line(aes(t, Object_Area_0)) +
        geom_point(aes(t, Object_Area_0, color = Classifier.Phenotype)) +
        scale_color_manual(values = col_values) +
        theme_light() +
        ggtitle(paste0("Area over time for track ", cell))
      ggsave(paste0("corr_", cell, "_", (max(p$t) + i), ".tiff"), units="px", width=1468, height=1100, dpi=300)
    }
  } else {
      for(i in 1:5) {
        ggplot(p) +
          geom_line(aes(t, Object_Area_0)) +
          geom_point(aes(t, Object_Area_0, color = Classifier.Phenotype)) +
          scale_color_manual(values = col_values) +
          theme_light() +
          ggtitle(paste0("Area over time for track ", cell))
        ggsave(paste0("corr_", cell, "_", (max(p$t) + i), ".tiff"), units="px", width=1468, height=1100, dpi=300)
      }
  }
}

############################### inter-division cells 

orig_cells <- allcells_merged_dat1 %>% 
  filter(t == 0)
orig_cells <- orig_cells$trackId

target_cells_df <- allcells_merged_dat1 %>% 
  filter(trackId %!in% orig_cells & parent == TRUE)
target_cells_list <- unique(target_cells_df$trackId)

cor_div_dat1 <- data.frame()
for(cell in sort(unique(target_cells_df$trackId))) {
  df_1 <- filter(target_cells_df, trackId == cell)
  if(max(df_1$lifetime) > 24) {
    df_2 <- df_1 %>%
      mutate(cor = cor(t, Object_Area_0, method = "pearson"))
    cor_div_dat1 <- rbind(cor_div_dat1, df_2)
    print(cell)
  }
}

cor_div_dat1_list <- cor_div_dat1[order(-cor_div_dat1$cor),]
cor_div_dat1_list <- unique(cor_div_dat1_list$trackId)
targets <- cor_div_dat1_list[4:13]

cor_div_dat1_alltime <- allcells_merged_dat1 %>% 
  filter(trackId %in% targets | parentTrackId %in% targets)

cor_div_dat2 <- data.frame()
for(cell in sort(unique(cor_div_dat1_alltime$trackId))) {
  df_1 <- filter(cor_div_dat1_alltime, trackId == cell | parentTrackId == cell)
  if(max(df_1$lifetime) > 24) {
    df_2 <- df_1 %>%
      mutate(cor = cor(t, Object_Area_0, method = "pearson"))
    cor_div_dat2 <- rbind(cor_div_dat2, df_2)
    print(cell)
  }
}

well <- "F6_1"
col_values <- c("Dead" = "cyan", "Unstained" = "green", "Alive" = "orangered", "Transitional" = "magenta")
fi <- list.files(paste0(root, filesep, well, "/HALO Markup/output/overlay_output"),full.names = T)
for(cell in targets) {
  filt <- cor_div_dat2 %>%
    filter(trackId == cell | parentTrackId == cell)
  for(i in min(filt$t):(max(filt$t)+5)) {
    tiff(paste0("~/corr_div_", well, "_", cell, "_", i, ".tif"), width = 1468, height = 1100)
    img=bioimagetools::readTIF(fi[i+1],as.is = T)
    img <- EBImage::resize(img, dim(img)[1]/1)
    plot(raster::as.raster(img[,,,1]))
    if(i %in% unique(filt$t)) {
      p <- filt %>% filter(t == i)
      p$Class[p$Classifier.Phenotype=="UnStained"]=2
      points(p$x, p$y+30, pch=(p$Class+2), col="black", cex=1.5, lwd = 3)
      text(p$x+50, p$y+30, labels=p$trackId, cex=1.5)
    }
    legend("topleft",
           legend = c("Dead", "Transitional", "Alive", "Unstained"),
           pch=c(1,2,3,4),
           cex = 1)
    mtext(fileparts(fi[i+1])$name, cex = 2)
    cor_num <- filt %>% 
      filter(trackId == cell)
    cor_num <- cor_num$cor
    mtext(paste0("Pearson correlation of track ", cell, ":  ", cor_num), side = 1, cex = 2)
    dev.off()
  }
}


for(cell in targets){
  p <- cor_div_dat2 %>% 
    filter(trackId == cell)
  for(i in unique(p$t)) {
    size <- p %>% 
      filter(t == i)
    ggplot(p) +
      geom_line(aes(t, Object_Area_0)) +
      geom_point(aes(t, Object_Area_0, color = Classifier.Phenotype)) +
      scale_color_manual(values = col_values) +
      geom_point(aes(x = i, y = size[1,43]), color = "royalblue", size = 3) +
      theme_light() +
      ggtitle(paste0("Area over time for track ", cell))
    ggsave(paste0("corr_div_", cell, "_", i, ".tiff"), units="px", width=1468, height=1100, dpi=300)
  }
  q <- cor_div_dat2 %>% 
    filter(parentTrackId == cell)
  if(nrow(q) == 2) {
    for(i in 1:6) {
      ggplot(p) +
        geom_line(aes(t, Object_Area_0)) +
        geom_point(aes(t, Object_Area_0, color = Classifier.Phenotype)) +
        scale_color_manual(values = col_values) +
        theme_light() +
        ggtitle(paste0("Area over time for track ", cell))
      ggsave(paste0("corr_div_", cell, "_", (max(p$t) + i), ".tiff"), units="px", width=1468, height=1100, dpi=300)
    }
  } else {
    for(i in 1:5) {
      ggplot(p) +
        geom_line(aes(t, Object_Area_0)) +
        geom_point(aes(t, Object_Area_0, color = Classifier.Phenotype)) +
        scale_color_manual(values = col_values) +
        theme_light() +
        ggtitle(paste0("Area over time for track ", cell))
      ggsave(paste0("corr_div_", cell, "_", (max(p$t) + i), ".tiff"), units="px", width=1468, height=1100, dpi=300)
    }
  }
}
