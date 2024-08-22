################################################################################
# In order to update this code for what you want to run, you will need to change
# ...the following lines (as of 8-15):
# 12, 13, 15, 181, 194, 227, 232, 234, 239-252, 319-330, 339-401, 
# This does not include changing file paths to match your layout.
# If I can't be reached, try jordan.elizabeth9452@gmail.com or +1 615-500-8301.
################################################################################

library(tidyverse)

# Project folder: all input and output goes through here
root <- "C:/Users/80027908/Desktop/K00_GemcitabineExposure_033023"
filesep <- "/"
# Set row to find row's tracking data and classifier data
row <- "E"

# Read in tracking results
tracking_csvs <- list.files(paste0(root, filesep, "Tracking_CSVs_", row), full.names = T)
# Read in the classifier file 
class <- read.csv(paste0(root, filesep, "Object_data", filesep, row, "_row_Object_Data.csv"))


################################################################################
# Merge tracking csv with classifier
# Location inputs can be row, well, or field
# Max time sets the number of timepoints to be processed
# Save over existing can be true or false
################################################################################

mergeData <- function(location, max_time, save_over_existing) {
  # If input is row letter...
  if(nchar(location) == 1) {
    print("Loading in row classifier and tracking csvs")
    row <- location
    class_1 <- class %>% 
      separate(col = Image.Location, 
               into = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "well_info", NA), 
               sep = "\\\\", 
               remove = FALSE)
    
    data <- list()
    # ...load in all data in row.
    for(j in 2:11) {
      for(i in 1:4) {
        dat <- read.csv(tracking_csvs[grepl(paste0(row, j, "_", i), tracking_csvs)]) %>% 
          mutate(well_info = paste0(row, j, "_", i))
        assign(paste0("dat_", row, j, "_", i), dat, .GlobalEnv)
        data <- append(list(dat), data)
      }
    }
  # If input is row and column (ie. well)...
  } else if(nchar(location) == 2) {
    print("Loading in well classifier and tracking csvs")
    well <- location
    class_1 <- class %>% 
      filter(grepl(well, Image.Location)) %>% 
      separate(col = Image.Location, 
               into = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "well_info", NA), 
               sep = "\\\\", 
               remove = FALSE)
    
    row <- strsplit(location, "")[[1]][1]
    data <- list()
    # ...load in all fields in well.
    for(i in 1:4) {
      dat <- read.csv(tracking_csvs[grepl(paste0(well, "_", i), tracking_csvs)]) %>% 
        mutate(well_info = paste0(well, "_", i))
      assign(paste0("dat_", well, "_", i), dat, .GlobalEnv)
      data <- append(list(dat), data)
    }
  # If input is row, column, and field...
  } else if(nchar(location) == 4) {
    print("Loading in field classifier and tracking csv")
    field <- location
    class_1 <- class %>% 
      filter(grepl(field, Image.Location)) %>% 
      separate(col = Image.Location, 
               into = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "well_info", NA), 
               sep = "\\\\", 
               remove = FALSE)
    # ...load the field.
    dat <- read.csv(tracking_csvs[grepl(field, tracking_csvs)]) %>% 
      mutate(well_info = field)
    assign(paste0("dat_", field), dat, .GlobalEnv)
    
    data <- list()
    data <- append(list(dat), data)
  } else {
    print("Input location does not match previous patterns. Location should be a row, well, or field.")
  }
  
  assign("data", data, .GlobalEnv)
  assign("class_1", class_1, .GlobalEnv)
  
  ##############################################################################
  print("Processing classifier")
  # Pull in the data for the well of interest
  class_1$time <- sapply(strsplit(class_1$Image.Location,"_"), function(x) x[length(x)])
  class_1$time <- gsub(".tif","",class_1$time)
  
  # Change times to frames
  map <- 0:(length(unique(class_1$time))-1)
  names(map) <- sort(unique(class_1$time))
  class_1$t=map[class_1$time]
  class_1 <- class_1[order(class_1$t),]
  
  #Create an "is not in" function
  '%!in%' <- function(x,y)!('%in%'(x,y))
  ##############################################################################
  
  print("Merging classifier and tracking csv")
  for(frame in 1:length(data)) {
    dat <- data[[frame]]
    # If the merged data does not exist or if save over existing is true...
    if(file.exists(paste0(root, filesep, "allcells_", dat[1,63], ".csv")) == FALSE | save_over_existing == TRUE) {
      print(paste0("Merging: ", dat[1,63]))
      # Rename columns
      ii=match(c("frame","Bounding_Box_Maximum_0","Bounding_Box_Maximum_1"),colnames(dat))
      colnames(dat)[ii]=c("t","x","y")
      
      seg <- class_1 %>% 
        filter(well_info == dat$well_info[1]) %>% 
        group_by(t) %>% 
        filter(t <= max_time)
      dat <- dat %>% 
        filter(t <= max_time)
      
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
      
      print(paste0("Marking: ", dat[1,63]))
      # Find end point class
      merged_1 <- filter(merged, trackId != -1)
      cells_list <- unique(merged_1$trackId)
      # Recreate parents list
      parents_list <- unique(merged$parentTrackId)
      parents_list <- parents_list[!parents_list %in% 0]
      
      # Add "parent" and "dividing" info
      allcells <- merged %>% 
        filter(trackId != -1) %>% 
        group_by(trackId) %>%
        mutate(lifetime = row_number(), 
               parent = case_when(trackId %in% parents_list ~ TRUE,
                                  .default = FALSE),
               dividing = case_when(parent == TRUE & t == max(t) ~ TRUE,
                                    .default = FALSE))
      
      # Save dataframe
      name <- paste("allcells", dat[1,63], sep = "_")
      assign(name, allcells, .GlobalEnv)
      write.csv(allcells, file = paste0(root, filesep, "allcells_", dat[1,63], ".csv"))
    }
  }
  print("Done")
}

################################################################################
# Run the function
mergeData("E4", 40, TRUE)
################################################################################

# All code below this point can only run with all fields in a well; no rows or
# ...individual fields without hardcoding those changes.

################################################################################
# Create track txt files for WGD model
################################################################################
dir.create(file.path(root, "Post-division"))
dir.create(file.path(root, "Inter-division"))
dir.create(file.path(root, "Orig-division"))

# Set which data to create files for
well <- "E4"
data_list <- c(paste0("allcells_", well, "_1"),
               paste0("allcells_", well, "_2"),
               paste0("allcells_", well, "_3"),
               paste0("allcells_", well, "_4"))

# For each field in list...
for(data in data_list) {
  df <- get(data)
  for(cell in unique(df$trackId)) {
    df <- get(data)
    row.names(df) <- NULL
    df <- filter(df, trackId == cell)
    if(min(df$t) != 0) {
      if(df$parent[1] == TRUE) {
        # ...get tracks that have and will divide.
        write.table(df_filt, paste0(root, "Inter-division/", data, "_", cell, "_inter-division.txt"), sep = "\t")
      }
      # ...get tracks that have divided (includes inter-division).
      write.table(df_filt, paste0(root, "Post-division/", data, "_", cell, "_post-division.txt"), sep = "\t")
    } else {
      # ...get original tracks (tracks that do not have parents and may or may not divide).
      write.table(df_filt, paste0(root, "Orig-division/", data, "_", cell, "_orig-division.txt"), sep = "\t")
    }
  }
}

################################################################################
# Create visuals of WGD results
################################################################################
library(patchwork)
library(matlab)

# t_val <- "t1"

row <- strsplit(well, "")[[1]][1]

# Load in WGD track assignment txt files
track_list <- list.files(paste0(root, filesep, "Misc", filesep, "track_IDS"), pattern = paste0(well, "_[0-9]_wgd_df"), full.names = T)
# If using t_val, use the line below
# track_list <- list.files(paste0(root, filesep, "Misc", filesep, "track_IDS"), pattern = paste0(well, "_[0-9]_wgd_df_", t_val), full.names = T)

# Reformat WGD data and merge
# Set number of dats according to how many WGD files there are... hopefully this 
# ...can be aged out once the column names of WGD0 are fixed.
dat1 <- read.table(track_list[1])
# colnames(dat1) <- c("track", "t", "x", "y", "Object_Area_0", "Classifier.Phenotype", "WGD")
# colnames(dat1) <- c("track", "t", "x", "y", "Size_in_pixels_0", "Classifier.Phenotype", "WGD")
# tracks <- dat1
dat2 <- read.table(track_list[2])
dat2 <- rownames_to_column(dat2, var = "track")
colnames(dat1) <- colnames(dat2)
dat3 <- read.table(track_list[3])
dat3 <- rownames_to_column(dat3, var = "track")
dat4 <- read.table(track_list[4])
dat4 <- rownames_to_column(dat4, var = "track")
tracks <- rbind(dat1, dat2, dat3, dat4)
# tracks <- rbind(dat1, dat2, dat3)
# tracks <- rbind(dat1, dat2)

# Pull out metadata
tracks <- tracks %>% 
  separate(col = track, 
           into = c(NA, "well", "field", "trackId", NA), 
           sep = "_", 
           remove = TRUE)
tracks$trackId <- as.numeric(tracks$trackId)
tracks$well_info <- paste0(tracks$well, "_", tracks$field)

# For visuals, generate correlations 
cor <- data.frame()
for(cell in sort(unique(tracks$trackId))) {
  df <- filter(tracks, trackId == cell) %>%
    mutate(cor = cor(t, Size_in_pixels_0, method = "pearson"))
  cor <- rbind(cor, df)
}

q <- rbind(get(paste0("allcells_", well, "_1")),
           get(paste0("allcells_", well, "_2")),
           get(paste0("allcells_", well, "_3")),
           get(paste0("allcells_", well, "_4")))
p <- left_join(cor, q)

p <- p %>% 
  group_by(trackId, well_info, WGD, parent) %>% 
  mutate(lifespan = length(t)) %>% 
  ungroup()

ggplot(p) +
  geom_histogram(aes(cor)) +
  facet_grid(cols = vars(WGD), rows = vars(parent))

ggplot(p) +
  geom_boxplot(aes(as.character(WGD), lifespan)) +
  facet_grid(cols = vars(parent))

ggplot(p) +
  geom_point(aes(lifespan, cor, color = as.character(WGD)))

# Create final output folders
dir.create(file.path(root, filesep, "Target_plots"))
dir.create(file.path(root, filesep, "Target_images"))
dir.create(file.path(root, filesep, "Target_images_halo"))
dir.create(file.path(root, filesep, "Target_csvs"))

# For every field in the WGD data...
for(fields in unique(tracks$field)) {
  # Get the data
  allcells_merged <- get(paste0("allcells_", well, "_", fields))

  # Until "cor", this code is kinda moot unless you want the parent and dividing
  # ...labels. Otherwise just use the whole "tracks" dataframe and not "allcells"
  tracks_filt <- tracks %>% 
    filter(field == fields) %>% 
    select("trackId", "t", "x", "y", "Size_in_pixels_0", "Classifier.Phenotype", "WGD")
  
  merged <- allcells_merged %>%
    filter(trackId %in% unique(tracks_filt$trackId))
  
  tracks_merged <- full_join(tracks_filt, merged, keep = FALSE)
  
  cor <- data.frame()
  for(cell in sort(unique(tracks_merged$trackId))) {
    df <- filter(tracks_merged, trackId == cell) %>%
      # Get the coordinates for filtering from edge
      mutate(cor = cor(t, Size_in_pixels_0, method = "pearson"),
             xmin = (min(x)),
             xmax = (max(x)),
             ymin = (min(y)),
             ymax = (max(y))) %>%
      # Filter cells that approach edge within lifetime
      filter(xmin > 50 &
               ymin > 50 &
               xmax < 1418 &
               ymax < 1050)
    # Exclude cells that divide or disappear before t30
    if(max(df$t) < 30) {
      next
    } else {
      cor <- rbind(cor, df)
    }
  }
  
  # This could be written much better
  # Numbers of "cor_list" must match the number of WGD types
  cor_list_0 <- cor %>%
    filter(WGD == 0) %>%
    group_by(trackId) %>%
    filter(t == max(t)) %>%
    reframe(trackId, cor = cor)
  
  cor_list_1 <- cor %>%
    filter(WGD == 1) %>%
    group_by(trackId) %>%
    filter(t == max(t)) %>%
    reframe(trackId, cor = cor)
  
  cor_list_2 <- cor %>%
    filter(WGD == 2) %>%
    group_by(trackId) %>%
    filter(t == max(t)) %>%
    reframe(trackId, cor = cor)
  
  cor_list_3 <- cor %>%
    filter(WGD == 3) %>%
    group_by(trackId) %>%
    filter(t == max(t)) %>%
    reframe(trackId, cor = cor)
  
  # Select tracks with correlations approaching the listed values
  cor_list_0 <- cor_list_0 %>%
    filter(abs(cor-0.3)==min(abs(cor-0.3)) |
             abs(cor-0.4)==min(abs(cor-0.4)) |
             abs(cor-0.5)==min(abs(cor-0.5)) |
             abs(cor-0.6)==min(abs(cor-0.6)) |
             abs(cor-0.7)==min(abs(cor-0.7)) |
             abs(cor-0.8)==min(abs(cor-0.8)) |
             abs(cor-0.9)==min(abs(cor-0.9)))
  
  cor_list_1 <- cor_list_1 %>%
    filter(abs(cor-0.3)==min(abs(cor-0.3)) |
             abs(cor-0.4)==min(abs(cor-0.4)) |
             abs(cor-0.5)==min(abs(cor-0.5)) |
             abs(cor-0.6)==min(abs(cor-0.6)) |
             abs(cor-0.7)==min(abs(cor-0.7)) |
             abs(cor-0.8)==min(abs(cor-0.8)) |
             abs(cor-0.9)==min(abs(cor-0.9)))
  
  cor_list_2 <- cor_list_2 %>%
    filter(abs(cor-0.3)==min(abs(cor-0.3)) |
             abs(cor-0.4)==min(abs(cor-0.4)) |
             abs(cor-0.5)==min(abs(cor-0.5)) |
             abs(cor-0.6)==min(abs(cor-0.6)) |
             abs(cor-0.7)==min(abs(cor-0.7)) |
             abs(cor-0.8)==min(abs(cor-0.8)) |
             abs(cor-0.9)==min(abs(cor-0.9)))
  
  cor_list_3 <- cor_list_3 %>%
    filter(abs(cor-0.3)==min(abs(cor-0.3)) |
             abs(cor-0.4)==min(abs(cor-0.4)) |
             abs(cor-0.5)==min(abs(cor-0.5)) |
             abs(cor-0.6)==min(abs(cor-0.6)) |
             abs(cor-0.7)==min(abs(cor-0.7)) |
             abs(cor-0.8)==min(abs(cor-0.8)) |
             abs(cor-0.9)==min(abs(cor-0.9)))
  
  # Get a list of names for these tracks
  cor_list <- c(unique(cor_list_0$trackId), unique(cor_list_1$trackId), unique(cor_list_2$trackId), unique(cor_list_3$trackId))
  
  # If there are no cells that meet the criteria...
  if(length(cor_list) == 0) next
  
  # Create info for import into ImageJ splicer
  for(cell in cor_list) {
    print(cell)
    cor_num <- cor %>%
      filter(trackId == cell)
    cor_num <- cor_num$cor
    p <- cor %>%
      filter(trackId == cell) %>%
      select("trackId", "t", "x", "y", "xmin", "xmax", "ymin", "ymax", "Size_in_pixels_0", "Classifier.Phenotype", "WGD") %>% 
      mutate(cor = cor_num)
    p$track_info <- paste0(p$well, "_", p$field, "_", p$trackId)
    rownames(p) <- NULL
    write.csv(p, file = paste0(root, filesep, "Target_csvs", filesep, well, "_", fields, "_", cell, ".csv", sep = ""))
  }
  
  # Load brightfield images
  fi <- list.files(paste0(root, row, "_row_images/", well, "_", fields), full.names = T)
  
  # For all the target cells in the list...
  for(cell in cor_list) {
    # Get the cell's track
    filt <- allcells_merged %>%
      filter(trackId == cell)
    start <- (min(filt$t) - 1)
    # Get the sibling track of the target cell
    filt_1 <- allcells_merged %>%
      filter(parentTrackId == filt$parentTrackId[1] & trackId != cell)
    # Get the parent track of the target cell
    filt <- allcells_merged %>% 
      filter(trackId == cell | parentTrackId == cell | trackId == filt$parentTrackId[1]) %>% 
      filter(t >= start) %>% 
      rbind(filt_1)
    
    # For all time points plus five more, plot the track on the image.
    for(j in min(filt$t):(max(filt$t)+5)) {
      cor_num <- cor %>%
        filter(trackId == cell)
      cor_num <- cor_num$cor
      tiff(paste0(root, filesep, "Target_images", filesep, well, "_", fields, "_", cell, "_", j, ".tif"), width = 1645, height = 1232)
      img=bioimagetools::readTIF(fi[j+1],as.is = T)
      img <- EBImage::resize(img, dim(img)[1]/1)
      plot(raster::as.raster(img[,,,1]))
      if(j %in% unique(filt$t)) {
        p <- filt %>% filter(t == j)
        p$Class[p$Classifier.Phenotype=="UnStained"]=2
        points(p$x, p$y+30, pch=(p$Class+2), col="black", cex=1.5, lwd = 3)
        text(p$x+50, p$y+30, labels=p$trackId, cex=1.5)
      }
      # Add a legend and info, if you are not cropping the images
      # legend("topleft",
      #       legend = c("Dead", "Transitional", "Alive", "Unstained"),
      #       pch=c(1,2,3,4),
      #       cex = 1)
      # mtext(fileparts(fi[j+1])$name, cex = 2)
      # mtext(paste0("Pearson correlation of track ", cell, ":  ", cor_num), side = 1, cex = 2)
      dev.off()
    }
    
    # Load HALO images
    fih <- list.files(paste0(root, filesep, well, "_", fields), full.names = T)
    
    # For all time points plus five more, plot the track on the HALO image.
    for(j in min(filt$t):(max(filt$t)+5)) {
      if(j+1 <= 40) {
        cor_num <- cor %>%
          filter(trackId == cell)
        cor_num <- cor_num$cor
        tiff(paste0(root, filesep, "Target_images_halo", filesep, well, "_", fields, "_", cell, "_", j, ".tif"), width = 1645, height = 1232)
        img=bioimagetools::readTIF(fih[j+1],as.is = T)
        img <- EBImage::resize(img, dim(img)[1]/1)
        plot(raster::as.raster(img[,,,1]))
        if(j %in% unique(filt$t)) {
          p <- filt %>% filter(t == j)
          p$Class[p$Classifier.Phenotype=="UnStained"]=2
          points(p$x, p$y+30, pch=(p$Class+2), col="white", cex=1.5, lwd = 3)
          text(p$x+50, p$y+30, labels=p$trackId, col="white", cex=1.5)
        }
        # Add a legend and info, if you are not cropping the images
        # legend("topleft",
        #       legend = c("Dead", "Transitional", "Alive", "Unstained"),
        #       pch=c(1,2,3,4),
        #       cex = 1)
        # mtext(fileparts(fi[j+1])$name, cex = 2)
        # mtext(paste0("Pearson correlation of track ", cell, ":  ", cor_num), side = 1, cex = 2)
        dev.off()
      } else {
        next
      }
    }
    
    # Set color scheme for cell class
    col_values <- c("Dead" = "cyan", "Unstained" = "green", "Alive" = "orangered", "Transitional" = "magenta")
    filt <- cor %>%
      filter(trackId == cell) 
    
    # I don't remember why this code right below this line is here so it stays...
    # filt <- filt[1:(max(filt$t) - min(filt$t)),]
   
    # For every time point in the track...
    for(i in unique(filt$t)) {
      size <- filt %>%
        filter(t == i)
      size <- unique(size)
      # Set the color of each WGD type
      if(size$WGD == 0) {
        col <- "red"
      }
      if(size$WGD == 1) {
        col <- "green"
      }
      if(size$WGD == 2) {
        col <- "blue"
      }
      if(size$WGD == 3) {
        col <- "purple"
      }
      
      # Load in a hand-cropped image of a legend 
      diy_legend <- jpeg::readJPEG(paste0(root, filesep, "diy_legend.jpg"), native = TRUE)
      
      # Graph plots
      ggplot(filt) +
        geom_line(aes(t, Object_Area_0)) +
        geom_line(aes(t, Size_in_pixels_0), color = "blue") +
        geom_point(aes(t, Object_Area_0, color = Classifier.Phenotype)) +
        scale_color_manual(values = col_values) +
        geom_point(aes(x = i, y = size$Object_Area_0), color = col, size = 3) +
        theme_light() +
        ggtitle(paste0("Area over time for track ", cell)) +
        inset_element(p = diy_legend, left = 0.8, bottom = 0.05, right = .99, top = 0.2)
      ggsave(path = paste0(root, filesep, "Target_plots", filesep), filename = paste0(well, "_", fields, "_", cell, "_", i, ".tiff"), units="px", width=1468, height=1100, dpi=300)
    }
    
    # Now graph extra plots... either 5 or 6 depending on if a track ends with
    # ...division (6) or terminates (5).
    q <- allcells_merged %>% 
      filter(parentTrackId == cell)
    if(nrow(q) == 2) {
      for(i in 1:6) {
        ggplot(filt) +
          geom_line(aes(t, Size_in_pixels_0), color = "blue") +
          geom_point(aes(t, Size_in_pixels_0, color = Classifier.Phenotype)) +
          scale_color_manual(values = col_values) +
          theme_light() +
          ggtitle(paste0("Area over time for track ", cell))
        ggsave(path = paste0(root, filesep, "Target_plots", filesep), filename = paste0(well, "_", fields, "_", cell, "_", (max(filt$t) + i), ".tiff"), units="px", width=1468, height=1100, dpi=300)
      }
    } else {
      for(i in 1:5) {
        ggplot(filt) +
          geom_line(aes(t, Size_in_pixels_0), color = "blue") +
          geom_point(aes(t, Size_in_pixels_0, color = Classifier.Phenotype)) +
          scale_color_manual(values = col_values) +
          theme_light() +
          ggtitle(paste0("Area over time for track ", cell))
        ggsave(path = paste0(root, filesep, "Target_plots", filesep), filename = paste0(well, "_", fields, "_", cell, "_", (max(filt$t) + i), ".tiff"), units="px", width=1468, height=1100, dpi=300)
      }
    }
    ggplot(filt) +
      geom_line(aes(t, Size_in_pixels_0), color = "blue") +
      geom_point(aes(t, Size_in_pixels_0, color = Classifier.Phenotype)) +
      scale_color_manual(values = col_values) +
      theme_light() +
      ggtitle(paste0("Area over time for track ", cell))
    ggsave(path = paste0(root, filesep, "Target_plots", filesep), filename = paste0(well, "_", fields, "_", cell, "_", (min(filt$t) - 1), ".tiff"), units="px", width=1468, height=1100, dpi=300)
  }
}

# All done! You should have generated:
# allcells csvs
# track txts
# target plots, images, halo images, and csvs

# The next step is to run the splicing script in ImageJ. To do so, open ImageJ
# ...and do the following:
# Plugins -> Macros -> Edit -> open the .ijm
# Then, hit the "run" button at the bottom!
```