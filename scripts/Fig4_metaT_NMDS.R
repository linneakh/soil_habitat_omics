##Linnea Honeker
##NMDS and PCA on metaT data from IMG
##2/26/24

library(ggfortify)
library(tidyverse)
library(factoextra)
library(ggnewscale)
library(ggrepel)
library(viridis)
library(vegan)
library(nlme)
library(ecodist)
library(ggplot2)
library(stats)
library(indicspecies)
library(cluster)
library(clusterProfiler)


source("./scripts/functions/NMDS_plotting_functions.R")


col_list = c("red", "orange", "green",  "blue",  "purple", "pink")
col_list_time = c("#F4A582", "#74ADD1",  "#D9EF8B")

Dir.f <- "./figures/Fig4_metaT_NMDS/"

# parameters for plots
h = 6
w = 8
res = 300
size = 16 #for poster

###inport data
#metaT de novo assembled and mapped to self (from IMG)
#bacteria
metaT.0 <- read.csv("./output/salazar_vst.csv", header = TRUE)
metaT <- metaT.0 %>%
  dplyr::select(-contains("216."),
           -contains("209."),
           -contains("203.")) %>%
  column_to_rownames(var = "X") %>%
  t() %>%
  as.data.frame()#remove outliers
 
#fungi
metaT.f.0 <- read.csv("./output/fungi_vst.csv", header = TRUE)
metaT.f <- metaT.f.0 %>%
  column_to_rownames(var = "X") %>%
  t() %>%
  as.data.frame()#remove outliers

#create metadata table for bacteria
metadata.b.0 <- as.data.frame(rownames(metaT)) 
colnames(metadata.b.0) <- "SampleID"

metadata.b <- metadata.b.0 %>%
  separate_wider_delim(SampleID,
                       delim = ".",
                       names = c("SampleID_num", "Treatment", "Zone", "Timepoint"), 
                       cols_remove = FALSE) %>%
  column_to_rownames(var = "SampleID") 

#create metadata table for fungi
metadata.f.0 <- as.data.frame(rownames(metaT.f)) 
colnames(metadata.f.0) <- "SampleID"

metadata.f <- metadata.f.0 %>%
  separate_wider_delim(SampleID,
                       delim = ".",
                       names = c("SampleID_num", "Treatment", "Zone", "Timepoint"), 
                       cols_remove = FALSE) %>%
  column_to_rownames(var = "SampleID") 




####-------bacterial nmds----------####
# Distance
dist <- vegdist(metaT, method = "euclidean")


# NMDS
set.seed(123)  # for reproducibility

nmds <- metaMDS(
  dist,
  k = 2,          # 2D solution
  trymax = 100,   # increase if convergence is difficult
  autotransform = FALSE
)


# Plot with your helpers (now site_scores carries Cluster from metadata_nmds)
nmds$stress

#extract coordinates
nmds_df <- as.data.frame(scores(nmds))
nmds_df$sample <- rownames(nmds_df)


# merge with metadata
nmds_df <- merge(nmds_df, metadata.b, by.x = "sample", by.y = "row.names") %>%
  select( -SampleID_num) %>%
  column_to_rownames(var = "sample")



# PERMANOVA (marginal effects)
distmat <- dist

adonis2(distmat ~ Zone+Treatment+Timepoint, metadata.b, permutations = 999, by = "margin")
# adonis2(formula = distmat ~ Zone + Treatment + Timepoint, data = metadata, permutations = 999, by = "margin")
# Df SumOfSqs      R2      F Pr(>F)    
# Zone       3    46563 0.22923 7.6807  0.001 ***
#   Treatment  1    13920 0.06853 6.8886  0.001 ***
#   Timepoint  2    11845 0.05831 2.9308  0.001 ***
#   Residual  64   129329 0.63670                  
# Total     70   203123 1.00000          




### nmds plot (timepoint and moisture)
w=10
h=8

#plot zone and treatment
make_nmds_plot(nmds_df, Group1 =Zone, Group2 =Treatment,col_list=col_list)

filename <- paste0(Dir.f, "Fig4_bacteria_nmds_zone_treatment.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)
filename <- paste0(Dir.f, "Fig4_bacteria_nmds_zone_treatment.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)

#plot zone and time
make_nmds_plot(nmds_df, Group1 =Zone, Group2 =Timepoint,col_list=col_list)

filename <- paste0(Dir.f, "bacteria_nmds_zone_timepoint.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)
filename <- paste0(Dir.f, "bacteria_nmds_zone_timepoint.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)


####-------fungal nmds----------####
# Distance
dist <- vegdist(metaT.f, method = "euclidean")


# NMDS
set.seed(123)  # for reproducibility

nmds <- metaMDS(
  dist,
  k = 2,          # 2D solution
  trymax = 100,   # increase if convergence is difficult
  autotransform = FALSE
)


# Plot with your helpers (now site_scores carries Cluster from metadata_nmds)
nmds$stress

#extract coordinates
nmds_df <- as.data.frame(scores(nmds))
nmds_df$sample <- rownames(nmds_df)


# merge with metadata
nmds_df <- merge(nmds_df, metadata.f, by.x = "sample", by.y = "row.names", all = TRUE) %>%
  select( -SampleID_num) %>%
  column_to_rownames(var = "sample")

metadata.f <- metadata.f %>%
  filter(rownames(.) %in% rownames(nmds_df))

# PERMANOVA (marginal effects)
distmat <- dist

adonis2(distmat ~ Zone+Treatment+Timepoint, metadata.f, permutations = 999, by = "margin")
# adonis2(formula = distmat ~ Zone + Treatment + Timepoint, data = metadata.f, permutations = 999, by = "margin")
# Df SumOfSqs      R2      F Pr(>F)    
# Zone       3    93141 0.14058 3.8404  0.001 ***
#   Treatment  1    14685 0.02216 1.8164  0.021 *  
#   Timepoint  2    20359 0.03073 1.2591  0.126    
# Residual  66   533573 0.80532                  
# Total     72   662563 1.00000   




### nmds plot (timepoint and moisture)
w=10
h=8

#plot zone and treatment
make_nmds_plot(nmds_df, Group1 =Zone, Group2 =Treatment,col_list=col_list)

filename <- paste0(Dir.f, "Fig4_fungi_nmds_zone_treatment.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)
filename <- paste0(Dir.f, "Fig4_fungi_nmds_zone_treatment.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)

#plot zone and time
make_nmds_plot(nmds_df, Group1 =Zone, Group2 =Timepoint,col_list=col_list)

filename <- paste0(Dir.f, "fungi_nmds_zone_timepoint.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)
filename <- paste0(Dir.f, "fungi_nmds_zone_timepoint.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)


