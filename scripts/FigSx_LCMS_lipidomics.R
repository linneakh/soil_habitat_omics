##Linnea Honeker
##honeker1@llnl.gov
##LC-MS analysis


##load libraries
library(readxl)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(factoextra)
library(ggfortify)
library(ggnewscale)
library(ggrepel)
library(viridis)
library(vegan)
library(remotes)
library(ggraph)
library(extrafont)
library(pmartR)
library(pals)
library(RColorBrewer)
library(colorspace)
library(clusterProfiler)

# # ── Color definitions ─────────────────────────────────────────────────────────
class_base_colors <- c(
  "Diacylglycerol"                       = "#A052B8",  # lighter deep purple
  "Triacylglycerol"                  = "#F0C060",  # lighter golden amber
  "Phosphatidylcholine"                = "#56A8D0",  # lighter steel blue
  #"Diacylglycercyl-N,N,N-trimethylhomoserine"   = "#E8707A",  # lighter raspberry
  "Phosphatidylethanolamine"                = "#7AB050",  # lighter olive green
  "Phosphatidylglycerol"                   = "#FF9A40",  # lighter deep orange
  "Phosphatidylinositol"                = "#40E0D0",  # lighter turquoise
  "Ceramide"                    = "#F8EE70"
)
  

zone_colors      <- c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red")
treatment_colors <- c(Drought = "brown", Untrt = "darkgreen")
timepoint_colors <- c(`4weeks` = "#F4A582", `8weeks` = "#74ADD1", `12weeks` = "#D9EF8B")



##load LC-MS data
file_path <- "./data/LC_MS_lipidomics.xlsx"
#kegg <- read.csv("./qSIP_output/LC-MS/HPOS_metaboanalyst_KEGGID_all.csv", header = TRUE)
#kegg$Query <- gsub(" ", "_", kegg$Query) 


##load extra functions
source("./scripts/functions/NMDS_plotting_functions.R")
source("./scripts/functions/lipidomics_functions.R")

# Get sheet names
sheet_names <- excel_sheets(file_path)
#"Groups"         "HILIC Positive" "HILIC Negative" "RP Positive"    "RP Negative"   

# Read each sheet into a data frame and store in a list
data_frames_list <- lapply(sheet_names, function(sheet_name) {
  read_excel(file_path, sheet = sheet_name, range = NULL)
})

# Extract the first data frame
peak_table_pos.0 = data_frames_list[[3]]

# Extract the second data frame
peak_table_neg.0 <- data_frames_list[[4]]

# ── Load and clean compound classifications ───────────────────────────────────
meta.comp <- rbind(
  peak_table_pos.0 %>% select(`Common Name`),
  peak_table_neg.0 %>% select(`Common Name`)
) %>%
  unique() %>%
  mutate(
    `Class` = case_when(
      grepl("Cer", `Common Name`) ~ "Ceramide",
      grepl("DG", `Common Name`) ~ "Diacylglycerol",
      grepl("TG", `Common Name`) ~ "Triacylglycerol",
      grepl("PC", `Common Name`) ~ "Phosphatidylcholine",
      grepl("DGTSA", `Common Name`) ~ "Diacylglycercyl-N,N,N-trimethylhomoserine",
      grepl("PE", `Common Name`) ~ "Phosphatidylethanolamine",
      grepl("PG", `Common Name`) ~ "Phosphatidylglycerol",
      grepl("PI", `Common Name`) ~ "Phosphatidylinositol"
    ),
  ) %>%
  column_to_rownames(var = "Common Name")
  



##HPOS
#reformat and filter pos table
colnames(peak_table_pos.0)[1] <- "mz"
colnames(peak_table_pos.0)[2] <- "rt"
colnames(peak_table_pos.0) = gsub("_Pos", "", colnames(peak_table_pos.0))
colnames(peak_table_pos.0) = gsub(" ", "", colnames(peak_table_pos.0))


#order column names to match neg table
peak_table_pos2 <-  peak_table_pos.0 %>%
  # unite(mz_name, c("mz", "CommonName"), sep = "_", remove = FALSE) %>%
  column_to_rownames(var = "CommonName") %>%
  select(order(colnames(.))) %>%
  rownames_to_column(var = "Name") %>%
  select(!contains("Pool")) %>%
  select(!contains("SRFA"))

colnames(peak_table_pos2) <- str_remove(colnames(peak_table_pos2), "_RP")


peak_table_mean <- peak_table_pos2 %>%
  #select(-mz, -rt) %>%
  group_by(Name) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) 

edata.0 <- peak_table_mean %>%
  select(-mz, -rt)


emeta.pos <-  peak_table_mean %>%
  select(Name, mz, rt) 

fdata.0 <- edata.0 %>%
  column_to_rownames(var = "Name") 


fdata <- as.data.frame(t(fdata.0)) %>%
  rownames_to_column(var = "SampleID") %>%
  dplyr::select (SampleID) %>%
  tidyr::separate_wider_delim(SampleID, delim = "_", 
                              names = c("DRIP", "ID", "Treatment0"),
                              cols_remove = FALSE) %>%
  mutate(Zone = case_when(
    startsWith(Treatment0, "RDD") ~ "RhizoDetritus",
    startsWith(Treatment0, "RDN") ~ "RhizoDetritus",
    startsWith(Treatment0, "R") ~ "Rhizo",
    startsWith(Treatment0, "B") ~ "Bulk",
    startsWith(Treatment0, "D") ~ "Detritus")) %>%
  mutate(Treatment = case_when(
    startsWith(Treatment0, "RDD") ~ "Drought",
    startsWith(Treatment0, "RDN") ~ "Untrt",
    startsWith(Treatment0, "RD") ~ "Drought",
    startsWith(Treatment0, "RN") ~ "Untrt",
    startsWith(Treatment0, "BD") ~ "Drought",
    startsWith(Treatment0, "BN") ~ "Untrt",
    startsWith(Treatment0, "DD") ~ "Drought",
    startsWith(Treatment0, "DN") ~ "Untrt")) %>%
  mutate(Timepoint = case_when(
    grepl("1", Treatment0) ~ "4weeks",
    grepl("2", Treatment0) ~ "8weeks",
    grepl("3", Treatment0) ~ "12weeks")) %>%
  #filter(!grepl("SRFA1", SampleID)) %>%
  select(-DRIP) %>%
  mutate(ID = str_remove(ID, "R")) %>%
  mutate(ID = str_remove(ID, "D")) %>%
  mutate(ID = str_remove(ID, "B")) %>%
  mutate(ID = str_remove(ID, "N")) %>%
  mutate(ID = paste("DRIP", ID, sep = "_")) %>%
  mutate(ID = case_when(
    Zone == "Bulk" ~ paste(ID, "b", sep = ""),
    TRUE ~ ID)) %>%
  unite(SampleID_new, c(ID, Treatment, Zone, Timepoint), sep = ".", remove = FALSE) %>%
  mutate(SampleType = case_when(
    grepl("ExBlk", SampleID) ~ "Blank",
    TRUE ~ "Sample"
  ))



edata <- edata.0[colnames(edata.0)[colnames(edata.0) %in% fdata$SampleID]]
edata <- cbind(edata.0[,1], edata)

#save meta-omics file with same names as other datasets

edata.t <- edata %>%
  column_to_rownames("Name") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%     # do this once
  left_join(fdata, by = "SampleID") %>%  # use dplyr join (keeps order, clearer)
  select(-c(ID, Treatment0, Zone, Treatment, Timepoint, SampleID, SampleType))  %>%# drop what you don’t need
  column_to_rownames(var="SampleID_new")

edata.pos <- edata.t %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var= "Name")

#remove features that are not at least 3x greater than median of blank
h_pos_clean <- blank_fold_filter_raw(
  edata.pos, 
  fdata,  
  edata_cname = "Name",
  fdata_cname = "SampleID_new",
  blank_label = "Blank",
  sample_type_cname = "SampleType",
  fold_threshold = 3,
  remove_blanks = TRUE
)
#no compound removed


#reformat and filter neg table
colnames(peak_table_neg.0)[1] <- "mz"
colnames(peak_table_neg.0)[2] <- "rt"
colnames(peak_table_neg.0) = gsub("_Pos", "", colnames(peak_table_neg.0))
colnames(peak_table_neg.0) = gsub(" ", "", colnames(peak_table_neg.0))


#order column names to match neg table
peak_table_neg2 <-  peak_table_neg.0 %>%
  # unite(mz_name, c("mz", "CommonName"), sep = "_", remove = FALSE) %>%
  column_to_rownames(var = "CommonName") %>%
  select(order(colnames(.))) %>%
  rownames_to_column(var = "Name") %>%
  select(!contains("Pool")) %>%
  select(!contains("SRFA"))

colnames(peak_table_neg2) <- str_remove(colnames(peak_table_neg2), "_RP_Neg")

peak_table_mean <- peak_table_neg2 %>%
  #select(-mz, -rt) %>%
  group_by(Name) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) 

edata.0 <- peak_table_mean %>%
  select(-mz, -rt)


emeta.neg <-  peak_table_mean %>%
  select(Name, mz, rt) 

fdata.0 <- edata.0 %>%
  column_to_rownames(var = "Name") 


fdata <- as.data.frame(t(fdata.0)) %>%
  rownames_to_column(var = "SampleID") %>%
  dplyr::select (SampleID) %>%
  filter(!grepl("Neg", SampleID)) %>%
  tidyr::separate_wider_delim(SampleID, delim = "_", 
                              names = c("DRIP", "ID", "Treatment0"),
                              cols_remove = FALSE) %>%
  mutate(Zone = case_when(
    startsWith(Treatment0, "RDD") ~ "RhizoDetritus",
    startsWith(Treatment0, "RDN") ~ "RhizoDetritus",
    startsWith(Treatment0, "R") ~ "Rhizo",
    startsWith(Treatment0, "B") ~ "Bulk",
    startsWith(Treatment0, "D") ~ "Detritus")) %>%
  mutate(Treatment = case_when(
    startsWith(Treatment0, "RDD") ~ "Drought",
    startsWith(Treatment0, "RDN") ~ "Untrt",
    startsWith(Treatment0, "RD") ~ "Drought",
    startsWith(Treatment0, "RN") ~ "Untrt",
    startsWith(Treatment0, "BD") ~ "Drought",
    startsWith(Treatment0, "BN") ~ "Untrt",
    startsWith(Treatment0, "DD") ~ "Drought",
    startsWith(Treatment0, "DN") ~ "Untrt")) %>%
  mutate(Timepoint = case_when(
    grepl("1", Treatment0) ~ "4weeks",
    grepl("2", Treatment0) ~ "8weeks",
    grepl("3", Treatment0) ~ "12weeks")) %>%
  #filter(!grepl("SRFA1", SampleID)) %>%
  select(-DRIP) %>%
  mutate(ID = str_remove(ID, "R")) %>%
  mutate(ID = str_remove(ID, "D")) %>%
  mutate(ID = str_remove(ID, "B")) %>%
  mutate(ID = str_remove(ID, "N")) %>%
  mutate(ID = paste("DRIP", ID, sep = "_")) %>%
  mutate(ID = case_when(
    Zone == "Bulk" ~ paste(ID, "b", sep = ""),
    TRUE ~ ID)) %>%
  unite(SampleID_new, c(ID, Treatment, Zone, Timepoint), sep = ".", remove = FALSE) %>%
  mutate(SampleType = case_when(
    grepl("ExBlk", SampleID) ~ "Blank",
    TRUE ~ "Sample"
  ))



edata <- edata.0[colnames(edata.0)[colnames(edata.0) %in% fdata$SampleID]]
edata <- cbind(edata.0[,1], edata)

#save meta-omics file with same names as other datasets

edata.t <- edata %>%
  column_to_rownames("Name") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%     # do this once
  left_join(fdata, by = "SampleID") %>%  # use dplyr join (keeps order, clearer)
  select(-c(ID, Treatment0, Zone, Treatment, Timepoint, SampleID, SampleType))  %>%# drop what you don’t need
  column_to_rownames(var="SampleID_new")

edata.neg <- edata.t %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var= "Name")


#remove features that are not at least 3x greater than median of blank
h_neg_clean <- blank_fold_filter_raw(
  edata.neg, 
  fdata,  
  edata_cname = "Name",
  fdata_cname = "SampleID_new",
  blank_label = "Blank",
  sample_type_cname = "SampleType",
  fold_threshold = 3,
  remove_blanks = TRUE
)
#no compound removed

#merge pos and neg and keep max of intensity for duplicates
h_clean_edata <- rbind(h_pos_clean$edata, h_neg_clean$edata)
h_clean_fdata <- rbind(h_pos_clean$fdata, h_neg_clean$fdata)
h_clean_emeta <- rbind(emeta.pos, emeta.neg)

h_clean_edata <- h_clean_edata %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE)

h_clean_emeta <- h_clean_emeta %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE)





####create omics data object
h_object <- as.metabData(
  e_data = h_clean_edata,
  f_data = h_clean_fdata,
  e_meta = h_clean_emeta,
  edata_cname = "Name",
  fdata_cname = "SampleID_new",
  emeta_cname = "Name",
  data_scale = "abundance",
  data_types = "Positive Ion"
)

class(h_object)
summary(h_object)


#plot transformed data
plot(edata_transform(h_object, data_scale = "log2"))

####filter biomolecules
mymolfilt <- molecule_filter(omicsData = h_object)
class(mymolfilt)
summary(mymolfilt) #all metabolites are found in all samples?

# coefficient of variation filter
h_object_groups <- group_designation(omicsData = h_object, main_effects = c("Timepoint",
                                                                            "Zone"))

mycvfilt <- cv_filter(omicsData = h_object_groups)
plot(mycvfilt, cv_threshold = 97)

h_object_groups_log2 <- edata_transform(h_object_groups, data_scale = "log2")
mycvfilt_log2 <- cv_filter(omicsData = h_object_groups_log2)
plot(mycvfilt_log2, cv_threshold = 97)

summary(h_object_groups_log2, cv_threshold = 97)
# Number Filtered Biomolecules: 21 


#apply the filter, this removes 8 biomolecules
h_object_groups_log2 <- applyFilt(filter_object = mycvfilt_log2, omicsData = h_object_groups_log2, cv_threshold = 97)

summary(h_object_groups_log2)

#imd-anova filter
myimdanovafilt <- imdanova_filter(omicsData = h_object_groups_log2)
summary(myimdanovafilt, min_nonmiss_anova = 2, min_nonmiss_gtest = 3) #no biomolecules filtered

####filter samples
myfilter <- rmd_filter(omicsData = h_object_groups_log2, metrics = c("Correlation", "MAD", "Skewness"))
plot(myfilter)

summary(myfilter, pvalue_threshold = 0.00000001)
#Filtered Samples: DRIP_15b.Drought.Bulk.4weeks, DRIP_165.Drought.RhizoDetritus.12weeks, DRIP_51.Drought.Rhizo.4weeks 

plot(myfilter, pvalue_threshold = 0.00000001, bw_theme = TRUE)

# get vector of potential outliers
potential_outliers <- summary(myfilter, pvalue_threshold = 0.00000001)$filtered_samples
# Filtered Samples: DRIP_51.Drought.Rhizo.4weeks 


# loop over potential outliers to generate plots
if (length(potential_outliers) > 0) {
  for (i in 1:length(potential_outliers)) {
    print(plot(myfilter, sampleID = potential_outliers[i]))
  }
}

mycor <- cor_result(omicsData = h_object_groups_log2)
plot(mycor, interactive = TRUE)

mypca <- dim_reduction(omicsData = h_object_groups_log2)

plot(mypca, interactive = TRUE)

#remove outliers found using metric above
h_object_groups_log2 <- applyFilt(filter_object = myfilter, 
                                  omicsData = h_object_groups_log2, pvalue_threshold = 0.00000001)
summary(h_object_groups_log2)

#numeric summary
edata_summary(h_object_groups_log2, by = "molecule", groupvar = NULL)

####normalization

# global median centering - apply the norm
h_object_groups_log2_norm <- normalize_global(
  omicsData = h_object_groups_log2,
  subset_fn = "all",
  norm_fn = "mean",
  apply_norm = TRUE,
  backtransform = TRUE
)
class(h_object_groups_log2_norm)

plot(h_object_groups_log2_norm)

####statstical analysis
#anova

h_object_groups_log2_norm_object <- group_designation(
  omicsData = h_object_groups_log2_norm,
  main_effects = c("Zone")
)

imd_filt <- imdanova_filter(omicsData = h_object_groups_log2_norm_object)
h_object_groups_log2_norm_object <- applyFilt(
  filter_object     = imd_filt,
  omicsData         = h_object_groups_log2_norm_object,
  min_nonmiss_anova = 2
)

# Step 2 — then run imd_anova
h_object_spans <- spans_procedure(h_object_groups_log2_norm_object)

h_object_pairwise <- imd_anova(
  omicsData              = h_object_groups_log2_norm_object,
  test_method            = "anova",
  pval_adjust_a_multcomp = "Tukey"
)

plot(h_object_pairwise, plot_type = "volcano")


###view PCA of normalized data
mycor <- cor_result(omicsData = h_object_groups_log2_norm_object)

mypca <- dim_reduction(omicsData = h_object_groups_log2_norm_object)

plot(mypca, interactive = TRUE)

#sig diff compounds (pairwise comparison)
p.cutoff = 0.2
rhizo.vs.bulk.up <- data.frame(compounds=h_object_pairwise$Name,
                               l2fc=h_object_pairwise$Fold_change_Rhizo_vs_Bulk, 
                               p.value=h_object_pairwise$P_value_A_Rhizo_vs_Bulk) %>%
  filter(p.value < p.cutoff &
           l2fc > 0)

rhizo.vs.bulk.down <- data.frame(compounds=h_object_pairwise$Name,
                                 l2fc=h_object_pairwise$Fold_change_Rhizo_vs_Bulk, 
                                 p.value=h_object_pairwise$P_value_A_Rhizo_vs_Bulk) %>%
  filter(p.value < p.cutoff &
           l2fc < 0)


rhizo.vs.rhizodetritus.up <- data.frame(compounds=h_object_pairwise$Name,
                                        l2fc=h_object_pairwise$Fold_change_Rhizo_vs_RhizoDetritus, 
                                        p.value=h_object_pairwise$P_value_A_Rhizo_vs_RhizoDetritus) %>%
  filter(p.value < p.cutoff &
  l2fc > 0)


rhizo.vs.rhizodetritus.down <- data.frame(compounds=h_object_pairwise$Name,
                                          l2fc=h_object_pairwise$Fold_change_Rhizo_vs_RhizoDetritus, 
                                          p.value=h_object_pairwise$P_value_A_Rhizo_vs_RhizoDetritus) %>%
  filter(p.value < p.cutoff &
           l2fc < 0)



rhizo.vs.detritus.up <- data.frame(compounds=h_object_pairwise$Name,
                                   l2fc=h_object_pairwise$Fold_change_Rhizo_vs_Detritus, 
                                   p.value=h_object_pairwise$P_value_A_Rhizo_vs_Detritus) %>%
  filter(p.value < p.cutoff &
  l2fc > 0)


rhizo.vs.detritus.down <- data.frame(compounds=h_object_pairwise$Name,
                                     l2fc=h_object_pairwise$Fold_change_Rhizo_vs_Detritus, 
                                     p.value=h_object_pairwise$P_value_A_Rhizo_vs_Detritus) %>%
  filter(p.value < p.cutoff &
           l2fc < 0)



rhizodet.vs.detritus.up <- data.frame(compounds=h_object_pairwise$Name,
                                      l2fc=h_object_pairwise$Fold_change_RhizoDetritus_vs_Detritus, 
                                      p.value=h_object_pairwise$P_value_A_RhizoDetritus_vs_Detritus) %>%
  filter(p.value < p.cutoff &
  l2fc > 0)


rhizodet.vs.detritus.down <- data.frame(compounds=h_object_pairwise$Name,
                                        l2fc=h_object_pairwise$Fold_change_RhizoDetritus_vs_Detritus, 
                                        p.value=h_object_pairwise$P_value_A_RhizoDetritus_vs_Detritus) %>%
  filter(p.value < p.cutoff &
           l2fc < 0)



rhizodet.vs.bulk.up <- data.frame(compounds=h_object_pairwise$Name,
                                  l2fc=h_object_pairwise$Fold_change_RhizoDetritus_vs_Bulk, 
                                  p.value=h_object_pairwise$P_value_A_RhizoDetritus_vs_Bulk) %>%
  filter(p.value < p.cutoff &
  l2fc > 0)


rhizodet.vs.bulk.down <- data.frame(compounds=h_object_pairwise$Name,
                                    l2fc=h_object_pairwise$Fold_change_RhizoDetritus_vs_Bulk, 
                                    p.value=h_object_pairwise$P_value_A_RhizoDetritus_vs_Bulk) %>%
  filter(p.value < p.cutoff &
           l2fc < 0)




bulk.vs.det.up <- data.frame(compounds=h_object_pairwise$Name,
                             l2fc=h_object_pairwise$Fold_change_Bulk_vs_Detritus, 
                             p.value=h_object_pairwise$P_value_A_Bulk_vs_Detritus) %>%
  filter(p.value < p.cutoff &
  l2fc > 0)


bulk.vs.det.down <- data.frame(compounds=h_object_pairwise$Name,
                               l2fc=h_object_pairwise$Fold_change_Bulk_vs_Detritus, 
                               p.value=h_object_pairwise$P_value_A_Bulk_vs_Detritus) %>%
  filter(p.value < p.cutoff &
           l2fc < 0)




#sig compounds to filter heatmap
sig.com <- rbind(rhizo.vs.bulk.up, rhizo.vs.bulk.down, rhizo.vs.rhizodetritus.up, 
                 rhizo.vs.rhizodetritus.down, rhizo.vs.detritus.up, rhizo.vs.detritus.down,
                 rhizodet.vs.bulk.up, rhizo.vs.bulk.down, bulk.vs.det.up, bulk.vs.det.down)

sig.com.unique <- sig.com %>%
  distinct(compounds)

##extract data and 
df.0 <- as.data.frame(h_object_groups_log2_norm_object$e_data)
df <- df.0 %>%
  rownames_to_column(var = "row") %>%
  column_to_rownames(var = "Name") %>%
  select(-row)

df.sig <- df %>%
  filter(rownames(.) %in% sig.com.unique$compounds)

df.z <- as.data.frame(t(scale(t(df))))


#save normalized and z-scaled dataframes
write.csv(df, "./output/lcms_lipids_norm.csv")

#load csv
df <- read.csv("./output/lcms_lipids_norm.csv", header = TRUE) %>%
  column_to_rownames(var = "X")



# Replace NAs with the minimum value in that feature's row
# (assumes below detection limit rather than truly absent)
df <- t(apply(df, 1, function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE)
  x
}))
df <- as.data.frame(df)

# Verify fixed
sum(is.na(df))  # should return 0
sum(!is.finite(as.matrix(df)))  # should return 0



# ══════════════════════════════════════════════════════════════════════════════
# SAMPLE METADATA
# ══════════════════════════════════════════════════════════════════════════════

meta <- fdata %>%
  column_to_rownames(var = "SampleID_new") %>%
  select(-Treatment0, -ID, -SampleID, -SampleType) %>%
  drop_na() %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Drought", "Untrt")),
    Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks"))
  )

# ══════════════════════════════════════════════════════════════════════════════
# Z-SCALE DATA (row-wise — each feature scaled across samples)
# ══════════════════════════════════════════════════════════════════════════════

df.z     <- as.data.frame(t(scale(t(df))))
df.z.sig <- as.data.frame(t(scale(t(df.sig))))

# ══════════════════════════════════════════════════════════════════════════════
# PRODUCE HEATMAPS
# ══════════════════════════════════════════════════════════════════════════════

outdir <- "./figures/Fig3_FigSx_LCMS"

# ── Compute global color scale across all datasets ────────────────────────────
# Do this BEFORE calling make_all_heatmaps
#filter meta
meta.f <- meta %>%
  filter(rownames(.) %in% colnames(df.z))


# compute once using all datasets (full + sig + all timepoints)
tp_samples_list <- lapply(c("4weeks", "8weeks", "12weeks"), function(tp) {
  rownames(meta.f)[meta.f$Timepoint == tp]
})


global_breaks <- compute_global_breaks(
  df.z, df.z.sig,
  df.z[, tp_samples_list[[1]]],
  df.z[, tp_samples_list[[2]]],
  df.z[, tp_samples_list[[3]]]
)

#global_breaks <- seq(-4, 4, length.out = 100)


# ── All timepoints ────────────────────────────────────────────────────────────
make_all_heatmaps(
  df_z     = df.z,
  df_z_sig = df.z.sig,
  meta     = meta.f,
  meta.comp = meta.comp,
  class_base_colors = class_base_colors,
  zone_colors            = zone_colors,
  treatment_colors       = treatment_colors,
  timepoint_colors       = timepoint_colors,
  outdir = outdir,
  #cluster_rows = TRUE,
  label  = "all-timepoints_",
  breaks = global_breaks
)

# ── Per timepoint ─────────────────────────────────────────────────────────────
timepoints <- c("4weeks", "8weeks", "12weeks")

# Define k values for each timepoint before the loop
k_values <- list("4weeks" = 6, "8weeks" = 4, "12weeks" = 4)  # adjust to match your visual clusters

for (tp in timepoints) {
  
  # samples belonging to this timepoint
  tp_samples <- rownames(meta)[meta$Timepoint == tp]
  tp_samples <- tp_samples[tp_samples %in% colnames(df.z)]
  
  if (length(tp_samples) == 0) {
    warning("No samples found for timepoint: ", tp)
    next
  }
  
  # subset data and metadata to this timepoint
  df_z_tp     <- df.z[,     tp_samples, drop = FALSE]
  df_z_sig_tp <- df.z.sig[, tp_samples, drop = FALSE]
  meta_tp     <- meta.f[tp_samples, , drop = FALSE]
  
  # drop Timepoint column from annotation since it's constant within a timepoint
  meta_tp <- meta_tp %>% select(-Timepoint)
  
  # drop timepoint from annotation colors
  tp_zone_colors      <- zone_colors
  tp_treatment_colors <- treatment_colors
  
  make_all_heatmaps(
    df_z     = df_z_tp,
    df_z_sig = df_z_sig_tp,
    meta     = meta_tp,
    meta.comp = meta.comp,
    class_base_colors = class_base_colors,
    zone_colors            = tp_zone_colors,
    treatment_colors       = tp_treatment_colors,
    timepoint_colors       = NULL,
    outdir = outdir,
    label  = paste0( tp, "_"),
    width = 15, height_all = 17, height_sig = 14,
    breaks = global_breaks
  )
  
  # Extract the hclust object from pheatmap
  heatmap_obj <- pheatmap(df_z_sig_tp, 
                          clustering_method = "complete",
                          silent = TRUE)
  
  # Get k value for this timepoint
  k <- k_values[[tp]]
  
  # Cut the row dendrogram using timepoint-specific k
  clusters <- cutree(heatmap_obj$tree_row, k = k)
  
  # See which lipids are in each cluster
  cluster_df <- data.frame(lipid = names(clusters), cluster = clusters) %>%
    arrange(cluster)
  
  # Add cluster membership back to your data
  df_z_sig_tp$cluster <- clusters
  
  df_z_sig_tp <- df_z_sig_tp %>%
    select(cluster)
  
  write.csv(df_z_sig_tp, paste0("./output/LCMS/",
                                tp, "_lipids_clusters.csv"))
}




#extract compounds from each cluster at each timepoint and associate with habitat
#4 weeks
week.4 <- read.csv("./output/LCMS/4weeks_lipids_clusters.csv")

rhizo.4 <- week.4 %>%
  filter(cluster==1 |
           cluster==2)
rhizo.4 <- rhizo.4$X

det.4 <- week.4 %>%
  filter(cluster==3 |
           cluster==4)
det.4 <- det.4$X

rhizodet.4 <- det.4

#8 weeks
week.8 <- read.csv("./output/LCMS/8weeks_lipids_clusters.csv")

rhizo.8 <- week.8 %>%
  filter(cluster==1)

rhizo.8 <- rhizo.8$X

det.8 <- week.8 %>%
  filter(cluster==2 |
           cluster==3 |
           cluster==4)
det.8 <- det.8$X

rhizodet.8 <- det.8

#12 weeks
week.12 <- read.csv("./output/LCMS/12weeks_lipids_clusters.csv")

rhizo.12 <- week.12 %>%
  filter(cluster==1 |
           cluster==2)

rhizo.12 <- rhizo.12$X

det.12 <- week.12 %>%
  filter(cluster==3 |
           cluster==4)
det.12 <- det.12$X

rhizodet.12 <- rhizo.12



#filter dataframe to above clusters
#12 weeks
to_plot_rhizo.12 <- meta.comp %>%
  rownames_to_column(var = "Compound") %>%
  filter(Compound %in% rhizo.12) %>%
  distinct() %>%
  mutate(habitat = "Rhizo") %>%
  mutate(Timepoint = "12weeks")

to_plot_det.12 <- meta.comp %>%
  rownames_to_column(var = "Compound") %>%
  filter(Compound %in% det.12) %>%
  distinct() %>%
  mutate(habitat = "Detritus") %>%
  mutate(Timepoint = "12weeks")

to_plot_rhizodet.12 <- to_plot_rhizo.12 %>%
  mutate(habitat = "RhizoDet")

#8 weeks
to_plot_rhizo.8 <- meta.comp %>%
  rownames_to_column(var = "Compound") %>%
  filter(Compound %in% rhizo.8) %>%
  distinct() %>%
  mutate(habitat = "Rhizo") %>%
  mutate(Timepoint = "8weeks")

to_plot_det.8 <- meta.comp %>%
  rownames_to_column(var = "Compound") %>%
  filter(Compound %in% det.8) %>%
  distinct() %>%
  mutate(habitat = "Detritus") %>%
  mutate(Timepoint = "8weeks")

to_plot_rhizodet.8 <- to_plot_det.8 %>%
  mutate(habitat = "RhizoDet")

#4 weeks
to_plot_rhizo.4 <- meta.comp %>%
  rownames_to_column(var = "Compound") %>%
  filter(Compound %in% rhizo.4) %>%
  distinct() %>%
  mutate(habitat = "Rhizo") %>%
  mutate(Timepoint = "4weeks")

to_plot_det.4 <- meta.comp %>%
  rownames_to_column(var = "Compound") %>%
  filter(Compound %in% det.4) %>%
  distinct() %>%
  mutate(habitat = "Detritus") %>%
  mutate(Timepoint = "4weeks")

to_plot_rhizodet.4 <- to_plot_det.4 %>%
  mutate(habitat = "RhizoDet")



#plot compound classes
#12 weeks

to_plot_combined <- rbind(to_plot_rhizo.12, to_plot_rhizodet.12, to_plot_det.12,
                           to_plot_rhizo.8, to_plot_rhizodet.8, to_plot_det.8,
                           to_plot_rhizo.4, to_plot_rhizodet.4, to_plot_det.4)
to_plot <- to_plot_combined %>%
  mutate(Class = factor(Class, levels = c("Diacylglycerol", "Triacylglycerol",
                                          "Phosphatidylcholine", "Phosphatidylethanolamine",
                                          "Phosphatidylglycerol", "Ceramide"))) %>%
  mutate(habitat = factor(habitat, levels = c("Rhizo", "RhizoDet", "Detritus"))) %>%
  group_by(Class, habitat, Timepoint) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks"))) %>%
  summarise(count = n()) %>%
  ggplot(aes(x=Timepoint, y=count, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~habitat, ncol = 3) +
  scale_fill_manual(values = class_base_colors) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("./figures/Fig3_FigSx_LCMS/FigSx_lipids_class_time.png", dpi = 300,
      height = 3, width = 5, units = "in")

ggsave("./figures/Fig3_FigSx_LCMS/FigSx_lipids_class_time.pdf", dpi = 300,
       height = 3, width = 5, units = "in")





##nmds
# Calculate nmds
#first transpose
df.t <- df %>%
  t() %>%
  as.data.frame() 

set.seed(123)
nmds = metaMDS(df.t, distance = "euclidean", autotransform = FALSE)
plot(nmds)

scores <- scores(nmds)
scores.samples <- scores$site

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  tidyr::separate_wider_delim(SampleID,
                              delim = ".",
                              names = c("SampleID", "Treatment", "Zone", "Time"),
                              cols_remove = FALSE) %>%
  mutate(Treatment = str_replace(Treatment, "Untrt", "Control")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Drought", "Control"))) %>%
  mutate(Zone = factor(Zone, levels = c("Rhizo", "RhizoDetritus", "Detritus", "Bulk")))

# envfit with aligned metadata
ef <- envfit(nmds ~ Zone + Treatment + Time, data = data.scores, permutations = 999)
ef

### nmds plot (zone and treatment)
w=12
h=8
nmds_plot <-  make_nmds_plot(data.scores, Zone, Treatment, col_list = zone_colors)
nmds_plot
filename <- paste0("./figures/Fig3_FigSx_LCMS/FigSx-lipids-nmds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0("./figures/Fig3_FigSx_LCMS/FigSx-lipids-nmds.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)



### nmds plot (treatment and zone)
nmds_plot <-  make_nmds_plot(data.scores, Time, Treatment, col_list = timepoint_colors)
nmds_plot
filename <- paste0("./figures/Fig3_FigSx_LCMS/other/lipids-nmds-timepoint-treatment.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

#permanova
distmat <- vegdist(t(as.data.frame(df)), method = "euclidean")


adonis2(distmat ~ Zone+Treatment+Time, data.scores, permutations = 999, by = "margin")
# adonis2(formula = distmat ~ Zone + Treatment + Time, data = data.scores, permutations = 999, by = "margin")
# Df SumOfSqs      R2       F Pr(>F)    
# Zone       3   1922.9 0.32529 15.7739  0.001 ***
#   Treatment  1    123.3 0.02085  3.0334  0.009 ** 
#   Time       2    296.6 0.05018  3.6497  0.002 ** 
#   Residual  88   3575.8 0.60491                   
# Total     94   5911.3 1.00000      


