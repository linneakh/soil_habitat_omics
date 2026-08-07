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

# ── Color definitions ─────────────────────────────────────────────────────────
superclass_base_colors <- c(
  "Lipids"                       = "#A052B8",  # lighter deep purple
  "Fatty Acyls"                  = "#F0C060",  # lighter golden amber
  "Organic acids"                = "#56A8D0",  # lighter steel blue
  "Organic nitrogen compounds"   = "#E8707A",  # lighter raspberry
  "Nucleic acids"                = "#7AB050",  # lighter olive green
  "Benzenoids"                   = "#FF9A40",  # lighter deep orange
  "Carbohydrates"                = "#40E0D0",  # lighter turquoise
  "Organic oxygen compounds"     = "#7A9AAA",  # lighter slate
  "Alkaloids"                    = "#F8EE70",  # lighter bright yellow
  "Polyketides"                  = "#D44080",  # lighter deep pink/magenta
  "Organoheterocyclic compounds" = "#C8FF60"   # lighter chartreuse
)

zone_colors      <- c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red")
treatment_colors <- c(Drought = "brown", Untrt = "darkgreen")
timepoint_colors <- c(`4weeks` = "#F4A582", `8weeks` = "#74ADD1", `12weeks` = "#D9EF8B")

class_hierarchy <- data.frame(
  superclass = c(
    "Lipids",                      "Lipids",           "Lipids",
    "Lipids",                      "Lipids",
    "Fatty Acyls",                 "Fatty Acyls",      "Fatty Acyls",      "Fatty Acyls",
    "Organic acids",               "Organic acids",    "Organic acids",
    "Organic nitrogen compounds",  "Organic nitrogen compounds", "Organic nitrogen compounds",
    "Nucleic acids",               "Nucleic acids",
    "Benzenoids",                  "Benzenoids",       "Benzenoids",
    "Alkaloids",                   "Alkaloids",
    "Polyketides",                 "Polyketides",      "Polyketides",
    "Carbohydrates",               "Carbohydrates",    "Carbohydrates", "Carbohydrates",
    "Organic oxygen compounds",    "Organic oxygen compounds",
    "Organoheterocyclic compounds"
  ),
  class = c(
    "Glycerophospholipids",        "Prenol lipids",    "Sphingolipids",
    "Sterol lipids",               "Glycerolipids",
    "Fatty acids",                 "Fatty esters",     "Fatty amides",     "Octadecanoids",
    "Amino acids and peptides",    "Carboxylic acids and derivatives", "Keto acids",
    "Amides",                      "Carnitines",       "Cholines",
    "Purines",                     "Pyrimidines",
    "Benzophenones",               "Benzamides",       "Anilines",
    "Tryptophan alkaloids",        "Nicotinic acid alkaloids",
    "Flavonoids",                  "Phenolic acids",   "Phenylpropanoids",
    "Monosaccharides",             "Disaccharides",    "Oligosaccharides", "Trisaccarides",
    "Organooxygen compounds",      "Alcohols and polyols",
    "Indolyl carboxylic acids"
  ),
  stringsAsFactors = FALSE
)


##load LC-MS data
file_path <- "./data/manuscript/DRIPSIP_LC_Metab_WK_CleanedUp.xlsx"
#kegg <- read.csv("./qSIP_output/LC-MS/HPOS_metaboanalyst_KEGGID_all.csv", header = TRUE)
#kegg$Query <- gsub(" ", "_", kegg$Query) 


##load extra functions
source("./scripts/manuscript/functions/extra_functions_PCA_NMDS_plotting.R")
source("./scripts/manuscript/functions/lc-ms-functions.R")

# Get sheet names
sheet_names <- excel_sheets(file_path)
#"Groups"         "HILIC Positive" "HILIC Negative" "RP Positive"    "RP Negative"   

# Read each sheet into a data frame and store in a list
data_frames_list <- lapply(sheet_names, function(sheet_name) {
  read_excel(file_path, sheet = sheet_name, range = NULL)
})

# Extract the first data frame
peak_table_pos.0 = data_frames_list[[2]]

# Extract the second data frame
peak_table_neg.0 <- data_frames_list[[3]]

# Extract the third data frame
peak_table_rpos.0 <- data_frames_list[[4]]

# Extract the fourth data fram
peak_table_rneg.0 <- data_frames_list[[5]]

# ── Load and clean compound classifications ───────────────────────────────────
meta.comp <- rbind(
  peak_table_pos.0 %>% select(Name, `Super class`, `Main class`, `Sub class`, Formula),
  peak_table_neg.0 %>% select(Name, `Super class`, `Main class`, `Sub class`, Formula)
) %>%
  unique() %>%
  column_to_rownames(var = "Name") %>%
  mutate(
    `Super class` = case_when(
      `Super class` == "Glycerophospholipids" ~ "Lipids",
      `Super class` == "Sphingolipids"        ~ "Lipids",
      TRUE ~ `Super class`
    ),
    `Main class` = case_when(
      `Main class` == "Prenol Lipids"                    ~ "Prenol lipids",
      `Main class` == "Sterol Lipids"                    ~ "Sterol lipids",
      `Main class` == "Carboxylic acids and derivitives" ~ "Carboxylic acids and derivatives",
      `Main class` == "Glycerophophocholines"            ~ "Glycerophosphocholines",
      TRUE ~ `Main class`
    ),
    `Sub class` = case_when(
      `Sub class` == "Sesquiterpeoids"       ~ "Sesquiterpenoids",
      `Sub class` == "Glycerophophocholines" ~ "Glycerophosphocholines",
      `Sub class` == "Oxolipins"             ~ "Oxylipins",
      TRUE ~ `Sub class`
    )
  )



##HPOS
#reformat and filter pos table
peak_table_pos.0$Name <- gsub(" ", "_", peak_table_pos.0$Name) 
colnames(peak_table_pos.0)[1] <- "mz"
colnames(peak_table_pos.0)[2] <- "rt"
colnames(peak_table_pos.0) = gsub("_pos", "", colnames(peak_table_pos.0))

peak_table_pos <- peak_table_pos.0 %>%
  filter(Tags == "Confirmed ID (HIgh Confidence)" |
           Tags == "Match with Mass and RT only" |
           Tags == "Not in Internal DB. MS/MS match only" |
           Tags == "Not in Internal DB. Potential ID (Low Confidence)" |
           Tags == "Potential ID (Low Confidence)")

#order column names to match neg table
peak_table_pos2 <-  peak_table_pos[,-c(4:12)] %>%
  unite(mz_name, c("mz", "Name"), sep = "_", remove = FALSE) %>%
  column_to_rownames(var = "mz_name") %>%
  select(order(colnames(.))) %>%
  rownames_to_column(var = "mz_name") 

colnames(peak_table_pos2) <- str_remove(colnames(peak_table_pos2), "_HPOS")


compound.class.pos <- peak_table_pos[,c(3, 5:10)]



edata.0 <- peak_table_pos2 %>%
  select(-mz, -rt) %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE) 


emeta.pos <-  peak_table_pos2 %>%
  select(Name, mz, rt) %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE) 

fdata.0 <- edata.0 %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE)  %>%
  column_to_rownames(var = "Name") 


fdata <- as.data.frame(t(fdata.0)) %>%
  rownames_to_column(var = "SampleID") %>%
  dplyr::select (SampleID) %>%
  tidyr::separate_wider_delim(SampleID,
                              delim = "_",
                              names = c("ID", "Treatment0"),
                              cols_remove = FALSE,
                              too_few = "align_start") %>%
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
  filter(!grepl("SRFA1", SampleID)) %>%
  #select(-DRIP) %>%
  mutate(ID = str_remove(ID, "R")) %>%
  mutate(ID = str_remove(ID, "D")) %>%
  mutate(ID = str_remove(ID, "B")) %>%
  mutate(ID = str_remove(ID, "N")) %>%
  mutate(ID = paste("DRIP", ID, sep = "_")) %>%
  mutate(ID = case_when(
    Zone == "Bulk" ~ paste(ID, "b", sep = ""),
    TRUE ~ ID)) %>%
  filter(SampleID != "SRFA1") %>%
  unite(SampleID_new, c(ID, Treatment, Zone, Timepoint), sep = ".", remove = FALSE) %>%
  mutate(SampleType = case_when(
    grepl("ExBlk", SampleID) ~ "Blank",
    TRUE ~ "Sample"
  ))

edata <- edata.0[colnames(edata.0)[colnames(edata.0) %in% fdata$SampleID]]
edata$Name <- edata.0$Name

#save meta-omics file with same names as other datasets

edata.t <- edata %>%
  column_to_rownames(var = "Name") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%     # do this once
  left_join(fdata, by = "SampleID") %>%  # use dplyr join (keeps order, clearer)
  select(-c(ID, Treatment0, Zone, Treatment, Timepoint, SampleID, SampleType))  %>%# drop what you don’t need
  column_to_rownames(var="SampleID_new")

edata.final <- edata.t %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var= "Name") 

#remove features that are not at least 3x greater than median of blank
h_pos_clean <- blank_fold_filter_raw(
  edata.final, 
  fdata,  
  edata_cname = "Name",
  fdata_cname = "SampleID_new",
  blank_label = "Blank",
  sample_type_cname = "SampleType",
  fold_threshold = 3,
  remove_blanks = TRUE
)

# Found 12 blank sample(s) and 97 non-blank sample(s).
# Fold-change filter (threshold = 3x): keeping 35 / 46 features, removing 11.
# Removed 12 blank sample(s).

#H NEG
#reformat and filter pos table
peak_table_neg.0$Name <- gsub(" ", "_", peak_table_neg.0$Name) 
colnames(peak_table_neg.0)[1] <- "mz"
colnames(peak_table_neg.0)[2] <- "rt"
colnames(peak_table_neg.0) = gsub("_HILIC", "", colnames(peak_table_neg.0))
colnames(peak_table_neg.0) = gsub("_neg", "", colnames(peak_table_neg.0))

peak_table_neg <- peak_table_neg.0 %>%
  filter(Tags == "Confirmed ID (HIgh Confidence)" |
           Tags == "Match with Mass and RT only" |
           Tags == "Not in Internal DB. MS/MS match only" |
           Tags == "Not in Internal DB. Potential ID (Low Confidence)" |
           Tags == "Potential ID (Low Confidence)")

#order column names to match neg table
peak_table_neg2 <-  peak_table_neg[,-c(4:12)] %>%
  unite(mz_name, c("mz", "Name"), sep = "_", remove = FALSE) %>%
  column_to_rownames(var = "mz_name") %>%
  select(order(colnames(.))) %>%
  rownames_to_column(var = "mz_name") 

colnames(peak_table_neg2) <- str_remove(colnames(peak_table_neg2), "_HILIC")


compound.class.neg <- peak_table_neg[,c(3, 5:10)]



edata.0 <- peak_table_neg2 %>%
  select(-mz, -rt) %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE) 


emeta.neg <-  peak_table_neg2 %>%
  select(Name, mz, rt) %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE) 

fdata.0 <- edata.0 %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE)  %>%
  column_to_rownames(var = "Name") 


fdata <- as.data.frame(t(fdata.0)) %>%
  rownames_to_column(var = "SampleID") %>%
  dplyr::select (SampleID) %>%
  separate_wider_delim(SampleID, 
                       delim = "_",
                       names = c("ID", "Treatment0"), 
                       cols_remove = FALSE,
                       too_few = "align_start") %>%
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
  filter(!grepl("SRFA1", SampleID)) %>%
  #select(-DRIP) %>%
  mutate(ID = str_remove(ID, "R")) %>%
  mutate(ID = str_remove(ID, "D")) %>%
  mutate(ID = str_remove(ID, "B")) %>%
  mutate(ID = str_remove(ID, "N")) %>%
  mutate(ID = paste("DRIP", ID, sep = "_")) %>%
  mutate(ID = case_when(
    Zone == "Bulk" ~ paste(ID, "b", sep = ""),
    TRUE ~ ID)) %>%
  filter(SampleID != "SRFA1") %>%
  unite(SampleID_new, c(ID, Treatment, Zone, Timepoint), sep = ".", remove = FALSE) %>%
  mutate(SampleType = case_when(
    grepl("ExBlk", SampleID) ~ "Blank",
    TRUE ~ "Sample"
  ))

edata <- edata.0[colnames(edata.0)[colnames(edata.0) %in% fdata$SampleID]]
edata$Name <- edata.0$Name

#save meta-omics file with same names as other datasets

edata.t <- edata %>%
  column_to_rownames(var = "Name") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%     # do this once
  left_join(fdata, by = "SampleID") %>%  # use dplyr join (keeps order, clearer)
  select(-c(ID, Treatment0, Zone, Treatment, Timepoint, SampleID, SampleType))  %>%# drop what you don’t need
  column_to_rownames(var="SampleID_new")

edata.final <- edata.t %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var= "Name") 


#remove features that are not at least 3x greater than median of blank
h_neg_clean <- blank_fold_filter_raw(
  edata.final, 
  fdata,  
  edata_cname = "Name",
  fdata_cname = "SampleID_new",
  blank_label = "Blank",
  sample_type_cname = "SampleType",
  fold_threshold = 3,
  remove_blanks = TRUE
)

# Found 12 blank sample(s) and 97 non-blank sample(s).
# Fold-change filter (threshold = 3x): keeping 21 / 29 features, removing 8.
# Removed 12 blank sample(s).

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
myfilter <- rmd_filter(omicsData = h_object_groups_log2, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness"))
plot(myfilter)

summary(myfilter, pvalue_threshold = 0.00000001)

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

#detect most different compounds for plotting (pairwise comparison)
rhizo.vs.bulk.up <- data.frame(compounds=h_object_pairwise$Name,
                               l2fc=h_object_pairwise$Fold_change_Rhizo_vs_Bulk, 
                               p.value=h_object_pairwise$P_value_A_Rhizo_vs_Bulk) %>%
  filter(p.value < 0.2 &
           l2fc > 0)

rhizo.vs.bulk.down <- data.frame(compounds=h_object_pairwise$Name,
                                 l2fc=h_object_pairwise$Fold_change_Rhizo_vs_Bulk, 
                                 p.value=h_object_pairwise$P_value_A_Rhizo_vs_Bulk) %>%
  filter(p.value < 0.2 &
           l2fc < 0)


rhizo.vs.rhizodetritus.up <- data.frame(compounds=h_object_pairwise$Name,
                                        l2fc=h_object_pairwise$Fold_change_Rhizo_vs_RhizoDetritus, 
                                        p.value=h_object_pairwise$P_value_A_Rhizo_vs_RhizoDetritus) %>%
  filter(p.value < 0.2 &
  l2fc > 0)


rhizo.vs.rhizodetritus.down <- data.frame(compounds=h_object_pairwise$Name,
                                          l2fc=h_object_pairwise$Fold_change_Rhizo_vs_RhizoDetritus, 
                                          p.value=h_object_pairwise$P_value_A_Rhizo_vs_RhizoDetritus) %>%
  filter(p.value < 0.2 &
           l2fc < 0)



rhizo.vs.detritus.up <- data.frame(compounds=h_object_pairwise$Name,
                                   l2fc=h_object_pairwise$Fold_change_Rhizo_vs_Detritus, 
                                   p.value=h_object_pairwise$P_value_A_Rhizo_vs_Detritus) %>%
  filter(p.value < 0.2 &
  l2fc > 0)


rhizo.vs.detritus.down <- data.frame(compounds=h_object_pairwise$Name,
                                     l2fc=h_object_pairwise$Fold_change_Rhizo_vs_Detritus, 
                                     p.value=h_object_pairwise$P_value_A_Rhizo_vs_Detritus) %>%
  filter(p.value < 0.2 &
           l2fc < 0)



rhizodet.vs.detritus.up <- data.frame(compounds=h_object_pairwise$Name,
                                      l2fc=h_object_pairwise$Fold_change_RhizoDetritus_vs_Detritus, 
                                      p.value=h_object_pairwise$P_value_A_RhizoDetritus_vs_Detritus) %>%
  filter(p.value < 0.2 &
  l2fc > 0)


rhizodet.vs.detritus.down <- data.frame(compounds=h_object_pairwise$Name,
                                        l2fc=h_object_pairwise$Fold_change_RhizoDetritus_vs_Detritus, 
                                        p.value=h_object_pairwise$P_value_A_RhizoDetritus_vs_Detritus) %>%
  filter(p.value < 0.2 &
           l2fc < 0)



rhizodet.vs.bulk.up <- data.frame(compounds=h_object_pairwise$Name,
                                  l2fc=h_object_pairwise$Fold_change_RhizoDetritus_vs_Bulk, 
                                  p.value=h_object_pairwise$P_value_A_RhizoDetritus_vs_Bulk) %>%
  filter(p.value < 0.2 &
  l2fc > 0)


rhizodet.vs.bulk.down <- data.frame(compounds=h_object_pairwise$Name,
                                    l2fc=h_object_pairwise$Fold_change_RhizoDetritus_vs_Bulk, 
                                    p.value=h_object_pairwise$P_value_A_RhizoDetritus_vs_Bulk) %>%
  filter(p.value < 0.2 &
           l2fc < 0)




bulk.vs.det.up <- data.frame(compounds=h_object_pairwise$Name,
                             l2fc=h_object_pairwise$Fold_change_Bulk_vs_Detritus, 
                             p.value=h_object_pairwise$P_value_A_Bulk_vs_Detritus) %>%
  filter(p.value < 0.2 &
  l2fc > 0)


bulk.vs.det.down <- data.frame(compounds=h_object_pairwise$Name,
                               l2fc=h_object_pairwise$Fold_change_Bulk_vs_Detritus, 
                               p.value=h_object_pairwise$P_value_A_Bulk_vs_Detritus) %>%
  filter(p.value < 0.2 &
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



#save normalized and z-scaled dataframes
# write.csv(df, "./output/manuscript/lcms_met_h_norm.csv")
# write.csv(df.z, "./output/manuscript/lcms_met_h_norm_z_score.csv")

#load csv
df <- read.csv("./output/manuscript/lcms_met_h_norm.csv", header = TRUE) %>%
  column_to_rownames(var = "X")



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

outdir <- "./manuscript/"

# ── Compute global color scale across all datasets ────────────────────────────
# Do this BEFORE calling make_all_heatmaps


# compute once using all datasets (full + sig + all timepoints)
tp_samples_list <- lapply(c("4weeks", "8weeks", "12weeks"), function(tp) {
  rownames(meta)[meta$Timepoint == tp]
})

# global_breaks <- compute_global_breaks(
#   df.z, df.z.sig,
#   df.z[, tp_samples_list[[1]]],
#   df.z[, tp_samples_list[[2]]],
#   df.z[, tp_samples_list[[3]]]
# )

global_breaks <- seq(-4, 4, length.out = 100)


# ── All timepoints ────────────────────────────────────────────────────────────
make_all_heatmaps(
  df_z     = df.z,
  df_z_sig = df.z.sig,
  meta     = meta,
  meta.comp = meta.comp,
  class_hierarchy        = class_hierarchy,
  superclass_base_colors = superclass_base_colors,
  zone_colors            = zone_colors,
  treatment_colors       = treatment_colors,
  timepoint_colors       = timepoint_colors,
  outdir = outdir,
  #cluster_rows = TRUE,
  label  = "lcms-h-all-timepoints",
  breaks = global_breaks
)

# ── Per timepoint ─────────────────────────────────────────────────────────────
timepoints <- c("4weeks", "8weeks", "12weeks")

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
  meta_tp     <- meta[tp_samples, , drop = FALSE]
  
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
    class_hierarchy        = class_hierarchy,
    superclass_base_colors = superclass_base_colors,
    zone_colors            = tp_zone_colors,
    treatment_colors       = tp_treatment_colors,
    timepoint_colors       = NULL,   # not shown within a single timepoint
    outdir = outdir,
    label  = paste0("Fig.3-lcms-h-heatmaps", tp),
    width = 15, height_all = 13, height_sig = 8,
    breaks = global_breaks
  )
}


###extract kegg ids for compounds associated with each habitat:
kegg_ids <- rbind(compound.class.neg, compound.class.pos) 

#define metabolites associated with each time period and habitat
rhizo.12 <- c("Guanine", "Hypoxanthine", "3-Hydroxybutyricacid", 
           "Mevalolactone", "L-Phenylalanyl-L-proline", "Leucylproline",
           "Adenine", "DL-Isoleucine", "Butanoylcarnitine", "Pantolactone",
           "2,4-Bis(2-methylbutan-2-yl)phenol",
           "9R-hydroxy-10E,12E-octadecadienoicacid,methylester",
           "Linoleicacid", "Oleicacid", "Leucyl-Valine",
           "2-Ethylhexanoicacid", "Indole-3-aldehyde")

det.12 <- c("PE-Nme", "PE2", 
           "Stachydrine", "Sugars-Disaccharides", "Adenosine",
           "L-Citrulline", "Sugars-Trisaccharides", "Choline",
           "Carnitine", "PE1", "PE3", "Prolylglycine", "Uracil",
           "DL-Carnitine", "Nalpha-acetyl-L-Lysine", "Acetyl-DL-carnitine")

rhizo.det.12 <- rhizo.12

rhizo.8 <- c("Pantolactone", "Hypoxanthine", 
             "Mevalolactone", "L-Phenylalanyl-L-proline", "Leucylproline",
             "Adenine","Guanine", "DL-Isoleucine", "3-Hydroxybutyricacid",
             "Butanoylcarnitine", "Indole-3-aldehyde","LinoleicAcid",
             "Oleicacid", "9R-hydroxy-10E,12E-octadecadienoicacid,methylester")
             

det.8 <- c("Acetyl-DL-carnitine", "Stachydrine", "Sugars-Disaccharides",
           "Sugars-Trisaccharides", 
           "L-Citrilline", "Leucyl-Valine", "Prolyglycine", "PE1", "PE3",
           "PE-Nme", "PE2", "2-Ethylhexanoicacid", "Choline", "Adenosine",
           "Carnitine", "DL-Carnitine", "Nalpha-aceytl-L-Lysine")

rhizo.det.8 <- det.8

rhizo.4 <- c("Adenosine","L-Citrulline","Sugars-Disaccharides","Sugars-Trisaccharides", 
             "PE2","PE-Nme","PE1", "PE2", "Prolyglycine", "Pantolactone", 
             "9R-hydroxy-10E,12E-octadecadienoicacid,methylester", "Mevalolactone",
             "L-Phenylalanyl-L-proline", "Leucylproline",
             "Hypoxanthine", "Adenine")

det.4 <- c("Stachydrine", "3-Hydroxybutyricacid", "3-Ethylhexanoicacid")

rhizo.det.4 <- det.4
#filter kegg ids to each list above
#12 weeks
kegg_ids_rhizo.12 <- kegg_ids %>%
  filter(Name %in% rhizo.12) %>%
  distinct() %>%
  mutate(habitat = "Rhizo") %>%
  mutate(Timepoint = "12weeks")

kegg_ids_rhizo.12.no.na <- kegg_ids_rhizo.12 %>%
  filter(KEGG != "NA") 

kegg_ids_det.12 <- kegg_ids %>%
  filter(Name %in% det.12) %>%
  distinct() %>%
  mutate(habitat = "Detritus") %>%
  mutate(Timepoint = "12weeks")

kegg_ids_det.12.no.na <- kegg_ids_det.12 %>%
  filter(KEGG != "NA") 

kegg_ids_rhizo.det.12 <- kegg_ids_rhizo.12 %>%
  mutate(habitat = "RhizoDet")

#8 weeks
kegg_ids_rhizo.8 <- kegg_ids %>%
  filter(Name %in% rhizo.8) %>%
  distinct() %>%
  mutate(habitat = "Rhizo") %>%
  mutate(Timepoint = "8weeks")

kegg_ids_rhizo.8.no.na <- kegg_ids_rhizo.8 %>%
  filter(KEGG != "NA") 

kegg_ids_det.8 <- kegg_ids %>%
  filter(Name %in% det.8) %>%
  distinct() %>%
  mutate(habitat = "Detritus") %>%
  mutate(Timepoint = "8weeks")

kegg_ids_det.8.no.na <- kegg_ids_det.8 %>%
  filter(KEGG != "NA") 

kegg_ids_rhizo.det.8 <- kegg_ids_det.8 %>%
  mutate(habitat = "RhizoDet")

#4 weeks
kegg_ids_rhizo.4 <- kegg_ids %>%
  filter(Name %in% rhizo.4) %>%
  distinct() %>%
  mutate(habitat = "Rhizo") %>%
  mutate(Timepoint = "4weeks")

kegg_ids_rhizo.4.no.na <- kegg_ids_rhizo.4 %>%
  filter(KEGG != "NA") 

kegg_ids_det.4 <- kegg_ids %>%
  filter(Name %in% det.4) %>%
  distinct() %>%
  mutate(habitat = "Detritus") %>%
  mutate(Timepoint = "4weeks")

kegg_ids_det.4.no.na <- kegg_ids_det.4 %>%
  filter(KEGG != "NA") 

kegg_ids_rhizo.det.4 <- kegg_ids_det.4 %>%
  mutate(habitat = "RhizoDet")


#run clusster profiler to find enriched pathways
#12 weeks
k_rhizo.12 <- enrichKEGG(gene = kegg_ids_rhizo.12.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_rhizo.12 <- filter_and_plot(k_rhizo.12,
                               min_count = 2,
                               title = "Rhizosphere, 12 weeks")
#nucleotide metabolism (adenine, guanine, hypoxanthine[all purines]), 
#biosynthesis of plant secondary metabolites (L-Isoleucine, 9Z-octadecenoic acid),
#cutein suberine and wax biosynthesis (9Z-octadecenoic acid)


m_rhizo.12 <- enrichMKEGG(gene = kegg_ids_rhizo.12.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)

k_det.12 <- enrichKEGG(gene = kegg_ids_det.12.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_det.12 <- filter_and_plot(k_det.12,
                             remove_patterns = c("Neuroactive ligand signaling",
                                                 "Bile secretion"),
                             min_count = 2,
                             title = "Detritusphere, 12 weeks")
#abc transporters (choline, adenosine), nucleotide metabolism (uracil [pyrimidine], adenosine),
#biosynthesis of secondary metabolites (L-citrulline)

m_det.12 <- enrichMKEGG(gene = kegg_ids_det.12.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_m_det.12 <- enrichplot::dotplot(m_det.12)

#8weeks
k_rhizo.8 <- enrichKEGG(gene = kegg_ids_rhizo.8.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_rhizo.8 <- filter_and_plot(k_rhizo.8,
                              min_count = 2,
                              title = "Rhizosphere, 8 weeks")

m_rhizo.8 <- enrichMKEGG(gene = kegg_ids_rhizo.8.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_m_rhizo.8 <- enrichplot::dotplot(m_rhizo.8)

k_det.8 <- enrichKEGG(gene = kegg_ids_det.8.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_det.8 <- filter_and_plot(k_det.8,
                            remove_patterns = c("Neuroactive ligand signaling",
                                                "Bile secretion"),
                            min_count = 2,
                            title = "Detritusphere, 8 weeks")

m_det.8 <- enrichMKEGG(gene = kegg_ids_det.8.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_m_rhizo.8 <- enrichplot::dotplot(m_det.8)


#4weeks
k_rhizo.4 <- enrichKEGG(gene = kegg_ids_rhizo.4.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_rhizo.4 <- filter_and_plot(k_rhizo.4, 
                              min_count = 2,
                              title = "Rhizosphere, 4weeks")

m_rhizo.4 <- enrichMKEGG(gene = kegg_ids_rhizo.4.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_m_rhizo.4 <- enrichplot::dotplot(m_rhizo.4)

k_det.4 <- enrichKEGG(gene = kegg_ids_det.4.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_det.4 <- filter_and_plot(k_det.4, 
                              min_count = 2,
                              title = "Detritusphere, 4weeks")

m_det.4 <- enrichMKEGG(gene = kegg_ids_det.4.no.na$KEGG, organism = 'cpd', pvalueCutoff = 0.05)
dp_m_det.4 <- enrichplot::dotplot(m_det.4)

#plot compound classes
#12 weeks

kegg_ids_combined <- rbind(kegg_ids_rhizo.12, kegg_ids_rhizo.det.12, kegg_ids_det.12,
                           kegg_ids_rhizo.8, kegg_ids_rhizo.det.8, kegg_ids_det.8,
                           kegg_ids_rhizo.4, kegg_ids_rhizo.det.4, kegg_ids_det.4)
kegg_ids.plot <- kegg_ids_combined %>%
  mutate(habitat = factor(habitat, levels = c("Rhizo", "RhizoDet", "Detritus"))) %>%
  group_by(`Super class`, habitat, Timepoint) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks"))) %>%
  summarise(count = n()) %>%
  ggplot(aes(x=Timepoint, y=count, fill = `Super class`)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~habitat, ncol = 3) +
  scale_fill_manual(values = superclass_base_colors) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("./figures/LC-MS/manuscript/Fig.3e-class-time.png", dpi = 300,
      height = 3, width = 6, units = "in")

ggsave("./figures/LC-MS/manuscript/Fig.3e-class-time.pdf", dpi = 300,
       height = 3, width = 6, units = "in")





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
  separate_wider_delim(SampleID, 
                       delim = ".",
                       names = c("SampleID", "Treatment", "Zone", "Time")) %>%
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
filename <- paste0("./figures/manuscript/Fig.3-h-nmds.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0("./figures/manuscript/Fig.3-h-nmds.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)



### nmds plot (treatment and zone)
nmds_plot <-  make_nmds_plot(data.scores, Time, Treatment, col_list = timepoint_colors)
nmds_plot
filename <- paste0("./figures/manuscript/h-nmds-timepoint-treatment.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

#permanova
distmat <- vegdist(t(as.data.frame(df)), dist = "euclidean")


adonis2(distmat ~ Zone+Treatment+Time, data.scores, permutations = 999, by = "margin")
# adonis2(formula = distmat ~ Zone + Treatment + Time, data = data.scores, permutations = 999, by = "margin")
# Df  SumOfSqs      R2       F Pr(>F)    
# Zone       3 0.0053509 0.29286 14.3181  0.001 ***
#   Treatment  1 0.0006873 0.03762  5.5176  0.002 ** 
#   Time       2 0.0010360 0.05670  4.1582  0.001 ***
#   Residual  90 0.0112114 0.61361                   
# Total     96 0.0182712 1.00000    


