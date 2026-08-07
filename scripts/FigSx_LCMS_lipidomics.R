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
file_path <- "./data/Metabolomics/2022_4_Replicates/DRIPSIP_LC_Metab_WK_CleanedUp.xlsx"
#kegg <- read.csv("./qSIP_output/LC-MS/HPOS_metaboanalyst_KEGGID_all.csv", header = TRUE)
#kegg$Query <- gsub(" ", "_", kegg$Query) 


##load extra functions
source("./scripts/extra_functions_PCA_NMDS_plotting.R")
source("./scripts/lc-ms-functions.R")

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

#sig diff compounds (pairwise comparison)
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

df.z <- as.data.frame(t(scale(t(df))))

df.z.sig <- as.data.frame(t(scale(t(df.sig))))

#save normalized and z-scaled dataframes
write.csv(df, "./output/omics_datasets/lcms_met_h_norm_blank_account.csv")
write.csv(df.z, "./output/omics_datasets/lcms_met_h_norm_z_blank_account.csv")

#load csv
df <- read.csv("./output/omics_datasets/lcms_met_h_norm_blank_account.csv", header = TRUE) %>%
  column_to_rownames(var = "X")
df.z <- read.csv("./output/omics_datasets/lcms_met_h_norm_z_blank_account.csv", header = TRUE) %>%
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

outdir <- "./Figures/LC-MS/blanks_removed"

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
    label  = paste0("lcms-h-", tp),
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
ggsave("./figures/LC-MS/blanks_removed/class-time.png", dpi = 300,
      height = 3, width = 8, units = "in")

ggsave("./figures/LC-MS/blanks_removed/class-time.pdf", dpi = 300,
       height = 3, width = 8, units = "in")




########################################## PCA 
###parameters for plots
h = 6.5
w = 9
res = 300
size = 13
color_group = c("red", "orange", "green",  "blue",  "purple", "pink")
color_moisture = c("white", "black")

list_of_shapes <- c(15,17,16)

#move rowname to column in metadata

meta.new <- meta %>%
  rownames_to_column(var = "SampleID")

# Calculate PCA with prcomp()
pca <- prcomp(t(df), scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot


scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
# filename <- paste0("./Figures/LC-MS/PCA_screeplot-log-norm-pos-neg_100.png")
# ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
# filename <- paste0("./Figures/LC-MS/PCA_cumulative_var-log-norm-pos-neg_100.png")
# ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- merge(pca_coordinates, meta.new, by = "SampleID")



# write.csv(pca_coordinates,file=paste0("./qSIP_output/LC-MS/PCA_individual_coordinates-log-norm-pos-neg_100.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
arrows.f <- arrows %>%
  filter(contrib > .9)
# write.csv(arrows,file=paste0("./qSIP_output/LC-MS/PCA_vector_coordinates-log-pos-neg_100.csv"), row.names=TRUE)



#pca time, moisture
pca_plot <-  #pca_biplot +
  ggplot() +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, color = Zone, shape = Treatment), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = color_group) +
  #scale_shape_manual(values= list_of_shapes) +
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/LC-MS/blanks_removed/PCA-time-moisture-pos-neg.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows.f, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows.f, aes(x=xend, y=yend), color = "black",
                  label=arrows.f$name, size=size/3, show.legend = FALSE)

pca_biplot 
filename <- paste0("./Figures/LC-MS/blanks_removed/PCA-time-moisture-pos-neg-biplot.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)


##nmds
# Calculate nmds
#first transpose
df.t <- df %>%
  t() %>%
  as.data.frame()

set.seed(123)
nmds = metaMDS(df.t, distance = "bray")
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
w=9
h=5
nmds_plot <-  make_nmds_plot(data.scores, Zone, Treatment, col_list = zone_colors)
nmds_plot
filename <- paste0("./figures/LC-MS/blanks_removed/h-nmds-zone-treatment.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)
filename <- paste0("./figures/LC-MS/blanks_removed/h-nmds-zone-treatment.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

### nmds plot (treatment and zone)
nmds_plot <-  make_nmds_plot(data.scores, Treatment, Zone, col_list = col_list_treatment)
nmds_plot
filename <- paste0("./figures/LC-MS/blanks_removed/h-nmds-treatment-zone.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)


### nmds plot (treatment and zone)
nmds_plot <-  make_nmds_plot(data.scores, Time, Treatment, col_list = timepoint_colors)
nmds_plot
filename <- paste0("./figures/LC-MS/blanks_removed/h-nmds-timepoint-treatment.png")
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


##rp-POS
#reformat and filter pos table
peak_table_rpos.0$Name <- gsub(" ", "_", peak_table_rpos.0$Name) 
colnames(peak_table_rpos.0)[1] <- "mz"
colnames(peak_table_rpos.0)[2] <- "rt"
colnames(peak_table_rpos.0) = gsub("_RP_Pos", "", colnames(peak_table_rpos.0))

peak_table_rpos <- peak_table_rpos.0 %>%
  filter(Tags == "Confirmed ID (HIgh Confidence)" |
           Tags == "Match with Mass and RT only" |
           Tags == "Not in Internal DB. MS/MS match only" |
           Tags == "Not in Internal DB. Potential ID (Low Confidence)" |
           Tags == "Potential ID (Low Confidence)")

#order column names to match neg table
peak_table_rpos2 <-  peak_table_rpos[,-c(4:11)] %>%
  unite(mz_name, c("mz", "Name"), sep = "_", remove = FALSE) %>%
  column_to_rownames(var = "mz_name") %>%
  select(order(colnames(.))) %>%
  rownames_to_column(var = "mz_name") 

colnames(peak_table_rpos2) <- str_remove(colnames(peak_table_rpos2), "_HPOS")


compound.class.rpos <- peak_table_rpos[,c(3, 5:9)]



edata.0 <- peak_table_rpos2 %>%
  select(-mz, -rt) %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE) 


emeta.pos <-  peak_table_rpos2 %>%
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
  separate(SampleID, c("ID", "Treatment0"), remove = FALSE) %>%
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
rp_pos_clean <- blank_fold_filter_raw(
  edata.final, 
  fdata,  
  edata_cname = "Name",
  fdata_cname = "SampleID_new",
  blank_label = "Blank",
  sample_type_cname = "SampleType",
  fold_threshold = 3,
  remove_blanks = TRUE
)

###rp-neg
#reformat and filter pos table
peak_table_rneg.0$Name <- gsub(" ", "_", peak_table_rneg.0$Name) 
colnames(peak_table_rneg.0)[1] <- "mz"
colnames(peak_table_rneg.0)[2] <- "rt"
colnames(peak_table_rneg.0) = gsub("_RP_Neg", "", colnames(peak_table_rneg.0))

peak_table_rneg <- peak_table_rneg.0 %>%
  filter(Tags == "Confirmed ID (HIgh Confidence)" |
           Tags == "Match with Mass and RT only" |
           Tags == "Not in Internal DB. MS/MS match only" |
           Tags == "Not in Internal DB. Potential ID (Low Confidence)" |
           Tags == "Potential ID (Low Confidence)")

#order column names to match neg table
peak_table_rneg2 <-  peak_table_rneg[,-c(4:11)] %>%
  unite(mz_name, c("mz", "Name"), sep = "_", remove = FALSE) %>%
  column_to_rownames(var = "mz_name") %>%
  select(order(colnames(.))) %>%
  select(-SRFA1) %>%
  rownames_to_column(var = "mz_name") 


compound.class.rneg <- peak_table_rneg[,c(3, 5:9)]



edata.0 <- peak_table_rneg2 %>%
  select(-mz, -rt) %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE) 


emeta.neg <-  peak_table_rneg2 %>%
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
  separate(SampleID, c("ID", "Treatment0"), remove = FALSE) %>%
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
rp_neg_clean <- blank_fold_filter_raw(
  edata.final, 
  fdata,  
  edata_cname = "Name",
  fdata_cname = "SampleID_new",
  blank_label = "Blank",
  sample_type_cname = "SampleType",
  fold_threshold = 3,
  remove_blanks = TRUE
)


#merge pos and neg and keep max of intensity for duplicates
rp_clean_edata <- rbind(rp_pos_clean$edata, rp_neg_clean$edata)
rp_clean_fdata <- rbind(rp_pos_clean$fdata, rp_neg_clean$fdata)
rp_clean_emeta <- rbind(emeta.pos, emeta.neg)

rp_clean_edata <- rp_clean_edata %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE)

rp_clean_emeta <- rp_clean_emeta %>%
  group_by(Name) %>%
  summarise_if(is.numeric, max, na.rm = TRUE) 


####create omics data object for rp
rp_object <- as.metabData(
  e_data = rp_clean_edata,
  f_data = rp_clean_fdata,
  e_meta = rp_clean_emeta,
  edata_cname = "Name",
  fdata_cname = "SampleID_new",
  emeta_cname = "Name",
  data_scale = "abundance"
)

class(rp_object)
summary(rp_object)


#plot transformed data
plot(edata_transform(rp_object, data_scale = "log2"))

####filter biomolecules
mymolfilt <- molecule_filter(omicsData = rp_object)
class(mymolfilt)
summary(mymolfilt) #all metabolites are found in all samples?

# coefficient of variation filter
rp_object_groups <- group_designation(omicsData = rp_object, main_effects = "Zone")

mycvfilt <- cv_filter(omicsData = rp_object_groups)
plot(mycvfilt, cv_threshold = 97)

rp_object_groups_log2 <- edata_transform(rp_object_groups, data_scale = "log2")
mycvfilt_log2 <- cv_filter(omicsData = rp_object_groups_log2)
plot(mycvfilt_log2, cv_threshold = 97)

summary(rp_object_groups_log2, cv_threshold = 97)
# Number Filtered Biomolecules: 21 


#apply the filter, this removes 8 biomolecules
rp_object_groups_log2 <- applyFilt(filter_object = mycvfilt_log2, omicsData = rp_object_groups_log2, cv_threshold = 97)
summary(rp_object_groups_log2)

#imd-anova filter
myimdanovafilt <- imdanova_filter(omicsData = rp_object_groups_log2)
summary(myimdanovafilt, min_nonmiss_anova = 2, min_nonmiss_gtest = 3) #no biomolecules filtered

####filter samples
myfilter <- rmd_filter(omicsData = rp_object_groups_log2, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness"))
plot(myfilter)

summary(myfilter, pvalue_threshold = 0.00000001)

plot(myfilter, pvalue_threshold = 0.0000001, bw_theme = TRUE)

# get vector of potential outliers
potential_outliers <- summary(myfilter, pvalue_threshold = 0.00000001)$filtered_samples
# Filtered Samples: DRIP_51.Drought.Rhizo.4weeks 


# loop over potential outliers to generate plots
if (length(potential_outliers) > 0) {
  for (i in 1:length(potential_outliers)) {
    print(plot(myfilter, sampleID = potential_outliers[i]))
  }
}

mycor <- cor_result(omicsData = rp_object_groups_log2)
plot(mycor, interactive = TRUE)

mypca <- dim_reduction(omicsData = rp_object_groups_log2)

plot(mypca, interactive = TRUE)

#remove outliers found using metric above
rp_object_groups_log2 <- applyFilt(filter_object = myfilter, omicsData = rp_object_groups_log2, pvalue_threshold = 0.00001)
summary(rp_object_groups_log2)

#numeric summary
edata_summary(rp_object_groups_log2, by = "molecule", groupvar = NULL)

####normalization

# global median centering - apply the norm
rp_object_groups_log2_norm <- normalize_global(
  omicsData = rp_object_groups_log2,
  subset_fn = "all",
  norm_fn = "mean",
  apply_norm = TRUE,
  backtransform = TRUE
)
class(rp_object_groups_log2_norm)

plot(rp_object_groups_log2_norm)

####statstical analysis
#anova
rp_object_groups_log2_norm_object <- group_designation(
  omicsData = rp_object_groups_log2_norm,
  main_effects = c("Zone")
)

all_pairwise_results <- imd_anova(
  omicsData = rp_object_groups_log2_norm_object,
  test_method = "anova"
)

rp_object_groups_log2_norm_object_pairwise <- imd_anova(
  omicsData = rp_object_groups_log2_norm_object,
  test_method = "anova",
  pval_adjust_a_multcomp = "Tukey"
)

plot(rp_object_groups_log2_norm_object_pairwise, plot_type = "volcano")



###view PCA of normalized data
mycor <- cor_result(omicsData = rp_object_groups_log2_norm_object)

mypca <- dim_reduction(omicsData = rp_object_groups_log2_norm_object)

plot(mypca, interactive = TRUE)

#sig diff compounds (pairwise comparison)
rhizo.vs.bulk.up <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                               l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_Rhizo_vs_Bulk, 
                               p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_Rhizo_vs_Bulk) %>%
  filter(p.value < 0.1 &
           l2fc > 0)

rhizo.vs.bulk.down <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                 l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_Rhizo_vs_Bulk, 
                                 p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_Rhizo_vs_Bulk) %>%
  filter(p.value < 0.1 &
           l2fc < 0)

rhizo.vs.rhizodetritus.up <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                        l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_Rhizo_vs_RhizoDetritus, 
                                        p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_Rhizo_vs_RhizoDetritus) %>%
  filter(p.value < 0.1 &
           l2fc > 0)

rhizo.vs.rhizodetritus.down <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                          l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_Rhizo_vs_RhizoDetritus, 
                                          p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_Rhizo_vs_RhizoDetritus) %>%
  filter(p.value < 0.1 &
           l2fc < 0)

rhizo.vs.detritus.up <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                   l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_Rhizo_vs_Detritus, 
                                   p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_Rhizo_vs_Detritus) %>%
  filter(p.value < 0.1 &
           l2fc > 0)

rhizo.vs.detritus.down <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                     l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_Rhizo_vs_Detritus, 
                                     p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_Rhizo_vs_Detritus) %>%
  filter(p.value < 0.1 &
           l2fc < 0)

rhizodet.vs.detritus.up <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                      l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_RhizoDetritus_vs_Detritus, 
                                      p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_RhizoDetritus_vs_Detritus) %>%
  filter(p.value < 0.1 &
           l2fc > 0)

rhizodet.vs.detritus.down <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                        l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_RhizoDetritus_vs_Detritus, 
                                        p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_RhizoDetritus_vs_Detritus) %>%
  filter(p.value < 0.1 &
           l2fc < 0)

rhizodet.vs.bulk.up <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                  l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_RhizoDetritus_vs_Bulk, 
                                  p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_RhizoDetritus_vs_Bulk) %>%
  filter(p.value < 0.1 &
           l2fc > 0)

rhizodet.vs.bulk.down <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                                    l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_RhizoDetritus_vs_Bulk, 
                                    p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_RhizoDetritus_vs_Bulk) %>%
  filter(p.value < 0.1 &
           l2fc < 0)


bulk.vs.det.up <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                             l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_Bulk_vs_Detritus, 
                             p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_Bulk_vs_Detritus) %>%
  filter(p.value < 0.1 &
           l2fc > 0)

bulk.vs.det.down <- data.frame(compounds=rp_object_groups_log2_norm_object_pairwise$Name,
                               l2fc=rp_object_groups_log2_norm_object_pairwise$Fold_change_Bulk_vs_Detritus, 
                               p.value=rp_object_groups_log2_norm_object_pairwise$P_value_A_Bulk_vs_Detritus) %>%
  filter(p.value < 0.1 &
           l2fc < 0)


#sig compounds to filter heatmap
sig.com <- rbind(rhizo.vs.bulk.up, rhizo.vs.bulk.down, rhizo.vs.rhizodetritus.up, 
                 rhizo.vs.rhizodetritus.down, rhizo.vs.detritus.up, rhizo.vs.detritus.down,
                 rhizodet.vs.bulk.up, rhizo.vs.bulk.down, bulk.vs.det.up, bulk.vs.det.down)

sig.com.unique <- sig.com %>%
  distinct(compounds)



##extract data and 
df.0 <- as.data.frame(rp_object_groups_log2_norm_object$e_data)
df <- df.0 %>%
  rownames_to_column(var = "row") %>%
  column_to_rownames(var = "Name") %>%
  select(-row)

df.sig <- df %>%
  filter(rownames(.) %in% sig.com.unique$compounds)

df.z <- as.data.frame(t(scale(t(df))))

df.z.sig <- as.data.frame(t(scale(t(df.sig))))

#save normalized and z-scaled dataframes
write.csv(df, "./output/omics_datasets/lcms_met_rp_norm_blank_account.csv")
write.csv(df.z, "./output/omics_datasets/lcms_met_rp_norm_z_blank_account.csv")

#load csv
df <- read.csv("./output/omics_datasets/lcms_met_rp_norm_blank_account.csv", header = TRUE) %>%
  column_to_rownames(var = "X")
df.z <- read.csv("./output/omics_datasets/lcms_met_rp_norm_z_blank_account.csv", header = TRUE) %>%
  column_to_rownames(var = "X") 

# ####classifyer
# #load compound classifications
# meta.comp <- read.csv("./output/LC-MS/classifications-hilic.csv") %>%
#   column_to_rownames(var = "X") 
# 
# meta.comp.sub <- meta.comp %>%
#   select(subclass)
# 
# meta.comp.class  <- meta.comp %>%
#   select(class)
# 
# meta.comp.supclass  <- meta.comp %>%
#   select(superclass)
# 
# #add classifications
# meta.comp.to.merge <- meta.comp %>%
#   rownames_to_column(var = "compound")
# 
# df.z.to.merge <- df.z %>%
#   rownames_to_column(var = "compound")
# 
# meta.comp.df <- meta.comp.to.merge %>%
#   merge(df.z.to.merge, by = "compound", all = TRUE) 


# ####compound classes in original file
# #load compound classifications
# meta.comp.pos <- peak_table_pos.0 %>%
#   select(Name, `Super class`, `Main class`, `Sub class`, Formula)
# 
# meta.comp.neg <- peak_table_neg.0 %>%
#   select(Name, `Super class`, `Main class`, `Sub class`, Formula)
# 
# meta.comp <- rbind(meta.comp.pos, meta.comp.neg) %>%
#   unique() %>%
#   column_to_rownames(var = "Name")
# 
# meta.comp.sup <- meta.comp %>%
#   select(`Super class`) %>%
#   rename(superclass = `Super class`) %>%
#   filter(!is.na(superclass)) %>%
#   filter(rownames(.) %in% rownames(df.z)) %>%
#   mutate(superclass = factor(superclass))
# 
# meta.comp.sup.sig <- meta.comp %>%
#   select(`Super class`) %>%
#   rename(superclass = `Super class`) %>%
#   filter(!is.na(superclass)) %>%
#   filter(rownames(.) %in% rownames(df.z.sig)) %>%
#   mutate(superclass = factor(superclass))
# 
# 
# 
# meta.comp.class  <- meta.comp %>%
#   select(`Main class`) %>%
#   rename(class = `Main class`) %>%
#   filter(!is.na(class)) %>%
#   filter(rownames(.) %in% rownames(df.z)) %>%
#   mutate(class = factor(class))
# 
# meta.comp.class.sig  <- meta.comp %>%
#   select(`Main class`) %>%
#   rename(class = `Main class`) %>%
#   filter(!is.na(class)) %>%
#   filter(rownames(.) %in% rownames(df.z.sig)) %>%
#   mutate(class = factor(class))
# 
# 
# meta.comp.subclass  <- meta.comp %>%
#   select(`Sub class`) %>%
#   rename(subclass = `Sub class`) %>%
#   filter(!is.na(subclass)) %>%
#   filter(rownames(.) %in% rownames(df.z)) %>%
#   mutate(subclass = factor(subclass))
# 
# meta.comp.subclass.sig  <- meta.comp %>%
#   select(`Sub class`) %>%
#   rename(subclass = `Sub class`) %>%
#   filter(!is.na(subclass)) %>%
#   filter(rownames(.) %in% rownames(df.z.sig)) %>%
#   mutate(subclass = factor(subclass)) 
# 
# #add classifications
# meta.comp.to.merge <- meta.comp %>%
#   rownames_to_column(var = "compound")
# 
# df.z.to.merge <- df.z %>%
#   rownames_to_column(var = "compound")
# 
# meta.comp.df <- meta.comp.to.merge %>%
#   merge(df.z.to.merge, by = "compound", all.y = TRUE) 
#   
# 
# 
# 
# 
# #heatmap
# meta <- fdata %>%
#   column_to_rownames(var = "SampleID_new") %>%
#   select(-Treatment0, -ID, -SampleID, -SampleType) %>%
#   drop_na()
# 
# #meta$Zone <- factor(meta$Zone, levels = c("Rhizo", "RhizoDetritus", "Detritus", "Bulk"))
# meta$Treatment <- factor(meta$Treatment, levels = c("Drought", "Untrt"))
# meta$Timepoint <- factor(meta$Timepoint, levels = c("4weeks", "8weeks", "12weeks"))
# 
# # Get the levels of your row annotation (subclass)
# subclass_levels <- levels(meta.comp.subclass$subclass) 
# subclass_levels.sig <- levels(meta.comp.subclass.sig$subclass)# adjust column name as needed
# # If meta.comp.subclass is just a factor, then:
# # subclass_levels <- levels(meta.comp.subclass)
# 
# # Build a Paired color set of the right length
# # subclass_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(subclass_levels))
# # names(subclass_colors) <- subclass_levels
# # 
# # subclass_colors.sig <- colorRampPalette(brewer.pal(12, "Paired"))(length(subclass_levels.sig))
# # names(subclass_colors.sig) <- subclass_levels.sig
# 
# subclass_colors <- setNames(
#   pals::polychrome(length(subclass_levels)),
#   subclass_levels
# )
# 
# subclass_colors.sig <- setNames(
#   pals::polychrome(length(subclass_levels.sig)),
#   subclass_levels.sig
# )
# 
# 
# # Get the levels of your row annotation (class)
# class_levels <- levels(meta.comp.class$class) 
# class_levels.sig <- levels(meta.comp.class.sig$class)# adjust column name as needed
# # If meta.comp.subclass is just a factor, then:
# # subclass_levels <- levels(meta.comp.subclass)
# 
# # Build a Paired color set of the right length
# # class_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(class_levels))
# # names(class_colors) <- class_levels
# # 
# # class_colors.sig <- colorRampPalette(brewer.pal(12, "Paired"))(length(class_levels.sig))
# # names(class_colors.sig) <- class_levels.sig
# 
# class_colors <- setNames(
#   pals::polychrome(length(class_levels)),
#   class_levels
# )
# 
# class_colors.sig <- setNames(
#   pals::polychrome(length(class_levels.sig)),
#   class_levels.sig
# )
# 
# 
# # Get the levels of your row annotation (superclass)
# superclass_levels <- levels(meta.comp.sup$superclass) 
# superclass_levels.sig <- levels(meta.comp.sup.sig$superclass)# adjust column name as needed
# # If meta.comp.subclass is just a factor, then:
# # subclass_levels <- levels(meta.comp.subclass)
# 
# # Build a Paired color set of the right length
# # superclass_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(superclass_levels))
# # names(superclass_colors) <- superclass_levels
# # 
# # superclass_colors.sig <- colorRampPalette(brewer.pal(12, "Paired"))(length(superclass_levels.sig))
# # names(superclass_colors.sig) <- superclass_levels.sig
# 
# superclass_colors <- setNames(
#   pals::polychrome(length(superclass_levels)),
#   superclass_levels
# )
# 
# superclass_colors.sig <- setNames(
#   pals::polychrome(length(superclass_levels.sig)),
#   superclass_levels.sig
# )
# 
# #annotation colors for subclass
# ann_colors = list(
#   Zone = c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red"),
#   Treatment = c(Drought = "brown", Untrt = "darkgreen"),
#   Timepoint = c(`4weeks` = "pink", `8weeks` = "lightblue", `12weeks` = "lightyellow"),
#   subclass = subclass_colors
# )
# 
# ann_colors.sig = list(
#   Zone = c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red"),
#   Treatment = c(Drought = "brown", Untrt = "darkgreen"),
#   Timepoint = c(`4weeks` = "pink", `8weeks` = "lightblue", `12weeks` = "lightyellow"),
#   subclass = subclass_colors.sig
# )
# 
# #annotation colors for class
# ann_colors_class = list(
#   Zone = c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red"),
#   Treatment = c(Drought = "brown", Untrt = "darkgreen"),
#   Timepoint = c(`4weeks` = "pink", `8weeks` = "lightblue", `12weeks` = "lightyellow"),
#   class = class_colors
# )
# 
# ann_colors.class.sig = list(
#   Zone = c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red"),
#   Treatment = c(Drought = "brown", Untrt = "darkgreen"),
#   Timepoint = c(`4weeks` = "pink", `8weeks` = "lightblue", `12weeks` = "lightyellow"),
#   class = class_colors.sig
# )
# 
# #annotation colors for superlass
# ann_colors_superclass = list(
#   Zone = c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red"),
#   Treatment = c(Drought = "brown", Untrt = "darkgreen"),
#   Timepoint = c(`4weeks` = "pink", `8weeks` = "lightblue", `12weeks` = "lightyellow"),
#   superclass = superclass_colors
# )
# 
# ann_colors.superclass.sig = list(
#   Zone = c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red"),
#   Treatment = c(Drought = "brown", Untrt = "darkgreen"),
#   Timepoint = c(`4weeks` = "pink", `8weeks` = "lightblue", `12weeks` = "lightyellow"),
#   superclass = superclass_colors.sig
# )
# 



####compound classes in original file
# load compound classifications
meta.comp.pos <- peak_table_rpos.0 %>% select(Name, `Super class`, `Main class`, `Sub class`, Formula)
meta.comp.neg <- peak_table_rneg.0 %>% select(Name, `Super class`, `Main class`, `Sub class`, Formula)

meta.comp <- rbind(meta.comp.pos, meta.comp.neg) %>%
  unique() %>%
  drop_na() %>%
  column_to_rownames(var = "Name") %>%
  mutate(
    # Fix superclass errors
    `Super class` = case_when(
      `Super class` == "Glycerophospholipids" ~ "Lipids",
      `Super class` == "Sphingolipids"        ~ "Lipids",
      TRUE ~ `Super class`
    ),
    # Fix capitalization inconsistencies
    `Main class` = case_when(
      `Main class` == "Prenol Lipids"                    ~ "Prenol lipids",
      `Main class` == "Sterol Lipids"                    ~ "Sterol lipids",
      `Main class` == "Carboxylic acids and derivitives" ~ "Carboxylic acids and derivatives",
      `Main class` == "Glycerophophocholines"            ~ "Glycerophosphocholines",
      TRUE ~ `Main class`
    ),
    # Fix subclass typos
    `Sub class` = case_when(
      `Sub class` == "Sesquiterpeoids"        ~ "Sesquiterpenoids",
      `Sub class` == "Glycerophophocholines"  ~ "Glycerophosphocholines",
      `Sub class` == "Oxolipins"              ~ "Oxylipins",
      TRUE ~ `Sub class`
    )
  )

# ── Define superclass base colors ─────────────────────────────────────────────
# ── Define superclass base colors ─────────────────────────────────────────────
# All superclasses present in your data:
# Lipids, Fatty Acyls, Organic acids, Organic nitrogen compounds,
# Nucleic acids, Benzenoids, Carbohydrates, Organic oxygen compounds,
# Alkaloids, Polyketides, Organoheterocyclic compounds

superclass_base_colors <- c(
  "Lipids"                       = "#E41A1C",  # red
  "Fatty Acyls"                  = "#FF7F00",  # orange
  "Organic acids"                = "#4DAF4A",  # green
  "Organic nitrogen compounds"   = "#377EB8",  # blue
  "Nucleic acids"                = "#984EA3",  # purple
  "Benzenoids"                   = "#A65628",  # brown
  "Carbohydrates"                = "#F768A1",  # pink
  "Organic oxygen compounds"     = "#999999",  # grey
  "Alkaloids"                    = "#FFED6F",  # yellow
  "Polyketides"                  = "#1B9E77",  # teal
  "Organoheterocyclic compounds" = "#B2DF8A"   # light green
)

# ── Define superclass → main class hierarchy ──────────────────────────────────
# Built directly from your actual data

class_hierarchy <- data.frame(
  superclass = c(
    # Lipids
    "Lipids",       "Lipids",         "Lipids",
    "Lipids",       "Lipids",
    # Fatty Acyls
    "Fatty Acyls",  "Fatty Acyls",    "Fatty Acyls",  "Fatty Acyls",
    # Organic acids
    "Organic acids", "Organic acids", "Organic acids",
    # Organic nitrogen compounds
    "Organic nitrogen compounds", "Organic nitrogen compounds",
    "Organic nitrogen compounds",
    # Nucleic acids
    "Nucleic acids", "Nucleic acids",
    # Benzenoids
    "Benzenoids",   "Benzenoids",     "Benzenoids",
    # Alkaloids
    "Alkaloids",    "Alkaloids",
    # Polyketides
    "Polyketides",  "Polyketides",    "Polyketides",
    # Carbohydrates
    "Carbohydrates", "Carbohydrates", "Carbohydrates",
    # Organic oxygen compounds
    "Organic oxygen compounds", "Organic oxygen compounds",
    # Organoheterocyclic compounds
    "Organoheterocyclic compounds"
  ),
  class = c(
    # Lipids — from your data:
    # Glycerophospholipids (1-Oleoyl-rac-Glycerol, 2-oleoyl-sn-glycero-3-phosphocholine)
    # Prenol lipids (diterpenoids, sesquiterpenoids, Retinal, Glycyrrhetinicacid)
    # Sphingolipids (Sphinganine)
    # Sterol Lipids (SimilaritytoDehydrocholicacid, DeoxycholicAcid)
    "Glycerophospholipids", "Prenol lipids",  "Sphingolipids",
    "Sterol Lipids",        "Prenol Lipids",  # note: both capitalizations appear in your data
    # Fatty Acyls — from your data:
    # Fatty acids, Fatty esters, Fatty amides, Octadecanoids
    "Fatty acids",          "Fatty esters",   "Fatty amides",   "Octadecanoids",
    # Organic acids — from your data:
    # Amino acids and peptides, Carboxylic acids and derivitives, Keto acids
    "Amino acids and peptides", "Carboxylic acids and derivitives", "Keto acids",
    # Organic nitrogen compounds — from your data:
    # Amides, Carnitines, Cholines
    "Amides",               "Carnitines",     "Cholines",
    # Nucleic acids — from your data:
    # Purines, Pyrimidines
    "Purines",              "Pyrimidines",
    # Benzenoids — from your data:
    # Benzophenones, Benzamides, Anilines
    "Benzophenones",        "Benzamides",     "Anilines",
    # Alkaloids — from your data:
    # Tryptophan alkaloids, Nicotinic acid alkaloids
    "Tryptophan alkaloids", "Nicotinic acid alkaloids",
    # Polyketides — from your data:
    # Flavonoids, Phenolic acids, Phenylpropanoids
    "Flavonoids",           "Phenolic acids", "Phenylpropanoids",
    # Carbohydrates — from your data:
    # Monosaccharides, Disaccharides, Oligosaccharides
    "Monosaccharides",      "Disaccharides",  "Oligosaccharides",
    # Organic oxygen compounds — from your data:
    # Organooxygen compounds, Alcohols and polyols
    "Organooxygen compounds", "Alcohols and polyols",
    # Organoheterocyclic compounds — from your data:
    # Indolyl carboxylic acids
    "Indolyl carboxylic acids"
  ),
  stringsAsFactors = FALSE
)
# ── Function to generate shades of a base color ───────────────────────────────
generate_shades <- function(base_color, n) {
  if (n == 1) return(setNames(base_color, NULL))
  colorspace::lighten(base_color, amount = seq(0.5, -0.3, length.out = n))
}

# ── Function to build class colors inheriting from superclass ─────────────────
build_class_colors <- function(hierarchy, base_colors, present_classes) {
  class_colors <- c()
  for (sc in names(base_colors)) {
    classes_in_sc <- hierarchy$class[hierarchy$superclass == sc]
    classes_in_sc <- classes_in_sc[classes_in_sc %in% present_classes]
    if (length(classes_in_sc) == 0) next
    shades <- generate_shades(base_colors[sc], length(classes_in_sc))
    names(shades) <- classes_in_sc
    class_colors <- c(class_colors, shades)
  }
  return(class_colors)
}

# ── Build ordered class levels (superclass order → class order within) ────────
build_ordered_classes <- function(hierarchy, base_colors, present_classes) {
  class_hierarchy %>%
    filter(class %in% present_classes) %>%
    arrange(match(superclass, names(base_colors))) %>%
    pull(class)
}

# ── Build annotation data frames ──────────────────────────────────────────────

# superclass
meta.comp.sup <- meta.comp %>%
  select(`Super class`) %>%
  rename(superclass = `Super class`) %>%
  filter(!is.na(superclass)) %>%
  filter(rownames(.) %in% rownames(df.z)) %>%
  mutate(superclass = factor(superclass, levels = names(superclass_base_colors)))

meta.comp.sup.sig <- meta.comp %>%
  select(`Super class`) %>%
  rename(superclass = `Super class`) %>%
  filter(!is.na(superclass)) %>%
  filter(rownames(.) %in% rownames(df.z.sig)) %>%
  mutate(superclass = factor(superclass, levels = names(superclass_base_colors)))

# class — ordered by superclass
class_levels     <- unique(na.omit(meta.comp$`Main class`[rownames(meta.comp) %in% rownames(df.z)]))
class_levels.sig <- unique(na.omit(meta.comp$`Main class`[rownames(meta.comp) %in% rownames(df.z.sig)]))

ordered_classes     <- build_ordered_classes(class_hierarchy, superclass_base_colors, class_levels)
ordered_classes.sig <- build_ordered_classes(class_hierarchy, superclass_base_colors, class_levels.sig)

meta.comp.class <- meta.comp %>%
  select(`Main class`) %>%
  rename(class = `Main class`) %>%
  filter(!is.na(class)) %>%
  filter(rownames(.) %in% rownames(df.z)) %>%
  mutate(class = factor(class, levels = ordered_classes)) %>%
  arrange(class)

meta.comp.class.sig <- meta.comp %>%
  select(`Main class`) %>%
  rename(class = `Main class`) %>%
  filter(!is.na(class)) %>%
  filter(rownames(.) %in% rownames(df.z.sig)) %>%
  mutate(class = factor(class, levels = ordered_classes.sig)) %>%
  arrange(class)

# subclass — use polychrome since too many to assign manually
subclass_levels     <- levels(factor(na.omit(meta.comp$`Sub class`[rownames(meta.comp) %in% rownames(df.z)])))
subclass_levels.sig <- levels(factor(na.omit(meta.comp$`Sub class`[rownames(meta.comp) %in% rownames(df.z.sig)])))

meta.comp.subclass <- meta.comp %>%
  select(`Sub class`) %>%
  rename(subclass = `Sub class`) %>%
  filter(!is.na(subclass)) %>%
  filter(rownames(.) %in% rownames(df.z)) %>%
  mutate(subclass = factor(subclass))

meta.comp.subclass.sig <- meta.comp %>%
  select(`Sub class`) %>%
  rename(subclass = `Sub class`) %>%
  filter(!is.na(subclass)) %>%
  filter(rownames(.) %in% rownames(df.z.sig)) %>%
  mutate(subclass = factor(subclass))

# ── Build color vectors ────────────────────────────────────────────────────────

# superclass — manually defined base colors
superclass_colors     <- superclass_base_colors[names(superclass_base_colors) %in% levels(meta.comp.sup$superclass)]
superclass_colors.sig <- superclass_base_colors[names(superclass_base_colors) %in% levels(meta.comp.sup.sig$superclass)]

# class — shades inheriting from superclass base color
class_colors     <- build_class_colors(class_hierarchy, superclass_base_colors, ordered_classes)
class_colors.sig <- build_class_colors(class_hierarchy, superclass_base_colors, ordered_classes.sig)

# subclass — polychrome since too many to assign manually
subclass_colors     <- setNames(pals::polychrome(length(subclass_levels)),     subclass_levels)
subclass_colors.sig <- setNames(pals::polychrome(length(subclass_levels.sig)), subclass_levels.sig)

# ── Sample metadata ────────────────────────────────────────────────────────────
meta <- fdata %>%
  column_to_rownames(var = "SampleID_new") %>%
  select(-Treatment0, -ID, -SampleID, -SampleType) %>%
  drop_na()

meta$Treatment <- factor(meta$Treatment, levels = c("Drought", "Untrt"))
meta$Timepoint <- factor(meta$Timepoint, levels = c("4weeks", "8weeks", "12weeks"))

# ── Annotation color lists ─────────────────────────────────────────────────────
zone_colors <- c(Rhizo = "green", RhizoDetritus = "blue", Detritus = "orange", Bulk = "red")
treatment_colors <- c(Drought = "brown", Untrt = "darkgreen")
timepoint_colors <- c(`4weeks` = "#F4A582", `8weeks` = "#74ADD1", `12weeks` = "#D9EF8B")

ann_colors <- list(
  Zone      = zone_colors,
  Treatment = treatment_colors,
  Timepoint = timepoint_colors,
  subclass  = subclass_colors
)

ann_colors.sig <- list(
  Zone      = zone_colors,
  Treatment = treatment_colors,
  Timepoint = timepoint_colors,
  subclass  = subclass_colors.sig
)

ann_colors_class <- list(
  Zone      = zone_colors,
  Treatment = treatment_colors,
  Timepoint = timepoint_colors,
  class     = class_colors
)

ann_colors_class.sig <- list(
  Zone      = zone_colors,
  Treatment = treatment_colors,
  Timepoint = timepoint_colors,
  class     = class_colors.sig
)

ann_colors_superclass <- list(
  Zone       = zone_colors,
  Treatment  = treatment_colors,
  Timepoint  = timepoint_colors,
  superclass = superclass_colors
)

ann_colors_superclass.sig <- list(
  Zone       = zone_colors,
  Treatment  = treatment_colors,
  Timepoint  = timepoint_colors,
  superclass = superclass_colors.sig
)

w=25
h=10

png("./Figures/LC-MS/blanks_removed/lcms-rp-pmart-pos-neg-heatmap-sublcass.png", width =w, height = h+3, res = 300, units = "in")
pheatmap(df.z, cluster_rows = TRUE, 
         cluster_cols = TRUE,
         #scale = "row",
         annotation_col = meta,
         annotation_row = meta.comp.subclass,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         drop_levels = FALSE)
dev.off()

png("./Figures/LC-MS/blanks_removed/lcms-rp-pmart-pos-neg-heatmap-subclass-sig.png", width =w, height = h, res = 300, units = "in")
pheatmap(df.z.sig, cluster_rows = TRUE, 
         cluster_cols = TRUE,
         #scale = "row",
         annotation_col = meta,
         annotation_colors = ann_colors.sig,
         annotation_row = meta.comp.subclass.sig,
         show_colnames = FALSE,
         drop_levels = FALSE,
         fontsize = 20)
dev.off()

png("./Figures/LC-MS/blanks_removed/lcms-rp-pmart-pos-neg-heatmap-mainclass.png", width =w, height = h+3,  res = 300, units = "in")
pheatmap(df.z, cluster_rows = TRUE, 
         cluster_cols = TRUE,
         #scale = "row",
         annotation_col = meta,
         annotation_colors = ann_colors_class,
         annotation_row = meta.comp.class,
         show_colnames = FALSE,
         drop_levels = FALSE,
         fontsize = 20)
dev.off()

png("./Figures/LC-MS/blanks_removed/lcms-rp-pmart-pos-neg-heatmap-mainclass-sig.png", width =w, height = h,  res = 300, units = "in")
pheatmap(df.z.sig, cluster_rows = TRUE, 
         cluster_cols = TRUE,
         #scale = "row",
         annotation_col = meta,
         annotation_colors = ann_colors_class.sig,
         annotation_row = meta.comp.class.sig,
         show_colnames = FALSE,
         drop_levels = FALSE,
         fontsize = 20)
dev.off()

png("./Figures/LC-MS/blanks_removed/lcms-rp-pmart-pos-neg-heatmap-superclass.png", width =w, height = h+3,  res = 300, units = "in")
pheatmap(df.z, cluster_rows = TRUE, 
         cluster_cols = TRUE,
         #scale = "row",
         annotation_col = meta,
         annotation_colors = ann_colors_superclass,
         annotation_row = meta.comp.sup,
         show_colnames = FALSE,
         drop_levels = FALSE,
         fontsize = 20)
dev.off()

png("./Figures/LC-MS/blanks_removed/lcms-rp-pmart-pos-neg-heatmap-superclass-sig.png", width =w, height = h, res = 300, units = "in")
pheatmap(df.z.sig, cluster_rows = TRUE, 
         cluster_cols = TRUE,
         #scale = "row",
         annotation_col = meta,
         annotation_colors = ann_colors_superclass.sig,
         annotation_row = meta.comp.sup.sig,
         show_colnames = FALSE,
         drop_levels = FALSE,
         fontsize = 20)
dev.off()



########################################## PCA 
###parameters for plots
h = 6.5
w = 9
res = 300
size = 13
color_group = c("red", "orange", "green",  "blue",  "purple", "pink")
color_moisture = c("white", "black")

list_of_shapes <- c(15,16,17)

#move rowname to column in metadata

meta.new <- meta %>%
  rownames_to_column(var = "SampleID")

# Calculate PCA with prcomp()
pca <- prcomp(t(df), scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot


scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
# filename <- paste0("./Figures/LC-MS/PCA_screeplot-log-norm-pos-neg_100.png")
# ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
# filename <- paste0("./Figures/LC-MS/PCA_cumulative_var-log-norm-pos-neg_100.png")
# ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- merge(pca_coordinates, meta.new, by = "SampleID")



# write.csv(pca_coordinates,file=paste0("./qSIP_output/LC-MS/PCA_individual_coordinates-log-norm-pos-neg_100.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
arrows.f <- arrows %>%
  filter(contrib > .9)
# write.csv(arrows,file=paste0("./qSIP_output/LC-MS/PCA_vector_coordinates-log-pos-neg_100.csv"), row.names=TRUE)



#pca time, moisture
pca_plot <-  #pca_biplot +
  ggplot() +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, color = Zone, shape = Treatment), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = color_group) +
  #scale_shape_manual(values= list_of_shapes) +
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/LC-MS/blanks_removed/PCA-time-moisture-rp.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows.f, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows.f, aes(x=xend, y=yend), color = "black",
                  label=arrows.f$name, size=size/3, show.legend = FALSE)

pca_biplot 
filename <- paste0("./Figures/LC-MS/blanks_removed/PCA-time-moisture-rp-biplot.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)


##nmds
# Calculate nmds
#first transpose
df.t <- df %>%
  t() %>%
  as.data.frame()

set.seed(123)
nmds = metaMDS(df.t, distance = "bray")
plot(nmds)

scores <- scores(nmds)
scores.samples <- scores$site

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)

data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
  separate(., SampleID, into = c("SampleID", "Treatment", "Zone", "Time"), sep = "\\.") %>%
  mutate(Treatment = str_replace(Treatment, "Untrt", "Control")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Drought", "Control")))

# envfit with aligned metadata
ef <- envfit(nmds ~ Zone + Treatment + Time, data = data.scores, permutations = 999)
ef

### nmds plot (zone and treatment)
w=8
h=6
nmds_plot <-  make_nmds_plot(data.scores, Zone, Treatment, col_list = col_list_zone)
nmds_plot
filename <- paste0("./figures/LC-MS/blanks_removed/rp-nmds-zone-treatment.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

### nmds plot (treatment and zone)
nmds_plot <-  make_nmds_plot(data.scores, Treatment, Zone, col_list = col_list_treatment)
nmds_plot
filename <- paste0("./figures/LC-MS/blanks_removed/rp-nmds-treatment-zone.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)


### nmds plot (treatment and zone)
nmds_plot <-  make_nmds_plot(data.scores, Time, Treatment, col_list = timepoint_colors)
nmds_plot
filename <- paste0("./figures/LC-MS/blanks_removed/rp-nmds-timepoint-treatment.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,nmds_plot)

#permanova
distmat <- dist(t(as.data.frame(df)))


adonis2(distmat ~ Zone*Treatment*Time, data.scores)
# adonis2(formula = distmat ~ Zone * Treatment * Timepoint, data = meta.new.f)
# Df SumOfSqs      R2      F Pr(>F)    
# Model    23   2037.5 0.48317 2.3575  0.001 ***
#   Residual 58   2179.5 0.51683                  
# Total    81   4217.1 1.00000                    
# ---          

adonis2(distmat ~ Zone, data.scores)
# adonis2(formula = distmat ~ Zone, data = meta.new.f)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     3    862.5 0.20452 6.6847  0.001 ***
#   Residual 78   3354.6 0.79548                  
# Total    81   4217.1 1.00000     

adonis2(distmat ~ Treatment, data.scores)
# adonis2(formula = distmat ~ Treatment, data = meta.new.f)
# Df SumOfSqs      R2      F Pr(>F)    
# Model     1    166.8 0.03954 3.2936  0.001 ***
#   Residual 80   4050.3 0.96046                  
# Total    81   4217.1 1.00000   

