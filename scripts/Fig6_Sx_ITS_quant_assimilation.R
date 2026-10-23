library(ggplot2)
library(ggbeeswarm)
library(tidyverse)
library(phyloseq)
library(nlme)
library(emmeans)

# ----------------------------
# plot parameters
# ----------------------------
source("./scripts/functions/EAF-functions.R")

col_list_det <- c("orange", "blue")
col_list_det_light <- c( "#FAD88E","lightblue")
col_list_all <- c("orange", "green", "blue")
col_list_all2 <- c("green", "orange", "blue")
col_list_rhizo <- c("green", "blue")
col_list_rhizo_light <- c("lightgreen", "lightblue")

col_list_treatment <- c("darkred", "darkgreen")
colors = c("black","gray","coral4","coral",
           "chartreuse3", "darkseagreen1","blue","lightblue",
           "yellow4", "yellow","darkorchid", "plum2",
           "darkred", "darksalmon","green4",
           "greenyellow","orange",
           "moccasin",   "hotpink4", "lightpink", "lightblue4", "lightcyan3",  "lightslateblue", "lightsteelblue1", "navy",
           "darkgray", "brown", "cornflowerblue","darkgoldenrod",
           "brown3","aliceblue","aquamarine",
           "beige","bisque","blue2","blueviolet","cyan","darkblue",
           "chocolate","aquamarine3","darkcyan","deeppink")

h <- 5
w <- 6.5
res <- 300
size <- 10

Dir.o <- "./output/qSIP/fungi/"
Dir.f <- "./figures/Fig6_Sx_qSIP/fungi/"

# ----------------------------
# assumptions for absolute abundance
# ----------------------------
#high range (low copy number and high mass)
copies_per_cell_high <- 15 
pgC_per_cell_high <- 50 

#mid range (mid copy number and mid mass)
copies_per_cell_mid <- 82 
pgC_per_cell_mid <- 16.5 

#low range (high copy number and low mass)
copies_per_cell_low <-  200
pgC_per_cell_low <- 7 

# ----------------------------
# load data
# ----------------------------
ps <- readRDS("./data/qSIP/DRIPITS_phyloseq.Rds")
eaf <- read.csv("./output/qSIP/ITS_EAF_13C_w_taxonomy.csv", header = TRUE)
q_copies <- read.csv("./output/qSIP/copies_per_g_dry_soil_for_quantitative_analysis.csv")

# ============================
# 1) subset to SIP 13C samples
# ============================
ps_s <- subset_samples(ps, SIP == "SIP" & Isotope == "13C")
ps_s <- prune_taxa(taxa_sums(ps_s) > 0, ps_s)

# fraction-level raw counts
otu_s <- as(otu_table(ps_s), "matrix")
if (taxa_are_rows(ps_s)) {
  otu_s <- t(otu_s)
}

# fraction-level relative abundance
ps_s_ra <- transform_sample_counts(ps_s, function(x) x / sum(x))
otu_s_ra <- as(otu_table(ps_s_ra), "matrix")
if (taxa_are_rows(ps_s_ra)) {
  otu_s_ra <- t(otu_s_ra)
}

# taxonomy
taxa.mat <- as(tax_table(ps_s), "matrix")

# sample data
sd <- as(sample_data(ps_s), "data.frame") 

# ============================
# 2) create SampleID
# ============================
sd2 <- sd %>%
  unite("SampleID", c("sample_id", "Habitat.Isotope", "Moisture"), remove = FALSE) %>%
  mutate(
    SampleID = str_remove_all(SampleID, " "),
    SampleID = str_replace(SampleID, "Detrtitus", "Detritus"),
    avg_ITS_g_soil = as.numeric(avg_ITS_g_soil)
  )

# sample lookup tables
sd_sampleid <- sd2 %>%
  select(SAMPLE, SampleID, Isotope)

sd_fraction <- sd2 %>%
  select(SAMPLE, SampleID, Fraction)

# metadata for reconstructed tube-level abundances
data.meta <- sd2 %>%
  select(SampleID, Fraction, avg_ITS_g_soil) %>%
  group_by(SampleID) %>%
  summarise(
    avg_ITS_g_soil = sum(avg_ITS_g_soil, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  separate(SampleID, c("Tube", "Habitat", "Moisture"), sep = "_", remove = FALSE) 

# ============================
# 3) fraction-level qPCR normalization
# ============================
# proportion of total tube 16S copies found in each fraction
copies_norm <- sd2 %>%
  select(SampleID, Fraction, avg_ITS_g_soil) %>%
  tidyr::pivot_wider(names_from = Fraction, values_from = avg_ITS_g_soil) %>%
  mutate(
    copies.tube = rowSums(across(starts_with("F")), na.rm = TRUE)
  ) %>%
  mutate(
    across(
      starts_with("F"),
      ~ .x / copies.tube,
      .names = "{.col}_norm"
    )
  ) %>%
  select(SampleID, ends_with("_norm"))

#check that columns add to 1
copies_norm_sum <- copies_norm %>%
  mutate(total_ra = rowSums(across(ends_with("_norm")), na.rm = TRUE))

# ============================
# 4) fraction-level relative abundance table
# ============================
otu_s_ra_df <- otu_s_ra %>%
  as.data.frame() %>%
  rownames_to_column("SAMPLE") %>%
  left_join(sd_fraction, by = "SAMPLE")
# columns: SAMPLE, SampleID, Fraction, DRIP...

# ============================
# 5) total reads per tube
# ============================
counts_total <- otu_s %>%
  as.data.frame() %>%
  rownames_to_column("SAMPLE") %>%
  left_join(sd_sampleid, by = "SAMPLE") %>%
  mutate(
    reads_in_fraction = rowSums(across(starts_with("DRIP")), na.rm = TRUE)
  ) %>%
  group_by(SampleID) %>%
  summarise(
    total_reads = sum(reads_in_fraction, na.rm = TRUE),
    .groups = "drop"
  )

# ============================
# 6) reconstruct fraction counts using qPCR proportions
# ============================
reads_new_total <- copies_norm %>%
  left_join(counts_total, by = "SampleID") %>%
  mutate(
    across(
      ends_with("_norm"),
      ~ .x * total_reads,
      .names = "{sub('_norm$', '', .col)}"
    )
  ) %>%
  select(-contains("norm")) %>%
  pivot_longer(
    cols = starts_with("F"),
    names_to = "Fraction",
    values_to = "new_count"
  )

# ============================
# 7) reconstruct ASV counts within each fraction, then sum to tube
# ============================
norm_ASV_table <- reads_new_total %>%
  left_join(otu_s_ra_df, by = c("SampleID", "Fraction")) %>%
  mutate(
    across(starts_with("DRIP"), ~ .x * new_count)
  ) %>%
  select(-new_count, -SAMPLE, -Fraction) %>%
  group_by(SampleID) %>%
  summarise(
    across(starts_with("DRIP"), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  column_to_rownames("SampleID")

# # ============================
# # 8) make reconstructed phyloseq object
# # ============================
# # rows = taxa, cols = samples for phyloseq
# data.otu.t <- t(norm_ASV_table)
# 
# # align taxonomy to reconstructed taxa
# taxa.mat <- taxa.mat[rownames(data.otu.t), , drop = FALSE]
# 
# otu <- otu_table(data.otu.t, taxa_are_rows = TRUE)
# map <- sample_data(data.meta)
# tax <- tax_table(taxa.mat, errorIfNULL = TRUE)
# 
# new_ps <- merge_phyloseq(otu, map, tax)

# ============================
# # 9) relative abundance from reconstructed tube-level counts
# # ============================
# new_ps_ra <- transform_sample_counts(new_ps, function(x) x / sum(x))
# #new_ps_ra <- prune_taxa(taxa_sums(new_ps_ra) > 0, new_ps_ra)
# 
# new_otu_ra <- as(otu_table(new_ps_ra), "matrix")
# if (taxa_are_rows(new_ps_ra)) {
#   new_otu_ra <- t(new_otu_ra)
# }
# 
# new_sd <- as(sample_data(new_ps_ra), "data.frame") %>%
#   rownames_to_column("SampleID")

## ==============================
## 9b) calcualte new relative abundance from new count table
##===========================
new_otu_ra <- norm_ASV_table %>%
  mutate(
    row.sums = rowSums(across(starts_with("D")), na.rm = TRUE)
  ) %>%
  mutate(
    across(
      starts_with("D"),
      ~ .x / row.sums,
      .names = "{.col}_ra"  
    )
  ) %>%
  select(contains("ra")) %>%
  rename_with(~ gsub("_ra$", "", .x))




new_otu_ra_df <- new_otu_ra %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(data.meta, by = "SampleID")

meta_cols <- c("SampleID", "Tube", "Habitat", "Moisture", "avg_ITS_g_soil")
asv_cols <- setdiff(colnames(new_otu_ra_df), meta_cols)

# ============================
# 10) absolute abundance and C per ASV HIGH
# ============================
otu.18S <- new_otu_ra_df %>%
  select(-c(Habitat, Moisture)) %>%
  merge(q_copies, by.x = "Tube", by.y = "Sample_ID") %>%
  mutate(
    X18S_copies_g_dry_soil = as.numeric(X18S_copies_g_dry_soil)
  ) %>%
  mutate(
    across(all_of(asv_cols), ~ .x * X18S_copies_g_dry_soil)
  )

otu.w.cell.high <- otu.18S %>%
  mutate(
    across(all_of(asv_cols), ~ .x / copies_per_cell_high)
  )

otu.cell.c.high <- otu.w.cell.high %>%
  mutate(
    across(all_of(asv_cols), ~ .x * pgC_per_cell_high)
  )

otu.w.cell.mid <- otu.18S %>%
  mutate(
    across(all_of(asv_cols), ~ .x / copies_per_cell_mid)
  )

otu.cell.c.mid <- otu.w.cell.mid %>%
  mutate(
    across(all_of(asv_cols), ~ .x * pgC_per_cell_mid)
  )

otu.w.cell.low <- otu.18S %>%
  mutate(
    across(all_of(asv_cols), ~ .x / copies_per_cell_low)
  )

otu.cell.c.low <- otu.w.cell.low %>%
  mutate(
    across(all_of(asv_cols), ~ .x * pgC_per_cell_low)
  )

# ============================
# 11) filter EAF table to active taxa
# ============================
eaf.f <- eaf %>%
  filter(lower > 0) %>%
  mutate(
    tube = as.character(tube),
    feature_id = as.character(feature_id),
    mean_resampled_EAF = as.numeric(mean_resampled_EAF)
  ) %>%
  select(-X)

active_asvs <- intersect(unique(eaf.f$feature_id), asv_cols)

otu.cell.c.f.high <- otu.cell.c.high %>%
  select(SampleID, Tube, Habitat, Moisture, avg_ITS_g_soil, any_of(active_asvs))

otu.cell.c.f.mid <- otu.cell.c.mid %>%
  select(SampleID, Tube, Habitat, Moisture, avg_ITS_g_soil, any_of(active_asvs))

otu.cell.c.f.low <- otu.cell.c.low %>%
  select(SampleID, Tube, Habitat, Moisture, avg_ITS_g_soil, any_of(active_asvs))

otu.cell.c.long.high <- otu.cell.c.f.high %>%
  pivot_longer(
    cols = -c(SampleID, Tube, Habitat, Moisture, avg_ITS_g_soil),
    names_to = "feature_id",
    values_to = "otu_cell_c_pg"
  ) %>%
  dplyr::rename(tube = Tube) %>%
  mutate(
    tube = as.character(tube),
    feature_id = as.character(feature_id)
  )

otu.cell.c.long.mid <- otu.cell.c.f.mid %>%
  pivot_longer(
    cols = -c(SampleID, Tube, Habitat, Moisture, avg_ITS_g_soil),
    names_to = "feature_id",
    values_to = "otu_cell_c_pg"
  ) %>%
  dplyr::rename(tube = Tube) %>%
  mutate(
    tube = as.character(tube),
    feature_id = as.character(feature_id)
  )

otu.cell.c.long.low <- otu.cell.c.f.low %>%
  pivot_longer(
    cols = -c(SampleID, Tube, Habitat, Moisture, avg_ITS_g_soil),
    names_to = "feature_id",
    values_to = "otu_cell_c_pg"
  ) %>%
  dplyr::rename(tube = Tube) %>%
  mutate(
    tube = as.character(tube),
    feature_id = as.character(feature_id)
  )

# ============================
# 12) join EAF and calculate labeled C 
# ============================
final_labeled_c_high <- otu.cell.c.long.high %>%
  left_join(eaf.f, by = c("tube", "feature_id")) %>%
  mutate(
    labeled_c_pg = otu_cell_c_pg * mean_resampled_EAF
  ) %>%
  filter(!is.na(mean_resampled_EAF), !is.na(labeled_c_pg))

final_labeled_c_mid <- otu.cell.c.long.mid %>%
  left_join(eaf.f, by = c("tube", "feature_id")) %>%
  mutate(
    labeled_c_pg = otu_cell_c_pg * mean_resampled_EAF
  ) %>%
  filter(!is.na(mean_resampled_EAF), !is.na(labeled_c_pg))

final_labeled_c_low <- otu.cell.c.long.low %>%
  left_join(eaf.f, by = c("tube", "feature_id")) %>%
  mutate(
    labeled_c_pg = otu_cell_c_pg * mean_resampled_EAF
  ) %>%
  filter(!is.na(mean_resampled_EAF), !is.na(labeled_c_pg))

#summary
#count ASVs with 0% relative abundance + 

#add labeled_c column
final_labeled_c_high <- final_labeled_c_high %>%
  mutate(labeled_c = case_when(
    grepl("13C-Rhizo", Habitat) ~ "Rhizo.13C",
    grepl("13C-Detritus", Habitat) ~ "Detritus.13C",
    Habitat == "Rhizo" ~ "Rhizo.13C",
    Habitat == "Detritus" ~ "Detritus.13C"
  )) %>%
  mutate(Habitat = case_when(
    Habitat == "Rhizo" ~ "13C-Rhizo",
    Habitat == "Detritus" ~ "13C-Detritus",
    TRUE ~ Habitat
  )) %>%
  mutate(Habitat = factor(Habitat, levels =c("13C-Rhizo", "13C-Detritus", 
                                             "12C-Rhizo + 13C-Detritus",
                                             "13C-Rhizo + 12C-Detritus"))) %>%
  mutate(Family = str_remove(Family, "f__"))

final_labeled_c_mid <- final_labeled_c_mid %>%
  mutate(labeled_c = case_when(
    grepl("13C-Rhizo", Habitat) ~ "Rhizo.13C",
    grepl("13C-Detritus", Habitat) ~ "Detritus.13C",
    Habitat == "Rhizo" ~ "Rhizo.13C",
    Habitat == "Detritus" ~ "Detritus.13C"
  )) %>%
  mutate(Habitat = case_when(
    Habitat == "Rhizo" ~ "13C-Rhizo",
    Habitat == "Detritus" ~ "13C-Detritus",
    TRUE ~ Habitat
  )) %>%
  mutate(Habitat = factor(Habitat, levels =c("13C-Rhizo", "13C-Detritus", 
                                             "12C-Rhizo + 13C-Detritus",
                                             "13C-Rhizo + 12C-Detritus"))) %>%
  mutate(Family = str_remove(Family, "f__"))

final_labeled_c_low <- final_labeled_c_low %>%
  mutate(labeled_c = case_when(
    grepl("13C-Rhizo", Habitat) ~ "Rhizo.13C",
    grepl("13C-Detritus", Habitat) ~ "Detritus.13C",
    Habitat == "Rhizo" ~ "Rhizo.13C",
    Habitat == "Detritus" ~ "Detritus.13C"
  )) %>%
  mutate(Habitat = case_when(
    Habitat == "Rhizo" ~ "13C-Rhizo",
    Habitat == "Detritus" ~ "13C-Detritus",
    TRUE ~ Habitat
  )) %>%
  mutate(Habitat = factor(Habitat, levels =c("13C-Rhizo", "13C-Detritus", 
                                             "12C-Rhizo + 13C-Detritus",
                                             "13C-Rhizo + 12C-Detritus"))) %>%
  mutate(Family = str_remove(Family, "f__"))
# ============================
# 13) save
# ============================
write.csv(
  final_labeled_c_high,
  paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_ITS_HIGH.csv"),
  row.names = FALSE
)

write.csv(
  final_labeled_c_mid,
  paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_ITS_MID.csv"),
  row.names = FALSE
)

write.csv(
  final_labeled_c_low,
  paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_ITS_LOW.csv"),
  row.names = FALSE
)

# ============================
# 14) plotting subsets high
# ============================
final_labeled_c_high <- read.csv(paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_ITS_HIGH.csv"),
                                                        header = TRUE)
final_labeled_c_rhizo <- final_labeled_c_high %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_pg))



final_labeled_c_p_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Phylum) %>%
  summarise(sum = sum(labeled_c_pg))

final_labeled_c_f_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Family) %>%
  summarise(sum = sum(labeled_c_pg) )%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_asv_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_det <- final_labeled_c_high %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_pg))%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_p_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Phylum) %>%
  summarise(sum = sum(labeled_c_pg))

final_labeled_c_f_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_asv_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_g_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Genus) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))






# ============================
# 15) plots
# ============================
#rhizosphere
plot.r.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_rhizo_total,
                                         tax_col = "Family",
                                         col_values = col_list_rhizo,
                                         col_values2 = col_list_rhizo_light,
                                         top_pct = 0.7,
                                         title = "13C assimilation (pg C g soil-1)",
                                         group1 = "13C-Rhizo",                    # <-- must match exact column name after spread
                                         group2 = "13C-Rhizo + 12C-Detritus"         # <-- must match exact column name after spread
)

ggsave(
  paste(Dir.f, "Fig6-rhizo13-assimilation-family-high.png"),plot.r.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  paste(Dir.f, "Fig6-rhizo13-assimilation-family-high.pdf"),plot.r.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)

#merge with table contining asvs to quanitfy asvs within families
r.asv.per.family <- plot.r.f.per$log2fc_data %>%
  merge(final_labeled_c_asv_rhizo_total, by = "Family")


r.thel <- r.asv.per.family %>%
  filter(Family == "Thelephoraceae")

unique(r.thel$feature_id) #4 unique ASVS

r.muc <- r.asv.per.family %>%
  filter(Family == "Mucoraceae")

unique(r.muc$feature_id) #3 unique ASVS

r.lasio <- r.asv.per.family %>%
  filter(Family == "Lasiosphaeriaceae")

unique(r.lasio$feature_id) #6 unique ASVS

#plot total rhizodeposition-derived 13C assimilation across total community
plot.r <- final_labeled_c_rhizo_total %>%
  ggplot(aes(x = Habitat, y = sum, fill = Habitat)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white",
               aes(group = Habitat), position = position_dodge(0.75)) + 
  scale_fill_manual(values = col_list_rhizo) +
  facet_wrap(~Moisture) +
  scale_y_continuous(
    limits = c(0, NA),           # force zero on y axis
    breaks = scales::pretty_breaks(n = 3),  # suggest ~6 break points
    labels = scales::scientific_format()      # optional: adds comma formatting to large numbers
  ) +
  labs(y="13C assimilation (pg C cell-1)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
ggsave(
  paste(Dir.f, "FigSx-rhizo13-assimilation-total-high.png"),
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  paste(Dir.f, "FigSx-rhizo13-assimilation-total-high.pdf"),
  height = h, width = w, dpi = 300, units = "in"
)


#detritusphere
plot.d.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_det_total,
                                         tax_col = "Family",
                                         col_values = col_list_det,
                                         col_values2 = col_list_det_light,
                                         top_pct = 0.70,
                                         title = "13C assimilation (pg C g soil-1)",
                                         group1 = "13C-Detritus",                    # <-- must match exact column name after spread
                                         group2 = "12C-Rhizo + 13C-Detritus"         # <-- must match exact column name after spread
)

ggsave(
  paste(Dir.f, "Fig6-det13-assimilation-family-high.png"),plot.d.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  paste(Dir.f, "Fig6-det13-assimilation-family-high.pdf"),plot.d.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)

#merge with table contining asvs to quanitfy asvs within families
d.asv.per.family <- plot.d.f.per$log2fc_data %>%
  merge(final_labeled_c_asv_det_total, by = "Family")

d.thel <- d.asv.per.family %>%
  filter(Family == "Thelephoraceae")

unique(d.thel$feature_id) #5 unique ASVS

d.muc <- d.asv.per.family %>%
  filter(Family == "Pleosporaceae")

unique(d.muc$feature_id) #14 unique ASVS

d.lasio <- d.asv.per.family %>%
  filter(Family == "Lasiosphaeriaceae")

unique(d.lasio$feature_id) #7 unique ASVS




plot.d <- final_labeled_c_det_total %>%
  ggplot(aes(x = Habitat, y =sum, fill = Habitat)) +
  geom_boxplot() +
  facet_wrap(~Moisture) +
  # stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white",
  #              aes(group = Habitat), position = position_dodge(0.75)) + 
  scale_fill_manual(values = col_list_det) +
  scale_y_continuous(
    limits = c(0, NA),           # force zero on y axis
    breaks = scales::pretty_breaks(n = 3),  # suggest ~3 break points
    labels = scales::scientific_format()      # optional: adds comma formatting to large numbers
  ) +
  labs(y="13C assimilation (pg C cell-1)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
   # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
ggsave(
  paste(Dir.f, "FigSx-det13-assimilation-total-high.png"),
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  paste(Dir.f, "FigSx-det13-assimilation-total-high.pdf"),
  height = h, width = w, dpi = 300, units = "in"
)


# ============================
# 16) statistics on labeled C
# ============================
# Detritus * phylum
m_det <- lm(sum ~ Habitat * Moisture * Phylum, data = final_labeled_c_p_det_total)
anova(m_det)
# Response: sum
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# Habitat                  1 1.2857e+15 1.2857e+15  3.7404   0.05817 .  
# Moisture                 1 3.8553e+14 3.8553e+14  1.1216   0.29413    
# Phylum                   8 6.6645e+16 8.3306e+15 24.2350 1.406e-15 ***
#   Habitat:Moisture         1 1.8930e+15 1.8930e+15  5.5071   0.02250 *  
#   Habitat:Phylum           7 3.0776e+15 4.3965e+14  1.2790   0.27746    
# Moisture:Phylum          8 3.5683e+15 4.4604e+14  1.2976   0.26371    
# Habitat:Moisture:Phylum  5 2.8075e+15 5.6150e+14  1.6335   0.16621    
# Residuals               56 1.9250e+16 3.4374e+14                     

em_det_hab <- emmeans(m_det, ~ Habitat | Moisture | Phylum)
pairs(em_det_hab, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
# Moisture = Normal, Phylum = p__Ascomycota:
#   contrast                                     estimate       SE df t.ratio p.value
# (12C-Rhizo + 13C-Detritus) - (13C-Detritus) -72078556 15138107 56  -4.761 <0.0001

# Detritus * family
m_det_f <- lm(sum ~ Habitat * Moisture * Family, data = final_labeled_c_f_det_total)
anova(m_det_f)
# Response: sum
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# Habitat                   1 1.4614e+14 1.4614e+14  5.2665  0.022218 *  
#   Moisture                  1 6.1494e+11 6.1494e+11  0.0222  0.881731    
# Family                   86 3.5222e+16 4.0956e+14 14.7592 < 2.2e-16 ***
#   Habitat:Moisture          1 2.7140e+14 2.7140e+14  9.7802  0.001883 ** 
#   Habitat:Family           60 8.9302e+14 1.4884e+13  0.5364  0.998151    
# Moisture:Family          60 1.7647e+15 2.9412e+13  1.0599  0.363398    
# Habitat:Moisture:Family  46 1.7069e+15 3.7106e+13  1.3372  0.075762 .  
# Residuals               432 1.1988e+16 2.7750e+13                     

em_det_hab_f <- emmeans(m_det_f, ~ Habitat | Moisture | Family)
pairs(em_det_hab_f, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
# Moisture = Normal, Family = Lasiosphaeriaceae:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  19094190 4301139 432   4.439 <0.0001
# 
# Moisture = Drought, Family = Lasiosphaeriaceae:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)   8129870 4023347 432   2.021  0.0439
# 
# Moisture = Normal, Family = Pleosporaceae:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  14599435 4301139 432   3.394  0.0008
# 
# Moisture = Normal, Family = Thelephoraceae:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  11434916 4301139 432   2.659  0.0081
# 
# Moisture = Drought, Family = Thelephoraceae:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  -9390840 4023347 432  -2.334  0.0200
# 
# Moisture = Normal, Family = unassigned:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  24893409 4301139 432   5.788 <0.0001
# 
# Moisture = Drought, Family = unassigned:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) -14888520 4023347 432  -3.701  0.0002

# Detritus
m_det <- lm(sum ~ Habitat * Moisture, data = final_labeled_c_det_total)
anova(m_det)

# Response: sum
# Df     Sum Sq    Mean Sq F value Pr(>F)
# Habitat           1 5.7895e+15 5.7895e+15  1.0558 0.3310
# Moisture          1 5.6818e+13 5.6818e+13  0.0104 0.9212
# Habitat:Moisture  1 1.3837e+16 1.3837e+16  2.5234 0.1466
# Residuals         9 4.9353e+16 5.4837e+15     

em_det_hab <- emmeans(m_det, ~ Habitat | Moisture)
pairs(em_det_hab, adjust = "tukey")
# Moisture = Normal:
#   contrast                                     estimate       SE df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) 112173725 60500000  9   1.855  0.0965
# 
# Moisture = Drought:
#   contrast                                     estimate       SE df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) -19342782 56600000  9  -0.342  0.7402

em_det_moist <- emmeans(m_det, ~ Moisture | Habitat)
pairs(em_det_moist, adjust = "tukey")
# Habitat = 13C-Detritus:
#   contrast          estimate       SE df t.ratio p.value
# Normal - Drought  57170011 56600000  9   1.011  0.3385
# 
# Habitat = 12C-Rhizo + 13C-Detritus:
#   contrast          estimate       SE df t.ratio p.value
# Normal - Drought -74346496 60500000  9  -1.230  0.2500

# Rhizosphere * phylum
m_rhizo <- lm(sum ~ Habitat * Moisture * Phylum, data = final_labeled_c_p_rhizo_total)
anova(m_rhizo)
# Response: sum
# Df     Sum Sq    Mean Sq F value   Pr(>F)   
# Habitat                  1 1.4543e+12 1.4543e+12  0.0070 0.933494   
# Moisture                 1 5.2833e+14 5.2833e+14  2.5501 0.115291   
# Phylum                   9 5.5983e+15 6.2203e+14  3.0024 0.004836 **
#   Habitat:Moisture         1 3.1028e+14 3.1028e+14  1.4976 0.225595   
# Habitat:Phylum           7 4.2341e+14 6.0487e+13  0.2920 0.954697   
# Moisture:Phylum          6 6.3614e+14 1.0602e+14  0.5117 0.797239   
# Habitat:Moisture:Phylum  4 3.9508e+14 9.8770e+13  0.4767 0.752644   
# Residuals               63 1.3052e+16 2.0718e+14         

em_rhizo_hab <- emmeans(m_rhizo, ~ Habitat | Moisture | Phylum)
pairs(em_rhizo_hab, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
#none

# Rhizosphere * family
m_rhizo_f <- lm(sum ~ Habitat * Moisture * Family, data = final_labeled_c_f_rhizo_total)
anova(m_rhizo_f)
# esponse: sum
# Df     Sum Sq    Mean Sq F value  Pr(>F)  
# Habitat                   1 3.0478e+11 3.0478e+11  0.0105 0.91834  
# Moisture                  1 6.0944e+13 6.0944e+13  2.1049 0.14767  
# Family                   90 3.4553e+15 3.8392e+13  1.3260 0.03769 *
#   Habitat:Moisture          1 4.7339e+13 4.7339e+13  1.6350 0.20180  
# Habitat:Family           63 5.3021e+14 8.4161e+12  0.2907 1.00000  
# Moisture:Family          59 1.0659e+15 1.8066e+13  0.6240 0.98593  
# Habitat:Moisture:Family  38 6.4898e+14 1.7078e+13  0.5899 0.97569  
# Residuals               372 1.0771e+16 2.8953e+13       

em_rhizo_hab_f <- emmeans(m_rhizo_f, ~ Habitat | Moisture | Family)
pairs(em_rhizo_hab_f, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)

# Moisture = Drought, Family = unassigned:
#   contrast                                 estimate      SE  df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) 20990685 3804801 372   5.517 <0.0001

# Rhizosphere
m_rhizo <- lm(sum ~ Habitat * Moisture, data = final_labeled_c_rhizo_total)
 anova(m_rhizo)
 # Response: sum
 # Df     Sum Sq    Mean Sq F value Pr(>F)
 # Habitat           1 2.8309e+14 2.8309e+14  0.0882 0.7720
 # Moisture          1 4.7159e+15 4.7159e+15  1.4693 0.2508
 # Habitat:Moisture  1 1.9594e+15 1.9594e+15  0.6105 0.4511
 # Residuals        11 3.5306e+16 3.2096e+15 
 
em_rhizo_hab <- emmeans(m_rhizo, ~ Habitat | Moisture)
pairs(em_rhizo_hab, adjust = "tukey") 
# Moisture = Normal:
#   contrast                                  estimate       SE df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) -18645827 43300000 11  -0.431  0.6748
# 
# Moisture = Drought:
#   contrast                                  estimate       SE df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus)  27427294 40100000 11   0.685  0.5077

# ============================
# 14) plotting subsets mid
# ============================
final_labeled_c <- read.csv(paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_ITS_MID.csv"),
                            header = TRUE)
final_labeled_c_rhizo <- final_labeled_c %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_pg))



final_labeled_c_p_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Phylum) %>%
  summarise(sum = sum(labeled_c_pg))

final_labeled_c_f_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Family) %>%
  summarise(sum = sum(labeled_c_pg) )%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_asv_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_det <- final_labeled_c %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_pg))%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_p_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Phylum) %>%
  summarise(sum = sum(labeled_c_pg))

final_labeled_c_f_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_g_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Genus) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_asv_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))






# ============================
# 15) plots
# ============================
#rhizosphere
plot.r.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_rhizo_total,
                                         tax_col = "Family",
                                         col_values = col_list_rhizo,
                                         col_values2 = col_list_rhizo_light,
                                         top_pct = 0.7,
                                         title = "13C assimilation (pg C g soil-1)",
                                         group1 = "13C-Rhizo",                    # <-- must match exact column name after spread
                                         group2 = "13C-Rhizo + 12C-Detritus"         # <-- must match exact column name after spread
)

ggsave(
  paste(Dir.f, "Fig6-rhizo13-assimilation-family-mid.png"),plot.r.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  paste(Dir.f, "Fig6-rhizo13-assimilation-family-mid.pdf"),plot.r.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)

#merge with table contining asvs to quanitfy asvs within families
r.asv.per.family <- plot.r.f.per$log2fc_data %>%
  merge(final_labeled_c_asv_rhizo_total, by = "Family")


r.thel <- r.asv.per.family %>%
  filter(Family == "Thelephoraceae")

unique(r.thel$feature_id) #4 unique ASVS

r.muc <- r.asv.per.family %>%
  filter(Family == "Mucoraceae")

unique(r.muc$feature_id) #3 unique ASVS

r.lasio <- r.asv.per.family %>%
  filter(Family == "Lasiosphaeriaceae")

unique(r.lasio$feature_id) #6 unique ASVS

#plot total rhizodeposition-derived 13C assimilation across total community
plot.r <- final_labeled_c_rhizo_total %>%
  ggplot(aes(x = Habitat, y = sum, fill = Habitat)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white",
               aes(group = Habitat), position = position_dodge(0.75)) + 
  scale_fill_manual(values = col_list_rhizo) +
  facet_wrap(~Moisture) +
  scale_y_continuous(
    limits = c(0, NA),           # force zero on y axis
    breaks = scales::pretty_breaks(n = 3),  # suggest ~6 break points
    labels = scales::scientific_format()      # optional: adds comma formatting to large numbers
  ) +
  labs(y="13C assimilation (pg C cell-1)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
ggsave(
  paste(Dir.f, "FigSx-rhizo13-assimilation-total-mid.png"),
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  paste(Dir.f, "FigSx-rhizo13-assimilation-total-mid.pdf"),
  height = h, width = w, dpi = 300, units = "in"
)


#detritusphere
plot.d.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_det_total,
                                         tax_col = "Family",
                                         col_values = col_list_det,
                                         col_values2 = col_list_det_light,
                                         top_pct = 0.70,
                                         title = "13C assimilation (pg C g soil-1)",
                                         group1 = "13C-Detritus",                    # <-- must match exact column name after spread
                                         group2 = "12C-Rhizo + 13C-Detritus"         # <-- must match exact column name after spread
)

ggsave(
  paste(Dir.f, "Fig6-det13-assimilation-family-mid.png"),plot.d.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  paste(Dir.f, "Fig6-det13-assimilation-family-mid.pdf"),plot.d.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)

#merge with table contining asvs to quanitfy asvs within families
d.asv.per.family <- plot.d.f.per$log2fc_data %>%
  merge(final_labeled_c_asv_det_total, by = "Family")

d.thel <- d.asv.per.family %>%
  filter(Family == "Thelephoraceae")

unique(d.thel$feature_id) #5 unique ASVS

d.muc <- d.asv.per.family %>%
  filter(Family == "Pleosporaceae")

unique(d.muc$feature_id) #14 unique ASVS

d.lasio <- d.asv.per.family %>%
  filter(Family == "Lasiosphaeriaceae")

unique(d.lasio$feature_id) #7 unique ASVS




plot.d <- final_labeled_c_det_total %>%
  ggplot(aes(x = Habitat, y =sum, fill = Habitat)) +
  geom_boxplot() +
  facet_wrap(~Moisture) +
  # stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white",
  #              aes(group = Habitat), position = position_dodge(0.75)) + 
  scale_fill_manual(values = col_list_det) +
  scale_y_continuous(
    limits = c(0, NA),           # force zero on y axis
    breaks = scales::pretty_breaks(n = 3),  # suggest ~3 break points
    labels = scales::scientific_format()      # optional: adds comma formatting to large numbers
  ) +
  labs(y="13C assimilation (pg C cell-1)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
ggsave(
  paste(Dir.f, "FigSx-det13-assimilation-total-mid.png"),
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  paste(Dir.f, "FigSx-det13-assimilation-total-mid.pdf"),
  height = h, width = w, dpi = 300, units = "in"
)

# ============================
# 16) statistics on labeled C
# ============================
#same as above

# ============================
# 14) plotting subsets low
# ============================
final_labeled_c <- read.csv(paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_ITS_LOW.csv"),
                            header = TRUE)
final_labeled_c_rhizo <- final_labeled_c %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_pg))



final_labeled_c_p_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Phylum) %>%
  summarise(sum = sum(labeled_c_pg))

final_labeled_c_f_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Family) %>%
  summarise(sum = sum(labeled_c_pg) )%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_asv_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_det <- final_labeled_c %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_pg))%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_p_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Phylum) %>%
  summarise(sum = sum(labeled_c_pg))

final_labeled_c_f_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_g_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, Genus) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_asv_det_total <- final_labeled_c_det %>%
  filter(labeled_c_pg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, Family) %>%
  summarise(sum = sum(labeled_c_pg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))






# ============================
# 15) plots
# ============================
#rhizosphere
plot.r.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_rhizo_total,
                                         tax_col = "Family",
                                         col_values = col_list_rhizo,
                                         col_values2 = col_list_rhizo_light,
                                         top_pct = 0.7,
                                         title = "13C assimilation (pg C g soil-1)",
                                         group1 = "13C-Rhizo",                    # <-- must match exact column name after spread
                                         group2 = "13C-Rhizo + 12C-Detritus"         # <-- must match exact column name after spread
)

ggsave(
  paste(Dir.f, "Fig6-rhizo13-assimilation-family-low.png"),plot.r.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  paste(Dir.f, "Fig6-rhizo13-assimilation-family-low.pdf"),plot.r.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)

#merge with table contining asvs to quanitfy asvs within families
r.asv.per.family <- plot.r.f.per$log2fc_data %>%
  merge(final_labeled_c_asv_rhizo_total, by = "Family")


r.thel <- r.asv.per.family %>%
  filter(Family == "Thelephoraceae")

unique(r.thel$feature_id) #4 unique ASVS

r.muc <- r.asv.per.family %>%
  filter(Family == "Mucoraceae")

unique(r.muc$feature_id) #3 unique ASVS

r.lasio <- r.asv.per.family %>%
  filter(Family == "Lasiosphaeriaceae")

unique(r.lasio$feature_id) #6 unique ASVS

#plot total rhizodeposition-derived 13C assimilation across total community
plot.r <- final_labeled_c_rhizo_total %>%
  ggplot(aes(x = Habitat, y = sum, fill = Habitat)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white",
               aes(group = Habitat), position = position_dodge(0.75)) + 
  scale_fill_manual(values = col_list_rhizo) +
  facet_wrap(~Moisture) +
  scale_y_continuous(
    limits = c(0, NA),           # force zero on y axis
    breaks = scales::pretty_breaks(n = 3),  # suggest ~6 break points
    labels = scales::scientific_format()      # optional: adds comma formatting to large numbers
  ) +
  labs(y="13C assimilation (pg C cell-1)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
   # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
ggsave(
  paste(Dir.f, "FigSx-rhizo13-assimilation-total-low.png"),
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  paste(Dir.f, "FigSx-rhizo13-assimilation-total-low.pdf"),
  height = h, width = w, dpi = 300, units = "in"
)


#detritusphere
plot.d.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_det_total,
                                         tax_col = "Family",
                                         col_values = col_list_det,
                                         col_values2 = col_list_det_light,
                                         top_pct = 0.70,
                                         title = "13C assimilation (pg C g soil-1)",
                                         group1 = "13C-Detritus",                    # <-- must match exact column name after spread
                                         group2 = "12C-Rhizo + 13C-Detritus"         # <-- must match exact column name after spread
)

ggsave(
  paste(Dir.f, "Fig6-det13-assimilation-family-low.png"),plot.d.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  paste(Dir.f, "Fig6-det13-assimilation-family-low.pdf"),plot.d.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)

#merge with table contining asvs to quanitfy asvs within families
d.asv.per.family <- plot.d.f.per$log2fc_data %>%
  merge(final_labeled_c_asv_det_total, by = "Family")

d.thel <- d.asv.per.family %>%
  filter(Family == "Thelephoraceae")

unique(d.thel$feature_id) #5 unique ASVS

d.muc <- d.asv.per.family %>%
  filter(Family == "Pleosporaceae")

unique(d.muc$feature_id) #14 unique ASVS

d.lasio <- d.asv.per.family %>%
  filter(Family == "Lasiosphaeriaceae")

unique(d.lasio$feature_id) #7 unique ASVS




plot.d <- final_labeled_c_det_total %>%
  ggplot(aes(x = Habitat, y =sum, fill = Habitat)) +
  geom_boxplot() +
  facet_wrap(~Moisture) +
  # stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white",
  #              aes(group = Habitat), position = position_dodge(0.75)) + 
  scale_fill_manual(values = col_list_det) +
  scale_y_continuous(
    limits = c(0, NA),           # force zero on y axis
    breaks = scales::pretty_breaks(n = 3),  # suggest ~3 break points
    labels = scales::scientific_format()      # optional: adds comma formatting to large numbers
  ) +
  labs(y="13C assimilation (pg C cell-1)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
ggsave(
  paste(Dir.f, "FigSx-det13-assimilation-total-low.png"),
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  paste(Dir.f, "FigSx-det13-assimilation-total-low.pdf"),
  height = h, width = w, dpi = 300, units = "in"
)

# ============================
# 16) statistics on labeled C
# ============================
#same as above
