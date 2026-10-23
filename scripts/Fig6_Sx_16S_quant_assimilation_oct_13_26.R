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

Dir.o <- "./output/qSIP/"
Dir.f <- "./figures/Fig6_Sx_qSIP/"

# ----------------------------
# assumptions for absolute abundance
# ----------------------------
#high range (low copy number and high mass)
copies_per_cell_high <- 1
fgC_per_cell_high <- 100 


#low range (high copy number and low mass)
copies_per_cell_low <-  6
fgC_per_cell_low <- 10 

# ----------------------------
# load data
# ----------------------------
ps <- readRDS("./data/qSIP/DRIP16S_phyloseq.Rds")
eaf <- read.csv("./output/qSIP/16S_EAF_13C_w_taxonomy.csv", header = TRUE)
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
    avg_16S_g_soil = as.numeric(avg_16S_g_soil)
  )

# sample lookup tables
sd_sampleid <- sd2 %>%
  select(SAMPLE, SampleID, Isotope)

sd_fraction <- sd2 %>%
  select(SAMPLE, SampleID, Fraction)

# metadata for reconstructed tube-level abundances
data.meta <- sd2 %>%
  select(SampleID, Fraction, avg_16S_g_soil) %>%
  group_by(SampleID) %>%
  summarise(
    avg_16S_g_soil = sum(avg_16S_g_soil, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  separate(SampleID, c("Tube", "Habitat", "Moisture"), sep = "_", remove = FALSE) 

# ============================
# 3) fraction-level qPCR normalization
# ============================
# proportion of total tube 16S copies found in each fraction
copies_norm <- sd2 %>%
  select(SampleID, Fraction, avg_16S_g_soil) %>%
  tidyr::pivot_wider(names_from = Fraction, values_from = avg_16S_g_soil) %>%
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

meta_cols <- c("SampleID", "Tube", "Habitat", "Moisture", "avg_16S_g_soil")
asv_cols <- setdiff(colnames(new_otu_ra_df), meta_cols)

# ============================
# 10) absolute abundance and C per ASV HIGH
# ============================
otu.16S <- new_otu_ra_df %>%
  select(-c(Habitat, Moisture)) %>%
  merge(q_copies, by.x = "Tube", by.y = "Sample_ID") %>%
  mutate(
    X16S_copies_g_dry_soil = as.numeric(X16S_copies_g_dry_soil)
  ) %>%
  mutate(
    across(all_of(asv_cols), ~ .x * X16S_copies_g_dry_soil)
  )

otu.w.cell.high <- otu.16S %>%
  mutate(
    across(all_of(asv_cols), ~ .x / copies_per_cell_high)
  )

otu.cell.c.high <- otu.w.cell.high %>%
  mutate(
    across(all_of(asv_cols), ~ .x * fgC_per_cell_high)
  )


otu.w.cell.low <- otu.16S %>%
  mutate(
    across(all_of(asv_cols), ~ .x / copies_per_cell_low)
  )

otu.cell.c.low <- otu.w.cell.low %>%
  mutate(
    across(all_of(asv_cols), ~ .x * fgC_per_cell_low)
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
  select(SampleID, Tube, Habitat, Moisture, avg_16S_g_soil, any_of(active_asvs))


otu.cell.c.f.low <- otu.cell.c.low %>%
  select(SampleID, Tube, Habitat, Moisture, avg_16S_g_soil, any_of(active_asvs))

otu.cell.c.long.high <- otu.cell.c.f.high %>%
  pivot_longer(
    cols = -c(SampleID, Tube, Habitat, Moisture, avg_16S_g_soil),
    names_to = "feature_id",
    values_to = "otu_cell_c_fg"
  ) %>%
  dplyr::rename(tube = Tube) %>%
  mutate(
    tube = as.character(tube),
    feature_id = as.character(feature_id)
  )


otu.cell.c.long.low <- otu.cell.c.f.low %>%
  pivot_longer(
    cols = -c(SampleID, Tube, Habitat, Moisture, avg_16S_g_soil),
    names_to = "feature_id",
    values_to = "otu_cell_c_fg"
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
    labeled_c_fg = otu_cell_c_fg * mean_resampled_EAF
  ) %>%
  filter(!is.na(mean_resampled_EAF), !is.na(labeled_c_fg))


final_labeled_c_low <- otu.cell.c.long.low %>%
  left_join(eaf.f, by = c("tube", "feature_id")) %>%
  mutate(
    labeled_c_fg = otu_cell_c_fg * mean_resampled_EAF
  ) %>%
  filter(!is.na(mean_resampled_EAF), !is.na(labeled_c_fg))

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
  mutate(family = str_remove(family, "f__"))


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
  mutate(family = str_remove(family, "f__"))
# ============================
# 13) save
# ============================
write.csv(
  final_labeled_c_high,
  paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_16S_HIGH.csv"),
  row.names = FALSE
)


write.csv(
  final_labeled_c_low,
  paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_16S_LOW.csv"),
  row.names = FALSE
)

# ============================
# 14) plotting subsets high
# ============================
final_labeled_c_high <- read.csv(paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_16S_HIGH.csv"),
                                                        header = TRUE)
final_labeled_c_rhizo <- final_labeled_c_high %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_fg))



final_labeled_c_p_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, phylum) %>%
  summarise(sum = sum(labeled_c_fg))

final_labeled_c_f_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, family) %>%
  summarise(sum = sum(labeled_c_fg) )%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_asv_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, family, genus) %>%
  summarise(sum = sum(labeled_c_fg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_det <- final_labeled_c_high %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_fg))%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_p_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, phylum) %>%
  summarise(sum = sum(labeled_c_fg))

final_labeled_c_f_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, family) %>%
  summarise(sum = sum(labeled_c_fg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_asv_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, family, genus) %>%
  summarise(sum = sum(labeled_c_fg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_g_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, genus) %>%
  summarise(sum = sum(labeled_c_fg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))






# ============================
# 15) plots
# ============================
#rhizosphere
plot.r.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_rhizo_total,
                                         tax_col = "family",
                                         col_values = col_list_rhizo,
                                         col_values2 = col_list_rhizo_light,
                                         top_pct = 0.7,
                                         title = "13C assimilation (fg C g soil-1)",
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
  merge(final_labeled_c_asv_rhizo_total, by = "family")


r.bac <- r.asv.per.family %>%
  filter(family == "Bacillaceae 1")

unique(r.bac$feature_id) #11 unique ASVS
unique(r.bac$genus) #11 unique ASVS
r.bac.sort <- r.bac %>%
  arrange(desc(sum)) #neobacillus and bacillus

r.strep <- r.asv.per.family %>%
  filter(family == "Streptomycetaceae")

unique(r.strep$feature_id) #17 unique ASVS
unique(r.strep$genus)
r.strep.sort <- r.strep %>%
  arrange(desc(sum)) #neobacillus and bacillus


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
  labs(y="13C assimilation (fg C cell-1)") +
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
                                         tax_col = "family",
                                         col_values = col_list_det,
                                         col_values2 = col_list_det_light,
                                         top_pct = 0.70,
                                         title = "13C assimilation (fg C g soil-1)",
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
  merge(final_labeled_c_asv_det_total, by = "family")

r.bac <- d.asv.per.family %>%
  filter(family == "Bacillaceae 1")

unique(r.bac$feature_id) #11 unique ASVS
unique(r.bac$genus) #11 unique ASVS
r.bac.sort <- r.bac %>%
  arrange(desc(sum)) #neobacillus and bacillus

r.strep <- d.asv.per.family %>%
  filter(family == "Streptomycetaceae")

unique(r.strep$feature_id) #17 unique ASVS
unique(r.strep$genus)
r.strep.sort <- r.strep %>%
  arrange(desc(sum)) #neobacillus and bacillus



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
  labs(y="13C assimilation (fg C cell-1)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.2),
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
m_det <- lm(sum ~ Habitat * Moisture * phylum, data = final_labeled_c_p_det_total)
anova(m_det)
# Response: sum
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# Habitat                   1 1.4426e+13 1.4426e+13  0.9114 0.3418033    
# Moisture                  1 1.4464e+13 1.4464e+13  0.9138 0.3411669    
# phylum                   18 1.1970e+16 6.6498e+14 42.0116 < 2.2e-16 ***
#   Habitat:Moisture          1 1.5671e+11 1.5671e+11  0.0099 0.9209190    
# Habitat:phylum           13 6.1595e+14 4.7381e+13  2.9934 0.0008522 ***
#   Moisture:phylum          13 3.3788e+13 2.5991e+12  0.1642 0.9995969    
# Habitat:Moisture:phylum  12 1.1789e+13 9.8246e+11  0.0621 0.9999966    
# Residuals               112 1.7728e+15 1.5828e+13                        

em_det_hab <- emmeans(m_det, ~ Habitat | Moisture | phylum)
pairs(em_det_hab, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
# Moisture = Normal, phylum = Actinobacteria:
#   contrast                                     estimate      SE  df t.ratio p.value
# (12C-Rhizo + 13C-Detritus) - (13C-Detritus) -10508309 3248426 112  -3.235  0.0016
# 
# Moisture = Drought, phylum = Actinobacteria:
#   contrast                                     estimate      SE  df t.ratio p.value
# (12C-Rhizo + 13C-Detritus) - (13C-Detritus) -12404235 3038624 112  -4.082 <0.0001
# 
# Moisture = Normal, phylum = Firmicutes:
#   contrast                                     estimate      SE  df t.ratio p.value
# (12C-Rhizo + 13C-Detritus) - (13C-Detritus)   8011973 3248426 112   2.466  0.0152
# 
# Moisture = Drought, phylum = Firmicutes:
#   contrast                                     estimate      SE  df t.ratio p.value
# (12C-Rhizo + 13C-Detritus) - (13C-Detritus)   7014407 3038624 112   2.308  0.0228

# Detritus * family
m_det_f <- lm(sum ~ Habitat * Moisture * family, data = final_labeled_c_f_det_total)
anova(m_det_f)
# Response: sum
# Df     Sum Sq    Mean Sq F value  Pr(>F)    
# Habitat                   1 1.2161e+12 1.2161e+12  0.9418 0.33210    
# Moisture                  1 6.8635e+12 6.8635e+12  5.3150 0.02138 *  
#   family                  160 6.6002e+15 4.1251e+13 31.9446 < 2e-16 ***
#   Habitat:Moisture          1 5.8174e+10 5.8174e+10  0.0450 0.83196    
# Habitat:family          128 4.9956e+14 3.9028e+12  3.0223 < 2e-16 ***
#   Moisture:family         118 5.4436e+13 4.6132e+11  0.3572 1.00000    
# Habitat:Moisture:family  80 8.1550e+12 1.0194e+11  0.0789 1.00000    
# Residuals               848 1.0951e+15 1.2913e+12                  

em_det_hab_f <- emmeans(m_det_f, ~ Habitat | Moisture | family)
pairs(em_det_hab_f, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
# Moisture = Normal, family = Bacillaceae 1:
#   contrast                                    estimate       SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) -7389589 927842.2 848  -7.964 <0.0001
# 
# Moisture = Drought, family = Bacillaceae 1:
#   contrast                                    estimate       SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) -7025244 867916.9 848  -8.094 <0.0001
# 
# Moisture = Normal, family = Enterobacteriaceae:
#   contrast                                    estimate       SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  1868501 927842.2 848   2.014  0.0443
# 
# Moisture = Normal, family = Streptomycetaceae:
#   contrast                                    estimate       SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  9875953 927842.2 848  10.644 <0.0001
# 
# Moisture = Drought, family = Streptomycetaceae:
#   contrast                                    estimate       SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) 10228448 867916.9 848  11.785 <0.0001

# Detritus
m_det <- lm(sum ~ Habitat * Moisture, data = final_labeled_c_det_total)
anova(m_det)

# Response: sum
# Df     Sum Sq    Mean Sq F value Pr(>F)
# Habitat           1 1.2648e+12 1.2648e+12  0.0031 0.9570
# Moisture          1 3.8077e+14 3.8077e+14  0.9243 0.3615
# Habitat:Moisture  1 8.6505e+11 8.6505e+11  0.0021 0.9645
# Residuals         9 3.7077e+15 4.1196e+14      

em_det_hab <- emmeans(m_det, ~ Habitat | Moisture)
pairs(em_det_hab, adjust = "tukey")
# Moisture = Normal:
#   contrast                                    estimate       SE df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  1957708 16600000  9   0.118  0.9086
# 
# Moisture = Drought:
#   contrast                                    estimate       SE df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)   917848 15500000  9   0.059  0.9541

em_det_moist <- emmeans(m_det, ~ Moisture | Habitat)
pairs(em_det_moist, adjust = "tukey")
# Habitat = 13C-Detritus:
#   contrast         estimate       SE df t.ratio p.value
# Normal - Drought 11369263 15500000  9   0.733  0.4820
# 
# Habitat = 12C-Rhizo + 13C-Detritus:
#   contrast         estimate       SE df t.ratio p.value
# Normal - Drought 10329402 16600000  9   0.623  0.5486

# Rhizosphere * phylum
m_rhizo <- lm(sum ~ Habitat * Moisture * phylum, data = final_labeled_c_p_rhizo_total)
anova(m_rhizo)
# Response: sum
# Df     Sum Sq    Mean Sq F value  Pr(>F)    
# Habitat                   1 1.7773e+11 1.7773e+11  0.0092 0.92372    
# Moisture                  1 1.0462e+13 1.0462e+13  0.5414 0.46294    
# phylum                   16 3.4379e+15 2.1487e+14 11.1184 < 2e-16 ***
#   Habitat:Moisture          1 8.2525e+13 8.2525e+13  4.2702 0.04039 *  
#   Habitat:phylum           16 3.6699e+13 2.2937e+12  0.1187 0.99999    
# Moisture:phylum          16 1.8366e+14 1.1479e+13  0.5940 0.88524    
# Habitat:Moisture:phylum  14 3.6806e+14 2.6290e+13  1.3604 0.17837    
# Residuals               161 3.1114e+15 1.9326e+13        

em_rhizo_hab <- emmeans(m_rhizo, ~ Habitat | Moisture | phylum)
pairs(em_rhizo_hab, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
# Moisture = Normal, phylum = Actinobacteria:
#   contrast                                 estimate      SE  df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) -8697045 3357565 161  -2.590  0.0105
# 
# Moisture = Drought, phylum = Actinobacteria:
#   contrast                                 estimate      SE  df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus)  9399467 3108501 161   3.024  0.0029

# Rhizosphere * family
m_rhizo_f <- lm(sum ~ Habitat * Moisture * family, data = final_labeled_c_f_rhizo_total)
anova(m_rhizo_f)
# Response: sum
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# Habitat                    1 1.4671e+11 1.4671e+11  0.0936  0.759688    
# Moisture                   1 1.0606e+12 1.0606e+12  0.6768  0.410896    
# family                   165 1.3531e+15 8.2004e+12  5.2328 < 2.2e-16 ***
#   Habitat:Moisture           1 1.5435e+13 1.5435e+13  9.8489  0.001748 ** 
#   Habitat:family           136 2.8415e+13 2.0893e+11  0.1333  1.000000    
# Moisture:family          127 1.4206e+14 1.1186e+12  0.7138  0.991577    
# Habitat:Moisture:family   93 3.7583e+14 4.0412e+12  2.5787 6.364e-13 ***
#   Residuals               1026 1.6079e+15 1.5671e+12    

em_rhizo_hab_f <- emmeans(m_rhizo_f, ~ Habitat | Moisture | family)
pairs(em_rhizo_hab_f, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)

# Moisture = Normal, family = Bacillaceae 1:
#   contrast                                 estimate       SE   df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) -6169061 956118.1 1026  -6.452 <0.0001
# 
# Moisture = Drought, family = Bacillaceae 1:
#   contrast                                 estimate       SE   df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus)  9215992 956118.1 1026   9.639 <0.0001
# 
# Moisture = Normal, family = Streptomycetaceae:
#   contrast                                 estimate       SE   df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) -7275989 956118.1 1026  -7.610 <0.0001
# 
# Moisture = Drought, family = Streptomycetaceae:
#   contrast                                 estimate       SE   df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus)  6519197 885193.4 1026   7.365 <0.0001

# Rhizosphere
m_rhizo <- lm(sum ~ Habitat * Moisture, data = final_labeled_c_rhizo_total)
 anova(m_rhizo)
 # Response: sum
 # Df     Sum Sq    Mean Sq F value Pr(>F)
 # Habitat           1 4.9483e+12 4.9483e+12  0.0058 0.9405
 # Moisture          1 8.2802e+13 8.2802e+13  0.0976 0.7606
 # Habitat:Moisture  1 1.1915e+15 1.1915e+15  1.4039 0.2610
 # Residuals        11 9.3359e+15 8.4872e+14   
 
em_rhizo_hab <- emmeans(m_rhizo, ~ Habitat | Moisture)
pairs(em_rhizo_hab, adjust = "tukey") 
# Moisture = Normal:
#   contrast                                  estimate       SE df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) -17856950 22300000 11  -0.803  0.4392
# 
# Moisture = Drought:
#   contrast                                  estimate       SE df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus)  18070375 20600000 11   0.877  0.3991



# ============================
# 14) plotting subsets low
# ============================
final_labeled_c <- read.csv(paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions_16S_LOW.csv"),
                            header = TRUE)
final_labeled_c_rhizo <- final_labeled_c %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_fg))



final_labeled_c_p_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, phylum) %>%
  summarise(sum = sum(labeled_c_fg))

final_labeled_c_f_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, family) %>%
  summarise(sum = sum(labeled_c_fg) )%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_asv_rhizo_total <- final_labeled_c_rhizo %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Rhizo.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, family) %>%
  summarise(sum = sum(labeled_c_fg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))

final_labeled_c_det <- final_labeled_c %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought")))

final_labeled_c_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture) %>%
  summarise(sum = sum(labeled_c_fg))%>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_p_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, phylum) %>%
  summarise(sum = sum(labeled_c_fg))

final_labeled_c_f_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, family) %>%
  summarise(sum = sum(labeled_c_fg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_g_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, genus) %>%
  summarise(sum = sum(labeled_c_fg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Detritus", 
                                              "12C-Rhizo + 13C-Detritus")))

final_labeled_c_asv_det_total <- final_labeled_c_det %>%
  filter(labeled_c_fg > 0) %>%
  filter(grepl("Detritus.13C", labeled_c)) %>%
  group_by(SampleID, Habitat, Moisture, feature_id, family) %>%
  summarise(sum = sum(labeled_c_fg)) %>%
  mutate(Habitat = factor(Habitat, levels = c("13C-Rhizo", 
                                              "13C-Rhizo + 12C-Detritus")))






# ============================
# 15) plots
# ============================
#rhizosphere
plot.r.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_rhizo_total,
                                         tax_col = "family",
                                         col_values = col_list_rhizo,
                                         col_values2 = col_list_rhizo_light,
                                         top_pct = 0.7,
                                         title = "13C assimilation (fg C g soil-1)",
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
  labs(y="13C assimilation (fg C cell-1)") +
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
  paste(Dir.f, "FigSx-rhizo13-assimilation-total-low.png"),
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  paste(Dir.f, "FigSx-rhizo13-assimilation-total-low.pdf"),
  height = h, width = w, dpi = 300, units = "in"
)


#detritusphere
plot.d.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_det_total,
                                         tax_col = "family",
                                         col_values = col_list_det,
                                         col_values2 = col_list_det_light,
                                         top_pct = 0.70,
                                         title = "13C assimilation (fg C g soil-1)",
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
  labs(y="13C assimilation (fg C cell-1)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.2),
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
