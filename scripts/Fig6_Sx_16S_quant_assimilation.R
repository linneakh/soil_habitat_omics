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
copies_per_cell <- 6
fgC_per_cell <- 10

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
# 10) absolute abundance and C per ASV
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

otu.w.cell <- otu.16S %>%
  mutate(
    across(all_of(asv_cols), ~ .x / copies_per_cell)
  )

otu.cell.c <- otu.w.cell %>%
  mutate(
    across(all_of(asv_cols), ~ .x * fgC_per_cell)
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

otu.cell.c.f <- otu.cell.c %>%
  select(SampleID, Tube, Habitat, Moisture, avg_16S_g_soil, any_of(active_asvs))

otu.cell.c.long <- otu.cell.c.f %>%
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
final_labeled_c <- otu.cell.c.long %>%
  left_join(eaf.f, by = c("tube", "feature_id")) %>%
  mutate(
    labeled_c_fg = otu_cell_c_fg * mean_resampled_EAF
  ) %>%
  filter(!is.na(mean_resampled_EAF), !is.na(labeled_c_fg))

#summary
#count ASVs with 0% relative abundance + 

#add labeled_c column
final_labeled_c <- final_labeled_c %>%
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
                                             "13C-Rhizo + 12C-Detritus")))

# ============================
# 13) save
# ============================
# write.csv(
#   final_labeled_c,
#   paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions.csv"),
#   row.names = FALSE
# )

# ============================
# 14) plotting subsets
# ============================
final_labeled_c <- read.csv(paste(Dir.o, "ASV_labeled_C_fg_per_sample_from_fractions.csv"),
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






# ============================
# 15) plots
# ============================
plot.r.p <- plot_top_boxplot_r_d(final_labeled_c_p_rhizo_total,
                                 tax_col = "phylum",
                                 col_values = col_list_rhizo,
                                 min_median = 100000,
                                 title="13C assimilation (fg C cell-1)")
ggsave(
  "./figures/qSIP/16S/rhizo13-assimilation-phyla-chat.png",
  height = h, width = w + 1, dpi = 300, units = "in"
)

plot.r.f <- plot_top_boxplot_r_d(final_labeled_c_f_rhizo_total,
                                 tax_col = "family",
                                 col_values = col_list_rhizo,
                                 min_median = 1000000,
                                 title="13C assimilation (fg C cell-1)")
ggsave(
  "./figures/qSIP/16S/rhizo13-assimilation-family-chat.png",
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  "./figures/qSIP/16S/rhizo13-assimilation-family-chat.pdf",
  height = h+2, width = w +1, dpi = 300, units = "in"
)

plot.r.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_rhizo_total,
                                         tax_col = "family",
                                         col_values = col_list_rhizo,
                                         col_values2 = col_list_rhizo_light,
                                         top_pct = 0.7,
                                         title = "13C assimilation (pg C g soil-1)",
                                         group1 = "13C-Rhizo",                    # <-- must match exact column name after spread
                                         group2 = "13C-Rhizo + 12C-Detritus"         # <-- must match exact column name after spread
)

ggsave(
  "./figures/qSIP/16S/rhizo13-assimilation-family-chat-per-70.png",plot.r.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  "./figures/qSIP/16S/rhizo13-assimilation-family-chat-per-70.pdf",plot.r.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)

#merge with table contining asvs to quanitfy asvs within families
r.asv.per.family <- plot.r.f.per$log2fc_data %>%
  merge(final_labeled_c_asv_det_total, by = "family")


r.strep<- r.asv.per.family %>%
  filter(family == "Streptomycetaceae")

r.firm <- r.asv.per.family %>%
  filter(family == "Bacillaceae 1")

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
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
ggsave(
  "./figures/qSIP/16S/rhizo13-assimilation-chat.png",
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  "./figures/qSIP/16S/rhizo13-assimilation-chat.pdf",
  height = h, width = w, dpi = 300, units = "in"
)

plot.d.p <- plot_top_boxplot_r_d(final_labeled_c_p_det_total,
                                 tax_col = "phylum",
                                 col_values = col_list_det,
                                 min_median = 100000,
                                 title="13C assimilation (fg C cell-1)")
ggsave(
  "./figures/qSIP/16S/det13-assimilation-phyla-chat.png",
  height = h, width = w + 1, dpi = 300, units = "in"
)

plot.d.f <- plot_top_boxplot_r_d(final_labeled_c_f_det_total,
                                 tax_col = "family",
                                 col_values = col_list_det,
                                 min_median = 1000000,
                                 title="13C assimilation (fg C cell-1)")
ggsave(
  "./figures/qSIP/16S/det13-assimilation-family-chat.png",
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  "./figures/qSIP/16S/det13-assimilation-family-chat.pdf",
  height =  h+2, width = w+1 , dpi = 300, units = "in"
)

plot.d.f.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_f_det_total,
                                         tax_col = "family",
                                         col_values = col_list_det,
                                         col_values2 = col_list_det_light,
                                         top_pct = 0.70,
                                         title = "13C assimilation (pg C g soil-1)",
                                         group1 = "13C-Detritus",                    # <-- must match exact column name after spread
                                         group2 = "12C-Rhizo + 13C-Detritus"         # <-- must match exact column name after spread
)

ggsave(
  "./figures/qSIP/16S/det13-assimilation-family-chat-per-70.png",plot.d.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)
ggsave(
  "./figures/qSIP/16S/det13-assimilation-family-chat-per-70.pdf",plot.d.f.per$boxplot,
  height = h+2, width = w +1, dpi = 300, units = "in"
)

#merge with table contining asvs to quanitfy asvs within families
d.asv.per.family <- plot.d.f.per$log2fc_data %>%
  merge(final_labeled_c_asv_det_total, by = "family")


d.strep<- d.asv.per.family %>%
  filter(family == "Streptomycetaceae")

d.firm <- d.asv.per.family %>%
  filter(family == "Bacillaceae 1")


# plot.d.asv.per <- plot_top_boxplot_r_d_per(data = final_labeled_c_asv_det_total,
#                                          tax_col = "feature_id",
#                                          col_values = col_list_det,
#                                          col_values2 = col_list_det_light,
#                                          top_pct = 0.70,
#                                          title = "13C assimilation (pg C g soil-1)",
#                                          group1 = "13C-Detritus",                    # <-- must match exact column name after spread
#                                          group2 = "12C-Rhizo + 13C-Detritus"         # <-- must match exact column name after spread
# )
# 
# #merge ASV data frame with family
# d.asv.per.family <- plot.d.asv.per$log2fc_data %>%
#   merge(final_labeled_c_asv_det_total, by = "feature_id")
# 
# d.strep<- d.asv.per.family %>%
#   filter(family == "Streptomycetaceae")
# 
# d.firm <- d.asv.per.family %>%
#   filter(family == "Bacillaceae 1")

#plot ratio pos to neg
plot.ratio.det <- log2fc_f_det %>%
  mutate(direction = case_when(
    log2fc > 0 ~ "pos",
    log2fc < 0 ~ "neg")) %>%
  mutate(family = reorder_within(family, log2fc, Moisture)) %>%  # order within each facet
  ggplot(aes(x = family, y = log2fc, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_x_reordered() +  # removes the "___Moisture" suffix from labels
  facet_wrap(~ Moisture, scales = "free_y")  # free

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
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
ggsave(
  "./figures/qSIP/16S/det13-assimilation-chat.png",
  height = h, width = w, dpi = 300, units = "in"
)

ggsave(
  "./figures/qSIP/16S/det13-assimilation-chat.pdf",
  height = h, width = w, dpi = 300, units = "in"
)


# ============================
# 16) statistics on labeled C
# ============================


# Detritus * phylum
m_det <- lm(sum ~ Habitat * Moisture * phylum, data = final_labeled_c_p_det_total)
anova(m_det)
# Response: sum
# Df     Sum Sq    Mean Sq F value  Pr(>F)    
# Habitat                   1 1.0716e+14 1.0716e+14  0.3208 0.57224    
# Moisture                  1 2.6579e+13 2.6579e+13  0.0796 0.77839    
# phylum                   18 3.5267e+17 1.9593e+16 58.6593 < 2e-16 ***
#   Habitat:Moisture          1 1.0495e+12 1.0495e+12  0.0031 0.95540    
# Habitat:phylum           13 6.8883e+15 5.2987e+14  1.5864 0.09937 .  
# Moisture:phylum          13 3.5522e+14 2.7324e+13  0.0818 0.99999    
# Habitat:Moisture:phylum  12 1.2339e+15 1.0282e+14  0.3078 0.98686    
# Residuals               112 3.7409e+16 3.3401e+14             

em_det_hab <- emmeans(m_det, ~ Habitat | Moisture | phylum)
pairs(em_det_hab, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
# Moisture = Normal, phylum = Actinobacteria:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  13135386 4060533 112   3.235  0.0016
# 
# Moisture = Drought, phylum = Actinobacteria:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  15505293 3798280 112   4.082 <0.0001
# 
# Moisture = Normal, phylum = Firmicutes:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) -10014967 4060533 112  -2.466  0.0152
# 
# Moisture = Drought, phylum = Firmicutes:
#   contrast                                     estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  -8768008 3798280 112  -2.308  0.0228

# Detritus * family
m_det_f <- lm(sum ~ Habitat * Moisture * family, data = final_labeled_c_f_det_total)
anova(m_det_f)
# Response: sum
# Df     Sum Sq    Mean Sq F value  Pr(>F)    
# Habitat                   1 1.9002e+12 1.9002e+12  0.9418 0.33210    
# Moisture                  1 1.0724e+13 1.0724e+13  5.3150 0.02138 *  
#   family                  160 1.0313e+16 6.4455e+13 31.9446 < 2e-16 ***
#   Habitat:Moisture          1 9.0897e+10 9.0897e+10  0.0450 0.83196    
# Habitat:family          128 7.8056e+14 6.0981e+12  3.0223 < 2e-16 ***
#   Moisture:family         118 8.5056e+13 7.2082e+11  0.3572 1.00000    
# Habitat:Moisture:family  80 1.2742e+13 1.5928e+11  0.0789 1.00000    
# Residuals               848 1.7110e+15 2.0177e+12                 

em_det_hab_f <- emmeans(m_det_f, ~ Habitat | Moisture | family)
pairs(em_det_hab_f, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
# Moisture = Normal, family = Bacillaceae 1:
#   contrast                                    estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) -9236987 1159803 848  -7.964 <0.0001
# 
# Moisture = Drought, family = Bacillaceae 1:
#   contrast                                    estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) -8781555 1084896 848  -8.094 <0.0001
# 
# Moisture = Normal, family = Enterobacteriaceae:
#   contrast                                    estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  2335627 1159803 848   2.014  0.0443
# 
# Moisture = Normal, family = Streptomycetaceae:
#   contrast                                    estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) 12344942 1159803 848  10.644 <0.0001
# 
# Moisture = Drought, family = Streptomycetaceae:
#   contrast                                    estimate      SE  df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus) 12785560 1084896 848  11.785 <0.0001

# Detritus
m_det <- lm(sum ~ Habitat * Moisture, data = final_labeled_c_det_total)
anova(m_det)

# Response: sum
# Df     Sum Sq    Mean Sq F value Pr(>F)
# Habitat           1 1.9762e+12 1.9762e+12  0.0031 0.9570
# Moisture          1 5.9495e+14 5.9495e+14  0.9243 0.3615
# Habitat:Moisture  1 1.3516e+12 1.3516e+12  0.0021 0.9645
# Residuals         9 5.7932e+15 6.4369e+14     

em_det_hab <- emmeans(m_det, ~ Habitat | Moisture)
pairs(em_det_hab, adjust = "tukey")
# Moisture = Normal:
#   contrast                                    estimate       SE df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  2447136 20700000  9   0.118  0.9086
# 
# Moisture = Drought:
#   contrast                                    estimate       SE df t.ratio p.value
# (13C-Detritus) - (12C-Rhizo + 13C-Detritus)  1147310 19400000  9   0.059  0.9541

# Moisture = Drought:
#   contrast                              estimate      SE df t.ratio
# Detritus.13C - Rhizo.12C.Detritus.13C  1772152 3250000  9   0.546
# p.value
# 0.5983

em_det_moist <- emmeans(m_det, ~ Moisture | Habitat)
pairs(em_det_moist, adjust = "tukey")

# Rhizosphere * phylum
m_rhizo <- lm(sum ~ Habitat * Moisture * phylum, data = final_labeled_c_p_rhizo_total)
anova(m_rhizo)
# Response: sum
# Df     Sum Sq    Mean Sq F value  Pr(>F)    
# Habitat                   1 2.7770e+11 2.7770e+11  0.0092 0.92372    
# Moisture                  1 1.6347e+13 1.6347e+13  0.5414 0.46294    
# phylum                   16 5.3717e+15 3.3573e+14 11.1184 < 2e-16 ***
#   Habitat:Moisture          1 1.2894e+14 1.2894e+14  4.2702 0.04039 *  
#   Habitat:phylum           16 5.7343e+13 3.5839e+12  0.1187 0.99999    
# Moisture:phylum          16 2.8697e+14 1.7936e+13  0.5940 0.88524    
# Habitat:Moisture:phylum  14 5.7509e+14 4.1078e+13  1.3604 0.17837    
# Residuals               161 4.8616e+15 3.0196e+13      

em_rhizo_hab <- emmeans(m_rhizo, ~ Habitat | Moisture | phylum)
pairs(em_rhizo_hab, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
# Moisture = Normal, phylum = Actinobacteria:
#   contrast                                  estimate      SE  df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) -10871307 4196956 161  -2.590  0.0105
# 
# Moisture = Drought, phylum = Actinobacteria:
#   contrast                                  estimate      SE  df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus)  11749334 3885626 161   3.024  0.0029

# Rhizosphere * family
m_rhizo_f <- lm(sum ~ Habitat * Moisture * family, data = final_labeled_c_f_rhizo_total)
anova(m_rhizo_f)
# Response: sum
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# Habitat                    1 2.2924e+11 2.2924e+11  0.0936  0.759688    
# Moisture                   1 1.6571e+12 1.6571e+12  0.6768  0.410896    
# family                   165 2.1142e+15 1.2813e+13  5.2328 < 2.2e-16 ***
#   Habitat:Moisture           1 2.4117e+13 2.4117e+13  9.8489  0.001748 ** 
#   Habitat:family           136 4.4398e+13 3.2645e+11  0.1333  1.000000    
# Moisture:family          127 2.2196e+14 1.7478e+12  0.7138  0.991577    
# Habitat:Moisture:family   93 5.8724e+14 6.3144e+12  2.5787 6.364e-13 ***
#   Residuals               1026 2.5123e+15 2.4486e+12     

em_rhizo_hab_f <- emmeans(m_rhizo_f, ~ Habitat | Moisture | family)
pairs(em_rhizo_hab_f, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)

# Moisture = Normal, family = Bacillaceae 1:
#   contrast                                 estimate      SE   df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) -7711326 1195148 1026  -6.452 <0.0001
# 
# Moisture = Drought, family = Bacillaceae 1:
#   contrast                                 estimate      SE   df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) 11519991 1195148 1026   9.639 <0.0001
# 
# Moisture = Normal, family = Streptomycetaceae:
#   contrast                                 estimate      SE   df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus) -9094986 1195148 1026  -7.610 <0.0001
# 
# Moisture = Drought, family = Streptomycetaceae:
#   contrast                                 estimate      SE   df t.ratio p.value
# (13C-Rhizo) - (13C-Rhizo + 12C-Detritus)  8148997 1106492 1026   7.365 <0.0001

# Rhizosphere
m_rhizo <- lm(sum ~ Habitat * Moisture, data = final_labeled_c_rhizo_total)
 anova(m_rhizo)
 # Response: sum
 # Df     Sum Sq    Mean Sq F value Pr(>F)
 # Habitat           1 7.7318e+12 7.7318e+12  0.0058 0.9405
 # Moisture          1 1.2938e+14 1.2938e+14  0.0976 0.7606
 # Habitat:Moisture  1 1.8617e+15 1.8617e+15  1.4039 0.2610
 # Residuals        11 1.4587e+16 1.3261e+15        

em_rhizo_hab <- emmeans(m_rhizo, ~ Habitat | Moisture)
pairs(em_rhizo_hab, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(p.value < 0.05)
#none

