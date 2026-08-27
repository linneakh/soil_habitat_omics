# ***************************************************************
# 
# MetaboDirect
# Data Exploration step
# MetaboDirect version 1.0
# by Christian Ayala
# Licensed under the MIT license. See LICENSE.md file.
#
# Modified by Linnea Hernandez
# 8/25/26
# ***************************************************************

# Loading libraries ----

library(ggvenn)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(vegan)
library(UpSetR)
library(lmerTest)
library(lme4)
library(emmeans)
library(afex)
library(tidyverse)


# Defining variables ----

output_dir <- './output/FTICR/'
figure_dir <- './figures/FigSxx_FTICR/'


# Loading custom functions ----
source("./scripts/functions/FTICR_functions.R")
source("./scripts/functions/NMDS_plotting_functions.R")

class_colors <- get_palette(palette = 'Set3', k = 9)
names(class_colors) <- c(classification$Class, 'Other')


# Import data ----

## Loading tables
df <- read.csv("./data/FTICR/FTICR_Report_processed_noNorm_MolecFormulas.csv", header = TRUE)
metadata <- read.csv("./data/FTICR/metadata.csv", header = TRUE) %>% select(-notes)

# Reformat data files ----

## Change all metadata columns to factors to avoid problems at plotting

  
## Intensity data file
df_longer_orig <- df %>%
  unite(Mass_formula_class, c("Mass", "MolecularFormula", "Class")) %>%
  pivot_longer(metadata$SampleID, names_to = "SampleID", values_to = "Intensity") %>%
  filter(Intensity > 0) %>% 
  left_join(metadata, by = 'SampleID') %>%
  mutate(Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks")))



#create new column called "presence" to indicate if compound found in both or are unique

present_compounds <- df_longer_orig %>%
  # extract replicate ID from sample name
  mutate(rep_id = sub("^X(\\d+)[A-Za-z0-9]*_.*", "\\1", SampleID)) %>%
  group_by(Mass_formula_class, Timepoint, Treatment, Zone) %>%
  mutate(Replicate = dense_rank(rep_id)) %>%
  ungroup() %>%
  distinct(Mass_formula_class, Replicate, Timepoint, Treatment, Zone) %>%
  # count how many replicates detected per compound x treatment combo
  group_by(Mass_formula_class, Treatment, Zone, Timepoint) %>%
  summarise(n_reps = n(), .groups = "drop") %>%
  # create one label per treatment combination
  mutate(combo = paste(Treatment, Timepoint, Zone, sep = "_")) %>%
  mutate(present = if_else(n_reps >= 2, "present", "not_present"))

# wide table: one row per compound, one column per treatment combo
presence_wide <- present_compounds %>%
  select(Mass_formula_class, combo, present) %>%
  pivot_wider(names_from = combo, values_from = present, values_fill = "not_present")

# long table for plotting: which combos is each compound present in
presence_long <- present_compounds %>%
  select(Mass_formula_class, combo, present)

df_group_moisture <- df_longer_orig %>%
  select(Mass_formula_class) %>% 
  distinct() %>% 
  merge(presence_long, by = "Mass_formula_class") %>%
  select(Mass_formula_class, combo, present) %>%
  filter(present == "present") %>%
  separate(combo, c("Treatment", "Timepoint", "Zone"), sep = "_") 

#filter df to keep only masses that are in at least 2 replicates from a treatment combination
df_longer_orig_f <- df_longer_orig %>%
  merge(df_group_moisture, by = c("Mass_formula_class", "Treatment", "Timepoint", "Zone"))
  

#class table
df_longer_class <- df_longer_orig_f %>%
  separate(Mass_formula_class, c("Mass", "Formula", "Class"), sep = "_") %>%
  group_by(SampleID, Treatment, Timepoint, Zone, Class) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(Class, nesting(SampleID, Treatment, Timepoint, Zone), fill = list(count = 0)) %>%
  group_by(SampleID) %>%
  mutate(rel_abundance = count / sum(count) * 100) %>%
  ungroup()

#####----- Plotting NOSC and size----####
NOSC_summary <- df_longer_orig %>%
  merge(df_group_moisture, by = c("Mass_formula_class", "Treatment", "Timepoint", "Zone")) %>%
  group_by(Zone, Treatment, Timepoint, SampleID) %>%
  summarise(mean_sample = mean(NOSC)) %>%
  group_by(Zone, Treatment, Timepoint) %>%
  summarise(
    mean = mean(mean_sample, na.rm = TRUE),
    se   = sd(mean_sample, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought")))

nosc_line_plot <- ggplot(NOSC_summary, 
                          aes(x = Timepoint, y =  mean, 
                              color = Zone, group = Zone)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = 0.1, linewidth = 1) +
  scale_color_manual(values = col_list_zone) +
  labs(y = "Nominal Oxidation State of Carbon (NOSC)", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Treatment) +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

nosc_line_plot
ggsave(paste(figure_dir, "/FigSx_NOSC.png", sep = ""), dpi = 300, height = 4, width = 8, unit = "in")
ggsave(paste(figure_dir, "/FigSx_NOSC.pdf", sep = ""), dpi = 300, height = 4, width = 8, unit = "in")


#plot size
size_summary <- df_longer_orig %>%
  separate(Mass_formula_class, c("Mass", "Formula", "Class"), sep = "_") %>%
  mutate(Mass = as.numeric(Mass)) %>%
  group_by(Zone, Treatment, Timepoint, SampleID) %>%
  summarise(mean_sample = mean(Mass)) %>%
  group_by(Zone, Treatment, Timepoint) %>%
  summarise(
    mean = mean(mean_sample, na.rm = TRUE),
    se   = sd(mean_sample, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought")))

size_line_plot <- ggplot(size_summary, 
                         aes(x = Timepoint, y =  mean, 
                             color = Zone, group = Zone)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = 0.1, linewidth = 1) +
  scale_color_manual(values = col_list_zone) +
  labs(y = "Mass (m/z)", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Treatment) +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
size_line_plot
ggsave(paste(figure_dir, "/FigSx_mass.png", sep = ""), dpi = 300, height = 4, width = 8, unit = "in")
ggsave(paste(figure_dir, "/FigSx_mass.pdf", sep = ""), dpi = 300, height = 4, width = 8, unit = "in")



####----------ANOVA on NOSC and size--------####
##NOSC normal
df_longer_NOSC_Norm <- df_longer_orig_f %>%
  filter(Treatment == "Normal") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(NOSC = mean(NOSC)) 

#lmer
m <- lm(NOSC ~ Zone * Timepoint,  data = df_longer_NOSC_Norm)
anova(m)

# Response: NOSC
# Df   Sum Sq   Mean Sq F value  Pr(>F)  
# Zone            3 0.036628 0.0122093  3.9480 0.01539 *
#   Timepoint       2 0.021498 0.0107488  3.4757 0.04137 *
#   Zone:Timepoint  6 0.060503 0.0100839  3.2607 0.01127 *
#   Residuals      37 0.114424 0.0030925      

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")

# Timepoint = 4weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus         0.01752 0.0393 37   0.445  0.9701
# Bulk - RhizoDet         0.06900 0.0393 37   1.755  0.3111
# Bulk - Rhizosphere      0.01249 0.0393 37   0.318  0.9887
# Detritus - RhizoDet     0.05148 0.0393 37   1.309  0.5631
# Detritus - Rhizosphere -0.00502 0.0393 37  -0.128  0.9992
# RhizoDet - Rhizosphere -0.05650 0.0393 37  -1.437  0.4851
# 
# Timepoint = 8weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus         0.03308 0.0393 37   0.841  0.8344
# Bulk - RhizoDet        -0.03708 0.0373 37  -0.994  0.7537
# Bulk - Rhizosphere      0.06436 0.0393 37   1.637  0.3712
# Detritus - RhizoDet    -0.07016 0.0373 37  -1.881  0.2537
# Detritus - Rhizosphere  0.03128 0.0393 37   0.795  0.8560
# RhizoDet - Rhizosphere  0.10144 0.0373 37   2.719  0.0467
# 
# Timepoint = 12weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus         0.05218 0.0393 37   1.327  0.5521
# Bulk - RhizoDet         0.13384 0.0393 37   3.404  0.0084
# Bulk - Rhizosphere      0.15249 0.0393 37   3.878  0.0023
# Detritus - RhizoDet     0.08166 0.0393 37   2.077  0.1796
# Detritus - Rhizosphere  0.10031 0.0393 37   2.551  0.0684
# RhizoDet - Rhizosphere  0.01865 0.0393 37   0.474  0.9643

##NOSC drought
df_longer_nosc_drought <- df_longer_orig_f %>%
  filter(Treatment == "Drought") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(NOSC = mean(NOSC)) 

#lmer
m <- lm(NOSC ~ Zone * Timepoint,  data = df_longer_nosc_drought)
anova(m)

# Response: NOSC
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Zone            3 0.060125 0.0200417  9.9329 6.542e-05 ***
#   Timepoint       2 0.006529 0.0032643  1.6178  0.212421    
# Zone:Timepoint  6 0.045405 0.0075674  3.7505  0.005317 ** 
#   Residuals      36 0.072638 0.0020177           

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")
# Timepoint = 4weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus          0.0488 0.0318 36   1.535  0.4279
# Bulk - RhizoDet          0.0911 0.0318 36   2.868  0.0332
# Bulk - Rhizosphere      -0.0580 0.0318 36  -1.826  0.2781
# Detritus - RhizoDet      0.0423 0.0318 36   1.333  0.5486
# Detritus - Rhizosphere  -0.1068 0.0318 36  -3.361  0.0096
# RhizoDet - Rhizosphere  -0.1491 0.0318 36  -4.694  0.0002
# 
# Timepoint = 8weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus          0.0427 0.0318 36   1.344  0.5415
# Bulk - RhizoDet          0.0178 0.0318 36   0.559  0.9434
# Bulk - Rhizosphere      -0.0629 0.0318 36  -1.981  0.2139
# Detritus - RhizoDet     -0.0249 0.0318 36  -0.785  0.8606
# Detritus - Rhizosphere  -0.1056 0.0318 36  -3.326  0.0105
# RhizoDet - Rhizosphere  -0.0807 0.0318 36  -2.540  0.0705
# 
# Timepoint = 12weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus          0.1012 0.0318 36   3.186  0.0151
# Bulk - RhizoDet         -0.0142 0.0318 36  -0.448  0.9697
# Bulk - Rhizosphere       0.0286 0.0318 36   0.900  0.8047
# Detritus - RhizoDet     -0.1154 0.0318 36  -3.634  0.0046
# Detritus - Rhizosphere  -0.0726 0.0318 36  -2.286  0.1204
# RhizoDet - Rhizosphere   0.0428 0.0318 36   1.348  0.5394
# #

##size normal
df_longer_mass_normal <- df_longer_orig_f %>%
  separate(Mass_formula_class, c("Mass", "Formula", "Class"), sep = "_") %>%
  mutate(Mass = as.numeric(Mass)) %>%
  filter(Treatment == "Normal") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(Mass = mean(Mass)) 

#lmer
m <- lm(Mass ~ Zone * Timepoint,  data = df_longer_mass_normal)
anova(m)

# Response: Mass
# Df  Sum Sq Mean Sq F value   Pr(>F)    
# Zone            3 1657.34  552.45  8.3901 0.000221 ***
#   Timepoint       2   67.04   33.52  0.5091 0.605185    
# Zone:Timepoint  6 1510.49  251.75  3.8234 0.004590 ** 
#   Residuals      37 2436.26   65.84       

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")
# Timepoint = 4weeks:
#   contrast               estimate   SE df t.ratio p.value
# Bulk - Detritus         -15.761 5.74 37  -2.747  0.0438
# Bulk - RhizoDet         -24.908 5.74 37  -4.341  0.0006
# Bulk - Rhizosphere      -10.418 5.74 37  -1.816  0.2824
# Detritus - RhizoDet      -9.147 5.74 37  -1.594  0.3943
# Detritus - Rhizosphere    5.343 5.74 37   0.931  0.7883
# RhizoDet - Rhizosphere   14.490 5.74 37   2.525  0.0724
# 
# Timepoint = 8weeks:
#   contrast               estimate   SE df t.ratio p.value
# Bulk - Detritus         -19.471 5.74 37  -3.393  0.0086
# Bulk - RhizoDet         -20.232 5.44 37  -3.717  0.0036
# Bulk - Rhizosphere      -20.639 5.74 37  -3.597  0.0050
# Detritus - RhizoDet      -0.762 5.44 37  -0.140  0.9990
# Detritus - Rhizosphere   -1.168 5.74 37  -0.204  0.9970
# RhizoDet - Rhizosphere   -0.407 5.44 37  -0.075  0.9998
# 
# Timepoint = 12weeks:
#   contrast               estimate   SE df t.ratio p.value
# Bulk - Detritus         -10.833 5.74 37  -1.888  0.2506
# Bulk - RhizoDet           6.214 5.74 37   1.083  0.7019
# Bulk - Rhizosphere        1.938 5.74 37   0.338  0.9865
# Detritus - RhizoDet      17.048 5.74 37   2.971  0.0256
# Detritus - Rhizosphere   12.771 5.74 37   2.226  0.1351
# RhizoDet - Rhizosphere   -4.276 5.74 37  -0.745  0.8781

##size drought
df_longer_mass_drought <- df_longer_orig_f %>%
  separate(Mass_formula_class, c("Mass", "Formula", "Class"), sep = "_") %>%
  mutate(Mass = as.numeric(Mass)) %>%
  filter(Treatment == "Drought") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(Mass = mean(Mass)) 

#lmer
m <- lm(Mass ~ Zone * Timepoint,  data = df_longer_mass_drought)
anova(m)

# Response: Mass
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Zone            3 3529.8 1176.61 55.8920 1.264e-13 ***
#   Timepoint       2  219.0  109.51  5.2019   0.01036 *  
#   Zone:Timepoint  6 1341.2  223.53 10.6182 8.947e-07 ***
#   Residuals      36  757.9   21.05   

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")

# Timepoint = 4weeks:
#   contrast               estimate   SE df t.ratio p.value
# Bulk - Detritus          -18.46 3.24 36  -5.691 <0.0001
# Bulk - RhizoDet          -21.94 3.24 36  -6.761 <0.0001
# Bulk - Rhizosphere       -28.68 3.24 36  -8.841 <0.0001
# Detritus - RhizoDet       -3.47 3.24 36  -1.071  0.7093
# Detritus - Rhizosphere   -10.22 3.24 36  -3.150  0.0165
# RhizoDet - Rhizosphere    -6.75 3.24 36  -2.080  0.1791
# 
# Timepoint = 8weeks:
#   contrast               estimate   SE df t.ratio p.value
# Bulk - Detritus          -18.26 3.24 36  -5.628 <0.0001
# Bulk - RhizoDet          -21.15 3.24 36  -6.518 <0.0001
# Bulk - Rhizosphere       -14.87 3.24 36  -4.582  0.0003
# Detritus - RhizoDet       -2.89 3.24 36  -0.889  0.8104
# Detritus - Rhizosphere     3.40 3.24 36   1.047  0.7235
# RhizoDet - Rhizosphere     6.28 3.24 36   1.936  0.2315
# 
# Timepoint = 12weeks:
#   contrast               estimate   SE df t.ratio p.value
# Bulk - Detritus          -28.33 3.24 36  -8.732 <0.0001
# Bulk - RhizoDet          -17.01 3.24 36  -5.243 <0.0001
# Bulk - Rhizosphere        -3.91 3.24 36  -1.206  0.6272
# Detritus - RhizoDet       11.32 3.24 36   3.489  0.0068
# Detritus - Rhizosphere    24.42 3.24 36   7.526 <0.0001
# RhizoDet - Rhizosphere    13.10 3.24 36   4.037  0.0015


####----------Plot classes--------####
#normal
class_line_time_normal <- df_longer_class %>%
  filter(Treatment == "Normal") %>%
  filter(Class != "Other") %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought"))) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks"))) %>%
  mutate(Class = factor(Class, levels = c("Amino sugar", "Lipid", "Protein",
                                          "Carbohydrate", "Cond. HC", "Unsat. HC",
                                          "Tannin", "Lignin"))) %>%
  group_by(Timepoint, Zone, Treatment, Class) %>%
  summarise(
    mean_ra = mean(rel_abundance, na.rm = TRUE),
    se_ra   = sd(rel_abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

class_line_plot_normal <- ggplot(class_line_time_normal, 
                          aes(x = Timepoint, y = mean_ra, 
                              color = Zone, group = Zone)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_ra - se_ra, 
                    ymax = mean_ra + se_ra),
                width = 0.2, linewidth = 0.4) +
  scale_color_manual(values = col_list_zone) +
  labs(y = "Relative abundance", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Class, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

class_line_plot_normal

filename <- file.path(figure_dir, 'FigSx-class-line-normal.png')
ggsave(filename, class_line_plot_normal, dpi = 300, width = 12, height = 5.5, units = "in")
filename <- file.path(figure_dir, 'FigSx-class-line-normal.pdf')
ggsave(filename, class_line_plot_normal, dpi = 300, width = 12, height = 5.5, units = "in")

#drought
class_line_time_drought <- df_longer_class %>%
  filter(Treatment == "Drought") %>%
  filter(Class != "Other") %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought"))) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks"))) %>%
  mutate(Class = factor(Class, levels = c("Amino sugar", "Lipid", "Protein",
                                          "Carbohydrate", "Cond. HC", "Unsat. HC",
                                          "Tannin", "Lignin"))) %>%
  group_by(Timepoint, Zone, Treatment, Class) %>%
  summarise(
    mean_ra = mean(rel_abundance, na.rm = TRUE),
    se_ra   = sd(rel_abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

class_line_plot_drought <- ggplot(class_line_time_drought, 
                                 aes(x = Timepoint, y = mean_ra, 
                                     color = Zone, group = Zone)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_ra - se_ra, 
                    ymax = mean_ra + se_ra),
                width = 0.2, linewidth = 0.4) +
  scale_color_manual(values = col_list_zone) +
  labs(y = "Relative abundance", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Class, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

class_line_plot_drought

filename <- file.path(figure_dir, 'FigSx-class-line-drought.png')
ggsave(filename, class_line_plot_drought, dpi = 300, width = 12, height = 5.5, units = "in")
filename <- file.path(figure_dir, 'FigSx-class-line-drought.pdf')
ggsave(filename, class_line_plot_drought, dpi = 300, width = 12, height = 5.5, units = "in")

#run ANOVA
#normal 
df_longer_class_normal <- df_longer_class %>%
  filter(Treatment == "Normal")

anova_class_normal <- run_lm_tukey (data = df_longer_class_normal,
              classes = unique(df_longer_class_normal$Class),
              treatment_filter = NULL,
              emmeans_formula = ~ Zone | Timepoint)

write.csv(anova_class_normal$anova, paste(output_dir, "FigSx-stats-anova-normal.csv"))
write.csv(anova_class_normal$tukey, paste(output_dir, "FigSx-stats-tukey-normal.csv"))

#drought 
df_longer_class_drought <- df_longer_class %>%
  filter(Treatment == "Drought")

anova_class_drought <- run_lm_tukey (data = df_longer_class_drought,
                                    classes = unique(data$Class),
                                    treatment_filter = NULL,
                                    emmeans_formula = ~ Zone | Timepoint)

write.csv(anova_class_drought$anova, paste(output_dir, "FigSx-stats-anova-drought.csv"))
write.csv(anova_class_drought$tukey, paste(output_dir, "FigSx-stats-tukey-drought.csv"))


####------ NMDS----------####
# parameters for plots
h = 6
w = 10
res = 300
size = 16


#widen fticr data frame
df_orig <- df_longer_orig_f %>%
  select(-c(C, H, O, N, S, P, NeutralMass, Error_ppm, El_comp,
            OC, HC, NOSC, GFE, DBE, DBE_O, AI, AI_mod, DBE_AI, present, Treatment,
            Timepoint, Zone, )) %>%
  spread(key = "SampleID", value = "Intensity") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  column_to_rownames(var = "Mass_formula_class")

#create a presence/absence table for jaccard beta diversity
df_orig_pres_ab.0 <- df_orig %>%
  mutate(across(where(is.numeric), ~ ifelse(.x > 0, 1, 0)))

#transpose table
df_orig_pres_ab <- df_orig_pres_ab.0 %>%
  t() %>%
  as.data.frame()

##NMDS
# Calculate nmds

set.seed(123)
nmds = metaMDS(df_orig_pres_ab, distance = "jaccard")
plot(nmds)

scores <- scores(nmds)
scores.samples <- scores$sites

#extract NMDS scores for x and y coordinates
data.scores = as.data.frame(scores.samples)
#data.scores = as.data.frame(scores)


data.scores$SampleID <- rownames(data.scores)
data.scores <- data.scores  %>% 
   merge(metadata, by = "SampleID") %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk")))
  

#determine environmetnal drivers (how do environmetnal variables drive patterns)

ef <- envfit(nmds ~ Zone + Treatment + Timepoint, data = data.scores,
             permutations = 999)
plot(nmds, type = "n")
points(nmds, pch = 21, bg = "lightgreen")
plot(ef, col = "blue", lwd = 2)  # draws arrows for each env var

# Goodness of fit:
#   r2 Pr(>r)    
# Zone      0.7312  0.001 ***
#   Treatment 0.0153  0.235    
# Timepoint 0.0151  0.563  

### nmds plot (timepoint and moisture)
make_nmds_plot(data.scores, Group1 =Zone, Group2 =Treatment, col_list_zone)

filename <- paste0(output_dir, "FigSx_zone_treatment.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)
filename <- paste0(output_dir, "FigSx_zone_treatment.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)

### nmds plot (timepoint and moisture)
make_nmds_plot(data.scores, Group1 =Zone, Group2 =Timepoint, col_list_zone)


filename <- paste0(output_dir, "zone_timepoint.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res)


##statistical tests
#df.df <- as.data.frame(df)
df.df <- as.data.frame(df_orig_pres_ab)

df.df <- df.df[sapply(df.df, is.numeric)]


##all time points##
dist <- vegdist(df.df, method = "jaccard")  # binary matrix

adonis2(dist ~ Zone + Treatment + Timepoint, data.scores, permutations = 999, by = "margin")
# Df SumOfSqs      R2       F Pr(>F)    
# Zone       3   3.8700 0.48150 31.3028  0.001 ***
#   Treatment  1   0.2195 0.02731  5.3257  0.001 ***
#   Timepoint  2   0.2444 0.03041  2.9651  0.006 ** 
#   Residual  90   3.7089 0.46147                   
# Total     96   8.0372 1.00000 













