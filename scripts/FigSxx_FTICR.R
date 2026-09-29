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

# class_colors <- get_palette(palette = 'Set3', k = 9)
# names(class_colors) <- c(classification$Class, 'Other')


# Import data ----

## Loading tables
#compounds with formulas
df <- read.csv("./data/FTICR/Report_processed_noNorm_MolecFormulas.csv", header = TRUE)

#all compounds
df.all <- read.csv("./data/FTICR/Report_processed_noNorm.csv", header = TRUE)

metadata <- read.csv("./data/FTICR/metadata.csv", header = TRUE) %>% select(-notes)

trans <- read.csv("./data/FTICR/transformations_summary_all.csv", header = TRUE) %>%
  merge(metadata, by = "SampleID") %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought"))) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks")))

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
#NOSC line
count_summary <- df_longer_orig_f %>%
  #merge(df_group_moisture, by = c("Mass_formula_class", "Treatment", "Timepoint", "Zone")) %>%
  group_by(Zone, Treatment, Timepoint, SampleID) %>%
  summarise(count = n()) %>%
  group_by(Zone, Treatment, Timepoint) %>%
  summarise(
    mean = mean(count, na.rm = TRUE),
    se   = sd(count, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought")))

count_line_plot <- ggplot(count_summary, 
                         aes(x = Timepoint, y =  mean, 
                             color = Zone, group = Zone)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = 0.1, linewidth = 1) +
  scale_color_manual(values = col_list_zone) +
  labs(y = "Number of compounds", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Treatment) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(text = element_text(size = 14))

# ,
# panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(),
# strip.background = element_blank())

count_line_plot
ggsave(paste(figure_dir, "/FigSx_count.png", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")
ggsave(paste(figure_dir, "/FigSx_count.pdf", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")

#NOSC line
NOSC_summary <- df_longer_orig_f %>%
  #merge(df_group_moisture, by = c("Mass_formula_class", "Treatment", "Timepoint", "Zone")) %>%
  group_by(Zone, Treatment, Timepoint, SampleID) %>%
  summarise(med_sample = median(NOSC)) %>%
  group_by(Zone, Treatment, Timepoint) %>%
  summarise(
    mean = mean(med_sample, na.rm = TRUE),
    se   = sd(med_sample, na.rm = TRUE) / sqrt(n()),
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
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(text = element_text(size = 14))
        
        # ,
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        # strip.background = element_blank())

nosc_line_plot
ggsave(paste(figure_dir, "/FigSx_NOSC.png", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")
ggsave(paste(figure_dir, "/FigSx_NOSC.pdf", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")

#NOSC violin

nosc_violin_plot <- df_longer_orig_f %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought"))) %>%
  ggplot(aes(x = Timepoint, y = NOSC, fill = Zone)) +
  geom_violin(alpha = 0.2, 
              position = position_dodge(width = 0.9), 
              show.legend = FALSE,
              linewidth = 0.2) +
  geom_boxplot(width = 0.15, 
               position = position_dodge(width = 0.9),
               outlier.shape = NA,
               outliers = FALSE,
               linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = col_list_zone) +
  facet_wrap(~Treatment, ncol = 1) +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())


ggsave(paste(figure_dir, "/FigSx_NOSC_violin.png", sep = ""), dpi = 300, height = 5, width = 7, unit = "in")
ggsave(paste(figure_dir, "/FigSx_NOSC_violin.pdf", sep = ""), dpi = 300, height = 5, width = 7, unit = "in")

#DBE line
DBE_summary <- df_longer_orig_f %>%
  drop_na(DBE) %>%
  #merge(df_group_moisture, by = c("Mass_formula_class", "Treatment", "Timepoint", "Zone")) %>%
  group_by(Zone, Treatment, Timepoint, SampleID) %>%
  summarise(med_sample = median(DBE)) %>%
  group_by(Zone, Treatment, Timepoint) %>%
  summarise(
    mean = mean(med_sample, na.rm = TRUE),
    se   = sd(med_sample, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought")))

DBE_line_plot <- ggplot(DBE_summary, 
                         aes(x = Timepoint, y =  mean, 
                             color = Zone, group = Zone)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = 0.1, linewidth = 1) +
  scale_color_manual(values = col_list_zone) +
  labs(y = "Double Bond Equivelent (DBE)", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Treatment) +
  theme_bw() +
  theme(text = element_text(size = 14))

# ,
# panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(),
# strip.background = element_blank())

DBE_line_plot
ggsave(paste(figure_dir, "/FigSx_DBE.png", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")
ggsave(paste(figure_dir, "/FigSx_DBE.pdf", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")

#dbe violin

dbe_violin_plot <- df_longer_orig_f %>%
  drop_na(DBE) %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought"))) %>%
  ggplot(aes(x = Timepoint, y = DBE, fill = Zone)) +
  geom_violin(alpha = 0.2, 
              position = position_dodge(width = 0.9), 
              show.legend = FALSE,
              linewidth = 0.2) +
  geom_boxplot(width = 0.15, 
               position = position_dodge(width = 0.9),
               outlier.shape = NA,
               outliers = FALSE,
               linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = col_list_zone) +
  facet_wrap(~Treatment, ncol = 1) +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
dbe_violin_plot
ggsave(paste(figure_dir, "/FigSx_DBE_violin.png", sep = ""), dpi = 300, height = 5, width = 7, unit = "in")
ggsave(paste(figure_dir, "/FigSx_DBE_violin.pdf", sep = ""), dpi = 300, height = 5, width = 7, unit = "in")

#AI line
AI_summary <- df_longer_orig_f %>%
  drop_na(AI) %>%
  #merge(df_group_moisture, by = c("Mass_formula_class", "Treatment", "Timepoint", "Zone")) %>%
  group_by(Zone, Treatment, Timepoint, SampleID) %>%
  summarise(med_sample = median(AI)) %>%
  group_by(Zone, Treatment, Timepoint) %>%
  summarise(
    mean = mean(med_sample, na.rm = TRUE),
    se   = sd(med_sample, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought")))

AI_line_plot <- ggplot(AI_summary, 
                        aes(x = Timepoint, y =  mean, 
                            color = Zone, group = Zone)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = 0.1, linewidth = 1) +
  scale_color_manual(values = col_list_zone) +
  labs(y = "Aromaticity Index (AI)", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Treatment) +
  theme_bw() +
  theme(text = element_text(size = 14))

# ,
# panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(),
# strip.background = element_blank())

AI_line_plot
ggsave(paste(figure_dir, "/FigSx_AI.png", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")
ggsave(paste(figure_dir, "/FigSx_AI.pdf", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")

#ai violin

ai_violin_plot <- df_longer_orig_f %>%
  drop_na(AI) %>%
  mutate(Zone = factor(Zone, levels = c("Rhizosphere", "RhizoDet", "Detritus", "Bulk"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Normal", "Drought"))) %>%
  ggplot(aes(x = Timepoint, y = AI, fill = Zone)) +
  geom_violin(alpha = 0.2, 
              position = position_dodge(width = 0.9), 
              show.legend = FALSE,
              linewidth = 0.2) +
  geom_boxplot(width = 0.15, 
               position = position_dodge(width = 0.9),
               outlier.shape = NA,
               outliers = FALSE,
               linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = col_list_zone) +
  facet_wrap(~Treatment, ncol = 1) +
  ylim(-2, 2.5) +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ai_violin_plot
ggsave(paste(figure_dir, "/FigSx_AI_violin.png", sep = ""), dpi = 300, height = 5, width = 7, unit = "in")
ggsave(paste(figure_dir, "/FigSx_AI_violin.pdf", sep = ""), dpi = 300, height = 5, width = 7, unit = "in")



####----------ANOVA on NOSC and DBE--------####
##NOSC normal
df_longer_NOSC_Norm <- df_longer_orig_f %>%
  filter(Treatment == "Normal") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(NOSC = median(NOSC)) 

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
  summarise(NOSC = median(NOSC)) 

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

##DBE normal
df_longer_dbe_normal <- df_longer_orig_f %>%
  drop_na(DBE) %>%
 # mutate(DBE = as.numeric(DBE)) %>%
  filter(Treatment == "Normal") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(DBE = median(DBE)) 

#lmer
m <- lm(DBE ~ Zone * Timepoint,  data = df_longer_dbe_normal)
anova(m)

# Response: DBE
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Zone            3  3.7433 1.24778  3.4198 0.0270680 *  
#   Timepoint       2  3.6950 1.84752  5.0636 0.0113819 *  
#   Zone:Timepoint  6 12.6127 2.10211  5.7613 0.0002604 ***
#   Residuals      37 13.5000 0.36486              

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")
# Timepoint = 4weeks:
#   contrast               estimate    SE df t.ratio p.value
# Bulk - Detritus           -1.00 0.427 37  -2.341  0.1070
# Bulk - RhizoDet           -0.25 0.427 37  -0.585  0.9359
# Bulk - Rhizosphere        -0.75 0.427 37  -1.756  0.3104
# Detritus - RhizoDet        0.75 0.427 37   1.756  0.3104
# Detritus - Rhizosphere     0.25 0.427 37   0.585  0.9359
# RhizoDet - Rhizosphere    -0.50 0.427 37  -1.171  0.6488
# 
# Timepoint = 8weeks:
#   contrast               estimate    SE df t.ratio p.value
# Bulk - Detritus           -0.25 0.427 37  -0.585  0.9359
# Bulk - RhizoDet           -0.75 0.405 37  -1.851  0.2666
# Bulk - Rhizosphere         0.50 0.427 37   1.171  0.6488
# Detritus - RhizoDet       -0.50 0.405 37  -1.234  0.6097
# Detritus - Rhizosphere     0.75 0.427 37   1.756  0.3104
# RhizoDet - Rhizosphere     1.25 0.405 37   3.085  0.0192
# 
# Timepoint = 12weeks:
#   contrast               estimate    SE df t.ratio p.value
# Bulk - Detritus            0.75 0.427 37   1.756  0.3104
# Bulk - RhizoDet            1.75 0.427 37   4.097  0.0012
# Bulk - Rhizosphere         2.00 0.427 37   4.683  0.0002
# Detritus - RhizoDet        1.00 0.427 37   2.341  0.1070
# Detritus - Rhizosphere     1.25 0.427 37   2.927  0.0285
# RhizoDet - Rhizosphere     0.25 0.427 37   0.585  0.9359

##DBE drought
df_longer_dbe_normal <- df_longer_orig_f %>%
  drop_na(DBE) %>%
  # mutate(DBE = as.numeric(DBE)) %>%
  filter(Treatment == "Drought") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(DBE = median(DBE)) 

#lmer
m <- lm(DBE ~ Zone * Timepoint,  data = df_longer_dbe_normal)
anova(m)

# Response: DBE
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Zone            3 2.7292 0.90972  3.5405 0.02406 *
#   Timepoint       2 1.5417 0.77083  3.0000 0.06237 .
# Zone:Timepoint  6 3.4583 0.57639  2.2432 0.06103 .
# Residuals      36 9.2500 0.25694    

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone)
pairs(em, adjust="tukey")
# 
# contrast               estimate    SE df t.ratio p.value
# Bulk - Detritus          0.1667 0.207 36   0.805  0.8515
# Bulk - RhizoDet         -0.3333 0.207 36  -1.611  0.3855
# Bulk - Rhizosphere      -0.4167 0.207 36  -2.013  0.2020
# Detritus - RhizoDet     -0.5000 0.207 36  -2.416  0.0920
# Detritus - Rhizosphere  -0.5833 0.207 36  -2.819  0.0373
# RhizoDet - Rhizosphere  -0.0833 0.207 36  -0.403  0.9776

##AI normal
df_longer_ai_normal <- df_longer_orig_f %>%
  drop_na(DBE) %>%
  # mutate(DBE = as.numeric(DBE)) %>%
  filter(Treatment == "Normal") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(AI = median(AI)) 

#lmer
m <- lm(AI ~ Zone * Timepoint,  data = df_longer_ai_normal)
anova(m)

# Response: AI
# Df   Sum Sq   Mean Sq F value   Pr(>F)   
# Zone            3 0.029892 0.0099641  6.6077 0.002061 **
#   Timepoint       2 0.008997 0.0044983  2.9830 0.069659 . 
# Zone:Timepoint  5 0.035909 0.0071819  4.7627 0.003653 **
#   Residuals      24 0.036191 0.0015080             

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")
# Timepoint = 4weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus        -0.06845 0.0336 24  -2.035  0.2033
# Bulk - RhizoDet        -0.01190 0.0336 24  -0.354  0.9844
# Bulk - Rhizosphere     -0.03750 0.0336 24  -1.115  0.6840
# Detritus - RhizoDet     0.05655 0.0388 24   1.456  0.4784
# Detritus - Rhizosphere  0.03095 0.0388 24   0.797  0.8551
# RhizoDet - Rhizosphere -0.02560 0.0388 24  -0.659  0.9113
# 
# Timepoint = 8weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus        -0.06892 0.0336 24  -2.050  0.1984
# Bulk - RhizoDet        -0.07566 0.0260 24  -2.904  0.0365
# Bulk - Rhizosphere      0.04079 0.0275 24   1.485  0.4614
# Detritus - RhizoDet    -0.00674 0.0325 24  -0.207  0.9968
# Detritus - Rhizosphere  0.10971 0.0336 24   3.262  0.0163
# RhizoDet - Rhizosphere  0.11645 0.0260 24   4.470  0.0009
# 
# Timepoint = 12weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus          nonEst     NA NA      NA      NA
# Bulk - RhizoDet         0.09212 0.0336 24   2.739  0.0296
# Bulk - Rhizosphere      0.10390 0.0336 24   3.090  0.0134
# Detritus - RhizoDet      nonEst     NA NA      NA      NA
# Detritus - Rhizosphere   nonEst     NA NA      NA      NA
# RhizoDet - Rhizosphere  0.01178 0.0275 24   0.429  0.9040

##DBE drought
df_longer_AI_drought <- df_longer_orig_f %>%
  drop_na(AI) %>%
  # mutate(DBE = as.numeric(DBE)) %>%
  filter(Treatment == "Drought") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(AI = median(AI)) 

#lmer
m <- lm(AI ~ Zone * Timepoint,  data = df_longer_AI_drought)
anova(m)

# Response: AI
# Df   Sum Sq   Mean Sq F value  Pr(>F)  
# Zone            3 0.014004 0.0046681  4.1162 0.01310 *
#   Timepoint       2 0.004457 0.0022286  1.9651 0.15488  
# Zone:Timepoint  6 0.016539 0.0027565  2.4306 0.04477 *
#   Residuals      36 0.040827 0.0011341      

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")
# Timepoint = 4weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus        -0.04106 0.0238 36  -1.724  0.3263
# Bulk - RhizoDet        -0.00449 0.0238 36  -0.188  0.9976
# Bulk - Rhizosphere     -0.02017 0.0238 36  -0.847  0.8315
# Detritus - RhizoDet     0.03657 0.0238 36   1.536  0.4274
# Detritus - Rhizosphere  0.02089 0.0238 36   0.877  0.8166
# RhizoDet - Rhizosphere -0.01569 0.0238 36  -0.659  0.9118
# 
# Timepoint = 8weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus        -0.03016 0.0238 36  -1.267  0.5895
# Bulk - RhizoDet        -0.06619 0.0238 36  -2.780  0.0409
# Bulk - Rhizosphere     -0.02707 0.0238 36  -1.137  0.6696
# Detritus - RhizoDet    -0.03602 0.0238 36  -1.513  0.4406
# Detritus - Rhizosphere  0.00310 0.0238 36   0.130  0.9992
# RhizoDet - Rhizosphere  0.03912 0.0238 36   1.643  0.3683
# 
# Timepoint = 12weeks:
#   contrast               estimate     SE df t.ratio p.value
# Bulk - Detritus         0.01325 0.0238 36   0.556  0.9442
# Bulk - RhizoDet        -0.07256 0.0238 36  -3.047  0.0214
# Bulk - Rhizosphere     -0.00716 0.0238 36  -0.301  0.9904
# Detritus - RhizoDet    -0.08581 0.0238 36  -3.604  0.0050
# Detritus - Rhizosphere -0.02041 0.0238 36  -0.857  0.8267
# RhizoDet - Rhizosphere  0.06540 0.0238 36   2.747  0.0442


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
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(y = "Relative abundance", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Class, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

class_line_plot_normal

filename <- file.path(figure_dir, 'FigSx-class-line-normal.png')
ggsave(filename, class_line_plot_normal, dpi = 300, width = 13, height = 5.5, units = "in")
filename <- file.path(figure_dir, 'FigSx-class-line-normal.pdf')
ggsave(filename, class_line_plot_normal, dpi = 300, width = 13, height = 5.5, units = "in")

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
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(y = "Relative abundance", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Class, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
class_line_plot_drought

filename <- file.path(figure_dir, 'FigSx-class-line-drought.png')
ggsave(filename, class_line_plot_drought, dpi = 300, width = 13, height = 5.5, units = "in")
filename <- file.path(figure_dir, 'FigSx-class-line-drought.pdf')
ggsave(filename, class_line_plot_drought, dpi = 300, width = 13, height = 5.5, units = "in")

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

####-------additional class anova----####
##class normal
df_longer_tannin_Norm <- df_longer_class %>%
  filter(Treatment == "Normal" &
           Class == "Tannin") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(rel_abundance = median(rel_abundance)) 

#lmer
m <- lm(rel_abundance ~ Zone * Timepoint,  data = df_longer_tannin_Norm)
anova(m)

# Response: rel_abundance
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Zone            3 27.648  9.2158  9.1184 0.0001192 ***
#   Timepoint       2  5.568  2.7841  2.7547 0.0766959 .  
# Zone:Timepoint  6 10.661  1.7768  1.7580 0.1349881    
# Residuals      37 37.395  1.0107  

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone )
pairs(em, adjust="tukey")

# contrast               estimate    SE df t.ratio p.value
# Bulk - Detritus          1.9046 0.410 37   4.640  0.0002
# Bulk - RhizoDet          1.8337 0.404 37   4.544  0.0003
# Bulk - Rhizosphere       1.3865 0.410 37   3.378  0.0090
# Detritus - RhizoDet     -0.0708 0.404 37  -0.175  0.9980
# Detritus - Rhizosphere  -0.5180 0.410 37  -1.262  0.5922
# RhizoDet - Rhizosphere  -0.4472 0.404 37  -1.108  0.6868

##class normal
df_longer_lignin_Norm <- df_longer_class %>%
  filter(Treatment == "Drought" &
           Class == "Lignin") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(rel_abundance = median(rel_abundance)) 

#lmer
m <- lm(rel_abundance ~ Zone * Timepoint,  data = df_longer_lignin_Norm)
anova(m)

# Response: rel_abundance
# Df  Sum Sq Mean Sq  F value    Pr(>F)    
# Zone            3 1179.80  393.27 110.6642 < 2.2e-16 ***
#   Timepoint       2   58.74   29.37   8.2644  0.001112 ** 
#   Zone:Timepoint  6   39.80    6.63   1.8665  0.113696    
# Residuals      36  127.93    3.55       

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone )
pairs(em, adjust="tukey")

# contrast               estimate   SE df t.ratio p.value
# Bulk - Detritus         -4.4415 0.77 36  -5.771 <0.0001
# Bulk - RhizoDet         -4.5061 0.77 36  -5.855 <0.0001
# Bulk - Rhizosphere       7.6614 0.77 36   9.955 <0.0001
# Detritus - RhizoDet     -0.0646 0.77 36  -0.084  0.9998
# Detritus - Rhizosphere  12.1030 0.77 36  15.726 <0.0001
# RhizoDet - Rhizosphere  12.1675 0.77 36  15.810 <0.0001

####----------Transformations-------####
#count line
trans_count_summary <- trans %>%
  group_by(Zone, Treatment, Timepoint, SampleID) %>%
  summarize(count = n()) %>%
  group_by(Zone, Treatment, Timepoint) %>%
  summarise(
    mean = mean(count, na.rm = TRUE),
    se   = sd(count, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) 

count_line_plot <- ggplot(trans_count_summary, 
                        aes(x = Timepoint, y =  mean, 
                            color = Zone, group = Zone)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = 0.1, linewidth = 0.4) +
  scale_color_manual(values = col_list_zone) +
  labs(y = "Transformation count", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Treatment) +
  theme_bw() +
  theme(text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

count_line_plot
ggsave(paste(figure_dir, "/FigSx_trans_count.png", sep = ""), dpi = 300, height = 5, width = 7, unit = "in")
ggsave(paste(figure_dir, "/FigSx_trans_count.pdf", sep = ""), dpi = 300, height = 5, width = 7, unit = "in")

#trans line
trans_group_summary <- trans %>%
  group_by(Zone, Treatment, Timepoint, SampleID, Transformation, Group) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(Zone, Treatment, Timepoint, SampleID) %>%           # regroup to just the sample level
  mutate(rel_abund = count / sum(count) * 100) %>%             # % of that sample's total transformations
  ungroup() %>%
  group_by(Zone, Treatment, Timepoint, Transformation, Group) %>%
  summarise(
    mean = mean(rel_abund, na.rm = TRUE),
    se   = sd(rel_abund, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

group_line_plot_normal <- ggplot(trans_group_summary %>% filter(Group == "Sugar" &
                                                                  Treatment == "Normal"), 
                          aes(x = Timepoint, y =  mean, 
                              color = Zone, group = Zone)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = 0.1, linewidth = 1) +
  scale_color_manual(values = col_list_zone) +
  scale_x_discrete(labels = c("4", "8", "12")) +
  labs(y = "Transformation count", x = "Timepoint", color = "Zone") +
  facet_wrap(~ Transformation, ncol = 10, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 14)) 

group_line_plot_normal
ggsave(paste(figure_dir, "/FigSx_trans_count.png", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")
ggsave(paste(figure_dir, "/FigSx_trans_count.pdf", sep = ""), dpi = 300, height = 4, width = 10, unit = "in")


##count normal anova
df_longer_count_Norm <- trans %>%
  filter(Treatment == "Normal") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(count = n()) 

#lm
m <- lm(count ~ Zone * Timepoint,  data = df_longer_count_Norm)
anova(m)

# Response: count
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# Zone            3 1.0811e+11 3.6037e+10 34.1724 9.287e-11 ***
#   Timepoint       2 3.3594e+09 1.6797e+09  1.5928 0.2169838    
# Zone:Timepoint  6 3.0909e+10 5.1515e+09  4.8850 0.0009122 ***
#   Residuals      37 3.9019e+10 1.0546e+09         

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")

# Timepoint = 4weeks:
#   contrast               estimate    SE df t.ratio p.value
# Rhizosphere - RhizoDet   -56540 23000 37  -2.462  0.0830
# Rhizosphere - Detritus   -59665 23000 37  -2.598  0.0615
# Rhizosphere - Bulk        94707 23000 37   4.124  0.0011
# RhizoDet - Detritus       -3124 23000 37  -0.136  0.9991
# RhizoDet - Bulk          151248 23000 37   6.587 <0.0001
# Detritus - Bulk          154372 23000 37   6.723 <0.0001
# 
# Timepoint = 8weeks:
#   contrast               estimate    SE df t.ratio p.value
# Rhizosphere - RhizoDet   -31038 21800 37  -1.425  0.4924
# Rhizosphere - Detritus   -43294 23000 37  -1.885  0.2517
# Rhizosphere - Bulk       105926 23000 37   4.613  0.0003
# RhizoDet - Detritus      -12256 21800 37  -0.563  0.9425
# RhizoDet - Bulk          136964 21800 37   6.287 <0.0001
# Detritus - Bulk          149220 23000 37   6.498 <0.0001
# 
# Timepoint = 12weeks:
#   contrast               estimate    SE df t.ratio p.value
# Rhizosphere - RhizoDet    -2404 23000 37  -0.105  0.9996
# Rhizosphere - Detritus   -84043 23000 37  -3.660  0.0042
# Rhizosphere - Bulk        -3253 23000 37  -0.142  0.9990
# RhizoDet - Detritus      -81638 23000 37  -3.555  0.0056
# RhizoDet - Bulk            -849 23000 37  -0.037  1.0000
# Detritus - Bulk           80790 23000 37   3.518  0.0062

##NOSC drought
df_longer_count_Drought <- trans %>%
  filter(Treatment == "Drought") %>%
  group_by(Zone,Timepoint, SampleID) %>%
  summarise(count = n()) 

#lm
m <- lm(count ~ Zone * Timepoint,  data = df_longer_count_Drought)
anova(m)

# Response: count
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# Zone            3 1.7741e+11 5.9137e+10 76.6938 1.053e-15 ***
#   Timepoint       2 7.9707e+09 3.9853e+09  5.1685   0.01063 *  
#   Zone:Timepoint  6 1.3876e+10 2.3127e+09  2.9993   0.01761 *  
#   Residuals      36 2.7759e+10 7.7108e+08         

# Post-hoc pairwise contrasts at each timepoint
em    <- emmeans(m, ~ Zone | Timepoint)
pairs(em, adjust="tukey")
# Timepoint = 4weeks:
#   contrast               estimate    SE df t.ratio p.value
# Rhizosphere - RhizoDet   -10825 19600 36  -0.551  0.9456
# Rhizosphere - Detritus    -7487 19600 36  -0.381  0.9808
# Rhizosphere - Bulk       116802 19600 36   5.949 <0.0001
# RhizoDet - Detritus        3338 19600 36   0.170  0.9982
# RhizoDet - Bulk          127627 19600 36   6.500 <0.0001
# Detritus - Bulk          124289 19600 36   6.330 <0.0001
# 
# Timepoint = 8weeks:
#   contrast               estimate    SE df t.ratio p.value
# Rhizosphere - RhizoDet   -38109 19600 36  -1.941  0.2295
# Rhizosphere - Detritus   -15653 19600 36  -0.797  0.8552
# Rhizosphere - Bulk       120127 19600 36   6.118 <0.0001
# RhizoDet - Detritus       22456 19600 36   1.144  0.6654
# RhizoDet - Bulk          158236 19600 36   8.059 <0.0001
# Detritus - Bulk          135780 19600 36   6.915 <0.0001
# 
# Timepoint = 12weeks:
#   contrast               estimate    SE df t.ratio p.value
# Rhizosphere - RhizoDet   -61248 19600 36  -3.119  0.0179
# Rhizosphere - Detritus  -107019 19600 36  -5.450 <0.0001
# Rhizosphere - Bulk        88306 19600 36   4.497  0.0004
# RhizoDet - Detritus      -45771 19600 36  -2.331  0.1098
# RhizoDet - Bulk          149554 19600 36   7.617 <0.0001
# Detritus - Bulk          195325 19600 36   9.948 <0.0001

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













