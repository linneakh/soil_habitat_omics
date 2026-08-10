##Linnea Honeker
##NMDS and PCA on metaT data from IMG
##2/26/24

library(ggfortify)
library(tidyverse)
library(ggrepel)
library(viridis)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)
library(stats)
library(emmeans)


# define colors
col_list_det = c("orange", "blue")
col_list_rhiz = c("green", "blue")
col_list_treatment = c("darkred", "darkgreen")


# parameters for plots
h = 3
w = 6
res = 300
size = 10

###inport data
maom.0 <- read.csv("./data/maom.csv", header = TRUE)

#change harvest times to timepoints and put treatments in order and rename columns
maom <- maom.0 %>%
  mutate(Timepoint = case_when(
    Harvest == "H1" ~ "4weeks",
    Harvest == "H2" ~ "8weeks",
    Harvest == "H3" ~ "12weeks"
  )) %>%
  mutate(Timepoint.num = case_when(
    Harvest == "H1" ~ 4,
    Harvest == "H2" ~ 8,
    Harvest == "H3" ~ 12
  )) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks"))) %>%
  mutate(Moisture = factor(Moisture, levels = c("Normal", "Drought"))) %>%
  rename(Zone = Treatment, Treatment = Moisture)
write.csv(maom, "./output/maom-num.csv")


##line plots
#plot 13C-detritus
maom.det <- maom %>%
  as.data.frame() %>%
  filter(Isotope.Loc == "Detritus.13C") %>%
  mutate(Treatment = factor(Treatment, c("Normal", "Drought")))

## 1. Summarise to get one value per Treatment × Zone × Timepoint.num
maom.sum <- maom.det %>%
  group_by(Treatment, Zone, Timepoint.num) %>%
  summarise(
    mean_Clay = mean(Clay.ug.13C, na.rm = TRUE),
    sd_Clay   = sd(Clay.ug.13C,   na.rm = TRUE),
    n         = dplyr::n(),
    se_Clay   = sd_Clay / sqrt(n)
  ) %>%
  ungroup() 

## 2. Plot: one point + error bars per group, lines for each Zone
plot.det <- ggplot(
  maom.sum,
  aes(x = Timepoint.num,
      y = mean_Clay,
      color = Zone,
      group = Zone)      # one line per zone
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_Clay - se_Clay,
        ymax = mean_Clay + se_Clay),
    width = 0.2
  ) +
  ylim(c(20,100)) +
  scale_color_manual(values = col_list_det) +
  facet_wrap(~ Treatment) +
  labs(x= "Harvest timepoint", 
       y = expression({}^{13}*C~"--"*MAOM~"(mg" * {}^{13}*C~g~soil^{-1} * ")")) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )

plot(plot.det)

ggsave("./figures/Fig.2-det-line.png", dpi = res, w=w, h=h, units = "in")
ggsave("./figures/Fig.2-det-line.pdf", dpi = res, w=w, h=h, units = "in")


#plot 13C-rhizo
maom.rhiz <- maom %>%
  as.data.frame() %>%
  filter(Isotope.Loc == "Rhizo.13C")

## 1. Summarise to get one value per Treatment × Zone × Timepoint.num
maom.sum <- maom.rhiz %>%
  group_by(Treatment, Zone, Timepoint.num) %>%
  summarise(
    mean_Clay = mean(Clay.ug.13C, na.rm = TRUE),
    sd_Clay   = sd(Clay.ug.13C,   na.rm = TRUE),
    n         = dplyr::n(),
    se_Clay   = sd_Clay / sqrt(n)
  ) %>%
  ungroup()

## 2. Plot: one point + error bars per group, lines for each Zone
plot.rhiz <- ggplot(
  maom.sum,
  aes(x = Timepoint.num,
      y = mean_Clay,
      color = Zone,
      group = Zone)      # one line per zone
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_Clay - se_Clay,
        ymax = mean_Clay + se_Clay),
    width = 0.2
  ) +
  scale_color_manual(values = col_list_rhiz) +
  facet_wrap(~ Treatment) +
  labs(x= "Harvest timepoint", 
       y = expression({}^{13}*C~"-"*MAOM~"(mg " * {}^{13}*C~g~soil^{-1} * ")")) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.002, 0.998),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )

plot(plot.rhiz)
ggsave("./figures/Fig.2-rhizo-line.png", dpi = res, w=w, h=h, units = "in")
ggsave("./figures/Fig.2-rhizo-line.pdf", dpi = res, w=w, h=h, units = "in")



#anova on linear model
#det normal only
maom.det.norm <- maom.det %>%
  filter(Treatment == "Normal")

m_lm <- lm(Clay.ug.13C ~ Zone * Timepoint, data = maom.det.norm)
anova(m_lm)
# esponse: Clay.ug.13C
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Zone            1 2656.5 2656.54 10.5697 0.002838 ** 
#   Timepoint       2 4647.0 2323.50  9.2446 0.000745 ***
#   Zone:Timepoint  2  214.1  107.06  0.4260 0.657028    
# Residuals      30 7540.0  251.33    

em <- emmeans(m_lm, ~ Timepoint)
pairs(em, adjust="tukey") 
# contrast         estimate   SE df t.ratio p.value
# 4weeks - 8weeks      2.63 6.47 30   0.407  0.9132
# 4weeks - 12weeks   -22.68 6.47 30  -3.504  0.0041
# 8weeks - 12weeks   -25.31 6.47 30  -3.910  0.0014

#det drought only
maom.det.dr <- maom.det %>%
  filter(Treatment == "Drought")

m_lm <- lm(Clay.ug.13C ~ Zone * Timepoint, data = maom.det.dr)
anova(m_lm)
# Response: Clay.ug.13C
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Zone            1  538.2  538.24  2.9212 0.097749 . 
# Timepoint       2 2942.8 1471.40  7.9858 0.001658 **
#   Zone:Timepoint  2 2927.1 1463.56  7.9432 0.001704 **
#   Residuals      30 5527.6  184.25        

em <- emmeans(m_lm, ~ Zone | Timepoint)
pairs(em, adjust="tukey") 
# Timepoint = 4weeks:
#   contrast                  estimate   SE df t.ratio p.value
# Detritus - Rhizo.Detritus    -1.12 7.84 30  -0.142  0.8876
# 
# Timepoint = 8weeks:
#   contrast                  estimate   SE df t.ratio p.value
# Detritus - Rhizo.Detritus    -8.56 7.84 30  -1.092  0.2836
# 
# Timepoint = 12weeks:
#   contrast                  estimate   SE df t.ratio p.value
# Detritus - Rhizo.Detritus    32.87 7.84 30   4.195  0.0002




#rhizo
#rhizo normal only
maom.rhiz.norm <- maom.rhiz %>%
  filter(Treatment == "Normal")

m_lm <- lm(Clay.ug.13C ~ Zone * Timepoint, data = maom.rhiz.norm)
anova(m_lm)
# Response: Clay.ug.13C
# Df  Sum Sq Mean Sq F value   Pr(>F)   
# Zone            1     0.1    0.08  0.0002 0.988826   
# Timepoint       2  5032.5 2516.25  6.2965 0.005071 **
#   Zone:Timepoint  2   135.2   67.59  0.1691 0.845174   
# Residuals      31 12388.3  399.62    

em <- emmeans(m_lm, ~ Timepoint)
pairs(em, adjust="tukey") 
# contrast         estimate   SE df t.ratio p.value
# 4weeks - 8weeks     -7.91 8.01 31  -0.987  0.5902
# 4weeks - 12weeks   -27.82 8.01 31  -3.471  0.0043
# 8weeks - 12weeks   -19.90 8.16 31  -2.439  0.0525

#rhizo drought only
maom.rhiz.dr <- maom.rhiz %>%
  filter(Treatment == "Drought")

m_lm <- lm(Clay.ug.13C ~ Zone * Timepoint, data = maom.rhiz.dr)
anova(m_lm)
# Response: Clay.ug.13C
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Zone            1  907.71  907.71 16.5204 0.0002589 ***
#   Timepoint       2   59.35   29.68  0.5401 0.5874723    
# Zone:Timepoint  2  281.30  140.65  2.5598 0.0917159 .  
# Residuals      35 1923.06   54.94              



