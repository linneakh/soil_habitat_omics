
# ============================================================
# Linnea Hernandez
# 8/15/26
# Hernandez163@llnl.gov
# 
# WGCNA analysis: Drought treatment
# ============================================================


# ----------------------------
# 1. Load packages and functions
# ----------------------------

library(WGCNA)
library(flashClust)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(tidyverse)
library(qvalue)
library(pheatmap)
library(clusterProfiler)
library(KEGGREST)
library(lmerTest)
library(nlme)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

source("./scripts/functions/wgcna_functions.R")


# ----------------------------
# 2. Define directories
# ----------------------------

Dir.f <- "./figures/Fig5_FigSx_WGCNA/drought/"
Dir.o <- "./output/WGCNA/drought/"

dir.create(Dir.f, recursive = TRUE, showWarnings = FALSE)
dir.create(Dir.o, recursive = TRUE, showWarnings = FALSE)


# ----------------------------
# 3. Define analysis parameters
# ----------------------------

#softPower        <- 9
minModuleSize    <- 40
deepSplitValue   <- 2
mergeCutHeight   <- 0.30
moduleSigAlpha   <- 0.05
kME_threshold    <- 0.7
minGenesFiltered <- 5
minSize          <- 5
cytoscapeThresh  <- 0.01
cytoscapeEdgeMin <- 0.20
dpi=300
h=5
w=5
size=10


# ----------------------------
# 4. Define factor levels and annotation colors
# ----------------------------

zone_levels <- c("Bulk", "Rhizo", "RhizoDetritus", "Detritus")
time_levels <- c("4weeks", "8weeks", "12weeks")

ann_color <- list(
  Zone = c(
    Bulk = "darkred",
    Rhizo = "green4",
    Detritus = "blue",
    RhizoDetritus = "purple"
  ),
  Timepoint = c(
    `4weeks` = "pink",
    `8weeks` = "orange",
    `12weeks` = "yellow"
  )
)


# ----------------------------
# 5. Load and preprocess expression data
# ----------------------------

datExpr.0 <- read.csv(
  "./output/salazar_drought_vst.csv",
  header = TRUE
) %>%
  column_to_rownames(var = "X")


# Filter features present in at least 20% of samples
limit <- 0.2 * nrow(datExpr.0)

datExpr <- datExpr.0 %>%
  arrange(rownames(.)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sum.count = rowSums(.)) %>%
  mutate(gr.zero = rowSums(. > 0)) %>%
  filter(gr.zero > limit) %>%
  dplyr::select(-c(gr.zero, sum.count))


# ----------------------------
# 6. Format sample trait data
# ----------------------------

datTraits <- datExpr %>%
  rownames_to_column(var = "SampleID") %>%
  dplyr::select(SampleID) %>%
  separate_wider_delim(
    SampleID,
    delim = ".",
    names = c("SampleID_num", "Treatment", "Zone", "Timepoint"),
    cols_remove = FALSE
  ) %>%
  mutate(
    Time = str_remove(Timepoint, "weeks"),
    Time = as.numeric(Time),
    RhizoDet = ifelse(Zone == "RhizoDetritus", 1, 0),
    Rhizo = ifelse(Zone == "Rhizo", 1, 0),
    Detritus = ifelse(Zone == "Detritus", 1, 0),
    Bulk = ifelse(Zone == "Bulk", 1, 0)
  ) %>%
  dplyr::select(-Timepoint, -Zone, -Treatment, -SampleID_num) %>%
  column_to_rownames(var = "SampleID")


# ----------------------------
# 7. Quality control
# ----------------------------

gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  datTraits <- datTraits[rownames(datExpr), , drop = FALSE]
}


# ----------------------------
# 8. Sample clustering and outlier filtering (step 1)
# ----------------------------

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)

png(
  file = paste0(
    Dir.f,
    "other/step_1_sample_clustering_outliers.png"
  ),
  width = 12,
  height = 9,
  units = "in",
  res = 300
)

par(cex = 0.6)
par(mar = c(0, 4, 2, 0))

plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = ""
)

abline(h = 110, col = "red")

dev.off()


clust <- cutreeStatic(
  sampleTree,
  cutHeight = 30000,
  minSize = 5
)

keepSamples <- (clust == 1)

datExpr <- datExpr[keepSamples, , drop = FALSE]
datTraits <- datTraits[rownames(datExpr), , drop = FALSE]

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)


# ----------------------------
# 9. Sample tree visualization and soft-threshold selection (step 2)
# ----------------------------

plot_sample_tree(
  datExpr,
  datTraits,
  prefix = "other/step_2"
)

sft <- plot_soft_threshold(
  datExpr,
  prefix = "other/step_2",
  save_plot = TRUE,
  networkType = "signed"
)

print(sft)


# ----------------------------
# 10. Construct the weighted gene co-expression network
# ----------------------------

softPower = 14

TOM <- TOMsimilarityFromExpr(
  datExpr,
  power = softPower,
  networkType = "signed"
)

dissTOM <- 1 - TOM

if (any(!is.finite(dissTOM))) {
  stop(
    "dissTOM contains NA/NaN/Inf values. ",
    "Check datExpr for missing or zero-variance genes."
  )
}


# ----------------------------
# 11. Gene clustering and dynamic module detection (step 3)
# ----------------------------

geneTree <- hclust(
  as.dist(dissTOM),
  method = "average"
)

plot_gene_tree(
  geneTree,
  prefix = "other/step_3"
)


dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = deepSplitValue,
  pamRespectsDendro = TRUE,
  minClusterSize = minModuleSize
)

dynamicColors <- labels2colors(dynamicMods)


# ----------------------------
# 12. Module eigengenes and module merging (step 4)
# ----------------------------

plot_module_tree(
  geneTree,
  colors_mat = dynamicColors,
  labels = "Dynamic Tree Cut",
  prefix = "other/step_4_dynamic",
  main_text = "Gene dendrogram and module colors"
)

MEList <- moduleEigengenes(
  datExpr,
  colors = dynamicColors
)

MEs <- orderMEs(MEList$eigengenes)

plot_me_clustering(
  MEs,
  cutHeight = mergeCutHeight,
  prefix = "other/step_4"
)

merge <- mergeCloseModules(
  datExpr,
  dynamicColors,
  cutHeight = mergeCutHeight,
  verbose = 3
)

mergedColors <- merge$colors
mergedMEs <- orderMEs(merge$newMEs)

write.csv(
  mergedMEs,
  paste0(Dir.o, "full_merged_mEs.csv")
)

write.csv(
  as.data.frame(scale(mergedMEs)),
  paste0(Dir.o, "full_merged_mEs_zscaled.csv")
)


# ----------------------------
# 13. Visualize and save merged modules (step 5)
# ----------------------------

plot_module_tree(
  geneTree,
  colors_mat = cbind(dynamicColors, mergedColors),
  labels = c("Dynamic Tree Cut", "Merged dynamic"),
  prefix = "other/step_5_merged",
  main_text = "Gene dendrogram and merged module colors"
)

filename.network <- paste0(
  Dir.o,
  "networkConstruction-auto-step-log.RData"
)

save(
  MEs,
  dynamicMods,
  dynamicColors,
  geneTree,
  datExpr,
  mergedMEs,
  mergedColors,
  TOM,
  datTraits,
  file = filename.network
)


# ----------------------------
# 14. Generate full analysis outputs
# ----------------------------

full_keys <- write_module_key(
  datExpr,
  mergedColors,
  prefix = "full"
)

plot_me_sample_heatmap(
  merge$oldMEs,
  datTraits,
  prefix = "other/full_oldMEs"
)

plot_me_sample_heatmap(
  mergedMEs,
  datTraits,
  prefix = "other/full_newMEs"
)

full_geneInfo <- compute_gene_membership_table(
  datExpr = datExpr,
  datTraits = datTraits,
  colors = mergedColors,
  MEs = mergedMEs,
  prefix = "full"
)


# ----------------------------
# 15. Define module ordering for figures
# ----------------------------

order_MEs <- c(
  "MEtan",
  "MEgreen",
  "MEblack",
  "MEpurple",
  "MEblue",
  "MEsalmon",
  "MEturquoise",
  "MEbrown"
)

order_MEs_colors <- c(
  "tan",
  "green",
  "black",
  "purple",
  "blue",
  "salmon",
  "turquoise",
  "brown"
)

order_MEs.fig <- c(
  "MEsalmon",
  "MEturquoise",
  "MEblue"
)

order_MEs_colors.fig <- c(
  "salmon",
  "turquoise",
  "blue"
)


# ----------------------------
# 16. Module-trait relationships
# ----------------------------

full_module_trait <- compute_module_trait_heatmap(
  MEs = mergedMEs,
  datTraits = datTraits,
  prefix = "full",
  module_order = order_MEs
)

full_module_trait_primary_fig <- plot_module_trait_heatmap_subset(
  full_module_trait,
  trait_subset = c(
    "Bulk",
    "Rhizo",
    "RhizoDet",
    "Detritus"
  ),
  module_subset = order_MEs.fig,
  prefix = "other/subset_primary_fig"
)

plot_module_trait_with_gene_bars(
  full_results = full_module_trait,
  df_key = full_keys$df_key,
  trait_subset = c(
    "Detritus",
    "RhizoDet",
    "Rhizo",
    "Bulk"
  ),
  module_subset = order_MEs_colors.fig,
  prefix = "other/full_combined"
)


#set exudate and litter components to plot:
common_exudates <- c(
  "Amino acids", "Aromatics/Phenolics",
  #"Flavonoids", 
  "Quaternary amines",
  "Mucilage sugars",
  "Organic acids",
  "Simple sugars"
)

common_litter <- c(
  "Cellulose", 
  "Chitin", 
  #"Cutin/Suberin", 
  "Hemicellulose", "Lignin",
  "Pectin", "Starch"
)


full_module_trait_figure_w_genes <- plot_module_trait_with_gene_bars(
  full_results = full_module_trait,
  df_key = full_keys$df_key,
  trait_subset = c(
    "Detritus",
    "RhizoDet",
    "Rhizo",
    "Bulk"
  ),
  module_subset = order_MEs_colors,
  prefix = "other/full_heatmap_genes_combined",
  exudate_levels = common_exudates,
  litter_levels = common_litter
)


subset_module_trait_figure_w_genes <- plot_module_trait_with_gene_bars(
  full_results = full_module_trait_primary_fig,
  df_key = full_keys$df_key,
  trait_subset = c(
    "Detritus",
    "RhizoDet",
    "Rhizo"
  ),
  module_subset = order_MEs_colors.fig,
  prefix = "Fig5",
  exudate_levels = common_exudates,
  litter_levels = common_litter
)

#save tables of root exudate and litter degradation genes for table Sx
write.csv(as.data.frame(subset_module_trait_figure_w_genes[1]), "./output/WGCNA/drought/TableSx_module_exudate.csv")
write.csv(as.data.frame(subset_module_trait_figure_w_genes[2]), "./output/WGCNA/drought/TableSx_litter.csv")


# ----------------------------
# 17. Eigengene visualization
# ----------------------------

plot_eigengenes_by_module(
  MEs = mergedMEs,
  datTraits = datTraits,
  prefix = "other/full",
  ncol = 4
)

w = w + 5
h = h + 10


# ----------------------------
# 18. Pathway analysis
# ----------------------------

run_pathway_plots(
  gene_table = full_geneInfo %>%
    dplyr::filter(ModuleColor != "grey"),
  prefix = "pathway/full"
)


# ----------------------------
# 19. Eigengene statistical analyses
# ----------------------------

full_eigengene_stats <- run_eigengene_zone_time_stats(
  MEs = mergedMEs,
  datTraits = datTraits,
  prefix = "full"
)

res_compare <- compare_modules_within_treatment(
  mergedMEs,
  datTraits,
  treatment_name = "Drought",
  prefix = "full"
)


# ----------------------------
# 20. Publication figure eigengene plots
# ----------------------------

plot_eigengenes_gg(
  mergedMEs[
    ,
    c("MEturquoise", "MEblue")
  ],
  datTraits,
  order_MEs = order_MEs_colors.fig,
  facet_by = "Module",
  color_by = "Module",
  prefix = "Fig5",
  ncol = 2,
  w = w,
  h = h - 2
)


