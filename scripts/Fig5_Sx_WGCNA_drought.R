# Linnea Honeker
# 2/15/22
# linneah@arizona.edu
# Objective: To perform WGCNA analysis on soil pyruvate metaT data paired with metabolomics/VOCs

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

source("./scripts/WGCNA/extra_functions_wgcna.R")

# ----------------------------
# Directories
# ----------------------------
Dir.f <- "./figures/all-kingdoms-taxize/bacteria/salazar/wgcna/drought/high_sp/"
Dir.o <- "./output/all_kingdoms/bacteria_only/salazar/wgcna/drought/high_sp/"

dir.create(Dir.f, recursive = TRUE, showWarnings = FALSE)
dir.create(Dir.o, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Parameters
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

source("./scripts/WGCNA/extra_functions_2.R")


# ----------------------------
# Load and preprocess data
# ----------------------------

datExpr.0 <- read.csv(
  "./output/all_kingdoms/bacteria_only/salazar/normalized_per_cell_transcript_drought_vst.csv",
  header = TRUE
) %>%
  column_to_rownames(var = "X")

# datExpr.0 <- read.csv(
#   "./output/all_kingdoms/bacteria_only/deseq2/vst.bacteria_new_database.csv",
#   header = TRUE) %>%
#   column_to_rownames(var = "X")

limit <- 0.2 * nrow(datExpr.0)

datExpr <- datExpr.0 %>%
  arrange(rownames(.)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sum.count = rowSums(.)) %>%
  mutate(gr.zero = rowSums(. > 0)) %>%
  filter(gr.zero > limit) %>%
  dplyr::select(-c(gr.zero, sum.count))

datTraits <- datExpr %>%
  rownames_to_column(var = "SampleID") %>%
  dplyr::select(SampleID) %>%
  separate_wider_delim(SampleID, 
                       delim = ".", 
                       names = c("SampleID_num", "Treatment", "Zone", "Timepoint"),
                       cols_remove = FALSE) %>%
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

# QC
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  datTraits <- datTraits[rownames(datExpr), , drop = FALSE]
}

# Sample clustering and outlier filtering
sampleTree <- hclust(dist(datExpr), method = "average")
png(file = paste0(Dir.f, "step_1_sample_clustering_outliers.png"),
    width = 12, height = 9, units = "in", res = 300)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
abline(h = 110, col = "red")
dev.off()

clust <- cutreeStatic(sampleTree, cutHeight = 30000, minSize = 5)
keepSamples <- (clust == 1 | clust == 2)

datExpr <- datExpr[keepSamples, , drop = FALSE]
datTraits <- datTraits[rownames(datExpr), , drop = FALSE]

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

plot_sample_tree(datExpr, datTraits, prefix = "step_2")
sft  <- plot_soft_threshold(datExpr, prefix = "step_2", save_plot = TRUE, networkType = "signed")
print(sft)
# ----------------------------
# Network construction
# ----------------------------
softPower = 14
TOM <- TOMsimilarityFromExpr(datExpr, power = softPower, networkType = "signed")
dissTOM <- 1 - TOM

if (any(!is.finite(dissTOM))) {
  stop("dissTOM contains NA/NaN/Inf values. Check datExpr for missing or zero-variance genes.")
}

geneTree <- hclust(as.dist(dissTOM), method = "average")
plot_gene_tree(geneTree, prefix = "step_3")

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = deepSplitValue,
  pamRespectsDendro = TRUE,
  minClusterSize = minModuleSize
)

dynamicColors <- labels2colors(dynamicMods)

plot_module_tree(
  geneTree,
  colors_mat = dynamicColors,
  labels = "Dynamic Tree Cut",
  prefix = "step_4_dynamic",
  main_text = "Gene dendrogram and module colors"
)

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- orderMEs(MEList$eigengenes)

plot_me_clustering(MEs, cutHeight = mergeCutHeight, prefix = "step_4")

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = mergeCutHeight, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- orderMEs(merge$newMEs)

write.csv(mergedMEs, paste0(Dir.o, "full_merged_mEs.csv"))
write.csv(as.data.frame(scale(mergedMEs)), paste0(Dir.o, "full_merged_mEs_zscaled.csv"))

plot_module_tree(
  geneTree,
  colors_mat = cbind(dynamicColors, mergedColors),
  labels = c("Dynamic Tree Cut", "Merged dynamic"),
  prefix = "step_5_merged",
  main_text = "Gene dendrogram and merged module colors"
)

filename.network <- paste0(Dir.o, "networkConstruction-auto-step-log.RData")
save(MEs, dynamicMods, dynamicColors, geneTree, datExpr, mergedMEs, mergedColors, TOM, datTraits,
     file = filename.network)

# ----------------------------
# FULL analysis outputs
# ----------------------------
full_keys <- write_module_key(datExpr, mergedColors, prefix = "full")
plot_me_sample_heatmap(merge$oldMEs, datTraits, prefix = "full_oldMEs")
plot_me_sample_heatmap(mergedMEs, datTraits, prefix = "full_newMEs")

full_geneInfo <- compute_gene_membership_table(
  datExpr = datExpr,
  datTraits = datTraits,
  colors = mergedColors,
  MEs = mergedMEs,
  prefix = "full"
)


#define order of MEs and colors for plotting
order_MEs <- c("MEtan", "MEgreen", "MEblack", "MEpurple", "MEblue", "MEsalmon", "MEturquoise", "MEbrown")

order_MEs_colors <-  c("tan", "green", "black", "purple", "blue", "salmon", "turquoise", "brown")

order_MEs.fig <- c("MEsalmon", "MEturquoise", "MEblue")

order_MEs_colors.fig <- c("salmon", "turquoise","blue" )

full_module_trait <- compute_module_trait_heatmap(
  MEs = mergedMEs,
  datTraits = datTraits,
  prefix = "full",
  module_order = order_MEs
)

full_module_trait_primary_fig <- plot_module_trait_heatmap_subset(
  full_module_trait,
  trait_subset = c("Bulk", "Rhizo", "RhizoDet", "Detritus"),
  module_subset = order_MEs.fig,
  prefix = "subset_primary_fig")

plot_module_trait_with_gene_bars(
  full_results  = full_module_trait,
  df_key        = full_keys$df_key,
  trait_subset  = c("Detritus", "RhizoDet","Rhizo", "Bulk"),
  module_subset = order_MEs_colors.fig,  
  prefix        = "full_heatmap_genes_combined"
)

plot_module_trait_with_gene_bars(
  full_results  = full_module_trait_primary_fig,
  df_key        = full_keys$df_key,
  trait_subset  = c("Detritus", "RhizoDet","Rhizo"),
  module_subset = order_MEs_colors.fig,  
  prefix        = "full_heatmap_genes_combined_subset"
)



plot_eigengenes_by_module(
  MEs = mergedMEs,
  datTraits = datTraits,
  prefix = "full",
  ncol=4
)

#redefine w and h
w = w+5
h = h+10
run_pathway_plots(
  gene_table = full_geneInfo %>% dplyr::filter(ModuleColor != "grey"),
  prefix = "full"
)


plot_exudate_litter_summary(
  df_key = full_keys$df_key %>% dplyr::filter(moduleColors != "grey"),
  prefix = "full"
)

export_cytoscape_filtered(
  TOM = TOM,
  datExpr = datExpr,
  colors = setNames(mergedColors, colnames(datExpr)),
  prefix = "full"
)

full_eigengene_stats <- run_eigengene_zone_time_stats(
  MEs = mergedMEs,
  datTraits = datTraits,
  prefix = "full"
)


res_compare <- compare_modules_within_treatment(
  mergedMEs, datTraits,
  treatment_name = "Drought",
  prefix = "full"
)

plot_eigengenes_gg(
  mergedMEs,
  datTraits,
  facet_by = c("Zone","Timepoint"),
  color_by = "Module",
  prefix = "full_grid_drought",
  ncol=3,
  w = w+10,
  h = h+5)

plot_eigengenes_gg(
  mergedMEs,
  datTraits,
  facet_by = c("Zone"),
  color_by = "Module",
  prefix = "full_grid_zone_drought",
  w = w+1,
  h = h+1)

plot_eigengenes_gg(
  mergedMEs[,c("MEturquoise", "MEblue")],
  datTraits,
  order_MEs = order_MEs_colors.fig,
  facet_by = "Module",
  color_by = "Module",
  prefix = "full_grid_drought_no_facet",
  ncol = 2,
  w = w,
  h = h-2)

plot_eigengenes_gg(
  mergedMEs[,c("MEturquoise", "MEblue")],
  datTraits,
  facet_by = c("Zone"),
  color_by = "Module",
  prefix = "full_grid_zone_drought_deg",
  ncol=4,
  w = w+1,
  h = h+1)


# ----------------------------
# FILTERED analysis outputs
# ----------------------------
filteredColors <- make_filtered_colors(datExpr, full_geneInfo)
filteredColors <- filteredColors[colnames(datExpr)]

MEList.filtered <- moduleEigengenes(datExpr, colors = filteredColors)
filteredMEs <- orderMEs(MEList.filtered$eigengenes)

write.csv(filteredMEs, paste0(Dir.o, "filtered_module_eigengenes.csv"))

plot_module_tree(
  geneTree,
  colors_mat = cbind(dynamicColors, mergedColors, filteredColors),
  labels = c("Dynamic Tree Cut", "Merged dynamic", "Filtered modules"),
  prefix = "step_6_filtered",
  main_text = "Gene dendrogram and filtered module colors"
)

filtered_keys <- write_module_key(datExpr, filteredColors, prefix = "filtered")
plot_me_sample_heatmap(filteredMEs, datTraits, prefix = "filtered_newMEs")

filtered_geneInfo <- compute_gene_membership_table(
  datExpr = datExpr,
  datTraits = datTraits,
  colors = filteredColors,
  MEs = filteredMEs,
  prefix = "filtered"
)



# order_MEs <- c("MEturquoise", "MEblack", "MEcyan", "MEmagenta", "MEgreen")

filtered_module_trait <- compute_module_trait_heatmap(
  MEs = filteredMEs,
  datTraits = datTraits,
  prefix = "filtered",
  module_order = order_MEs
)

plot_eigengenes_by_module(
  MEs = filteredMEs,
  datTraits = datTraits,
  prefix = "filtered"
)

# run_pathway_plots(
#   gene_table = filtered_geneInfo %>%
#     dplyr::filter(ModuleColor != "grey", MM.assigned.significant),
#   prefix = "filtered"
# )

plot_exudate_litter_summary(
  df_key = filtered_keys$df_key %>% dplyr::filter(moduleColors != "grey"),
  prefix = "filtered"
)

export_cytoscape_filtered(
  TOM = TOM,
  datExpr = datExpr,
  colors = filteredColors,
  prefix = "filtered"
)

filtered_eigengene_stats <- run_eigengene_zone_time_stats(
  MEs = filteredMEs,
  datTraits = datTraits,
  prefix = "filtered"
)


res_compare <- compare_modules_within_treatment(
  filteredMEs, datTraits,
  treatment_name = "Drought",
  prefix = "filtered"
)

plot_eigengenes_gg(
  filteredMEs[,order_MEs],
  datTraits,
  facet_by = c("Zone","Timepoint"),
  color_by = "Module",
  prefix = "filtered_grid",
  ncol=3,
  w = w+10,
  h = h+5)

plot_eigengenes_gg(
  filteredMEs[,order_MEs],
  datTraits,
  facet_by = c("Zone"),
  color_by = "Module",
  prefix = "filtered_grid_zone",
  w = w+1,
  h = h+1)

plot_eigengenes_gg(
  filteredMEs[,c("MEblue", "MEturquoise")],
  datTraits,
  facet_by = c("Zone"),
  color_by = "Module",
  prefix = "filtered_grid_zone_drought_deg",
  ncol=4,
  w = w+1,
  h = h+1)


# ----------------------------
# Save key objects
# ----------------------------
save(
  datExpr, datTraits, TOM, dissTOM, geneTree,
  dynamicMods, dynamicColors,
  mergedColors, mergedMEs,
  filteredColors, filteredMEs,
  full_geneInfo, filtered_geneInfo,
  full_eigengene_stats, filtered_eigengene_stats,
  file = paste0(Dir.o, "wgcna_full_and_filtered_workspace.RData")
)

