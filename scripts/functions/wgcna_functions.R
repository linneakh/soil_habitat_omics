#define colors
colors <- c(
  "black", "gray", "coral4", "coral",
  "chartreuse3", "darkseagreen1", "blue", "lightblue",
  "yellow4", "yellow", "darkorchid", "plum2",
  "darkred", "darksalmon", "green4",
  "greenyellow", "orange", "moccasin",
  "hotpink4", "lightpink", "lightblue4", "lightcyan3",
  "lightslateblue", "lightsteelblue1", "navy",
  "darkgray", "brown", "cornflowerblue", "darkgoldenrod",
  "brown3", "aliceblue", "aquamarine",
  "beige", "bisque", "blue2", "blueviolet", "cyan", "darkblue",
  "chocolate", "aquamarine3", "darkcyan", "deeppink", "black", "gray", "coral4", "coral",
  "chartreuse3", "darkseagreen1", "blue", "lightblue",
  "yellow4", "yellow", "darkorchid", "plum2",
  "darkred", "darksalmon", "green4",
  "greenyellow", "orange", "moccasin",
  "hotpink4", "lightpink", "lightblue4", "lightcyan3",
  "lightslateblue", "lightsteelblue1", "navy",
  "darkgray", "brown", "cornflowerblue", "darkgoldenrod",
  "brown3", "aliceblue", "aquamarine",
  "beige", "bisque", "blue2", "blueviolet", "cyan", "darkblue",
  "chocolate", "aquamarine3", "darkcyan", "deeppink"
)


# =============================================================================
# WGCNA & CLUSTERPROFILE FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# 1. DATA PREPARATION & HELPER FUNCTIONS
# -----------------------------------------------------------------------------

safe_qvalues <- function(pvec) {
  qobj <- try(qvalue::qvalue(p = pvec), silent = TRUE)
  if (inherits(qobj, "try-error")) {
    return(p.adjust(pvec, method = "BH"))
  }
  qobj$qvalues
}

make_sample_annotation <- function(datTraits) {
  datTraits %>%
    rownames_to_column(var = "SampleID") %>%
    separate_wider_delim(SampleID,
                         delim = ".",
                         names = c("ID", "Treatment", "Zone", "Timepoint"),
                         cols_remove = FALSE) %>%
    mutate(
      Zone = factor(Zone, levels = zone_levels),
      Timepoint = factor(Timepoint, levels = time_levels)
    ) %>%
    column_to_rownames(var = "SampleID") %>%
    dplyr::select(Zone, Timepoint)
}

make_filtered_colors <- function(datExpr, geneInfo) {
  filteredColors <- setNames(geneInfo$AssignedModule, geneInfo$Gene)
  filteredColors <- filteredColors[colnames(datExpr)]
  
  keepGenes <- geneInfo$Gene[geneInfo$MM.assigned.sig.and.strong]
  filteredColors[!(names(filteredColors) %in% keepGenes)] <- "grey"
  
  tabFilt <- table(filteredColors)
  smallMods <- names(tabFilt)[tabFilt < minGenesFiltered & names(tabFilt) != "grey"]
  filteredColors[filteredColors %in% smallMods] <- "grey"
  
  filteredColors
}

# -----------------------------------------------------------------------------
# 2. NETWORK CONSTRUCTION — SOFT THRESHOLD & GENE TREE
# -----------------------------------------------------------------------------

plot_soft_threshold <- function(datExpr, prefix = NULL, save_plot = TRUE, networkType = "signed") {
  powers <- c(1:10, seq(from = 12, to = 20, by = 2))
  sft <- pickSoftThreshold(datExpr,
                           powerVector = powers,
                           verbose = 5,
                           networkType = networkType
  )
  
  if (save_plot && !is.null(prefix)) {
    png(paste0(Dir.f, prefix, "_soft_threshold.png"),
        width = 9, height = 5, units = "in", res = 300)
  }
  
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit, signed R^2",
    type = "n",
    main = "Scale independence"
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  abline(h = 0.80, col = "red")
  
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = "Mean connectivity"
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  
  if (save_plot && !is.null(prefix)) {
    dev.off()
  }
  
  return(sft)
}

plot_gene_tree <- function(geneTree, prefix) {
  filename <- paste0(Dir.f, prefix, "_gene_tree.pdf")
  pdf(file = filename, width = 12, height = 9)
  plot(
    geneTree,
    xlab = "",
    sub = "",
    main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE,
    hang = 0.04
  )
  dev.off()
}

plot_module_tree <- function(geneTree, colors_mat, labels, prefix, main_text) {
  filename <- paste0(Dir.f, prefix, "_module_tree.png")
  png(file = filename, width = 8, height = 6, units = "in", res = 300)
  plotDendroAndColors(
    geneTree,
    colors_mat,
    labels,
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = main_text
  )
  dev.off()
}

plot_me_clustering <- function(MEs, cutHeight, prefix) {
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")
  
  filename <- paste0(Dir.f, prefix, "_eigengene_clustering.png")
  png(file = filename, width = 7, height = 6, units = "in", res = 300)
  plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  abline(h = cutHeight, col = "red")
  dev.off()
  
  invisible(METree)
}

# -----------------------------------------------------------------------------
# 3. GENE & MODULE MEMBERSHIP
# -----------------------------------------------------------------------------

gene_info <- function(traits, variable, MEs, data, mergedColor = mergedColors, Dir = Dir.o) {
  var.df = as.data.frame(get(variable, traits))
  rownames(var.df) <- rownames(traits) 
  names(var.df) = variable
  
  modNames = substring(names(MEs), 3)
  geneModuleMembership = cor(data, MEs, use = "p") 
  nGenes = ncol(data);
  nSamples = nrow(data);
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  table(rownames(data) == rownames(var.df))
  geneTraitSignificance = cor(data, var.df, use = "p");
  GSPvalue = as.data.frame(corPvalueStudent(geneTraitSignificance, nSamples));
  GSPvalue.0 = as.vector(GSPvalue[,variable])
  GSQvalue.0 = qvalue(p=GSPvalue.0, lambda = 0, fdr.level = 0.05)
  GSQvalue <- as.data.frame(GSQvalue.0$qvalues, row.names=names(data))
  summary(GSQvalue.0)
  pi0 <- GSQvalue.0$pi0
  lfdr <- GSQvalue.0$lfdr 
  
  probes=names(data)
  geneInfo0 = data.frame(datExpr = probes, ModuleColor = mergedColors,
                         geneTraitSignificance,                      
                         GSQvalue)
  
  filename = paste(Dir, variable, "_not_ordered.csv", sep = "")
  write.csv(geneInfo0, filename)
  
  var.vec = get(variable, traits)
  NS1=networkScreening(y=var.vec, datME=MEs, datExpr=data,
                       oddPower=3, blockSize=1000, minimumSampleSize=4,
                       addMEy=TRUE, removeDiag=FALSE, weightESy=0.5, corOptions = "use = 'p', method = 'spearman'")
  filename2 = paste(Dir, variable, "_gene_screening.csv", sep = "")
  write.csv(NS1, filename2)
}

gene_info_no_trait <- function(MEs, data, Dir=Dir.o) {
  modNames = substring(names(MEs), 3)
  geneModuleMembership = cor(data, MEs, use = "p")
  
  nGenes = ncol(data);
  nSamples = nrow(data);
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  probes=names(data)
  geneInfo0 = data.frame(datExpr = probes, ModuleColor = mergedColors)
  
  filename = paste(Dir.o, "_not_ordered.csv", sep = "")
  write.csv(geneInfo0, filename)
  return(geneInfo0)
}

gene_info_w_KO_flux <- function(traits, variable, MEs, data, module_colors, Dir = Dir.o, Kegg = kegg) {
  var.df = as.data.frame(traits[,variable]) 
  rownames(var.df) <- rownames(traits) 
  names(var.df) = variable
  
  modNames = substring(names(MEs), 3)
  
  geneModuleMembership = cor(data, MEs, use = "p")
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  MMPvector = as.vector(MMPvalue)
  MMQvalue.0 = qvalue(p=MMPvector, lambda = 0, fdr.level = 0.05)
  MMQvalue <- as.data.frame(MMQvalue.0$qvalues, row.names=names(data))
  names(MMQvalue) = paste("q.MM", modNames, sep="");
  
  geneTraitSignificance = cor(data, var.df, use = "p");
  GSPvalue = as.data.frame(corPvalueStudent(geneTraitSignificance, nSamples));
  GSPvalue.0 = as.vector(GSPvalue[,variable])
  GSQvalue.0 = qvalue(p=GSPvalue.0, lambda = 0, fdr.level = 0.05)
  GSQvalue <- as.data.frame(GSQvalue.0$qvalues, row.names=names(data))
  summary(GSQvalue.0)
  pi0 <- GSQvalue.0$pi0
  lfdr <- GSQvalue.0$lfdr   
  
  probes=names(data)
  geneInfo0 = data.frame(Gene = probes, ModuleColor = module_colors,
                         GS = geneTraitSignificance, GS.Q = GSQvalue)
  
  modOrder = order(abs(cor(MEs, var.df, use = "p")));
  
  for (mod in 1:ncol(geneModuleMembership)) {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]], MMQvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""),                        
                         paste("q.MM.", modNames[modOrder[mod]], sep=""))
  }
  
  geneInfoK <- geneInfo0 %>%
    merge(Kegg, by = "Gene", all.x=TRUE, all.y=FALSE) 
  
  filename = paste(Dir, variable, "_ordered.csv", sep = "")
  write.csv(geneInfoK, filename)
}

compute_gene_membership_table <- function(datExpr, datTraits, colors, MEs, prefix) {
  MEs <- orderMEs(MEs)
  modNames <- substring(names(MEs), 3)
  genes <- colnames(datExpr)
  nSamples <- nrow(datExpr)
  
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  names(geneModuleMembership) <- paste0("MM.", modNames)
  
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(MMPvalue) <- paste0("p.MM.", modNames)
  
  MMQvalue <- MMPvalue
  for (i in seq_len(ncol(MMPvalue))) {
    MMQvalue[[i]] <- safe_qvalues(MMPvalue[[i]])
  }
  names(MMQvalue) <- paste0("q.MM.", modNames)
  
  Rhizo <- as.data.frame(datTraits$Rhizo)
  names(Rhizo) <- "Rhizo"
  
  geneTraitSignificance <- as.data.frame(cor(datExpr, Rhizo, use = "p"))
  names(geneTraitSignificance) <- "GS.Rhizo"
  
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(GSPvalue) <- "p.GS.Rhizo"
  
  GSQvalue <- data.frame(
    q.GS.Rhizo = p.adjust(GSPvalue[[1]], method = "BH"),
    row.names = genes
  )
  
  geneInfo <- data.frame(
    Gene = genes,
    ModuleColor = colors,
    geneTraitSignificance,
    GSPvalue,
    GSQvalue,
    stringsAsFactors = FALSE
  )
  
  modOrder <- order(-abs(cor(MEs, Rhizo, use = "p")))
  
  for (mod in seq_len(ncol(geneModuleMembership))) {
    idx <- modOrder[mod]
    geneInfo[[paste0("MM.", modNames[idx])]]   <- geneModuleMembership[[idx]]
    geneInfo[[paste0("p.MM.", modNames[idx])]] <- MMPvalue[[idx]]
    geneInfo[[paste0("q.MM.", modNames[idx])]] <- MMQvalue[[idx]]
  }
  
  geneInfo$AssignedModule <- geneInfo$ModuleColor
  geneInfo$MM.assigned <- NA_real_
  geneInfo$p.MM.assigned <- NA_real_
  geneInfo$q.MM.assigned <- NA_real_
  
  for (m in modNames) {
    inMod <- geneInfo$AssignedModule == m
    geneInfo$MM.assigned[inMod]   <- geneInfo[inMod, paste0("MM.", m)]
    geneInfo$p.MM.assigned[inMod] <- geneInfo[inMod, paste0("p.MM.", m)]
    geneInfo$q.MM.assigned[inMod] <- geneInfo[inMod, paste0("q.MM.", m)]
  }
  
  grey_idx <- geneInfo$AssignedModule == "grey"
  geneInfo$MM.assigned[grey_idx] <- NA
  geneInfo$p.MM.assigned[grey_idx] <- NA
  geneInfo$q.MM.assigned[grey_idx] <- NA
  
  geneInfo$MM.assigned.significant <-
    geneInfo$q.MM.assigned < moduleSigAlpha
  
  geneInfo$MM.assigned.sig.and.strong <-
    geneInfo$q.MM.assigned < moduleSigAlpha &
    abs(geneInfo$MM.assigned) >= kME_threshold
  
  write.csv(geneInfo,
            file = paste0(Dir.o, prefix, "_geneInfo_with_MM_pq_assigned.csv"), row.names = FALSE)
  write.csv(subset(geneInfo, MM.assigned.significant),
            file = paste0(Dir.o, prefix, "_geneInfo_significant_within_assigned_module.csv"), row.names = FALSE)
  write.csv(subset(geneInfo, MM.assigned.sig.and.strong),
            file = paste0(Dir.o, prefix, "_geneInfo_significant_and_strong_within_assigned_module.csv"), row.names = FALSE)
  
  geneInfo
}

write_module_key <- function(datExpr, colors, prefix) {
  df_key <- data.frame(
    KO = colnames(datExpr),
    moduleColors = colors,
    stringsAsFactors = FALSE
  )
  write.csv(df_key, paste0(Dir.o, prefix, "_Module_key.csv"), row.names = FALSE)
  
  kegg.brite <- read.csv("./output/KEGG/Kegg_brite_all.csv")
  if ("X" %in% colnames(kegg.brite)) kegg.brite$X <- NULL
  
  df_key_pathways <- merge(df_key, kegg.brite, by = "KO", all.x = TRUE)
  readr::write_delim(df_key_pathways,
                     paste0(Dir.o, prefix, "_df_key_pathways.txt"),
                     delim = "\t")
  
  list(df_key = df_key, df_key_pathways = df_key_pathways)
}

write_module_key_family <- function(datExpr, colors, prefix) {
  df_key <- data.frame(
    KO = colnames(datExpr),
    moduleColors = colors,
    stringsAsFactors = FALSE
  )
  
  df_key <- df_key %>%
    separate(KO, c("KO", "Family"), sep = "\\.", remove = FALSE) 
  write.csv(df_key, paste0(Dir.o, prefix, "_Module_key.csv"), row.names = FALSE)
  
  kegg.brite <- read.csv("./output/KEGG/Kegg_brite_all.csv")
  if ("X" %in% colnames(kegg.brite)) kegg.brite$X <- NULL
  
  df_key_pathways <- merge(df_key, kegg.brite, by = "KO", all.x = TRUE)
  readr::write_delim(df_key_pathways,
                     paste0(Dir.o, prefix, "_df_key_pathways.txt"),
                     delim = "\t")
  
  list(df_key = df_key, df_key_pathways = df_key_pathways)
}

# -----------------------------------------------------------------------------
# 4. MODULE-TRAIT RELATIONSHIPS
# -----------------------------------------------------------------------------

plot_sample_tree <- function(datExpr, datTraits, prefix) {
  sampleTree <- flashClust(dist(datExpr), method = "average")
  traitColors <- numbers2colors(datTraits, signed = TRUE, naColor = "grey")
  
  filename <- paste0(Dir.f, prefix, "_sample_dendrogram_trait_heatmap.png")
  png(filename, width = 10, height = 7, units = "in", res = 300)
  plotDendroAndColors(
    sampleTree, traitColors,
    groupLabels = names(datTraits),
    main = "Sample dendrogram and trait heatmap"
  )
  dev.off()
  
  invisible(sampleTree)
}

compute_module_trait_heatmap <- function(MEs, datTraits, prefix,
                                         trait_order = c("Bulk", "Rhizo", "RhizoDet", "Detritus"),
                                         module_order = colnames(orderMEs(MEs))) {
  nSamples <- nrow(datTraits)
  
  datTraits.sub <- datTraits[, trait_order, drop = FALSE]
  MEs.sub <- MEs[, module_order, drop = FALSE] 
  datTraits.sub <- datTraits.sub[match(rownames(MEs.sub), rownames(datTraits.sub)), , drop = FALSE]
  
  moduleTraitCor <- t(cor(MEs.sub, datTraits.sub, use = "p"))
  moduleTraitPvalue <- corPvalueStudent(as.matrix(moduleTraitCor), nSamples)
  moduleTraitQvalue <- matrix(
    p.adjust(as.vector(moduleTraitPvalue), method = "fdr"),
    nrow = nrow(moduleTraitPvalue),
    ncol = ncol(moduleTraitPvalue),
    byrow = FALSE,
    dimnames = dimnames(moduleTraitPvalue)
  )
  
  write.csv(moduleTraitCor, paste0(Dir.o, prefix, "_module_trait_cor.csv"))
  write.csv(moduleTraitPvalue, paste0(Dir.o, prefix, "_module_trait_pvalue.csv"))
  write.csv(moduleTraitQvalue, paste0(Dir.o, prefix, "_module_trait_qvalue.csv"))
  
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitQvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  filename <- paste0(Dir.f, prefix, "_module_trait_heatmap.png")
  png(file = filename, width = 9, height = 7, units = "in", res = 300)
  par(mar = c(10, 12, 5, 5))
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = colnames(moduleTraitCor),
    yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.7,
    cex.lab.y = 1,
    cex.lab.x = 1,
    zlim = c(-1, 1),
    main = "Module-trait relationships"
  )
  dev.off()
  
  list(cor = moduleTraitCor, p = moduleTraitPvalue, q = moduleTraitQvalue)
}

plot_module_trait_heatmap_subset <- function(full_results, 
                                             trait_subset = NULL,
                                             module_subset = NULL,
                                             prefix = "subset") {
  # Pull from already-computed full results
  moduleTraitCor    <- full_results$cor
  moduleTraitQvalue <- full_results$q
  
  # Subset traits (rows) if specified
  if (!is.null(trait_subset)) {
    moduleTraitCor    <- moduleTraitCor[trait_subset, , drop = FALSE]
    moduleTraitQvalue <- moduleTraitQvalue[trait_subset, , drop = FALSE]
  }
  
  # Subset modules (columns) if specified
  if (!is.null(module_subset)) {
    moduleTraitCor    <- moduleTraitCor[, module_subset, drop = FALSE]
    moduleTraitQvalue <- moduleTraitQvalue[, module_subset, drop = FALSE]
  }
  
  # Build text matrix
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitQvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  # Plot
  filename <- paste0(Dir.f, prefix, "_module_trait_heatmap_subset.png")
  png(file = filename, width = 9, height = 7, units = "in", res = 300)
  par(mar = c(10, 12, 5, 5))
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = colnames(moduleTraitCor),
    yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.7,
    cex.lab.y = 1,
    cex.lab.x = 1,
    zlim = c(-1, 1),
    main = "Module-trait relationships"
  )
  dev.off()
  
  invisible(list(cor = moduleTraitCor, q = moduleTraitQvalue))
}

plot_module_trait_with_gene_bars <- function(full_results,
                                             df_key,
                                             trait_subset = NULL,
                                             module_subset = NULL,
                                             prefix = "combined") {
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  
  # -------------------------
  # 1. Subset correlation results (same as subset function)
  # -------------------------
  moduleTraitCor    <- full_results$cor
  moduleTraitQvalue <- full_results$q
  
  # Strip ME prefix from column names
  # Strip ME prefix from column names
  colnames(moduleTraitCor)    <- str_remove(colnames(moduleTraitCor), "^ME")
  colnames(moduleTraitQvalue) <- str_remove(colnames(moduleTraitQvalue), "^ME")
  
  # Also strip ME from module_subset if present, and save as the ordering reference
  if (!is.null(module_subset)) {
    moduleTraitCor    <- moduleTraitCor[, module_subset, drop = FALSE]
    moduleTraitQvalue <- moduleTraitQvalue[, module_subset, drop = FALSE]
  }
  
  if (!is.null(trait_subset)) {
    moduleTraitCor    <- moduleTraitCor[trait_subset, ,drop = FALSE]
    moduleTraitQvalue <- moduleTraitQvalue[trait_subset, ,drop = FALSE]
  }
  
  # Clean ordered vector to use as factor levels everywhere
  module_order <- str_remove(colnames(moduleTraitCor), "^ME")
  
  # Use module_order to set x-axis order, fall back to column order
  if (is.null(module_subset)) {
    module_subset <- colnames(moduleTraitCor)
  }
  
  # -------------------------
  # 2. Heatmap as ggplot 
  # (converts labeledHeatmap logic to ggplot so patchwork can align it)
  # -------------------------
  cor_df <- as.data.frame(moduleTraitCor) %>%
    rownames_to_column(var = "Trait") %>%
    pivot_longer(-Trait, names_to = "Module", values_to = "Correlation")
  
  q_df <- as.data.frame(moduleTraitQvalue) %>%
    rownames_to_column(var = "Trait") %>%
    pivot_longer(-Trait, names_to = "Module", values_to = "Qvalue")
  
  heatmap_df <- left_join(cor_df, q_df, by = c("Trait", "Module")) %>%
    mutate(
      Module = factor(Module, levels = module_subset),
      Trait = factor(Trait, levels = trait_subset),
      label  = paste0(signif(Correlation, 2), "\n(", signif(Qvalue, 1), ")")
    ) %>%
    mutate(Module = str_remove(Module, "ME"))
  
  p_heat <- ggplot(heatmap_df, aes(x = Module, y = Trait, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 2.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Correlation") +
    scale_x_discrete(limits = module_subset) +  # force order here
    theme_bw() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 12),
      panel.grid   = element_blank(),
      strip.background = element_rect(fill = "grey90")
    ) +
    labs(y = "Trait")
  # -------------------------
  # 3. Load exudate & litter gene lists
  # -------------------------
  exudate.list <- read.csv("./data/root_exudate_KOs_5_21_26_Claude_v2.csv", header = TRUE) %>%
    # filter(
    #   Role != "Amino acid catabolism",
    #   Role != "Amino acid metabolism",
    #   Role != "Amino acid uptake",
    #   Role != "Compatible solutes",
    #   Role != "Glyoxylate cycle",
    #   Role != "Ring cleavage",
    #   Role != "Small organics",
    #   Role != "Organic acid catabolism",
    #   Role != "Organic acid uptake"
    # )
       filter(
         Root_Exudate_Class != "Xylose",
         Root_Exudate_Class != "Alcohols",
         Root_Exudate_Class != "Trehalose",
         Root_Exudate_Class != "Amino acids",
         Root_Exudate_Class != "Galactose",
         Root_Exudate_Class != "Arabinose",
         Root_Exudate_Class != "Gluconate",
         Root_Exudate_Class != "GABA",
         Root_Exudate_Class != "Aromatic",
         Root_Exudate_Class != "Flavonoids",
         Role != "Acetate catabolism"
       )
  
  litter.list <- read.csv("./data/litter_degradation_KOs_5_21_26_Claude.csv", header = TRUE) %>%
    filter(
      Litter_Component != "Chitin",
      Litter_Component != "Arabinogalactan"
    )
  
  # -------------------------
  # 4. Merge with module assignments
  # -------------------------
  
  exudate.modules <- merge(df_key, exudate.list, by = "KO") %>%
    mutate(moduleColors = factor(moduleColors, levels = module_subset))
  
  litter.modules <- merge(df_key, litter.list, by = "KO") %>%
    mutate(moduleColors = factor(moduleColors, levels = module_subset))
  
  # -------------------------
  # 5. Count genes per module, with ALL modules present even if zero
  # This is the key step - expand to full module list first
  # -------------------------
  all_modules <- data.frame(moduleColors = factor(module_subset, levels = module_subset))
  
  exudate_counts <- exudate.modules %>%
    group_by(moduleColors, Role) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(moduleColors = factor(module_subset, levels = module_subset),
             Role,
             fill = list(n = 0)) %>%
    drop_na() %>%
    filter(n >0) %>%
    mutate(moduleColors = factor(moduleColors, levels = module_subset))  # re-level after complete
  
  litter_counts <- litter.modules %>%
    group_by(moduleColors, Litter_Component) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(moduleColors = factor(module_subset, levels = module_subset), 
             Litter_Component, 
             fill = list(n = 0))%>%
    drop_na() %>%
    filter(n >0) %>%
    mutate(moduleColors = factor(moduleColors, levels = module_subset))  # re-level after complete
  
  # -------------------------
  # 6. Color palettes
  # -------------------------
  gray_palette_e <- gray.colors(length(unique(exudate_counts$Role)), start = 0.8, end = 0.25)
  gray_palette_l <- gray.colors(length(unique(litter_counts$Litter_Component)), start = 0.8, end = 0.25)
  
  # -------------------------
  # 7. Exudate bar plot — x axis hidden, aligned to heatmap
  # -------------------------
  p_exudate <- ggplot(exudate_counts, aes(x = moduleColors, y = n, fill = Role)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = gray_palette_e) +
    scale_x_discrete(drop = FALSE) +
    theme_bw() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid   = element_blank(),
      legend.key.width = unit(0.5, "cm"),
      legend.text = element_text(size = 12)
    ) +
    labs(y = "Exudate\ngenes", fill = "Role")
  
  # -------------------------
  # 8. Litter bar plot — x axis shown for module labels
  # -------------------------
  p_litter <- ggplot(litter_counts, aes(x = moduleColors, y = n, fill = Litter_Component)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = gray_palette_l) +
    scale_x_discrete(drop = FALSE) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
      panel.grid   = element_blank(),
      legend.key.width = unit(0.5, "cm"),
      legend.text = element_text(size = 12)
    ) +
    labs(x = "Module", y = "Litter\ngenes", fill = "Litter substrate")
  
  # -------------------------
  # 9. Stitch together with patchwork
  # -------------------------
  combined <- p_heat / p_exudate / p_litter +
    plot_layout(heights = c(3, 2, 2)) # heatmap taller than bars
  
  # -------------------------
  # 10. Save
  # -------------------------
  filename <- paste0(Dir.f, prefix, "_heatmap_with_gene_bars.pdf")
  ggsave(filename, combined,
         width = max(8, length(module_subset) * 0.6),
         height = 10,
         units = "in",
         dpi = 300)
  
  invisible(combined)
}

# -----------------------------------------------------------------------------
# 5. STATISTICAL TESTING — EIGENGENE ~ ZONE/TIMEPOINT
# -----------------------------------------------------------------------------

run_eigengene_zone_stats <- function(MEs, datTraits, prefix) {
  library(emmeans)
  
  MEs <- as.data.frame(MEs)
  
  trait_df <- datTraits %>%
    rownames_to_column(var = "SampleID") %>%
    separate_wider_delim(SampleID,
                         delim = ".",
                         names = c("ID", "Treatment", "Zone", "Timepoint"),
                         cols_remove = FALSE) %>%
    mutate(
      Zone = factor(Zone),
      Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks"))
    )
  
  me_long <- MEs %>%
    rownames_to_column(var = "SampleID") %>%
    pivot_longer(cols = starts_with("ME"),
                 names_to = "Module",
                 values_to = "Eigengene")
  
  stat_df <- me_long %>%
    left_join(trait_df, by = "SampleID")
  
  modules <- unique(stat_df$Module)
  anova_out <- list()
  posthoc_out <- list()
  
  for (m in modules) {
    dfm <- stat_df %>% filter(Module == m)
    fit <- lm(Eigengene ~ Zone, data = dfm)
    
    aov_tab <- anova(fit)
    aov_df <- data.frame(
      Module = m,
      Term = rownames(aov_tab),
      DF = aov_tab$Df,
      Fvalue = aov_tab$`F value`,
      Pvalue = aov_tab$`Pr(>F)`,
      row.names = NULL
    ) %>%
      filter(Term %in% c("Zone"))
    anova_out[[m]] <- aov_df
    
    emm <- emmeans(fit, ~ Zone)
    contrasts <- pairs(emm, adjust = "BH")
    contrasts_df <- as.data.frame(contrasts) %>%
      mutate(Module = m)
    posthoc_out[[m]] <- contrasts_df
  }
  
  anova_tbl <- bind_rows(anova_out)
  anova_tbl$Qvalue <- p.adjust(anova_tbl$Pvalue, method = "BH")
  posthoc_tbl <- bind_rows(posthoc_out)
  
  write.csv(anova_tbl, paste0(Dir.o, prefix, "_eigengene_zone_anova.csv"), row.names = FALSE)
  write.csv(posthoc_tbl, paste0(Dir.o, prefix, "_eigengene_zone_within_posthoc.csv"), row.names = FALSE)
  
  return(list(anova = anova_tbl, posthoc = posthoc_tbl))
}

run_eigengene_zone_time_stats <- function(MEs, datTraits, prefix) {
  library(emmeans)
  
  MEs <- as.data.frame(MEs)
  
  trait_df <- datTraits %>%
    rownames_to_column(var = "SampleID") %>%
    separate_wider_delim(SampleID,
                         delim = ".",
                         names = c("ID", "Treatment", "Zone", "Timepoint"),
                         cols_remove = FALSE) %>%
    mutate(
      Zone = factor(Zone),
      Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks"))
    )
  
  me_long <- MEs %>%
    rownames_to_column(var = "SampleID") %>%
    pivot_longer(cols = starts_with("ME"),
                 names_to = "Module",
                 values_to = "Eigengene")
  
  stat_df <- me_long %>%
    left_join(trait_df, by = "SampleID")
  
  modules <- unique(stat_df$Module)
  anova_out <- list()
  posthoc_out <- list()
  
  for (m in modules) {
    dfm <- stat_df %>% filter(Module == m)
    fit <- lm(Eigengene ~ Zone * Timepoint, data = dfm)
    
    aov_tab <- anova(fit)
    aov_df <- data.frame(
      Module = m,
      Term = rownames(aov_tab),
      DF = aov_tab$Df,
      Fvalue = aov_tab$`F value`,
      Pvalue = aov_tab$`Pr(>F)`,
      row.names = NULL
    ) %>%
      filter(Term %in% c("Zone", "Timepoint", "Zone:Timepoint"))
    anova_out[[m]] <- aov_df
    
    emm <- emmeans(fit, ~ Zone | Timepoint)
    contrasts <- pairs(emm, adjust = "BH")
    contrasts_df <- as.data.frame(contrasts) %>%
      mutate(Module = m)
    posthoc_out[[m]] <- contrasts_df
  }
  
  anova_tbl <- bind_rows(anova_out)
  anova_tbl$Qvalue <- p.adjust(anova_tbl$Pvalue, method = "BH")
  posthoc_tbl <- bind_rows(posthoc_out)
  
  write.csv(anova_tbl, paste0(Dir.o, prefix, "_eigengene_zone_time_anova.csv"), row.names = FALSE)
  write.csv(posthoc_tbl, paste0(Dir.o, prefix, "_eigengene_zone_within_time_posthoc.csv"), row.names = FALSE)
  
  return(list(anova = anova_tbl, posthoc = posthoc_tbl))
}

compare_modules_within_treatment <- function(MEs, datTraits, treatment_name, prefix) {
  library(emmeans)
  
  trait_df <- datTraits %>%
    rownames_to_column(var = "SampleID") %>%
    separate_wider_delim(SampleID,
                         delim = ".",
                         names = c("ID", "Treatment", "Zone", "Timepoint"),
                         cols_remove = FALSE)
  
  me_long <- MEs %>%
    as.data.frame() %>%
    rownames_to_column(var = "SampleID") %>%
    pivot_longer(cols = starts_with("ME"),
                 names_to = "Module",
                 values_to = "Eigengene")
  
  df <- me_long %>%
    left_join(trait_df, by = "SampleID") 
  
  fit <- lm(Eigengene ~ Module * Zone, data = df)
  
  emm <- emmeans(fit, ~ Module | Zone)
  contrasts <- pairs(emm, adjust = "BH")
  res <- as.data.frame(contrasts)
  write.csv(res, paste0(Dir.o, prefix, "_module_comparison_within_", treatment_name, ".csv"), row.names = FALSE)
  
  emm <- emmeans(fit, ~ Zone | Module)
  contrasts <- pairs(emm, adjust = "BH")
  res <- as.data.frame(contrasts)
  write.csv(res, paste0(Dir.o, prefix, "_module_comparison_within_", treatment_name, "_zone.csv"), row.names = FALSE)
  
  return(res)
}

# -----------------------------------------------------------------------------
# 6. EIGENGENE VISUALIZATION
# -----------------------------------------------------------------------------

plot_me_sample_heatmap <- function(MEs, datTraits, prefix) {
  col_ann <- make_sample_annotation(datTraits)
  
  MEs_use <- as.data.frame(MEs)
  if ("MEgrey" %in% colnames(MEs_use)) {
    MEs_use <- dplyr::select(MEs_use, -MEgrey)
  }
  
  MEs_use <- MEs_use[order(match(rownames(MEs_use), rownames(col_ann))), , drop = FALSE]
  col_ann <- col_ann[match(rownames(MEs_use), rownames(col_ann)), , drop = FALSE]
  
  filename <- paste0(Dir.f, prefix, "_eigengene_heatmap.pdf")
  pdf(file = filename, height = 8, width = 6)
  pheatmap(
    MEs_use,
    cluster_col = TRUE,
    cluster_row = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    fontsize = 6,
    annotation_row = col_ann,
    annotation_colors = ann_color
  )
  dev.off()
}

plot_modules <- function(module_eigengenes, module_color) {
  ME = module_eigengenes[, paste("ME",module_color, sep="")]
  ME <- as.data.frame(ME)
  ME$SampleID <- rownames(module_eigengenes)
  ME <- ME %>%
    mutate(Treat = case_when(
      grepl("Drought",SampleID) ~ "Drought",
      grepl( "Untrt",SampleID) ~ "Control")) %>%
    mutate(Time = case_when(
      grepl("4weeks", SampleID) ~ "4weeks",
      grepl("8weeks", SampleID) ~ "8weeks",
      grepl("12weeks", SampleID) ~ "12weeks"
    ))  
  ME$Time <- factor(ME$Time, levels = c("4weeks", "8weeks", "12weeks"))
  
  barplot_time <- ggplot(ME, aes(x=Treat, y= ME)) +
    geom_point() + geom_boxplot(fill = module_color) +
    facet_wrap(~Time) +
    labs(title = paste(module_color, "module", sep = "_"), y = "Eigengene Expression", x= "Treatment") +
    theme_bw() +
    theme(text = element_text(size = size, color = "black"),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
          axis.text.y = element_text(color = "black"))
  print(barplot_time)
}

plot_modules_drought_no_time <- function(module_eigengenes, metadata, module_color) {
  metadata <- metadata %>%
    rownames_to_column(var = "SampleID") %>%
    separate_wider_delim(SampleID,
                         delim = ".",
                         names = c("ID", "Treatment", "Zone", "Timepoint"),
                         cols_remove = FALSE)
  
  ME = module_eigengenes[, paste("ME",module_color, sep="")]
  ME <- as.data.frame(ME)
  ME$SampleID <- rownames(module_eigengenes)
  ME <- ME %>%
    merge(metadata, by = "SampleID") %>%
    mutate(Zone = factor(Zone, levels = c("Bulk", "Rhizo", "RhizoDetritus", "Detritus")))
  
  barplot_time <- ggplot(ME, aes(x=Zone, y= ME)) +
    geom_point() + geom_boxplot(fill = module_color) +
    labs(title = paste(module_color, "module", sep = "_"), y = "Eigengene Expression", x= "Treatment") +
    theme_bw() +
    theme(text = element_text(size = size, color = "black"),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"))
  print(barplot_time)
}

plot_modules_by_zone <- function(module_eigengenes, metadata, module_color, ncol=2) {
  metadata <- metadata %>%
    rownames_to_column(var = "SampleID") %>%
    separate_wider_delim(SampleID,
                         delim = ".",
                         names = c("ID", "Treatment", "Zone", "Timepoint"),
                         cols_remove = FALSE)
  
  ME = module_eigengenes[, paste("ME",module_color, sep="")]
  ME <- as.data.frame(ME)
  ME$SampleID <- rownames(module_eigengenes)
  ME <- ME %>%
    merge(metadata, by = "SampleID") %>%
    mutate(Zone = factor(Zone, levels = c("Bulk", "Rhizo", "RhizoDetritus", "Detritus")))
  
  ME$Timepoint <- factor(ME$Timepoint, levels = c("4weeks", "8weeks", "12weeks"))
  
  barplot_time <- ggplot(ME, aes(x=Timepoint, y= ME)) +
    geom_point() + geom_boxplot(fill = module_color) +
    facet_wrap(~Zone, ncol = ncol) +
    labs(title = paste(module_color, "module", sep = "_"), y = "Eigengene Expression", x= "Treatment") +
    theme_bw() +
    theme(text = element_text(size = size, color = "black"),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"))
  print(barplot_time)
}

plot_modules_by_time <- function(module_eigengenes, metadata, module_color, ncol=2) {
  metadata <- metadata %>%
    rownames_to_column(var = "SampleID") %>%
    separate_wider_delim(SampleID,
                         delim = ".",
                         names = c("ID", "Treatment", "Zone", "Timepoint"),
                         cols_remove = FALSE)
  
  ME = module_eigengenes[, paste("ME",module_color, sep="")]
  ME <- as.data.frame(ME)
  ME$SampleID <- rownames(module_eigengenes)
  ME <- ME %>%
    merge(metadata, by = "SampleID") %>%
    mutate(Zone = factor(Zone, levels = c("Bulk", "Rhizo", "RhizoDetritus", "Detritus")))
  
  ME$Timepoint <- factor(ME$Timepoint, levels = c("4weeks", "8weeks", "12weeks"))
  
  barplot_time <- ggplot(ME, aes(x=Zone, y= ME)) +
    geom_point() + geom_boxplot(fill = module_color) +
    facet_wrap(~Timepoint, ncol = ncol) +
    labs(title = paste(module_color, "module", sep = "_"), y = "Eigengene Expression", x= "Treatment") +
    theme_bw() +
    theme(text = element_text(size = size, color = "black"),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"))
  print(barplot_time)
}

plot_eigengenes_gg <- function(MEs, datTraits,
                               facet_by = c("Timepoint"),
                               color_by = "Zone",
                               shape_by = NULL,
                               prefix = "plot",
                               save = TRUE,
                               order_MEs=order_MEs_colors,
                               ncol=1,
                               w = w,
                               h = h) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  trait_df <- datTraits %>%
    rownames_to_column(var = "SampleID") %>%
    separate_wider_delim(SampleID,
                         delim = ".",
                         names=c("ID", "Treatment", "Zone", "Timepoint"),
                         cols_remove = FALSE) %>%
    mutate(
      Zone = factor(Zone),
      Timepoint = factor(Timepoint, levels = c("4weeks", "8weeks", "12weeks")),
      Treatment = factor(Treatment)
    )
  
  me_long <- as.data.frame(MEs) %>%
    rownames_to_column(var = "SampleID") %>%
    pivot_longer(cols = starts_with("ME"),
                 names_to = "Module",
                 values_to = "Eigengene") %>%
    filter(Module != "MEgrey") %>%
    mutate(Module = substring(Module, 3)) %>%
    mutate(Module = factor(Module, levels = c(order_MEs)))
  
  df <- me_long %>%
    left_join(trait_df, by = "SampleID") %>%
    mutate(Zone = factor(Zone, c("Bulk", "Rhizo", "RhizoDetritus", "Detritus")))
  
  p <- ggplot(df, aes(x = Module, y = Eigengene, color = Module)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7,
                aes(shape = if (!is.null(shape_by)) .data[[shape_by]] else NULL)) +
    scale_color_identity(guide = "legend") +
    geom_hline(yintercept=0, linetype = "dashed", color = "grey50")
  
  p2 <- ggplot(df, aes(x = Zone, y = Eigengene, color = Module)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7,
                aes(shape = if (!is.null(shape_by)) .data[[shape_by]] else NULL)) +
    scale_color_identity(guide = "legend") +
    geom_hline(yintercept=0, linetype = "dashed", color = "grey50")
  
  if (length(facet_by) == 1 && facet_by == "Module") {
    p <- p2 + facet_wrap(as.formula(paste("~", facet_by)), ncol = ncol) +
      labs(x = "Zone",y = "Eigengene expression",color = color_by)
  } else if (length(facet_by) == 1) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)), ncol = ncol) +
      labs(x = "Module",y = "Eigengene expression",color = color_by)
  } else if (length(facet_by) == 2) {
    p <- p + facet_wrap(as.formula(paste(facet_by[1], "~", facet_by[2])), ncol = ncol) +
      labs(x = "Module",y = "Eigengene expression",color = color_by)
  } else {
    warning("facet_by must be 'Module', a single variable, or two variables")
  }
  
  p <- p +
    ylim(c(-0.4, 0.4)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid   = element_blank(),
      strip.background = element_rect(fill = "grey90")
    ) 
  
  if (save) {
    ggsave(
      filename = paste0(Dir.f, prefix, "_eigengene_plot.pdf"),
      plot = p,
      width = w,
      height = h,
      units = "in",
      dpi = 300
    )
  }
  
  return(p)
}

plot_eigengenes_by_module <- function(MEs, datTraits, prefix, ncol=4) {
  modules <- substring(colnames(MEs), 3)
  
  for (m in modules) {
    try({
      plot_modules_by_time(MEs, datTraits, m, ncol=ncol)
      ggsave(
        filename = paste0(Dir.f, prefix, "_", m, "_EGexp_zone_time.png"),
        width = 5, height = 3, units = "in", dpi = 300
      )
    }, silent = TRUE)
  }
}

# -----------------------------------------------------------------------------
# 7. PATHWAY ENRICHMENT — KEGG DOTPLOTS & NETWORKS
# -----------------------------------------------------------------------------

plot_dp_k <- function(data, module_color, categories=40, Dir=Dir.f) {
  data.mod <- data %>%
    filter(ModuleColor == module_color)
  data.ko <- data.mod$KO
  k <- enrichKEGG(gene = data.ko, organism = 'ko', pvalueCutoff = 0.05)
  dp <- enrichplot::dotplot(k, showCategory = categories)
  print(dp)
}

plot_dp_m <- function(data, module_color, categories=40, Dir=Dir.f) {
  data.mod <- data %>%
    filter(ModuleColor == module_color)
  data.ko <- data.mod$KO
  m <- enrichMKEGG(gene = data.ko, organism = 'ko', pvalueCutoff = 0.2, qvalueCutoff = 0.2)
  dp <- enrichplot::dotplot(m, showCategory = categories)
  print(dp)
}

plot_net_k <- function(data, module_color, categories=40, Dir=Dir.f) {
  data.mod <- data %>%
    filter(ModuleColor == module_color) 
  data.ko <- data.mod$KO
  k <- enrichKEGG(gene = data.ko, organism = 'ko', pvalueCutoff = 0.05)
  
  if (is.null(k) || nrow(k@result) == 0 || 
      all(k@result$geneID == "") ||
      length(k@geneSets) == 0) {
    message("Skipping cnetplot for: ", module_color)
    return(invisible(NULL))
  }
  
  nw <- tryCatch({
    cnetplot(k, showCategory = categories) +
      scale_color_manual(values = c(module_color, module_color)) +
      guides(color = FALSE)
  }, error = function(e) {
    message("cnetplot failed for ", module_color, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (!is.null(nw)) print(nw)
  return(invisible(nw))
}

plot_net_m <- function(data, module_color, categories=40, Dir=Dir.f) {
  data.mod <- data %>%
    filter(ModuleColor == module_color) %>%
    drop_na()
  data.ko <- data.mod$KO
  m <- enrichMKEGG(gene = data.ko, organism = 'ko', pvalueCutoff = 1, qvalueCutoff = 1)
  
  if (is.null(m) || nrow(m@result) == 0 || 
      all(m@result$geneID == "") ||
      length(m@geneSets) == 0) {
    message("Skipping cnetplot for: ", module_color)
    return(invisible(NULL))
  }
  
  nw <- tryCatch({
    cnetplot(m, showCategory = categories) +
      scale_color_manual(values = c(module_color, module_color)) +
      guides(color = FALSE)
  }, error = function(e) {
    message("cnetplot failed for ", module_color, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (!is.null(nw)) print(nw)
  return(invisible(nw))
}

run_pathway_plots <- function(gene_table, prefix) {
  data.m.df <- gene_table %>%
    dplyr::select(Gene, AssignedModule) %>%
    dplyr::rename(KO = Gene, ModuleColor = AssignedModule) %>%
    dplyr::filter(ModuleColor != "grey") %>%
    as.data.frame()
  
  module.counts <- table(data.m.df$ModuleColor)
  keep.modules <- names(module.counts[module.counts >= minGenesFiltered])
  data.m.df <- data.m.df %>% dplyr::filter(ModuleColor %in% keep.modules)
  modules <- unique(data.m.df$ModuleColor)
  
  for (m in modules) {
    plot_dp_k(data.m.df, m)
    ggsave(paste0(Dir.f, prefix, "_", m, "_dotplot_kegg.png"),
           width = w, height = h, units = "in", dpi = 300)
    
    nw_k <- plot_net_k(data.m.df, m)
    if (!is.null(nw_k)) {
      ggsave(paste0(Dir.f, prefix, "_", m, "_network_kegg.png"),
             width = w, height = h, units = "in", dpi = 300)
    }
    
    plot_dp_m(data.m.df, m)
    ggsave(paste0(Dir.f, prefix, "_", m, "_dotplot_module.png"),
           width = w, height = h, units = "in", dpi = 300)
    
    nw_m <- plot_net_m(data.m.df, m)
    if (!is.null(nw_m)) {
      ggsave(paste0(Dir.f, prefix, "_", m, "_network_module.png"),
             width = w, height = h, units = "in", dpi = 300)
    }
  }
}

# -----------------------------------------------------------------------------
# 8. BIOLOGICAL INTERPRETATION — EXUDATE/LITTER & CYTOSCAPE EXPORT
# -----------------------------------------------------------------------------

plot_exudate_litter_summary <- function(df_key, prefix) {
  exudate.list <- read.csv("./data/root_exudate_KOs_5_21_26_Claude.csv", header = TRUE) %>%
    filter(
      Role != "Amino acid catabolism",
      Role != "Amino acid metabolism",
      Role != "Amino acid uptake",
      Role != "Compatible solutes",
      Role != "Glyoxylate cycle",
      Role != "Ring cleavage",
      Role != "Small organics",
      Role != "Organic acid catabolism",
      Role != "Organic acid uptake"
    )
  
  litter.list <- read.csv("./data/litter_degradation_KOs_5_21_26_Claude.csv", header = TRUE) %>%
    filter(Litter_Component != "Pectin",
           Litter_Component != "Arabinogalactan",
           Litter_Component != "Mannan")
  
  exudate.modules <- merge(df_key, exudate.list, by = "KO")
  litter.modules  <- merge(df_key, litter.list, by = "KO")
  
  module.order <- substring(colnames(mergedMEs), 3)
  module.order <- module.order[module.order != "grey"]
  
  exudate.modules$moduleColors <- factor(exudate.modules$moduleColors, levels = module.order)
  litter.modules$moduleColors  <- factor(litter.modules$moduleColors, levels = module.order)
  
  gray_palette_e <- gray.colors(length(unique(exudate.modules$Role)), start = 0.8, end = 0.25)
  gray_palette_l <- gray.colors(length(unique(litter.modules$Litter_Component)), start = 0.8, end = 0.25)
  
  plot.exudate <- exudate.modules %>%
    ggplot(aes(x = moduleColors, fill = Role)) +
    geom_bar() +
    scale_fill_manual(values = gray_palette_e) +
    theme_bw() +
    theme(panel.grid = element_blank())
  ggsave(paste0(Dir.f, prefix, "_exudate_genes.png"),
         plot = plot.exudate, dpi = 300, height = 4, width = 5, units = "in")
  
  plot.litter <- litter.modules %>%
    ggplot(aes(x = moduleColors, fill = Litter_Component)) +
    geom_bar() +
    scale_fill_manual(values = gray_palette_l) +
    theme_bw() +
    theme(panel.grid = element_blank())
  ggsave(paste0(Dir.f, prefix, "_litter_genes.png"),
         plot = plot.litter, dpi = 300, height = 4, width = 5, units = "in")
  
  write.csv(exudate.modules, paste0(Dir.o, prefix, "_exudate_modules.csv"), row.names = FALSE)
  write.csv(litter.modules,  paste0(Dir.o, prefix, "_litter_modules.csv"),  row.names = FALSE)
}

plot_exudate_litter_summary_family <- function(df_key, prefix) {
  exudate.list <- read.csv("./data/root_exudate_consumption_KOs_UPDATED_LH.csv", header = TRUE) %>%
    filter(
      Role != "Amino acid catabolism",
      Role != "Amino acid metabolism",
      Role != "Amino acid uptake",
      Role != "Compatible solutes",
      Role != "Glyoxylate cycle",
      Role != "Ring cleavage",
      Role != "Small organics",
      Role != "Organic acid catabolism",
      Role != "Organic acid uptake"
    )
  
  litter.list <- read.csv("./data/litter_degradation_KOs_FIXED_LH.csv", header = TRUE) %>%
    filter(Litter_Component != "Pectin",
           Litter_Component != "Arabinogalactan",
           Litter_Component != "Mannan")
  
  exudate.modules <- merge(df_key, exudate.list, by = "KO")
  litter.modules  <- merge(df_key, litter.list, by = "KO")
  
  module.order <- substring(colnames(mergedMEs), 3)
  module.order <- module.order[module.order != "grey"]
  
  exudate.modules$moduleColors <- factor(exudate.modules$moduleColors, levels = module.order)
  litter.modules$moduleColors  <- factor(litter.modules$moduleColors, levels = module.order)
  
  gray_palette_e <- gray.colors(length(unique(exudate.modules$Role)), start = 0.8, end = 0.25)
  gray_palette_l <- gray.colors(length(unique(litter.modules$Litter_Component)), start = 0.8, end = 0.25)
  
  gray_palette_ef <- gray.colors(length(unique(exudate.modules$Family)), start = 0.8, end = 0.25)
  gray_palette_lf <- gray.colors(length(unique(litter.modules$Family)), start = 0.8, end = 0.25)
  
  plot.exudate <- exudate.modules %>%
    ggplot(aes(x = moduleColors, fill = Role)) +
    geom_bar() +
    scale_fill_manual(values = gray_palette_e) +
    theme_bw() +
    theme(panel.grid = element_blank())
  ggsave(paste0(Dir.f, prefix, "_exudate_genes.png"),
         plot = plot.exudate, dpi = 300, height = 4, width = 5, units = "in")
  
  plot.exudate.family <- exudate.modules %>%
    ggplot(aes(x = moduleColors, fill = Family)) +
    geom_bar() +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(panel.grid = element_blank())
  ggsave(paste0(Dir.f, prefix, "_exudate_genes_family.png"),
         plot = plot.exudate.family, dpi = 300, height = 4, width = 10, units = "in")
  
  plot.litter <- litter.modules %>%
    ggplot(aes(x = moduleColors, fill = Litter_Component)) +
    geom_bar() +
    scale_fill_manual(values = gray_palette_l) +
    theme_bw() +
    theme(panel.grid = element_blank())
  ggsave(paste0(Dir.f, prefix, "_litter_genes.png"),
         plot = plot.litter, dpi = 300, height = 4, width = 5, units = "in")
  
  plot.litter.family <- litter.modules %>%
    ggplot(aes(x = moduleColors, fill = Family)) +
    geom_bar() +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(panel.grid = element_blank())
  ggsave(paste0(Dir.f, prefix, "_litter_genes_family.png"),
         plot = plot.litter.family, dpi = 300, height = 4, width = 10, units = "in")
  
  write.csv(exudate.modules, paste0(Dir.o, prefix, "_exudate_modules.csv"), row.names = FALSE)
  write.csv(litter.modules,  paste0(Dir.o, prefix, "_litter_modules.csv"),  row.names = FALSE)
}

export_cytoscape_filtered <- function(TOM, datExpr, colors, prefix) {
  probes <- colnames(datExpr)
  colors <- colors[probes]
  
  keepGenes <- probes[colors != "grey"]
  inModule <- probes %in% keepGenes
  modProbes <- probes[inModule]
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modProbes, modProbes)
  modColors <- colors[inModule]
  
  cyt <- exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste0(Dir.o, prefix, "_CytoscapeInput_edges.txt"),
    nodeFile = paste0(Dir.o, prefix, "_CytoscapeInput_nodes.txt"),
    weighted = TRUE,
    threshold = cytoscapeThresh,
    nodeNames = modProbes,
    nodeAttr = modColors
  )
  
  edge.0 <- as.data.frame(cyt$edgeData)
  node.0 <- as.data.frame(cyt$nodeData)
  
  edge <- edge.0 %>%
    dplyr::filter(weight > cytoscapeEdgeMin) %>%
    dplyr::select(-c(fromAltName, toAltName))
  
  node.1 <- node.0 %>%
    merge(., edge, by.x = "nodeName", by.y = "fromNode") %>%
    dplyr::select(-toNode)
  
  node.2 <- node.0 %>%
    merge(., edge, by.x = "nodeName", by.y = "toNode") %>%
    dplyr::select(-fromNode)
  
  node <- rbind(node.1, node.2) %>%
    dplyr::distinct(nodeName, .keep_all = TRUE) %>%
    dplyr::select(-c(weight, direction))
  
  write.table(edge, file = paste0(Dir.o, prefix, "_CytoscapeInput_edges_filt.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(node, file = paste0(Dir.o, prefix, "_CytoscapeInput_nodes_filt.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}