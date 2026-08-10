#functions for lc-ms analysis with Pmart

#' Filter features where sample signal is NOT at least fold_threshold x blank signal
#'
#' @param omicsData  A pmartR metabData (or lipidData) object
#' @param blank_label  The value in your SampleType column that identifies blanks
#' @param sample_type_cname  Column name in f_data that contains sample type labels
#' @param fold_threshold  Minimum fold-change of median(sample) / median(blank) to KEEP a feature
#' @param remove_blanks  If TRUE, blank samples are removed from the returned object
#'
#' @return A filtered metabData object


#filter features from raw data 
blank_fold_filter_raw <- function(edata,
                                  fdata,
                                  edata_cname,
                                  fdata_cname,
                                  blank_label       = "Blank",
                                  sample_type_cname = "SampleType",
                                  fold_threshold    = 3,
                                  remove_blanks     = TRUE) {
  
  # --- 1. Validate inputs ---
  if (!sample_type_cname %in% colnames(fdata)) {
    stop(paste0("Column '", sample_type_cname, "' not found in fdata. ",
                "Available columns: ", paste(colnames(fdata), collapse = ", ")))
  }
  
  blank_samps  <- fdata[[fdata_cname]][fdata[[sample_type_cname]] == blank_label]
  sample_samps <- fdata[[fdata_cname]][fdata[[sample_type_cname]] != blank_label]
  
  if (length(blank_samps) == 0) {
    stop(paste0("No samples found with '", sample_type_cname, "' == '", blank_label, "'. ",
                "Check your blank_label argument."))
  }
  
  message(sprintf("Found %d blank sample(s) and %d non-blank sample(s).",
                  length(blank_samps), length(sample_samps)))
  
  # --- 2. Compute median signals (raw abundance scale) ---
  feature_col <- edata[[edata_cname]]
  blank_mat   <- as.matrix(edata[, blank_samps,  drop = FALSE])
  sample_mat  <- as.matrix(edata[, sample_samps, drop = FALSE])
  
  median_blank  <- apply(blank_mat,  1, median, na.rm = TRUE)
  median_sample <- apply(sample_mat, 1, median, na.rm = TRUE)
  
  # --- 3. Compute fold-change and decide what to keep ---
  fold_change <- median_sample / median_blank
  
  # Keep if fold_change >= threshold OR blank is 0/NA (feature not in blank)
  keep_flag <- is.na(fold_change) | is.nan(fold_change) | fold_change >= fold_threshold
  
  n_total   <- nrow(edata)
  n_kept    <- sum(keep_flag)
  n_removed <- n_total - n_kept
  
  message(sprintf(
    "Fold-change filter (threshold = %sx): keeping %d / %d features, removing %d.",
    fold_threshold, n_kept, n_total, n_removed
  ))
  
  # --- 4. Build summary table ---
  filter_summary <- data.frame(
    feature       = feature_col,
    median_blank  = median_blank,
    median_sample = median_sample,
    fold_change   = fold_change,
    kept          = keep_flag
  )
  
  # --- 5. Filter features from edata ---
  edata <- edata[keep_flag, ]
  
  # --- 6. Optionally remove blank samples from both edata and fdata ---
  if (remove_blanks) {
    edata <- edata[, !colnames(edata) %in% blank_samps, drop = FALSE]
    fdata <- fdata[fdata[[sample_type_cname]] != blank_label, ]
    message(sprintf("Removed %d blank sample(s).", length(blank_samps)))
  }
  
  return(list(
    edata          = edata,
    fdata          = fdata,
    filter_summary = filter_summary
  ))
}





# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

#' Generate n shades of a base color from light to dark
generate_shades <- function(base_color, n) {
  if (n == 1) return(setNames(base_color, NULL))
  colorspace::lighten(base_color, amount = seq(0.5, -0.3, length.out = n))
}

#' Build class colors inheriting shades from superclass base colors
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

#' Build ordered class levels grouped by superclass
build_ordered_classes <- function(hierarchy, base_colors, present_classes) {
  hierarchy %>%
    filter(class %in% present_classes) %>%
    arrange(match(superclass, names(base_colors))) %>%
    pull(class)
}

#' Build all annotation data frames and color lists for a given df.z
build_annotations <- function(df_z, meta.comp, class_hierarchy,
                              superclass_base_colors,
                              zone_colors, treatment_colors, timepoint_colors) {
  
  present_rows <- rownames(df_z)
  
  # superclass
  meta_sup <- meta.comp %>%
    select(`Super class`) %>%
    rename(superclass = `Super class`) %>%
    filter(!is.na(superclass), rownames(.) %in% present_rows) %>%
    mutate(superclass = factor(superclass, levels = names(superclass_base_colors)))
  
  # class
  present_classes  <- unique(na.omit(meta.comp$`Main class`[rownames(meta.comp) %in% present_rows]))
  ordered_cls      <- build_ordered_classes(class_hierarchy, superclass_base_colors, present_classes)
  
  meta_class <- meta.comp %>%
    select(`Main class`) %>%
    rename(class = `Main class`) %>%
    filter(!is.na(class), rownames(.) %in% present_rows) %>%
    mutate(class = factor(class, levels = ordered_cls)) %>%
    arrange(class)
  
  # subclass
  present_subclass <- levels(factor(na.omit(
    meta.comp$`Sub class`[rownames(meta.comp) %in% present_rows]
  )))
  
  meta_subclass <- meta.comp %>%
    select(`Sub class`) %>%
    rename(subclass = `Sub class`) %>%
    filter(!is.na(subclass), rownames(.) %in% present_rows) %>%
    mutate(subclass = factor(subclass))
  
  # colors
  sup_colors     <- superclass_base_colors[names(superclass_base_colors) %in% levels(meta_sup$superclass)]
  cls_colors     <- build_class_colors(class_hierarchy, superclass_base_colors, ordered_cls)
  subcls_colors  <- setNames(pals::polychrome(length(present_subclass)), present_subclass)
  
  # annotation color lists
  ann_sup <- list(
    Zone = zone_colors, Treatment = treatment_colors, Timepoint = timepoint_colors,
    superclass = sup_colors
  )
  ann_cls <- list(
    Zone = zone_colors, Treatment = treatment_colors, Timepoint = timepoint_colors,
    class = cls_colors
  )
  ann_sub <- list(
    Zone = zone_colors, Treatment = treatment_colors, Timepoint = timepoint_colors,
    subclass = subcls_colors
  )
  
  list(
    meta_sup      = meta_sup,
    meta_class    = meta_class,
    meta_subclass = meta_subclass,
    ann_sup       = ann_sup,
    ann_cls       = ann_cls,
    ann_sub       = ann_sub,
    row_order_cls = rownames(meta_class)  # rows sorted by class
  )
}

compute_global_breaks <- function(..., n_breaks = 100) {
  # pass all df_z datasets you want to share a scale
  all_data <- list(...)
  all_vals <- unlist(lapply(all_data, function(d) as.vector(as.matrix(d))))
  all_vals <- all_vals[!is.na(all_vals)]
  
  # use a symmetric range so 0 is always the midpoint (white)
  limit <- max(abs(quantile(all_vals, c(0.01, 0.99))))  # trim extreme outliers
  seq(-limit, limit, length.out = n_breaks)
}

compute_global_breaks <- function(..., n_breaks = 100) {
  all_data <- list(...)
  all_vals <- unlist(lapply(all_data, function(d) as.vector(as.matrix(d))))
  all_vals <- all_vals[!is.na(all_vals)]
  
  # use absolute max to guarantee ALL values fall within break range
  limit <- max(abs(all_vals))
  seq(-limit, limit, length.out = n_breaks)
}


#' Save a single pheatmap to png
save_heatmap <- function(filepath, data, annotation_col, annotation_row,
                         annotation_colors, row_order = NULL,
                         width = 25, height = 13, fontsize = 14,
                         cluster_rows = TRUE,
                         breaks = breaks) {          # ← add breaks argument
  
  annotation_col <- annotation_col[colnames(data), , drop = FALSE]
  
  if (!is.null(row_order)) {
    row_order      <- row_order[row_order %in% rownames(data)]
    data           <- data[row_order, ]
    annotation_row <- annotation_row[row_order, , drop = FALSE]
    cluster_rows   <- FALSE
  }
  
  annotation_row <- annotation_row[rownames(data), , drop = FALSE]
  
  # build color palette matched to breaks
  if (!is.null(breaks)) {
    n_colors <- length(breaks) - 1
    colors   <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(n_colors)
    # blue = low, white = 0, red = high
  } else {
    colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  }
  
  if (!is.null(breaks)) {
    n_colors <- length(breaks) - 1
    colors   <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(n_colors)
    
    # Assert they match — this will error loudly if something is wrong
    stopifnot(length(colors) == length(breaks) - 1)
    
    message(sprintf("breaks: %d values [%.2f, %.2f] | colors: %d",
                    length(breaks), min(breaks), max(breaks), length(colors)))
  }
  
  if (!is.null(breaks)) {
    data[data < min(breaks)] <- min(breaks)
    data[data > max(breaks)] <- max(breaks)
  }
  
  pdf(filepath, width = width, height = height)
  pheatmap(
    data,
    cluster_rows      = cluster_rows,
    cluster_cols      = TRUE,
    annotation_col    = annotation_col,
    annotation_row    = annotation_row,
    annotation_colors = annotation_colors,
    show_colnames     = FALSE,
    drop_levels       = FALSE,
    fontsize          = fontsize,
    breaks            = breaks,       # ← pass to pheatmap
    color             = colors        # ← pass matched colors
  )
  dev.off()
  message("Saved: ", filepath)
}



#' Produce all 6 heatmaps (superclass/class/subclass × all/sig) for a dataset
make_all_heatmaps <- function(df_z, df_z_sig, meta, meta.comp,
                              class_hierarchy, superclass_base_colors,
                              zone_colors, treatment_colors, timepoint_colors,
                              outdir, label,
                              width = 15, height_all = 13, height_sig = 10,
                              breaks = breaks) {
  
  # build annotations for full and sig datasets
  ann     <- build_annotations(df_z,     meta.comp, class_hierarchy,
                               superclass_base_colors,
                               zone_colors, treatment_colors, timepoint_colors)
  ann_sig <- build_annotations(df_z_sig, meta.comp, class_hierarchy,
                               superclass_base_colors,
                               zone_colors, treatment_colors, timepoint_colors)
  
  # align meta to data columns
  meta_aligned     <- meta[colnames(df_z),     , drop = FALSE]
  meta_aligned_sig <- meta[colnames(df_z_sig), , drop = FALSE]
  
  plots <- list(
    list(fp = file.path(outdir, paste0("other/", label, "_met-superclass.pdf")),
         data = df_z,     ann_col = meta_aligned,     ann_row = ann$meta_sup,
         ann_colors = ann$ann_sup,     h = height_all, sig = FALSE),
    
    list(fp = file.path(outdir, paste0("Fig3_", label, "_met-superclass-sig.pdf")),
         data = df_z_sig, ann_col = meta_aligned_sig, ann_row = ann_sig$meta_sup,
         ann_colors = ann_sig$ann_sup, h = height_sig, sig = TRUE),
    
    list(fp = file.path(outdir, paste0("other/",label, "_met-mainclass.pdf")),
         data = df_z,     ann_col = meta_aligned,     ann_row = ann$meta_class,
         ann_colors = ann$ann_cls,     h = height_all, sig = FALSE),
    
    list(fp = file.path(outdir, paste0("other/",label, "_met-mainclass-sig.pdf")),
         data = df_z_sig, ann_col = meta_aligned_sig, ann_row = ann_sig$meta_class,
         ann_colors = ann_sig$ann_cls, h = height_sig, sig = TRUE),
    
    list(fp = file.path(outdir, paste0("other/",label, "_met-subclass.pdf")),
         data = df_z,     ann_col = meta_aligned,     ann_row = ann$meta_subclass,
         ann_colors = ann$ann_sub,     h = height_all, sig = FALSE),
    
    list(fp = file.path(outdir, paste0("other/",label, "_met-subclass-sig.pdf")),
         data = df_z_sig, ann_col = meta_aligned_sig, ann_row = ann_sig$meta_subclass,
         ann_colors = ann_sig$ann_sub, h = height_sig, sig = TRUE)
  )
  
  for (p in plots) {
    save_heatmap(
      filepath          = p$fp,
      data              = p$data,
      annotation_col    = p$ann_col,
      annotation_row    = p$ann_row,
      annotation_colors = p$ann_colors,
      row_order         = p$row_order %||% NULL,
      width             = width,
      height            = p$h,
      breaks            = breaks
    )
  }
}

# null coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

#filter dotplot
filter_and_plot <- function(kegg_result, 
                            remove_patterns = NULL,
                            min_count = 2,
                            title = "") {
  
  filtered <- kegg_result
  filtered@result <- kegg_result@result %>%
    filter(Count >= min_count)
  
  if (!is.null(remove_patterns)) {
    pattern <- paste(remove_patterns, collapse = "|")
    filtered@result <- filtered@result %>%
      filter(!grepl(pattern, Description, ignore.case = TRUE))
  }
  
  enrichplot::dotplot(filtered) +
    ggtitle(title) +
    theme_bw()
}
