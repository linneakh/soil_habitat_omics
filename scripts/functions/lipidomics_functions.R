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


# 
# compute_global_breaks <- function(..., n_breaks = 100) {
#   # pass all df_z datasets you want to share a scale
#   all_data <- list(...)
#   all_vals <- unlist(lapply(all_data, function(d) as.vector(as.matrix(d))))
#   all_vals <- all_vals[!is.na(all_vals)]
#   
#   # use a symmetric range so 0 is always the midpoint (white)
#   limit <- max(abs(quantile(all_vals, c(0.01, 0.99))))  # trim extreme outliers
#   seq(-limit, limit, length.out = n_breaks)
# }

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
  
  #annotation_colors <- class_base_colors[names(class_base_colors) %in% levels(meta_sup$superclass)]
  
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
                              class_base_colors,
                              zone_colors, treatment_colors, timepoint_colors,
                              outdir, label,
                              width = 15, height_all = 13, height_sig = 10,
                              breaks = NULL) {
  
  # ── align sample metadata to data columns ──────────────────────────────────
  meta_aligned     <- meta[colnames(df_z),     , drop = FALSE]
  meta_aligned_sig <- meta[colnames(df_z_sig), , drop = FALSE]
  
  # ── align compound metadata to data rows ───────────────────────────────────
  meta_comp_aligned     <- meta.comp[rownames(df_z),     , drop = FALSE]
  meta_comp_aligned_sig <- meta.comp[rownames(df_z_sig), , drop = FALSE]
  
  # ── screen for missing class colors — full dataset ─────────────────────────
  present_classes     <- unique(na.omit(meta_comp_aligned$Class))
  present_classes_sig <- unique(na.omit(meta_comp_aligned_sig$Class))
  all_present_classes <- unique(c(present_classes, present_classes_sig))
  
  missing_classes <- all_present_classes[!all_present_classes %in% names(class_base_colors)]
  
  if (length(missing_classes) > 0) {
    message("── Missing class colors detected ──────────────────────────────")
    message("The following classes have no color in class_base_colors:")
    message(paste(" •", missing_classes, collapse = "\n"))
    message("Adding grey (#CCCCCC) for unmatched classes — update class_base_colors to fix.")
    
    extra_colors <- setNames(
      rep("#CCCCCC", length(missing_classes)),
      missing_classes
    )
    class_base_colors <- c(class_base_colors, extra_colors)
  } else {
    message("── Class color check passed — all classes have colors ──────────")
  }
  
  # ── handle NAs in class column ─────────────────────────────────────────────
  na_count     <- sum(is.na(meta_comp_aligned$Class))
  na_count_sig <- sum(is.na(meta_comp_aligned_sig$Class))
  
  if (na_count > 0 || na_count_sig > 0) {
    message(sprintf("NA classes found: %d in full dataset, %d in sig dataset — replacing with 'Unknown'",
                    na_count, na_count_sig))
    meta_comp_aligned$Class[is.na(meta_comp_aligned$Class)]         <- "Unknown"
    meta_comp_aligned_sig$Class[is.na(meta_comp_aligned_sig$Class)] <- "Unknown"
    
    if (!"Unknown" %in% names(class_base_colors)) {
      class_base_colors["Unknown"] <- "#CCCCCC"
    }
  }
  
}



