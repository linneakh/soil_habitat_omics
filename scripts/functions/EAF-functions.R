# ============================
# Functions for EAF analysis
# ============================

#' Subset phyloseq to SIP 13C samples and return OTU + taxonomy matrices
prep_sip_ps <- function(ps) {
  ps_s <- subset_samples(ps, SIP == "SIP" & Isotope == "13C")
  ps_s <- prune_taxa(taxa_sums(ps_s) > 0, ps_s)
  
  otu_s <- as(otu_table(ps_s), "matrix")
  if (taxa_are_rows(ps_s)) otu_s <- t(otu_s)
  
  ps_s_ra <- transform_sample_counts(ps_s, function(x) x / sum(x))
  otu_s_ra <- as(otu_table(ps_s_ra), "matrix")
  if (taxa_are_rows(ps_s_ra)) otu_s_ra <- t(otu_s_ra)
  
  list(
    ps       = ps_s,
    otu      = otu_s,
    otu_ra   = otu_s_ra,
    taxa     = as(tax_table(ps_s), "matrix"),
    metadata = as(sample_data(ps_s), "data.frame")
  )
}

#' Build summarised EAF x family data for a given habitat filter
make_family_sum <- function(data, habitat_pattern, tax_col = "family",
                            moisture_levels = c("Normal", "Drought")) {
  data %>%
    filter(grepl(habitat_pattern, Habitat)) %>%
    mutate(Moisture = factor(Moisture, levels = moisture_levels)) %>%
    group_by(Habitat, Moisture, Tube, .data[[tax_col]]) %>%
    summarise(sum = sum(weighted_eaf), .groups = "drop")
}

#' Boxplot of weighted EAF by family for top families
library(tidytext)

plot_top_boxplot <- function(data, tax_col = "family",
                             col_values = col_list_rhizo,
                             min_median = 0.001) {
  top_fams <- data %>%
    group_by(.data[[tax_col]], Habitat, Moisture) %>%
    summarise(median_sum = median(sum, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[tax_col]]) %>%
    summarise(max_median = max(median_sum, na.rm = TRUE), .groups = "drop") %>%
    filter(max_median > min_median) %>%
    arrange(desc(max_median)) %>%
    pull(.data[[tax_col]])
  
  data %>%
    filter(.data[[tax_col]] %in% top_fams) %>%
    group_by(Moisture) %>%
    mutate(
      median_per_group = median(sum[match(.data[[tax_col]], .data[[tax_col]])], na.rm = TRUE),
      "{tax_col}" := reorder_within(.data[[tax_col]], -sum, Moisture, fun = median)
    ) %>%
    ungroup() %>%
    ggplot(aes(x = .data[[tax_col]], y = sum, fill = Habitat)) +
    geom_boxplot(position = position_dodge(0.75)) +
    
    facet_wrap(~Moisture, scales = "free_x") +  # free_x allows different x orders per facet
    scale_x_reordered() +  
    stat_summary(fun = mean, geom = "point", 
                 shape = 23, size = 1, 
                 position = position_dodge(0.75),
                 color = "black", fill = "white",
                 aes(group = Habitat)) +
    labs(y = "Weighted EAF (sum)") +
    scale_fill_manual(values = col_values) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


plot_top_boxplot_r_d <- function(data, tax_col = "family",
                                 col_values = col_list,
                                 min_median = 0.001,
                                 title = title) {
  
  # Step 1: Identify top taxa by their max median across groups
  top_fams <- data %>%
    group_by(.data[[tax_col]], Habitat, Moisture) %>%
    summarise(median_sum = median(sum, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[tax_col]]) %>%
    summarise(max_median = max(median_sum, na.rm = TRUE), .groups = "drop") %>%
    filter(max_median > min_median) %>%
    arrange(desc(max_median)) %>%
    pull(.data[[tax_col]])
  
  # Step 2: Compute mean per tax_col × Moisture for ordering within each facet
  order_vals <- data %>%
    filter(.data[[tax_col]] %in% top_fams) %>%
    group_by(.data[[tax_col]], Moisture) %>%
    summarise(order_stat = median(sum, na.rm = TRUE), .groups = "drop")
  
  # Compute the correct factor levels outside the pipe
  fct_levels <- order_vals %>%
    arrange(Moisture, desc(order_stat)) %>%
    mutate(facet_label = paste0(.data[[tax_col]], "___", Moisture)) %>%
    pull(facet_label)
  
  # Step 3: Join ordering values and apply reorder_within
  data %>%
    filter(.data[[tax_col]] %in% top_fams) %>%
    mutate("{tax_col}" := as.character(.data[[tax_col]])) %>%
    left_join(order_vals, by = c(tax_col, "Moisture")) %>%
    mutate(
      "{tax_col}" := factor(
        paste0(.data[[tax_col]], "___", Moisture),
        levels = fct_levels
      )
    ) %>%
    ggplot(aes(x = .data[[tax_col]], y = sum, fill = Habitat)) +
    geom_boxplot(position = position_dodge(0.75),
                 linewidth = 0.2) +
    facet_wrap(~Moisture, scales = "free_x") +
    scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
    # stat_summary(fun = mean, geom = "point",
    #              shape = 23, size = 1,
    #              position = position_dodge(0.75),
    #              color = "black", fill = "white",
    #              aes(group = Habitat)) +
    labs(y = title) +
    scale_fill_manual(values = col_values) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = c(0.002, 0.998),
      legend.justification = c("left", "top"),
      legend.title = element_blank()
    )
}

plot_top_boxplot_r_d_per <- function(data, tax_col = "family",
                                     col_values = col_list,
                                     col_values2 = col_list_2,
                                     top_pct = 0.95,
                                     title = title,
                                     group1 = NULL,  # single habitat e.g. "13C-Detritus"
                                     group2 = NULL) { # combined habitat e.g. "12C-Rhizo + 13C-Detritus"
  
  # Step 1: Identify top taxa by cumulative % of max median across groups
  top_fams <- data %>%
    group_by(.data[[tax_col]], Habitat, Moisture) %>%
    summarise(median_sum = median(sum, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[tax_col]]) %>%
    summarise(max_median = max(median_sum, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_median)) %>%
    mutate(cumulative_pct = cumsum(max_median) / sum(max_median)) %>%
    filter(cumulative_pct <= top_pct) %>%
    filter(.data[[tax_col]] != "unassigned" &
             .data[[tax_col]] != "Chloroplast") %>%
    pull(.data[[tax_col]])
  
  # Step 2: Compute median per tax_col x Moisture for ordering within each facet
  order_vals <- data %>%
    filter(.data[[tax_col]] %in% top_fams) %>%
    group_by(.data[[tax_col]], Moisture) %>%
    summarise(order_stat = median(sum, na.rm = TRUE), .groups = "drop")
  
  # Compute the correct factor levels outside the pipe
  fct_levels <- order_vals %>%
    arrange(Moisture, desc(order_stat)) %>%
    mutate(facet_label = paste0(.data[[tax_col]], "___", Moisture)) %>%
    pull(facet_label)
  
  # Step 3: Build plot and save to object
  result_plot <- data %>%
    filter(.data[[tax_col]] %in% top_fams) %>%
    mutate("{tax_col}" := as.character(.data[[tax_col]])) %>%
    left_join(order_vals, by = c(tax_col, "Moisture")) %>%
    mutate(
      "{tax_col}" := factor(
        paste0(.data[[tax_col]], "___", Moisture),
        levels = fct_levels
      )
    ) %>%
    ggplot(aes(x = .data[[tax_col]], y = sum, fill = Habitat)) +
    geom_boxplot(position = position_dodge(0.75), linewidth = 0.2) +
    facet_wrap(~Moisture, scales = "free_x") +
    scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
    labs(y = title) +
    scale_fill_manual(values = col_values) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = c(0.002, 0.998),
      legend.justification = c("left", "top"),
      legend.title = element_blank()
    )
  
  # Step 4: Calculate log2fc dataframe
  
  
  result_df <- data %>%
    filter(.data[[tax_col]] %in% top_fams) %>%
    group_by(Habitat, Moisture, .data[[tax_col]]) %>%
    summarise(med_sum = median(sum, na.rm = TRUE), .groups = "drop") %>%
    spread(key = Habitat, value = med_sum) %>%
    drop_na() %>%
    mutate(log2fc = log2(.data[[group2]]) - log2(.data[[group1]])) %>%
    arrange(Moisture, desc(log2fc))
  

  # Step 5: Build log2fc plot ordered by log2fc within each Moisture facet
  log2fc_plot <- result_df %>%
    mutate(direction = case_when(
      log2fc > 0 ~ "pos",
      log2fc < 0 ~ "neg")) %>%
    mutate("{tax_col}" := reorder_within(.data[[tax_col]], log2fc, Moisture)) %>%
    ggplot(aes(x = .data[[tax_col]], y = log2fc, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_x_reordered() +
    scale_fill_manual(values = col_values2) +
    facet_wrap(~Moisture, scales = "free_y") +
    labs(y = paste0("log2fc (", group2, " / ", group1, ")"),
         x = tax_col) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank()
    )

  combined_plot <- result_plot | log2fc_plot  # | places side by side

  return(list(
    combined_plot = combined_plot,
    boxplot = result_plot,
    log2fc_plot = log2fc_plot,
    log2fc_data = result_df
  ))
  

}



#' Stacked proportional barplot colored by top N families per group
plot_prop_barplot <- function(data, tax_col = "family", n_top = 4,
                              col_values, x_var = "Habitat",
                              facet_var = "Moisture") {
  # Summarise to mean per group
  plot_data <- data %>%
    ungroup() %>%
    select(all_of(c(tax_col, x_var, facet_var, "sum"))) %>%
    group_by(across(all_of(c(tax_col, x_var, facet_var)))) %>%
    summarise(sum = mean(sum), .groups = "drop") %>%
    filter(sum > 0)
  
  # Top families per group combination
  top_families <- plot_data %>%
    group_by(across(all_of(c(x_var, facet_var)))) %>%
    slice_max(sum, n = n_top) %>%
    ungroup() %>%
    distinct(.data[[tax_col]]) %>%
    pull(.data[[tax_col]])
  
  n_colors <- length(top_families)
  cat("Number of top families:", n_colors, "\n")
  cat("Top families:", paste(top_families, collapse = ", "), "\n")
  
  # Normalize and order within each bar
  plot_data_norm <- plot_data %>%
    group_by(across(all_of(c(x_var, facet_var)))) %>%
    mutate(prop = sum / sum(sum)) %>%
    ungroup() %>%
    mutate(is_top = .data[[tax_col]] %in% top_families) %>%
    group_by(across(all_of(c(x_var, facet_var)))) %>%
    arrange(across(all_of(c(x_var, facet_var))), is_top, prop) %>%
    mutate(order_id = row_number()) %>%
    ungroup()
  
  # Color mapping
  family_colors <- c(
    setNames(col_values[1:n_colors], top_families),
    "non_top" = "grey90",
    "Other"   = "grey90"
  )
  
  # Dummy row for "Other" legend entry
  dummy_row <- plot_data_norm %>%
    slice(1) %>%
    mutate("{tax_col}" := "Other", sum = 0, prop = 0,
           is_top = FALSE, order_id = 0)
  
  plot_data_norm %>%
    bind_rows(dummy_row) %>%
    mutate(
      fill_family = case_when(
        .data[[tax_col]] == "Other"              ~ "Other",
        .data[[tax_col]] %in% top_families       ~ as.character(.data[[tax_col]]),
        TRUE                                     ~ "non_top"
      ),
      fill_family = factor(fill_family)
    ) %>%
    ggplot(aes(
      x     = .data[[x_var]],
      y     = prop,
      fill  = fill_family,
      group = interaction(.data[[x_var]], .data[[facet_var]], order_id)
    )) +
    geom_bar(stat = "identity", position = "stack",
             color = "black", linewidth = 0.3) +
    facet_wrap(as.formula(paste("~", facet_var))) +
    scale_fill_manual(
      values = family_colors,
      breaks = c(top_families, "Other"),
      labels = c(top_families, "Other")
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = x_var, y = "Proportional Assimilation", fill = tax_col) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#' Run lm + emmeans for habitat and moisture contrasts
run_emmeans <- function(data, formula = sum ~ Habitat * Moisture) {
  m <- lm(formula, data = data)
  cat("=== ANOVA ===\n")
  print(anova(m))
  
  em_hab   <- emmeans(m, ~ Habitat | Moisture)
  em_moist <- emmeans(m, ~ Moisture | Habitat)
  
  cat("\n=== Habitat contrasts ===\n")
  print(pairs(em_hab, adjust = "tukey"))
  cat("\n=== Moisture contrasts ===\n")
  print(pairs(em_moist, adjust = "tukey"))
  
  invisible(list(model = m, em_habitat = em_hab, em_moisture = em_moist))
}

#' Run t-tests per moisture level
run_ttests <- function(data, moisture_levels = c("Normal", "Drought")) {
  results <- lapply(moisture_levels, function(m) {
    d <- data %>% filter(Moisture == m)
    t <- t.test(sum ~ Habitat, data = d)
    cat("\n=== Moisture:", m, "===\n")
    print(t)
    t
  })
  setNames(results, moisture_levels)
}
