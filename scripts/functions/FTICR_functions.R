# ***************************************************************
# Linnea Hernandez
# Functions for fticr
# 
# ***************************************************************



#Linnea's additional functions
#run lm on data
run_lm_tukey <- function(data,
                         classes = unique(data$Class),
                         treatment_filter = NULL,
                         emmeans_formula = ~ Zone | Timepoint) {
  
  # Filter by treatment if specified
  if (!is.null(treatment_filter)) {
    data <- data %>% filter(Treatment == treatment_filter)
  }
  
  # Store results
  anova_results <- list()
  tukey_results <- list()
  
  for (cls in classes) {
    
    # Filter to class
    df_cls <- data %>% filter(Class == cls)
    
    # Skip if too few observations
    if (nrow(df_cls) < 10) {
      warning("Skipping ", cls, " - too few observations")
      next
    }
    
    # Fit linear model
    m <- tryCatch(
      lm(rel_abundance ~ Zone * Timepoint, data = df_cls),
      error = function(e) {
        warning("Model failed for ", cls, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(m)) next
    
    # ANOVA table
    aov_table <- broom::tidy(anova(m)) %>%
      mutate(Class = cls) %>%
      select(Class, term, df, sumsq, meansq, statistic, p.value)
    
    anova_results[[cls]] <- aov_table
    
    # Tukey pairwise comparisons using specified formula
    em <- tryCatch(
      emmeans(m, emmeans_formula),
      error = function(e) {
        warning("emmeans failed for ", cls, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(em)) next
  
    pairs_result <- pairs(em, adjust = "tukey")
    
    tukey_table <- broom::tidy(pairs_result) %>%
      mutate(Class = cls)
    
    tukey_results[[cls]] <- tukey_table
  }
  
  # Combine all results
  anova_combined <- bind_rows(anova_results) %>%
    group_by(term) %>%
    mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>%
    ungroup()
  
  tukey_combined <- bind_rows(tukey_results) %>%
    mutate(p.fdr = p.adjust(adj.p.value, method = "fdr"))
  
  return(list(
    anova = anova_combined,
    tukey = tukey_combined
  ))
}
