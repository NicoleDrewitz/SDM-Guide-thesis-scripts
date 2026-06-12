# =============================================================================
# OPTIONAL: Fine-Tuning Parameter Analysis - 31 Feature Classes
# Run AFTER the main AICc comparison script
# Created with Claude.ai on 18-10-2025 by Nicole Drewitz
# =============================================================================

cat("\n=== Fine-Tuning Parameter Analysis ===\n")

# Verify aicc_results exists from main script
if (!exists("aicc_results")) {
  stop("Please run the main AICc comparison script first!")
}

if (!exists("best_model_full")) {
  stop("best_model_full not found. Run the main script first.")
}

if (!exists("default_model")) {
  stop("default_model not found. Run the main script first.")
}

# =============================================================================
# PLOT 1: FINE-TUNING WITH BEST/DEFAULT HIGHLIGHTED
# =============================================================================

cat("\n→ Creating fine-tuning plot with best and default models highlighted\n")

# Add highlighting variables
aicc_results_plot <- aicc_results %>%
  mutate(
    highlight = case_when(
      features == best_model_full$features ~ "Best Model",
      features == default_model$features ~ "Default Model",
      TRUE ~ "Other"
    ),
    alpha_level = case_when(
      highlight != "Other" ~ 1.0,
      TRUE ~ 0.3
    ),
    line_width = case_when(
      highlight != "Other" ~ 1.2,
      TRUE ~ 0.6
    )
  )

# Identify top 10 feature classes by lowest mean AICc
feature_performance_main <- aicc_results %>%
  group_by(features) %>%
  summarise(mean_aicc = mean(AICc), .groups = "drop") %>%
  arrange(mean_aicc)

top_10_features <- head(feature_performance_main$features, 10)

# Colorblind-friendly palette for top 10
top_10_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                   "#D55E00", "#CC79A7", "#000000", "#E1BE6C", "#B35806")

# Create color mapping
color_map <- setNames(top_10_colors, top_10_features)

# Find default model with specific criteria
default_model <- aicc_results %>%
  filter(reg_mult == 1.0, features == "LQPH") %>%
  slice(1)

if (nrow(default_model) == 0) {
  default_model <- aicc_results %>%
    filter(reg_mult == 1.0) %>%
    arrange(AICc) %>%
    slice(1)
}

fine_tuning_plot <- ggplot(aicc_results, aes(x = reg_mult, y = AICc, group = features)) +
  # Grey lines for other features
  geom_line(data = aicc_results %>% filter(!(features %in% top_10_features)),
            color = "gray70", linewidth = 0.6, alpha = 0.5) +
  # Colored lines for top 10
  geom_line(data = aicc_results %>% filter(features %in% top_10_features),
            aes(color = features), linewidth = 1.0, alpha = 0.8) +
  # Points
  geom_point(data = aicc_results %>% filter(features %in% top_10_features),
             aes(color = features), size = 2.5, alpha = 0.7) +
  geom_point(data = aicc_results %>% filter(!(features %in% top_10_features)),
             color = "gray70", size = 1.5, alpha = 0.4) +
  
  # Circle best model
  geom_point(data = best_model_full, 
             aes(x = reg_mult, y = AICc), 
             color = "black", size = 8, shape = 1, stroke = 2.5, inherit.aes = FALSE) +
  geom_point(data = best_model_full, 
             aes(x = reg_mult, y = AICc), 
             fill = color_map[best_model_full$features], color = color_map[best_model_full$features],
             size = 4, inherit.aes = FALSE) +
  
  # Circle default model
  geom_point(data = default_model, 
             aes(x = reg_mult, y = AICc), 
             color = "black", size = 8, shape = 1, stroke = 2.5, inherit.aes = FALSE) +
  geom_point(data = default_model, 
             aes(x = reg_mult, y = AICc), 
             fill = color_map[default_model$features], color = color_map[default_model$features],
             size = 4, inherit.aes = FALSE) +
  
  # Color scale
  scale_color_manual(values = top_10_colors, name = "Top 10 Feature Classes", breaks = top_10_features) +
  
  # Labels
  annotate("label", 
           x = best_model_full$reg_mult, 
           y = best_model_full$AICc, 
           label = sprintf("Best\n%s\nRM: %.1f", best_model_full$features, best_model_full$reg_mult),
           hjust = -0.1, vjust = 0.5, 
           fontface = "bold", size = 3, fill = alpha("white", 0.9), 
           label.size = 0.5, label.padding = unit(0.3, "lines")) +
  
  annotate("label", 
           x = default_model$reg_mult, 
           y = default_model$AICc, 
           label = sprintf("Default\n%s\nRM: %.1f", default_model$features, default_model$reg_mult),
           hjust = -0.1, vjust = -0.5, 
           fontface = "bold", size = 3, fill = alpha("white", 0.9), 
           label.size = 0.5, label.padding = unit(0.3, "lines")) +
  
  labs(
    title = "Model Selection: AICc across Regularization Multipliers and Features",
    subtitle = sprintf("Comparing %d feature class combinations (top 10 in color, others in grey)", 
                       length(unique(aicc_results$features))),
    x = "Regularization Multiplier",
    y = "AICc (lower is better)"
  ) +
  
  # Theme
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(fine_tuning_plot)

# =============================================================================
# PLOT 2: FEATURE CLASS PERFORMANCE SUMMARY
# =============================================================================

cat("→ Creating feature class performance summary\n")

# Calculate summary statistics for each feature class
feature_performance <- aicc_results %>%
  group_by(features) %>%
  summarise(
    mean_aicc = mean(AICc),
    min_aicc = min(AICc),
    max_aicc = max(AICc),
    range_aicc = max_aicc - min_aicc,
    sd_aicc = sd(AICc),
    .groups = "drop"
  ) %>%
  arrange(mean_aicc)

# Create feature performance plot
feature_summary_plot <- ggplot(feature_performance, aes(x = reorder(features, mean_aicc), 
                                                        y = mean_aicc)) +
  geom_segment(aes(xend = features, y = min_aicc, yend = max_aicc), 
               color = "gray60", linewidth = 1) +
  geom_point(aes(y = min_aicc), color = "#009E73", size = 3, alpha = 0.7) +
  geom_point(aes(y = mean_aicc), color = "black", size = 2.5) +
  geom_point(aes(y = max_aicc), color = "#CC79A7", size = 3, alpha = 0.7) +
  
  # Highlight best model's feature class
  geom_point(data = feature_performance %>% filter(features == best_model_full$features),
             aes(y = mean_aicc), color = "#D55E00", size = 5, shape = 18) +
  
  coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Feature Class Performance Summary",
    subtitle = "Range of AICc values across all regularization multipliers",
    x = "Feature Class",
    y = "AICc",
    caption = "Black dot = mean | Green = best (min) | Pink = worst (max) | Orange diamond = overall best model"
  )

print(feature_summary_plot)

# =============================================================================
# PLOT 3: REGULARIZATION MULTIPLIER COMPARISON
# =============================================================================

cat("→ Creating regularization multiplier comparison\n")

# Calculate summary statistics for each regularization multiplier
regmult_performance <- aicc_results %>%
  group_by(reg_mult) %>%
  summarise(
    mean_aicc = mean(AICc),
    median_aicc = median(AICc),
    min_aicc = min(AICc),
    max_aicc = max(AICc),
    sd_aicc = sd(AICc),
    q25 = quantile(AICc, 0.25),
    q75 = quantile(AICc, 0.75),
    n_models = n(),
    .groups = "drop"
  )

# Violin and boxplot comparison
regmult_comparison <- ggplot(aicc_results, aes(x = factor(reg_mult), y = AICc)) +
  geom_violin(fill = "gray85", color = "gray60", alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = AICc), width = 0.15, alpha = 0.4, size = 1.5) +
  
  # Mean line
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, 
               color = "#E69F00", stroke = 1.5) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), 
               color = "#E69F00", linewidth = 1, linetype = "dashed") +
  
  # Highlight best model
  geom_point(data = best_model_full, 
             aes(x = factor(reg_mult), y = AICc), 
             color = "#D55E00", size = 6, shape = 18, 
             stroke = 1.5, inherit.aes = FALSE) +
  
  # Highlight default model
  geom_point(data = default_model, 
             aes(x = factor(reg_mult), y = AICc), 
             color = "#0072B2", size = 6, shape = 17, 
             stroke = 1.5, inherit.aes = FALSE) +
  
  scale_color_viridis_c(option = "viridis", end = 0.9) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray50")
  ) +
  
  labs(
    title = "Regularization Multiplier Performance Comparison",
    subtitle = sprintf("Distribution of AICc across %d feature classes", 
                       length(unique(aicc_results$features))),
    x = "Regularization Multiplier",
    y = "AICc (lower is better)",
    color = "AICc",
    caption = "Orange diamonds = mean | Large diamond = best model | Large triangle = default model"
  )

print(regmult_comparison)

# =============================================================================
# PLOT 4: REGULARIZATION STATISTICS
# =============================================================================

cat("→ Creating regularization statistics summary\n")

regmult_stats_plot <- regmult_performance %>%
  pivot_longer(cols = c(mean_aicc, median_aicc, min_aicc, sd_aicc), 
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric_label = case_when(
      metric == "mean_aicc" ~ "Mean AICc",
      metric == "median_aicc" ~ "Median AICc",
      metric == "min_aicc" ~ "Best (Min) AICc",
      metric == "sd_aicc" ~ "Std Dev of AICc"
    ),
    metric_label = factor(metric_label, 
                          levels = c("Mean AICc", "Median AICc", 
                                     "Best (Min) AICc", "Std Dev of AICc"))
  ) %>%
  ggplot(aes(x = reg_mult, y = value)) +
  geom_line(linewidth = 1.2, color = "#0072B2") +
  geom_point(size = 3, color = "#0072B2") +
  
  geom_vline(xintercept = best_model_full$reg_mult, 
             linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  geom_vline(xintercept = default_model$reg_mult, 
             linetype = "dashed", color = "#0072B2", linewidth = 0.8, alpha = 0.5) +
  
  facet_wrap(~metric_label, scales = "free_y", ncol = 2) +
  
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray50")
  ) +
  
  labs(
    title = "Regularization Multiplier: Summary Statistics",
    x = "Regularization Multiplier",
    y = "Value",
    caption = "Orange dashed line = best model | Blue dashed line = default model"
  )

print(regmult_stats_plot)

# =============================================================================
# SUMMARY STATISTICS AND TABLES
# =============================================================================

cat("\n=== Top 5 Feature Classes (by mean AICc) ===\n")
top_5 <- head(feature_performance, 5)
print(top_5)

cat("\n=== Bottom 5 Feature Classes (by mean AICc) ===\n")
bottom_5 <- tail(feature_performance, 5)
print(bottom_5)

cat("\n=== Feature Classes with Least Sensitivity to Regularization ===\n")
least_sensitive <- head(feature_performance %>% arrange(range_aicc), 5)
print(least_sensitive)

cat("\n=== Regularization Multiplier Performance Summary ===\n")
print(regmult_performance %>% arrange(mean_aicc) %>%
        select(reg_mult, mean_aicc, median_aicc, min_aicc, sd_aicc))

cat("\n=== Improvement over untuned settings ===\n")
improvement_pct <- ((default_model$AICc - best_model_full$AICc) / default_model$AICc) * 100
cat(sprintf("Default model AICc: %.2f\n", default_model$AICc))
cat(sprintf("Best model AICc: %.2f\n", best_model_full$AICc))
cat(sprintf("Improvement: %.2f%%\n", improvement_pct))

# Correlation analysis
cat("\n=== Trend Analysis ===\n")
cor_mean <- cor(regmult_performance$reg_mult, regmult_performance$mean_aicc)
cor_sd <- cor(regmult_performance$reg_mult, regmult_performance$sd_aicc)

cat(sprintf("Correlation (reg_mult vs mean_aicc): %.3f\n", cor_mean))
cat(sprintf("Correlation (reg_mult vs std_dev): %.3f\n", cor_sd))

if (cor_mean > 0) {
  cat("→ Higher regularization tends to increase AICc (worse performance)\n")
} else {
  cat("→ Higher regularization tends to decrease AICc (better performance)\n")
}

if (cor_sd > 0) {
  cat("→ Higher regularization increases variability across features\n")
} else {
  cat("→ Higher regularization decreases variability across features\n")
}

# Save plots and results
cat("\n=== Saving Results ===\n")

ggsave(file.path(output_dir, "01_fine_tuning_best_vs_default.png"), fine_tuning_plot, 
       width = 12, height = 7, dpi = 300, bg = "white")
cat("✓ Saved: 01_fine_tuning_best_vs_default.png\n")

ggsave(file.path(output_dir, "02_feature_class_performance.png"), feature_summary_plot, 
       width = 10, height = 12, dpi = 300, bg = "white")
cat("✓ Saved: 02_feature_class_performance.png\n")

ggsave(file.path(output_dir, "03_regmult_comparison.png"), regmult_comparison, 
       width = 12, height = 8, dpi = 300, bg = "white")
cat("✓ Saved: 03_regmult_comparison.png\n")

ggsave(file.path(output_dir, "04_regmult_statistics.png"), regmult_stats_plot, 
       width = 12, height = 8, dpi = 300, bg = "white")
cat("✓ Saved: 04_regmult_statistics.png\n")

# Save feature performance table
write.csv(feature_performance, 
          file.path(output_dir, "feature_class_performance_summary.csv"), 
          row.names = FALSE)
cat("✓ Saved: feature_class_performance_summary.csv\n")

# Save regularization multiplier performance table
write.csv(regmult_performance, 
          file.path(output_dir, "regmult_performance_summary.csv"), 
          row.names = FALSE)
cat("✓ Saved: regmult_performance_summary.csv\n")

cat("\n=== Fine-Tuning Analysis Complete ===\n")
