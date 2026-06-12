# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q2
# Selects best fine-tuning settings by comparing models using AICc
# Create with Claude.ai on 04-10-2025 by Nicole Drewitz
# Last updated 2026-01-19

# Runs selection of best model parameters based on suitability values extracted from sampling locations reserved for model selection.
# Results: Best model = RegMult: 0.5, LTH (Linear, Threshold, Hinge)
# =============================================================================

# =============================================================================
# WORKSPACE SETUP
# =============================================================================

# Load required libraries
required_packages <- c("ggplot2", "here", "dplyr", "tictoc", "tidyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

save_plots <- FALSE

# paths relative to the RProject root
input_dir <- here::here("Input_files/MaxEnt_files")
fig_dir  <- here::here("figures")
output_dir <- here::here ("analysis_output_data")

parameters_output_csv <- file.path(input_dir, "Q2-a-model_predictions_comprehensive.csv")
stepwise_aicc_csv <- file.path (input_dir, "Q2-a-habitat_suitability_all_stepwise_runs.csv")

data <- read.csv(parameters_output_csv)
stepwise_data <- read.csv(stepwise_aicc_csv)

tic() # start script runtime counter

cat("Available columns for parameter data:\n")
print(names(data))
cat("\n")

cat("Available columns for stepwise AICc data:\n")
print(names(stepwise_data))
cat("\n")

# =============================================================================
# AICc COMPARISON OF HABITAT SUITABILITY MODELS WITH DIFFERENT PARAMETERS
# =============================================================================

cat("\n=== AICc Model Comparison ===\n")

# Load the comprehensive prediction data
model_data <- read.csv(parameters_output_csv, stringsAsFactors = FALSE)

# Extract habitat suitability columns
hs_cols <- grep("^habitat_suitability_", names(model_data), value = TRUE)
cat("Found", length(hs_cols), "habitat suitability predictions for comparison\n\n")

# Initialize results dataframe
aicc_results <- data.frame(
  Model = character(),
  reg_mult = numeric(),
  features = character(),
  K = numeric(),
  AICc = numeric(),
  Delta_AICc = numeric(),
  Weight = numeric(),
  stringsAsFactors = FALSE
)

# Calculate AICc for each model using logistic regression
for (i in seq_along(hs_cols)) {
  col <- hs_cols[i]
  model_name <- sub("^habitat_suitability_", "", col)

  # Parse regularization multiplier and features from model name
  parts <- strsplit(model_name, "_")[[1]]

  if (grepl("regmult", parts[1], ignore.case = TRUE)) {
    if (length(parts) >= 3 && grepl("^[0-9]+$", parts[2]) && grepl("^[0-9]+$", parts[3])) {
      reg_mult <- as.numeric(paste0(parts[2], ".", parts[3]))
      features <- paste(parts[4:length(parts)], collapse = "_")
    } else if (length(parts) >= 2 && grepl("^[0-9]+$", parts[2])) {
      reg_mult <- as.numeric(parts[2])
      features <- paste(parts[3:length(parts)], collapse = "_")
    } else {
      reg_mult <- NA
      features <- paste(parts[-1], collapse = "_")
    }
  } else {
    reg_mult <- NA
    features <- paste(parts, collapse = "_")
  }

  # Fit logistic regression model with habitat suitability as predictor
  fit <- glm(presence.absence ~ model_data[[col]],
             family = binomial(link = "logit"),
             data = model_data)

  # Extract AIC and number of parameters
  aic <- AIC(fit)
  k <- length(coef(fit))
  n <- nrow(model_data)

  # Calculate AICc (corrected for small sample sizes)
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)

  # Store results
  aicc_results <- rbind(aicc_results, data.frame(
    Model = model_name,
    reg_mult = reg_mult,
    features = features,
    K = k,
    AICc = aicc,
    Delta_AICc = NA,
    Weight = NA,
    stringsAsFactors = FALSE
  ))
}

# Calculate Delta AICc and Akaike weights
min_aicc <- min(aicc_results$AICc)
aicc_results$Delta_AICc <- aicc_results$AICc - min_aicc

# Calculate Akaike weights
aicc_results$Weight <- exp(-0.5 * aicc_results$Delta_AICc) /
  sum(exp(-0.5 * aicc_results$Delta_AICc))

# Sort by regularization multiplier then by AICc
aicc_results <- aicc_results[order(aicc_results$reg_mult, aicc_results$AICc), ]
rownames(aicc_results) <- NULL

# Display results
cat("Model Comparison Results (sorted by Regularization Multiplier):\n")
cat(strrep("=", 100), "\n")
print(aicc_results, digits = 4)
cat("\n")

# =============================================================================
# IDENTIFY BEST MODEL PARAMETERS
# =============================================================================

best_idx <- which.min(aicc_results$AICc)
best_model_full <- aicc_results[best_idx, ]
cat("Best Model (lowest AICc):", best_model_full$Model, "\n")
cat("Regularization Multiplier:", best_model_full$reg_mult, "\n")
cat("Feature Classes:", best_model_full$features, "\n")
cat("AICc =", round(best_model_full$AICc, 2), "\n")
cat("Akaike Weight =", round(best_model_full$Weight, 4), "\n\n")

# Define default model (typically regmult_1_LQH or similar)
# Adjust this based on your MaxEnt default settings
default_model <- aicc_results %>%
  filter(reg_mult == 1.0, features == "LQH") %>%
  slice(1)

if (nrow(default_model) == 0) {
  # If no exact match, find closest to default
  default_model <- aicc_results %>%
    filter(reg_mult == min(aicc_results$reg_mult)) %>%
    filter(features %in% c("LQH", "LQ", "LH")) %>%
    slice(1)
}

cat("Default Model (MaxEnt default settings):", default_model$Model, "\n")
cat("AICc =", round(default_model$AICc, 2), "\n\n")

# SAVE PARAMETER TUNING AICC RESULTS===========================================
output_aicc_csv <- file.path(output_dir, "Q2-a-AICc_parameter_tuning_comparison_results.csv")
write.csv(aicc_results, output_aicc_csv, row.names = FALSE)
cat("Parameter tuning AICc results saved to:", output_aicc_csv, "\n")

# =============================================================================
# AICc ANALYSIS FOR STEPWISE VARIABLE REMOVAL
# =============================================================================

cat("\n=== AICc Analysis for Stepwise Variable Removal ===\n")

# Function to calculate AICc from a logistic regression model
calc_AICc <- function(model) {
  n <- length(model$fitted.values)
  k <- length(coef(model))
  AIC_val <- AIC(model)
  AICc_val <- AIC_val + (2 * k * (k + 1)) / (n - k - 1)
  return(AICc_val)
}

# Calculate AICc for each habitat suitability column
habitat_columns <- names(stepwise_data)[5:ncol(stepwise_data)]
stepwise_aicc_results <- data.frame(
  model_name = character(),
  AICc = numeric(),
  stringsAsFactors = FALSE
)

for (col_name in habitat_columns) {
  model_subset <- stepwise_data[!is.na(stepwise_data[[col_name]]), c("presence.absence", col_name)]

  if (nrow(model_subset) > 0) {
    tryCatch({
      glm_model <- glm(as.formula(paste("presence.absence ~", col_name)),
                       data = model_subset, family = binomial())
      stepwise_aicc_results <- rbind(stepwise_aicc_results,
                                     data.frame(model_name = col_name,
                                                AICc = calc_AICc(glm_model)))
    }, error = function(e) {
      cat("  Error calculating AICc for", col_name, ":", conditionMessage(e), "\n")
    })
  }
}

# Add number of variables and sort by AICc
stepwise_aicc_results <- stepwise_aicc_results %>%
  arrange(AICc) %>%
  mutate(number_of_variables = 11 - as.integer(sub("RUN(\\d+)", "\\1", model_name)))

# Identify best model
best_stepwise_full <- stepwise_aicc_results[1, ]

# PRINT RESULTS ===============================================================
cat("\n=== Stepwise AICc Comparison Results ===\n")
print(stepwise_aicc_results[, c("model_name", "number_of_variables", "AICc")])
cat("\n=== Best Stepwise Model (Lowest AICc) ===\n")
print(best_stepwise_full)

# Save results
output_stepwise_csv <- file.path(output_dir, "Q2-a-AICc_stepwise_comparison_results.csv")
write.csv(stepwise_aicc_results, output_stepwise_csv, row.names = FALSE)
cat("\nStepwise AICc results saved to:", output_stepwise_csv, "\n")

# =============================================================================
# FORMATTING FOR PARAMETER VISUALS
# =============================================================================

# Identify top 10 feature classes by lowest mean AICc
feature_performance_summary <- aicc_results %>%
  group_by(features) %>%
  summarise(mean_aicc = mean(AICc), .groups = "drop") %>%
  arrange(mean_aicc)

top_10_features <- head(feature_performance_summary$features, 10)

# Colorblind-friendly palette for top 10
top_10_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                   "#D55E00", "#CC79A7", "#000000", "#E1BE6C", "#B35806")

# Create color mapping
color_map <- setNames(top_10_colors, top_10_features)

# Add color grouping to results
aicc_results <- aicc_results %>%
  mutate(
    feature_group = ifelse(features %in% top_10_features, features, "Other"),
    feature_color = ifelse(features %in% top_10_features, color_map[features], "gray70")
  )

# Find default model
default_model <- aicc_results %>%
  filter(reg_mult == 1.0, features == "LQPH") %>%
  slice(1)

if (nrow(default_model) == 0) {
  default_model <- aicc_results %>%
    filter(reg_mult == 1.0) %>%
    arrange(AICc) %>%
    slice(1)
}

# =============================================================================
# PLOT 1: MODEL PARAMETERS SELECTION PLOTS
# =============================================================================
cat ("\n=== Model Parameters Selection Plot ===\n")

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
           label = sprintf(" Best\n %s\nRM: %.1f", best_model_full$features, best_model_full$reg_mult),
           hjust = -0.5, vjust = 0.2,
           fontface = "bold", size = 3, fill = alpha("white", 0.9),
           label.size = 0.5, label.padding = unit(0.5, "lines")) +

  annotate("label",
           x = default_model$reg_mult,
           y = default_model$AICc,
           label = sprintf("Default\n%s\nRM: %.1f", default_model$features, default_model$reg_mult),
           hjust = 1.25, vjust = -0.5,
           fontface = "bold", size = 3, fill = alpha("white", 0.9),
           label.size = 0.5, label.padding = unit(0.3, "lines")) +

  labs(
    title = "Model Selection: AICc across Regularization Multipliers and Features",
    x = "Regularization Multiplier",
    y = "AICc", # lower is better, >2 difference to be significant
    caption = "Top 10 feature classes (by mean AICc) shown in color; others in grey"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 9, hjust = 0),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white")
  )

fine_tuning_plot

toc()
cat("\n=== AICc Analysis Complete ===\n")

# =============================================================================
# ADDITIONAL PLOTS
# =============================================================================

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
# PLOT 4: STEPWISE VARIABLE REMOVAL
# =============================================================================

cat("=== Creating Stepwise Variable Removal Plot ===\n")

stepwise_plot <- ggplot(stepwise_aicc_results, aes(x = number_of_variables, y = AICc)) +
  geom_line(linewidth = 1, color = "black") +
  geom_point(size = 3, color = "black") +
  geom_point(data = best_stepwise_full, aes(x = number_of_variables, y = AICc),
             color = "red", size = 4, shape = 18) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Stepwise Variable Removal: AICc vs Number of Variables",
    x = "Number of Variables",
    y = "AICc" # lower is better, needs >2 difference
  ) +
  scale_x_continuous(breaks = seq(min(stepwise_aicc_results$number_of_variables),
                                  max(stepwise_aicc_results$number_of_variables),
                                  by = 1)) +
  # Best model line
  geom_hline(yintercept = best_stepwise_full$AICc, linetype = "dashed",
             color = "red", linewidth = 0.8) +
  # Best model annotation
  annotate("text",
           x = max(stepwise_aicc_results$number_of_variables),
           y = best_stepwise_full$AICc,
           label = paste0("Best model: ", best_stepwise_full$model_name,
                          " (", best_stepwise_full$number_of_variables, " variables)"),
           hjust = 1.2, vjust = -0.5, color = "red", fontface = "bold", size = 3.5)

print(stepwise_plot)

# =============================================================================
# SUMMARY STATISTICS AND TABLES
# =============================================================================

cat("\n=== Top 5 Feature Classes (by mean AICc) ===\n")
top_5 <- head(feature_performance, 5)
print(top_5)

cat("\n=== Bottom 5 Feature Classes (by mean AICc) ===\n")
bottom_5 <- tail(feature_performance, 5)
print(bottom_5)

cat("\n=== Improvement over Default Settings ===\n")
improvement_pct <- ((default_model$AICc - best_model_full$AICc) / default_model$AICc) * 100
cat(sprintf("Default model AICc: %.2f\n", default_model$AICc))
cat(sprintf("Best model AICc: %.2f\n", best_model_full$AICc))
cat(sprintf("Improvement: %.2f%%\n", improvement_pct))

# =============================================================================
# SAVE PLOTS
# =============================================================================
# name plots (use *_file for paths)
fine_tuning_output_file <- file.path(fig_dir, "Q2-a-01_fine_tuning_aicc_comparison.png")
feature_summary_file    <- file.path(fig_dir, "Q2-a-02_feature_class_performance.png")
regmult_comparison_file <- file.path(fig_dir, "Q2-a-03_regmult_comparison.png")
stepwise_plot_file <- file.path(fig_dir, "Q2-a-04_stepwise_plot.png")

if (save_plots) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(fine_tuning_output_file, plot = fine_tuning_plot, width = 7, height = 7, dpi = 300)
  cat("✓ Saved: 01_fine_tuning_output_plot.png\n")

  ggsave(feature_summary_file, plot = feature_summary_plot, width = 10, height = 12, dpi = 300)
  cat("✓ Saved: 02_feature_class_performance.png\n")

  ggsave(regmult_comparison_file, plot = regmult_comparison, width = 12, height = 8, dpi = 300)
  cat("✓ Saved: 03_regmult_comparison.png\n")

  ggsave(stepwise_plot_file, plot = stepwise_plot, width = 8, height = 6, dpi = 300)
  cat("✓ Saved: 04_stepwise_plot.png\n")

}
