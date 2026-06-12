# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q2: Final model performance summary in caret package
# Created with claude.ai (sonnet 4.5) on 2025-10-28 by Nicole Drewitz
# Last updated 2026_01-14

# Use this to assess model performance at sampling location reserved for analysis
# =============================================================================

# =============================================================================
# WORKSPACE SETUP
# =============================================================================

# Load required libraries
required_packages <- c("caret", "pROC", "ggplot2", "here")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# MaxSSS threshold calculated by MaxEnt
MaxSSS_threshold <- 0.5092 # set from maxentResults.csv (Maximum training sensitivity plus specificity Cloglog threshold)

# paths relative to the RProject root
data_dir <- here::here("Data")
fig_dir  <- here::here("figures")

input_csv <- file.path(here::here("Data/cloudberry_suitability_at_analysis_locations.csv"))

data <- read.csv(input_csv)

cat("Available columns:\n")
print(names(data))
cat("\n")

save_plots <- TRUE
Climate_Scenario <- "current"

# =============================================================================
# FORMAT DATA FOR ROC
# =============================================================================

# Use presence.absence as observed and run_number column as predicted values
observed  <- data$presence.absence
predicted <- data[[Climate_Scenario]]

# Remove any NA values
valid_idx <- !is.na(observed) & !is.na(predicted)
observed  <- observed[valid_idx]
predicted <- predicted[valid_idx]

# =============================================================================
# PLOT SUITABILITY AT PRESENCE/ABSENCE LOCATIONS
# =============================================================================

# Read combined results data
combined_results <- read.csv(input_csv)
# Convert presence.absence to factor for better grouping in plots
combined_results$presence.absence <- factor(combined_results$presence.absence, levels = c(0, 1), labels = c("Absence", "Presence"))

# Create boxplot comparing habitat suitability for presence vs absence
p_box <- ggplot(combined_results, aes(x = presence.absence, y = .data[[Climate_Scenario]], fill = presence.absence)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 16, outlier.size = 2) +
  scale_fill_manual(values = c("Absence" = "#7986CB", "Presence" = "#E57373")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = paste("Habitat Suitability by Presence/Absence -", Climate_Scenario),
       x = "Presence/Absence",
       y = "Habitat Suitability") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11))

print(p_box)

# Plot option 2 ================================================================
# Create density plot comparing habitat suitability for presence vs absence
p_density <- ggplot(combined_results, aes(x = .data[[Climate_Scenario]], fill = presence.absence)) +
  geom_density(alpha = 0.6, adjust = 0.4) + # adjust, adjusts the density curve smoothing
  scale_fill_manual(values = c("Absence" = "#7986CB", "Presence" = "#E57273"),
                    labels = c("True Absences", "True Presences")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "Presence and absence density \nacross threshold habitat suitability values",
       x = "Modelled Habitat suitability",
       y = "Density",
       fill = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 13),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.key.size = unit(0.8, "cm"),
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank())

print(p_density)

# =============================================================================
# ROC ANALYSIS FOR EXISTING EXTRACTED CSV
# =============================================================================

# Create ROC object
roc_obj <- roc(observed, predicted, levels = c(0, 1), direction = "<")
# Calculate Youden's index to find optimal threshold
coords_result <- coords(
  roc_obj, "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")
)

# Get AUC value
auc_value <- auc(roc_obj)
# Prepare ROC data frame
roc_data <- data.frame(
  sensitivity  = roc_obj$sensitivities,
  specificity  = roc_obj$specificities
)

# Plot
roc_plot <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "gray") +
  geom_point(
    aes(x = 1 - coords_result$specificity,
        y = coords_result$sensitivity),
    color = "red", size = 4, shape = 19
  ) +
  annotate(
    "text",
    x = 1 - coords_result$specificity,
    y = coords_result$sensitivity - 0.08,
    label = sprintf(
      "Youden Optimal Threshold: %.3f\nSens: %.3f, Spec: %.3f",
      coords_result$threshold,
      coords_result$sensitivity,
      coords_result$specificity
    ),
    hjust = 0.05, size = 5, color = "red"
  ) +
  annotate(
    "text",
    x = 0.7, y = 0.6,
    label = sprintf("Test AUC: %.4f", auc_value),
    size = 6, fontface = "bold"
  ) +
  labs(
    title = sprintf("ROC Curve at True Presence/Absence Locations"),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title  = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

print(roc_plot)

# =============================================================================
# MAXENT OPTIMAL THRESHOLD (MAXSS= MAX TRAINING SENSITIVY + SPECIFICITY)
# =============================================================================

# Turn probabilities into class predictions (0/1 to match observed)
pred_class <- ifelse(predicted >= MaxSSS_threshold, 1, 0)

# Make sure they are factors with the same level order
obs_factor  <- factor(observed,  levels = c(0, 1))
pred_factor <- factor(pred_class, levels = c(0, 1))

cmatrix <- confusionMatrix(pred_factor, obs_factor)
cmatrix

# convert to table figure =====================================================
cm_table <- as.data.frame(cmatrix$table)
# cm$table has rows = Prediction, columns = Reference

# use this label for counts only
tp <- cmatrix$table["1", "1"]
fn <- cmatrix$table["0", "1"]
fp <- cmatrix$table["1", "0"]
tn <- cmatrix$table["0", "0"]
cm_label <- sprintf("TN: %d   FP: %d\nFN: %d   TP: %d", tn, fp, fn, tp)

n_absence <- sum(data$presence.absence == 0)
n_presence <- sum(data$presence.absence == 1)

# use this label for % of presences or % of absences
perc_tn <- round(100 * tn / n_absence, 1)  # % of absences correct
perc_fp <- round(100 * fp / n_absence, 1)  # % of absences predicted present
perc_fn <- round(100 * fn / n_presence, 1) # % of presences predicted absent
perc_tp <- round(100 * tp / n_presence, 1) # % of presences correct
cm_label_perc <- sprintf("T Absence: %d (%.1f%% abs)   F Presence: %d (%.1f%% abs)\nF Absence: %d (%.1f%% pres)   T Presence: %d (%.1f%% pres)",
                         tn, perc_tn, fp, perc_fp, fn, perc_fn, tp, perc_tp)

# add to suitability density plot =============================================
p_density_cm <- p_density +
  annotate(
    "text",
    x = max(predicted, na.rm = TRUE) * 0.01,  # adjust position as needed
    y = max(ggplot_build(p_density)$data[[1]]$density, na.rm = TRUE) * 0.6,
    label = cm_label_perc,
    hjust = -0.8,
    vjust = 8.5,
    size  = 3.5
  )

print(p_density_cm)

# =============================================================================
# SAVE PLOTS
# =============================================================================
if (save_plots) {
  # ensure figures directory exists
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  #density_file <- file.path(fig_dir, sprintf("%s_density_plot.png", run_number))
  density_file_cmatrix <- file.path(fig_dir, "Q2-b-Model_performance_density_with_confusion_matrix.png")
  density_file <- file.path(fig_dir, "Q2-b-Model_performance_density.png")
  box_file     <- file.path(fig_dir, "Q2-b-Model_performance_boxplot.png")
  roc_file     <- file.path(fig_dir, "Q2-b-Model_performance_ROC_Curve.png")

  #ggsave(density_file, plot = p_density, width = 8, height = 6, dpi = 300)
  #cat("Density plot saved as:", density_file, "\n")

  ggsave(density_file_cmatrix, plot = p_density_cm, width = 10, height = 6, dpi = 300)
  cat("Density plot with confusion matrix saved as:", density_file_cmatrix, "\n")
  
  ggsave(density_file, plot = p_density, width = 10, height = 6, dpi = 300)
  cat("Density plot with confusion matrix saved as:", density_file, "\n")

  ggsave(roc_file,     plot = roc_plot,  width = 8, height = 6, dpi = 400)
  cat("\nROC plot saved as:", roc_file, "\n")

  ggsave(box_file,     plot = p_box,     width = 8, height = 6, dpi = 300)
  cat("Boxplot saved as:", box_file, "\n")
}
