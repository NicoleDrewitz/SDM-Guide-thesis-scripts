# =============================================================================
# Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Selects best fine-tuning settings bu comparing models using AICc
# Created with Claude.ai on 04-10-2025 by Nicole Drewitz
# Last updated 2025-10-18

# Results: Best model = RegMult: 0.5, LTH (Linear, Threshold, Hinge)
# =============================================================================

required_packages <- c("dismo", "terra", "raster", "ggplot2", "dplyr",
                       "sf", "gridExtra", "rJava", "BBmisc", "MuMIn", "tictoc", "here")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

tic()
# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# File paths macOS
base_dir <- here::here("Output_files")
setwd(base_dir)

# MaxEnt results directory from previous run
maxent_results_dir <- file.path(base_dir, "MaxEnt_Results_parameter_tuning")

# Presence and absence data files reserved for model selection
presence_file_AICc <- here::here("Input_files/V2_cloudberry_presence_select.csv")
absence_file_AICc <- here::here("Input_files/V2_cloudberry_absence_select.csv")

output_dir <- file.path(maxent_results_dir, "Model_Comparison_AICc")
dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)

# =============================================================================
# LOAD SPECIES DATA AND PREDICTIONS
# =============================================================================

cat("=== Loading Presence and Absence Data ===\n")

# Load presence data
presence_data <- read.csv(presence_file_AICc, stringsAsFactors = FALSE)
presence_points <- presence_data[, c("longitude", "latitude")]
presence_points <- na.omit(presence_points)

cat("Presence points:", nrow(presence_points), "\n")

# Load absence data
absence_data <- read.csv(absence_file_AICc, stringsAsFactors = FALSE)
absence_points <- absence_data[, c("longitude", "latitude")]
absence_points <- na.omit(absence_points)

cat("Absence points:", nrow(absence_points), "\n")

# =============================================================================
# EXTRACT PREDICTIONS FROM MAXENT RUNS
# =============================================================================

cat("\n=== Extracting Predictions from MaxEnt Runs ===\n")

# Find all model run directories
run_dirs <- list.dirs(maxent_results_dir, recursive = FALSE, full.names = TRUE)

if (length(run_dirs) == 0) {
  stop("No model run directories found in:", maxent_results_dir)
}

cat("Found", length(run_dirs), "model runs\n")

# Initialize results list
model_predictions <- list()
model_metadata <- list()

# Extract predictions and metadata from each run
for (run_dir in run_dirs) {

  run_name <- basename(run_dir)
  cat("Processing:", run_name, "\n")

  # Find ASC prediction files
  asc_files <- list.files(run_dir, pattern = "\\.asc$", full.names = TRUE, ignore.case = TRUE)
  asc_files <- asc_files[!grepl("species|clamp|mess", tolower(asc_files))]

  # Find maxent model file
  model_files <- list.files(run_dir, pattern = "\\.lambdas$", full.names = TRUE)

  if (length(asc_files) > 0 && length(model_files) > 0) {

    # Load prediction raster
    pred_raster <- raster::raster(asc_files[1])

    # Load model results
    results_file <- file.path(run_dir, "maxentResults.csv")

    if (file.exists(results_file)) {
      results <- read.csv(results_file, stringsAsFactors = FALSE)

      # Extract key metrics
      train_auc <- NA
      test_auc <- NA
      n_params <- NA

      train_cols <- grep("Training.AUC|AUC.training", names(results), ignore.case = TRUE, value = TRUE)
      test_cols <- grep("Test.AUC|AUC.test", names(results), ignore.case = TRUE, value = TRUE)
      param_cols <- grep("X.Parameters|parameters", names(results), ignore.case = TRUE, value = TRUE)

      if (length(train_cols) > 0) train_auc <- as.numeric(results[[train_cols[1]]][1])
      if (length(test_cols) > 0) test_auc <- as.numeric(results[[test_cols[1]]][1])
      if (length(param_cols) > 0) n_params <- as.numeric(results[[param_cols[1]]][1])

      # --- Replace your current parsing block with this one ---
      parts <- strsplit(run_name, "_")[[1]]

      if (grepl("regmult", parts[1], ignore.case = TRUE)) {

        # handle decimals like regmult_0_5_LPH
        if (length(parts) >= 3 && grepl("^[0-9]+$", parts[2]) && grepl("^[0-9]+$", parts[3])) {
          reg_mult <- as.numeric(paste0(parts[2], ".", parts[3]))
          features <- paste(parts[4:length(parts)], collapse = "_")

          # handle integers like regmult_1_LPH
        } else if (length(parts) >= 2 && grepl("^[0-9]+$", parts[2])) {
          reg_mult <- as.numeric(parts[2])
          features <- paste(parts[3:length(parts)], collapse = "_")

          # fallback when regex fails (no valid number)
        } else {
          reg_mult <- NA
          features <- paste(parts[-1], collapse = "_")
        }

      } else {
        reg_mult <- NA
        features <- paste(parts, collapse = "_")
      }

      model_predictions[[run_name]] <- pred_raster
      model_metadata[[run_name]] <- list(
        reg_mult = reg_mult,
        features = features,
        train_auc = train_auc,
        test_auc = test_auc,
        n_params = n_params,
        file = asc_files[1]
      )

      cat("  Reg Mult:", reg_mult, "| Features:", features,
          "| Train AUC:", round(train_auc, 3), "| Test AUC:", round(test_auc, 3), "\n")
    }
  }
}

cat("Successfully loaded", length(model_predictions), "models\n")

# =============================================================================
# BUILD COMPREHENSIVE DATA TABLE WITH PREDICTIONS
# =============================================================================

cat("\n=== Building Data Table with Model Predictions ===\n")

# Combine presence and absence points
all_points <- rbind(
  data.frame(presence_points, presence.absence = 1),
  data.frame(absence_points, presence.absence = 0)
)

# Extract habitat suitability values at all points
extraction_data <- data.frame(
  sample.No = 1:nrow(all_points),
  presence.absence = all_points$presence.absence,
  longitude = all_points[, 1],
  latitude = all_points[, 2],
  stringsAsFactors = FALSE
)

cat("Total evaluation points:", nrow(extraction_data),
    "(", sum(extraction_data$presence.absence == 1), "presence,",
    sum(extraction_data$presence.absence == 0), "absence)\n")

# Extract values from each model prediction
for (model_name in names(model_predictions)) {

  pred_raster <- model_predictions[[model_name]]

  # Create point matrix for extraction
  points_matrix <- as.matrix(extraction_data[, c("longitude", "latitude")])
  extracted_values <- raster::extract(pred_raster, points_matrix)

  # Clean column name for habitat suitability
  col_name <- paste0("habitat_suitability_", model_name)
  extraction_data[[col_name]] <- extracted_values
}

# Remove rows with any NA values
extraction_data_clean <- extraction_data[complete.cases(extraction_data), ]

cat("Data table with", nrow(extraction_data_clean), "complete presence records\n")
cat("Habitat suitability predictions from", length(model_predictions), "models\n")

# Save to CSV
output_csv <- file.path(output_dir, "model_predictions_comprehensive.csv")
write.csv(extraction_data_clean, output_csv, row.names = FALSE)
cat("Saved to:", output_csv, "\n")

# =============================================================================
# AICc COMPARISON OF HABITAT SUITABILITY MODELS
# =============================================================================

cat("\n=== AICc Model Comparison ===\n")

# Load the comprehensive prediction data
model_data <- read.csv(output_csv, stringsAsFactors = FALSE)

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

# Identify best model
best_idx <- which.min(aicc_results$AICc)
best_model_full <- aicc_results[best_idx, ]
cat("Best Model (lowest AICc):", best_model_full$Model, "\n")
cat("Regularization Multiplier:", best_model_full$reg_mult, "\n")
cat("Feature Classes:", best_model_full$features, "\n")
cat("AICc =", round(best_model_full$AICc, 2), "\n")
cat("Akaike Weight =", round(best_model_full$Weight, 4), "\n\n")

# Define default model (typically regmult_1_LQH or similar)
# Adjust this based on your MaxEnt untuned settings
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

cat("Default Model (MaxEnt untuned settings):", default_model$Model, "\n")
cat("AICc =", round(default_model$AICc, 2), "\n\n")

# Models within 2 Delta AICc (substantial support)
substantial <- aicc_results[aicc_results$Delta_AICc <= 2, ]
if (nrow(substantial) > 1) {
  cat("Models with substantial support (Delta AICc <= 2):\n")
  print(substantial[, c("Model", "reg_mult", "features", "AICc", "Delta_AICc", "Weight")], digits = 4)
  cat("\n")
}

# Save results
output_aicc_csv <- file.path(output_dir, "AICc_comparison_results.csv")
write.csv(aicc_results, output_aicc_csv, row.names = FALSE)
cat("AICc results saved to:", output_aicc_csv, "\n")

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

# Create line plot
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
    x = "Regularization Multiplier",
    y = "AICc (lower is better)",
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
output_plot <- file.path(output_dir, "AICc_comparison_plot.png")
ggsave(output_plot, fine_tuning_plot, width = 14, height = 7, dpi = 300)
cat("Plot saved to:", output_plot, "\n")

toc()
cat("\n=== AICc Analysis Complete ===\n")
cat("You can now run the optional fine-tuning analysis with: source('path_to_optional_script.R')\n")
