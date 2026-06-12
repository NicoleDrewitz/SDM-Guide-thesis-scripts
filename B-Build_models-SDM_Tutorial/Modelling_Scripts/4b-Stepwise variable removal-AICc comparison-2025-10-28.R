# =============================================================================
# Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Selects best model with AICc by sampling habitat suitability from presence/absence locations (or fruiting/non-fruiting locations)
# for selecting best model from stepwise variable removal models
# Create by Nicole Drewitz with Claude.ai October 4, 2025
# last updated: 28 October

# best model with 10 variables from regmult14.5_threshold
# =============================================================================

#_____________________________________________________________________________
# USER GUIDE
#  Data required:
#   Folders of MaxEnt models to compare (including ASC of habitat suitability)
#   Species presence locations (CSV) - not used in MaxEnt models
#   Species absence locations (if selecting best tuning parameters from presence/absence AICc ranking)
# Used to select best tuning parameters
#   For stepwise variable removal using the selected tuning parameters afterwards,
#       it may be better to create models manually in MaxEnt to have access to the latest software version
# line 55 update folder name (best settings used in stepwise variable removal)

# =============================================================================
# INSTALL AND LOAD REQUIRED LIBRARIES
# =============================================================================

required_packages <- c("dismo", "terra", "ggplot2", "dplyr", "sf", "gridExtra", "rJava", "raster", "here")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# File path
base_dir <- here::here("Output_files")
setwd(base_dir)

# Presence and absence data files reserved for model selection
presence_file_AICc <- here::here("Input_files/V2_cloudberry_presence_select.csv")
absence_file_AICc <- here::here("Input_files/V2_cloudberry_absence_select.csv")

# stepwise variable removal path
output_dir <- file.path(base_dir, "MaxEnt_Results_Stepwise_Variable_removal")

file_dir <- file.path(output_dir, "Model_Comparison_AICc")
dir.create(file_dir, recursive = TRUE, showWarnings = TRUE)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

cat("=== Loading Species Occurrence Data ===\n")

# Load and prepare presence data
presence_data_AICc <- read.csv(presence_file_AICc, stringsAsFactors = FALSE)
presence_points <- na.omit(presence_data_AICc[, c("longitude", "latitude")])

# Load and prepare absence data
absence_data <- read.csv(absence_file_AICc, stringsAsFactors = FALSE)
absence_points <- na.omit(absence_data[, c("longitude", "latitude")])

cat("Presence points for AICc:", nrow(presence_points), "\n")
cat("Absence points for AICc:", nrow(absence_points), "\n")

# =============================================================================
# HABITAT SUITABILITY SAMPLING - SINGLE COMBINED CSV
# =============================================================================

cat("\n=== Sampling Habitat Suitability at Evaluation Locations ===\n")

# Initialize base dataframe with location info
combined_results <- data.frame(
  sample.No = 1:(nrow(presence_points) + nrow(absence_points)),
  presence.absence = c(rep(1, nrow(presence_points)), rep(0, nrow(absence_points))),
  longitude = c(presence_points$longitude, absence_points$longitude),
  latitude = c(presence_points$latitude, absence_points$latitude)
)

# Combine all location points for extraction
all_points <- rbind(presence_points, absence_points)

# Find all model directories with RUN prefix
model_dirs <- list.dirs(output_dir, recursive = FALSE)
model_dirs <- model_dirs[grepl("^RUN", basename(model_dirs))]

if (length(model_dirs) == 0) {
  stop("No model directories with 'RUN' prefix found for habitat suitability sampling.\n")
}

cat("Found", length(model_dirs), "model directories to process.\n")

# =============================================================================
# HABITAT SUITABILITY SAMPLING - SINGLE COMBINED CSV
# =============================================================================

cat("\n=== Sampling Habitat Suitability at Evaluation Locations ===\n")

# Initialize base dataframe with location info
combined_results <- data.frame(
  sample.No = 1:(nrow(presence_points) + nrow(absence_points)),
  presence.absence = c(rep(1, nrow(presence_points)), rep(0, nrow(absence_points))),
  longitude = c(presence_points$longitude, absence_points$longitude),
  latitude = c(presence_points$latitude, absence_points$latitude)
)

# Combine all location points for extraction
all_points <- rbind(presence_points, absence_points)

# =============================================================================
# PROCESS MODELS
# =============================================================================

cat("\n=== Processing Models ===\n")

# Process each model directory
for (model_dir in model_dirs) {
  model_name <- basename(model_dir)
  cat("Processing model:", model_name, "\n")
  
  # Look for ASC prediction files
  asc_files <- list.files(model_dir, pattern = "\\.asc$", full.names = TRUE)
  
  if (length(asc_files) == 0) {
    cat("  No ASC files found in", model_name, "\n")
    next
  }
  
  for (asc_file in asc_files) {
    tryCatch({
      # Load the habitat suitability raster
      habitat_raster <- rast(asc_file)
      
      # Extract habitat suitability values at all locations
      habitat_values <- terra::extract(habitat_raster, all_points)
      
      # Create column name based on folder name
      column_name <- model_name
      
      # Get the actual raster layer name (first column after ID)
      raster_col_name <- names(habitat_values)[2]
      
      # Add to combined results using the folder name as column name
      combined_results[[column_name]] <- habitat_values[[raster_col_name]]
      
      cat("  Added column:", column_name, "\n")
      
    }, error = function(e) {
      cat("  Error processing", basename(asc_file), ":", conditionMessage(e), "\n")
    })
  }
}

# =============================================================================
# SAVE COMBINED RESULTS
# =============================================================================

cat("\n=== Saving Combined Results ===\n")

# Save the combined results to CSV
output_csv <- file.path(file_dir, "habitat_suitability_all_stepwise_runs.csv")
write.csv(combined_results, output_csv, row.names = FALSE)

cat("Combined results saved to:", output_csv, "\n")
cat("Total columns:", ncol(combined_results), "\n")
cat("Total rows:", nrow(combined_results), "\n")
#print(names(combined_results))
#print(head(habitat_values))
#print(names(habitat_values))

# =============================================================================
# CALCULATE AICc FOR ALL MODELS
# =============================================================================

cat("\n=== Calculating AICc for All Models ===\n")

# Function to calculate AICc from a logistic regression model
calc_AICc <- function(model) {
  n <- length(model$fitted.values)
  k <- length(coef(model))
  AIC_val <- AIC(model)
  AICc_val <- AIC_val + (2 * k * (k + 1)) / (n - k - 1)
  return(AICc_val)
}

# Initialize results storage
aicc_results <- data.frame(
  model_name = character(),
  AICc = numeric(),
  stringsAsFactors = FALSE
)

# Calculate AICc for each habitat suitability column
habitat_columns <- names(combined_results)[5:ncol(combined_results)]

for (col_name in habitat_columns) {
  # Remove rows with NA values for this model
  model_data <- combined_results[!is.na(combined_results[[col_name]]), c("presence.absence", col_name)]
  
  if (nrow(model_data) > 0) {
    tryCatch({
      # Fit logistic regression model
      formula_str <- paste("presence.absence ~", col_name)
      glm_model <- glm(as.formula(formula_str), data = model_data, family = binomial())
      
      # Calculate AICc
      aicc <- calc_AICc(glm_model)
      
      # Store results
      aicc_results <- rbind(aicc_results, data.frame(model_name = col_name, AICc = aicc))
      
    }, error = function(e) {
      cat("  Error calculating AICc for", col_name, ":", conditionMessage(e), "\n")
    })
  }
}

# Order by lowest AICc
aicc_results <- aicc_results %>% arrange(AICc)

cat("\n=== AICc Comparison Results ===\n")
print(aicc_results)

cat("\n=== Best Model (Lowest AICc) ===\n")
best_model <- aicc_results[1, ]
print(best_model)

# =============================================================================
# PREPARE DATA FOR PLOTTING
# =============================================================================

cat("\n=== Preparing Plot Data ===\n")

# Print model names
cat("Model names found:\n")
print(aicc_results$model_name)

# Load stepwise removal log to get number of variables
stepwise_log_path <- file.path(output_dir, "stepwise_removal_log.csv")

if (file.exists(stepwise_log_path)) {
  stepwise_log <- read.csv(stepwise_log_path, stringsAsFactors = FALSE)
  
  cat("\nStepwise removal log loaded:\n")
  print(stepwise_log)
  
  # Merge with aicc_results based on model name (iteration = RUN)
  # Convert iteration to match RUN format (e.g., 1 becomes "RUN01")
  aicc_results <- aicc_results %>%
    left_join(stepwise_log %>% 
                select(iteration, variables_remaining) %>%
                mutate(
                  model_name = paste0("RUN", sprintf("%02d", iteration)),
                  number_of_variables = variables_remaining + 1
                ) %>%
                select(model_name, number_of_variables),
              by = "model_name")
  
  cat("\nAICc results with number of variables:\n")
  print(aicc_results[, c("model_name", "number_of_variables", "AICc")])
  
} else {
  cat("\nWarning: stepwise_removal_log.csv not found at:", stepwise_log_path, "\n")
  cat("Creating plot with run numbers instead.\n")
  aicc_results <- aicc_results %>%
    mutate(number_of_variables = 1:n())
}

# Identify best model
best_model_full <- aicc_results[1, ]

cat("\nBest model identified:\n")
print(best_model_full)

# =============================================================================
# CREATE PLOT
# =============================================================================

cat("=== Creating Stepwise Variable Removal Plot ===\n")

stepwise_plot <- ggplot(aicc_results, aes(x = number_of_variables, y = AICc)) +
  geom_line(linewidth = 1, color = "#440154FF") +
  geom_point(size = 3, color = "#440154FF") + 
  geom_point(data = best_model_full, aes(x = number_of_variables, y = AICc), 
             color = "red", size = 4, shape = 18) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Stepwise Variable Removal: AICc vs Number of Variables",
    x = "Number of Variables",
    y = "AICc (lower is better)"
  ) +
  scale_x_continuous(breaks = seq(min(aicc_results$number_of_variables), 
                                  max(aicc_results$number_of_variables), 
                                  by = 1)) +
  # Best model line
  geom_hline(yintercept = best_model_full$AICc, linetype = "dashed", 
             color = "red", linewidth = 0.8) +
  # Best model annotation
  annotate("text", 
           x = max(aicc_results$number_of_variables), 
           y = best_model_full$AICc, 
           label = paste0("Best model: ", best_model_full$model_name, 
                          " (", best_model_full$number_of_variables, " variables)"),
           hjust = 1, vjust = -0.5, color = "red", fontface = "bold", size = 3.5)

stepwise_plot
ggsave(file.path(file_dir, "AICc_comparison_plot.png"), stepwise_plot, width = 8, height = 6, dpi = 300)
# Display plot

# =============================================================================
# PLOT WITH VARIABLES OF BEST MODEL
# =============================================================================

library(xml2)

best_model_dir <- file.path(output_dir, best_model_full$model_name)

html_path <- file.path(best_model_dir, "maxent.html")
html_content <- read_html(html_path)
# This will vary, but find the table with percent contributions
percent_table <- xml_find_first(html_content, "//table[.//th[contains(text(),'Percent contribution')]]")
percent_rows <- xml_find_all(percent_table, ".//tr[position()>1]")
variables <- xml_text(xml_find_all(percent_rows, ".//td[1]"))
percentages <- as.numeric(xml_text(xml_find_all(percent_rows, ".//td[2]")))
percent_df <- data.frame(Variable = variables, Percent = percentages)

