# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q3

# R Script to Extract Values from ASC Files for Multiple Climate Scenarios
# Extracted from presence location reservced for model analysis
# Author: Nicole Drewitz with claude.ai
# Date: 2026-02-16
# =============================================================================

# Load required libraries
required_packages <- c("raster", "dplyr", "tidyr", "tictoc", "here")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

tic()

# Define base paths
base_path <- here::here("Input_files/Environmental_variables") 
#base_path <- "/Users/nicoledrewitz/Library/CloudStorage/OneDrive-Menntaský/Thesis Research Files - Mac Connection/A - Thesis Research Files - Mac Connection/1 - GIS files/c - Environmental variables - 2025-08-21/A_1_c- Environmental variables - modelled area/A_1_c - Environmental Variables - Final Model - 2026-01-15"

locations_path <- file.path(here::here("Input_files/V2_cloudberry_presence_analysis.csv"))

output_path <- here::here("Data/Q3-extracted_environmental_values-analysis_locations.csv")

# Define climate scenario folders
climate_folders <- list(
  current = "Environmental Variables - Final - 1 - Current",
  ssp126 = "Environmental Variables - Final - 2 - ssp126",
  ssp370 = "Environmental Variables - Final - 3 - ssp370",
  ssp585 = "Environmental Variables - Final - 4 - ssp585"
)

# Define variables that should ONLY be extracted from current climate folder
# (e.g., topographic variables, soil variables, or other static variables)
current_only_variables <- c(
  "soil_phh2o_0-5-30cm_mean_averaged",
  "soil_soc_0-5-30cm_mean_averaged",
  "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged",
  "Topo_northness"
  # Add other variable names or patterns that should only come from current climate
)
current_only_variables

# Read location data
cat("Reading location data...\n")
locations <- read.csv(locations_path)

# Check required columns
if (!all(c("longitude", "latitude") %in% names(locations))) {
  stop("Location file must contain 'longitude' and 'latitude' columns")
}

# Add sample.No if not present
if (!"sample.No" %in% names(locations)) {
  locations$sample.No <- 1:nrow(locations)
}

# Initialize results dataframe
results_list <- list()

# Function to check if a variable should only be extracted from current climate
is_current_only <- function(filename) {
  # Check if filename contains any of the current_only_variables patterns
  any(sapply(current_only_variables, function(pattern) {
    grepl(pattern, filename, ignore.case = TRUE)
  }))
}

# Loop through each climate scenario
for (scenario_name in names(climate_folders)) {
  cat(paste0("\nProcessing climate scenario: ", scenario_name, "\n"))

  folder_path <- file.path(base_path, climate_folders[[scenario_name]])

  # Check if folder exists
  if (!dir.exists(folder_path)) {
    warning(paste0("Folder not found: ", folder_path))
    next
  }

  # Get all ASC files in the folder
  asc_files <- list.files(folder_path, pattern = "\\.asc$", full.names = TRUE, ignore.case = TRUE)

  if (length(asc_files) == 0) {
    warning(paste0("No ASC files found in: ", folder_path))
    next
  }

  cat(paste0("Found ", length(asc_files), " ASC files\n"))

  # Create a temporary dataframe for this scenario
  scenario_data <- data.frame(
    sample.No = locations$sample.No,
    longitude = locations$longitude,
    latitude = locations$latitude,
    climateScenario = scenario_name
  )

  # Extract values from each ASC file
  for (asc_file in asc_files) {
    filename <- basename(asc_file)
    variable_name <- tools::file_path_sans_ext(filename)

    # Skip current-only variables for future scenarios
    if (scenario_name != "current" && is_current_only(filename)) {
      cat(paste0("  Skipping '", filename, "' (current climate only)\n"))
      next
    }

    cat(paste0("  Extracting: ", filename, "\n"))

    tryCatch({
      # Read raster
      r <- raster::raster(asc_file)

      # Extract values at locations
      extracted_values <- raster::extract(r, locations[, c("longitude", "latitude")])

      # Add to scenario data
      scenario_data[[variable_name]] <- extracted_values

    }, error = function(e) {
      warning(paste0("Error processing ", filename, ": ", e$message))
      scenario_data[[variable_name]] <- NA
    })
  }

  # Add to results list
  results_list[[scenario_name]] <- scenario_data
}

# Combine all scenarios into one dataframe
cat("\nCombining all scenarios...\n")
final_results <- bind_rows(results_list)

# Reorder columns to have sample.No, longitude, latitude, climateScenario first
first_cols <- c("sample.No", "longitude", "latitude", "climateScenario")
other_cols <- setdiff(names(final_results), first_cols)
final_results <- final_results[, c(first_cols, other_cols)]

# Save to CSV
cat(paste0("\nSaving results to: ", output_path, "\n"))
write.csv(final_results, output_path, row.names = FALSE)

# Print summary
cat("\n=== EXTRACTION COMPLETE ===\n")
cat(paste0("Total rows: ", nrow(final_results), "\n"))
cat(paste0("Total columns: ", ncol(final_results), "\n"))
cat(paste0("Climate scenarios: ", paste(unique(final_results$climateScenario), collapse = ", "), "\n"))
cat(paste0("Output file: ", output_path, "\n"))

toc()
# Display first few rows
cat("\nFirst few rows:\n")
print(head(final_results))
