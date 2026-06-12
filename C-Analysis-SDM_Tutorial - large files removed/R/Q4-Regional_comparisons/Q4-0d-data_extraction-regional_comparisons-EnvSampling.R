# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q4:
# R Script to Extract Environmental Values from ASC Files for Multiple Climate Scenarios
# AT RANDOM SAMPLING LOCATIONS - WITH REGION COLUMN

# Author: Nicole Drewitz
# Date: 2026-02-25 with Perplexity AI
# =============================================================================

# Load required libraries
required_packages <- c("here", "terra", "dplyr", "readr", "tictoc", "tools")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

tic()

# Define paths relative to RProject root
sampling_locations_dir <- here::here("Data/Q4-random_sampling_locations")
output_dir <- here::here("Data/Q4-random_sampling_locations-with_values")
crs_laea <- "EPSG:3035"

# Define climate scenario folders (matching your structure)
climate_folders <- list(
  current = "Environmental Variables - Final - 1 - Current",
  ssp126 = "Environmental Variables - Final - 2 - ssp126",
  ssp370 = "Environmental Variables - Final - 3 - ssp370",
  ssp585 = "Environmental Variables - Final - 4 - ssp585"
)

# ASC files directory for environmental variables
asc_dir <- #here::here("Input_files/Environmental_variables")
  "/Users/nicoledrewitz/Library/CloudStorage/OneDrive-Menntaský/Thesis Research Files - Mac Connection/A - Thesis Research Files - Mac Connection/1 - GIS files/c - Environmental variables - 2025-08-21/A_1_c- Environmental variables - modelled area/A_1_c - Environmental Variables - Final Model - 2026-01-15"

# Define variables that should ONLY be extracted from current climate folder
current_only_variables <- c(
  "soil_phh2o_0-5-30cm_mean_averaged",
  "soil_soc_0-5-30cm_mean_averaged",
  "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged",
  "Topo_northness"
)
cat("Current-only variables:", paste(current_only_variables, collapse = ", "), "\n")

# Function to check if a variable should only be extracted from current climate
is_current_only <- function(filename) {
  any(sapply(current_only_variables, function(pattern) {
    grepl(pattern, filename, ignore.case = TRUE)
  }))
}


# Get all sampling location CSV files
location_files <- list.files(sampling_locations_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
cat("Found", length(location_files), "sampling location files\n")

# Process each sampling location file separately
for (points_file in location_files) {
  cat(paste0("\n=== Processing: ", basename(points_file), " ===\n"))

  # Read sampling locations
  sampling_pts <- read_csv(points_file, show_col_types = FALSE)

  # **Per-row region handling**
  region_col <- names(sampling_pts)[1]
  cat("  Region column ('", region_col, "'): ", length(unique(sampling_pts[[region_col]])), " unique values\\n")
  sampling_pts <- sampling_pts %>%
    rename(region = !!sym(region_col))

  if (!all(c("x", "y") %in% names(sampling_pts))) {
    stop("File must contain 'x' and 'y' columns: ", basename(points_file))
  }

  # Convert to SpatVector (region preserved in data frame)
  sp_pts <- vect(sampling_pts, geom = c("x", "y"), crs = crs_laea)

  # Initialize results list for this sampling file
  results_list <- list()

  # Loop through each climate scenario
  for (scenario_name in names(climate_folders)) {
    cat(paste0("Processing climate scenario: ", scenario_name, "\n"))

    scenario_folder <- file.path(asc_dir, climate_folders[[scenario_name]])

    if (!dir.exists(scenario_folder)) {
      warning(paste0("Folder not found: ", scenario_folder))
      next
    }

    asc_files <- list.files(scenario_folder, pattern = "\\.asc$", full.names = TRUE, ignore.case = TRUE)
    if (length(asc_files) == 0) {
      warning(paste0("No ASC files found in: ", scenario_folder))
      next
    }

    cat(paste0("  Found ", length(asc_files), " ASC files\n"))

    # Create scenario data WITH region column
    scenario_data <- data.frame(
      region = sampling_pts$region,  # ← ADDED REGION COLUMN
      x = sampling_pts$x,
      y = sampling_pts$y,
      climateScenario = scenario_name
    )

    # Extract values from each ASC file
    for (asc_file in asc_files) {
      filename <- basename(asc_file)
      variable_name <- tools::file_path_sans_ext(filename)

      # Skip current-only variables for future scenarios
      if (scenario_name != "current" && is_current_only(filename)) {
        cat(paste0("    Skipping '", filename, "' (current climate only)\n"))
        next
      }

      cat(paste0("    Extracting: ", filename, "\n"))

      tryCatch({
        r <- rast(asc_file)
        crs(r) <- crs_laea
        extracted_values <- terra::extract(r, sp_pts)[, 2]
        scenario_data[[variable_name]] <- extracted_values
      }, error = function(e) {
        warning(paste0("Error processing ", filename, ": ", e$message))
        scenario_data[[variable_name]] <- NA
      })
    }

    results_list[[scenario_name]] <- scenario_data
  }

  # Combine all scenarios for this sampling file
  cat("Combining scenarios...\n")
  final_results <- bind_rows(results_list)

  # Reorder columns: region, x, y, climateScenario first
  first_cols <- c("region", "x", "y", "climateScenario")  # ← REGION FIRST
  other_cols <- setdiff(names(final_results), first_cols)
  final_results <- final_results[, c(first_cols, other_cols)]

  # Create output filename
  base_name <- tools::file_path_sans_ext(basename(points_file))
  out_file <- file.path(output_dir, paste0("Q4-0d-", base_name, "_EnvVariables.csv"))

  # Save
  write_csv(final_results, out_file)
  cat(paste0("Saved: ", basename(out_file), "\n"))
  cat(paste0("  Rows: ", nrow(final_results), ", Columns: ", ncol(final_results), "\n"))
}

cat("\n=== ALL FILES PROCESSED ===\n")
toc()