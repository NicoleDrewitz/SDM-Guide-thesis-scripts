# ==============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Script to fill in missing Snow water equivalent values (swe)
    #  swe = 0, where swe= no data AND snow cover days =0
    # this is an imperfect approximation, as snow cover days = 0,
        #is 0 annual days estimated with >= 1kg/m^2 of snow

# Author: Nicole Drewitz with Claude AI
# Last updated: 2025-12-1
# ==============================================================================

library(terra)
library(tictoc)
library(here)

tic() # record runtime

# Define input and output directories
input_folder <- here::here("Environmental_Variables/Processed continuous variables")
output_folder <- here::here("Environmental_Variables/Processed continuous variables/corrected")

# Create output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Find all SWE .asc files in the input folder
swe_files <- list.files(input_folder, pattern = "swe_.*\\.asc$", full.names = TRUE)

if (length(swe_files) == 0) {
  stop("No SWE .asc files found in the input folder.")
}

cat("Found", length(swe_files), "SWE files to process\n\n")

# Process each SWE file
for (swe_file in swe_files) {
  
  # Extract the base filename
  swe_basename <- basename(swe_file)
  
  # Create corresponding SCD filename by replacing "swe" with "scd"
  scd_basename <- sub("swe_", "scd_", swe_basename)
  scd_file <- file.path(input_folder, scd_basename)
  
  cat("Processing:", swe_basename, "\n")
  
  # Check if corresponding SCD file exists
  if (!file.exists(scd_file)) {
    cat("  WARNING: Corresponding SCD file not found:", scd_basename, "\n")
    cat("  Skipping this file.\n\n")
    next
  }
  
  # Load rasters
  swe <- rast(swe_file)
  scd <- rast(scd_file)
  
  # Diagnostic: Check for NA values and specific values
  n_swe_na <- global(is.na(swe), "sum", na.rm = TRUE)[1,1]
  n_scd_0 <- global(scd == 0, "sum", na.rm = TRUE)[1,1]
  n_both <- global(is.na(swe) & scd == 0, "sum", na.rm = TRUE)[1,1]
  
  cat("  SWE pixels with NA (NoData):", n_swe_na, "\n")
  cat("  SCD pixels with value 0:", n_scd_0, "\n")
  cat("  Pixels meeting both conditions (NA & 0):", n_both, "\n")
  
  # Count pixels that will be changed (before correction)
  n_changed <- n_both
  
  # Perform the conditional replacement
  # Where SWE = NA (9999 NoData) AND SCD = 0, replace SWE with SCD (which is 0)
  swe_corrected <- ifel(is.na(swe) & scd == 0, scd, swe)
  
  # Create output filename (keep same name as input)
  output_file <- file.path(output_folder, swe_basename)
  
  # Save the corrected raster
  writeRaster(swe_corrected, output_file, overwrite = TRUE, NAflag=9999)
  
  cat("  Corrected file saved to:", output_file, "\n")
  cat("  Number of pixels changed:", n_changed, "\n\n")
}

toc() # end runtime recording

cat("All files processed successfully!\n")