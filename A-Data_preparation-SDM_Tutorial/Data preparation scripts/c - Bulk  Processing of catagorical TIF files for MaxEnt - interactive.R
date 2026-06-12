# ==============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Bulk processing of TIF files for MaxEnt
# nearest neighbour resampling
# Created on 18-09-2025 with Claude.ai by Nicole Drewitz
# Working for categorical variables

# 1. Extract desired files (i.e. environmental variables)
# 2. Reproject (WG84 to 3035 Europe if not already in 3035)
# 3. Clip by study area
# 4. Save as ASC file with 1000m^2 resolution and nodata = 9999

# USER GUIDE:
# 1. choose input and output folder
# 2. choose between option 1, 2, or 3
# 3. answer script prompts
# 4. optional - paste in terminal window to keep r running (caffeinate -i) or (caffeinate -i -t 7200) to set a timer for 7200 seconds = 2 hours
          # If files are too large to proccess, change order of code to clip by study area first, then reproject
# ==============================================================================
# Load necessary libraries ====
library(terra)
library(sf)
library(here)

# 1. Define categorical variable directories ====
input_folder <- here::here("Input files/Other categorical variables download")
output_folder <- here::here("Environmental_Variables/Downloads/Processed Categorical variables")

# Create output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

extract_dir <- file.path(output_folder, "zip_extraction")
if (!dir.exists(extract_dir)) {
  dir.create(extract_dir, recursive = TRUE)
}

# Load clipping vector ====
clip_vector <- st_read(here::here("Input files/Modelling_area-regional_divisions.gpkg"))

# Convert to terra vector
clip_vect <- vect(clip_vector)

# 1. Extract desired files (i.e. environmental variables)
file_list <- list.files(input_folder, pattern = "\\.(tif|tiff|geotiff)$", full.names = TRUE, ignore.case = TRUE)

# Global variables to store user preferences
global_nodata_option <- NULL
global_additional_nodata <- NULL

#________________________________________________________________________________
# MANUAL SETUP: Choose your option by uncommenting ONE of the following sections: ====

# OPTION 1: No additional nodata values for any file
# global_nodata_option <- "none"
# global_additional_nodata <- NULL

# OPTION 2: Same additional nodata value(s) for all files  
# global_nodata_option <- "same_for_all"
# global_additional_nodata <- c(-9999, 5)  # Replace with your actual nodata values

# OPTION 3: Prompt for each file individually
global_nodata_option <- "individual"
global_additional_nodata <- NULL
#________________________________________________________________________________

# Function to get nodata handling preferences (now just displays current settings) ====
get_nodata_preferences <- function() {
  cat("\n=== NODATA HANDLING SETUP ===\n")
  cat("You have", length(file_list), "files to process.\n")
  cat("Standard nodata value will be set to 9999 for all files.\n\n")
  
  if (global_nodata_option == "none") {
    cat("??? Current setting: No additional nodata values for any file\n")
  } else if (global_nodata_option == "same_for_all") {
    cat("??? Current setting: Same additional nodata values for all files\n")
    if (!is.null(global_additional_nodata)) {
      cat("  Additional nodata values:", paste(global_additional_nodata, collapse = ", "), "\n")
    }
  } else {
    cat("??? Current setting: Individual prompts for each file\n")
  }
  
  cat("===============================\n")
}

# Function to get additional nodata values based on user preference ====
get_additional_nodata <- function(filename) {
  cat("\nProcessing file:", basename(filename), "\n")
  
  if (global_nodata_option == "none") {
    return(NULL)
    
  } else if (global_nodata_option == "same_for_all") {
    if (!is.null(global_additional_nodata) && length(global_additional_nodata) > 0) {
      cat("Using global additional nodata values:", paste(global_additional_nodata, collapse = ", "), "\n")
    }
    return(global_additional_nodata)
    
  } else { # individual prompting
    cat("Standard nodata value will be set to 9999\n")
    
    # Ask user if there are additional nodata values
    response <- readline(prompt = "Are there any additional integer codes that should be treated as nodata? (y/n): ")
    
    additional_nodata <- NULL
    if (tolower(response) %in% c("y", "yes")) {
      nodata_input <- readline(prompt = "Enter the integer code(s) for additional nodata values (separate multiple values with commas): ")
      
      # Parse the input
      if (nzchar(nodata_input)) {
        additional_nodata <- as.numeric(unlist(strsplit(nodata_input, ",")))
        # Remove any NA values that might result from invalid input
        additional_nodata <- additional_nodata[!is.na(additional_nodata)]
        
        if (length(additional_nodata) > 0) {
          cat("Additional nodata values:", paste(additional_nodata, collapse = ", "), "\n")
        }
      }
    }
    
    return(additional_nodata)
  }
}

# 2. Reproject function (WGS84 to EPSG 3035 Europe if not already 3035) ====
reproject_raster <- function(r) {
  target_crs <- "EPSG:3035"
  if (crs(r) != target_crs) {
    # Use "near" method for categorical data to preserve integer values
    r <- project(r, target_crs, method = "near")
  }
  return(r)
}

# Get user preferences for nodata handling before processing files
get_nodata_preferences()

# Process each file ====
for (file in file_list) {
  # Load raster
  r <- rast(file)
  
  # Get additional nodata values from user
  additional_nodata <- get_additional_nodata(file)
  
  # 2. Reproject raster if needed
  r <- reproject_raster(r)
  
  # 3. Clip by study area (crop + mask)
  r_cropped <- crop(r, ext(clip_vect))
  r_masked <- mask(r_cropped, clip_vect)
  
  # 4. Save as ASC file with 1000m x 1000m resolution and nodata=9999
  # FIXED: Create template with exact 1000m resolution and proper grid alignment
  
  # Get the extent from the vector (clipping vector)
  shapefile_extent <- ext(clip_vect)
  
  # Create template raster with exactly 1000m x 1000m resolution
  # Round extent to align with 1000m grid
  xmin_aligned <- floor(shapefile_extent[1] / 1000) * 1000
  ymin_aligned <- floor(shapefile_extent[3] / 1000) * 1000
  xmax_aligned <- ceiling(shapefile_extent[2] / 1000) * 1000
  ymax_aligned <- ceiling(shapefile_extent[4] / 1000) * 1000
  
  # Create extent with aligned coordinates
  aligned_extent <- ext(xmin_aligned, xmax_aligned, ymin_aligned, ymax_aligned)
  
  # Create template raster with exact 1000m resolution
  r_template <- rast(aligned_extent, resolution = c(1000, 1000), crs = crs(clip_vect))
  
  # Project the input raster to the template grid
  r_projected <- project(r, r_template, method = "near") #for categorical variables or ranked variables (possible unequal rank)
  
  # Now mask with the vector
  r_resampled <- mask(r_projected, clip_vect)
  
  # Handle nodata values
  vals <- values(r_resampled)
  
  # Convert additional nodata values to NA first
  if (!is.null(additional_nodata) && length(additional_nodata) > 0) {
    for (nodata_val in additional_nodata) {
      vals[vals == nodata_val & !is.na(vals)] <- NA
      cat("Converted", sum(vals == nodata_val & !is.na(vals), na.rm = TRUE), "pixels with value", nodata_val, "to nodata\n")
    }
  }
  
  # Set NA value to 9999
  NAflag(r_resampled) <- 9999
  vals[is.na(vals)] <- 9999
  values(r_resampled) <- vals
  
  # Output filename with ASC extension
  output_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(file)), ".asc"))
  
  # Write raster as ASCII file with nodata flag 9999 ====
  writeRaster(r_resampled, filename = output_file, filetype = "AAIGrid", overwrite = TRUE, NAflag = 9999)
  
  cat("Saved:", basename(output_file), "\n")
  cat("Resolution:", res(r_resampled), "meters\n")
  cat("Extent:", as.vector(ext(r_resampled)), "\n")
  cat("-----------------------------------\n")
}
print("Processing complete")
