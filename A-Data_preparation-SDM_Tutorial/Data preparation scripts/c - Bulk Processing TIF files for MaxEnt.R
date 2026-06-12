# ==============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Bulk processing of TIF files for MaxEnt
# Created on 17-09-2025 with Claude.ai by Nicole Drewitz
# Nearest neighbour resampling used for reprojecting, but bilinear resampling for changing pixel size
# Code does not remove all nodata values

# Extract desired files (i.e. environmental variables)
# Reproject (WG84 to 3035 LAEA Europe if not already in 3035)
# Clip by study area
# Save as ASC file with 1000m^2 resolution and nodata = 9999

# USER GUIDE:
  # 0. create separate folders for categorical and continuous variables
  # 1. choose between option 1 or 2
  # 2. choose input and output folder (input folder can remain zipped)
  # 3. set vector (vector) boundary
  # 4. set target_crs (target projection)
  # 5. select resampling method from line 169
# ==============================================================================
# Load necessary libraries ====
library(terra)
library(sf)
library(tictoc)
library(here)

tic() # to record script run time
# Global variables to store user preferences
global_nodata_option <- NULL
global_additional_nodata <- NULL

#________________________________________________________________________________
# 1. MANUAL SETUP: Choose your option by uncommenting ONE of the following sections: ====

# OPTION 1: No additional nodata values for any file
 global_nodata_option <- "none" # will be set to uniform value compatible for MaxEnt at end of script
 global_additional_nodata <- NULL

# OPTION 2: Setting additional nodata value(s) (same for all files) 
 #global_nodata_option <- "same_for_all"
 #global_additional_nodata <- c(-9999, 255)  # Replace with your actual nodata values
 
# OPTION 3: Setting additional nodata value(s) (different for all files) 
 # use script "Bulk processing of categorical TIF files for MaxEnt"
#________________________________________________________________________________
# 2. Define if continuous of categorical/ranked variables will be used ====
# Only one set of folder can be use at a time (also check the resampling method used near end of script is correct)

# Define continuous variable directories
input_folder <- "Environmental_Variables/Continuous variables" # you need to manually create separate folders to continuous and categorical variables.
output_folder <-"Environmental_Variables/Processed continuous variables"
 
# Define categorical variable directories
#input_folder <-here::here("Environmental_Variables/Categorical variables")
#output_folder <- here::here("Environmental_Variables/Processed categorical variables")

#________________________________________________________________________________

# Create output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

extract_dir <- file.path(output_folder, "zip_extraction")
if (!dir.exists(extract_dir)) {
  dir.create(extract_dir, recursive = TRUE)
}

# 3. Load clipping vector (set path to folder) ====
clip_vector <- st_read(here::here("Input files/Modelling_area-regional_divisions.gpkg"))

# Convert to terra vector
clip_vect <- vect(clip_vector)

# 1. Extract desired files (i.e. environmental variables) ====
file_list <- list.files(input_folder, pattern = "\\.(tif|tiff|geotiff|asc|ascii)$", full.names = TRUE, ignore.case = TRUE)

# Function to get nodata handling preferences (now just displays current settings)
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
  }
  
  cat("===============================\n")
}

# Function to get additional nodata values based on user preference
get_additional_nodata <- function(filename) {
  cat("\nProcessing file:", basename(filename), "\n")
  
  if (global_nodata_option == "none") {
    return(NULL)
    
  } else if (global_nodata_option == "same_for_all") {
    if (!is.null(global_additional_nodata) && length(global_additional_nodata) > 0) {
      cat("Using global additional nodata values:", paste(global_additional_nodata, collapse = ", "), "\n")
    }
    return(global_additional_nodata)
  }
}

# 4. Reproject function (WGS84 to EPSG 3035 Europe if not already 3035) ====
reproject_raster <- function(r) {
  target_crs <- "EPSG:3035"
  if (crs(r) != target_crs) {
    r <- project(r, target_crs, method = "near")
  }
  return(r)
}

# Display nodata handling preferences before processing
get_nodata_preferences()

# Main processing loop ====
for (file in file_list) {
  # Get additional nodata values for this file
  additional_nodata <- get_additional_nodata(file)
  
  # Load raster
  r <- rast(file)
  
  # Handle additional nodata values if specified
  if (!is.null(additional_nodata) && length(additional_nodata) > 0) {
    vals <- values(r)
    for (nodata_val in additional_nodata) {
      vals[vals == nodata_val] <- NA
      cat("Converted", nodata_val, "values to NA\n")
    }
    values(r) <- vals
  }
  
  # 3. Reproject raster if needed ====
  r <- reproject_raster(r)
  
  # 2. Clip by study area (crop + mask) ====
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
  
  # Create template raster with desired resolution
  r_extent <- ext(xmin_aligned, xmax_aligned, ymin_aligned, ymax_aligned)
  r_template <- rast(r_extent, resolution = c(1000, 1000), crs = crs(r_masked)) #1000m resolution set here
  
  # 5. Resample to desired resolution
  #r_resampled <- resample(r_masked, r_template, method = "near") # for categorical variables (nearest neighbour)
  r_resampled <- resample(r_masked, r_template, method = "average") #best for continuous environmental data
  
  # Set NA value to 9999
  NAflag(r_resampled) <- 9999
  vals <- values(r_resampled)
  vals[is.na(vals)] <- 9999
  values(r_resampled) <- vals
  
  # Output filename with ASC extension
  output_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(file)), ".asc"))
  
  # Write raster as ASCII file with nodata flag 9999
  writeRaster(r_resampled, filename = output_file, filetype = "AAIGrid", overwrite = TRUE, NAflag = 9999)
  
  cat("Completed processing:", basename(file), "\n")
}

toc()
print("Processing complete")
