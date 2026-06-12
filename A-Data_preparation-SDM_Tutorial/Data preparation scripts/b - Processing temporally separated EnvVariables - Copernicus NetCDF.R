# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Copernicus NetCDF Raster Extraction
# Merging tif files separated temporally
# By Nicole Drewitz with Claude.ai
# Last updated 15-9-2025

# Load required libraries
library(terra)
library(sf)
library(here)


# Set working directory and file paths
# Adjust these paths to match your file locations
nc_folder <- here::here("Environmental_Variables/NDVI")  # Folder containing your .nc files
modelling_area_path <- here::here("Input files/Modelling_area-regional_divisions.gpkg")  # Path to a modelling area vector file
output_folder <- here::here("Environmental_Variables/NDVI - clipped means")  # Output folder for TIFF files

#____________________________________________
# Check file usability
#___________________________________

# DEBUGGING STEPS

# 1. Check if the folder exists and list all files
cat("=== STEP 1: Checking folder contents ===\n")
cat("NC folder exists:", dir.exists(nc_folder), "\n")
all_files <- list.files(nc_folder, full.names = TRUE)
cat("Total files in folder:", length(all_files), "\n")
cat("All files:\n")
for(i in 1:min(10, length(all_files))) {  # Show first 10 files
  cat(" -", basename(all_files[i]), "\n")
}
if(length(all_files) > 10) cat("... and", length(all_files) - 10, "more files\n")

# 2. Check specifically for .nc files
cat("\n=== STEP 2: Looking for .nc files ===\n")
nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
cat("NetCDF files found:", length(nc_files), "\n")
if(length(nc_files) > 0) {
  cat("NC files:\n")
  for(file in nc_files) {
    cat(" -", basename(file), "\n")
  }
} else {
  # Try different patterns
  nc_files_alt1 <- list.files(nc_folder, pattern = "\\.NC$", full.names = TRUE)
  nc_files_alt2 <- list.files(nc_folder, pattern = "\\.[nN][cC]$", full.names = TRUE)
  cat("Trying .NC (uppercase):", length(nc_files_alt1), "files\n")
  cat("Trying case-insensitive:", length(nc_files_alt2), "files\n")
  
  if(length(nc_files_alt2) > 0) {
    nc_files <- nc_files_alt2
    cat("Using case-insensitive pattern. Files found:\n")
    for(file in nc_files) {
      cat(" -", basename(file), "\n")
    }
  }
}

# 3. If we found NC files, examine the first one
if(length(nc_files) > 0) {
  cat("\n=== STEP 3: Examining first NetCDF file ===\n")
  first_nc <- nc_files[1]
  cat("Examining:", basename(first_nc), "\n")
  
  tryCatch({
    # Load the first NetCDF file
    nc_rast <- rast(first_nc)
    
    cat("Number of layers:", nlyr(nc_rast), "\n")
    cat("Layer names:\n")
    layer_names <- names(nc_rast)
    for(i in 1:length(layer_names)) {
      cat(" ", i, ":", layer_names[i], "\n")
    }
    
    # Look for "mean" in layer names
    mean_layers <- grep("mean", layer_names, ignore.case = TRUE, value = TRUE)
    cat("\nLayers containing 'mean':", length(mean_layers), "\n")
    if(length(mean_layers) > 0) {
      for(layer in mean_layers) {
        cat(" -", layer, "\n")
      }
    } else {
      cat("No layers found with 'mean' in the name.\n")
      cat("Looking for other common patterns...\n")
      
      # Check for other common patterns
      avg_layers <- grep("avg|average", layer_names, ignore.case = TRUE, value = TRUE)
      if(length(avg_layers) > 0) {
        cat("Layers with 'avg' or 'average':\n")
        for(layer in avg_layers) cat(" -", layer, "\n")
      }
      
      # Maybe the layer is just the variable itself
      ndvi_layers <- grep("ndvi", layer_names, ignore.case = TRUE, value = TRUE)
      if(length(ndvi_layers) > 0) {
        cat("Layers with 'ndvi':\n")
        for(layer in ndvi_layers) cat(" -", layer, "\n")
      }
    }
    
    # Show raster info
    cat("\nRaster information:\n")
    cat("Dimensions:", nrow(nc_rast), "rows,", ncol(nc_rast), "cols\n")
    cat("Resolution:", res(nc_rast)[1], "x", res(nc_rast)[2], "\n")
    cat("CRS:", as.character(crs(nc_rast)), "\n")
    
  }, error = function(e) {
    cat("Error loading NetCDF file:", e$message, "\n")
  })
}

# 4. Check modelling area file
cat("\n=== STEP 4: Checking vector ===\n")
cat("Modelling area file exists:", file.exists(modelling_area_path), "\n")

if(file.exists(modelling_area_path)) {
  tryCatch({
    clip_shape <- vect(modelling_area_path)
    cat("modelling area loaded successfully\n")
    cat("Number of features:", nrow(clip_shape), "\n")
    cat("CRS:", as.character(crs(clip_shape)), "\n")
  }, error = function(e) {
    cat("Error loading modelling area file:", e$message, "\n")
  })
} else {
  cat("Modelling area file not found at specified path\n")
  # Check if folder exists
  shapefile_folder <- dirname(modelling_area_path)
  cat("Vector folder exists:", dir.exists(shapefile_folder), "\n")
  if(dir.exists(shapefile_folder)) {
    shp_files <- list.files(shapefile_folder, pattern = "\\.shp$", full.names = TRUE)
    cat("Shapefiles in folder:\n")
    for(file in shp_files) {
      cat(" -", basename(file), "\n")
    }
  }
}

#______________________________________
# Clipped means
#________________________________

# Create output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

print("Starting processing...")
print(paste("Looking for NetCDF files in:", nc_folder))

# Get list of all .nc files
nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
print(paste("Found", length(nc_files), "NetCDF files"))

if(length(nc_files) == 0) {
  stop("No NetCDF files found!")
}

# Load vector
print("Loading vector...")
clip_shape <- vect(modelling_area_path)
print(paste("vector loaded with", nrow(clip_shape), "features"))

# Initialize list to store mean rasters
mean_rasters <- list()

print("Processing NetCDF files...")

# Loop through each NetCDF file to extract mean layers
for(i in 1:length(nc_files)) {
  nc_file <- nc_files[i]
  print(paste("Processing file", i, "of", length(nc_files), ":", basename(nc_file)))
  
  # Load NetCDF file
  nc_rast <- rast(nc_file)
  
  # Extract mean layer
  mean_rast <- nc_rast[["mean"]]
  
  # Crop to vector extent to speed up processing
  cropped_rast <- crop(mean_rast, clip_shape)
  
  # Store in list
  mean_rasters[[i]] <- cropped_rast
  
  print(paste("  ??? Extracted mean layer from", basename(nc_file)))
}

print("Combining all mean rasters...")

# Create a raster stack from all mean layers
all_means_stack <- rast(mean_rasters)
print(paste("Created stack with", nlyr(all_means_stack), "layers"))

# Calculate the average across all time periods
print("Calculating average across all time periods...")
averaged_raster <- app(all_means_stack, fun = mean, na.rm = TRUE)

# Clip the averaged raster to the vector
print("Clipping to vector boundary...")
final_clipped <- mask(averaged_raster, clip_shape)

# Create output filename
output_file <- file.path(output_folder, "NDVI_LTS_1999-2019_averaged_mean_clipped.tif")

# Save the final averaged and clipped raster
print("Saving final TIFF file...")
writeRaster(final_clipped, output_file, overwrite = TRUE)

# Print summary information
print("=== PROCESSING COMPLETE ===")
print(paste("Input files processed:", length(nc_files)))
print(paste("Output file:", output_file))
print(paste("Output file size:", round(file.info(output_file)$size / 1024^2, 2), "MB"))

# Display some stats about the final raster
print("=== RASTER STATISTICS ===")
print(paste("Dimensions:", nrow(final_clipped), "x", ncol(final_clipped), "cells"))
print(paste("Resolution:", round(res(final_clipped)[1], 6), "x", round(res(final_clipped)[2], 6), "degrees"))
print(paste("Min value:", round(minmax(final_clipped)[1], 4)))
print(paste("Max value:", round(minmax(final_clipped)[2], 4)))

print("Done! Check your output folder for the averaged NDVI file.")

#_____________________________________________
# Processing range in values of each pixel
#______________________________________________

print("Starting NDVI range processing...")
print(paste("Looking for NetCDF files in:", nc_folder))

# Get list of all .nc files
nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
print(paste("Found", length(nc_files), "NetCDF files"))

if(length(nc_files) == 0) {
  stop("No NetCDF files found!")
}

# Load vector
print("Loading vector...")
clip_shape <- vect(modelling_area_path)
print(paste("vector loaded with", nrow(clip_shape), "features"))

# Initialize lists to store max and min rasters
max_rasters <- list()
min_rasters <- list()

print("Processing NetCDF files to extract max and min layers...")

# Loop through each NetCDF file to extract max and min layers
for(i in 1:length(nc_files)) {
  nc_file <- nc_files[i]
  print(paste("Processing file", i, "of", length(nc_files), ":", basename(nc_file)))
  
  # Load NetCDF file
  nc_rast <- rast(nc_file)
  
  # Extract max and min layers
  max_rast <- nc_rast[["max"]]
  min_rast <- nc_rast[["min"]]
  
  # Crop to vector extent to speed up processing
  cropped_max <- crop(max_rast, clip_shape)
  cropped_min <- crop(min_rast, clip_shape)
  
  # Store in lists
  max_rasters[[i]] <- cropped_max
  min_rasters[[i]] <- cropped_min
  
  print(paste("  ??? Extracted max and min layers from", basename(nc_file)))
}

print("Combining all max and min rasters...")

# Create raster stacks from all max and min layers
all_max_stack <- rast(max_rasters)
all_min_stack <- rast(min_rasters)
print(paste("Created max stack with", nlyr(all_max_stack), "layers"))
print(paste("Created min stack with", nlyr(all_min_stack), "layers"))

# Find the highest max value across all time periods for each pixel
print("Finding highest max values across all time periods...")
global_max <- app(all_max_stack, fun = max, na.rm = TRUE)

# Find the lowest min value across all time periods for each pixel
print("Finding lowest min values across all time periods...")
global_min <- app(all_min_stack, fun = min, na.rm = TRUE)

# Calculate the range (max - min) for each pixel
print("Calculating range (max - min) for each pixel...")
range_raster <- global_max - global_min

# Clip the range raster to the vector
print("Clipping to vector boundary...")
final_range_clipped <- mask(range_raster, clip_shape)

# Create output filename
output_file <- file.path(output_folder, "NDVI_LTS_1999-2019_range_clipped.tif")

# Save the final range raster
print("Saving final range TIFF file...")
writeRaster(final_range_clipped, output_file, overwrite = TRUE)

# Print summary information
print("=== PROCESSING COMPLETE ===")
print(paste("Input files processed:", length(nc_files)))
print(paste("Output file:", output_file))
print(paste("Output file size:", round(file.info(output_file)$size / 1024^2, 2), "MB"))

# Display some stats about the final raster
print("=== RASTER STATISTICS ===")
print(paste("Dimensions:", nrow(final_range_clipped), "x", ncol(final_range_clipped), "cells"))
print(paste("Resolution:", round(res(final_range_clipped)[1], 6), "x", round(res(final_range_clipped)[2], 6), "degrees"))
print(paste("Min range value:", round(minmax(final_range_clipped)[1], 4)))
print(paste("Max range value:", round(minmax(final_range_clipped)[2], 4)))

print("Done! The range raster shows NDVI variability (highest max - lowest min) for each pixel across all time periods.")