# ==============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Merging and processing spatially separated tiff tiles downloaded from WEkEO
# Accessing Copernicus data often results in multiple files for download
# Created on 16-09-2025 with Perplexity by Nicole Drewitz

# Change folder and merged file name at end of script in part I (Folder of compressed folders with tif files)
                                                            # or part II (folder of tif files)
# ==============================================================================
library(terra)
library(tools)
library(here)

#_____________________________________________________________________________
#PART I: Folder of compressed folders with TIF files ====
#_____________________________________________________________________________
# Set your folder path and run the merge
main_folder <- here::here("Input files/Tree_Cover_Density_more") #holding compressed folders of tif files
merged_raster <- merge_zip_rasters(main_folder, "Tree_Cover_Density_more.tif", cleanup = FALSE)

# Complete function to merge zip files
merge_zip_rasters <- function(main_folder, output_file = "merged.tif", cleanup = FALSE) {
  
  # Set up permanent directory
  output_directory <- here::here("Environmental_Variables/Merged spatially separate tiles")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  extract_dir <- file.path(output_directory, "zip_extraction")
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir, recursive = TRUE)
  }
  
  # Find all zip files
  zip_files <- list.files(main_folder, pattern = "\\.zip$", full.names = TRUE, ignore.case = TRUE)
  
  if (length(zip_files) == 0) {
    stop("No zip files found in: ", main_folder)
  }
  
  cat("Found", length(zip_files), "zip files\n")
  
  tif_files <- character()
  
  # Extract and find TIFs
  for (zip_file in zip_files) {
    cat("Processing:", basename(zip_file), "\n")
    extract_path <- file.path(extract_dir, tools::file_path_sans_ext(basename(zip_file)))
    
    # Create extraction directory
    if (!dir.exists(extract_path)) {
      dir.create(extract_path, recursive = TRUE)
    }
    
    # Extract zip file
    utils::unzip(zip_file, exdir = extract_path)
    
    # Find TIF files
    tifs <- list.files(extract_path, pattern = "\\.(tif|tiff)$", 
                       full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    
    if (length(tifs) > 0) {
      tif_files <- c(tif_files, tifs)
      cat("  Found", length(tifs), "TIF file(s)\n")
    } else {
      cat("  No TIF files found in this archive\n")
    }
  }
  
  if (length(tif_files) == 0) {
    stop("No TIF files found in any zip archives")
  }
  
  cat("\nTotal TIF files found:", length(tif_files), "\n")
  cat("Loading and merging rasters...\n")
  
  # Load and merge rasters ====
  rasters <- lapply(tif_files, function(x) {
    tryCatch({
      terra::rast(x)
    }, error = function(e) {
      warning("Could not load: ", basename(x))
      return(NULL)
    })
  })
  
  # Remove NULL entries ====
  rasters <- rasters[!sapply(rasters, is.null)]
  
  if (length(rasters) == 0) {
    stop("No rasters could be loaded successfully")
  }
  
  cat("Successfully loaded", length(rasters), "rasters\n")
  
  # Merge using mosaic
  if (length(rasters) == 1) {
    merged <- rasters[[1]]
  } else {
    merged <- do.call(terra::mosaic, rasters)
  }
  
  # Save to permanent directory ====
  output_path <- file.path(output_directory, output_file)
  cat("Saving merged raster to:", output_path, "\n")
  
  terra::writeRaster(merged, output_path, overwrite = TRUE,
                     gdal = c("COMPRESS=LZW", "TILED=YES"))
  
  # Clean up only if requested
  if (cleanup) {
    cat("Cleaning up extracted files...\n")
    unlink(extract_dir, recursive = TRUE)
  } else {
    cat("Extracted files retained in:", extract_dir, "\n")
  }
  
  cat("Merge complete!\n")
  cat("Output dimensions:", dim(merged), "\n")
  
  return(merged)
}

#Completed: Forest_type_       north / IceFar / South / East / EastFin
#Completed: Tree_Cover_Density_north / IceFar / South / East / EastFin

#_____________________________________________________________________________
#PART II: Folder of TIF files ====
#_____________________________________________________________________________

# Load required libraries
library(terra)
library(tools)

# Set folder path
folder_path <- here::here("Input files/Merged forest variables")

# Function to merge TIFF files by prefix
merge_tiff_by_prefix <- function(folder_path, prefix, output_name) {
  # Get all TIFF files in the folder
  all_files <- list.files(folder_path, pattern = "\\.tif$|\\.tiff$", full.names = TRUE, ignore.case = TRUE)
  
  # Filter files that start with the specified prefix
  target_files <- all_files[grepl(paste0("^", prefix), basename(all_files))]
  
  if (length(target_files) == 0) {
    cat("No files found with prefix:", prefix, "\n")
    return(NULL)
  }
  
  cat("Found", length(target_files), "files with prefix:", prefix, "\n")
  cat("Files to merge:\n")
  for (f in target_files) {
    cat(" -", basename(f), "\n")
  }
  
  # Load all rasters ====
  cat("Loading rasters...\n")
  raster_list <- lapply(target_files, function(x) {
    tryCatch({
      rast(x)
    }, error = function(e) {
      cat("Error loading", basename(x), ":", e$message, "\n")
      return(NULL)
    })
  })
  
  # Remove any NULL entries (failed loads) ====
  raster_list <- raster_list[!sapply(raster_list, is.null)]
  
  if (length(raster_list) == 0) {
    cat("No rasters could be loaded successfully\n")
    return(NULL)
  }
  
  if (length(raster_list) == 1) {
    cat("Only one raster found, copying to output name\n")
    merged_raster <- raster_list[[1]]
  } else {
    # Merge rasters - this will handle overlapping areas
    cat("Merging rasters...\n")
    merged_raster <- do.call(mosaic, c(raster_list, fun = "mean"))
  }
  
  # Create output file path ====
  output_path <- file.path(folder_path, paste0(output_name, ".tif"))
  
  # Write merged raster ====
  cat("Writing merged raster to:", output_path, "\n")
  writeRaster(merged_raster, output_path, overwrite = TRUE)
  
  cat("Successfully merged", length(raster_list), "rasters into", output_name, ".tif\n\n")
  
  return(merged_raster)
}

# Main execution ====
cat("Starting TIFF merging process...\n")
cat("Working directory:", folder_path, "\n\n")

# Check if folder exists
if (!dir.exists(folder_path)) {
  stop("Folder does not exist: ", folder_path)
}

# Merge Tree_Cover_Density files ====
cat("=== Merging Tree_Cover_Density files ===\n")
tree_cover <- merge_tiff_by_prefix(folder_path, "Tree_Cover_Density", "Tree_Cover_Density")

# Merge Forest_type files ====
cat("=== Merging Forest_type files ===\n")
forest_type <- merge_tiff_by_prefix(folder_path, "Forest_type", "Forest_type")

cat("Merging process completed!\n")

# Optional: Display summary information ====
if (!is.null(tree_cover)) {
  cat("\nTree_Cover_Density summary:\n")
  cat("Dimensions:", dim(tree_cover)[1], "x", dim(tree_cover)[2], "cells\n")
  cat("Extent:", as.vector(ext(tree_cover)), "\n")
  cat("Data range:", range(values(tree_cover), na.rm = TRUE), "\n")
}

if (!is.null(forest_type)) {
  cat("\nForest_type summary:\n")
  cat("Dimensions:", dim(forest_type)[1], "x", dim(forest_type)[2], "cells\n")
  cat("Extent:", as.vector(ext(forest_type)), "\n")
  cat("Data range:", range(values(forest_type), na.rm = TRUE), "\n")
}
