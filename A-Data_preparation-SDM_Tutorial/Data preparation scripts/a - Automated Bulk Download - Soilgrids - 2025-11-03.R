# ==============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# SoilGrids Data Download Script
# Author: Nicole Drewitz, Drafted with Claude AI, 10-2025
# Last updated: 3-11-2025

# Purpose: Downloads global soil files first, then clips them in a separate step. Adapted from WebDAV terra documentation. https://docs.isric.org/globaldata/soilgrids/SoilGrids_faqs_01.html
# currently this script averages the mean value for measurements at 0-30cm depths.
# Downloads variables from https://files.isric.org/soilgrids/latest/data/ 
# Usage: Update your ROI (region of interest), soil variables/depths of interest

# Data citation:
# ISRIC. (2025, April 28). SoilGrids Documentation. https://docs.isric.org/globaldata/soilgrids/ 
# 2020 May, ISRIC, SoilGrids Documentation 
# Visualized at SoilGrids250m 2.0 

#(Common soil chemical and physical properties: 
#    Poggio, L., de Sousa, L. M., Batjes, N. H., Heuvelink, G. B. M., Kempen, B., Ribeiro, E., and Rossiter, D.: SoilGrids 2.0: producing soil information for the globe with quantified spatial uncertainty, SOIL, 7, 217–240, 2021. https://doi.org/10.5194/soil-7-217-2021  
# Soil water content at different pressure heads: 
    #Turek, M.E., Poggio, L., Batjes, N. H., Armindo, R. A., de Jong van Lier, Q., de Sousa, L.M., Heuvelink, G. B. M. : Global mapping of volumetric water retention at 100, 330 and 15 000 cm suction using the WoSIS database, International Soil and Water Conservation Research, 11-2, 225-239, 2023. https://doi.org/10.1016/j.iswcr.2022.08.001)
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

# Load required libraries
library(Rcpp)
library(terra)
library(sf)
library(here)

# Your desired output directory
output_dir <- here::here("Environmental_Variables/Soilgrids_data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: DEFINE LAYER OF INTEREST
# ==============================================================================

# Base URLs for different resolutions
base_urls <- list(
  "250m" = "https://files.isric.org/soilgrids/latest/data/",
  "1000m" = "https://files.isric.org/soilgrids/latest/data_aggregated/1000m/",
  "5000m" = "https://files.isric.org/soilgrids/latest/data_aggregated/5000m/"
)

# EDIT THESE VARIABLES AS NEEDED:
# Available variables: "bdod", "cec", "cfvo", "clay", "nitrogen", "ocd", "ocs", 
#                      "soc", "phh2o", "sand", "silt", "wv0010", "wv1500", "wv0033"
vois <- c("bdod", "cec", "cfvo", "clay", "nitrogen", "ocd", "ocs", 
          "soc", "phh2o", "sand", "silt", "wv0010", "wv1500", "wv0033")

# Available depths: "0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm"
depths <- c("0-5cm", "5-15cm", "15-30cm")

# Available quantiles: "mean", "Q0.05", "Q0.5", "Q0.95", "uncertainty"
quantiles <- c("mean")

# Available resolutions: "250m", "1000m", "5000m"
# Note: 1000m and 5000m only available for mean quantile
resolution <- "1000m"

# NEW OPTIONS: How to aggregate multiple depths?
# aggregate_depths: Set to "none", "sum", or "average"
#   - "none": Keep each depth as a separate file
#   - "sum": Sum all selected depths (useful for stocks like carbon)
#   - "average": Average all selected depths (useful for concentrations/properties)
aggregate_depths <- "average"

# ==============================================================================
# STEP 2: DEFINE REGION OF INTEREST (ROI)
# ==============================================================================

# Choose ONE of the following methods:

# METHOD 2a: Manual bounding box coordinates
# Uncomment and edit these lines if using manual coordinates:
#in_crs <- "EPSG:4326"  # CRS of your bounding box
#xmin <- -73.24
#xmax <- 33.5
#ymin <- 54.56
#ymax <- 45
#bb <- ext(c(xmin, xmax, ymin, ymax))

# METHOD 2b: Extract from ROI file
# Uncomment and edit these lines if using a ROI file:
roi_path <- st_read(here::here("Input files/Modelling_area-regional_divisions.gpkg")) # or .gdb, .tif, etc to use a file with the modelling area to set the ROI
roi <- vect(roi_path)  # use vect() for vector or rast() for raster
bb <- ext(roi)
in_crs <- crs(roi)

# METHOD 2c: Global data (no ROI)
# Leave bb and in_crs undefined to download global data
# Note: This may take considerable time!

# ==============================================================================
# STEP 3: BUILD FILE PATHS AND DOWNLOAD DATA
# ==============================================================================

# Define CRS for SoilGrids layers
rst_crs <- "ESRI:54052"  # Interrupted Goode Homolosine projection

# Loop through all selected variables and download
for (voi in vois) {
  # Adjust depths for OCS (organic carbon stock)
  # OCS is only available as 0-30cm, so force single depth and disable aggregation for this variable
  if (voi == "ocs") {
    cur_depths <- "0-30cm"
    cur_aggregate <- "none"  # Can't aggregate a single depth
  } else {
    cur_depths <- depths
    cur_aggregate <- aggregate_depths
  }
  
  for (quantile in quantiles) {
    
    # Initialize list to store rasters for aggregation
    depth_rasters <- list()
    
    for (depth in cur_depths) {
      
      # Build file path based on resolution
      if (resolution == "250m") {
        rstFile <- paste0("/vsicurl/", base_urls[[resolution]], 
                          voi, "/", voi, "_", depth, "_", quantile, ".vrt")
      } else if (resolution == "1000m") {
        rstFile <- paste0("/vsicurl/", base_urls[[resolution]], 
                          voi, "/", voi, "_", depth, "_", quantile, "_1000.tif")
      } else {
        rstFile <- paste0("/vsicurl/", base_urls[[resolution]], 
                          voi, "/", voi, "_", depth, "_quantile_5000.tif")
      }
      
      cat("Processing:", rstFile, "\n")
      
      # Load raster
      rst <- rast(rstFile)
      crs(rst) <- rst_crs
      
      # Subset to ROI if bb is defined
      if (exists("bb")) {
        bb_proj <- project(bb, from = crs(in_crs), to = crs(rst_crs))
        window(rst) <- bb_proj
      }
      
      # Optional: reproject to ROI CRS
      # Uncomment if you want output in original ROI CRS:
      # if (exists("in_crs")) {
      #   rst <- project(rst, in_crs)
      # }
      
      if (cur_aggregate != "none" && length(cur_depths) > 1) {
        # Store raster for later aggregation
        depth_rasters[[depth]] <- rst
        cat("  Stored for", cur_aggregate, "\n")
      } else {
        # Save individual depth file
        output_filename <- paste0(voi, "_", depth, "_", quantile, ".tif")
        output_path <- file.path(output_dir, output_filename)
        
        cat("Saving to:", output_path, "\n")
        writeRaster(rst, output_path, overwrite = TRUE)
        
        # Plot for visual check
        plot(rst, main = paste(voi, depth, quantile))
      }
    }
    
    # Aggregate depths if option is enabled and multiple depths selected
    if (cur_aggregate != "none" && length(depth_rasters) > 1) {
      cat("Aggregating", length(depth_rasters), "depth layers for", voi, quantile, "using", cur_aggregate, "\n")
      
      # Perform aggregation based on method
      if (cur_aggregate == "sum") {
        rst_agg <- Reduce("+", depth_rasters)
        agg_label <- "summed"
      } else if (cur_aggregate == "average") {
        rst_agg <- Reduce("+", depth_rasters) / length(depth_rasters)
        agg_label <- "averaged"
      }
      
      # Create output filename for aggregated raster
      depth_range <- paste0(gsub("cm", "", cur_depths[1]), "-", 
                            gsub(".*-", "", cur_depths[length(cur_depths)]))
      output_filename <- paste0("soil_", voi, "_", depth_range, "_", quantile, "_", agg_label, ".tif")
      output_path <- file.path(output_dir, output_filename)
      
      cat("Saving", agg_label, "raster to:", output_path, "\n")
      writeRaster(rst_agg, output_path, overwrite = TRUE)
      
      # Plot for visual check
      plot(rst_agg, main = paste(voi, toupper(cur_aggregate), "of depths", quantile))
    }
  }
}

cat("\nAll downloads complete! Files saved to:", output_dir, "\n")

# ==============================================================================
# OPTIONAL: DOWNLOAD WRB SOIL CLASSIFICATION DATA
# ==============================================================================

# Uncomment this section if you want WRB soil classification layers:

# Available classes: "Acrisols", "Histosols", "Podzols", etc.
wrb_class <- "Histosols"

# WRB layers use different CRS
wrb_crs <- "EPSG:4326"

# Base URL for WRB data
wrb_base_url <- "https://files.isric.org/soilgrids/latest/data/wrb/"

# List of all possible tile numbers (you may need to expand this list)
# These were found by browsing the directory structure
all_possible_tiles <- c("7", "8", "9", "10", "11", "12", "13", "14", "15", 
                        as.character(16:100), as.character(101:144), 
                        as.character(149:150), as.character(151:200), 
                        as.character(201:250), as.character(251:300), "303", "304", 
                        as.character(310:340), "342", as.character(354:360), 
                        as.character(363:375), as.character(377:400), as.character(401:450), 
                        as.character(452:470), as.character(472:490), as.character(494:510))

cat("\n=== Finding WRB tiles that intersect with your ROI ===\n")

# Check if bounding box is defined
if (!exists("bb")) {
  stop("Please define a bounding box (bb) in Step 2 before running WRB download")
}

# Project bounding box to WRB CRS
bb_wrb <- project(bb, from = crs(in_crs), to = crs(wrb_crs))

cat("Your ROI extent (EPSG:4326):\n")
cat("  xmin:", bb_wrb[1], "xmax:", bb_wrb[2], "\n")
cat("  ymin:", bb_wrb[3], "ymax:", bb_wrb[4], "\n\n")

# Find tiles that intersect with ROI
tiles_to_download <- c()
tiles_with_data <- c()

cat("Testing tiles for intersection with your ROI...\n")

for (tile in all_possible_tiles) {
  rstFile_wrb <- paste0("/vsicurl/", wrb_base_url, wrb_class, "/", tile, ".tif")
  
  tryCatch({
    # Try to load tile
    rst_wrb <- rast(rstFile_wrb)
    crs(rst_wrb) <- wrb_crs
    
    # Get tile extent
    tile_ext <- ext(rst_wrb)
    
    # Check if tile intersects with ROI
    if (tile_ext[1] <= bb_wrb[2] && tile_ext[2] >= bb_wrb[1] &&
        tile_ext[3] <= bb_wrb[4] && tile_ext[4] >= bb_wrb[3]) {
      tiles_to_download <- c(tiles_to_download, tile)
      cat("  ✓ Tile", tile, "intersects ROI\n")
    }
    
    tiles_with_data <- c(tiles_with_data, tile)
    
  }, error = function(e) {
    # Tile doesn't exist or can't be accessed - skip silently
  })
}

if (length(tiles_to_download) == 0) {
  cat("\n⚠ No tiles found that intersect with your ROI.\n")
  cat("This could mean:\n")
  cat("  1. Your ROI doesn't contain any", wrb_class, "\n")
  cat("  2. The tile list needs to be updated\n")
  cat("  3. There's an issue with the data access\n")
} else {
  cat("\n=== Downloading", length(tiles_to_download), "tiles ===\n")
  cat("Tiles to download:", paste(tiles_to_download, collapse = ", "), "\n\n")
  
  # Download and process each tile
  for (tile in tiles_to_download) {
    rstFile_wrb <- paste0("/vsicurl/", wrb_base_url, wrb_class, "/", tile, ".tif")
    
    cat("Processing tile:", tile, "\n")
    
    tryCatch({
      rst_wrb <- rast(rstFile_wrb)
      crs(rst_wrb) <- wrb_crs
      
      # Crop to ROI
      rst_wrb_crop <- crop(rst_wrb, bb_wrb)
      
      # Optional: reproject to ROI CRS
      # if (exists("in_crs") && in_crs != wrb_crs) {
      #   rst_wrb_crop <- project(rst_wrb_crop, in_crs)
      # }
      
      # Save tile
      wrb_filename <- paste0("soil_", wrb_class, "_tile_", tile, ".tif")
      wrb_path <- file.path(output_dir, wrb_filename)
      writeRaster(rst_wrb_crop, wrb_path, overwrite = TRUE)
      
      cat("  Saved:", wrb_filename, "\n")
      
    }, error = function(e) {
      cat("  Error processing tile", tile, ":", conditionMessage(e), "\n")
    })
  }
  
  cat("\n=== Merging tiles into single file ===\n")
  
  # Load all downloaded tiles
  tile_files <- list.files(output_dir, pattern = paste0("wrb_", wrb_class, "_tile_.*\\.tif$"), 
                           full.names = TRUE)
  
  if (length(tile_files) > 1) {
    tile_rasters <- lapply(tile_files, rast)
    
    # Merge tiles
    merged_wrb <- do.call(merge, tile_rasters)
    
    # Save merged file
    merged_filename <- paste0("wrb_", wrb_class, "_merged.tif")
    merged_path <- file.path(output_dir, merged_filename)
    writeRaster(merged_wrb, merged_path, overwrite = TRUE)
    
    cat("Merged file saved:", merged_filename, "\n")
    
    # Plot merged result
    plot(merged_wrb, main = paste("WRB", wrb_class, "(merged)"))
    
  } else if (length(tile_files) == 1) {
    cat("Only one tile downloaded - no merging needed\n")
    plot(rast(tile_files[1]), main = paste("WRB", wrb_class))
  }
}

# Delete individual tile files after merging
file.remove(tile_files)
cat("Individual tile files have been deleted.\n")

cat("\nWRB processing complete!\n")
