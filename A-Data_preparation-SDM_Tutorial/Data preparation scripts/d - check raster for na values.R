# ==============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Checking for NA raster values within modelling area
  # identifies areas that may be excluded from model if the variable is important
  # variable dependent a manual fix could be possible
      # For example, missing low snow water equivalent values
          # set to 0 where annual snow cover days is also 0
# Create 2025-12-1 with claude.ai by Nicole Drewitz
# ==============================================================================

# Setup ====
library(terra)
library(sf)
library(here)

# Define folder to check
folder_path <- here::here("Environmental_Variables/Environmental variable dataset for modelling")

# Load vector for clipping
clip_vector <- st_read(here::here("Input files/Modelling_area-regional_divisions.gpkg"))

# Find all raster files (common extensions)
raster_files <- list.files(folder_path, 
                           pattern = "\\.(tif|asc|img|grd)$", 
                           full.names = TRUE,
                           ignore.case = TRUE)

if (length(raster_files) == 0) {
  stop("No raster files found in the folder.")
}

cat("Checking", length(raster_files), "raster files for NA values within vector extent...\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Create output folder for NA maps ====
na_maps_folder <- file.path(folder_path, "NA_location_maps")
if (!dir.exists(na_maps_folder)) {
  dir.create(na_maps_folder, recursive = TRUE)
}

# Create a summary data frame ====
results <- data.frame(
  filename = character(),
  total_cells = numeric(),
  na_cells = numeric(),
  na_percentage = numeric(),
  has_na = logical(),
  stringsAsFactors = FALSE
)

# Check each raster ====
for (raster_file in raster_files) {
  
  filename <- basename(raster_file)
  cat("Checking:", filename, "\n")
  
  # Load raster
  r <- rast(raster_file)
  
  # First, identify which cells in the ORIGINAL raster are NA (before masking)
  original_na <- is.na(r)
  
  # Clip and mask both the original raster and the NA map to vector
  r_clipped <- mask(crop(r, vect(clip_vector)), vect(clip_vector))
  original_na_clipped <- mask(crop(original_na, vect(clip_vector)), vect(clip_vector))
  
  # Count cells: only count NAs that existed BEFORE masking (real data gaps)
  na_cells <- global(original_na_clipped == 1, "sum", na.rm = TRUE)[1,1]
  
  # Count valid cells (had data before masking)
  valid_cells <- global(original_na_clipped == 0, "sum", na.rm = TRUE)[1,1]
  
  # Total cells within vector
  total_cells <- valid_cells + na_cells
  
  # Calculate percentage
  na_percentage <- (na_cells / total_cells) * 100
  
  # Check if has any NA
  has_na <- na_cells > 0
  
  # Add to results
  results <- rbind(results, data.frame(
    filename = filename,
    total_cells = total_cells,
    na_cells = na_cells,
    na_percentage = round(na_percentage, 2),
    has_na = has_na
  ))
  
  # Print summary ====
  if (has_na) {
    cat("  ??? Contains NA values:", na_cells, "cells (", 
        round(na_percentage, 2), "%)\n")
    
    # Create NA location map - only showing original NAs within vector
    # Mask the original_na_clipped so ocean areas are not displayed
    na_map <- original_na_clipped
    
    # Save as TIF ====
    output_tif <- file.path(na_maps_folder, 
                            paste0(tools::file_path_sans_ext(filename), "_NA_map.tif"))
    writeRaster(na_map, output_tif, overwrite = TRUE)
    cat("  NA map saved to:", output_tif, "\n")
    
    # Plot
    par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
    
    # Plot original data
    plot(r_clipped, main = paste("Original:", filename), 
         col = terrain.colors(100))
    plot(vect(clip_vector), add = TRUE)
    
    # Plot NA locations (1 = NA, 0 = has data) ====
    plot(na_map, main = "NA Locations (red = NA)", 
         col = c("lightblue", "red"), legend = FALSE)
    plot(vect(clip_vector), add = TRUE, border = "black", lwd = 1)
    
    cat("  Press [Enter] to continue to next raster...")
    readline()
    
  } else {
    cat("  ??? No NA values found\n")
  }
  
  cat("\n")
}

cat(paste(rep("=", 70), collapse = ""), "\n") 
cat("SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Print full results table ====
print(results, row.names = FALSE)

# Summary statistics
cat("\n")
cat("Files with NA values:", sum(results$has_na), "/", nrow(results), "\n")
cat("Files without NA values:", sum(!results$has_na), "/", nrow(results), "\n")

# Save results to CSV ====
#output_file <- file.path(folder_path, "NA_check_results.csv")
#write.csv(results, output_file, row.names = FALSE)
#cat("\nResults saved to:", output_file, "\n")
cat("NA location maps saved to:", na_maps_folder, "\n")
