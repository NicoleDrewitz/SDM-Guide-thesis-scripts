# =============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Check error message for "CRS does not match"
# paste after receiving error message, or open files in GIS software
# Created with claude.ai on 18 October, 2025 by Nicole Drewitz
# =============================================================================

library(terra)

crs_list <- lapply(raster_files, function(f) crs(rast(f)))
# Extract the character representation of each CRS
crs_texts <- sapply(crs_list, function(x) terra::crs(x, proj=TRUE))
unique_crs <- unique(crs_texts)

if (length(unique_crs) > 1) {
  message("Warning: Not all rasters have identical CRS definitions.")
  message("Different CRS found:")
  # Print each unique CRS separately
  for (i in seq_along(unique_crs)) {
    message(paste0("[", i, "] ", unique_crs[i]))
  }
} else {
  message("All rasters have the same CRS.")
}

# Warning: Not all rasters have identical CRS definitions.
# Different CRS found:
#  [1] +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs
#  [2] 

# proj=laea is EPGS:3035 Lambert Azimuthal Equal Area (Europe)
# [2] blank = unassigned (file may match format for correct projection)
