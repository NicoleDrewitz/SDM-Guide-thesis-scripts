# ==============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Assign 0 values where peat is is not found (assumption)
# Author: Nicole Drewitz with Perplexity AI
# Lasted updated: 23/09/2025
# ==============================================================================

library(terra)
library(here)

# Load the ASC raster file
r <- rast(here::here("Environmental_Variables/Peat.asc")) #example file
# Identify the NoData value of the raster
nodata_value <- terra::NAflag(r)

# Replace NoData (NA) values with 0
r[is.na(r[])] <- 0

# Optional: explicitly set NAflag to a value if needed (e.g., 9999)
# terra::NAflag(r) <- 9999 

# Save the raster as new ASCII file
writeRaster(r, here::here("Environmental_Variables/Processed categorial variables/Peat.asc"))