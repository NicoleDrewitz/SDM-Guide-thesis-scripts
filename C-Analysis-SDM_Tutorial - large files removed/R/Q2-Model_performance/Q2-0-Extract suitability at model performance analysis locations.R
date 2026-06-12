# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q2
# Create CSV for final model performance summary
# Created with claude.ai (sonnet 4.5) on 2025-10-28 by Nicole Drewitz
# Last updated 2025-10-29

# must assign files paths
# =============================================================================

# =============================================================================
# INSTALL AND LOAD REQUIRED LIBRARIES
# =============================================================================

library(raster)
library(here)

# 1. Load rasters

current <- raster(here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_1 - 1981-2010.asc"))
ssp126 <- raster(here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_2 - ssp126.asc"))
ssp370 <- raster(here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_3 - ssp370.asc"))
ssp585 <- raster(here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_4 - ssp585.asc"))

# 2. Read your sampling location files
absence_data  <- read.csv(here::here("Input_files/V2_cloudberry_absence_analysis.csv"),
                          stringsAsFactors = FALSE)
presence_data <- read.csv(here::here("Input_files/V2_cloudberry_presence_analysis.csv"),
                          stringsAsFactors = FALSE)

# If needed, rename coordinate columns to longitude/latitude to match your later code
# names(presence_data)[names(presence_data) == "x"] <- "longitude"
# names(presence_data)[names(presence_data) == "y"] <- "latitude"
# names(absence_data)[names(absence_data) == "x"] <- "longitude"
# names(absence_data)[names(absence_data) == "y"] <- "latitude"

# 3. Extract raster values at those coordinates
presence_vals <- data.frame(
  current = extract(current, presence_data[, c("longitude", "latitude")]),
  ssp126  = extract(ssp126,  presence_data[, c("longitude", "latitude")]),
  ssp370  = extract(ssp370,  presence_data[, c("longitude", "latitude")]),
  ssp585  = extract(ssp585,  presence_data[, c("longitude", "latitude")])
)

absence_vals <- data.frame(
  current = extract(current, absence_data[, c("longitude", "latitude")]),
  ssp126  = extract(ssp126,  absence_data[, c("longitude", "latitude")]),
  ssp370  = extract(ssp370,  absence_data[, c("longitude", "latitude")]),
  ssp585  = extract(ssp585,  absence_data[, c("longitude", "latitude")])
)

# 4. Build combined results dataframe (using your structure)
combined_results <- data.frame(
  sample.No = 1:(nrow(presence_data) + nrow(absence_data)),
  presence.absence = c(rep(1, nrow(presence_data)), rep(0, nrow(absence_data))),
  longitude = c(presence_data$longitude, absence_data$longitude),
  latitude  = c(presence_data$latitude,  absence_data$latitude),
  current   = c(presence_vals$current,   absence_vals$current),
  ssp126    = c(presence_vals$ssp126,    absence_vals$ssp126),
  ssp370    = c(presence_vals$ssp370,    absence_vals$ssp370),
  ssp585    = c(presence_vals$ssp585,    absence_vals$ssp585)
)

# 5. Write to CSV
write.csv(combined_results,
          file = here::here("Data/cloudberry_suitability_at_analysis_locations.csv"),
          row.names = FALSE)
