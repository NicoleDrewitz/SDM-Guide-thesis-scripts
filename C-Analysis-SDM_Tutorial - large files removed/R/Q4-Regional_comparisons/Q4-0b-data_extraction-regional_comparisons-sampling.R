# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q4:
# Extract Habitat Suitability Values at Random Sampling Locations
# created on 20 feb, 2026 by Nicole Drewitz with Claude AI
# =============================================================================
# Extracts values from current and future scenario ASC rasters at sampling
# points and saves results with a 'habitat_suitability' column.
# change sampling location file and file name as needed
# hard code field name on line 104
# =============================================================================

library(here)
library(terra)
library(dplyr)
library(readr)
library(terra)

# Paths relative to the RProject root
sampling_locations_dir <- here::here("Data/Q4-random_sampling_locations")
output_dir <- here::here("Data/Q4-random_sampling_locations-with_values")
output_summary <- here::here("analysis_output_data")

# -----------------------------------------------------------------------------
# 1. Load rasters
# -----------------------------------------------------------------------------
asc_dir <- here::here("Input_files/MaxEnt_files")

current <- rast(file.path(asc_dir, "Rubus_chamaemorus-Scenario_1 - 1981-2010.asc"))
ssp126 <- rast(file.path(asc_dir, "Rubus_chamaemorus-Scenario_2 - ssp126.asc"))
ssp370 <- rast(file.path(asc_dir, "Rubus_chamaemorus-Scenario_3 - ssp370.asc"))
ssp585 <- rast(file.path(asc_dir, "Rubus_chamaemorus-Scenario_4 - ssp585.asc"))

crs_laea <- "EPSG:3035"
crs(current) <- crs_laea
crs(ssp126)  <- crs_laea
crs(ssp370)  <- crs_laea
crs(ssp585)  <- crs_laea

message("--- Raster diagnostics ---")
message("current extent: ", paste(as.vector(ext(current)), collapse = ", "))
message("current CRS: ", terra::crs(current, describe = TRUE)$name)

# -----------------------------------------------------------------------------
# 2. Load sampling locations
# -----------------------------------------------------------------------------
FIELD         <- "ECO_Cor" # Inland, ECO_Cor, CoastCNT

points_file <- file.path(
  sampling_locations_dir, # update file name to match!
  #"Q4-0a-random_points_for_sampling-Inland-n1000each-mindist10000m.csv" # greenland and svalbard excluded in previous script
  #"Q4-0a-random_points_for_sampling-CoastCNT-n200each-mindist2000m.csv"
  "Q4-0a-random_points_for_sampling-ECO_Cor-n300each-mindist5000m.csv"
)

sampling_pts <- read_csv(points_file, show_col_types = FALSE)

if (!all(c("x", "y") %in% names(sampling_pts))) {
  stop("Sampling locations CSV must contain columns named 'x' and 'y'.")
}

# Diagnostics — check coordinate ranges match raster extent
message("--- Point coordinate ranges ---")
message("x range: ", min(sampling_pts$x), " to ", max(sampling_pts$x))
message("y range: ", min(sampling_pts$y), " to ", max(sampling_pts$y))

# Convert to SpatVector in EPSG:3035
sp_pts <- vect(sampling_pts, geom = c("x", "y"), crs = crs_laea)

# -----------------------------------------------------------------------------
# 3. Extract raster values
# -----------------------------------------------------------------------------
out_df <- sampling_pts %>%
  mutate(
    suitability_current = terra::extract(current, sp_pts)[, 2],
    suitability_ssp126  = terra::extract(ssp126,  sp_pts)[, 2],
    suitability_ssp370  = terra::extract(ssp370,  sp_pts)[, 2],
    suitability_ssp585  = terra::extract(ssp585,  sp_pts)[, 2]
  )

# Diagnostics — check how many NAs were produced
message("--- NA counts in extracted columns ---")
message("suitability_current NAs: ", sum(is.na(out_df$suitability_current)), " / ", nrow(out_df))
message("suitability_ssp126  NAs: ", sum(is.na(out_df$suitability_ssp126)),  " / ", nrow(out_df))
message("suitability_ssp370  NAs: ", sum(is.na(out_df$suitability_ssp370)),  " / ", nrow(out_df))
message("suitability_ssp585  NAs: ", sum(is.na(out_df$suitability_ssp585)),  " / ", nrow(out_df))

# -----------------------------------------------------------------------------
# 4. Save
# -----------------------------------------------------------------------------
out_file <- file.path(output_dir,
#"Q4-0b-1-random_points_suitability-Inland-n1000each-mindist10000m.csv" # Svalbard and Greenland exclude in previous script
#"Q4-0b-2-random_points_suitability-CoastCNT-n200each-mindist2000m.csv"
"Q4-0b-3-random_points_suitability-ECO_Cor-n300each-mindist5000m.csv"
)
write_csv(out_df, out_file)
out_df
message("Done! File saved to: ", out_file)

# -----------------------------------------------------------------------------
# 5. Summary statistics by category
# -----------------------------------------------------------------------------
summary_df <- out_df %>%
  group_by(ECO_Cor) %>% # hard code column name!!! ECO_Cor, Inland, ECO_Cor, CoastCNT
  summarise(
    n_points       = n(),
    mean_Current      = mean(suitability_current, na.rm = TRUE),
    `mean_SSP1-2.6` = mean(suitability_ssp126,  na.rm = TRUE),
    `mean_SSP3-7.0` = mean(suitability_ssp370,  na.rm = TRUE),
    `mean_SSP5-8.5` = mean(suitability_ssp585,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    change_126 = `mean_SSP1-2.6` - mean_Current,
    change_370 = `mean_SSP3-7.0` - mean_Current,
    change_585 = `mean_SSP5-8.5` - mean_Current
  )

# Extract label from source filename (between "suitability-" and "-n")
source_filename <- basename(points_file)
label <- stringr::str_extract(source_filename, "(?<=suitability-)(.+?)(?=-n)")
# If source uses the input naming pattern instead, fall back to out_file name
if (is.na(label)) {
  label <- stringr::str_extract(basename(out_file), "(?<=suitability-)(.+?)(?=-n)")
}

summary_file <- file.path(output_summary,
                          paste0("Q4-0b-suitability_summary-", label, "-SUMMARY.csv")
)
write_csv(summary_df, summary_file)

message("Summary saved to: ", summary_file)
print(summary_df)
