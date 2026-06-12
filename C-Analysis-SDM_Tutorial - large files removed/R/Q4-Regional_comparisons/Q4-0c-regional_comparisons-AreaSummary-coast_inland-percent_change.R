# ============================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q4:
# Extracting area above thresholds summary with percent area from each category
# e.g. % habitat above threshold of total coastal area
# created by Nicole Drewitz on 21 Feb, 2026 with Claude AI
# ============================================================
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(here)

# ── Paths ────────────────────────────────────────────────────
asc_dir <- here::here("Input_files/MaxEnt_files")
VECTOR_PATH <- here::here("Input_files/Modelling_area-regional_divisions.gpkg")
output_data <- here::here("analysis_output_data")

# ── Settings ─────────────────────────────────────────────────
THRESHOLDS    <- c(0.5092, 0.75, 0.90)
FIELD         <- "CoastCNT" # Inland, ECO_Cor, CoastCNT
Category_NAME <- "Coasts_by_country" # used to name file. change when FIELD is changed. e.g. "Inland_coastal", "Ecoregions", Coasts_by_country"
CRS_LAEA      <- "EPSG:3035"

options(scipen = 999) # suppress scientific notation

# ── Load rasters ─────────────────────────────────────────────
message("Loading rasters...")
current <- rast(file.path(asc_dir, "Rubus_chamaemorus-Scenario_1 - 1981-2010.asc"))
ssp126 <- rast(file.path(asc_dir, "Rubus_chamaemorus-Scenario_2 - ssp126.asc"))
ssp370 <- rast(file.path(asc_dir, "Rubus_chamaemorus-Scenario_3 - ssp370.asc"))
ssp585 <- rast(file.path(asc_dir, "Rubus_chamaemorus-Scenario_4 - ssp585.asc"))

crs(current) <- CRS_LAEA
crs(ssp126)  <- CRS_LAEA
crs(ssp370)  <- CRS_LAEA
crs(ssp585)  <- CRS_LAEA

# ── Load & reproject vector ───────────────────────────────────
message("Loading vector layer...")
regions <- st_read(VECTOR_PATH, quiet = TRUE) |>
  st_transform(CRS_LAEA)

if (!FIELD %in% names(regions))
  stop("Field '", FIELD, "' not found. Available: ", paste(names(regions), collapse = ", "))

message("Unique values in '", FIELD, "': ",
        paste(sort(unique(regions[[FIELD]])), collapse = ", "))

# ── Helper: standardize category labels ──────────────────────
label_category <- function(x) {
  trimws(as.character(x))
}

# ── Helper: total cells per category (denominator) ───────────
# Mask the region raster by the suitability raster so only cells with a
# valid suitability value are counted — ensures all bands sum to 100%.
get_total_cells <- function(r_suit, regions, field) {
  region_r <- rasterize(vect(regions), r_suit, field = field)
  region_r <- mask(region_r, r_suit)
  as.data.frame(region_r, na.rm = TRUE) |>
    setNames("category_raw") |>
    mutate(category = label_category(category_raw)) |>
    filter(!is.na(category), category != "") |>
    count(category, name = "total_cells")
}

# ── Helper: count cells per band ─────────────────────────────
# Bands (non-overlapping):
#   "poor"  : suitability <  THRESHOLDS[1]
#   band i  : suitability >= THRESHOLDS[i] and < THRESHOLDS[i+1]
#   top band: suitability >= THRESHOLDS[n]
count_bands <- function(r, regions, field) {
  region_r <- rasterize(vect(regions), r, field = field)
  n        <- length(THRESHOLDS)

  extract_band <- function(mask_raster, threshold_val) {
    stk <- c(mask_raster, region_r)
    names(stk) <- c("in_band", "category_raw")
    as.data.frame(stk, na.rm = TRUE) |>
      filter(in_band == 1) |>
      mutate(category = label_category(category_raw)) |>
      filter(!is.na(category), category != "") |>
      group_by(category) |>
      summarise(area_km2 = n(), .groups = "drop") |>
      mutate(threshold = threshold_val)
  }

  # Poor habitat band (below lowest threshold)
  poor <- extract_band(r < THRESHOLDS[1], threshold_val = "poor")

  # Suitability bands
  suit <- imap_dfr(THRESHOLDS, function(lo, i) {
    hi <- if (i < n) THRESHOLDS[i + 1] else Inf
    extract_band((r >= lo) & (r < hi), threshold_val = as.character(lo))
  })

  all_combos <- expand.grid(
    category  = all_categories,
    threshold = c("poor", as.character(THRESHOLDS)),
    stringsAsFactors = FALSE
  )

  bind_rows(poor, suit) |>
    right_join(all_combos, by = c("category", "threshold")) |>
    mutate(area_km2 = replace_na(area_km2, 0)) # no cells with values = no suitable habitat at 1km2 resolution
}

# ── Compute denominator ───────────────────────────────────────
message("Computing total cells per category...")
total_cells <- get_total_cells(current, regions, FIELD)
all_categories <- total_cells$category
message("  ", paste(apply(total_cells, 1, paste, collapse = ": "), collapse = " | "))

# ── Scenarios ─────────────────────────────────────────────────
scenarios <- list(
  current = current,
  ssp126  = ssp126,
  ssp370  = ssp370,
  ssp585  = ssp585
)

# ── Extract counts ────────────────────────────────────────────
message("Counting cells per scenario / category / band...")
results <- imap_dfr(scenarios, function(r, nm) {
  message("  ", nm)
  count_bands(r, regions, FIELD) |> mutate(scenario = nm)
})

# ── Normalize ─────────────────────────────────────────────────
results <- results |>
  left_join(total_cells, by = "category") |>
  mutate(pct_area = round((area_km2 / total_cells) * 100, 4)) |>
  dplyr::select(-total_cells)

message("\nPreview:")
print(results |> arrange(scenario, category, threshold))
#print(results, n = Inf)

# ── Save ──────────────────────────────────────────────────────
out_csv <- file.path(output_data, paste0("Q4-0c-area_data-", Category_NAME, "-percent_change_total.csv"))
write.csv(results, out_csv, row.names = FALSE)
message("Saved: ", out_csv)
