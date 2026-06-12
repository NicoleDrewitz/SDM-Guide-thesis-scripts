# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q4: Point Generation for Regional Comparisons - Data Extraction
# Random Point Generation per Polygon - Excluding Modelled Raster Pixels

# Purpose: Generate points for MaxEnt regional comparisons by randomly placing
# points within raster cells that have suitability predictions, but
# do not already have another point within the cell. Min distance between points can be set.
# Set max No. of points for each polygon category (Inland/Coast/ECO_cor/EcoCoast divisions)

# Created: 19 Feb 2026 by Nicole Drewitz (with Claude AI and Perplexity AI for editing)
# Last updated: 20 Feb 2026
# =============================================================================

# ---- CORE LIBRARIES ----
library(sf)        # Vector polygon operations (st_read, st_union)
library(terra)     # Raster processing (rast, crop, mask, xyFromCell)
library(dplyr)     # Data manipulation (bind_rows)
library(here)      # Project-relative paths

# ---- PROJECT PATHS ----
input_dir <- here::here("Data/Q4-random_sampling_locations")

# ---- USER CONFIGURATION ----
FIELD_NAME   <- "Inland" # Inland, CoastCNT, ECO_Cor, EcoCoast
N_POINTS     <- 1000
MIN_DIST_M   <- 10000       # Minimum distance between points (meters)
# field name needs to be hardcoded on line 162!!!!!!!!!!!!
    # " data.frame(
    # CoastCNT = poly_id,
    #x = coords[kept_idx, 1],
    #y = coords[kept_idx, 2]) "

# ANALYSIS CONFIGURATIONS:
# FIELD_NAME              | Inland | CoastCNT                     | ECO_Cor                                 | EcoCoast
# N_POINTS (per category) | 1000   | 200                          | 300                                     | 150
# MIN_DIST_M              | 10000  | 2000 = 2km                   | 5000 = 5 km                             | 2000
# no. categories          | 2      | 29 coasts + 7 inland         | 10                                      | 66
# raster contents         | binary | coast area/inland by country | ecoregions in modelling area            | ECO_cor + CoastCNT
# notes                   | below-1| analyze coastal areas only   | https://ecoregions.appspot.com/         | combined analysis (not conducted for thesis)
#                         | 1-Greenland/Svalbard points removed with QGIS after to eliminate over-representation in suitability values where data is limited. points increased to account for removal.
#                         |

SET_SEED     <- 42
NODATA_VAL   <- -9999
CRS_CODE     <- "EPSG:3035" # prevents false CRS mismatch warnings

# ── Optional geographic exclusion ──────────────────────────────────────────
# Uncomment to exclude regions by LAEA (EPSG:3035) coordinate bounds.
# Check printed ranges below to calibrate thresholds.
# Common use: remove Greenland (far west) and/or Svalbard (far north)
# to avoid over-representation where occurrence data is sparse. e.g. Svalbard and Greenland are removed from coastal/inland comparison because they would over-represent extreme coastal habitat.

# ── Optional country exclusion ─────────────────────────────────
EXCL_ISO3 <- c("GRL", "SJM")  # Greenland, Svalbard
# ───────────────────────────────────────────────────────────────────────────

# ---- INPUT FILES ----
VECTOR_PATH  <- here::here("Input_files/Modelling_area-regional_divisions.gpkg")

RASTER_PATH  <- here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_1 - 1981-2010.asc")

OUTPUT_CSV   <- file.path(input_dir, paste0("Q4-0a-random_points_for_sampling-", FIELD_NAME, "-n", N_POINTS, "each-mindist", MIN_DIST_M, "m.csv"))

# ---- STEP 1: LOAD SPATIAL DATA ----
cat("Loading vector...\n")
polys <- st_read(VECTOR_PATH, quiet = TRUE)

cat("Column names:\n") # check file and available fields
print(colnames(polys))
cat("\nUnique values in ECO_cor (or similar):\n")
if ("ECO_cor" %in% colnames(polys)) {
  print(unique(polys$ECO_cor))
  cat("Any non-NA values:", sum(!is.na(polys$ECO_cor)), "\n")
} else {
  cat("Column 'ECO_cor' not found. Similar names:\n")
  print(colnames(polys)[grep("eco|cor", colnames(polys), ignore.case = TRUE)])
}
print(st_layers(VECTOR_PATH))  # Check GPKG layers

cat("Loading raster...\n")
ras <- rast(RASTER_PATH)

NAflag(ras) <- NODATA_VAL
ras[ras == NODATA_VAL] <- NA
crs(ras) <- CRS_CODE

# ---- STEP 2: CREATE OCCUPIED PIXEL MASK ----
cat("🔄 Building occupied pixel mask...\n")
occupied <- !is.na(ras)
cat(sprintf("   Free cells available: %.0f (%.1f%% of total)\n",
            sum(!values(occupied, na.rm = TRUE)),
            100 * sum(!values(occupied, na.rm = TRUE)) / ncell(occupied)))

# ---- STEP 3: GREEDY MINIMUM-DISTANCE FILTER ----
# Retains points in random order, dropping any candidate closer than
# MIN_DIST_M to an already-accepted point. O(n²) but fine for n ≤ few thousand.

filter_min_distance <- function(coords_matrix, min_dist) {
  if (nrow(coords_matrix) == 0) return(integer(0))

  kept <- integer(nrow(coords_matrix))  # pre-allocate index store
  n_kept <- 0L

  for (i in seq_len(nrow(coords_matrix))) {
    if (n_kept == 0L) {
      n_kept <- n_kept + 1L
      kept[n_kept] <- i
    } else {
      # Euclidean distances to all accepted points (CRS is metric LAEA)
      dx <- coords_matrix[kept[seq_len(n_kept)], 1] - coords_matrix[i, 1]
      dy <- coords_matrix[kept[seq_len(n_kept)], 2] - coords_matrix[i, 2]
      dists <- sqrt(dx^2 + dy^2)
      if (all(dists >= min_dist)) {
        n_kept <- n_kept + 1L
        kept[n_kept] <- i
      }
    }
  }
  kept[seq_len(n_kept)]
}

# ---- STEP 4: POINT SAMPLING FUNCTION ----
# Samples all free cells, shuffles them, then greedily keeps spatially spread points.

sample_polygon <- function(poly, poly_id, n_max, occupied_mask, min_dist) {
  poly_vect  <- vect(poly)
  poly_vect  <- project(poly_vect, CRS_CODE)
  ras_crop   <- crop(occupied_mask, poly_vect)
  ras_mask   <- mask(ras_crop, poly_vect)

  free_cells <- which(values(ras_mask) == TRUE)
  coords <- xyFromCell(ras_mask, free_cells)

  if (length(free_cells) == 0) {
    warning(sprintf("⚠️  No free cells in polygon '%s' - skipping", poly_id))
    return(NULL)
  }

  # Shuffle all candidates first so greedy selection isn't spatially biased
  free_cells <- sample(free_cells)
  coords     <- xyFromCell(ras_mask, free_cells)   # all candidate coordinates
  
  poly_sub <- polys[polys[[FIELD_NAME]] == uid, ]
  
  cat(sprintf("      %d free cells → applying %.0fm minimum distance filter...\n",
              nrow(coords), min_dist))

  kept_idx <- filter_min_distance(coords, min_dist)

  if (length(kept_idx) == 0) {
    warning(sprintf("⚠️  No points survived distance filter in '%s' - try reducing MIN_DIST_M", poly_id))
    return(NULL)
  }

  # Honour n_max cap after distance filtering
  if (length(kept_idx) > n_max) kept_idx <- kept_idx[seq_len(n_max)]

  cat(sprintf("      → kept %d points (cap: %d)\n", length(kept_idx), n_max))

  data.frame(
    Inland = poly_id, # CHANGE FIELD NAME HERE!!! (e.g. ECO_Cor, Inland, CoastCNT, ECO_Cor, EcoCoast)
    x = coords[kept_idx, 1],
    y = coords[kept_idx, 2]
  )
}

# ---- STEP 5: GENERATE POINTS PER POLYGON ----
cat("🔄 Sampling random points per", FIELD_NAME, "category...\n")
unique_ids <- unique(polys[[FIELD_NAME]])

if (length(unique_ids) == 0) {
  stop("❌ No unique values found in FIELD_NAME column")
} else {
  cat(sprintf("   Processing %d unique %s categories\n", length(unique_ids), FIELD_NAME))

  set.seed(SET_SEED)
  results <- list()

  for (i in seq_along(unique_ids)) {
    uid <- unique_ids[i]
    cat(sprintf("   [%d/%d] Processing: %s\n", i, length(unique_ids), uid))

    poly_sub   <- polys[polys[[FIELD_NAME]] == uid, ]
    
    # ── Apply country exclusion before union ──────────────────────────────
    if (exists("EXCL_ISO3")) {
      n_before <- nrow(poly_sub)
      poly_sub <- poly_sub[!poly_sub$ISO3_CODE %in% EXCL_ISO3, ]
      cat(sprintf("      Country filter removed %d polygons (%d remain)\n",
                  n_before - nrow(poly_sub), nrow(poly_sub)))
    }
    # ─────────────────────────────────────────────────────────────────────
    
    poly_union <- st_union(poly_sub)
    poly_sf    <- st_sf(geometry = poly_union)

    pts <- tryCatch(
      sample_polygon(poly_sf, uid, N_POINTS, occupied, MIN_DIST_M),
      error = function(e) {
        warning(sprintf("❌ Error in %s: %s", uid, e$message))
        return(NULL)
      }
    )

    if (!is.null(pts)) results[[as.character(uid)]] <- pts
  }
}

# ---- STEP 6: COMBINE AND VALIDATE OUTPUT ----
cat("🔄 Combining and validating results...\n")
all_points <- bind_rows(results)

if (nrow(all_points) == 0) {
  stop("❌ No points generated - check polygon-raster overlap, free cells, and MIN_DIST_M")
}

cat(sprintf("✅ SUCCESS: Generated %d points total\n", nrow(all_points)))
cat(sprintf("   Points per %s category:\n", FIELD_NAME))
print(table(all_points[[FIELD_NAME]]))
cat(sprintf("   Points per category: min=%d, max=%d, mean=%.0f\n",
            min(table(all_points[[FIELD_NAME]])),
            max(table(all_points[[FIELD_NAME]])),
            mean(table(all_points[[FIELD_NAME]]))))

# ---- STEP 7: EXPORT ----
write.csv(all_points, OUTPUT_CSV, row.names = FALSE)
cat(sprintf("💾 Saved %d points to: %s\n", nrow(all_points), OUTPUT_CSV))
cat("\n📍 Next steps:\n")
cat("   1. Extract suitability values with these points using Q4-0b-data_extraction-regional_comparisons-sampling.R\n")
cat("   2. Plot MaxEnt regional comparisons with Q4-a-....R scripts\n")
cat("   3. Save and discussion mean relevent change\n\n")
