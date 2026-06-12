#===============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q2:
# Spatial autocorrelation check at sampling locations
# Created by Nicole Drewitz with Claude.ai on 18 Feb, 2026
#===============================================================================

# ── Spatial Autocorrelation Check for Cloudberry Presence/Absence Data ──────────────
# CRS: LAEA EPSG:3035 — coordinates in meters

library(ncf)
library(here)
fig_dir  <- here::here("figures")

#── 1. Load presence and absence data ──────────────────────────────────────────────
absence_file_select  <- read.csv(here::here("Input_files/V2_cloudberry_absence_analysis.csv"),
                          stringsAsFactors = FALSE)
presence_file_select <- read.csv(here::here("Input_files/V2_cloudberry_presence_analysis.csv"),
                          stringsAsFactors = FALSE)

cloudberry_presence <- presence_file_select
cloudberry_absence <- absence_file_select

# Disable scientific notation globally (affects plots, console, printing)
options(scipen = 999)

# Combine datasets
cloudberry <- rbind(
  transform(cloudberry_presence, presence = 1),
  transform(cloudberry_absence, presence = 0)
)

#── 1b. Exclude Greenland and Svalbard ─────────────────────────────────────────────
# LAEA EPSG:3035 approximate bounding box for mainland Europe + Iceland
# Greenland sits west of ~-1,500,000 m (LAEA x), Svalbard north of ~5,200,000 m (LAEA y)
# Adjust these thresholds based on your actual coordinate ranges printed below
cat("Longitude (x) range BEFORE filter:", range(cloudberry$longitude), "\n")
cat("Latitude  (y) range BEFORE filter:", range(cloudberry$latitude),  "\n")

cloudberry <- cloudberry[
  cloudberry$longitude > -1500000 &   # exclude Greenland (far west)
    cloudberry$latitude  <  5600000,    # exclude Svalbard (far north)
]

cat("Records after geographic filter:", nrow(cloudberry), "\n")
cat("Longitude (x) range AFTER filter:", range(cloudberry$longitude), "\n")
cat("Latitude  (y) range AFTER filter:", range(cloudberry$latitude),  "\n")
head(cloudberry)
str(cloudberry)
cat("Total records:", nrow(cloudberry), "\n")
cat("Presence records:", sum(cloudberry$presence), "\n")
cat("Absence records:", sum(1 - cloudberry$presence), "\n")
cat("Longitude range:", range(cloudberry$longitude), "\n")
cat("Latitude range:", range(cloudberry$latitude), "\n")

#write.csv(cloudberry, here::here("analysis_output_data/V2_cloudberry_clipped_mainland.csv"), row.names = FALSE)

#── 2. Run spline correlogram ─────────────────────────────────────────────────────
# Study area spans ~4,000 km west-east (Greenland to Finland) and
# ~3,000 km north-south. xmax of 100,000 m (100 km) captures ecologically
# relevant autocorrelation distances.
# resamp = 99 for quick check; use 999 for publication.

set.seed(42)
Correlog <- spline.correlog(
  x = cloudberry[, "longitude"],
  y = cloudberry[, "latitude"],
  z = cloudberry[, "presence"],
  xmax = 100000,  # 100 km in metres
  resamp = 99 # bootstrap resampling
)

#── 3. Plot (NO SCIENTIFIC NOTATION) ──────────────────────────────────────────────
plot(Correlog,
     main = "Spline Correlogram\nCloudberry model selection presence/absence locations",
     xlab = "Distance (km)",  # Changed label
     ylab = "Moran's I",
     xaxt = "n")  # Suppress default x-axis

# Custom x-axis in kilometers
x_ticks_m <- pretty(par("usr")[1:2], 10)  # 10 = axis ticks
x_ticks_km <- x_ticks_m / 1000           # Convert to km
axis(1, at = x_ticks_m, labels = format(round(x_ticks_km), big.mark = ","))  # Use original positions, km labels

abline(h = 0, lty = 2, col = "grey50") # 95% confidence interval from resampling

#── 4. Extract thinning distance ───────────────────────────────────────────────────
cat("\n--- Suggested minimum thinning distance ---\n")
intercept_km <- round(Correlog$x.intercept / 1000, 1)
cat("x-intercept (autocorrelation ~ 0):", formatC(Correlog$x.intercept, format = "d", big.mark = ","), "metres\n")
cat(" =", formatC(intercept_km, format = "d", big.mark = ","), "km\n")

#===============================================================================
# NEAREST NEIGHBOUR DISTANCE IN MODELLING PRESENCES
#===============================================================================
library(sf)

# Load and convert to sf object
# Use filtered presence points only
cloudberry_presence_filtered <- cloudberry[cloudberry$presence == 1, ]

pts <- st_as_sf(cloudberry_presence_filtered, coords = c("longitude", "latitude"), crs = 3035)

summary(nn_min_km)
hist(nn_min_km, breaks = 50)  # no xlim first, to see the full range

# Calculate nearest-neighbour distances
nn_dists <- st_distance(pts)
diag(nn_dists) <- NA
nn_min <- apply(nn_dists, 1, min, na.rm = TRUE)
nn_min_km <- as.numeric(nn_min) / 1000

#===============================================================================
# PLOT
#===============================================================================
# Plot filtered distribution
hist(nn_min_km,
     breaks = 50,
     main = "Nearest Neighbour Distances in Modelling Locations",
     xlab = "Distance to nearest neighbour (km)",
     xlim = c(0, 30),
     col = "steelblue", border = "white")
axis(1, at = seq(0, 30, by = 5)) # add ticks "by" every...
abline(v = 1, lty = 2, col = "black", lwd = 2)
abline(v = median(nn_min_km), lty = 2, col = "red", lwd = 2)
abline(v = mean(nn_min_km), lty = 2, col = "orange", lwd = 2)

legend("top",
       legend = c("1 km resolution",
                  paste0("Median: ", round(median(nn_min_km), 1), " km"),
                  paste0("Mean: ", round(mean(nn_min_km), 1), " km")),  # ← ) closes c()
       col = c("black", "red", "orange"),
       lty = 2, lwd = 2)
