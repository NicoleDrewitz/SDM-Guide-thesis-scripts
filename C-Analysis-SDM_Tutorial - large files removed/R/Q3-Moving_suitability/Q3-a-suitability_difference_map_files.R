# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q3
# created rasters of the difference in suitability between current and future climates
# Create with Claude.ai on 27-01-2026 by Nicole Drewitz

# =============================================================================

# =============================================================================
# WORKSPACE SETUP
# =============================================================================
# Load required library
library(raster)
library(here)

# Set working directory to your folder
setwd(base_path <- here::here("Input_files/Environmental_variables"))

# Load the rasters
current <- raster("Rubus_chamaemorus - current.asc")
ssp126 <- raster("Rubus_chamaemorus_Environmental Variables - Final - 2 - ssp126.asc")
ssp370 <- raster("Rubus_chamaemorus_Environmental Variables - Final - 3 - ssp370.asc")
ssp585 <- raster("Rubus_chamaemorus_Environmental Variables - Final - 4 - ssp585.asc")

# Calculate change (future - current) for each scenario
change_ssp126 <- ssp126 - current
change_ssp370 <- ssp370 - current
change_ssp585 <- ssp585 - current

# Write the change rasters to ASC files
writeRaster(change_ssp126, "Rubus_chamaemorus_change_ssp126.asc",
            format = "ascii", overwrite = TRUE)
writeRaster(change_ssp370, "Rubus_chamaemorus_change_ssp370.asc",
            format = "ascii", overwrite = TRUE)
writeRaster(change_ssp585, "Rubus_chamaemorus_change_ssp585.asc",
            format = "ascii", overwrite = TRUE)

cat("\nChange rasters created successfully!\n")

# =============================================================================
# SUITABILITY AREA CHANGE PLOT
# =============================================================================

# Load required libraries for plotting
library(ggplot2)
library(raster)  # for cell area if needed

# Get cell area in km² (assuming projected CRS; adjust if lat/long)
cell_area <- area(current)[1,1]  # km² per cell, assuming uniform

# Function to compute cumulative area (km²) up to each quantile (= habitat suitability range)
compute_cum_area <- function(rast) {
  vals <- values(rast)
  vals <- sort(vals[!is.na(vals)])
  n <- length(vals)
  quants <- quantile(vals, c(0.05, 0.25, 0.50, 0.75, 0.95))
  cum_cells <- numeric(5)
  for(i in 1:5) {
    cum_cells[i] <- sum(vals <= quants[i])
  }
  return(cum_cells * cell_area)
}

# Current baseline areas
current_areas <- compute_cum_area(current)

# Future areas
ssp126_areas <- compute_cum_area(ssp126)
ssp370_areas <- compute_cum_area(ssp370)
ssp585_areas <- compute_cum_area(ssp585)

# Change in area (km²)
change_areas <- rbind(
  Current = 0,
  ssp126 = ssp126_areas - current_areas,
  ssp370 = ssp370_areas - current_areas,
  ssp585 = ssp585_areas - current_areas
)

# Data frame for plotting (stacked changes)
quants_data <- data.frame(
  scenario = rep(rownames(change_areas), each = 5),
  quant = rep(c("5th", "25th", "50th", "75th", "95th"), nrow(change_areas)),
  area_km2 = as.vector(t(change_areas))
)

# Print for check
cat("Change in cumulative area (km²) up to each quantile:\n")
print(quants_data)

# Plot
p <- ggplot(quants_data, aes(x = scenario, y = area_km2, fill = quant)) +
  geom_area(position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("5th" = "#2166AC", "25th" = "#67A9CF", "50th" = "#F7F7F7",
                               "75th" = "#FDAE61", "95th" = "darkred"), guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Climate scenario", y = "Area (km²)",
       title = "Change in habitat suitability range") +
  annotate("text", x = 1.5, y = min(quants_data$area_km2)*0.95,
           label = "Climatic\nbaseline", size = 3.5, hjust = 0.5) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))

print(p)

cat("\nPlot created! Y-axis: stacked change in km² up to each suitability quantile vs current.\n")
