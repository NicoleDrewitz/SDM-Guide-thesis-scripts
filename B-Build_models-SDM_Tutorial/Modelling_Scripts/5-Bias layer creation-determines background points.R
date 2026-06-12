# =============================================================================
# Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# MaxEnt Bias Layer Creation Script
# This script creates a bias layer based on background point sampling
# Created by Nicole Drewitz with Claude.ai (Sonnet 4.5) on 2026-01-15

# Background points = pseudo-absences (i.e. random locations to be compared with presences)
# This bias layer matches the background sampling used in the variable and model selection process.
# This script does not allow background points to be selected within the same cells as presence locations used in modelling.
# It also encourages more background points to be sampled near presence cells (target-group background sampling)
# =============================================================================

# paste this script after generating your background points from your modelling process script.
# Only change made was to increased the number of background points
#         (20 000 possible points increased to 82 000 points) to allow the possible background points selected to be increased from 10 000 to 50 000.

# =============================================================================
# CREATE BIAS LAYER
# =============================================================================

# Improved bias layer creation (replace your existing CREATE BIAS LAYER section)
cat("\n=== Creating Improved Bias Layer ===\n")

# Get template raster
template_raster <- initial_env_stack[[1]]

# Rasterize background points to create density surface
point_density <- rasterize(background_points_clean, template_raster, fun = "count", background = 0)

# Smooth with focal window (mean)
bias_raw <- focal(point_density, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)

# Robust normalization: add tiny epsilon first, then scale to 0-1
epsilon <- 1e-10
bias_layer <- bias_raw + epsilon
bias_min <- cellStats(bias_layer, min)
bias_max <- cellStats(bias_layer, max)
bias_layer <- (bias_layer - bias_min) / (bias_max - bias_min)

# Final check: force all values > 0.001, cap at 1
bias_layer[bias_layer < 0.001] <- 0.001
bias_layer[bias_layer > 1] <- 1
bias_layer[is.na(bias_layer)] <- 0.001

# Verify
cat("Min value:", cellStats(bias_layer, min), "\n")
cat("Max value:", cellStats(bias_layer, max), "\n")
cat("Zeros present?", any(getValues(bias_layer) <= 0), "\n")  # Extracts values first

cat("Improved bias layer created successfully\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

# Save bias layer as .asc for MaxEnt
output_asc <- file.path(base_dir, "maxent_bias_layer.asc")
writeRaster(bias_layer, output_asc, overwrite = TRUE, format = "ascii")
cat("\nBias layer saved to:", output_asc, "\n")

# Save background points
bg_output <- file.path(base_dir, "background_points.csv")
write.csv(data.frame(longitude = background_points_clean[,1], 
                     latitude = background_points_clean[,2]),
          bg_output, row.names = FALSE)
cat("Background points saved to:", bg_output, "\n")

# =============================================================================
# CREATE VISUALIZATION
# =============================================================================

cat("\n=== Creating Visualization ===\n")

# Optional save
save_plot <- TRUE  # Set to FALSE to skip saving
if (save_plot) {
  png(file.path(base_dir, "bias_layer_plot.png"), width = 1000, height = 600)
}
par(mfrow = c(1, 2))

# Plot 1: Bias layer only (sampling density)
my_colors <- colorRampPalette(c("white", "blue", "black"))(100)  # Option 3: Blue-White-Red
plot(bias_layer, main = "MaxEnt Bias Layer\n(Sampling Density)", 
     col = my_colors)

# Plot 2: Presence and background points only (no raster/colors)
plot(0, type = "n", xlim = par("usr")[1:2], ylim = par("usr")[3:4],
     main = "Presence and Background Points", xlab = "", ylab = "", axes = FALSE)
points(background_points_clean, pch = ".", cex = 0.3, col = "#7986CB")
points(clean_occurrence_points, pch = 20, cex = 0.3, col = "#E57373")
legend("topright", 
       legend = c("Presence (creating model)", "Background"),
       col = c("#E57373", "#7986CB"), 
       pch = c(20, 20),
       cex = 0.8,
       pt.cex = 1.5,
       bg = "white")

if (save_plot) {
  dev.off()
  cat("Visualization saved to: bias_layer_plot.png\n")
} else {
  cat("Visualization printed to screen (no save).\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat("Background points generated:", nrow(background_points_clean), "\n")
cat("Bias layer range:", round(cellStats(bias_layer, min), 4), "to", round(cellStats(bias_layer, max), 4), "\n")
cat("\nTo use in MaxEnt: Specify 'maxent_bias_layer.asc' as the bias file parameter\n")
cat("\nDone!\n")