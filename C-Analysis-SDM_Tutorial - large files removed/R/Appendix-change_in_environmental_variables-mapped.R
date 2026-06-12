# =============================================================================
# Integrated MaxEnt Species Distribution Model
# Nordic Cloudberry (Rubus chamaemorus) Distribution Modeling
# Appendix - change in environmental variable mapped
# Created with claude.ai (sonnet 4.5) on 2026-01-27 by Nicole Drewitz
# edited with Perplexity.ai

# =============================================================================

# ============================================
# VARIABLES TO EXCLUDE (edit this list)
# ============================================
exclude_variables <- c(
  "soil_phh2o_0-5-30cm_mean_averaged.asc",
  "soil_soc_0-5-30cm_mean_averaged.asc",
  "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged.asc",
  "Topo_northness.asc"
)
# ============================================

# Load required libraries
library(raster)
library(RColorBrewer)
library(terra)  # Modern raster package

# Set working directory
main_folder <- "C:/Users/nidre4892/OneDrive - UiT Office 365/Documents/Temporary ASC files for QGIS/A_1_c - Environmental Variables - Final Model - 2026-01-15"
setwd(main_folder)

# Define scenario folders
current_folder <- "Environmental Variables - Final - 1 - Current"
ssp126_folder <- "Environmental Variables - Final - 2 - ssp126"
ssp370_folder <- "Environmental Variables - Final - 3 - ssp370"
ssp585_folder <- "Environmental Variables - Final - 4 - ssp585"

# Get list of ASC files from current scenario
current_files <- list.files(current_folder, pattern = "\\.asc$", full.names = TRUE)
variable_names <- basename(current_files)

# Filter out excluded variables
if (length(exclude_variables) > 0) {
  keep_indices <- !(variable_names %in% exclude_variables)
  current_files <- current_files[keep_indices]
  variable_names <- variable_names[keep_indices]

  cat("Excluded variables:", paste(exclude_variables, collapse = ", "), "\n")
  cat("Remaining variables:", length(variable_names), "\n\n")
}
variable_names

# Function to calculate change (future - current)
calculate_change <- function(current_file, future_folder, var_name) {
  future_file <- file.path(future_folder, var_name)
  if (!file.exists(future_file)) {
    warning(paste("Future file not found:", future_file))
    return(NULL)
  }
  current_raster <- rast(current_file)
  future_raster <- rast(future_file)
  # Calculate change
  change <- future_raster - current_raster
  return(change)
}

# Define scenarios
scenarios <- list(
  ssp126 = ssp126_folder,
  ssp370 = ssp370_folder,
  ssp585 = ssp585_folder
)

# Create custom color palette: blue (negative), middle (0), red (positive)
custom_palette <- colorRampPalette(c("#0571b0", "lightgrey", "#ca0020"))(101)

# Loop through each variable and scenario
cat("Displaying change maps...\n")

# Store all change rasters for later plotting
all_changes <- list()

for (i in seq_along(current_files)) {
  var_name <- variable_names[i]
  cat("\nProcessing variable:", var_name, "\n")

  for (scenario_name in names(scenarios)) {
    scenario_folder <- scenarios[[scenario_name]]
    cat("  Scenario:", scenario_name, "\n")

    # Calculate change
    change_raster <- calculate_change(current_files[i], scenario_folder, var_name)

    if (!is.null(change_raster)) {
      # Store for later use
      plot_key <- paste(var_name, scenario_name, sep = "_")
      all_changes[[plot_key]] <- list(
        raster = change_raster,
        var_name = var_name,
        scenario = scenario_name
      )

      # Set up plotting window
      par(mfrow = c(1, 1), mar = c(4, 4, 3, 6))

      # Get the range of values to center the color scale at zero
      raster_range <- range(values(change_raster), na.rm = TRUE)
      max_abs <- max(abs(raster_range))

      # Plot with custom palette centered at zero
      plot(change_raster,
           main = paste("Change in", gsub("\\.asc$", "", var_name),
                        "\n(", scenario_name, "- Current)"),
           col = custom_palette,  # Blue=negative, Pale yellow=0, Red=positive
           range = c(-max_abs, max_abs),  # Center at zero
           axes = TRUE,
           box = TRUE,
           cex.main = 1.2,
           cex.axis = 0.9)
    }
  }
}

cat("\n=== All maps displayed! ===\n")

# Automatically save as matrices of plots (3x2 landscape layout)
n_plots <- length(all_changes)
plots_per_file <- 6  # 2 rows x 3 columns (landscape)
n_files <- ceiling(n_plots / plots_per_file)

if (n_plots >= 1) {
  cat(paste0("\nSaving plots as 3x2 landscape matrices across ", n_files, " file(s)...\n"))

  plot_names <- names(all_changes)

  for (file_idx in 1:n_files) {
    start_idx <- (file_idx - 1) * plots_per_file + 1
    end_idx <- min(file_idx * plots_per_file, n_plots)

    filename <- paste0("variable_change_map_matrix", file_idx, ".png")
    cat(paste0("Creating file ", file_idx, " (plots ", start_idx, "-", end_idx, ")... "))

    png(filename, width = 3000, height = 2000, res = 200)  # Landscape: wide x short
    par(mfrow = c(2, 3), mar = c(3, 3, 2, 4))  # 2 rows x 3 columns

    for (j in start_idx:end_idx) {
      change_data <- all_changes[[plot_names[j]]]

      # Get the range of values to center the color scale at zero
      raster_range <- range(values(change_data$raster), na.rm = TRUE)
      max_abs <- max(abs(raster_range))

      plot(change_data$raster,
           main = paste(gsub("\\.asc$", "", change_data$var_name),
                        "-", change_data$scenario),
           col = custom_palette,
           range = c(-max_abs, max_abs),
           axes = TRUE,
           cex.main = 0.9,
           cex.axis = 0.7)
    }
    dev.off()
    cat("Saved as '", filename, "'\n", sep = "")
  }
} else {
  cat("\nNo plots to save.\n")
}

# Reset plotting parameters
par(mfrow = c(1, 1))
cat("\n=== Script complete! ===\n")
