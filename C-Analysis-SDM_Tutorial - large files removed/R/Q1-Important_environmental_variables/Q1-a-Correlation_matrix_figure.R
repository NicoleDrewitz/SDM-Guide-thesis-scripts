# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q1

# Environmental Variables Correlation Heatmap
# Created by Nicole Drewitz with Claude.ai (4.5 Sonnet) on 2026-01-23
# =============================================================================

# =============================================================================
# CONFIGURATION SECTION - MODIFY THESE SETTINGS
# =============================================================================

fig_dir  <- here::here("figures")
output_data <- here::here("analysis_output_data")
save_plots <- TRUE

CONFIG <- list(
  # Number of cells to sample for correlation analysis
  sample_size = 50000,
  # Random seed for reproducibility
  random_seed = 123
)

# Directory paths
PATHS <- list(
  env_raster_dir = here::here("Input_files/Environmental_variables") 
  #env_raster_dir = "/Users/nicoledrewitz/Library/CloudStorage/OneDrive-Menntaský/Thesis Research Files - Mac Connection/A - Thesis Research Files - Mac Connection/1 - GIS files/c - Environmental variables - 2025-08-21/A_1_c- Environmental variables - modelled area/A_1_c - Environmental Variables - Final Model - 2026-01-15/Environmental Variables - Final - 1 - Current"
  )

# Variables to INCLUDE in analysis (list only the variables you want)==========
# Leave empty c() to include all variables found in directory

INCLUDE_VARIABLES <- c("bio11",
                       "bio13",
                       "bio19",
                       "fcf",
                       "swe",
                       "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged",
                       "soil_phh2o_0-5-30cm_mean_averaged",
                       "soil_soc_0-5-30cm_mean_averaged",
                       "Topo_northness")
####"kg_0", excluded because it is categorical

# =============================================================================
# LOAD REQUIRED PACKAGES
# =============================================================================

cat("=== Loading Required Packages ===\n")
required_packages <- c("corrplot", "terra", "ggplot2", "here", "usdm")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# DISPLAY CONFIGURATION
# =============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("=== CONFIGURATION SETTINGS ===\n")
cat(rep("=", 80), "\n\n", sep = "")
cat("Sample size:              ", CONFIG$sample_size, "\n")
cat("Random seed:              ", CONFIG$random_seed, "\n")
cat("Correlation method:       ", CONFIG$correlation_method, "\n")
cat("\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Shorten variable names at underscore for plotting
shorten_var_names <- function(var_names) {
  sapply(var_names, function(name) {
    parts <- strsplit(name, "_")[[1]]
    if (length(parts) > 1) {
      paste(parts[1], parts[2], sep = "_")
    } else {
      name
    }
  }, USE.NAMES = FALSE)
}

# =============================================================================
# LOAD ENVIRONMENTAL VARIABLES
# =============================================================================

cat("=== Loading Environmental Variables ===\n")
raster_files <- list.files(PATHS$env_raster_dir, pattern = "\\.asc$",
                           full.names = TRUE, ignore.case = TRUE)
if (length(raster_files) == 0) {
  stop("No .asc files found in environmental variables directory!")
}

# Load all rasters
env_raster <- rast(raster_files)

# Assign names
initial_var_names <- tools::file_path_sans_ext(basename(raster_files))
names(env_raster) <- initial_var_names

# Filter to include only specified variables (if any specified)
if (length(INCLUDE_VARIABLES) > 0) {
  # Check which requested variables exist
  available_vars <- INCLUDE_VARIABLES[INCLUDE_VARIABLES %in% names(env_raster)]
  missing_vars <- setdiff(INCLUDE_VARIABLES, available_vars)

  if (length(missing_vars) > 0) {
    cat("WARNING: These requested variables were not found:\n")
    cat(paste("  -", missing_vars, collapse = "\n"), "\n\n")
  }

  if (length(available_vars) == 0) {
    stop("None of the specified variables were found in the directory!")
  }

  env_raster <- env_raster[[available_vars]]
  cat("Using", nlyr(env_raster), "specified environmental layers\n")
} else {
  cat("Using all", nlyr(env_raster), "environmental layers found\n")
}

cat("Variables included:\n", paste("  -", names(env_raster), collapse = "\n"), "\n\n")

set.seed(CONFIG$random_seed)
# sample raster cells (rows = cells, cols = variables)
env_sample <- terra::spatSample(
  env_raster,
  size   = CONFIG$sample_size,
  method = "random",
  na.rm  = TRUE,
  as.data.frame = TRUE
)

# drop coordinates if present
env_sample <- env_sample[ , names(env_raster), drop = FALSE]

# =============================================================================
# SAMPLE RASTER CELLS FOR CORRELATION ANALYSIS
# =============================================================================
cat("=== Sampling Raster Cells for Correlation Analysis ===\n")
# Determine correlation method
if (is.null(CONFIG$correlation_method)) {
  normality_test <- apply(env_sample, 2, function(x) {
    if(length(x) > 3) {
      test_data <- if(length(x) > 5000) sample(x, 1000) else x
      shapiro.test(test_data)$p.value > 0.05
    } else TRUE
  })
  cor_method <- if (all(normality_test)) "pearson" else "spearman"
  cat("Auto-detected correlation method:", cor_method, "\n")
} else {
  cor_method <- CONFIG$correlation_method
  cat("Using specified correlation method:", cor_method, "\n")

  # env_sample is your data.frame of numeric columns
  # Choose layout, e.g. 3 rows x 3 cols for 9 variables
  par(mfrow = c(3, 3),        # matrix of panels
      mar   = c(4, 4, 2, 1))  # margins: bottom, left, top, right
  invisible(lapply(names(env_sample), function(v) {
    boxplot(
      env_sample[[v]],
      main = v,              # column name as title
      ylab = "",             # or a nicer label if you want
      col  = "lightgray",
      border = "black"
    )
  }))
}

# =============================================================================
# CALCULATE AND PLOT CORRELATION MATRIX
# =============================================================================
cat("=== Calculating Correlation Matrix ===\n")
cor_matrix <- cor(env_sample, use = "complete.obs", method = cor_method)
# Calculate total pairs
n_vars <- ncol(cor_matrix)
total_pairs <- (n_vars * (n_vars - 1)) / 2
cat("Total variable pairs evaluated:", total_pairs, "\n\n")

cat("=== Generating Correlation Heatmap ===\n")
# Create shortened names for plotting
short_names <- shorten_var_names(rownames(cor_matrix))
cor_matrix_plot <- cor_matrix
rownames(cor_matrix_plot) <- short_names
colnames(cor_matrix_plot) <- short_names

# =============================================================================
# PLOT WITH CORRELATION VALUES AND SCATTERPLOTS
# =============================================================================
# file path
corr_plot_file <- file.path(fig_dir, "Q1-a-corr_plot.png")
if (save_plots) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  png(filename = corr_plot_file, width = 12, height = 10, units = "in", res = 400)

  # Use corrplot with correlation coefficients added
  corrplot(
    cor_matrix_plot,
    method = "color",
    type = "upper",
    tl.col = "black",
    tl.srt = 45,
    tl.cex = 0.8,
    col = colorRampPalette(c("#0571b0", "white", "#ca0020"))(200),
    cl.cex = 0.8,
    title = paste0("Correlation Matrix (", cor_method, ")"),
    mar = c(0, 0, 2, 0),
    addCoef.col = "black",  # Add correlation coefficients in black
    number.cex = 1,       # Size of correlation numbers
    number.digits = 2,      # Number of decimal places
    diag = FALSE            # Don't show diagonal
  )

  dev.off()
  cat("✓ Saved: corr_plot.png (with correlation values)\n")
}

# =============================================================================
# CREATE FULL SCATTERPLOT MATRIX WITH CORRELATION HEATMAP
# =============================================================================
if (save_plots) {
  scatter_plot_file <- file.path(fig_dir, "Q1-a-corr_scatterplot_matrix.png")
  png(filename = scatter_plot_file, width = 14, height = 14, units = "in", res = 400)

  # Rename columns for the scatterplot matrix
  env_sample_plot <- env_sample
  colnames(env_sample_plot) <- short_names

  # Create scatterplot matrix with correlation values
  pairs(env_sample_plot,
        lower.panel = function(x, y, ...) {
          points(x, y, pch = 19, col = rgb(0, 0, 0, 0.2), cex = 0.5)
          lines(lowess(x, y), col = "#ca0020", lwd = 2)
        },
        upper.panel = function(x, y, ...) {
          usr <- par("usr")
          on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          # Calculate correlation
          r <- cor(x, y, use = "complete.obs", method = cor_method)
          # Color based on correlation strength
          col <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(200)[round((r + 1) * 99.5) + 1]
          rect(0, 0, 1, 1, col = col, border = NA)
          # Add correlation value
          txt <- format(r, digits = 2)
          text(0.5, 0.5, txt, cex = 1.5, font = 2)
        },
        diag.panel = function(x, ...) {
          usr <- par("usr")
          on.exit(par(usr))
          par(usr = c(usr[1:2], 0, 1.5))
          h <- hist(x, plot = FALSE)
          breaks <- h$breaks
          nB <- length(breaks)
          y <- h$counts
          y <- y/max(y)
          rect(breaks[-nB], 0, breaks[-1], y, col = "lightgray", border = "black")
        },
        gap = 0.2,
        cex.labels = 1.6,# variable label
        main = paste0("Correlation Scatterplot Matrix (", cor_method, ")")
  )

  dev.off()
  cat("✓ Saved: corr_scatterplot_matrix.png (full scatterplot matrix)\n")
}

# =============================================================================
# CALCULATE AND PRINT VIF TABLE FOR ALL VARIABLES (using usdm)
# =============================================================================

cat("=== Calculating VIF for ALL Variables ===\n")

# Assign short names
colnames(env_sample) <- short_names

# Compute VIF for all vars
vif_all <- vif(env_sample)

# Results as data.frame
vif_table_all <- data.frame(
  Variable = vif_all$Variables,
  VIF = round(as.numeric(vif_all$VIF), 3)
)
rownames(vif_table_all) <- NULL

print(vif_table_all, row.names = FALSE)

cat("\nVIF > 5-10 indicates multicollinearity issues for MaxEnt modeling.\n")

# Save the VIF data to CSV
write.csv(vif_table_all, file.path(output_data, "Q1-a-VIF.csv"), row.names = FALSE)
cat("Data saved")
