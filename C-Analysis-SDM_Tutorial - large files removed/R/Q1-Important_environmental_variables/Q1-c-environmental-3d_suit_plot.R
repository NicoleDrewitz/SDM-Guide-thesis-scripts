# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q1

# 3D plot visualization relationship between 2 variables and habitat suitability
# Create with Claude.ai on 30-01-2026 by Nicole Drewitz
# =============================================================================

# Load required libraries
required_packages <- c("raster", "plotly", "here")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

fig_dir  <- here::here("figures")

# Set working directory and load data
setwd(here::here("Input_files/MaxEnt_files"))
current <- raster("Rubus_chamaemorus-Scenario_1 - 1981-2010.asc")

# Directory paths
PATHS <- list(
  env_raster_dir = here::here("Input_files/Environmental_variables")
)

# Variables to include
INCLUDE_VARIABLES <- c("bio11", "soil_soc_0-5-30cm_mean_averaged")

# Variable display names (optional, for plot labels)
VAR_NAMES <- c(
  "bio11" = "Winter temp (daily)",
  "soil_soc_0-5-30cm_mean_averaged" = "Soil Organic Carbon"
)

# Load environmental variables automatically
env_vars <- lapply(INCLUDE_VARIABLES, function(var) {
  raster(file.path(PATHS$env_raster_dir, paste0(var, ".asc")))
})
names(env_vars) <- INCLUDE_VARIABLES

# Function to create 3D scatter plot
plot_3d_suitability <- function(var1_name, var2_name, suitability,
                                var1_label = NULL,
                                var2_label = NULL,
                                max_points = 10000) {

  # Get rasters from loaded variables
  env_var1 <- env_vars[[var1_name]]
  env_var2 <- env_vars[[var2_name]]

  # Use custom labels or default to variable names
  if(is.null(var1_label)) var1_label <- ifelse(var1_name %in% names(VAR_NAMES), VAR_NAMES[var1_name], var1_name)
  if(is.null(var2_label)) var2_label <- ifelse(var2_name %in% names(VAR_NAMES), VAR_NAMES[var2_name], var2_name)

  # Stack rasters and convert to data frame
  stacked <- stack(env_var1, env_var2, suitability)
  df <- as.data.frame(stacked, xy = TRUE, na.rm = TRUE)
  colnames(df) <- c("x_coord", "y_coord", "env1", "env2", "suitability")

  # Remove any rows with NA values
  df <- na.omit(df)

  cat("Total points:", nrow(df), "\n")

  # Sample data to manageable size
  if(nrow(df) > max_points) {
    set.seed(123)
    df <- df[sample(nrow(df), max_points), ]
    cat("Sampled to:", max_points, "points\n")
  }

  # Create 3D scatter plot
  p <- plot_ly(df, x = ~env1, y = ~env2, z = ~suitability,
               type = "scatter3d", mode = "markers",
               marker = list(size = 2, color = ~suitability,
                             colorscale = "Viridis", showscale = TRUE,
                             colorbar = list(title = "Suitability"))) %>%
    layout(
      scene = list(
        xaxis = list(title = var1_label),
        yaxis = list(title = var2_label),
        zaxis = list(title = "Suitability", range = c(0, 1))
      ),
      title = ""
    )

  return(p)
}

# Create plot using variable names from INCLUDE_VARIABLES
plot_no_smooth <- plot_3d_suitability(
  var1_name = "bio11",
  var2_name = "soil_soc_0-5-30cm_mean_averaged",
  suitability = current,
  max_points = 10000
)

# Display plot
plot_no_smooth

