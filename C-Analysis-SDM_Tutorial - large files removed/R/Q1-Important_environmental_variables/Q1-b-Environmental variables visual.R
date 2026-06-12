# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q1

# MaxEnt Species Distribution Model - Response Curves
# (Habitat suitability vs each Environmental variable)
# create 23 jan 2026 by Nicole Drewitz with claude.ai
# Last updated: 18 May, 2026

# Response curves:
# Others at mean/average = The partial effect of the variable when only that variable is changed (i.e. determined by how much suitability changes when only one variable changes in the final model)
# Single variable = The unadjusted relationship with species presence when all other variables are ignored
# If these values are very different, check if the environmental variable is correlated with another within that range.
# e.g. if winter temp just above freezing causes suitability to drop (sole predictor model), is it logical that the other the variables would remain the same (multi-variable model)?
# More here: https://naturalis.github.io/mebioda/doc/week2/w2d5/maxent_ouput_pres.pdf
# ===============================================================================

# Install and load required packages
required_packages <- c("here", "ggplot2", "dplyr", "tidyr", "patchwork", "stringr", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# Paths
input_dir <- here::here("Input_files/MaxEnt_files")
fig_dir  <- here::here("figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

#===============================================================================
# 1. SETUP & CONFIGURATION
#===============================================================================

# Model settings
THRESHOLD      <- 0.5092   # MaxSSS threshold - horizontal dashed line (under threshold is outside of predicted habitat range)
EXCLUDE_ZERO   <- TRUE     # Exclude y == 0 rows (except Topo_northness) (Can also change this to FALSE and remove y=0 from CSV) y=0, typically this means the value was outside of sampled range.
CAT_VAR        <- "kg0"    # Categorical variable name
DATASET_TYPE   <- "both"   # "only" | "average" | "both"

# Variable display labels
VAR_TITLES <- c(
  "bio11"             = "Temperature of Coldest Quarter",
  "bio13"             = "Precipitation in Wettest Month",
  "bio19"             = "Precipitation of Coldest Quarter",
  "fcf"               = "Frost Change Frequency",
  "swe"               = "Snowpack",
  "soil_phh2o"        = "Soil pH (0-30 cm depth)",
  "soil_soc"          = "Soil Organic Carbon",
  "soil_waterV33"     = "Soil Water Content",
  "Topo_northness"    = "Aspect",
  "kg0"               = "Köppen-Geiger Climate Class"
)

VAR_UNITS <- c(
  "bio11"               = "Daily mean (°C)",
  "bio13"               = "Mean (kg/m²)",
  "bio19"               = "Monthly mean (kg/m²)",
  "fcf"                 = "Annual frequency",
  "swe"                 = "kg/m²/year",
  "soil_phh2o"          = "Mean pH (in soil water)",
  "soil_soc"            = "Mean (dg/kg)",
  "soil_waterV33"       = "% (at 33 kPa)",
  "Topo_northness"      = "Northness (index))",
  "kg0"                 = "Class"
)

get_title <- function(var) ifelse(var %in% names(VAR_TITLES), VAR_TITLES[[var]], var)
get_unit  <- function(var) ifelse(var %in% names(VAR_UNITS),  VAR_UNITS[[var]],  "")

# rescale x-axis for 2 variables
get_x_scale <- function(var) {
  if (var == "soil_waterV33") {
    scale_x_continuous(
      name   = get_unit(var),
      labels = function(x) x / 10,  # Divide tick values by 10 when printing
      breaks = scales::breaks_width(width = 20)   # fewer labels than 10
    )
  } else if (var == "soil_phh2o") {
    scale_x_continuous(
      name   = get_unit(var),
      labels = function(x) x / 10,
      breaks = scales::breaks_width(width = 10)
    )
  } else {
    scale_x_continuous(name = get_unit(var))  # Default for all other continuous vars
  }
}

#===============================================================================
# 2. LOAD DATA
#===============================================================================

data_path <- function(...) here::here("Input_files/MaxEnt_files", ...)

if (DATASET_TYPE == "only") {

  data_only    <- read.csv(data_path("Q1-b-plots_only_dat_files_from_MaxEnt.csv"), sep = ";")
  grid_filename <- "Q1-b-cloudberry_response_curves-only_1_variable_models-2x5.png"

} else if (DATASET_TYPE == "average") {

  data_average  <- read.csv(data_path("Q1-b-plots_with_others_average_dat_files_from_MaxEnt.csv"), sep = ";")
  grid_filename  <- "Q1-b-cloudberry_response_curves-other-variables-mean_models-2x5.png"

} else if (DATASET_TYPE == "both") {

  data_only     <- read.csv(data_path("Q1-b-plots_only_dat_files_from_MaxEnt.csv"), sep = ";")
  data_average  <- read.csv(data_path("Q1-b-plots_with_others_average_dat_files_from_MaxEnt.csv"), sep = ";")
  grid_filename  <- "Q1-b-cloudberry_response_curves-comparison-2x5.png"

} else {
  stop("DATASET_TYPE must be 'only', 'average', or 'both'.")
}

# Derive variable names from column headers (strip x_ / y_ prefix)
ref_data  <- if (exists("data_only")) data_only else data_average
var_names <- unique(gsub("^[xy]_", "", names(ref_data)))
var_names <- var_names[var_names != ""]   # drop any empty strings

cat("Variables found:", paste(var_names, collapse = ", "), "\n")
cat("Total panels:   ", length(var_names), "\n")
var_names

#===============================================================================
# 3. SHARED THEME & SCALE HELPERS
#===============================================================================

# Fixed y limits so the threshold line never appears to "move"
Y_LIMITS <- c(0, 1)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(face = "bold", size = 14),
    axis.text         = element_text(size = 14),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.margin       = margin(t = 5, r = 5, b = 25, l = 5),
    legend.position   = "bottom"
  )

fixed_y_scale <- scale_y_continuous(
  name   = "Habitat suitability",
  labels = number_format(accuracy = 0.1),
  limits = Y_LIMITS,       # <-- KEY FIX: threshold stays at same position in every panel
  expand = expansion(mult = c(0, 0.02))
)

#===============================================================================
# 4. PLOT FUNCTIONS
#===============================================================================

# --- 4a. Single-dataset plot (gradient style) ---------------------------------

make_single_plot <- function(data, var, threshold) {

  x_col <- paste0("x_", var)
  y_col <- paste0("y_", var)

  pd <- data.frame(x = data[[x_col]], y = data[[y_col]]) %>%
    arrange(x)

  if (EXCLUDE_ZERO && var != "Topo_northness") pd <- filter(pd, y > 0)

  if (var == CAT_VAR) {
    # Categorical bar chart
    ggplot(pd, aes(x = as.factor(x), y = y)) +
      geom_col(fill = "#2166ac", color = "black", linewidth = 0.3) +
      geom_hline(yintercept = threshold, linetype = "dashed", linewidth = 0.6) +
      scale_x_discrete(
        name   = get_unit(var),
        labels = function(x) str_wrap(as.character(x), 8)
      ) +
      fixed_y_scale +
      labs(title = get_title(var)) +
      base_theme +
      theme(legend.position = "none")

  } else {
    # Continuous line + ribbon
    ggplot(pd, aes(x = x, y = y)) +
      geom_ribbon(
        aes(ymin = threshold, ymax = pmax(y, threshold)),
        fill = "#2166ac", alpha = 0.2
      ) +
      geom_line(color = "#2166ac", linewidth = 1.2) +
      geom_hline(yintercept = threshold, linetype = "dashed", linewidth = 0.6) +
      get_x_scale(var) +
      fixed_y_scale +
      labs(title = get_title(var)) +
      base_theme +
      theme(legend.position = "none")
  }
}

# --- 4b. Comparison plot (single-var vs others-mean) -------------------------

COLORS <- c("only" = "#2166ac", "average" = "#b2182b")
LABELS <- c("only" = "Single variable model", "average" = "Only 1 variable changed in multi-variable model")
#LABELS <- c("only" = "Single variable", "average" = "Others at mean")

make_comparison_plot <- function(data_only, data_average, var, threshold) {

  x_col <- paste0("x_", var)
  y_col <- paste0("y_", var)

  pd <- data.frame(
    x       = data_only[[x_col]],
    only    = data_only[[y_col]],
    average = data_average[[y_col]]
  ) %>%
    arrange(x) %>%
    pivot_longer(c(only, average), names_to = "model", values_to = "y")

  if (EXCLUDE_ZERO && var != "Topo_northness") pd <- filter(pd, y > 0)

  if (var == CAT_VAR) {
    ggplot(pd, aes(x = as.factor(x), y = y, fill = model)) +
      geom_col(position = "dodge", color = "black", linewidth = 0.3) +
      geom_hline(yintercept = threshold, linetype = "dashed", linewidth = 0.6) +
      scale_fill_manual(values = COLORS, labels = LABELS, name = NULL) +
      scale_x_discrete(
        name   = get_unit(var),
        labels = function(x) str_wrap(as.character(x), 8)
      ) +
      fixed_y_scale +
      labs(title = get_title(var)) +
      base_theme

  } else {
    ggplot(pd, aes(x = x, y = y, color = model, fill = model, group = model)) +
      geom_ribbon(
        aes(ymin = threshold, ymax = pmax(y, threshold)),
        alpha = 0.2, color = NA          # color = NA avoids ribbon outline
      ) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = threshold, linetype = "dashed", linewidth = 0.6) +
      scale_color_manual(values = COLORS, labels = LABELS, name = NULL) +
      scale_fill_manual(values  = COLORS, labels = LABELS, name = NULL) +
      get_x_scale(var) +
      fixed_y_scale +
      labs(title = get_title(var)) +
      base_theme
  }
}

#===============================================================================
# 5. BUILD PANEL GRID
#===============================================================================

cat("Building 2x5 grid...\n")

plots <- if (DATASET_TYPE == "both") {
  lapply(var_names, make_comparison_plot,
         data_only = data_only, data_average = data_average,
         threshold = THRESHOLD)
} else if (DATASET_TYPE == "only") {
  lapply(var_names, make_single_plot, data = data_only,    threshold = THRESHOLD)
} else {
  lapply(var_names, make_single_plot, data = data_average, threshold = THRESHOLD)
}

# Strip redundant y-axis labels from panels that aren't in the left column
# (columns 1 and 6 in a 5-column, 2-row layout)
left_col_idx <- c(1, 6)
for (i in seq_along(plots)) {
  if (!i %in% left_col_idx) {
    plots[[i]] <- plots[[i]] +
      theme(
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
}

# Remove legends from all panels
plots_no_legend <- lapply(plots, function(p) p + theme(legend.position = "none"))

# Create a dummy plot just for the legend
legend_plot <- ggplot(
  data.frame(model = c("only", "average"), x = 1, y = 1),
  aes(x, y, color = model)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = COLORS, labels = LABELS, name = NULL) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 14, face = "bold")
  )

# Combine panels + legend row
panel_letters <- letters[1:length(plots_no_legend)]

plots_tagged <- Map(function(p, lab) {
  p +
    labs(tag = lab) +
    theme(
      plot.tag = element_text(face = "bold"),
      plot.tag.position = "bottomleft",
      plot.margin = margin(t = 5, r = 5, b = 18, l = 5)
    )
}, plots_no_legend, panel_letters)

combined <- wrap_plots(plots_tagged, ncol = 5, nrow = 2) /
  wrap_plots(legend_plot) +
  plot_layout(heights = c(20, 1))
#===============================================================================
# 6. SAVE OUTPUT
#===============================================================================

out_path <- file.path(fig_dir, grid_filename)

png(out_path, width = 20 * 600, height = 10 * 600, res = 700)
print(combined)
dev.off()

cat("Saved:", out_path, "\n")
