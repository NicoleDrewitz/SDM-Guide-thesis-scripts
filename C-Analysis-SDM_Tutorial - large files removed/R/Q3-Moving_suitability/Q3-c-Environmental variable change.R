# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q3
# R Script to Plot Climate Scenario Comparisons
# Author: Nicole Drewitz with Claude.ai
# Date: 2026-02-16
# Last updated: 19 Feb, 2026 with Perplexity AI
# Description: Creates box-and-whisker plot matrices comparing environmental variables
#              across climate scenarios
# =============================================================================

# Load required libraries
required_packages <- c("here", "dplyr", "tidyr", "ggplot2", "patchwork")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Define paths relative to the RProject root
data_dir <- here::here("data")
fig_dir  <- here::here("figures")

# Read the extracted data
cat("Reading extracted environmental data...\n")
data <- read.csv(file.path(data_dir, "Q3-extracted_environmental_values-analysis_locations.csv"))

# =============================================================================
# DEFINE VARIABLES
# =============================================================================
# Define current-only variables (static across climate scenarios)
# After scaling data, RE-run detection with correct names
current_only_variables <- c(
  "soil_phh2o_0.5.30cm_mean_averaged",
  "soil_soc_0.5.30cm_mean_averaged",
  "soil_waterV33_LoamClay_0033kPa_0.5.30cm_mean_averaged",
  "Topo_northness"
)

categorical_variables <- "kg0"

# Identify which variables are climate-dependent vs current-only
all_variables <- setdiff(names(data), c("sample.No", "longitude", "latitude", "climateScenario"))

# Check which current_only patterns match actual column names
# Note: Using gsub to normalize hyphens/periods for matching
climate_dependent_vars <- all_variables[!sapply(all_variables, function(var) {
  var_normalized <- gsub("-", ".", var, fixed = TRUE)
  any(sapply(current_only_variables, function(pattern) {
    pattern_normalized <- gsub("-", ".", pattern, fixed = TRUE)
    grepl(pattern_normalized, var_normalized, fixed = TRUE)
  }))
})]
# exclude categorical variable(s)
climate_dependent_vars <- setdiff(climate_dependent_vars, categorical_variables)
climate_dependent_vars

current_only_vars <- all_variables[sapply(all_variables, function(var) {
  var_normalized <- gsub("-", ".", var, fixed = TRUE)
  any(sapply(current_only_variables, function(pattern) {
    pattern_normalized <- gsub("-", ".", pattern, fixed = TRUE)
    grepl(pattern_normalized, var_normalized, fixed = TRUE)
  }))
})]
current_only_vars

cat(paste0("\nFound ", length(climate_dependent_vars), " climate-dependent variables\n"))
cat(paste0("Found ", length(current_only_vars), " current-only variables\n"))

# =============================================================================
# VARIABLE LABELS AND UNITS
# =============================================================================

# Title labels - FIXED to match actual column names with DOTS
var_titles <- c(
  # Climate-dependent (these work)
  "bio11" = "Temperature of Coldest Quarter",
  "bio13" = "Precipitation in wettest month",
  "bio19" = "Precipitation of Coldest Quarter",
  "fcf"   = "Frost change frequency",
  "swe"   = "Snowpack",

  # Current-only - FIXED (dots instead of dashes)
  "soil_phh2o_0.5.30cm_mean_averaged" = "Soil pH in top 30 cm",
  "soil_soc_0.5.30cm_mean_averaged" = "Soil Organic Carbon",
  "soil_waterV33_LoamClay_0033kPa_0.5.30cm_mean_averaged" = "Soil water content",
  "Topo_northness" = "Aspect",

  "kg0" = "Köppen-Geiger climate class"
)

# Y-axis unit labels - FIXED (dots instead of dashes)
var_units <- c(
  # Climate-dependent (these work)
  "bio11" = "Daily mean (°C)",
  "bio13" = "Mean (kg/m²)",
  "bio19" = "Monthly mean (kg/m²)",
  "fcf"   = "Annual frequency",
  "swe"   = "kg/m²/year",

  # Current-only - FIXED (dots instead of dashes)
  "soil_phh2o_0.5.30cm_mean_averaged" = "mean pH (in soil water)",
  "soil_soc_0.5.30cm_mean_averaged" = "mean (dg/kg)",
  "soil_waterV33_LoamClay_0033kPa_0.5.30cm_mean_averaged" = "% (at 33 kPa)",
  "Topo_northness" = "(northness)",

  "kg0" = "class"
)

# Helper functions - FIXED (no magrittr dependency)
get_var_title <- function(var) {
  if (var %in% names(var_titles)) {
    return(var_titles[[var]])
  } else {
    return(var)  # fallback to raw name
  }
}

get_var_unit <- function(var) {
  if (var %in% names(var_units)) {
    return(var_units[[var]])
  } else {
    return("")  # no units
  }
}

# Read the extracted data
cat("Reading extracted environmental data...\n")
data <- read.csv(file.path(data_dir, "Q3-extracted_environmental_values-analysis_locations.csv"))

# [Rest of your variable detection code stays the same...]
current_only_variables <- c(
  "soil_phh2o_0-5-30cm_mean_averaged",
  "soil_soc_0-5-30cm_mean_averaged",
  "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged",
  "Topo_northness"
)

scenario_colours <- c(
  "current" = "#feebe2",
  "ssp126" = "#b3cde3",
  "ssp370" = "#8c96c6",
  "ssp585" = "#88419d"
)

# =============================================================================
# PLOTTING FORMATTING
# =============================================================================

# Function to create a box plot for a single variable
create_boxplot <- function(data, variable_name, show_legend = FALSE) {
  plot_data <- data %>%
    filter(!is.na(.data[[variable_name]]))

  p <- ggplot(plot_data, aes(x = climateScenario, y = .data[[variable_name]], fill = climateScenario)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(
      values = scenario_colours,
      labels = c("1981-2010" = "Current", "ssp126" = "SSP1-2.6",
                 "ssp370" = "SSP3-7.0", "ssp585" = "SSP5-8.5")
    ) +
    labs(
      title = get_var_title(variable_name),           # Short title
      x = NULL,
      y = paste0(get_var_unit(variable_name)),
      fill = "Climate Scenario"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = if(show_legend) "bottom" else "none"
    )
  return(p)
}

# Function to create a box plot for current-only variables
create_current_only_boxplot <- function(data, variable_name) {
  plot_data <- data %>%
    filter(climateScenario == "current", !is.na(.data[[variable_name]]))

  p <- ggplot(plot_data, aes(x = "", y = .data[[variable_name]])) +
    geom_boxplot(fill = "#feebe2", outlier.size = 0.5) +
    # scale by dividing by 10
    scale_y_continuous(
      name = get_var_unit(variable_name),
      labels = if(variable_name %in% c("soil_phh2o_0.5.30cm_mean_averaged",
                                       "soil_waterV33_LoamClay_0033kPa_0.5.30cm_mean_averaged")) {
        function(x) paste0(round(x/10, 1))
      } else {
        waiver()
      }
    ) +

    labs(
      title = get_var_title(variable_name),           # Short title
      x = NULL,
      y = paste0(get_var_unit(variable_name))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  return(p)
}

# ============================================================================
# CLIMATE-DEPENDENT VARIABLES PLOT MATRIX
# ============================================================================

if (length(climate_dependent_vars) > 0) {
  cat("\nCreating plot matrix for climate-dependent variables...\n")

  # Create individual plots
  plot_list <- lapply(seq_along(climate_dependent_vars), function(i) {
    var <- climate_dependent_vars[i]
    # Show legend only on the last plot
    show_legend <- (i == length(climate_dependent_vars))
    create_boxplot(data, var, show_legend = show_legend)
  })

  # Calculate grid dimensions (aim for roughly square layout)
  n_plots <- length(plot_list)
  n_cols <- ceiling(sqrt(n_plots))
  n_rows <- ceiling(n_plots / n_cols)

  # Combine plots using patchwork
  combined_plot <- wrap_plots(plot_list, ncol = n_cols) +
    plot_annotation(
      title = "Climate-Dependent Environmental Variables Across Climate Scenarios",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )

  # Save the plot
  output_file_climate <- file.path(fig_dir, "Q3-c-Change_in_variables-continuous.png")

  # Calculate dynamic height based on number of rows
  plot_height <- max(8, n_rows * 3)
  plot_width <- max(10, n_cols * 3)

  ggsave(
    output_file_climate,
    combined_plot,
    width = plot_width,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )

  cat(paste0("Saved climate-dependent variables plot to: ", output_file_climate, "\n"))
  cat(paste0("  Dimensions: ", plot_width, " x ", plot_height, " inches\n"))
}

# ============================================================================
# CURRENT-ONLY VARIABLES PLOT MATRIX
# ============================================================================

if (length(current_only_vars) > 0) {
  cat("\nCreating plot matrix for current-only variables...\n")

  # Create individual plots
  plot_list_current <- lapply(current_only_vars, function(var) {
    create_current_only_boxplot(data, var)
  })

  # Calculate grid dimensions
  n_plots <- length(plot_list_current)
  n_cols <- ceiling(sqrt(n_plots))
  n_rows <- ceiling(n_plots / n_cols)

  # Combine plots using patchwork
  combined_plot_current <- wrap_plots(plot_list_current, ncol = n_cols) +
    plot_annotation(
      title = "Static Environmental Variables (Current Climate Only)",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )

  # Save the plot
  output_file_current <- file.path(fig_dir, "Q3-c-Soil_and_northness_variables.png")

  # Calculate dynamic height based on number of rows
  plot_height <- max(8, n_rows * 3)
  plot_width <- max(10, n_cols * 3)

  ggsave(
    output_file_current,
    combined_plot_current,
    width = plot_width,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )

  cat(paste0("Saved current-only variables plot to: ", output_file_current, "\n"))
  cat(paste0("  Dimensions: ", plot_width, " x ", plot_height, " inches\n"))
}

# ============================================================================
# CATEGORICAL VARIABLES STACKED BAR CHART
# ============================================================================

if (length(categorical_variables) > 0) {
  cat("\nCreating stacked bar charts for categorical variables...\n")

  # Create individual stacked bar charts for each categorical variable
  plot_list_categorical <- lapply(categorical_variables, function(var) {
    # Count occurrences of each category by climate scenario
    plot_data <- data %>%
      filter(!is.na(.data[[var]])) %>%
      group_by(climateScenario, .data[[var]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(climateScenario) %>%
      mutate(
        total = sum(count),
        percentage = (count / total) * 100
      )

    # Limit to top 7 categories if more exist
    top_categories <- plot_data %>%
      group_by(.data[[var]]) %>%
      summarise(total_count = sum(count), .groups = "drop") %>%
      arrange(desc(total_count)) %>%
      slice_head(n = 7) %>%
      pull(.data[[var]])

    plot_data_filtered <- plot_data %>%
      filter(.data[[var]] %in% top_categories)

    # Create the stacked bar chart
    p <- ggplot(plot_data_filtered,
                aes(x = climateScenario, y = percentage,
                    fill = factor(.data[[var]]))) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Thinner bars
      scale_fill_brewer(palette = "Set3", name = "Class") + # or name = var
      scale_x_discrete(
        labels = c(
          "current" = "Current",
          "ssp126" = "SSP1-2.6",
          "ssp370" = "SSP3-7.0",
          "ssp585" = "SSP5-8.5"
        )
      ) +
      labs(
        x = "Climate Scenario",
        y = "Area (%)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Bigger x labels
        axis.text.y = element_text(size = 14),  # Bigger y labels
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12)
      )
    return(p)
  })

  # Calculate grid dimensions
  n_plots <- length(plot_list_categorical)
  if (n_plots == 1) {
    n_cols <- 1
    n_rows <- 1
  } else {
    n_cols <- ceiling(sqrt(n_plots))
    n_rows <- ceiling(n_plots / n_cols)
  }

  # Combine plots using patchwork
  combined_plot_categorical <- wrap_plots(plot_list_categorical, ncol = n_cols) +
    plot_annotation(
      title = "Koppen-Geiger classification",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )

  # Save the plot
  output_file_categorical <- file.path(fig_dir, "Q3-c-change_in_variables-kg0.png")

  # Calculate dynamic dimensions
  plot_height <- max(6, n_rows * 4)
  plot_width <- max(8, n_cols * 5)

  ggsave(
    output_file_categorical,
    combined_plot_categorical,
    width = plot_width,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )

  cat(paste0("Saved categorical variables plot to: ", output_file_categorical, "\n"))
  cat(paste0("  Dimensions: ", plot_width, " x ", plot_height, " inches\n"))

  # Print category counts
  cat("\nCategory counts by climate scenario:\n")
  for (var in categorical_variables) {
    cat(paste0("\n", var, ":\n"))
    category_summary <- data %>%
      filter(!is.na(.data[[var]])) %>%
      group_by(climateScenario, .data[[var]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(climateScenario, desc(count))
    print(category_summary)
  }
}

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")

# Calculate summary stats for climate-dependent variables
if (length(climate_dependent_vars) > 0) {
  cat("\nClimate-dependent variables by scenario:\n")
  summary_stats <- data %>%
    select(climateScenario, all_of(climate_dependent_vars)) %>%
    group_by(climateScenario) %>%
    summarise(across(everything(), list(
      mean = ~mean(.x, na.rm = TRUE),
      median = ~median(.x, na.rm = TRUE),
      sd = ~sd(.x, na.rm = TRUE)
    ), .names = "{.col}_{.fn}"))

  print(summary_stats)
}

cat("\n=== PLOTTING COMPLETE ===\n")
cat(paste0("Figures saved to: ", fig_dir, "\n"))

# ============================================================================
# ONE-WAY ANOVA + FDR CORRECTION (climateScenario effect)
# ============================================================================

if (length(climate_dependent_vars) > 0) {
  cat("\n=== One-way ANOVA with FDR correction ===\n")

  # Run one-way ANOVA for each variable: variable ~ climateScenario
  anova_results <- lapply(climate_dependent_vars, function(var) {
    f <- as.formula(paste(var, "~ climateScenario"))
    fit <- aov(f, data = data)

    tab <- summary(fit)[[1]]
    p_val <- tab["climateScenario", "Pr(>F)"]

    data.frame(
      variable = var,
      p_value = p_val,
      stringsAsFactors = FALSE
    )
  })

  anova_results <- do.call(rbind, anova_results)

  # FDR correction (Benjamini–Hochberg)
  anova_results$p_fdr <- p.adjust(anova_results$p_value, method = "BH")

  # Optional: significance flags
  alpha <- 0.05
  anova_results$significant_raw  <- anova_results$p_value  < alpha
  anova_results$significant_fdr  <- anova_results$p_fdr    < alpha

  print(anova_results)
}

#=== One-way ANOVA with FDR correction ===
# variable      p_value        p_fdr significant_raw significant_fdr
#1    bio11 0.000000e+00 0.000000e+00            TRUE            TRUE
#2    bio13 4.296778e-16 7.161297e-16            TRUE            TRUE
#3    bio19 1.365176e-05 1.365176e-05            TRUE            TRUE
#4      fcf 4.244489e-12 5.305611e-12            TRUE            TRUE
#5      swe 1.241330e-85 3.103325e-85            TRUE            TRUE

