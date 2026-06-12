# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q4: Region comparisons
# Clustered Box-and-Whisker style plots of Habitat Suitability by Category
# created by Nicole Drewitz on 20 feb, 2026 with claude.ai
# =============================================================================

# =============================================================================
# Guide
# For each change in field changes need to be made (commented/uncommented) for:
#     a. Paths (input_file, output_file, output_CSV)
#     b. category colours
#     c. comment/uncomment plot combinations (e.g.coast/inland used violinplot, boxplot, summmary_stats mean dot)
#     d. plot subtitle if using boxplots
#     e. plot 3 starts with a line excluding categories with the prefix "Inland"
# =============================================================================

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis) # for generating colours automatically
library(ggridges) # only for ridge plot

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
data_dir   <- here::here("Data/Q4-random_sampling_locations-with_values")
fig_dir    <- here::here("figures")
output_data <- here::here("analysis_output_data")

input_file  <- file.path(data_dir,
                         #"Q4-0b-1-random_points_suitability-Inland-n1000each-mindist10000m.csv" # Greenland and Svalbard excluded
                         "Q4-0b-3-random_points_suitability-ECO_Cor-n300each-mindist5000m.csv"
                         #"Q4-0b-2-random_points_suitability-CoastCNT-n200each-mindist2000m.csv"
                         )
output_file <- file.path(fig_dir,
                         #"Q4-b-boxplots_suitability_by_scenario-coast_inland.png"
                         "Q4-b-boxplots_suitability_by_ecoregion.png"
                         #"Q4-b-ridgeplots_suitability_by_CoastCNT.png"
                         )
output_CSV <- file.path(output_data, #               # assign end of file name
                        #"coast_inland.csv"
                        "ecoregions.csv"
                        #"CoastCNT"
                        )

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------
df <- readr::read_csv(input_file, show_col_types = FALSE)

# The first column holds the category (e.g. "coast" / "inland").
names(df)[1] <- "category"

# -----------------------------------------------------------------------------
# 3. Reshape to long format
# -----------------------------------------------------------------------------
df_long <- df %>%
  dplyr::select(category, dplyr::starts_with("suitability_")) %>%
  tidyr::pivot_longer(
    cols      = dplyr::starts_with("suitability_"),
    names_to  = "scenario",
    values_to = "suitability"
  ) %>%
  dplyr::mutate(
    # Clean up scenario labels for the plot
    scenario = dplyr::recode(scenario,
                             "suitability_current" = "Current",
                             "suitability_ssp126"  = "SSP1-2.6",
                             "suitability_ssp370"  = "SSP3-7.0",
                             "suitability_ssp585"  = "SSP5-8.5"
    ),
    # Fix factor order so scenarios run chronologically
    scenario = factor(scenario,
                      levels = c("Current", "SSP1-2.6", "SSP3-7.0", "SSP5-8.5"))
  )

# -----------------------------------------------------------------------------
# 4. Set category colours - 3 options
# -----------------------------------------------------------------------------
# generating category colours
#n_cats <- length(unique(df_long$category))  # Counts categories dynamically
#category_colours <- viridis(n_cats, option = "D")  # "D" = viridis; try "A"(magma),"B"(inferno),"C"(plasma)
#names(category_colours) <- unique(df_long$category)

# or use manual colour setting option
category_colours <- c(
  #"coast"  = "#3d85c8", "inland" = "#e69500")

  # ecoregions below (Eco_Cor field)
  "Baltic mixed forests" = "#225555", "European Atlantic mixed forests" = "#DDCC77", "Kalaallit Nunaat Arctic steppe" = "#88CCEE", "Kola Peninsula tundra" = "#5237dc", #"Rock and Ice" = "#BBBBBB", # more colours "#AA4499", "#CC6677"
  "Russian Arctic desert"= "#5D6D7E", "Sarmatic mixed forests" = "#F0E442", "Scandinavian and Russian taiga" = "#E69F00", "Scandinavian coastal conifer forests" = "#0072B2", "Scandinavian Montane Birch forest and grasslands" = "#ab2b6b")

# -----------------------------------------------------------------------------
# 5. Option a - Bar plots
# -----------------------------------------------------------------------------
# Plot 1

# adjust plot width spacing between plot in same climate scenario
width_total <- 0.8  # total width to spread across per scenario

# Create a numeric x + offset column
df_plot <- df_long %>%
  mutate(
    x_num  = as.numeric(factor(scenario)),
    cat_idx = as.numeric(factor(category)),  # Simple: just get numeric index of categories
    x_nudge = x_num + (cat_idx / (length(unique(category)) + 1) - 0.5) * width_total
  )

p <- ggplot(df_plot, aes(x = x_nudge, y = suitability, fill   = category, colour = category, # uncomment this line to include categories starting with "Inland"
                group  = interaction(scenario, category))) +

  #geom_violin(width = 0.5, alpha = 0.75, linewidth = 0.55, trim = FALSE) +
  #geom_boxplot(aes(group = interaction(scenario, category)), width = 0.07, outlier.shape = 16, outlier.size  = 1.5, alpha = 0.9, colour = "black") +

  #stat_summary(fun = mean, geom = "point", size = 2.5, colour = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.5, linewidth = 0.8) + # 95% confidence interval

  scale_x_continuous(breaks = unique(df_plot$x_num), labels = levels(factor(df_plot$scenario))
  ) +
  scale_fill_manual(values = category_colours, name = "Category") +
  scale_colour_manual(values = category_colours, name = "Category") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    title    = "Habitat Suitability by Climate Scenario",
    subtitle = "change in mean with 95% confidence interval",
    #subtitle = "dot= mean, points randomly sampled. points from Greenland removed",
    x = "Climate scenario", y = "Habitat suitability"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, colour = "grey40", margin = margin(b = 10)),
    axis.title         = element_text(size = 12),
    axis.text          = element_text(size = 11),
    legend.title       = element_text(size = 11, face = "bold"),
    legend.text        = element_text(size = 10),
    legend.position    = "right",
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.4),
    plot.margin        = margin(12, 16, 12, 12)
  )
p
# -----------------------------------------------------------------------------
# Calculate means per scenario × category
# -----------------------------------------------------------------------------
means <- df_plot %>%
  group_by(scenario, category) %>%
  summarise(mean_y  = mean(suitability, na.rm = TRUE),x_nudge = mean(x_nudge), .groups = "drop"
  )

# -----------------------------------------------------------------------------
# Add line to plot
# -----------------------------------------------------------------------------
p2 <- p +
  geom_line(
    data     = means, #medians,
    aes(x    = x_nudge, y    = mean_y, #median_y,
        colour = category, group  = category), linewidth = 0.8, linetype  = "solid"
  ) +
  geom_point(
    data   = means, #medians,
    aes(x  = x_nudge,
        y  = mean_y, #median_y,
        colour = category), size   = 3, shape  = 18   # diamond shape to distinguish from mean dots
  )
p2
# -----------------------------------------------------------------------------
# Complete summary table (means + changes + sample sizes)
# -----------------------------------------------------------------------------
df_means <- df_long %>%
  group_by(category, scenario) %>%
  summarise(mean_suit = mean(suitability, na.rm = TRUE), .groups = "drop")

# Join means back onto df_long so each row knows its category-scenario mean
df_plot <- df_long %>%
  left_join(df_means, by = c("category", "scenario"))

mean_changes <- df_long %>%
  filter(scenario == "Current") %>%
  count(category, name = "n_points") %>%
  left_join(
    df_means %>%
      tidyr::pivot_wider(names_from = scenario, values_from = mean_suit, names_prefix = "mean_") %>%
      mutate(
        change_126  = get("mean_SSP1-2.6") - mean_Current,
        change_370  = `mean_SSP3-7.0` - mean_Current,
        change_585  = `mean_SSP5-8.5` - mean_Current
      ),
    by = "category"
  ) %>%
  dplyr::select(category, n_points, starts_with("mean_"), starts_with("change_")) %>%
  print()

# -----------------------------------------------------------------------------
# SAVE CSV - already made with Q4-0b script
# -----------------------------------------------------------------------------
#base_name <- tools::file_path_sans_ext(basename(output_CSV))  # Gets name set at start of script
#out_file <- file.path(output_data, paste0("Q4-b-suitability_summary-", base_name, ".csv"))
#write_csv(mean_changes, out_file)
#message("CSV saved to: ", out_file)

# -----------------------------------------------------------------------------
# 5. Option b - Plot 3 - ridge plot
# -----------------------------------------------------------------------------

# Colour ramp extracted from QGIS QML file
# Below-threshold (< 0.5092): grey; above-threshold: blue → white → red

qgis_colours <- c(         # colours match suitability maps
  "0"       = "black",
  "0.50919" = "#818181",
  "0.5092"  = "#023858",   # deep navy blue
  "0.55828" = "#0571b0",   # strong blue
  "0.60736" = "#2c89c4",   # mid blue
  "0.65644" = "#74b4d4",   # light blue
  "0.70552" = "#b8d5e8",   # pale blue
  "0.7546"  = "#e8c9bb",   # pale red-warm (replaces near-white)
  "0.80368" = "#e8896a",   # salmon
  "0.85276" = "#d94f3a",   # strong red-orange
  "0.90184" = "#b81c26",   # deep red
  "0.95092" = "#8b0000",   # dark red
  "1"       = "#580000"    # very dark red
)

# Build a colorRamp function from the QGIS stops
qgis_stops     <- as.numeric(names(qgis_colours))
qgis_hex       <- unname(qgis_colours)

# Interpolation function: given a value 0-1, return a hex colour
qgis_ramp <- colorRamp(qgis_hex)   # equally spaced by default, so we map via stops
suit_to_colour <- function(x) {
  # map x through the unevenly-spaced stops
  clipped <- pmin(pmax(x, 0), 1)
  # find position in [0,1] relative to the stop positions
  pos <- approx(qgis_stops, seq(0, 1, length.out = length(qgis_stops)), xout = clipped)$y
  rgb_vals <- qgis_ramp(pos)
  rgb(rgb_vals[, 1], rgb_vals[, 2], rgb_vals[, 3], maxColorValue = 255)
}

# -----------------------------------------------------------------------------
# Plot

p3 <- ggplot(df_plot %>% filter(!grepl("^Inland", category)), # use this line to exclude categories starting with "Inland" (and comment out next line)
#p3 <- ggplot(df_plot,
                aes(x = suitability, y  = category,
                fill = mean_suit,   # fill = mean suitability of that ridge
                group = interaction(category, scenario))) +
  geom_density_ridges(
    alpha          = 0.85,
    scale          = 0.9,
    rel_min_height = 0.01,
    colour         = "grey30",
    linewidth      = 0.3
  ) +
  facet_wrap(~ scenario, nrow = 1) +
  # Colour gradient matching QGIS ramp
  scale_fill_gradientn(
    colours = qgis_hex,
    values  = scales::rescale(qgis_stops),   # map unevenly-spaced stops correctly
    limits  = c(0, 1),
    name    = "Mean\nsuitability"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Habitat Suitability by Climate Scenario",
    x     = "Habitat suitability",
    y     = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    axis.title         = element_text(size = 12),
    axis.text.x        = element_text(size = 10),
    axis.text.y        = element_text(size = 12),
    legend.title       = element_text(size = 13, face = "bold"),
    legend.text        = element_text(size = 11),
    legend.position    = "right",
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4),
    strip.text         = element_text(size = 12, face = "bold")
  )
p3

# -----------------------------------------------------------------------------
# 6. Save
# -----------------------------------------------------------------------------
#ggsave(filename = output_file, plot = p, width    = 9, height   = 5.5, dpi      = 300)
ggsave(filename = output_file, plot  = p2, width    = 9, height   = 5.5, dpi      = 300)
#ggsave(filename = output_file, plot  = p3, width    = 9, height   = 5.5, dpi      = 300)
#ggsave(filename = output_file, plot  = p3, width    = 14, height   = 12, dpi      = 300) # large plot
message("Plot saved to: ", output_file)
