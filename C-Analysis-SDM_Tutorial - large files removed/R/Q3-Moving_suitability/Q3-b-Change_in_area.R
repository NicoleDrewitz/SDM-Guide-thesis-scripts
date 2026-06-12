# =============================================================================
# Nordic Cloudberry (Rubus chamaemorus) Distribution Modeling
# Area above specific suitability thresholds across climate scenarios
#
# Modified version: Tracks area above thresholds of 0.5092, 0.75, and 0.9
# Create with Claude.ai on 29-01-2026 by Nicole Drewitz
# last update 16-02-2026
# =============================================================================

required_packages <- c("raster", "dplyr", "tidyr", "ggplot2", "here")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

fig_dir  <- here::here("figures")
output_data <- here::here("analysis_output_data")

current <- raster(here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_1 - 1981-2010.asc"))
ssp126 <- raster(here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_2 - ssp126.asc"))
ssp370 <- raster(here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_3 - ssp370.asc"))
ssp585 <- raster(here::here("Input_files/MaxEnt_files/Rubus_chamaemorus-Scenario_4 - ssp585.asc"))

# Define suitability thresholds
thresholds <- c(0.5092, 0.75, 0.9)

# =============================================================================
# Function to calculate area above each threshold
# =============================================================================
calculate_area_above_thresholds <- function(raster_obj, scenario_name) {
  vals <- values(raster_obj)
  vals <- vals[!is.na(vals)]

  # Calculate area above each threshold
  results <- data.frame(
    scenario = scenario_name,
    threshold = thresholds,
    threshold_label = paste0(">", thresholds * 100, "%"),
    area_km2 = sapply(thresholds, function(t) sum(vals >= t))
  )

  return(results)
}

# Calculate for all scenarios
current_data <- calculate_area_above_thresholds(current, "Current")
ssp126_data <- calculate_area_above_thresholds(ssp126, "SSP1-2.6")
ssp370_data <- calculate_area_above_thresholds(ssp370, "SSP3-7.0")
ssp585_data <- calculate_area_above_thresholds(ssp585, "SSP5-8.5")

all_data <- bind_rows(current_data, ssp126_data, ssp370_data, ssp585_data)

# Assign numeric values to scenarios for plotting
all_data <- all_data %>%
  mutate(
    scenario_numeric = case_when(
      scenario == "Current" ~ 1,
      scenario == "SSP1-2.6" ~ 2,
      scenario == "SSP3-7.0" ~ 3,
      scenario == "SSP5-8.5" ~ 4
    ),
    scenario = factor(scenario, levels = c("Current", "SSP1-2.6", "SSP3-7.0", "SSP5-8.5"))
  )

# Calculate percent change from current for each threshold
change_data <- all_data %>%
  group_by(threshold) %>%
  mutate(
    current_area = area_km2[scenario == "Current"],
    percent_change = ((area_km2 - current_area) / current_area) * 100,
    absolute_change = area_km2 - current_area
  ) %>%
  ungroup()

# =============================================================================
# Save cell counts to CSV before plotting
# =============================================================================
write.csv(all_data, file.path(output_data, "Q3-b-area_above_thresholds_Nordic.csv"), row.names = FALSE)

# =============================================================================
# PLOTTING
# =============================================================================

# Create line plot showing area above each threshold
p1 <- ggplot(all_data, aes(x = scenario_numeric, y = area_km2,
                           color = threshold_label, group = threshold_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(">50.92%" = "#0571b0",
               ">75%" = "#643076",
               ">90%" = "#ca0020"),
    name = "Suitability\nThreshold"
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("Current", "SSP1-2.6", "SSP3-7.0", "SSP5-8.5")
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  labs(
    x = "Climate Scenario",
    y = "Habitat Area (km²)",
    title = "Area Above Suitability Thresholds for Rubus chamaemorus"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11)
  )

# Create bar plot showing percent change from current
p2 <- ggplot(change_data %>% filter(scenario != "Current"),
             aes(x = scenario, y = percent_change, fill = threshold_label)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(
    values = c(">=50.92%" = "#0571b0",
               ">=75%" = "#643076",
               ">=90%" = "#ca0020"),
    name = "Suitability\nThreshold"
  ) +
  labs(
    x = "Climate Scenario",
    y = "Change from Current (%)",
    title = "Percent Change in Area Above Suitability Thresholds"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11)
  )

# Save the plots
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# =============================================================================
# SAVE
# =============================================================================

output_file1 <- file.path(fig_dir, "Q3-b-area_above_thresholds_Nordic.png")
ggsave(output_file1, p1, width = 10, height = 6, dpi = 300)
cat(paste("\nLine plot saved to:", output_file1, "\n"))

output_file2 <- file.path(fig_dir, "Q3-b-percent_change_thresholds_Nordic.png")
ggsave(output_file2, p2, width = 10, height = 6, dpi = 300)
cat(paste("Bar plot saved to:", output_file2, "\n"))

# Save the data to CSV
if (!dir.exists(output_data)) {
  dir.create(output_data, recursive = TRUE)
}

write.csv(all_data, file.path(output_data, "Q3-b-area_above_thresholds_Nordic.csv"), row.names = FALSE)
write.csv(change_data, file.path(output_data, "Q3-b-threshold_changes_Nordic.csv"), row.names = FALSE)
cat("Data saved to output directory\n")
