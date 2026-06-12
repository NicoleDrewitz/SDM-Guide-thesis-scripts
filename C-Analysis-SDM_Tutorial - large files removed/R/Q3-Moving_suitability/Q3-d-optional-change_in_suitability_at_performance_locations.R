# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q3: optional

# difference in suitability between climates scenarios at true presence/absence locations used to assess model performance
# ANOVA included
# Create with Claude.ai on 27-01-2026 by Nicole Drewitz
# last updated 2026-01-29
# =============================================================================

# Load required libraries
required_packages <- c("dplyr", "tidyr", "ggplot2", "here")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# paths relative to the RProject root
data_dir <- here::here("Data")
output_data <- here::here("analysis_output_data")
fig_dir  <- here::here("figures")
input_csv <- file.path(
  data_dir, "cloudberry_suitability_at_analysis_locations.csv"
)

data <- read.csv(input_csv)
cat("Available columns:\n")
print(names(data))
cat("\n")

save_plots <- FALSE

# =============================================================================
# STATISTICAL ANALYSIS - ANOVA
# =============================================================================

# Reshape data for ANOVA
plot_data <- data %>%
  dplyr::select(presence.absence, current, ssp126, ssp370, ssp585) %>%
  tidyr::pivot_longer(cols = c(current, ssp126, ssp370, ssp585),
               names_to = "scenario",
               values_to = "suitability") %>%
  dplyr::mutate(scenario = factor(scenario, levels = c("current", "ssp126", "ssp370", "ssp585")),
         presence.absence = factor(presence.absence, levels = c(0, 1),
                                   labels = c("Absence", "Presence")))

# Two-way ANOVA: suitability ~ scenario * presence.absence
anova_model <- aov(suitability ~ scenario * presence.absence, data = plot_data)
anova_results <- summary(anova_model)

cat("\n=============================================================================\n")
cat("TWO-WAY ANOVA RESULTS\n")
cat("=============================================================================\n\n")
print(anova_results)

# Extract p-values for quarto
p_scenario <- anova_results[[1]]["scenario", "Pr(>F)"]
p_presence <- anova_results[[1]]["presence.absence", "Pr(>F)"]
p_interaction <- anova_results[[1]]["scenario:presence.absence", "Pr(>F)"]

# After running the ANOVA, save results to an RData file
anova_results_list <- list(
  anova_summary = anova_results,  # Contains the formatted "< 2.2e-16"
  anova_model = anova_model
)

# Save to data directory
saveRDS(anova_results_list, file.path(output_data, "3b-Change_at-Performance_locations-anova_results.rds"))

# =============================================================================
# PLOTTING
# =============================================================================

# Reshape data for plotting (long format)
plot_data <- data %>%
  dplyr::select(presence.absence, current, ssp126, ssp370, ssp585) %>%
  tidyr::pivot_longer(cols = c(current, ssp126, ssp370, ssp585),
               names_to = "scenario",
               values_to = "suitability") %>%
  dplyr::mutate(scenario = factor(scenario, levels = c("current", "ssp126", "ssp370", "ssp585")),
         presence.absence = factor(presence.absence, levels = c(0, 1),
                                   labels = c("Absence", "Presence")))

# Calculate medians for dashed lines
medians <- plot_data %>%
  group_by(scenario, presence.absence) %>%
  summarise(median = median(suitability, na.rm = TRUE), .groups = 'drop')

# Create violin plots
violin_plot <- ggplot(plot_data, aes(x = scenario, y = suitability, fill = presence.absence)) +
  geom_violin(trim = FALSE, alpha = 0.8,
              position = position_identity(),
              scale = "width") +
  # Add colored dashed lines at medians
  geom_segment(data = medians,
               aes(x = as.numeric(scenario) - 0.4,
                   xend = as.numeric(scenario) + 0.4,
                   y = median, yend = median,
                   color = presence.absence),
               linetype = "dashed", linewidth = 0.5,
               inherit.aes = FALSE) +
  scale_fill_manual(values = c("Absence" = "#7986CB", "Presence" = "#E57273"),
                    labels = c("True absence", "True presence")) +
  scale_color_manual(values = c("Absence" = "blue", "Presence" = "red"), # for median lines
                     guide = "none") +  # Hide color legend
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1),
                     labels = seq(0, 100, 10)) +
  scale_x_discrete(labels = c("current" = "Current",
                              "ssp126" = "SSP126",
                              "ssp370" = "SSP370",
                              "ssp585" = "SSP585")) +
  labs(x = NULL, y = "Suitability Score (%)", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position = c(0.19, 0.13),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.5), color = "black", linewidth = 0.5),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = 10),
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(violin_plot)

# =============================================================================
# SAVE PLOTS
# =============================================================================

if (save_plots) {
  # ensure figures directory exists
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  violin_plot_file <- file.path(fig_dir, "change_in_suitability-violin_plot.png")
  ggsave(violin_plot_file, plot = violin_plot, width = 8, height = 6, dpi = 300)
  cat("Violin plot saved as: ", violin_plot_file, "\n")
}
