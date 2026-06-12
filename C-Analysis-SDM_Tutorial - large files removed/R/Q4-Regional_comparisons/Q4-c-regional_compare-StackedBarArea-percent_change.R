# ============================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q4: Plotting area above thresholds summary with percent area from each category
# e.g. % habitat above threshold of total coastal area
# Stacked bars in % of total area
# created by Nicole Drewitz on 21 Feb, 2026 with Claude AI
# comment/uncomment out line 97-98 to exclude categories with prefix Inland
# ============================================================
library(dplyr)
library(ggplot2)
library(here)

# ── Config ────────────────────────────────────────────────────
Category_NAME  <-"Inland_Coastal"  # used in output file name & plot title. e.g. "Inland_Coastal", "Ecoregions", "Coasts_by_country"
HORIZONTAL  <- TRUE    # TRUE = horizontal bars; FALSE = vertical
WRAP_COLS   <- 1       # number of columns in the facet grid. Inland_Coastal= 1, Ecoregions= 2,Coasts_by_country= 4
# choose to comment out line at start of plot that excludes plotting categories with prefix "Inland"

# add file created with 0_4-part3-regional_comparisons-AreaSummary-coast_inland-percent_change.R
csv_path    <- here::here("analysis_output_data",
                          "Q4-0c-area_data-Inland_coastal-percent_change_total.csv" # use horizontal = FALSE
                          #"Q4-0c-area_data-Ecoregions-percent_change_total.csv"
                          #"Q4-0c-area_data-Coasts_by_country-percent_change_total.csv"
                          )

# ── Paths ─────────────────────────────────────────────────────
fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Read CSV ──────────────────────────────────────────────────
raw <- read.csv(csv_path, stringsAsFactors = FALSE)
raw$category <- trimws(raw$category)
raw$scenario <- trimws(raw$scenario)

message("CSV rows: ", nrow(raw))
message("Scenarios : ", paste(unique(raw$scenario),  collapse = ", "))
message("Categories: ", paste(unique(raw$category),  collapse = ", "))
message("Thresholds: ", paste(unique(raw$threshold), collapse = ", "))

# ── Factor levels (display order) ────────────────────────────
scenario_levels <- c("ssp585", "ssp370", "ssp126", "current")
scenario_labels <- c("SSP5-8.5", "SSP3-7.0", "SSP1-2.6", "Current")

# threshold column contains "poor" and numeric strings like "0.5092"
band_levels <- c("0.9", "0.75", "0.5092", "poor")
band_labels <- c("\u2265 0.90", "\u2265 0.75", "\u2265 0.51", "Poor habitat")

# ── Derive category levels & labels from CSV ──────────────────
csv_categories  <- sort(unique(raw$category))          # e.g. c("coastal","inland")
category_labels <- tools::toTitleCase(csv_categories)  # e.g. c("Coastal","Inland")

# ── Prepare plot data ─────────────────────────────────────────
plot_data <- raw |>
  mutate(
    scenario  = factor(scenario,
                       levels = scenario_levels,
                       labels = scenario_labels),
    threshold = factor(as.character(threshold),
                       levels = band_levels,
                       labels = band_labels),
    category  = factor(category,
                       levels = csv_categories,
                       labels = category_labels)
  ) |>
  filter(!is.na(scenario), !is.na(threshold), !is.na(category))

# Sanity check: pct per scenario × category should sum to ~100
check <- plot_data |>
  group_by(scenario, category) |>
  summarise(total_pct = sum(pct_area), .groups = "drop")
message("pct sums (should all be ~100):")
print(check)

# ── Colours ───────────────────────────────────────────────────
fill_colours <- setNames(
  c("#1A5C26", "#3A9E4A", "#A8D5A2", "#C4C4C4"),
  band_labels
)

# ── Plot ──────────────────────────────────────────────────────

# axis setup differs by orientation
if (HORIZONTAL) {
  scale_primary   <- scale_x_discrete(expand = expansion(add = 0.6))
  scale_secondary <- scale_y_continuous(labels = function(x) paste0(x, "%"),
                                        limits = c(0, 100.001),
                                        expand = expansion(mult = c(0.02, 0)))
  axis_labs <- labs(x = "% of Total Category Area", y = NULL)
} else {
  scale_primary   <- scale_x_discrete(expand = expansion(add = 0.6))
  scale_secondary <- scale_y_continuous(labels = function(x) paste0(x, "%"),
                                        limits = c(0, 100.001),
                                        expand = expansion(mult = c(0, 0.02)))
  axis_labs <- labs(x = NULL, y = "% of Total Category Area")
}

p <- ggplot(plot_data,
#p <- ggplot(plot_data |> filter(!grepl("^Inland", category, ignore.case = TRUE)), # use this line to exclude plotting categories starting with "Inland"
            aes(x = scenario, y = pct_area, fill = threshold)) +
  geom_col(position = position_stack(reverse = TRUE),
           alpha = 0.95, colour = NA) +
  { if (HORIZONTAL) coord_flip() } +
  scale_fill_manual(values = fill_colours, name = NULL) +
  scale_primary +
  scale_secondary +
  facet_wrap(~ category, ncol = WRAP_COLS, labeller = label_wrap_gen(width = 35)) + # width = characters allowed before line break
  axis_labs +
  #labs(title = paste0(tools::toTitleCase(Category_NAME), ": All Suitability \u2014 % of Category Area")) +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = 8, face = "bold", margin = margin(b = 1, t = 1)), # b = bottom margin space, t = top
    axis.text          = element_text(size = 9, margin = margin(r = 20)),
    axis.text.x        = element_text(
      angle = if (HORIZONTAL) 0 else 30,
      hjust = if (HORIZONTAL) 0.5 else 1),
    axis.title         = element_text(size = 12),
    axis.title.x       = element_text(margin = margin(t = 10)),
    axis.title.y       = element_text(margin = margin(r = 9)),
    plot.title         = element_text(size = 9, hjust = 0.5),
    legend.position    = "bottom",
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.5, "cm"),
    plot.margin        = margin(16, 20, 10, 16)
  )
p

# ── Save ──────────────────────────────────────────────────────
out_path <- file.path(fig_dir,
                      paste0("Q4-c-area-stackedBar-", Category_NAME,"-percent_change_of_area.png"))
#ggsave(out_path, plot = p, width = 9, height = 7.5, dpi = 300, bg = "#1C1C1E")

ggsave(out_path, plot = p, width = 9.72, height = 6.3, dpi = 500) # to print on full A4 paper landscape
message("Plot saved: ", out_path)
