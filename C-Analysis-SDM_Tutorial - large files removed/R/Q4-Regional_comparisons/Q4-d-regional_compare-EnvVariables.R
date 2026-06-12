# =============================================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q4:
# Script to Create Summary Statistics Table for Environmental Variables
# From Random Sampling Location Extraction Output Files (WITH REGION COLUMN)
# Adds Δ (future - current) columns for each scenario
# Author: Nicole Drewitz
# Date: 2026-02-25 (with Perplexity AI)
# Last updated: 2026-02-26

# remember that kg0 is categorical! - separate CSV produced at end
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tictoc)

tic()

# ── Paths ──────────────────────────────────────────────────────────────────────
data_dir   <- here::here("Data/Q4-random_sampling_locations-with_values")
output_dir <- here::here("analysis_output_data")

# ── Variable metadata ──────────────────────────────────────────────────────────
# Human-readable names and units for each environmental variable
var_titles <- c(
  "bio11" = "Temperature of Coldest Quarter",
  "bio13" = "Precipitation in Wettest Month",
  "bio19" = "Precipitation of Coldest Quarter",
  "fcf"   = "Frost Change Frequency",
  "swe"   = "Snowpack",
  "soil_phh2o_0-5-30cm_mean_averaged"              = "Soil pH in Top 30 cm",
  "soil_soc_0-5-30cm_mean_averaged"                = "Soil Organic Carbon",
  "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged" = "Soil Water Content (33 kPa)",
  "Topo_northness" = "Aspect",
  "kg0"            = "Köppen-Geiger Climate Class"
)

var_units <- c(
  "bio11" = "Daily mean (°C)",
  "bio13" = "Mean (kg/m²)",
  "bio19" = "Monthly mean (kg/m²)",
  "fcf"   = "Annual count of days crossing 0℃",
  "swe"   = "kg/m²/year",
  "soil_phh2o_0-5-30cm_mean_averaged"              = "mean pH (in soil water)",
  "soil_soc_0-5-30cm_mean_averaged"                = "mean (dg/kg)",
  "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged" = "% (at 33 kPa)",
  "Topo_northness" = "northness",
  "kg0"            = "class"
)

# Lookup table joining titles and units to variable names (used in joins later)
lookup_table <- tibble(
  variable = names(var_titles),
  title    = unname(var_titles),
  unit     = unname(var_units)
)

# ── Load data ──────────────────────────────────────────────────────────────────
# Read all per-region extraction CSVs into one combined dataframe
output_files <- list.files(data_dir, pattern = "^Q4-0d-.*_EnvVariables\\.csv$", full.names = TRUE)
cat("Found", length(output_files), "extraction output files\n")
if (length(output_files) == 0) stop("No extraction output files found.")

all_data <- map_dfr(output_files, \(f) {
  read_csv(f, show_col_types = FALSE) |> mutate(file_source = basename(f))
})

cat("Rows:", nrow(all_data), "| Regions:", paste(unique(all_data$region), collapse = ", "), "\n")

# Only keep variables that are both defined in metadata and present in the data
#available_vars <- intersect(names(var_titles), names(all_data))
available_vars <- intersect(names(var_titles), names(all_data)) |> setdiff("kg0")
missing_vars   <- setdiff(names(var_titles), available_vars)
if (length(missing_vars) > 0) cat("⚠️ Missing:", paste(missing_vars, collapse = ", "), "\n")

# ── Summary statistics ─────────────────────────────────────────────────────────
# Pivot to long format and compute n, mean, min, max per region × scenario × variable
summary_stats <- all_data |>
  select(region, climateScenario, all_of(available_vars)) |>
  pivot_longer(all_of(available_vars), names_to = "variable", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(region, variable, climateScenario) |>
  summarise(
    n_points = n(),
    mean_val = round(mean(value), 3),
    min_val  = round(min(value),  3),
    max_val  = round(max(value),  3),
    .groups  = "drop"
  )

# ── Climate deltas (future − current) ─────────────────────────────────────────
# For each region × variable, subtract current-climate mean from each future scenario
delta_table <- summary_stats |>
  select(region, variable, climateScenario, mean_val, min_val, max_val) |>
  pivot_wider(names_from = climateScenario, values_from = c(mean_val, min_val, max_val)) |>
  mutate(
    across(c(mean_val_ssp126, mean_val_ssp370, mean_val_ssp585),
           \(x) round(x - mean_val_current, 3), .names = "delta_{.col}"),
    across(c(min_val_ssp126,  min_val_ssp370,  min_val_ssp585),
           \(x) round(x - min_val_current,  3), .names = "delta_{.col}"),
    across(c(max_val_ssp126,  max_val_ssp370,  max_val_ssp585),
           \(x) round(x - max_val_current,  3), .names = "delta_{.col}")
  ) |>
  select(region, variable, starts_with("delta"))

# ── Export one CSV per environmental variable ──────────────────────────────────
# Attach human-readable titles, then split and save
summary_stats |>
  left_join(lookup_table, by = "variable") |>
  rename(`Environmental variable` = title) |>
  group_by(`Environmental variable`) |>
  group_walk(\(df, key) {
    safe_name <- gsub("[^a-zA-Z0-9_-]", "_", key$`Environmental variable`)
    out_file  <- file.path(output_dir, paste0("4b-env_summary_", safe_name, ".csv"))
    #write_csv(df, out_file, na = "")
    cat("Saved:", basename(out_file), "\n")
  })

toc()

# ── Build publication-style summary table ─────────────────────────────────────
# Target structure (matching screenshot):
#   Environmental variable | Unit | Region | Current climate (mean, min-max) |
#   Low emissions (Δ) | Medium emissions (Δ) | High emissions (Δ)

final_table <- summary_stats |>
  left_join(lookup_table, by = "variable") |>

  # Rescale pH and soil water content (stored as ×10 in source data)
  mutate(across(
    c(mean_val, min_val, max_val),
    \(x) if_else(
      variable %in% c("soil_phh2o_0-5-30cm_mean_averaged",
                      "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged"),
      round(x / 10, 4),
      x
    )
  )) |>
  # Round to 0 decimals for all variables except for aspect (2 decimal places)
  mutate(across(
    c(mean_val, min_val, max_val),
    \(x) if_else(
      variable == "Topo_northness",
      round(x, 2),
      round(x, 0)
    )
  )) |>

# Compute deltas vs current within each region × variable
  group_by(region, variable) |>
  mutate(
    current_mean = mean_val[climateScenario == "current"],
    current_min  = min_val[climateScenario  == "current"],
    current_max  = max_val[climateScenario  == "current"],
    delta_mean   = round(mean_val - current_mean, 3),
    #delta_min    = round(min_val  - current_min,  3),
    #delta_max    = round(max_val  - current_max,  3),
    delta_fmt = if_else(
      climateScenario != "current",
      #paste0(delta_mean, " (", delta_min, "–", delta_max, ")"),
      paste0(delta_mean),
      NA_character_
    )
  ) |>
  ungroup() |>

  # Keep only the columns we need before pivoting
  select(title, unit, region, climateScenario, mean_val, min_val, max_val, delta_fmt) |>
  
  # Pivot scenarios to wide
  pivot_wider(
    names_from  = climateScenario,
    values_from = c(mean_val, min_val, max_val, delta_fmt)
  ) |>
  
  # Rename to match headers
  rename(
    `Environmental variable` = title,
    `Unit`                   = unit,
    `Region`                 = region,
    `Current climate (mean)` = mean_val_current,
    `Current climate (min)`  = min_val_current,
    `Current climate (max)`  = max_val_current,
    `Change with low emissions (mean)`    = delta_fmt_ssp126,
    `Change with medium emissions (mean)` = delta_fmt_ssp370,
    `Change with high emissions (mean)`   = delta_fmt_ssp585
  ) |>
  
  select(-delta_fmt_current, -mean_val_ssp126, -mean_val_ssp370, -mean_val_ssp585,
         -min_val_ssp126,  -min_val_ssp370,  -min_val_ssp585,
         -max_val_ssp126,  -max_val_ssp370,  -max_val_ssp585) |>
  
  # Sort to match variable order in var_titles
  arrange(match(`Environmental variable`, unname(var_titles)), Region)

# ── Export ─────────────────────────────────────────────────────────────────────
final_table |>
  left_join(distinct(all_data, file_source, region), by = c("Region" = "region")) |>
  group_by(file_source) |>
  group_walk(\(df, key) {
    safe_name <- regmatches(key$file_source, regexpr("(?<=sampling-)(.+?)(?=-n)", key$file_source, perl = TRUE))
    out_file  <- file.path(output_dir, paste0("Q4-d-Regional-", safe_name, "-environmental_variables_summary_stats.csv"))
    write_csv(df, out_file, na = "")
    cat("Saved:", basename(out_file), "\n")
  })

# ==============================================================================
# ── Köppen-Geiger categorical summary ─────────────────────────────────────────
# For kg0 (categorical), compute most/least frequent class per region × scenario
# ==============================================================================

kg0_freq <- all_data |>
  select(region, climateScenario, kg0) |>
  filter(!is.na(kg0)) |>
  count(region, climateScenario, kg0) |>
  group_by(region, climateScenario) |>
  summarise(
    most_frequent_class  = kg0[which.max(n)],
    least_frequent_class = kg0[which.min(n)],
    .groups = "drop"
  ) |>
  arrange(region, climateScenario) |>
  left_join(distinct(all_data, file_source, region), by = "region") |>
  group_by(file_source) |>
  group_walk(\(df, key) {
    safe_name <- regmatches(key$file_source, regexpr("(?<=sampling-)(.+?)(?=-n)", key$file_source, perl = TRUE))
    out_file  <- file.path(output_dir, paste0("Q4-d-Regional-", safe_name, "-kg0_koppen_geiger_frequency_summary_stats.csv"))
    write_csv(df, out_file, na = "")
    cat("Saved:", basename(out_file), "\n")
  })