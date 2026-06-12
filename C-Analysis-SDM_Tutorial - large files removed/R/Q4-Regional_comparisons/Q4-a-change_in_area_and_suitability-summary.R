# ============================================================
# Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt
# Q4: build area summary tables and pie charts

# Merges "Q4-0c-area_data-*" and "Q4-0b-suitability_summary-*" CSVs
# from output_data into summary table CSVs matching
# the format: Region | Current suitability (mean) |
#   avg change ssp126/370/585 | Suitability sample size |
#   Current Area (km2) | Change in Area | Difference from Nordic average
# Created on 23 Feb, 2026 by Nicole Drewitz with Claude AI
# Last updated 24 Feb, 2026 with Perplexity AI
# ============================================================

library(here)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

data_dir    <- here::here("Data/Q4-random_sampling_locations-with_values")
fig_dir     <- here::here("figures")
output_data <- here::here("analysis_output_data")

# ----------------------------------------------------------
# USER OPTION: Select threshold for area data
# Available options: "0.5092", "0.75", "0.9"
# ----------------------------------------------------------
selected_threshold <- "0.75" # threshold for highly suitable habitat

# ----------------------------------------------------------
# USER OPTION: Manual mapping between area file tags and
#              suitability file tags.
#   area_tag     = suffix of Q4-0c-area_data-<tag>.csv
#   suit_tag     = suffix of Q4-0b-suitability_summary-<tag>.csv
#   exclude_prefix = character vector of Region name prefixes to drop
#                    from BOTH area and suitability data for this pair.
#                    Set to NA (no quotes) to keep all regions.
#   Output CSV will be named after the area_tag (left side).
# ----------------------------------------------------------
tag_pairs <- tribble(
  ~area_tag,                                     ~suit_tag,       ~exclude_prefix,
  "Coasts_by_country-percent_change_total",     "CoastCNT",       list("Inland"),
  "Ecoregions-percent_change_total",            "ECO_Cor",        list(NA),
  "Inland_coastal-percent_change_total",        "Inland",         list(NA)
  )

# ----------------------------------------------------------
# 1. Discover files and extract region tag from filename
#    Convention: Q4-0c-area_data-<tag>.csv / Q4-0b-suitability_summary-<tag>.csv
# ----------------------------------------------------------

area_files <- list.files(output_data, pattern = "^Q4-0c-area_data", full.names = TRUE)
suit_files <- list.files(output_data, pattern = "^Q4-0b-suitability_summary", full.names = TRUE)

if (length(area_files) == 0) stop("No area_data* files found.")
if (length(suit_files) == 0) stop("No suitability_summary* files found.")

extract_region_tag <- function(paths, type = c("area", "suitability")) {
  type <- match.arg(type)
  basename(paths) %>%
    str_remove("\\.csv$") %>%
    str_remove("^Q4-[^-]+-") %>%              # remove "Q4-0c-" or "Q4-0b-"
    str_remove("^area_data-") %>%             # remove "area_data-" for area files
    str_remove("^suitability_summary-") %>%   # remove "suitability_summary-" for suit files
    str_remove("-SUMMARY$")                   # remove trailing "-SUMMARY" for suit files
}

area_tags <- extract_region_tag(area_files, "area")
suit_tags  <- extract_region_tag(suit_files, "suitability")

message("Area file region tags found:        ", paste(area_tags, collapse = ", "))
message("Suitability file region tags found: ", paste(suit_tags,  collapse = ", "))

# Validate that every tag in the mapping actually exists on disk
missing_area <- setdiff(tag_pairs$area_tag, area_tags)
missing_suit <- setdiff(tag_pairs$suit_tag,  suit_tags)
if (length(missing_area) > 0) stop("tag_pairs area_tag(s) not found on disk: ", paste(missing_area, collapse = ", "))
if (length(missing_suit) > 0) stop("tag_pairs suit_tag(s) not found on disk: ",  paste(missing_suit,  collapse = ", "))

message("Processing ", nrow(tag_pairs), " matched pair(s).")

# ----------------------------------------------------------
# Helper: normalise category names to Title Case
# ----------------------------------------------------------
to_title <- function(x) str_to_title(str_trim(x))

# ----------------------------------------------------------
# 2. Process each matched pair
# ----------------------------------------------------------

process_pair <- function(area_tag, suit_tag, exclude_prefix) {

  # unwrap the list() wrapper used to allow NA alongside character vectors
  exclude_prefix <- unlist(exclude_prefix)
  do_exclude     <- !all(is.na(exclude_prefix))

  message("\n========== Processing: area='", area_tag, "'  suit='", suit_tag, "' ==========")
  if (do_exclude) {
    message("  Excluding regions with prefix(es): ", paste(exclude_prefix, collapse = ", "))
  }

  area_file <- area_files[area_tags == area_tag]
  suit_file <- suit_files[suit_tags  == suit_tag]

  area_raw <- read_csv(area_file, show_col_types = FALSE) %>%
    mutate(
      category  = to_title(category),
      threshold = as.character(as.numeric(threshold))
    )

  suit_raw <- read_csv(suit_file, show_col_types = FALSE)
  
  # Rename the first column to `category` (to match your existing logic)
  first_col <- names(suit_raw)[1]
  suit_raw <- suit_raw %>%
    rename(category = all_of(first_col)) %>%
    mutate(category = to_title(category))

  # ---- Exclude regions by prefix (applied to both files) ----
  if (do_exclude) {
    prefix_pattern <- paste0("^(", paste(exclude_prefix, collapse = "|"), ")")
    area_raw <- area_raw %>% filter(!str_detect(category, regex(prefix_pattern, ignore_case = TRUE)))
    suit_raw <- suit_raw %>% filter(!str_detect(category, regex(prefix_pattern, ignore_case = TRUE)))
    message("  Rows after exclusion — area: ", nrow(area_raw), "  suitability: ", nrow(suit_raw))
  }

  # ---- Filter area to selected threshold and above ----------
  available_thresholds <- unique(area_raw$threshold)
  message("  Available thresholds: ", paste(available_thresholds, collapse = ", "))

  if (!selected_threshold %in% available_thresholds) {
    stop(paste0(
      "selected_threshold '", selected_threshold, "' not found for area tag '", area_tag, "'.\n",
      "Available: ", paste(available_thresholds, collapse = ", ")
    ))
  }

  area_filtered <- area_raw %>%
    filter(as.numeric(threshold) >= as.numeric(selected_threshold))
  message("  Rows retained after threshold filter: ", nrow(area_filtered))
  
  # ---- Total current area (all thresholds, for denominator) ----
  total_area_all_thresholds <- area_raw %>%
    filter(str_detect(tolower(scenario), "current|baseline|historic")) %>%
    group_by(category) %>%
    summarise(total_area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")
  
  # ---- Current area -----------------------------------------
  current_area <- area_filtered %>%
    filter(str_detect(tolower(scenario), "current|baseline|historic")) %>%
    group_by(category) %>%
    summarise(current_area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")

  # ---- Change in area ---------------------------------------
  area_change <- area_filtered %>%
    filter(!str_detect(tolower(scenario), "current|baseline|historic")) %>%
    mutate(scenario_col = case_when(
      str_detect(tolower(scenario), "126|ssp1|1-2.6") ~ "area_ssp126",
      str_detect(tolower(scenario), "370|ssp3|3-7.0") ~ "area_ssp370",
      str_detect(tolower(scenario), "585|ssp5|5-8.5") ~ "area_ssp585",
      TRUE ~ paste0("area_", scenario)
    )) %>%
    group_by(category, scenario_col) %>%
    summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = scenario_col, values_from = area_km2) %>%
    left_join(current_area, by = "category") %>%
    mutate(across(starts_with("area_"),
                  ~ ifelse(current_area_km2 == 0, NA_real_,
                           (. - current_area_km2) / current_area_km2),
                  .names = "change_{.col}")) %>%
    # optionally drop the raw area_* columns if you don’t need them:
    select(-starts_with("area_"))

  # ---- Suitability summary ----------------------------------
  suit_summary <- suit_raw %>%
    rename(
      Region            = category,
      current_suit_mean = mean_Current,
      ssp126_suit       = `mean_SSP1-2.6`,
      ssp370_suit       = `mean_SSP3-7.0`,
      ssp585_suit       = `mean_SSP5-8.5`,
      avg_change_ssp126 = change_126,
      avg_change_ssp370 = change_370,
      avg_change_ssp585 = change_585,
      suitability_n     = n_points
    ) %>%
    select(Region, current_suit_mean, avg_change_ssp126, avg_change_ssp370,
           avg_change_ssp585, suitability_n)

  # ---- Join -------------------------------------------------
  combined <- suit_summary %>%
    left_join(area_change, by = c("Region" = "category")) %>%
    left_join(total_area_all_thresholds, by = c("Region" = "category")) %>%
    rename(
      `Current suitability (mean)` = current_suit_mean,
      `Average change ssp126`      = avg_change_ssp126,
      `Average change ssp370`      = avg_change_ssp370,
      `Average change ssp585`      = avg_change_ssp585,
      `Suitability sample size`    = suitability_n,
      `Current Area (km2)`         = current_area_km2,
      `Change in Area (ssp126)`    = change_area_ssp126,
      `Change in Area (ssp370)`    = change_area_ssp370,
      `Change in Area (ssp585)`    = change_area_ssp585
    ) %>%
    mutate(
      `Proportion of total area (current, highly suitable)` =
        `Current Area (km2)` / total_area_km2
    ) %>%
    select(-total_area_km2)

  # ---- Difference from Nordic average (within this file) ----
  nordic_avg_current    <- mean(suit_raw$mean_Current, na.rm = TRUE)
  nordic_avg_change_370 <- mean(suit_raw$change_370,   na.rm = TRUE)

  combined <- combined %>%
    mutate(
      `Difference from Nordic average suitability (current)`    = `Current suitability (mean)` - nordic_avg_current,
      `Difference from Nordic average suitability (change into ssp370)` = `Average change ssp370`      - nordic_avg_change_370
    )

  # ---- Round ------------------------------------------------
  combined <- combined %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  # ---- NA diagnostics ---------------------------------------
  cat("\n  -- NA COUNTS for area tag:", area_tag, "--\n")
  print(colSums(is.na(combined)))
  na_rows <- combined %>% filter(if_any(everything(), is.na))
  if (nrow(na_rows) > 0) {
    cat("  Rows with NAs:\n"); print(na_rows)
  } else {
    cat("  OK: No NAs.\n")
  }

  # ---- Write output (named after area_tag) ------------------
  out_path <- file.path(
    output_data,
    paste0("Q4-a-summary_table-", area_tag, ".csv")
  )
  write_csv(combined, out_path)
  message("  Written: ", out_path)
  combined  # return tibble for use after the loop
}

# ----------------------------------------------------------
# 3. Run over all pairs; results in a named list of tibbles
# ----------------------------------------------------------

combined_list <- pmap(tag_pairs, process_pair) %>%
  set_names(tag_pairs$area_tag)
message("\nDone. combined_list contains ", length(combined_list),
        " tibble(s): ", paste(names(combined_list), collapse = ", "))

# ============================================================
# PIE CHARTS: Current Area (km2) per region, per source CSV
# Colourblind-friendly palette (Okabe-Ito)
# ============================================================

library(ggplot2)
library(scales)

# Primary distinguishable colours (colourblind-safe, ~8 major categories)
main_palette <- c(
  "#E69F00", "#ab2b6b", "#0072B2", "#F0E442", "#56B4E9", "#CC79A7", "#009E73", "#D55E00" # palette for ecoregions (so coastal ecoregions are blue)
  # "#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7","#F0E442" # palette for coasts by country
  #"#E69F00","#56B4E9" palette for coastal vs inland
)

# Minor slices — greys of varying lightness
minor_palette <- c(
  "#AAAAAA","#888888","#666666","#444444"
)

# ----------------------------------------------------------
# Re-read area data and suitability files to rebuild area_filtered
# ----------------------------------------------------------

area_files <- list.files(output_data, pattern = "^Q4-0c-area_data", full.names = TRUE)
suit_files <- list.files(output_data, pattern = "^Q4-0b-suitability_summary", full.names = TRUE)

extract_region_tag <- function(paths, type = c("area", "suitability")) {
  type <- match.arg(type)
  basename(paths) %>%
    str_remove("\\.csv$") %>%
    str_remove("^Q4-[^-]+-") %>%
    str_remove("^area_data-") %>%
    str_remove("^suitability_summary-") %>%
    str_remove("-SUMMARY$")
}

area_tags <- extract_region_tag(area_files, "area")
suit_tags  <- extract_region_tag(suit_files, "suitability")

# Build a tibble with all area data and file source
area_all <- map2_dfr(area_files, area_tags, ~{
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(.source_file = basename(.y))
})

# Rename first column to `category` if needed
if (names(area_all)[1] != "category") {
  area_all <- area_all %>% rename(category = names(area_all)[1])
}

area_all <- area_all %>%
  mutate(category = to_title(category))

# Filter to selected threshold and current scenario
area_filtered <- area_all %>%
  mutate(threshold_num = as.numeric(str_extract(threshold, "\\d+\\.?\\d*"))) %>%
  filter(!is.na(threshold_num), threshold_num >= as.numeric(selected_threshold))

# Build suit_with_source
suit_with_source <- map2_dfr(suit_files, suit_tags, ~{
  df <- read_csv(.x, show_col_types = FALSE)
  first_col <- names(df)[1]
  df %>%
    rename(category = all_of(first_col)) %>%
    mutate(category = to_title(category),
           .source_file = basename(.x))
})

# ----------------------------------------------------------
# Now run the pie chart loop
# ----------------------------------------------------------

pie_plots <- list()

for (i in seq_len(nrow(tag_pairs))) {
  suit_tag       <- tag_pairs$suit_tag[i]
  exclude_prefix <- unlist(tag_pairs$exclude_prefix[[i]])
  do_exclude     <- !all(is.na(exclude_prefix))
  src_file       <- suit_files[suit_tags == suit_tag]
  
  regions_in_file <- suit_with_source %>%
    filter(.source_file == basename(src_file)) %>%
    pull(category) %>%
    unique()
  
  if (do_exclude) {
    prefix_pattern  <- paste0("^(", paste(exclude_prefix, collapse = "|"), ")")
    regions_in_file <- regions_in_file[
      !str_detect(regions_in_file, regex(prefix_pattern, ignore_case = TRUE))
    ]
  }
  
  # Subset area data for these regions at selected threshold, current scenario
  pie_data_raw <- area_filtered %>%
    filter(
      str_detect(tolower(scenario), "current|baseline|historic"),
      category %in% regions_in_file
    ) %>%
    group_by(category) %>%
    summarise(current_area_km2 = mean(area_km2, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(current_area_km2)) %>%
    mutate(pct = current_area_km2 / sum(current_area_km2))
  
  # Separate major (>=3%) and minor (<3%) slices
  major_data <- pie_data_raw %>% filter(pct >= 0.03)
  minor_data <- pie_data_raw %>% filter(pct <  0.03)
  has_minor  <- nrow(minor_data) > 0
  
  # Assign main palette colours to major slices
  major_colours <- setNames(
    main_palette[((seq_len(nrow(major_data)) - 1) %% length(main_palette)) + 1],
    major_data$category
  )
  
  # Assign individual grey shades to minor slices
  if (has_minor) {
    minor_colours_individual <- setNames(
      minor_palette[((seq_len(nrow(minor_data)) - 1) %% length(minor_palette)) + 1],
      minor_data$category
    )
  } else {
    minor_colours_individual <- c()
  }
  
  all_colours <- c(major_colours, minor_colours_individual)
  
  # For the legend: collapse all minor slice names to a single "< 3% area" entry
  if (has_minor) {
    legend_breaks <- c(names(major_colours), names(minor_colours_individual)[1])
    legend_labels <- c(names(major_colours), "< 3% area")
    override_fills <- c(unname(major_colours), "#888888")
  } else {
    legend_breaks  <- names(major_colours)
    legend_labels  <- names(major_colours)
    override_fills <- unname(major_colours)
  }
  
  pie_data <- pie_data_raw %>%
    mutate(
      pct_label = if_else(pct >= 0.04, paste0(round(pct * 100), "%"), ""),
      category  = factor(category, levels = category)
    )
  
  cat("\n", strrep("=", 60), "\n")
  cat("Source file:", basename(src_file), "\n")
  cat(strrep("=", 60), "\n")
  print(
    pie_data_raw %>%
      select(Region = category, `Current Area (km2)` = current_area_km2, `% of Total` = pct) %>%
      mutate(`% of Total` = round(`% of Total` * 100, 2))
  )
  
  p <- ggplot(pie_data, aes(x = "", y = current_area_km2, fill = category)) +
    geom_col(width = 1, colour = "white", linewidth = 0.5) +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(
      values = all_colours,
      breaks = legend_breaks,
      labels = legend_labels,
      name   = "Region",
      guide  = guide_legend(
        override.aes = list(fill = override_fills)
      )
    ) +
    geom_text(
      aes(label = pct_label),
      position = position_stack(vjust = 0.5),
      size     = 3.5,
      colour   = "black",
      fontface = "bold"
    ) +
    labs(
      title    = "Current Area (km²) by Region",
      subtitle = paste0(
        "Source: ", basename(src_file),
        "  |  Threshold: ", selected_threshold
      ),
      caption  = paste0("Total area: ", comma(round(sum(pie_data$current_area_km2))), " km²")
    ) +
    theme_void(base_size = 13) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle   = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      plot.caption    = element_text(hjust = 0.5, size = 9, colour = "grey40"),
      legend.position = "right",
      legend.title    = element_text(face = "bold"),
      plot.margin     = margin(10, 10, 10, 10)
    )
  
  print(p)
  
  file_stem <- tools::file_path_sans_ext(basename(src_file))
  out_fig <- file.path(fig_dir, paste0("Q4-a-pie_current_area-", file_stem, ".png"))
  ggsave(out_fig, plot = p, width = 8, height = 6, dpi = 200)
  pie_plots[[file_stem]] <- p
}

# See available plots
names(pie_plots)