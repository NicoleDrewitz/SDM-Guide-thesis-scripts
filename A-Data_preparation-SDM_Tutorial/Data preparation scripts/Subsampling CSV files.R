# =============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Subsampling CSV files by proportions
# Author: Nicole Drewitz with Claude AI
# Last updated 18-9-2025

# Usage: add your own CSV files to the "Input files" folder with columns 'species', 'longitude' and 'latitude'
# after running this script use GIS software to limit locations to 1 max/raster pixel (of environmental variables), this is a default setting build into MaxEnt.
# or write script to ensure each point is a less 1km apart (if environmental variables are 1km2 resolution)

# Data citation:
# GBIF.org User. (2025b). Occurrence Download—Present [Data set]. The Global Biodiversity Information Facility. https://doi.org/10.15468/DL.5CM8KM 
# Norwegian Environment Agency. (2024). ANO [Geodatabase]. Natur i Norge (NiN). https://kartkatalog.miljodirektoratet.no/Dataset/Details/2054?lang=en-us 
# =============================================================================
library(here) # sets path names within this R Project
library(readr)
# =============================================================================
# CONFIGURATION SECTION
# Edit the variables below to switch between presence and absence data runs,
# adjust split proportions, or change input/output paths.
# =============================================================================

# --- Absence data (commented out; uncomment to process absence records) ---
#split_percentages <- c(0.4, 0.4, 0.2)
#split_names <- c("select", "analysis", "remainder")
#input_csv <- here::here("Input files/A_1_d_Absence_all_2025-09-11.csv")
#output_directory <- here::here("Location subsets CSV/Locations-absence-cloudberry-DATE")
#file_base_name <- "cloudberry_absence"

# --- Presence data (active) --- (file provided as example)
split_percentages <- c(0.1, 0.1, 0.8) # 10% selecting model settings, 10% model performance analysis, 80% for modelling in MaxEnt
split_names <- c("select", "analysis", "MaxEnt")
input_csv <- here::here("Input files/A_1_d_Presence_all_2025-09-11.csv")
output_directory <- here::here("Location subsets CSV")
file_base_name <- "cloudberry_presence"

random_seed <- NULL # Set to an integer (e.g. 42) for a reproducible shuffle; NULL = random each run

# =============================================================================
# FUNCTION: split_csv
# Randomly shuffles and partitions a CSV into named subsets, then writes each
# subset to its own file in the specified output directory.
#
# Parameters:
#   input_file        – path to the source CSV
#   output_dir        – directory where output CSVs are written (created if absent)
#   base_name         – filename prefix shared by all output files
#   split_percentages – numeric vector of proportions that must sum to 1.0
#   split_names       – character vector of labels; one per proportion
#   seed              – optional integer for set.seed(); NULL means non-reproducible
# =============================================================================
split_csv <- function(input_file, output_dir, base_name, 
                      split_percentages, split_names,
                      seed = NULL) {
  # Validate inputs before doing any work
  if (abs(sum(split_percentages) - 1.0) > 0.001) stop("Split percentages must sum to 1.0")
  if (length(split_percentages) != length(split_names)) stop("Mismatch in split counts")
  
  # Fix the random shuffle if a seed was supplied
  if (!is.null(seed)) set.seed(seed)
  
  # Load data and record row count
  data <- read_csv(input_file)
  n_rows <- nrow(data)
  
  # Generate a random permutation of row indices (i.e. shuffle the dataset)
  indices <- sample(n_rows)
  
  # Convert proportions to integer row counts; floor() prevents over-allocation
  split_n <- floor(split_percentages * n_rows)
  split_n[length(split_n)] <- n_rows - sum(split_n[-length(split_n)])  # Assign any rows lost to rounding into the last split so totals always equal n_rows
  
  # Slice the shuffled index vector into contiguous blocks, one per split
  splits <- vector("list", length(split_n))
  start_idx <- 1
  for (i in seq_along(split_n)) {
    end_idx <- start_idx + split_n[i] - 1
    splits[[i]] <- data[indices[start_idx:end_idx], ]
    start_idx <- end_idx + 1
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Write each split to disk and report row count + path to the console
  for (i in seq_along(splits)) {
    file_path <- file.path(output_dir, paste0(base_name, "_", split_names[i], ".csv"))
    write_csv(splits[[i]], file_path)
    cat(sprintf("%s: %d rows -> %s\n", split_names[i], nrow(splits[[i]]), file_path))
  }
}

# =============================================================================
# Call the function using the configuration variables above
# =============================================================================
split_csv(input_csv, output_directory, file_base_name, split_percentages, split_names, random_seed)
