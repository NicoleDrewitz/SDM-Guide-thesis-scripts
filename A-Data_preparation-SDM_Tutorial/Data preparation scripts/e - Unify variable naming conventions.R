# =============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Batch file renaming script
# Create by Nicole Drewitz with Claude.ai (Sonnet 4.5)

# Removes specific prefixes and suffixes from file names
# Needed for MaxEnt to recognize climate ssp scenario projections as the same environmental variable to make future predictions
# Same method could be used for hindcasting habitat suitability (could be another method for model validation if historic species range is known)
# Last updated: 2025-01-15
# =============================================================================

library(here)

# Set the main folder path
main_folder <- here::here("Environmental_Variables/Environmental variable dataset for modelling")

# Define prefixes and suffixes to remove ====
prefixes_to_remove <- c("gfdl-esm4_ssp126_", "gfdl-esm4_ssp370_", "gfdl-esm4_ssp585_")
suffixes_to_remove <- c("_2071-2100_eu", "_1981-2010_eu")

# Get all files recursively from the folder and subfolders ====
all_files <- list.files(main_folder, full.names = TRUE, recursive = TRUE)

# Counter for renamed files
renamed_count <- 0

# Loop through each file ====
for (file_path in all_files) {
  # Get the directory and file name
  dir_name <- dirname(file_path)
  file_name <- basename(file_path)
  
  # Store original name for comparison
  original_name <- file_name
  
  # Remove prefixes ====
  for (prefix in prefixes_to_remove) {
    if (startsWith(file_name, prefix)) {
      file_name <- sub(paste0("^", prefix), "", file_name)
    }
  }
  
  # Remove suffixes (before file extension) ====
  file_ext <- tools::file_ext(file_name)
  file_base <- tools::file_path_sans_ext(file_name)
  
  for (suffix in suffixes_to_remove) {
    if (endsWith(file_base, suffix)) {
      file_base <- sub(paste0(suffix, "$"), "", file_base)
    }
  }
  
  # Reconstruct filename with extension ====
  if (nchar(file_ext) > 0) {
    new_file_name <- paste0(file_base, ".", file_ext)
  } else {
    new_file_name <- file_base
  }
  
  # Rename if the name has changed
  if (original_name != new_file_name) {
    new_file_path <- file.path(dir_name, new_file_name)
    
    # Check if the new file name already exists
    if (!file.exists(new_file_path)) {
      file.rename(file_path, new_file_path)
      cat("Renamed:", original_name, "->", new_file_name, "\n")
      renamed_count <- renamed_count + 1
    } else {
      cat("Warning: Cannot rename", original_name, "- target file already exists\n")
    }
  }
}

cat("\nTotal files renamed:", renamed_count, "\n")