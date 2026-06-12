#===============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Automated WEKEO Downloading
# Accessing Copernicus data often results in multiple files for download (e.g. spatially separated tiles)
# Created on 16-09-2025 with Perplexity by Nicole Drewitz
#===============================================================================

# Install and load required packages
if(!require("hdar")) install.packages("hdar")
library(hdar)
if(!require("rstudioapi")) install.packages("rstudioapi")
library(rstudioapi)

# Log into WEkEO ====
username <- askForPassword("WEkEO Username")
password <- askForPassword("WEkEO Password")
client <- Client$new(username, password, save_credentials = TRUE)
#If error occurs, enter login again, try later, or reset your password
client$get_token()
client$show_terms()
client$terms_and_conditions()
client$terms_and_conditions(term_id = "all")

# Choose extents ====
# East area (over WEkEO download limit, only downloads northern Fennoscandia)
bbox_whole = c(-73.73951958338158,54.05662350026431,33.21129872285927,83.83000823746343)
# Iceland and Faroes
bbox_IceFar = c(-24.78477358310494, 61.09766101464476, -5.80369378052789, 66.63516657822314)

# South area (southern Fennoscandia + Estonia)
bbox_south <- c(4.203523078137387, 54.0844134626607, 31.77775648761179, 63.63127872571266)
# North area (north of N53)

bbox_NW <- c(-74, -20.5, 69, 84)
bbox_NE <- c(-20.5, 33, 69, 84)
bbox_SW <- c(-74, -20.5, 54, 69)
bbox_SE <- c(-20.5, 33, 54, 69)

# API request (forest cover density) ====
query <- list(
 dataset_id = "EO:EEA:DAT:HRL:TCF",
 bbox = bbox_south,
  productType = "Tree Cover Density",
  resolution = "100m",
 year = "2021",
  itemsPerPage = 200,
 startIndex = 0)

# API request (forest type)
#query <- list(
#  dataset_id = "EO:EEA:DAT:HRL:TCF",
#  bbox = bbox_NW,
#  productType = "Forest Type",
#  resolution = "100m",
#  year = "2021",
#  itemsPerPage =200,
#  startIndex =0)

# Convert query to JSON (if needed) ====
library(jsonlite)
json_query <- toJSON(query, auto_unbox = TRUE)

# Search HDA service
matches <- client$search(json_query, limit=100)
# Show found files and IDs
cat("Found", length(matches$results), "files\n")
print(sapply(matches$results, function(x) x$id))

# Extract all file IDs ====
all_file_ids <- sapply(matches$results, function(x) x$id)
# Extract unique IDs from the results
unique_ids <- unique(sapply(matches$results, function(x) x$id))
cat("Number of unique file IDs:", length(unique_ids), "\n")

# Bulk download to directory ====
output_directory <- here::here("Environmental_variables/Tree_Cover_Density_South")
#output_directory <- here::here("Environmental_Variables/Forest_Type_NW")
# Create download directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)}

#matches$download(output_directory, force = FALSE)
# Download each unique file by its ID
for(file_id in unique_ids) {
  cat("Downloading file ID:", file_id, "\n")
  client$download(id = file_id, path = output_directory)
}

# Replace your problematic download loop with this: ====
matches$download(output_directory, force = FALSE)

# List downloaded files
print(list.files(output_directory))

# List already downloaded files in output directory
downloaded_files <- list.files(output_directory)

#_________________________________________________
# Identify which files are missing (not downloaded or failed previously) ====
# Assuming downloaded file names contain the IDs, adjust if necessary

missing_files_ids <- setdiff(all_file_ids, downloaded_files)

if (length(missing_files_ids) == 0) {
  cat("All files already downloaded.\n")
} else {
  cat("Retrying download for", length(missing_files_ids), "missing files...\n")
  missing_matches <- Filter(function(x) x$id %in% missing_files_ids, matches$results)
  
  for (item in missing_matches) {
    file_id <- item$id
    cat("Downloading file ID:", file_id, "...\n")
    tryCatch({
      item$download(output_directory, force = FALSE)
      cat("Downloaded:", file_id, "\n")
    }, error = function(e) {
      cat("Failed to download", file_id, ":", e$message, "\n")
    })
  }
}

cat("Finished downloads. Files in directory:\n")
print(list.files(output_directory))

# Re-authenticate to refresh connection ====
client <- Client$new(username, password)
client$get_token()
client$terms_and_conditions(term_id = "all")