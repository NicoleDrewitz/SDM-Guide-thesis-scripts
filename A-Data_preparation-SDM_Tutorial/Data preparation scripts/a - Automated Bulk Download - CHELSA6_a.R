#===============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# CHELSA Bioclim+ Bulk Download
# Author: Nicole Drewitz, Drafted with Perplexity AI, 09-2025
# Last updated: 15-05-2026 with Claude AI

# Purpose: Downloads CHELSA climate rasters from a URL list, then clips to extent.
# Usage: You can create a new WGET file from this download source: https://www.chelsa-climate.org/datasets/chelsa_bioclim and put in "Input files" folder
# Slow loading? change download_dir file path to your desktop
# If the download was interrupted another script was created to continue the CHELSA download.

# Data citation:
# Brun, P., Zimmermann, N. E., Hari, C., Pellissier, L., Karger, D. N. (2022). 
# CHELSA-BIOCLIM+ A novel set of global climate-related predictors at kilometre-resolution. 
# EnviDat. https://www.doi.org/10.16904/envidat.332.
#===============================================================================

library(httr)    # HTTP downloads
library(tictoc)  # Timing
library(terra)   # Raster operations (clip/crop)
library(here)    # Relative file paths

tic()  # Start timer

# --- Configuration -----------------------------------------------------------
wget_file_path <- here("Input files/CHELSA6-bioclim download list - envidatS3paths - 2025-21-1.txt")
download_dir   <- here("Environmental_Variables/Bioclim_CHELSA") # if you are using this R project from cloud storage, files will download easier if you assign a file path to your desktop.

# Bounding box for Europe (includes Iceland, Greenland, Svalbard, Europe west of Urals)
study_area_extent <- c(xmin = -75, xmax = 60, ymin = 35, ymax = 85)

clip_files <- TRUE  # Set FALSE to skip the clipping step

# --- Setup -------------------------------------------------------------------
dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)
setwd(download_dir)

# --- Read & clean URL list ---------------------------------------------------
cat("Reading URLs...\n")
urls <- readLines(wget_file_path, warn = FALSE)
urls <- trimws(urls)
urls <- urls[urls != "" & !startsWith(urls, "#")]  # drop blanks and comments
urls <- urls[grepl("^https?://", urls)]            # keep only valid URLs
cat("Valid URLs found:", length(urls), "\n")

if (length(urls) == 0) {
  cat("First 10 lines of file:\n")
  print(readLines(wget_file_path, n = 10, warn = FALSE))
  stop("No valid URLs found — check file format.")
}

# =============================================================================
# DOWNLOAD (and optionally clip) FILES
# =============================================================================
cat("\n=== DOWNLOADING FILES ===\n")

# Counters
n_success <- 0
n_skipped <- 0
n_failed  <- 0
failed_files <- list()

n_clipped <- n_clip_skipped <- n_clip_failed <- 0
study_area_ext <- ext(study_area_extent)

# Clip a single downloaded file to European extent, delete the global original
clip_to_europe <- function(input_path) {
  input_file <- basename(input_path)
  
  # Build a shorter output name: drop "CHELSA_" prefix and "_V.2.1" version tag
  output_file <- input_file |>
    sub("^CHELSA_", "", x = _) |>
    sub("_V\\.2\\.1", "", x = _) |>
    sub("\\.tif$", "_eu.tif", x = _)
  output_path <- file.path(download_dir, output_file)
  
  if (file.exists(output_path)) {
    cat("  CLIP SKIP (exists):", output_file, "\n")
    return("skipped")
  }
  
  cat("  Clipping to Europe ...")
  
  tryCatch({
    r_crop <- crop(rast(input_path), study_area_ext)
    writeRaster(r_crop, output_path, overwrite = TRUE,
                gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
    
    orig_mb <- file.size(input_path)  / 1024^2
    clip_mb <- file.size(output_path) / 1024^2
    cat(sprintf(" OK (%.1f → %.1f MB, %.0f%% smaller)\n",
                orig_mb, clip_mb, (1 - clip_mb / orig_mb) * 100))
    
    file.remove(input_path)  # delete global file to save space
    return("success")
    
  }, error = function(e) {
    cat(" FAILED:", e$message, "\n")
    return("failed")
  })
}

# Download a single file; clips immediately if clip_files = TRUE
download_file <- function(url, timeout_seconds = 300) {
  filename   <- basename(url)
  local_path <- file.path(download_dir, filename)
  
  # Derive the clipped output name to check if already fully processed
  clipped_name <- filename |>
    sub("^CHELSA_", "", x = _) |>
    sub("_V\\.2\\.1", "", x = _) |>
    sub("\\.tif$", "_eu.tif", x = _)
  clipped_path <- file.path(download_dir, clipped_name)
  
  # Skip if clipped file already exists, or (if not clipping) raw file > 1 MB
  already_clipped <- clip_files && file.exists(clipped_path)
  already_raw     <- !clip_files && file.exists(local_path) && file.size(local_path) > 1e6
  
  if (already_clipped || already_raw) {
    label <- if (already_clipped) clipped_name else filename
    cat("SKIP:", label, "\n")
    return(list(status = "skipped", filename = filename))
  }
  
  # Remove any leftover partial file
  if (file.exists(local_path)) file.remove(local_path)
  
  tryCatch({
    cat("Downloading:", filename, "...")
    resp <- GET(url, timeout(timeout_seconds))
    
    if (status_code(resp) != 200) {
      cat(" FAILED (HTTP", status_code(resp), ")\n")
      return(list(status = "failed", filename = filename,
                  reason = paste("HTTP", status_code(resp))))
    }
    
    # Write to disk
    writeBin(content(resp, "raw"), local_path)
    actual_size   <- file.size(local_path)
    expected_size <- as.numeric(headers(resp)$`content-length`)  # may be NA
    
    # Validate: must be > 1 MB and match Content-Length if provided
    if (actual_size <= 1e6) {
      cat(" FAILED (too small:", actual_size, "bytes)\n")
      file.remove(local_path)
      return(list(status = "failed", filename = filename, reason = "file_too_small"))
    }
    
    if (!is.na(expected_size) && actual_size != expected_size) {
      cat(" FAILED (incomplete:", actual_size, "vs", expected_size, "bytes)\n")
      file.remove(local_path)
      return(list(status = "failed", filename = filename, reason = "incomplete_download"))
    }
    
    cat(sprintf(" OK (%.1f MB)\n", actual_size / 1024^2))
    
    # Clip immediately after successful download
    if (clip_files) clip_to_europe(local_path)
    
    return(list(status = "success", filename = filename, size_mb = actual_size / 1024^2))
    
  }, error = function(e) {
    cat(" ERROR:", e$message, "\n")
    if (file.exists(local_path)) file.remove(local_path)
    return(list(status = "failed", filename = filename, reason = e$message))
  })
}

# Download loop
for (i in seq_along(urls)) {
  cat(sprintf("\n[%d/%d  %.1f%%] ", i, length(urls), i / length(urls) * 100))
  
  result <- download_file(urls[[i]])
  
  if      (result$status == "success") { n_success <- n_success + 1; Sys.sleep(2) }
  else if (result$status == "skipped")   n_skipped <- n_skipped + 1
  else {
    n_failed <- n_failed + 1
    failed_files[[length(failed_files) + 1]] <- result
  }
  
  # Progress report every 10 files
  if (i %% 10 == 0) {
    cat(sprintf("\n--- Progress: %d done, %d skipped, %d failed, %d remaining ---\n",
                n_success, n_skipped, n_failed, length(urls) - i))
  }
}

cat(sprintf("\nDownload summary: %d OK | %d skipped | %d failed\n",
            n_success, n_skipped, n_failed))

# Save failed-download list for retry script
if (n_failed > 0) {
  writeLines(
    c("CHELSA files that failed to download:",
      strrep("=", 50), "",
      sapply(failed_files, function(x) paste(x$filename, "-", x$reason))),
    file.path(download_dir, "failed_downloads.txt")
  )
  cat("Failures written to failed_downloads.txt\n")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n=== FINAL SUMMARY ===\n")

pattern      <- if (clip_files) "_eu\\.tif$" else "\\.tif$"
final_files  <- list.files(download_dir, pattern = pattern)
total_size_gb <- sum(file.size(file.path(download_dir, final_files)), na.rm = TRUE) / 1024^3

cat("Files in directory:", length(final_files), "\n")
cat(sprintf("Total size: %.2f GB\n", total_size_gb))

# Show up to 5 sample files
cat("\nSample files:\n")
for (f in head(final_files, 5))
  cat(sprintf("  %s  (%.1f MB)\n", f, file.size(file.path(download_dir, f)) / 1024^2))
if (length(final_files) > 5) cat("  ... and", length(final_files) - 5, "more\n")

cat("\nOutput directory:", download_dir, "\n")

if (n_failed > 0)
  cat("\nNOTE: Some downloads failed — see failed_downloads.txt\n")

cat(if (clip_files) "\nDone: files downloaded and clipped.\n"
    else            "\nDone: files downloaded.\n")

toc()
