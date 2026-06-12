# ============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# RETRY FAILED DOWNLOADS / CHELSA Bioclim+ Bulk Download
# Paste this at the end of your script or run separately
# ============================================================================

library(httr)
library (here)

# Configuration - SAME AS YOUR MAIN SCRIPT
download_dir <- here::here("Environmental_Variables/Bioclim_CHELSA")
setwd(download_dir)

cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
cat("RETRYING FAILED DOWNLOADS\n")
cat(paste(rep("=", 60), collapse=""), "\n", sep="")

# Check if failed_downloads.txt exists
failed_list_file <- file.path(download_dir, "failed_downloads.txt")

if (!file.exists(failed_list_file)) {
  cat("No failed_downloads.txt file found.\n")
  cat("Either all downloads succeeded or the file was deleted.\n")
} else {
  # Read the failed downloads file
  failed_lines <- readLines(failed_list_file, warn = FALSE)
  
  # Extract only lines that start with "CHELSA_" (these are the filenames)
  # Skip header lines, empty lines, and error message continuation lines
  failed_lines <- failed_lines[grepl("^CHELSA_", failed_lines)]
  
  # Extract just the filename part (before the " - " reason)
  failed_filenames <- sapply(strsplit(failed_lines, " - "), function(x) trimws(x[1]))
  
  # Remove any duplicates
  failed_filenames <- unique(failed_filenames)
  
  cat("Found", length(failed_filenames), "failed downloads to retry\n\n")
  
  if (length(failed_filenames) > 0) {
    # Reconstruct URLs from filenames
    # CHELSA URLs follow a predictable pattern
    retry_urls <- character(0)
    
    for (filename in failed_filenames) {
      # Parse filename to reconstruct URL
      # Example: CHELSA_bio01_1981-2010_V.2.1.tif
      # Example: CHELSA_gfdl-esm4_ssp126_bio01_2071-2100_V.2.1.tif
      
      if (grepl("_gfdl-esm4_", filename)) { #GFDL-ESM4 is a model of bioclimatic prediction variables
        # Future projection file
        parts <- strsplit(filename, "_")[[1]]
        model <- "GFDL-ESM4"
        ssp <- parts[3]  # ssp126, ssp370, or ssp585
        variable <- parts[4]
        period <- parts[5]
        
        url <- paste0("https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/",
                      variable, "/", period, "/", model, "/", ssp, "/", filename)
      } else {
        # Historical file (1981-2010)
        parts <- strsplit(filename, "_")[[1]]
        variable <- parts[2]
        period <- parts[3]
        
        url <- paste0("https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/",
                      variable, "/", period, "/", filename)
      }
      
      retry_urls <- c(retry_urls, url)
    }
    
    cat("Reconstructed", length(retry_urls), "URLs\n")
    cat("Starting retry...\n\n")
    
    # Retry counters
    retry_success <- 0
    retry_failed <- 0
    still_failed <- list()
    
    # Function to download a single file (same as main script)
    download_file <- function(url, timeout_seconds = 300) {
      filename <- basename(url)
      local_path <- file.path(download_dir, filename)
      
      # Check if file now exists
      if (file.exists(local_path) && file.size(local_path) > 1000000) {
        file_size_mb <- round(file.size(local_path) / 1024^2, 2)
        cat("SKIP:", filename, "(now exists -", file_size_mb, "MB)\n")
        return(list(status = "skipped", filename = filename))
      }
      
      # Remove any existing partial file
      if (file.exists(local_path)) {
        file.remove(local_path)
      }
      
      # Try download
      tryCatch({
        cat("Downloading:", filename, "...")
        
        response <- GET(url, timeout(timeout_seconds))
        
        if (status_code(response) == 200) {
          content_length <- headers(response)$`content-length`
          expected_size <- if(!is.null(content_length)) as.numeric(content_length) else NULL
          
          file_content <- content(response, "raw")
          writeBin(file_content, local_path)
          
          actual_size <- file.size(local_path)
          
          if (file.exists(local_path) && actual_size > 1000000) {
            if (!is.null(expected_size) && actual_size != expected_size) {
              cat(" FAILED (incomplete download)\n")
              file.remove(local_path)
              return(list(status = "failed", filename = filename, reason = "incomplete_download"))
            }
            
            file_size_mb <- round(actual_size / 1024^2, 2)
            cat(" SUCCESS (", file_size_mb, "MB)\n")
            return(list(status = "success", filename = filename, size_mb = file_size_mb))
          } else {
            cat(" FAILED (file too small)\n")
            if (file.exists(local_path)) file.remove(local_path)
            return(list(status = "failed", filename = filename, reason = "file_too_small"))
          }
        } else {
          cat(" FAILED (HTTP", status_code(response), ")\n")
          return(list(status = "failed", filename = filename, reason = paste("HTTP", status_code(response))))
        }
        
      }, error = function(e) {
        cat(" ERROR:", e$message, "\n")
        if (file.exists(local_path)) file.remove(local_path)
        return(list(status = "failed", filename = filename, reason = e$message))
      })
    }
    
    # Retry each failed download
    for (i in 1:length(retry_urls)) {
      url <- retry_urls[i]
      
      progress_pct <- round(i/length(retry_urls)*100, 1)
      cat("\n[", i, "/", length(retry_urls), " - ", progress_pct, "%] ")
      
      result <- download_file(url, timeout_seconds = 300)
      
      if (result$status == "success") {
        retry_success <- retry_success + 1
      } else if (result$status == "failed") {
        retry_failed <- retry_failed + 1
        still_failed[[length(still_failed) + 1]] <- result
      }
      
      if (result$status == "success") {
        Sys.sleep(2)
      }
    }
    
    # Summary
    cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
    cat("RETRY SUMMARY\n")
    cat(paste(rep("=", 60), collapse=""), "\n", sep="")
    cat("Attempted:", length(retry_urls), "\n")
    cat("Successfully downloaded:", retry_success, "\n")
    cat("Still failed:", retry_failed, "\n")
    
    # Update failed downloads file
    if (retry_failed > 0) {
      failed_df <- data.frame(
        filename = sapply(still_failed, function(x) x$filename),
        reason = sapply(still_failed, function(x) x$reason),
        stringsAsFactors = FALSE
      )
      
      writeLines(
        c("CHELSA Files That Could Not Be Downloaded (After Retry):",
          paste(rep("=", 50), collapse=""),
          "",
          paste(failed_df$filename, "-", failed_df$reason)),
        failed_list_file
      )
      cat("\nUpdated failed_downloads.txt with remaining", retry_failed, "failed files\n")
    } else {
      # All succeeded - delete the failed file
      file.remove(failed_list_file)
      cat("\n✓ All failed downloads now succeeded! Deleted failed_downloads.txt\n")
    }
  }
}

cat("\nRetry complete!\n")