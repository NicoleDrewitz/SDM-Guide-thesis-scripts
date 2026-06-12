# =============================================================================
# RECREATE PARAMETER_TUNING_RESULTS.CSV FROM MAXENT OUTPUT FILES
# =============================================================================

# Set your base directory (same as original script)
base_dir <- "/Users/nicoledrewitz/Downloads"
# base_dir <- "C:/Users/nidre4892/Downloads"
setwd(base_dir)

output_dir <- file.path(base_dir, "MaxEnt_Results_parameter_tuning")

cat("=== Recreating Parameter Tuning Results ===\n")
cat("Searching for MaxEnt run directories in:", output_dir, "\n\n")

# Initialize results data frame
tuning_results <- data.frame(
  reg_mult = numeric(),
  features = character(),
  train_auc = numeric(),
  test_auc = numeric(),
  stringsAsFactors = FALSE
)

# Get all subdirectories in the output folder
run_dirs <- list.dirs(output_dir, full.names = TRUE, recursive = FALSE)

if (length(run_dirs) == 0) {
  stop("No run directories found in: ", output_dir)
}

cat("Found", length(run_dirs), "run directories\n\n")

# Process each directory
for (run_dir in run_dirs) {
  dir_name <- basename(run_dir)
  
  # Expected format: regmult_X_X_FEATURES (e.g., regmult_0_5_LQ or regmult_10_LQPTH)
  parts <- strsplit(dir_name, "_")[[1]]
  
  if (length(parts) >= 3 && parts[1] == "regmult") {
    # Construct reg_mult: combine parts[2] and parts[3] if fractional, e.g., "0_5" → 0.5
    if (grepl("^[0-9]+_[0-9]+$", paste(parts[2], parts[3], sep="_"))) {
      reg_mult_str <- paste(parts[2], parts[3], sep=".")
      reg_mult <- as.numeric(reg_mult_str)
      feature_index <- 4
    } else {
      reg_mult <- as.numeric(parts[2])
      feature_index <- 3
    }
    
    # Extract feature combination (everything after regmult)
    features <- paste(parts[feature_index:length(parts)], collapse = "_")
    
    cat("Processing:", dir_name, "\n")
    cat("  RegMult:", reg_mult, "Features:", features, "\n")
    
    results_file <- file.path(run_dir, "maxentResults.csv")
    
    if (file.exists(results_file)) {
      tryCatch({
        results_table <- read.csv(results_file, stringsAsFactors = FALSE)
        
        if (nrow(results_table) > 0) {
          train_cols <- grep("Training.AUC|AUC.training|Training AUC", 
                             names(results_table), ignore.case = TRUE, value = TRUE)
          test_cols <- grep("Test.AUC|AUC.test|Test AUC", 
                            names(results_table), ignore.case = TRUE, value = TRUE)
          
          train_auc <- if (length(train_cols) > 0) as.numeric(results_table[[train_cols[1]]][1]) else NA
          test_auc <- if (length(test_cols) > 0) as.numeric(results_table[[test_cols[1]]][1]) else NA
          
          tuning_results <- rbind(tuning_results, data.frame(
            reg_mult = reg_mult,
            features = features,
            train_auc = train_auc,
            test_auc = test_auc,
            stringsAsFactors = FALSE
          ))
          
          cat("  SUCCESS - Train AUC:", round(train_auc, 3),
              "Test AUC:", round(test_auc, 3), "\n\n")
        } else {
          cat("  WARNING: Results file is empty\n\n")
        }
      }, error = function(e) {
        cat("  ERROR reading results file:", conditionMessage(e), "\n\n")
      })
    } else {
      cat("  WARNING: maxentResults.csv not found\n\n")
    }
  } else {
    cat("WARNING: Could not parse directory name:", dir_name, "\n\n")
  }
}


# Sort results by regularization multiplier and features
tuning_results <- tuning_results[order(tuning_results$reg_mult, tuning_results$features), ]

# Save recreated results
output_file <- file.path(output_dir, "parameter_tuning_results_recreated.csv")
write.csv(tuning_results, output_file, row.names = FALSE)

cat("\n=== Summary ===\n")
cat("Total runs recovered:", nrow(tuning_results), "\n")
cat("Regularization multipliers:", length(unique(tuning_results$reg_mult)), "\n")
cat("Feature combinations:", length(unique(tuning_results$features)), "\n")
cat("\nResults saved to:", output_file, "\n")

# Display top models by test AUC
if (nrow(tuning_results) > 0) {
  cat("\n=== Top 10 Models by Test AUC ===\n")
  top_models <- tuning_results[order(-tuning_results$test_auc), ]
  print(head(top_models, 10))
  
  cat("\n=== All Results ===\n")
  print(tuning_results)
  
  # Check for any missing values
  missing_train <- sum(is.na(tuning_results$train_auc))
  missing_test <- sum(is.na(tuning_results$test_auc))
  
  if (missing_train > 0 || missing_test > 0) {
    cat("\n=== WARNING ===\n")
    cat("Missing Train AUC values:", missing_train, "\n")
    cat("Missing Test AUC values:", missing_test, "\n")
  }
}