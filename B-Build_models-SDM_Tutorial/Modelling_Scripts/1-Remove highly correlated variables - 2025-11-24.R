# =============================================================================
# Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Configurable Multicollinearity Analysis for Nordic Cloudberry (Rubus chamaemorus)
# by Nicole Drewitz with Claude AI in October 2025
# last updated 24 November 2025 (with lowest_importance_first option)
# =============================================================================
library(here)
# =============================================================================
# CONFIGURATION SECTION - MODIFY THESE SETTINGS
# =============================================================================

CONFIG <- list(
  # Variable importance metric to use
  # Options: "contribution", "permutation.importance"
  importance_metric = "permutation.importance",

  # Correlation pair priority strategy
  # Options: "highest_correlation_first" (address most strongly correlated pairs first)
  #          "lowest_correlation_first" (address weakly correlated pairs first)
  #          "lowest_importance_first" (address pairs whose removal candidate has lowest importance first)
  # Note: Within each pair, the variable with LOWER importance is always removed
  pair_priority = "lowest_importance_first",

  # Correlation threshold for identifying problematic pairs
  correlation_threshold = 0.75,

  # Number of cells to sample for correlation analysis
  sample_size = 50000,

  # Random seed for reproducibility
  random_seed = 123,

  # VIF threshold for warnings
  vif_warning_threshold = 5,

  # Correlation method (will be auto-selected based on normality test if NULL)
  # Options: "pearson", "spearman", NULL (auto-detect)
  correlation_method = NULL
)

# Directory paths
PATHS <- list(
  env_raster_dir = here::here("Input_files/Environmental_Variables-processed-current_climate"),
  base_dir = here::here("Output_files/Reducing_environmental_variable_set"),
  maxent_output_dir = here::here("MaxEnt-TrialRUNS-untuned_settings/MaxEnt_Results_Single_Run") # assign your MaxEnt untuned settings RUN folder
)

# OPTIONAL: Variables to exclude from analysis
# Examples: c("bio1", "bio2") or c() for no exclusions
# Can use prefixes (check spelling)
EXCLUDE_VARIABLES <- c("soil_silt") # silt+clay+sand =100% so all 3 are unnecessary

# =============================================================================
# DISPLAY CONFIGURATION
# =============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("=== CONFIGURATION SETTINGS ===\n")
cat(rep("=", 80), "\n\n", sep = "")
cat("Importance metric:        ", CONFIG$importance_metric, "\n")
cat("Pair priority:            ", CONFIG$pair_priority, "\n")
cat("Correlation threshold:    ", CONFIG$correlation_threshold, "\n")
cat("Sample size:              ", CONFIG$sample_size, "\n")
cat("Random seed:              ", CONFIG$random_seed, "\n")
cat("VIF warning threshold:    ", CONFIG$vif_warning_threshold, "\n")
cat("Correlation method:       ", ifelse(is.null(CONFIG$correlation_method), "Auto-detect", CONFIG$correlation_method), "\n")
cat("\n")

# =============================================================================
# LOAD REQUIRED PACKAGES
# =============================================================================

cat("=== Loading Required Packages ===\n")
if (!requireNamespace("corrplot", quietly = TRUE)) install.packages("corrplot")
if (!requireNamespace("usdm", quietly = TRUE)) install.packages("usdm")
library(corrplot)
library(usdm)
library(raster)
library(terra)
cat("Packages loaded successfully\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Convert variable names (replace all hyphens with periods to match MaxEnt format)
convert_var_name <- function(var_name) {
  gsub("-", ".", var_name)
}

# Shorten variable names at underscore for plotting
shorten_var_names <- function(var_names) {
  sapply(var_names, function(name) {
    parts <- strsplit(name, "_")[[1]]
    if (length(parts) > 1) {
      paste(parts[1], parts[2], sep = "_")
    } else {
      name
    }
  }, USE.NAMES = FALSE)
}

# =============================================================================
# LOAD ENVIRONMENTAL VARIABLES
# =============================================================================

cat("=== Loading Environmental Variables ===\n")
raster_files <- list.files(PATHS$env_raster_dir, pattern = "\\.asc$",
                           full.names = TRUE, ignore.case = TRUE)
if (length(raster_files) == 0) {
  stop("No .asc files found in environmental variables directory!")
}

initial_env_raster <- rast(raster_files)

# Assign names immediately after loading
initial_var_names <- tools::file_path_sans_ext(basename(raster_files))
names(initial_env_raster) <- initial_var_names

# Exclude variables by prefix from your exclusion list
if (length(EXCLUDE_VARIABLES) > 0) {
  exclude_pattern <- paste0("^", EXCLUDE_VARIABLES, collapse = "|")
  keep_layers <- !grepl(exclude_pattern, names(initial_env_raster))
  initial_env_raster <- initial_env_raster[[keep_layers]]
  cat("Excluded variables by prefix:", paste(EXCLUDE_VARIABLES, collapse = ", "), "\n")
}

cat("Loaded", nlyr(initial_env_raster), "environmental layers\n\n")

# =============================================================================
# CATEGORICAL VARIABLE IDENTIFICATION
# =============================================================================

cat("=== Identifying Categorical Variables ===\n")
categorical_prefixes <- c(
  "kg0", "Forest_type", "Peat", "Water", "Soil_Depth",
  "Topo_geomflat", "Topo_geomhollow", "Topo_geompit",
  "Topo_geomslope", "Topo_geomvalley"
)

categorical_vars <- unlist(
  lapply(categorical_prefixes, function(prefix) {
    initial_var_names[startsWith(initial_var_names, prefix)]
  })
)

cat("Identified", length(categorical_vars), "categorical variables\n")
if (length(categorical_vars) > 0) {
  cat("Categorical variables:\n", paste("  -", categorical_vars, collapse = "\n"), "\n")
}
cat("\n")

# =============================================================================
# DUPLICATE FILENAME CHECK
# =============================================================================

cat("=== Checking for Duplicate Filenames ===\n")
file_basenames <- tools::file_path_sans_ext(basename(raster_files))
if (any(duplicated(file_basenames))) {
  cat("ERROR: Duplicate .asc filenames exist. Please remove duplicates.\n")
  stop("Cannot proceed with duplicates.")
}
cat("No duplicate filenames found\n\n")

# Separate continuous and categorical variables
env_vars <- initial_env_raster
categorical_vars <- categorical_vars[categorical_vars %in% names(env_vars)]
continuous_vars <- setdiff(names(env_vars), categorical_vars)
cat("Continuous variables:", length(continuous_vars), "\n")
cat("Categorical variables:", length(categorical_vars), "\n\n")
env_vars_continuous <- env_vars[[continuous_vars]]

# =============================================================================
# EXCLUDE VARIABLES (IF ASSIGNED)
# =============================================================================

# Function to exclude variables by prefix
exclude_by_prefix <- function(var_names, exclude_prefixes) {
  pattern <- paste0("^", exclude_prefixes, collapse = "|")
  keep <- !grepl(pattern, var_names)
  return(var_names[keep])
}

# Exclude from variable lists
continuous_vars <- exclude_by_prefix(continuous_vars, EXCLUDE_VARIABLES)
categorical_vars <- exclude_by_prefix(categorical_vars, EXCLUDE_VARIABLES)

# =============================================================================
# SAMPLE RASTER CELLS FOR CORRELATION ANALYSIS
# =============================================================================

cat("=== Sampling Raster Cells for Correlation Analysis ===\n")
set.seed(CONFIG$random_seed)
n_cells <- min(CONFIG$sample_size, ncell(env_vars_continuous))
sample_cells <- sample(1:ncell(env_vars_continuous), n_cells)
env_sample <- as.data.frame(env_vars_continuous[sample_cells])
env_sample <- na.omit(env_sample)
cat("Using", nrow(env_sample), "cells for correlation analysis\n\n")

# =============================================================================
# CALCULATE CORRELATION MATRIX
# =============================================================================

cat("=== Calculating Correlation Matrix ===\n")

# Determine correlation method
if (is.null(CONFIG$correlation_method)) {
  normality_test <- apply(env_sample, 2, function(x) {
    if(length(x) > 3) {
      test_data <- if(length(x) > 5000) sample(x, 1000) else x
      shapiro.test(test_data)$p.value > 0.05
    } else TRUE
  })
  cor_method <- if (all(normality_test)) "pearson" else "spearman"
  cat("Auto-detected correlation method:", cor_method, "\n")
} else {
  cor_method <- CONFIG$correlation_method
  cat("Using specified correlation method:", cor_method, "\n")
}

cor_matrix <- cor(env_sample, use = "complete.obs", method = cor_method)

# Calculate total possible pairs
n_vars <- ncol(cor_matrix)
total_pairs <- (n_vars * (n_vars - 1)) / 2
cat("Total variable pairs evaluated:", total_pairs, "\n\n")

cat("=== Generating Initial Correlation Heatmap ===\n")

# Plot correlation matrix as heatmap
plot_correlation_heatmap <- function(cor_matrix, title, method) {
  short_names <- shorten_var_names(rownames(cor_matrix))
  cor_matrix_plot <- cor_matrix
  rownames(cor_matrix_plot) <- short_names
  colnames(cor_matrix_plot) <- short_names

  corrplot(cor_matrix_plot,
           method = "color",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.8,
           col = colorRampPalette(c("#0571b0", "white", "#ca0020"))(200),
           #addCoef.col = "black", # uncomment to overlay values on heatmap
           #number.cex = 0.5,
           cl.cex = 0.8,
           title = paste0(title, " (", method, ")"),
           mar = c(0, 0, 2, 0))
}
plot_correlation_heatmap(cor_matrix, "Initial Correlation Matrix", cor_method)

# =============================================================================
# IDENTIFY HIGH CORRELATION PAIRS
# =============================================================================

cat("\n=== Identifying High Correlation Pairs ===\n")
high_cor_pairs <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(ncol(cor_matrix)-1)) {
  for (j in (i+1):ncol(cor_matrix)) {
    if (abs(cor_matrix[i, j]) > CONFIG$correlation_threshold) {
      high_cor_pairs <- rbind(high_cor_pairs, data.frame(
        Variable1 = rownames(cor_matrix)[i],
        Variable2 = colnames(cor_matrix)[j],
        Correlation = cor_matrix[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("Found", nrow(high_cor_pairs), "pairs with |r| >", CONFIG$correlation_threshold, "\n\n")

# =============================================================================
# LOAD MAXENT VARIABLE IMPORTANCE
# =============================================================================

contributions <- NULL
maxent_loaded <- FALSE

if (!is.null(PATHS$maxent_output_dir) && dir.exists(PATHS$maxent_output_dir)) {
  cat("=== Loading MaxEnt Variable Importance ===\n")
  cat("MaxEnt output directory:", PATHS$maxent_output_dir, "\n")

  results_file <- list.files(PATHS$maxent_output_dir,
                             pattern = "maxentResults.csv",
                             full.names = TRUE,
                             ignore.case = TRUE,
                             recursive = TRUE)

  if (length(results_file) > 0) {
    cat("Found maxentResults.csv at:", results_file[1], "\n")
    maxent_results <- read.csv(results_file[1], stringsAsFactors = FALSE)

    # Build column pattern based on config
    col_pattern <- paste0("\\.", gsub("\\.", "\\\\.", CONFIG$importance_metric), "$")
    importance_cols <- grep(col_pattern, names(maxent_results), value = TRUE)

    if (length(importance_cols) > 0) {
      contributions <- data.frame(
        Variable = gsub(col_pattern, "", importance_cols),
        Importance = as.numeric(maxent_results[1, importance_cols]),
        stringsAsFactors = FALSE
      )
      contributions <- contributions[!is.na(contributions$Importance), ]
      contributions$original_order <- 1:nrow(contributions)

      cat("Successfully loaded", CONFIG$importance_metric, "for",
          nrow(contributions), "variables\n")

      contributions_sorted <- contributions[order(contributions$Importance,
                                                  decreasing = TRUE), ]
      cat("\nTop 10 variables by", CONFIG$importance_metric, ":\n")
      print(head(contributions_sorted[, c("Variable", "Importance")], 10),
            row.names = FALSE)
      maxent_loaded <- TRUE
    } else {
      cat("WARNING: Could not find columns matching pattern:", col_pattern, "\n")
    }
  } else {
    cat("WARNING: maxentResults.csv not found\n")
  }
}
cat("\n")

# =============================================================================
# VARIABLE REMOVAL ALGORITHM
# =============================================================================

if (maxent_loaded && nrow(high_cor_pairs) > 0) {

  cat(rep("=", 80), "\n", sep = "")
  cat("=== ANALYZING VARIABLES IN CORRELATED PAIRS ===\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Get all variables involved in high correlations
  all_vars_in_pairs <- unique(c(high_cor_pairs$Variable1, high_cor_pairs$Variable2))
  cat("Variables involved in high correlations:", length(all_vars_in_pairs), "\n")

  # Add importance data
  vars_with_importance <- data.frame(
    Variable = all_vars_in_pairs,
    stringsAsFactors = FALSE
  )

  vars_with_importance$Importance <- sapply(vars_with_importance$Variable, function(v) {
    v_converted <- convert_var_name(v)
    idx <- which(contributions$Variable == v_converted)
    if (length(idx) > 0) contributions$Importance[idx[1]] else NA
  })

  vars_with_importance$Order <- sapply(vars_with_importance$Variable, function(v) {
    v_converted <- convert_var_name(v)
    idx <- which(contributions$Variable == v_converted)
    if (length(idx) > 0) contributions$original_order[idx[1]] else NA
  })

  # Sort by importance for display (lowest to highest)
  vars_with_importance <- vars_with_importance[order(
    vars_with_importance$Importance,
    -vars_with_importance$Order
  ), ]

  cat("\n=== Variables Ranked by", CONFIG$importance_metric, "===\n")
  cat("(For reference - actual removal is pair-based)\n\n")
  print(vars_with_importance, row.names = FALSE)

  # =============================================================================
  # GREEDY REMOVAL ALGORITHM - ADDRESS PAIRS BY CONFIGURED PRIORITY
  # =============================================================================

  cat("\n", rep("=", 80), "\n", sep = "")
  cat("=== APPLYING VARIABLE REMOVAL ALGORITHM ===\n")
  cat("Pair priority strategy:", CONFIG$pair_priority, "\n")
  cat("Within each pair: Variable with LOWER", CONFIG$importance_metric, "is removed\n")
  cat(rep("=", 80), "\n\n", sep = "")

  vars_to_remove <- character(0)
  remaining_pairs <- high_cor_pairs
  removal_log <- data.frame(
    Variable = character(),
    Importance = numeric(),
    PairCorrelation = numeric(),
    N_Pairs_Eliminated = integer(),
    Reason = character(),
    stringsAsFactors = FALSE
  )

  step <- 0
  while (nrow(remaining_pairs) > 0 && step < 100) {
    step <- step + 1

    # Compute auxiliary columns for sorting and removal decisions
    remaining_pairs$AbsCorr <- abs(remaining_pairs$Correlation)

    if (CONFIG$pair_priority == "lowest_importance_first") {
      # Find importance per variable in the pairs
      remaining_pairs$Var1Importance <- sapply(remaining_pairs$Variable1, function(v) {
        v_converted <- convert_var_name(v)
        idx <- which(contributions$Variable == v_converted)
        if (length(idx) > 0) contributions$Importance[idx[1]] else Inf
      })
      remaining_pairs$Var2Importance <- sapply(remaining_pairs$Variable2, function(v) {
        v_converted <- convert_var_name(v)
        idx <- which(contributions$Variable == v_converted)
        if (length(idx) > 0) contributions$Importance[idx[1]] else Inf
      })

      # Variable with lower importance in the pair defined here
      remaining_pairs$LowerImportance <- pmin(remaining_pairs$Var1Importance, remaining_pairs$Var2Importance)

      # Sort pairs by ascending lower importance (lowest importance variables get removed first)
      # Tie-break with descending absolute correlation to remove highest correlation pairs with low importance first
      remaining_pairs <- remaining_pairs[order(remaining_pairs$LowerImportance, -remaining_pairs$AbsCorr), ]

    } else if (CONFIG$pair_priority == "highest_correlation_first") {
      remaining_pairs <- remaining_pairs[order(-remaining_pairs$AbsCorr), ]
    } else if (CONFIG$pair_priority == "lowest_correlation_first") {
      remaining_pairs <- remaining_pairs[order(remaining_pairs$AbsCorr), ]
    } else {
      # Default fallback to highest correlation first
      remaining_pairs <- remaining_pairs[order(-remaining_pairs$AbsCorr), ]
    }

    # Select the first pair for processing
    target_pair <- remaining_pairs[1, ]
    var1 <- target_pair$Variable1
    var2 <- target_pair$Variable2
    pair_corr <- target_pair$Correlation

    # Retrieve importance and order
    var1_converted <- convert_var_name(var1)
    var2_converted <- convert_var_name(var2)

    idx1 <- which(contributions$Variable == var1_converted)
    idx2 <- which(contributions$Variable == var2_converted)

    importance1 <- if (length(idx1) > 0) contributions$Importance[idx1[1]] else NA
    importance2 <- if (length(idx2) > 0) contributions$Importance[idx2[1]] else NA
    order1 <- if (length(idx1) > 0) contributions$original_order[idx1[1]] else NA
    order2 <- if (length(idx2) > 0) contributions$original_order[idx2[1]] else NA

    # Determine which variable to remove (lowest importance, then later order)
    if (is.na(importance1) && is.na(importance2)) {
      var_to_remove <- var1
      reason <- sprintf("Both missing importance (pair r=%.3f)", pair_corr)
      var_importance <- NA
    } else if (is.na(importance1)) {
      var_to_remove <- var1
      reason <- sprintf("Missing importance vs %.2f%% (pair r=%.3f)", importance2, pair_corr)
      var_importance <- NA
    } else if (is.na(importance2)) {
      var_to_remove <- var2
      reason <- sprintf("Missing importance vs %.2f%% (pair r=%.3f)", importance1, pair_corr)
      var_importance <- NA
    } else {
      if (importance1 < importance2) {
        var_to_remove <- var1
        var_importance <- importance1
        reason <- sprintf("Lower %s: %.2f%% vs %.2f%% (pair r=%.3f)", CONFIG$importance_metric, importance1, importance2, pair_corr)
      } else if (importance2 < importance1) {
        var_to_remove <- var2
        var_importance <- importance2
        reason <- sprintf("Lower %s: %.2f%% vs %.2f%% (pair r=%.3f)", CONFIG$importance_metric, importance2, importance1, pair_corr)
      } else {
        if (order1 > order2) {
          var_to_remove <- var1
          var_importance <- importance1
          reason <- sprintf("Tied %s (%.2f%%), later order (pair r=%.3f)", CONFIG$importance_metric, importance1, pair_corr)
        } else {
          var_to_remove <- var2
          var_importance <- importance2
          reason <- sprintf("Tied %s (%.2f%%), later order (pair r=%.3f)", CONFIG$importance_metric, importance2, pair_corr)
        }
      }
    }

    # Count how many pairs this removal eliminates
    n_pairs_removed <- sum(remaining_pairs$Variable1 == var_to_remove | remaining_pairs$Variable2 == var_to_remove)

    removal_log <- rbind(removal_log, data.frame(
      Variable = var_to_remove,
      Importance = var_importance,
      PairCorrelation = pair_corr,
      N_Pairs_Eliminated = n_pairs_removed,
      Reason = reason,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("Step %d: Pair [%s - %s], r=%.3f\n", step, var1, var2, pair_corr))
    cat(sprintf("        Removing %-40s (%s)\n", var_to_remove, reason))
    cat(sprintf("        Eliminates %d total pairs\n\n", n_pairs_removed))

    # Remove variable and update pairs
    vars_to_remove <- c(vars_to_remove, var_to_remove)
    remaining_pairs <- remaining_pairs[
      remaining_pairs$Variable1 != var_to_remove &
        remaining_pairs$Variable2 != var_to_remove,
    ]
  }

  cat("\n=== REMOVAL COMPLETE ===\n")
  cat("Total variables removed:", length(vars_to_remove), "\n")
  cat("Total steps:", step, "\n\n")

  # =============================================================================
  # FINAL RESULTS
  # =============================================================================

  vars_to_keep <- setdiff(continuous_vars, vars_to_remove)

  cat(rep("=", 80), "\n", sep = "")
  cat("=== FINAL REMOVAL SUMMARY ===\n")
  cat(rep("=", 80), "\n\n", sep = "")

  cat("Variables removed (in order):\n\n")
  print(removal_log, row.names = FALSE)

  cat("\n=== FINAL VARIABLE SET ===\n")
  cat("Continuous variables retained:", length(vars_to_keep), "\n")
  cat("Continuous variables removed:", length(vars_to_remove), "\n")
  cat("Categorical variables (unchanged):", length(categorical_vars), "\n")
  cat("Total variables for modeling:", length(vars_to_keep) + length(categorical_vars), "\n\n")

  cat("Retained continuous variables:\n")
  cat(paste(vars_to_keep, collapse = "\n"), "\n\n")

  # =============================================================================
  # VERIFY WITH FINAL CORRELATION PLOT
  # =============================================================================

  if (length(vars_to_keep) > 1) {
    cat("=== Generating Final Correlation Heatmap ===\n")
    env_sample_filtered <- env_sample[, vars_to_keep, drop = FALSE]
    cor_matrix_filtered <- cor(env_sample_filtered, use = "complete.obs", method = cor_method)

    plot_correlation_heatmap(cor_matrix_filtered,
                             paste0("Final Correlation Matrix - ", length(vars_to_keep), " variables"),
                             cor_method)

    # Count remaining high correlations
    n_high_cor_remaining <- 0
    for (i in 1:(ncol(cor_matrix_filtered)-1)) {
      for (j in (i+1):ncol(cor_matrix_filtered)) {
        if (abs(cor_matrix_filtered[i, j]) > CONFIG$correlation_threshold) {
          n_high_cor_remaining <- n_high_cor_remaining + 1
        }
      }
    }

    cat("\nRemaining pairs with |r| >", CONFIG$correlation_threshold, ":", n_high_cor_remaining, "\n")

    if (n_high_cor_remaining > 0) {
      cat("WARNING: Some correlations remain above threshold.\n")
    } else {
      cat("SUCCESS: All correlations above", CONFIG$correlation_threshold, "have been eliminated!\n")
    }
  }

} else if (nrow(high_cor_pairs) == 0) {
  cat("No highly correlated variable pairs found (threshold =", CONFIG$correlation_threshold, ")\n")
  vars_to_keep <- continuous_vars
} else {
  cat("Skipping removal (no MaxEnt data)\n")
  vars_to_keep <- continuous_vars
}

# =============================================================================
# VERIFY WITH FINAL CORRELATION PLOT
# =============================================================================

if (length(vars_to_keep) > 1) {
  cat("=== Generating Final Correlation Heatmap ===\n")
  env_sample_filtered <- env_sample[, vars_to_keep, drop = FALSE]
  cor_matrix_filtered <- cor(env_sample_filtered,
                             use = "complete.obs",
                             method = cor_method)

  plot_correlation_heatmap(cor_matrix_filtered,
                           paste0("Final Correlation Matrix - ",
                                  length(vars_to_keep), " variables"),
                           cor_method)

  # Count remaining high correlations
  n_high_cor_remaining <- 0
  for (i in 1:(ncol(cor_matrix_filtered)-1)) {
    for (j in (i+1):ncol(cor_matrix_filtered)) {
      if (abs(cor_matrix_filtered[i, j]) > CONFIG$correlation_threshold) {
        n_high_cor_remaining <- n_high_cor_remaining + 1
      }
    }
  }

  cat("\nRemaining pairs with |r| >", CONFIG$correlation_threshold, ":",
      n_high_cor_remaining, "\n")

  if (n_high_cor_remaining > 0) {
    cat("WARNING: Some correlations remain above threshold.\n")
  } else {
    cat("SUCCESS: All correlations above", CONFIG$correlation_threshold,
        "have been eliminated!\n")
  }

} else if (nrow(high_cor_pairs) == 0) {
  cat("No highly correlated variable pairs found (threshold =",
      CONFIG$correlation_threshold, ")\n")
  vars_to_keep <- continuous_vars
} else {
  cat("Skipping removal (no MaxEnt data)\n")
  vars_to_keep <- continuous_vars
}

# =============================================================================
# VIF ANALYSIS
# =============================================================================

if (exists("vars_to_keep") && length(vars_to_keep) > 1) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("=== VIF ANALYSIS ===\n")
  cat(rep("=", 80), "\n\n", sep = "")

  env_vars_retained <- env_vars[[vars_to_keep]]

  tryCatch({
    vif_result <- vif(env_vars_retained)
    vif_df <- data.frame(
      Variable = vif_result$Variables,
      VIF = vif_result$VIF,
      stringsAsFactors = FALSE
    )
    vif_df <- vif_df[order(vif_df$VIF, decreasing = TRUE), ]

    cat("VIF Results (sorted by VIF):\n\n")
    print(vif_df, row.names = FALSE)

    high_vif <- vif_df[vif_df$VIF > CONFIG$vif_warning_threshold, ]
    if (nrow(high_vif) > 0) {
      cat("\nWARNING: Variables with VIF >", CONFIG$vif_warning_threshold, ":\n")
      print(high_vif, row.names = FALSE)
    } else {
      cat("\nAll VIF values <", CONFIG$vif_warning_threshold,
          ". Multicollinearity well-controlled.\n")
    }
  }, error = function(e) {
    cat("Error calculating VIF:", e$message, "\n")
  })
}

cat("\n", rep("=", 80), "\n", sep = "")
cat("=== ANALYSIS COMPLETE ===\n")
cat(rep("=", 80), "\n", sep = "")

#VIF = 1: No multicollinearity. The predictor is not correlated with other predictors.
#VIF 1 - 3: Generally considered safe; not a cause for concern.
#VIF 3 - 5: Indicates moderate multicollinearity. This level may warrant attention, but is not always problematic depending on the context.
#VIF > 5: High multicollinearity. The predictor's coefficient estimates may be unstable and unreliable. Consider corrective actions such as removing or combining predictors.
#VIF > 10: Serious multicollinearity. The variance of the coefficient is highly inflated, and the model's results are likely to be misleading. Immediate action is recommended.

# =============================================================================
# SAVE FINAL REMOVAL SUMMARY WITH VIF
# =============================================================================

# Build base summary: all continuous variables, their status, and importance
all_continuous_summary <- data.frame(
  Variable = continuous_vars,
  Status = ifelse(continuous_vars %in% vars_to_remove, "REMOVED", "RETAINED"),
  stringsAsFactors = FALSE
)

if (maxent_loaded) {
  all_continuous_summary$Importance <- sapply(all_continuous_summary$Variable, function(v) {
    v_converted <- convert_var_name(v)
    idx <- which(contributions$Variable == v_converted)
    if (length(idx) > 0) contributions$Importance[idx[1]] else NA
  })
} else {
  all_continuous_summary$Importance <- NA_real_
}

# Add pairwise-removal info from removal_log (only for removed variables)
if (exists("removal_log") && nrow(removal_log) > 0) {
  removal_log$Variable <- as.character(removal_log$Variable)
  all_continuous_summary <- merge(
    all_continuous_summary,
    removal_log,
    by = "Variable",
    all.x = TRUE,
    suffixes = c("", "_removal")
  )
} else {
  all_continuous_summary$PairCorrelation <- NA_real_
  all_continuous_summary$N_Pairs_Eliminated <- NA_integer_
  all_continuous_summary$Reason <- NA_character_
}

# Add VIF values (only available for retained variables)
if (exists("vif_df")) {
  all_continuous_summary <- merge(
    all_continuous_summary,
    vif_df,
    by = "Variable",
    all.x = TRUE
  )
} else {
  all_continuous_summary$VIF <- NA_real_
}

# Optional: sort for readability (retained first, then by VIF or importance)
all_continuous_summary <- all_continuous_summary[
  order(all_continuous_summary$Status,
        -ifelse(is.na(all_continuous_summary$Importance), -Inf, all_continuous_summary$Importance)),
]

# Write CSV to disk
output_csv_path <- file.path(PATHS$base_dir, "Final_Removal_by_Correlation_Summary.csv")
write.csv(all_continuous_summary, output_csv_path, row.names = FALSE)

cat("\nFinal removal summary with VIF saved to:\n", output_csv_path, "\n")

# =============================================================================
# CORRELATION CHECK FOR SPECIFIC VARIABLE
# =============================================================================

# --- CONFIGURATION: Specify the variable to analyze ---
target_variable <- "swe_1981-2010_eu"  # CHANGE THIS to your variable of interest
target_variable <- "npp_1981-2010_eu"

cat("\n", rep("=", 80), "\n", sep = "")
cat("=== CORRELATION ANALYSIS FOR SPECIFIC VARIABLE ===\n")
cat(rep("=", 80), "\n\n", sep = "")

# Check if target variable exists
if (!target_variable %in% continuous_vars) {
  cat("ERROR: Variable '", target_variable, "' not found in continuous variables.\n", sep = "")
  cat("\nAvailable continuous variables:\n")
  cat(paste(sort(continuous_vars), collapse = "\n"), "\n")
} else {
  cat("Target variable:", target_variable, "\n")
  cat("Correlation threshold:", CONFIG$correlation_threshold, "\n")
  cat("Correlation method:", cor_method, "\n\n")

  # Get correlations with target variable from original matrix
  if (target_variable %in% rownames(cor_matrix)) {
    target_correlations <- cor_matrix[target_variable, ]
    target_correlations <- target_correlations[names(target_correlations) != target_variable]

    # Find variables above threshold
    high_cor_vars <- target_correlations[abs(target_correlations) > CONFIG$correlation_threshold]
    high_cor_vars <- sort(high_cor_vars, decreasing = TRUE)

    cat("=== Variables Highly Correlated with", target_variable, "===\n")
    cat("(|r| >", CONFIG$correlation_threshold, ")\n\n")

    if (length(high_cor_vars) > 0) {
      correlation_df <- data.frame(
        Variable = names(high_cor_vars),
        Correlation = as.numeric(high_cor_vars),
        Abs_Correlation = abs(high_cor_vars),
        Status = ifelse(names(high_cor_vars) %in% vars_to_remove,
                        "REMOVED", "RETAINED"),
        stringsAsFactors = FALSE
      )

      # Add importance if available
      if (maxent_loaded) {
        correlation_df$Importance <- sapply(correlation_df$Variable, function(v) {
          v_converted <- convert_var_name(v)
          idx <- which(contributions$Variable == v_converted)
          if (length(idx) > 0) {
            sprintf("%.2f%%", contributions$Importance[idx[1]])
          } else {
            "N/A"
          }
        })
      }

      # Sort by absolute correlation
      correlation_df <- correlation_df[order(-correlation_df$Abs_Correlation), ]

      cat("Found", nrow(correlation_df), "variables with |r| >",
          CONFIG$correlation_threshold, "\n\n")

      print(correlation_df[, !names(correlation_df) %in% "Abs_Correlation"],
            row.names = FALSE)

      # Summary statistics
      cat("\n=== Summary ===\n")
      cat("Total correlated variables:", nrow(correlation_df), "\n")
      cat("  - Removed:", sum(correlation_df$Status == "REMOVED"), "\n")
      cat("  - Retained:", sum(correlation_df$Status == "RETAINED"), "\n")
      cat("Correlation range: [",
          sprintf("%.3f", min(correlation_df$Correlation)), ", ",
          sprintf("%.3f", max(correlation_df$Correlation)), "]\n", sep = "")

      # Check if target variable was removed
      if (target_variable %in% vars_to_remove) {
        cat("\nNOTE: Target variable '", target_variable,
            "' was REMOVED during multicollinearity analysis.\n", sep = "")
      } else if (exists("vars_to_keep")) {
        cat("\nNOTE: Target variable '", target_variable,
            "' was RETAINED for modeling.\n", sep = "")
      }

    } else {
      cat("No variables found with |r| >", CONFIG$correlation_threshold, "\n")
      cat("\nHighest correlations with", target_variable, "(top 10):\n\n")

      all_cors <- sort(abs(target_correlations), decreasing = TRUE)[1:min(10, length(target_correlations))]
      top_df <- data.frame(
        Variable = names(all_cors),
        Correlation = target_correlations[names(all_cors)],
        Abs_Correlation = all_cors,
        stringsAsFactors = FALSE
      )
      print(top_df[, c("Variable", "Correlation")], row.names = FALSE)
    }

  } else {
    cat("ERROR: Variable '", target_variable,
        "' not found in correlation matrix.\n", sep = "")
  }
}

cat("\n", rep("=", 80), "\n", sep = "")

