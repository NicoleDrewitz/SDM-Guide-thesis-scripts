# =============================================================================
# Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Stepwise Variable Removal MaxEnt Analysis - VECTORIZED VERSION
# Removes lowest permutation importance variables iteratively based on MaxEnt HTML output
# Created with Claude AI October 2025 by Nicole Drewitz
# Last updated: November 18, 2025

# only 10 000 background points (MaxEnt default) to speed processing.
# Optimized with vectorized operations for faster processing
# Variable update 17 Jan, 2026, RM=0.5 and LPTH
# =============================================================================

# =============================================================================
# INSTALL AND LOAD REQUIRED LIBRARIES
# =============================================================================

required_packages <- c("dismo", "terra", "ggplot2", "dplyr", "sf", "rJava", "xml2", "rvest", "raster", "tictoc", "here")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

 tic()
# =============================================================================
# USER CONFIGURATION - MODIFY THESE SETTINGS
# =============================================================================

FIXED_REG_MULTIPLIER <- 0.5
FIXED_FEATURES <- c("linear", "product", "threshold", "hinge") 
MINIMUM_VARIABLES <- 3
FEATURE_NAME <- paste0(substr(FIXED_FEATURES, 1, 1), collapse = "")
base_dir <- here::here("Output_files")
setwd(base_dir)
species_file_MaxEnt <- here::here("Input_files/V2_cloudberry_presence_MaxEnt.csv")
env_vars_dir <- here::here("Input_files/Environmental_Variables-processed-current_climate")
maxent_jar <- "/Applications/maxent.jar"
output_dir <- file.path(base_dir, "MaxEnt_Results_Stepwise_Variable_removal")
dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)

# =============================================================================
# HELPER FUNCTION: Normalize variable names for matching/removal
# =============================================================================

normalize_varname <- function(x) {
  gsub("[.-]", "", x)
}

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

cat("=== Loading Species Occurrence Data ===\n")
species_data_MaxEnt <- read.csv(species_file_MaxEnt, stringsAsFactors = FALSE)
presence_points <- species_data_MaxEnt[, c("longitude", "latitude")]
presence_points <- na.omit(presence_points)
cat("MaxEnt training points:", nrow(presence_points), "\n")

cat("\n=== Loading Environmental Variables ===\n")
raster_files <- list.files(env_vars_dir, pattern = "\\.asc$", full.names = TRUE, ignore.case = TRUE)
if (length(raster_files) == 0) stop("No .asc files found!")

env_raster <- rast(raster_files)
var_names <- tools::file_path_sans_ext(basename(raster_files))
names(env_raster) <- var_names
cat("Loaded", nlyr(env_raster), "total environmental layers\n")

cat("\n=== Identifying Categorical Variables ===\n")
categorical_prefixes <- c("kg0"
  #, "Water","Forest_type","Peat","Soil_Depth",
  #"Topo_geomflat","Topo_geomhollow","Topo_geompit","Topo_geomslope","Topo_geomvalley"
)
categorical_pattern <- paste0("^(", paste(categorical_prefixes, collapse = "|"), ")")
categorical_vars <- names(env_raster)[grepl(categorical_pattern, names(env_raster))]

if(length(categorical_vars) > 0) {
  cat("Found", length(categorical_vars), "categorical variables:\n")
  cat(paste("  ", categorical_vars, collapse = "\n"), "\n")
} else cat("No categorical variables found\n")

cat("\n=== Identifying Continuous Variables ===\n")
continuous_prefixes <- c(
 "fcf", "bio11", "bio13", "bio19", "soil_waterV33", "soil_phh2o", "soil_soc", "swe", "Topo_northness"
)
continuous_pattern <- paste0("^(", paste(continuous_prefixes, collapse = "|"), ")")
continuous_vars <- names(env_raster)[grepl(continuous_pattern, names(env_raster))]

if(length(continuous_vars) > 0) {
  cat("Found", length(continuous_vars), "continuous variables:\n")
  cat(paste("  ", continuous_vars, collapse = "\n"), "\n")
} else cat("No continuous variables found\n")

cat("\n=== Filtering Environmental Variables ===\n")
# Filter raster stack to selected variables
selected_vars <- c(categorical_vars, continuous_vars)
if(length(selected_vars) == 0) stop("ERROR: No variables identified!")
env_raster <- subset(env_raster, selected_vars)
initial_env_raster <- env_raster
initial_var_names <- names(env_raster)
initial_categorical_vars <- categorical_vars

convert_terra_to_raster <- function(terra_raster) {
  raster_list <- lapply(1:nlyr(terra_raster), function(i) raster::raster(terra_raster[[i]]))
  raster::stack(raster_list)
}

cat("\n=== Data Preparation ===\n")
set.seed(123)
env_stack <- convert_terra_to_raster(env_raster)
n_total <- nrow(presence_points)
training_indices <- sample(1:n_total, round(0.75 * n_total))
training_points <- presence_points[training_indices, ]
testing_points <- presence_points[-training_indices, ]

all_values <- raster::extract(env_stack, training_points)
clean_occurrence_points <- training_points[complete.cases(all_values), ]

test_values <- raster::extract(env_stack, testing_points)
testing_points <- testing_points[complete.cases(test_values), ]

cat("Clean training points:", nrow(clean_occurrence_points), "\n")
cat("Clean testing points:", nrow(testing_points), "\n")

cat("\nInitializing Java for MaxEnt...\n")
tryCatch({
  .jinit(parameters = "-Xmx4g")
  cat("Java initialized successfully with 4GB memory\n")
}, error = function(e) {
  tryCatch({
    .jinit(parameters = "-Xmx2g")
    cat("Java initialized with 2GB memory\n")
  }, error = function(e2) {
    tryCatch({
      .jinit()
      cat("Java initialized with untuned settings\n")
    }, error = function(e3) {
      stop("CRITICAL ERROR: Could not initialize Java: ", conditionMessage(e3))
    })
  })
})

# =============================================================================
# Extended Extract Variable Importance to parse percent contribution as well
# =============================================================================
extract_variable_importance <- function(html_file_path) {
  cat("\n--- Debugging HTML Parsing ---\nParsing file:", html_file_path, "\n")
  tryCatch({
    html_content <- read_html(html_file_path)
    tables <- html_content %>% html_nodes("table")
    cat("Found", length(tables), "tables in HTML\n")
    importance_data <- NULL
    for (i in seq_along(tables)) {
      table_data <- tables[[i]] %>% html_table(fill = TRUE)
      cat("Table", i, "- Dimensions:", nrow(table_data), "rows x", ncol(table_data), "cols\n")
      if (ncol(table_data) > 0) cat("  Column names:", paste(names(table_data), collapse = ", "), "\n")
      
      if (ncol(table_data) >= 3) {
        col_names <- tolower(names(table_data))
        var_col_idx <- grep("variable", col_names, ignore.case = TRUE)
        perm_col_idx <- grep("permutation.*importance", col_names, ignore.case = TRUE)
        contrib_col_idx <- grep("percent.*contribution", col_names, ignore.case = TRUE)
        
        if (length(var_col_idx) == 0) var_col_idx <- which(col_names == "variable")
        if (length(perm_col_idx) == 0) perm_col_idx <- grep("permutation|importance", col_names, ignore.case = TRUE)
        if (length(contrib_col_idx) == 0) contrib_col_idx <- grep("percent|contribution", col_names, ignore.case = TRUE)
        
        if (length(var_col_idx) > 0 && length(perm_col_idx) > 0 && length(contrib_col_idx) > 0) {
          var_col_idx <- var_col_idx[1]
          perm_col_idx <- perm_col_idx[1]
          contrib_col_idx <- contrib_col_idx[1]
          
          cat("  Using columns:", var_col_idx, "(variable),", perm_col_idx, "(permutation importance), and", contrib_col_idx, "(percent contribution)\n")
          importance_data <- data.frame(
            variable = as.character(table_data[[var_col_idx]]),
            permutation_importance = as.numeric(as.character(table_data[[perm_col_idx]])),
            percent_contribution = as.numeric(as.character(table_data[[contrib_col_idx]])),
            stringsAsFactors = FALSE
          )
          valid_rows <- !is.na(importance_data$variable) &
            importance_data$variable != "" &
            !is.na(importance_data$permutation_importance) &
            !is.na(importance_data$percent_contribution) &
            importance_data$permutation_importance >= 0
          importance_data <- importance_data[valid_rows, ]
          cat("  Valid rows after filtering:", nrow(importance_data), "\n")
          if (nrow(importance_data) > 0) {
            importance_data$table_order <- seq_len(nrow(importance_data))
            cat("  SUCCESS: Found importance table with", nrow(importance_data), "variables\n")
            break
          }
        }
      }
    }
    if (is.null(importance_data) || nrow(importance_data) == 0) {
      cat("\n!!! Table parsing failed. Attempting text-based parsing... !!!\n")
      html_text <- html_content %>% html_text2()
      lines <- strsplit(html_text, "\n")[[1]]
      importance_lines <- grep("^[A-Za-z_][A-Za-z0-9_.-]*\\s+[0-9.]+", lines, value = TRUE)
      cat("Found", length(importance_lines), "potential importance lines\n")
      if (length(importance_lines) > 0) {
        parsed_data <- lapply(importance_lines, function(line) {
          parts <- strsplit(trimws(line), "\\s+")[[1]]
          if (length(parts) >= 2) {
            var_name <- parts[1]
            numeric_parts <- parts[sapply(parts, function(x) !is.na(as.numeric(x)))]
            if (length(numeric_parts) > 0) {
              importance_val <- as.numeric(numeric_parts[length(numeric_parts)])
              return(data.frame(variable = var_name, 
                                permutation_importance = importance_val,
                                percent_contribution = NA,
                                stringsAsFactors = FALSE))
            }
          }
          return(NULL)
        })
        parsed_data <- do.call(rbind, parsed_data[!sapply(parsed_data, is.null)])
        if (!is.null(parsed_data) && nrow(parsed_data) > 0) {
          parsed_data$table_order <- seq_len(nrow(parsed_data))
          importance_data <- parsed_data
          cat("Text parsing successful:", nrow(importance_data), "variables found\n")
        }
      }
    }
    cat("--- End HTML Parsing Debug ---\n\n")
    return(importance_data)
  }, error = function(e) {
    cat("!!! ERROR in extract_variable_importance:", e$message, "!!!\n")
    cat("Call stack:\n")
    print(sys.calls())
    return(NULL)
  })
}

# =============================================================================
# GENERATE BACKGROUND POINTS (ONCE)
# =============================================================================

cat("\n=== Generating Background Points ===\n")
initial_env_stack <- convert_terra_to_raster(initial_env_raster)
set.seed(123)
background_points <- randomPoints(
  initial_env_stack,
  n = 20000,
  p = clean_occurrence_points,
  excludep = TRUE, # cells with presence points will not be selected for background points
  prob = TRUE # more background sampling in cells closer to presence cells
)
bg_env_values <- raster::extract(initial_env_stack, background_points)
background_points_clean <- background_points[complete.cases(bg_env_values), ]
cat("Generated", nrow(background_points_clean), "clean background points\n")
cat("These same background points will be used for all iterations\n")

# =============================================================================
# STEPWISE VARIABLE REMOVAL LOOP
# =============================================================================

current_variables <- initial_var_names
iteration <- 1
removal_log <- data.frame(
  iteration = integer(),
  variables_remaining = integer(),
  variables_remaining_listed = character(),
  runtime = numeric(),
  removed_variable = character(),
  removed_var_importance = numeric(),
  removed_var_contribution = numeric(),
  removed_by = character(),
  tied_variables = character(),
  variables_used = character(),
  train_auc = numeric(),
  test_auc = numeric(),
  stringsAsFactors = FALSE
)
MAX_ITERATIONS <- length(current_variables) - MINIMUM_VARIABLES

cat("\n=== Starting Stepwise Variable Removal ===\n")
cat("Starting with", length(current_variables), "filtered variables\n")
cat("Minimum variables to retain:", MINIMUM_VARIABLES, "\n")
cat("Maximum possible iterations:", MAX_ITERATIONS, "\n\n")

while (length(current_variables) > MINIMUM_VARIABLES && iteration <= MAX_ITERATIONS) {
  
  cat("\n========================================\n")
  cat("--- Iteration", iteration, "of max", MAX_ITERATIONS, "---\n")
  cat("========================================\n")
  cat("Variables remaining:", length(current_variables), "\n")
  
  iter_start_time <- Sys.time()
  iter_dir <- file.path(output_dir, paste0("RUN", sprintf("%02d", iteration)))
  dir.create(iter_dir, recursive = TRUE, showWarnings = FALSE)
  
  current_indices <- which(initial_var_names %in% current_variables)
  current_env_raster <- initial_env_raster[[current_indices]]
  current_env_stack <- convert_terra_to_raster(current_env_raster)
  current_categorical_vars <- intersect(initial_categorical_vars, current_variables)
  current_factors_vec <- NULL
  if (length(current_categorical_vars) > 0) {
    current_factors_vec <- which(names(current_env_stack) %in% current_categorical_vars)
  }
  
  cat("Using", nrow(background_points_clean), "pre-generated background points\n")
  
  maxent_args <- c(
    paste0("betamultiplier=", FIXED_REG_MULTIPLIER),
    "randomtestpoints=25",
    "maximumbackground=10000",
    "removeduplicates=true",
    "allowpartialdata=true",
    "autofeature=false",
    "outputformat=logistic",
    "askoverwrite=false",
    "skipifexists=false",
    "responsecurves=true",
    "pictures=false",
    "plots=false",
    "writeclampgrid=false",
    "writemess=false",
    "randomseed=true",
    "addsamplestobackground=true",
    "addallsamplestobackground=false",
    "fadebyclamping=false",
    "extrapolate=true",
    "doclamp=true",
    "outputgrids=true",
    "maximumiterations=500",
    "convergencethreshold=0.00001",
    "threads=1"
  )
  
  all_features <- c("linear", "quadratic", "product", "threshold", "hinge")
  feature_args <- sapply(all_features, function(feat) {
    paste0(feat, "=", tolower(as.character(feat %in% FIXED_FEATURES)))
  }, USE.NAMES = FALSE)
  maxent_args <- c(maxent_args, feature_args)
  
  cat("Running MaxEnt with", length(current_variables), "variables...\n")
  
  tryCatch({
    
    model <- maxent(
      x = current_env_stack,
      p = clean_occurrence_points,
      path = iter_dir,
      args = maxent_args,
      factors = current_factors_vec
    )
    
    if (!is.null(model)) {
      cat("Creating habitat suitability prediction...\n")
      prediction <- predict(model, current_env_stack)
      asc_output <- file.path(iter_dir, paste0("habitat_suitability_", FEATURE_NAME, "_",
                                               FIXED_REG_MULTIPLIER, ".asc"))
      raster::writeRaster(prediction, asc_output, format = "ascii", overwrite = TRUE)
      cat("Prediction files created successfully\n")
    }
    
    train_auc <- NA
    test_auc <- NA
    results_files <- list.files(iter_dir, pattern = "maxentResults\\.csv$", full.names = TRUE)
    if (length(results_files) > 0) {
      results_table <- read.csv(results_files[1], stringsAsFactors = FALSE)
      if (nrow(results_table) > 0) {
        train_cols <- grep("Training.AUC|AUC.training", names(results_table), ignore.case = TRUE, value = TRUE)
        test_cols <- grep("Test.AUC|AUC.test", names(results_table), ignore.case = TRUE, value = TRUE)
        if (length(train_cols) > 0) train_auc <- as.numeric(results_table[[train_cols[1]]][1])
        if (length(test_cols) > 0) test_auc <- as.numeric(results_table[[test_cols[1]]][1])
      }
    }
    
    cat("Model results - Train AUC:", round(train_auc, 3), "Test AUC:", round(test_auc, 3), "\n")
    
    html_files <- list.files(iter_dir, pattern = "\\.html$", full.names = TRUE)
    
    if (length(html_files) > 0) {
      cat("\nParsing variable permutation importance from HTML output...\n")
      importance <- extract_variable_importance(html_files[1])
      
      if (!is.null(importance) && nrow(importance) > 0) {
        cat("\nVariable permutation importance (sorted):\n")
        importance_sorted <- importance[order(importance$permutation_importance, -importance$table_order), ]
        print(importance_sorted)
        
        if (length(current_variables) > MINIMUM_VARIABLES) {
          min_importance_value <- min(importance$permutation_importance)
          min_importance_rows <- importance[importance$permutation_importance == min_importance_value, ]
          
          removed_by <- "permutation_importance"
          tied_vars <- ""
          
          if (nrow(min_importance_rows) > 1) {
            tied_vars <- paste(min_importance_rows$variable, collapse = "; ")
            # Fix: choose highest table_order among ties (max = least important)
            lowest_importance_var <- min_importance_rows$variable[which.max(min_importance_rows$table_order)]
            removed_by <- "table_order"
            cat("\nMultiple variables tied with importance", round(min_importance_value, 2), "- using table order\n")
            cat("Tied variables:", tied_vars, "\n")
          } else {
            lowest_importance_var <- min_importance_rows$variable[1]
          }
          lowest_importance_value <- min_importance_value
          lowest_importance_contribution <- importance$percent_contribution[importance$variable == lowest_importance_var]
          
          # Normalize variable names for matching/removal
          lowest_importance_var_norm <- normalize_varname(lowest_importance_var)
          current_vars_norm <- sapply(current_variables, normalize_varname)
          
          # Find exact match in normalized names
          match_idx <- which(current_vars_norm == lowest_importance_var_norm)
          
          cat("\n*** REMOVING VARIABLE ***\n")
          cat("Variable to remove:", lowest_importance_var, "\n")
          cat("Permutation importance:", round(lowest_importance_value, 2), "\n")
          cat("Percent contribution:", round(lowest_importance_contribution, 2), "\n")
          cat("Removed by:", removed_by, "\n")
          
          iter_end_time <- Sys.time()
          iter_runtime <- as.numeric(difftime(iter_end_time, iter_start_time, units = "mins"))
          previous_count <- length(current_variables)
          
          if(length(match_idx) == 1) {
            var_to_remove <- current_variables[match_idx]
            current_variables <- current_variables[-match_idx]
            new_count <- length(current_variables)
            
            cat("\n*** VARIABLE REMOVAL CHECK ***\n")
            cat("Variables before removal:", previous_count, "\n")
            cat("Variables after removal:", new_count, "\n")
            cat("SUCCESS: Variable successfully removed\n")
            
            vars_remaining_list <- paste(current_variables, collapse = ", ")
            vars_used_str <- paste(current_variables, collapse = ", ")
            
            removal_log <- rbind(removal_log, data.frame(
              iteration = iteration,
              variables_remaining = new_count,
              variables_remaining_listed = vars_remaining_list,
              runtime = round(iter_runtime, 2),
              removed_variable = var_to_remove,
              removed_var_importance = lowest_importance_value,
              removed_var_contribution = lowest_importance_contribution,
              removed_by = removed_by,
              tied_variables = tied_vars,
              variables_used = vars_used_str,
              train_auc = train_auc,
              test_auc = test_auc,
              stringsAsFactors = FALSE
            ))
          } else {
            cat("\n!!! WARNING: Variable was NOT removed! Name mismatch detected. !!!\n")
            cat("Attempted to remove normalized variable:", lowest_importance_var_norm, "\n")
            cat("No exact match found after normalization in current variables\n")
            cat("Searching for approximate matches...\n")
            similar <- grep(lowest_importance_var_norm, current_vars_norm, value = TRUE)
            if(length(similar) > 0) {
              cat("Possible matches found:\n")
              cat(paste("  ", similar, collapse = "\n"), "\n")
            }
            cat("\n!!! STOPPING TO PREVENT INFINITE LOOP !!!\n")
            break
          }
        } else {
          cat("\n*** Reached minimum number of variables. Stopping. ***\n")
          break
        }
      } else {
        cat("\n!!! WARNING: Could not parse variable permutation importance from HTML output !!!\n")
        cat("Available HTML files:", html_files, "\n")
        cat("!!! STOPPING !!!\n")
        break
      }
    } else {
      cat("\n!!! ERROR: No HTML output files found !!!\n")
      cat("!!! STOPPING !!!\n")
      break
    }
    
  }, error = function(e) {
    cat("\n!!! ERROR in iteration", iteration, ": !!!\n")
    cat(conditionMessage(e), "\n")
    cat("!!! STOPPING !!!\n")
    break
  })
  
  iteration <- iteration + 1
  gc()
}

cat("\n========================================\n")
cat("=== STEPWISE REMOVAL COMPLETE ===\n")
cat("========================================\n")
cat("Total iterations completed:", iteration - 1, "\n")
cat("Final number of variables:", length(current_variables), "\n")
cat("Final variables remaining:\n")
cat(paste("  ", current_variables, collapse = "\n"), "\n\n")

log_file <- file.path(output_dir, "variable_removal_log.csv")
write.csv(removal_log, log_file, row.names = FALSE)
cat("Removal log saved to:", log_file, "\n")

cat("\nVariable Removal Summary:\n")
if (nrow(removal_log) > 0) {
  summary_cols <- c("iteration", "variables_remaining", "removed_variable",
                    "removed_by", "removed_var_importance", "removed_var_contribution",
                    "variables_used", "train_auc", "test_auc")
  print(removal_log[, summary_cols])
}

cat("\n=== Analysis Complete ===\n")

# =============================================================================
# FINAL RESULTS AND SUMMARY
# =============================================================================

cat("\n=== Stepwise Variable Removal Complete ===\n")
cat("Final variables remaining:", paste(current_variables, collapse = ", "), "\n")
cat("Total iterations completed:", nrow(removal_log), "\n")

# Create summary visualization
if (nrow(removal_log) > 0) {
  cat("\nCreating summary plots...\n")
  
  # Plot AUC trends
  p1 <- ggplot(removal_log, aes(x = variables_remaining)) +
    geom_line(aes(y = train_auc, color = "Training AUC"), linewidth = 1) +
    geom_line(aes(y = test_auc, color = "Test AUC"), linewidth = 1) +
    geom_point(aes(y = train_auc, color = "Training AUC"), size = 2) +
    geom_point(aes(y = test_auc, color = "Test AUC"), size = 2) +
    scale_color_manual(values = c("Training AUC" = "blue", "Test AUC" = "red")) +
    labs(title = "Model Performance vs Number of Variables",
         x = "Number of Variables Remaining",
         y = "AUC (from MaxEnt)",
         color = "Metric") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, "maxent_auc_trends.png"), p1, width = 8, height = 6, dpi = 300)
  
  cat("Summary plots saved to:", output_dir, "\n")
}

toc()
cat("\nAnalysis complete!\n")
