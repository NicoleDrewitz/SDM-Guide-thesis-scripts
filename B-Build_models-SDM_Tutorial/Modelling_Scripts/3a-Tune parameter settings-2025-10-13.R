# =============================================================================
# Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Iterate through parameter tuning settings

# Compatible with MaxEnt version 3.4.4 - Includes ASC file generation
# dismo package actually runs version 3.4.3
# only 10 000 background points (maxent default) to speed tuning
# Terminal command to keep computer on for 2 hours: caffeinate -i -t 3600
# Author: Nicole Drewitz with Claude AI
# Last updated: October 13, 2025
# variable set updated 16 Jan, 2026 - 10 variables
# =============================================================================

# =============================================================================
# INSTALL AND LOAD REQUIRED LIBRARIES
# =============================================================================

required_packages <- c("dismo", "terra", "ggplot2", "dplyr", "sf", "gridExtra", "rJava", "raster", "tictoc", "here")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

tic() # time script runtime

# =============================================================================
# HELPER FUNCTION: Convert terra to raster for MaxEnt
# =============================================================================

convert_terra_to_raster <- function(terra_raster) {
  " Convert terra SpatRaster to raster stack for use with dismo::maxent"
  raster_list <- list()
  for (i in 1:nlyr(terra_raster)) {
    raster_list[[i]] <- raster::raster(terra_raster[[i]])
  }
  return(raster::stack(raster_list))
}

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# File paths
base_dir <- here::here("Output_files")
setwd(base_dir)

# File paths
species_file_MaxEnt <- here::here("Input_files/V2_cloudberry_presence_MaxEnt.csv")
env_vars_dir <- here::here("Input_files/Environmental_Variables-processed-current_climate")
maxent_jar <- "/Applications/maxent.jar"

output_dir <- file.path(base_dir, "MaxEnt_Results_parameter_tuning")
dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

cat("=== Loading Species Occurrence Data ===\n")

# Load and prepare presence data for MaxEnt
species_data_MaxEnt <- read.csv(species_file_MaxEnt, stringsAsFactors = FALSE)
presence_points <- species_data_MaxEnt[, c("longitude", "latitude")]
presence_points <- na.omit(presence_points)

cat("MaxEnt training points:", nrow(presence_points), "\n")

# Load environmental variables
cat("\n=== Loading Environmental Variables ===\n")
raster_files <- list.files(env_vars_dir, pattern = "\\.asc$", full.names = TRUE, ignore.case = TRUE)

if (length(raster_files) == 0) {
  stop("No .asc files found in environmental variables directory!")
}

env_raster <- rast(raster_files)
var_names <- tools::file_path_sans_ext(basename(raster_files))
names(env_raster) <- var_names

cat("Loaded", nlyr(env_raster), "total environmental layers\n")

# =============================================================================
# IDENTIFY CATEGORICAL VARIABLES
# =============================================================================

cat("\n=== Identifying Categorical Variables ===\n")

categorical_prefixes <- c("kg0") #add , within brackets to list more variables

categorical_vars <- character(0)
for (prefix in categorical_prefixes) {
  matching_vars <- names(env_raster)[startsWith(names(env_raster), prefix)]
  if (length(matching_vars) > 0) {
    categorical_vars <- c(categorical_vars, matching_vars)
  }
}

if (length(categorical_vars) > 0) {
  cat("Found", length(categorical_vars), "categorical variables:\n")
  for (var in categorical_vars) {
    cat("  ", var, "\n")
  }
} else {
  cat("No categorical variables found\n")
}

# =============================================================================
# IDENTIFY CONTINUOUS VARIABLES
# =============================================================================

cat("\n=== Identifying Continuous Variables ===\n")

# Define continuous variable prefixes to match the correlation subset
continuous_prefixes <- c("fcf", "bio11", "bio13", "bio19", "soil_waterV33", "soil_phh2o", "soil_soc", "swe", "Topo_northness")

# Function to find variable names starting with any prefix in the list
find_vars_by_prefix <- function(prefixes, var_names) {
  matched_vars <- character(0)
  for (prefix in prefixes) {
    matching <- var_names[startsWith(var_names, prefix)]
    if (length(matching) > 0) {
      matched_vars <- c(matched_vars, matching)
    }
  }
  unique(matched_vars)
}

continuous_vars <- find_vars_by_prefix(continuous_prefixes, names(env_raster))

if (length(continuous_vars) > 0) {
  cat("Found", length(continuous_vars), "continuous variables:\n")
  for (var in continuous_vars) {
    cat("  ", var, "\n")
  }
} else {
  cat("No continuous variables found\n")
}

# =============================================================================
# FILTER RASTER STACK TO SELECTED VARIABLES ONLY
# =============================================================================

cat("\n=== Filtering Environmental Variables ===\n")

# Combine identified categorical and continuous variables
selected_vars <- c(categorical_vars, continuous_vars)

if (length(selected_vars) == 0) {
  stop("ERROR: No variables were identified as categorical or continuous!")
}

cat("Total selected variables:", length(selected_vars), "\n")
cat("  Categorical:", length(categorical_vars), "\n")
cat("  Continuous:", length(continuous_vars), "\n")

# Subset the raster stack to keep ONLY these variables
env_raster <- subset(env_raster, selected_vars)

cat("\nFiltered to", nlyr(env_raster), "environmental layers for modeling\n")
cat("Variables included:\n")
for (var in names(env_raster)) {
  cat("  ", var, "\n")
}

# Store initial variables for iteration tracking
initial_env_raster <- env_raster
initial_var_names <- names(env_raster)
initial_categorical_vars <- categorical_vars

# =============================================================================
# DATA SPLITTING AND CLEANING
# =============================================================================

cat("\n=== Data Preparation ===\n")
set.seed(123)

# Split presence data for training/testing
n_total <- nrow(presence_points)
training_indices <- sample(1:n_total, round(0.75 * n_total))
training_points <- presence_points[training_indices, ]
testing_points <- presence_points[-training_indices, ]

# Convert terra to raster for extraction (needed for compatibility)
env_stack <- convert_terra_to_raster(env_raster)

# Clean data - remove points with NA values in any environmental layer
all_values <- raster::extract(env_stack, training_points)
complete_cases <- complete.cases(all_values)
training_points <- training_points[complete_cases, ]

test_values <- raster::extract(env_stack, testing_points)
test_complete_cases <- complete.cases(test_values)
testing_points <- testing_points[test_complete_cases, ]

cat("Clean training points:", nrow(training_points), "\n")
cat("Clean testing points:", nrow(testing_points), "\n")

# =============================================================================
# IMPROVED MAXENT PARAMETER TUNING
# =============================================================================

cat("\n=== Improved MaxEnt Parameter Tuning ===\n")

# Check MaxEnt availability first
if (!file.exists(maxent_jar)) {
  cat("ERROR: MaxEnt jar file not found at:", maxent_jar, "\n")
  cat("Trying to find MaxEnt in dismo package...\n")
  
  alternative_paths <- c(
    system.file("java", "maxent.jar", package = "dismo"),
    file.path(system.file(package = "dismo"), "java", "maxent.jar")
  )
  
  maxent_found <- FALSE
  for (path in alternative_paths) {
    if (file.exists(path)) {
      maxent_jar <- path
      maxent_found <- TRUE
      cat("Found MaxEnt at:", path, "\n")
      break
    }
  }
  
  if (!maxent_found) {
    stop("CRITICAL ERROR: MaxEnt jar file not found! Please install MaxEnt and update the path.")
  }
}

# Initialize Java with better error handling
cat("Initializing Java for MaxEnt...\n")
tryCatch({
  .jinit(parameters = "-Xmx4g")
  cat("Java initialized successfully with 4GB memory\n")
}, error = function(e) {
  cat("Java initialization with 4GB failed, trying 2GB:", conditionMessage(e), "\n")
  tryCatch({
    .jinit(parameters = "-Xmx2g")
    cat("Java initialized with 2GB memory\n")
  }, error = function(e2) {
    cat("Java initialization with 2GB failed, trying default:", conditionMessage(e2), "\n")
    tryCatch({
      .jinit()
      cat("Java initialized with untuned settings\n")
    }, error = function(e3) {
      stop("CRITICAL ERROR: Could not initialize Java: ", conditionMessage(e3))
    })
  })
})

# More robust data preprocessing using terra package
cat("Preprocessing occurrence data...\n")

# Convert points to SpatVector
training_vect <- vect(training_points, geom = c("longitude", "latitude"))

# Get cell numbers for each point
cell_numbers <- cells(env_raster[[1]], training_vect)[, "cell"]

# Remove duplicate points in same cells
unique_cells <- !duplicated(cell_numbers)
training_points_clean <- training_points[unique_cells, ]

cat("Removed", sum(!unique_cells), "duplicate points within same pixels\n")

# Ensure no NA values with more robust checking
env_values <- raster::extract(env_stack, training_points_clean)
#env_values <- extract(env_raster, training_points_clean)
complete_rows <- complete.cases(env_values)
training_points_final <- training_points_clean[complete_rows, ]

cat("Final clean training points:", nrow(training_points_final), "\n")

# Generate background points more systematically
cat("Generating background points...\n")
set.seed(123)  # Ensure reproducibility
background_points <- dismo::randomPoints(
  env_stack,  # Use the converted raster stack
  n = 20000, #generate extra since many will be removed as NAs
  p = training_points_final,
  excludep = TRUE,  # Excludes pixels with presence points
  prob = TRUE       # Use probability-based sampling
)

# Remove background points with NA values
bg_env_values <- raster::extract(env_stack, background_points)
bg_complete <- complete.cases(bg_env_values)
background_points_clean <- background_points[bg_complete, ]

cat("Clean background points:", nrow(background_points_clean), "\n")

# Create factors vector for categorical variables
if (length(categorical_vars) > 0) {
  factors_vec <- which(names(env_stack) %in% categorical_vars)
} else {
  factors_vec <- NULL
}

# Enhanced parameter combinations
reg_mult_values <- seq(0.5, 5, by = 0.5)
feature_combinations <- list(
  "L" = c("linear"),
  "Q" = c("quadratic"),
  "P" = c("product"),
  "T" = c("threshold"),
  "H" = c("hinge"),
  "LQ" = c("linear", "quadratic"),
  "LP" = c("linear", "product"),
  "LT" = c("linear", "threshold"),
  "LH" = c("linear", "hinge"),
  "QP" = c("quadratic", "product"),
  "QT" = c("quadratic", "threshold"),
  "QH" = c("quadratic", "hinge"),
  "PT" = c("product", "threshold"),
  "PH" = c("product", "hinge"),
  "TH" = c("threshold", "hinge"),
  "LQP" = c("linear", "quadratic", "product"),
  "LQT" = c("linear", "quadratic", "threshold"),
  "LQH" = c("linear", "quadratic", "hinge"),
  "LPT" = c("linear", "product", "threshold"),
  "LPH" = c("linear", "product", "hinge"),
  "LTH" = c("linear", "threshold", "hinge"),
  "QPT" = c("quadratic", "product", "threshold"),
  "QPH" = c("quadratic", "product", "hinge"),
  "QTH" = c("quadratic", "threshold", "hinge"),
  "PTH" = c("product", "threshold", "hinge"),
  "LQPT" = c("linear", "quadratic", "product", "threshold"),
  "LQPH" = c("linear", "quadratic", "product", "hinge"),
  "LQTH" = c("linear", "quadratic", "threshold", "hinge"),
  "LPTH" = c("linear", "product", "threshold", "hinge"),
  "QPTH" = c("quadratic", "product", "threshold", "hinge"),
  "LQPTH" = c("linear", "quadratic", "product", "threshold", "hinge")
)

tuning_results <- data.frame(
  reg_mult = numeric(),
  features = character(),
  train_auc = numeric(),
  test_auc = numeric(),
 # n_parameters = numeric(),
  stringsAsFactors = FALSE
)

# Run enhanced parameter tuning
total_runs <- length(reg_mult_values) * length(feature_combinations)
run_count <- 0

for (reg_mult in reg_mult_values) {
  for (feat_name in names(feature_combinations)) {
    run_count <- run_count + 1
    features <- feature_combinations[[feat_name]]
    
    cat("Run", run_count, "of", total_runs, "- RegMult:", reg_mult, "Features:", feat_name, "\n")
    
    run_dir <- file.path(output_dir, paste0("regmult_", gsub("\\.", "_", reg_mult), "_", feat_name))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    
    tryCatch({
      # Comprehensive MaxEnt arguments matching GUI defaults
      maxent_args <- c(
        paste0("betamultiplier=", reg_mult),
        "randomtestpoints=25",
        "maximumbackground=10000",
        "removeduplicates=true",
        "allowpartialdata=true",
        "autofeature=false",
        "outputformat=logistic",
        "askoverwrite=false",
        "skipifexists=false",
        "responsecurves=false",
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
      
      # Explicit feature control
      all_features <- c("linear", "quadratic", "product", "threshold", "hinge")
      for (feat in all_features) {
        if (feat %in% features) {
          maxent_args <- c(maxent_args, paste0(feat, "=true"))
        } else {
          maxent_args <- c(maxent_args, paste0(feat, "=false"))
        }
      }
      
      # Run MaxEnt with explicit background points and categorical variables
      model <- maxent(
        x = env_stack,  # Use converted raster stack
        p = training_points_final,
        a = background_points_clean,
        path = run_dir,
        args = maxent_args,
        factors = factors_vec
      )
      
      # Create habitat suitability prediction
      if (!is.null(model)) {
        tryCatch({
          cat("  Creating habitat suitability prediction...\n")
          
          # Make prediction using raster stack
          prediction <- predict(model, env_stack)
          
          # Save as ASC file
          asc_output <- file.path(run_dir, paste0("habitat_suitability_", feat_name, "_", reg_mult, ".asc"))
          raster::writeRaster(prediction, asc_output, format = "ascii", overwrite = TRUE)
          
          cat("    SUCCESS: Prediction files created\n")
        }, error = function(pred_err) {
          cat("    Prediction failed:", conditionMessage(pred_err), "\n")
        })
      }
      
      # Extract comprehensive results
      train_auc <- NA
      test_auc <- NA
      #n_parameters <- NA
      
      # Look for results files
      results_files <- list.files(run_dir, pattern = "maxentResults\\.csv$", full.names = TRUE)
      
      if (length(results_files) > 0) {
        results_table <- read.csv(results_files[1], stringsAsFactors = FALSE)
        
        if (nrow(results_table) > 0) {
          # Extract AUC values
          train_cols <- grep("Training.AUC|AUC.training", names(results_table), ignore.case = TRUE, value = TRUE)
          test_cols <- grep("Test.AUC|AUC.test", names(results_table), ignore.case = TRUE, value = TRUE)
          #param_cols <- grep("X.Parameters|parameters", names(results_table), ignore.case = TRUE, value = TRUE)
          
          if (length(train_cols) > 0) train_auc <- as.numeric(results_table[[train_cols[1]]][1])
          if (length(test_cols) > 0) test_auc <- as.numeric(results_table[[test_cols[1]]][1])
          #if (length(param_cols) > 0) n_parameters <- as.numeric(results_table[[param_cols[1]]][1])
        }
      }
      
      # Store results
      tuning_results <- rbind(tuning_results, data.frame(
        reg_mult = reg_mult,
        features = feat_name,
        train_auc = train_auc,
        test_auc = test_auc,
        #n_parameters = n_parameters,
        stringsAsFactors = FALSE
      ))
      
      cat("  Results - Train AUC:", round(train_auc, 3),
          "Test AUC:", round(test_auc, 3),
          #"Parameters:", n_parameters,
          "\n")
      
    }, error = function(e) {
      cat("  Error in run:", conditionMessage(e), "\n")
    })
    
    # Clean up memory
    gc()
  }
}

# Save tuning results
write.csv(tuning_results, file.path(output_dir, "parameter_tuning_results.csv"), row.names = FALSE)

toc() # end script stopwatch
cat("\n=== Parameter Tuning Complete ===\n")
cat("Results saved to:", file.path(output_dir, "parameter_tuning_results.csv"), "\n")
