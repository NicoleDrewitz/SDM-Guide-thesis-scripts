# =============================================================================
# Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# Single Run MaxEnt Model for Analysis
# Created with Claude AI by Nicole Drewitz on 14 October, 2025
# Last updated structurally: 24 November, 2025
# Last updated with variables: 16 January, 2026

# ==============================================================================
# produce models manually in MaxEnt to have access to the latest software version (3.4.4)
#   this code defaults to version 3.4.3
# max iterations increased to 5000

# ==============================================================================
# USER GUIDE
# ==============================================================================
# Downloads required
#   https://biodiversityinformatics.amnh.org/open_source/maxent/
#   https://www.java.com/en/download/manual.jsp
#   Environmental variables (ASC)
#   Species presence locations (CSV)
#   Species absence locations (if selecting best tuning parameters from presence/absence AICc ranking)
# Used to select best tuning parameters
#   For stepwise variable removal using the selected tuning parameters afterwards,
#       it may be better to create models manually in MaxEnt to have access to the latest software version

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

tic() # start timer

# =============================================================================
# USER CONFIGURATION - MODIFY THESE SETTINGS
# =============================================================================
# response curves can be set to false for first run

# Set your preferred parameters (use same for all runs)
FIXED_REG_MULTIPLIER <- 1  # Set your preferred regularization multiplier
FIXED_FEATURES <- c("linear", "quadratic", "product", "hinge")  # untuned settings with RegMult=1

# Feature combination name for output folders
FEATURE_NAME <- paste0(substr(FIXED_FEATURES, 1, 1), collapse = "")  # e.g., "LQH"

# OPTIONAL: Variables to exclude from analysis
# Examples: c("bio1", "bio2") or c() for no exclusions
# File prefixes can be used. bio1 must "bio1_" or else bio11, bio12, etc are also excluded.

#EXCLUDE_VARIABLES <- c("soil_silt") # silt+clay+sand =100% so all are not needed, likely to be removed anyway at such coarse scale

# Variables to remove after high correlation filtering
# Bioclim variable starting with CHELSA_ are from before they were redownloaded from the new website
#EXCLUDE_VARIABLES <- c("soil_silt", "CHELSA_gdgfgd0", "CHELSA_gdgfgd5", "CHELSA_bio6", "CHELSA_ngd5", # silt+sand+clay= 100% so not all are necessary 
#                        "CHELSA_bio7", "CHELSA_bio14", "CHELSA_bio1_", "CHELSA_gddlgd0", "NDVI_mean", "CHELSA_bio3", # bio3 manually replaced with fcf
#                        "CHELSA_bio18", "CHELSA_bio11", "CHELSA_gdd5", "CHELSA_bio8", "soil_ocs", "CHELSA_bio9", # swe not removed
#                        "CHELSA_bio5", "CHELSA_lgd", "CHELSA_scd", "CHELSA_bio17", "CHELSA_bio12", "CHELSA_gdgfgd10",
#                        "CHELSA_gst", "CHELSA_bio10", "CHELSA_bio16", "CHELSA_fgd", "CHELSA_gdd0", "CHELSA_gdd10",
#                        "CHELSA_gddlgd10", "CHELSA_gddlgd5", "CHELSA_gsl", "CHELSA_ngd0", "CHELSA_ngd10", "soil_bdod", "Topo_slope")

# Variables for top 15 (predicted) best variables
      # Created after viewing PCA and response curve from above reduced variable model.
#EXCLUDE_VARIABLES <- c("soil_silt", "CHELSA_gdgfgd0", "CHELSA_gdgfgd5", "CHELSA_bio6", "CHELSA_ngd5", # silt+sand+clay= 100% so not all are necessary 
 #                      "CHELSA_bio7", "CHELSA_bio14", "CHELSA_bio1_", "CHELSA_gddlgd0", "NDVI_mean", "CHELSA_bio3", # bio3 manually replaced with fcf
  #                     "CHELSA_bio18", "CHELSA_bio11", "CHELSA_gdd5", "CHELSA_bio8", "soil_ocs", "CHELSA_bio9", # swe not removed
   #                    "CHELSA_bio5", "CHELSA_lgd", "CHELSA_scd", "CHELSA_bio17", "CHELSA_bio12", "CHELSA_gdgfgd10",
    #                   "CHELSA_gst", "CHELSA_bio10", "CHELSA_bio16", "CHELSA_fgd", "CHELSA_gdd0", "CHELSA_gdd10",
     #                  "CHELSA_gddlgd10", "CHELSA_gddlgd5", "CHELSA_gsl", "CHELSA_ngd0", "CHELSA_ngd10", "soil_bdod", "Topo_slope", # after this row are the new removed variables (for 15 variable model)
      #                 "CHELSA_bio13", "soil_cfvo", "Topo_elevation", "CHELSA_bio15", "CHELSA_bio4", "NDVI_range", "soil_cec",
       #                "soil_Histosols", "CHELSA_bio2", "Tree_Cover", "soil_ocd", "Topo_northness", "Topo_tpi", "Forest_type",
        #               "Peat", "Soil_Depth", "Topo_geomflat", "Topo_geomhollow", "Topo_geompit", "Topo_geomslope", "Topo_geomvalley", "Water")

# Variables to remove after high correlation filtering (2025-12-4)
EXCLUDE_VARIABLES <- c("soil_silt", "bio04", "bio07", "gdd10", "gdgfgd5", "bio14", "bio18", "npp", "bio06", "bio08", "bio09", "bio05",
                       "scd", "bio16", "soil_bdod", "Topo_tri", "bio17", "soil_ocd", "bio01", "bio03", "bio10", "bio12", "fgd", "gdd0", "gdd5",
                       "gddlgd10", "gdgfgd10", "gsl", "gst", "lgd", "ngd0", "ngd10", "ngd5" # next variables were excluded after PCA1 filtering (18 variable remaining), updated 16 Jan
                       , "bio02", "bio15", "gddlgd0", "gsp", "soil_cfvo", "soil_cec", "soil_nitrogen", "soil_sand",
                       "soil_waterV1500", "Water", "Topo_tpi", "Tree_Cover", "Forest_type", "Soil_Depth", "Topo_geom")

# File paths - UPDATE THESE TO MATCH YOUR SETUP
base_dir <- here::here("MaxEnt-TrialRUNS-untuned_settings")
setwd(base_dir)

# File paths
species_file_MaxEnt <- here::here("Input_files/V2_cloudberry_presence_MaxEnt.csv")
env_vars_dir <- here::here("Input_files/Environmental_Variables-processed-current_climate")
maxent_jar <- "/Applications/maxent.jar" # set to macOS path, update to match your file path

#output_dir <- file.path(base_dir, "MaxEnt_Results_Single_Run1_Before_Cor_Filter") # before correlation filtering
#output_dir <- file.path(base_dir, "MaxEnt_Results_Single_Run2_After_Cor_Filter") # after correlation filtering
output_dir <- file.path(base_dir, "MaxEnt_Results_Single_Run3_After_PC1_filter") # after PCA1 filtering (2. Reduce multi-collinearity)
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

# Assign names once, immediately after loading
var_names <- tools::file_path_sans_ext(basename(raster_files))
names(env_raster) <- var_names

# Exclude specified variables by prefix
if (length(EXCLUDE_VARIABLES) > 0) {
  exclude_pattern <- paste0("^", EXCLUDE_VARIABLES, collapse = "|")
  keep_layers <- !grepl(exclude_pattern, names(env_raster))
  env_raster <- env_raster[[keep_layers]]
  cat("Excluded variables by prefix:", paste(EXCLUDE_VARIABLES, collapse = ", "), "\n")
}

cat("Loaded", nlyr(env_raster), "total environmental layers\n")

# =============================================================================
# IDENTIFY CATEGORICAL VARIABLES
# =============================================================================

cat("\n=== Identifying Categorical Variables ===\n")

categorical_prefixes <- c(
  "kg0", "Forest_type", "Peat", "Water", "Soil_Depth",
  "Topo_geomflat", "Topo_geomhollow", "Topo_geompit", "Topo_geomslope", "Topo_geomvalley"
)

# =============================================================================
# PRINT VARIABLES
# =============================================================================

# Print categorical variables
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

# Identify continuous variables as those not found in categorical_vars
continuous_vars <- setdiff(names(env_raster), categorical_vars)

# Print the continuous variables
cat("\n=== Continuous Variables ===\n")
if (length(continuous_vars) > 0) {
  cat("Found", length(continuous_vars), "continuous variables:\n")
  for (var in continuous_vars) {
    cat("  ", var, "\n")
  }
} else {
  cat("No continuous variables found\n")
}

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
# DATA SPLITTING AND CLEANING
# =============================================================================

cat("\n=== Data Preparation ===\n")
set.seed(123)

# Convert to raster stack for extraction
env_stack <- convert_terra_to_raster(env_raster)

# Split presence data for training/testing
n_total <- nrow(presence_points)
training_indices <- sample(1:n_total, round(0.75 * n_total))
training_points <- presence_points[training_indices, ]
testing_points <- presence_points[-training_indices, ]

# Clean data - remove points with NA values in any environmental layer
all_values <- raster::extract(env_stack, training_points)
complete_cases <- complete.cases(all_values)
clean_occurrence_points <- training_points[complete_cases, ]

test_values <- raster::extract(env_stack, testing_points)
test_complete_cases <- complete.cases(test_values)
testing_points <- testing_points[test_complete_cases, ]

cat("Clean training points:", nrow(clean_occurrence_points), "\n")
cat("Clean testing points:", nrow(testing_points), "\n")

# =============================================================================
# JAVA INITIALIZATION
# =============================================================================

# Initialize Java with error handling
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
# GENERATE BACKGROUND POINTS
# =============================================================================

cat("\n=== Generating Background Points ===\n")

set.seed(123)  # Ensure reproducibility
background_points <- randomPoints(
  env_stack,
  n = 88000, # extra points to ensure at least 50000 remain after NAs are removed
  p = clean_occurrence_points,
  excludep = TRUE, # excludes pixels with presence points
  prob = TRUE # more background sampling in cells closer to presence cells
)

# Remove background points with NA values
bg_env_values <- raster::extract(env_stack, background_points)
bg_complete <- complete.cases(bg_env_values)
background_points_clean <- background_points[bg_complete, ]
cat("Generated", nrow(background_points_clean), "clean background points\n")

# =============================================================================
# RUN MAXENT MODEL
# =============================================================================

cat("\n=== Running MaxEnt Model ===\n")
cat("Variables:", length(names(env_raster)), "\n")

# Update categorical variables
current_factors_vec <- NULL
if (length(categorical_vars) > 0) {
  current_factors_vec <- which(names(env_stack) %in% categorical_vars)
}

cat("Using", nrow(background_points_clean), "background points\n")

# Prepare MaxEnt arguments
maxent_args <- c(
  paste0("betamultiplier=", FIXED_REG_MULTIPLIER),
  "randomtestpoints=25",
  "maximumbackground=50000",
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
  "maximumiterations=5000", #increased to 5000
  "convergencethreshold=0.00001",
  "threads=1"
)

# Explicit feature control
all_features <- c("linear", "quadratic", "product", "threshold", "hinge")
for (feat in all_features) {
  if (feat %in% FIXED_FEATURES) {
    maxent_args <- c(maxent_args, paste0(feat, "=true"))
  } else {
    maxent_args <- c(maxent_args, paste0(feat, "=false"))
  }
}

# Run MaxEnt
cat("Running MaxEnt...\n")

tryCatch({
  model <- maxent(
    x = env_stack,
    p = clean_occurrence_points,
    path = output_dir,
    args = maxent_args,
    factors = current_factors_vec
  )
  
  # Create habitat suitability prediction
  if (!is.null(model)) {
    cat("Creating habitat suitability prediction...\n")
    prediction <- predict(model, env_stack)
    
    # Save prediction files
    asc_output <- file.path(output_dir, paste0("habitat_suitability_", FEATURE_NAME, "_",
                                               FIXED_REG_MULTIPLIER, ".asc"))
    
    raster::writeRaster(prediction, asc_output, format = "ascii", overwrite = TRUE)
    cat("Prediction files created successfully\n")
  }
  
  # Extract AUC values from results
  train_auc <- NA
  test_auc <- NA
  
  results_files <- list.files(output_dir, pattern = "maxentResults\\.csv$", full.names = TRUE)
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
  
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})

toc() # end timer
cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", output_dir, "\n")
