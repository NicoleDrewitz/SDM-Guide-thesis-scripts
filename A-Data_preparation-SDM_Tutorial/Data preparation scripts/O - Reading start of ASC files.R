# =============================================================================
# A-Data Processing / Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

#reading first rows in ASC. files
#edited from claude.ai on 25-08-2025 by Nicole Drewitz

library(here)
ASC_file <- here::here("Environmental_Variables/Environmental variable dataset for modelling/reprojected_bio12.asc") # example variable

if (!file.exists(ASC_file)) {
  stop("File does not exist. Check the file path.")
}

data <- read.table(ASC_file, header = TRUE, sep = "\t")
print(head(data, 6))