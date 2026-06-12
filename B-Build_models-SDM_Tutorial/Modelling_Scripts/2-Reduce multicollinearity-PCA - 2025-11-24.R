# =============================================================================
# Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

# PCA-Based Correlation Visualization for Environmental Variables
# by Nicole Drewitz, 30 October 2025
# Provides multiple visualizations to understand variable relationships
# Last updated 13 Jan, 2026

# lines 67-75 defines variables excluded from PCA (pretends they are all categorical)

# =============================================================================
# LOAD REQUIRED PACKAGES
# =============================================================================

required_packages <- c("FactoMineR", "factoextra", "corrplot", "ggplot2", 
                       "RColorBrewer", "gridExtra", "reshape2", "car", "terra", "raster", "here")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# CONFIGURATION
# =============================================================================

env_raster_dir <- here::here("Input_files/Environmental_Variables-processed-current_climate")
output_dir <- here::here("Output_files/Reducing_environmental_variable_set/PCA_Analysis")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Set high-quality graphics parameters
pdf_width <- 12
pdf_height <- 10

# =============================================================================
# LOAD ENVIRONMENTAL VARIABLES
# =============================================================================

cat("\n=== Loading Environmental Variables ===\n")
raster_files <- list.files(env_raster_dir, pattern = "\\.asc$", full.names = TRUE, ignore.case = TRUE)
if (length(raster_files) == 0) stop("No .asc files found!")

# Load with terra (handles .asc natively)
env_raster <- rast(raster_files)
cat("Loaded", length(raster_files), "raster files\n")
names(env_raster) <- tools::file_path_sans_ext(basename(raster_files))  # Clean names
print(env_raster)

# Extract variable names for later use
var_names <- names(env_raster)
cat("Variable names extracted:", length(var_names), "variables\n")
print(head(var_names))

# =============================================================================
# IDENTIFY CATEGORICAL VARIABLES AND ANY OTHERS NOT IN USE (EXCLUDE FROM PCA)
# =============================================================================
# Categorical variables or other excluded variables. Ensure bio1 is "bio1_" or bio11, bio12, etc will also be excluded.

cat("\n=== Identifying Variables to Exclude ===\n")
categorical_prefixes <- c( # January 2025
    "kg0", "Forest_type", "Peat", "Water", "Soil_Depth",
    "Topo_geomflat", "Topo_geomhollow", "Topo_geompit", "Topo_geomslope", "Topo_geomvalley", # continuous variables listed after this line
    "soil_silt",
    "bio04", "bio07", "gdd10", "gdgfgd5", "bio14", "bio18", "npp", "bio06", "bio08", "bio09", "bio05", #swe was not excluded
    "scd", "bio16", "soil_bdod", "Topo_tri", "bio17", "soil_ocd", "bio01", "bio03", "bio10", "bio12", "fgd", "gdd0", "gdd5",
    "gddlgd10", "gdgfgd10", "gsl", "gst", "lgd", "ngd0", "ngd10", "ngd5" # continuous variables after this line were removed after PCA_1 review
    , "bio02", "bio15", "gddlgd0", "gsp", "soil_cfvo", "soil_cec", "soil_nitrogen", "soil_sand",
    "soil_waterV1500", "Topo_tpi", "Tree_Cover"
    )

categorical_vars <- unlist(
  lapply(categorical_prefixes, function(prefix) var_names[startsWith(var_names, prefix)])
)
categorical_vars <- categorical_vars[categorical_vars %in% names(env_raster)]
continuous_vars <- setdiff(var_names, categorical_vars)

cat("Continuous variables:", length(continuous_vars), "\n")
cat("Categorical variables (excluded from PCA):", length(categorical_vars), "\n")
continuous_vars

# =============================================================================
# SAMPLE AND PREPARE DATA FOR PCA
# =============================================================================

cat("\n=== Sampling Raster Data ===\n")
env_continuous <- env_raster[[continuous_vars]]

set.seed(123)
n_sample <- min(60000, ncell(env_continuous))
sample_cells <- sample(1:ncell(env_continuous), n_sample)
env_data <- as.data.frame(env_continuous[sample_cells])
env_data <- na.omit(env_data)

cat("Sample size for PCA:", nrow(env_data), "cells\n")
cat("Variables in analysis:", ncol(env_data), "\n")

# =============================================================================
# SHORTEN VARIABLE NAMES FOR READABILITY
# =============================================================================

cat("\n=== Shortening Variable Names ===\n")

# Shorten variable names - remove everything after second underscore and replace hyphens
shorten_var_names <- function(var_names) {
  sapply(var_names, function(name) {
    # Replace hyphens with underscores first
    name <- gsub("-", "_", name)
    
    # Split on underscore and keep first two parts
    parts <- strsplit(name, "_")[[1]]
    if (length(parts) >= 2) {
      paste(parts[1], parts[2], sep = "_")
    } else {
      name
    }
  }, USE.NAMES = FALSE)
}

colnames(env_data) <- shorten_var_names(colnames(env_data))
cat("Variable names shortened for plotting\n")

# =============================================================================
# CALCULATE CORRELATION MATRIX
# =============================================================================

cat("\n=== Calculating Correlation Matrix ===\n")

# Test for normality
normality_test <- apply(env_data, 2, function(x) {
  if(length(x) > 3) {
    test_data <- if(length(x) > 5000) sample(x, 1000) else x
    shapiro.test(test_data)$p.value > 0.05
  } else TRUE
})

cor_method <- if (all(normality_test)) "pearson" else "spearman"
cat("Using", cor_method, "correlation\n")

cor_matrix <- cor(env_data, method = cor_method, use = "complete.obs")

# =============================================================================
# START SINGLE PDF OUTPUT AND DISPLAY IN PLOT WINDOW
# =============================================================================

cat("\n=== Creating Combined PDF Report ===\n")

# Open PDF device for combined output
pdf(file.path(output_dir, "PCA_Complete_Report.pdf"), width = pdf_width, height = pdf_height)

# =============================================================================
# VISUALIZATION 1: ENHANCED CORRELATION MATRIX
# =============================================================================

cat("\n=== Creating Enhanced Correlation Matrix ===\n")

corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.7,
         addCoef.col = "black",
         number.cex = 0.4,
         col = colorRampPalette(c("#1f78b4", "white", "#E46726"))(200),
         title = paste0("Correlation Matrix (", cor_method, ") - Hierarchical Clustering"),
         mar = c(0, 0, 2, 0))

# Create version with only high correlations labeled
cor_matrix_display <- cor_matrix
cor_matrix_display[abs(cor_matrix_display) < 0.75] <- NA

corrplot(cor_matrix_display,
         method = "color",
         type = "upper",
         order = "original",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.7,
         addCoef.col = "black",
         number.cex = 0.5,
         col = colorRampPalette(c("#1f78b4", "white", "#E46726"))(200),
         title = paste0("High Correlations Only (|r| > 0.75)"),
         mar = c(0, 0, 2, 0),
         insig = "blank")

cat("Added correlation matrices\n")

# =============================================================================
# VISUALIZATION 2: CORRELATION HEATMAP WITH DENDROGRAMS
# =============================================================================

cat("\n=== Creating Correlation Heatmap with Clustering ===\n")

# Hierarchical clustering
dist_matrix <- as.dist(1 - abs(cor_matrix))
hc <- hclust(dist_matrix, method = "ward.D2")

# Reorder correlation matrix
cor_matrix_ordered <- cor_matrix[hc$order, hc$order]

# Melt for ggplot
cor_melted <- melt(cor_matrix_ordered)
colnames(cor_melted) <- c("Var1", "Var2", "Correlation")

# Create heatmap
p_heatmap <- ggplot(cor_melted, aes(Var1, Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "#1f78b4", mid = "white", high = "#E46726",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = paste0(cor_method, "\nCorrelation")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right") +
  coord_fixed() +
  ggtitle("Correlation Heatmap (Hierarchically Clustered)")

print(p_heatmap)
cat("Added correlation heatmap\n")

# =============================================================================
# PERFORM PCA
# =============================================================================

cat("\n=== Performing PCA ===\n")

# Standardize data
env_data_scaled <- scale(env_data)

# Perform PCA
pca_result <- PCA(env_data_scaled, scale.unit = FALSE, graph = FALSE)

# Get summary
pca_summary <- summary(pca_result)
cat("\nVariance explained by first 10 components:\n")
print(pca_result$eig[1:min(10, nrow(pca_result$eig)), ])

# =============================================================================
# VISUALIZATION 3: SCREE PLOT
# =============================================================================

cat("\n=== Creating Scree Plot ===\n")

print(fviz_eig(pca_result, 
               addlabels = TRUE,
               ylim = c(0, 50),
               ncp = 20,
               main = "Scree Plot - Variance Explained by Principal Components",
               xlab = "Principal Component",
               ylab = "Percentage of Variance Explained (%)"))

cat("Added scree plot\n")

# =============================================================================
# VISUALIZATION 4: BIPLOT (VARIABLES AND INDIVIDUALS)
# =============================================================================

cat("\n=== Creating Biplot ===\n")

print(fviz_pca_biplot(pca_result,
                      repel = TRUE,
                      col.var = "#1f78b4",
                      col.ind = "#E46726",
                      alpha.ind = 0.3,
                      pointsize = 1,
                      labelsize = 3,
                      label = "var",
                      title = "PCA Biplot - PC1 vs PC2",
                      select.var = list(contrib = 30)))

cat("Added biplot PC1 vs PC2\n")

# =============================================================================
# VISUALIZATION 5: VARIABLE CONTRIBUTIONS TO PCs
# =============================================================================

cat("\n=== Creating Variable Contribution Plots ===\n")

# PC1
print(fviz_contrib(pca_result, 
                   choice = "var", 
                   axes = 1, 
                   top = 10,
                   title = "Top 10 Variable Contributions to PC1",
                   fill = "#1f78b4"))

# PC2
print(fviz_contrib(pca_result, 
                   choice = "var", 
                   axes = 2, 
                   top = 10,
                   title = "Top 10 Variable Contributions to PC2",
                   fill = "#E46726"))

# Combined PC1 and PC2
print(fviz_contrib(pca_result, 
                   choice = "var", 
                   axes = 1:2, 
                   top = 10,
                   title = "Top 10 Variable Contributions to PC1 & PC2",
                   fill = "#1b9e77"))

cat("Added contribution plots\n")

# =============================================================================
# VISUALIZATION 6: VARIABLE CORRELATION CIRCLE
# =============================================================================

cat("\n=== Creating Variable Correlation Circles ===\n")

# PC1 vs PC2
print(fviz_pca_var(pca_result,
                   col.var = "contrib",
                   gradient.cols = c("black", "#E46726", "blue"),
                   repel = TRUE,
                   labelsize = 3,
                   title = "Variable Correlation Circle - PC1 vs PC2"))

# PC3 vs PC4
print(fviz_pca_var(pca_result,
                   axes = c(3, 4),
                   col.var = "contrib",
                   gradient.cols = c("black", "#E46726", "blue"),
                   repel = TRUE,
                   labelsize = 3,
                   title = "Variable Correlation Circle - PC3 vs PC4"))

cat("Added variable correlation circles\n")

# =============================================================================
# VISUALIZATION 7: VARIABLE QUALITY OF REPRESENTATION (COS2)
# =============================================================================

cat("\n=== Creating Quality of Representation Plots ===\n")

print(fviz_pca_var(pca_result,
                   col.var = "cos2",
                   gradient.cols = c("black", "#E46726", "blue"),
                   repel = TRUE,
                   labelsize = 3,
                   title = "Variable Quality of Representation (cos2) - PC1 vs PC2"))

cat("Added quality of representation plot\n")

# =============================================================================
# VISUALIZATION 8: CORRELATION GROUPS/CLUSTERS
# =============================================================================

cat("\n=== Identifying Variable Groups ===\n")

# Use hierarchical clustering to identify groups
n_clusters <- 5
var_clusters <- cutree(hc, k = n_clusters)

# Colorblind friendly palette (Okabe-Ito), 5 colors
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

# Create dendrogram
plot(hc, 
     main = paste("Hierarchical Clustering of Variables (", n_clusters, "groups)"),
     xlab = "Variables",
     ylab = "Height (1 - |correlation|)",
     cex = 0.6)
rect.hclust(hc, k = n_clusters, border = cb_palette)

# Show cluster assignments
cluster_df <- data.frame(
  Variable = names(var_clusters),
  Cluster = var_clusters
)
cluster_df <- cluster_df[order(cluster_df$Cluster), ]

cat("\nVariable cluster assignments:\n")
for (i in 1:n_clusters) {
  cluster_vars <- cluster_df$Variable[cluster_df$Cluster == i]
  cat("\nCluster", i, "(", length(cluster_vars), "variables):\n")
  cat(paste("  -", cluster_vars, collapse = "\n"), "\n")
}

# =============================================================================
# VISUALIZATION 9: PCA WITH CLUSTER COLORING
# =============================================================================

cat("\n=== Creating PCA Plot with Cluster Colors ===\n")

print(fviz_pca_var(pca_result,
                   col.var = factor(var_clusters),
                   palette = brewer.pal(n_clusters, "Set2"),
                   legend.title = "Variable\nCluster",
                   repel = TRUE,
                   labelsize = 3,
                   title = "PCA Variables Colored by Correlation Cluster"))

cat("Added PCA plot colored by clusters\n")

# =============================================================================
# VISUALIZATION 10: CORRELATION NETWORK WITHIN TOP VARIABLES
# =============================================================================

cat("\n=== Creating Correlation Network for Top Variables ===\n")

# Get top contributing variables
var_contrib <- pca_result$var$contrib[, 1:2]
top_var_idx <- order(rowSums(var_contrib), decreasing = TRUE)[1:15]
top_vars <- rownames(var_contrib)[top_var_idx]

# Subset correlation matrix
top_vars_short <- intersect(top_vars, rownames(cor_matrix))
cor_matrix_top <- cor_matrix[top_vars_short, top_vars_short]
top_vars <- shorten_var_names(rownames(var_contrib)[top_var_idx]) # comfirm shortened names

corrplot(cor_matrix_top,
         method = "circle",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.8,
         addCoef.col = "black",
         number.cex = 0.5,
         col = colorRampPalette(c("#1f78b4", "white", "#E46726"))(200),
         title = "Correlation Network - Top 30 Contributing Variables",
         mar = c(0, 0, 2, 0))

cat("Added correlation network for top variables\n")

# =============================================================================
# VISUALIZATION 11: VIF and Corrplots for Each Cluster
# =============================================================================

cat("\n=== VIF and Corrplots for Each Variable Cluster ===\n")

# Helper function to calculate VIFs from a data frame (no categorical vars!)
calculate_vif <- function(df) {
  vifs <- sapply(names(df), function(var) {
    tryCatch({
      lm_formula <- as.formula(paste(var, "~ ."))
      model <- lm(lm_formula, data = df)
      r2 <- summary(model)$r.squared
      
      if (is.na(r2)) return(NA)
      if (r2 >= 0.9999) return(Inf)  # Near-perfect collinearity
      
      vif <- 1 / (1 - r2)
      return(vif)
    }, error = function(e) {
      warning(paste("VIF calculation failed for variable:", var, 
                    "\n  Error:", e$message))
      return(NA)
    })
  })
  vifs
}

# For each cluster, print title, VIF stats, and corrplot
for (cl in sort(unique(cluster_df$Cluster))) {
  cluster_vars <- cluster_df$Variable[cluster_df$Cluster == cl]
  cluster_dat <- env_data[, cluster_vars, drop = FALSE]
  
  # Skip clusters with fewer than 2 variables
  if (ncol(cluster_dat) < 2) next
  
  # Add new summary page in PDF
  plot.new()
  text(0.5, 0.92, paste0("VARIABLE CLUSTER ", cl, " (", ncol(cluster_dat), " vars)"), cex = 1.4, font = 2)
  text(0.5, 0.86, "Variance Inflation Factors (VIF)", cex = 1.1, font = 2)
  
  # Try/catch for VIF in case of collinearity or degenerate cases
  vif_vals <- tryCatch({ calculate_vif(cluster_dat) }, error=function(e) rep(NA, ncol(cluster_dat)))
  vif_text <- capture.output(print(round(vif_vals, 2)))
  y_pos <- 0.80
  for (i in seq_along(vif_text)) {
    text(0.05, y_pos, vif_text[i], cex = 0.85, family="mono", adj=0)
    y_pos <- y_pos-0.03
  }
  
  # Add correlation plot
  cor_mat_cluster <- cor(cluster_dat, use = "complete.obs", method = cor_method)
  corrplot(cor_mat_cluster, 
           method = "color",
           type = "upper",
           order = "hclust",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.95,
           addCoef.col = "black",
           number.cex = 0.55,
           col = colorRampPalette(c("#1f78b4", "white", "#E46726"))(200),
           title = paste0("Cluster ", cl, ": Within-Cluster Correlation Matrix"),
           mar = c(0, 0, 2, 0))
}

cat("Added VIF and cluster corrplots to the PDF\n")

# =============================================================================
# CREATE & EXPORT COMBINED VARIABLE SUMMARY + PCA LOADINGS
# =============================================================================

cat("\n=== Creating Combined Variable Summary + PCA Loadings ===\n")

# Base table with all variables
summary_table <- data.frame(
  Cluster = var_clusters[names(env_data)],
  Variable = names(env_data),
  stringsAsFactors = FALSE
)

# PCA Top 10 rankings
top10_pc1_pc2 <- names(sort(rowSums(pca_result$var$contrib[, 1:2]), decreasing = TRUE)[1:10])
top10_pc1 <- names(sort(pca_result$var$contrib[, 1], decreasing = TRUE)[1:10])
top10_pc2 <- names(sort(pca_result$var$contrib[, 2], decreasing = TRUE)[1:10])

summary_table$Top10_PC1_PC2 <- match(summary_table$Variable, top10_pc1_pc2)
summary_table$Top10_PC1 <- match(summary_table$Variable, top10_pc1)
summary_table$Top10_PC2 <- match(summary_table$Variable, top10_pc2)

# Calculate VIF per cluster (robust loop)
summary_table$VIF_within_cluster <- NA_real_
for (cl in unique(cluster_df$Cluster)) {
  cluster_vars <- cluster_df$Variable[cluster_df$Cluster == cl]
  if (length(cluster_vars) >= 2) {
    cluster_dat <- env_data[, cluster_vars, drop = FALSE]
    vif_vals <- tryCatch(calculate_vif(cluster_dat), error = function(e) NULL)
    if (!is.null(vif_vals) && length(vif_vals) == length(cluster_vars)) {
      cluster_rows <- summary_table$Variable %in% cluster_vars & summary_table$Cluster == cl
      summary_table$VIF_within_cluster[cluster_rows] <- vif_vals[match(
        summary_table$Variable[cluster_rows], cluster_vars)]
    }
  }
}

# Add PCA coordinates, contributions, cos2
pc_cols <- paste0("PC", 1:5)
summary_table[pc_cols] <- pca_result$var$coord[summary_table$Variable, 1:5]
summary_table$Contribution_PC1_PC2 <- rowSums(pca_result$var$contrib[summary_table$Variable, 1:2])
summary_table$Cos2_PC1_PC2 <- rowSums(pca_result$var$cos2[summary_table$Variable, 1:2])

# Final sort and format
summary_table <- summary_table[order(summary_table$Cluster, -summary_table$Contribution_PC1_PC2), ]
summary_table$VIF_within_cluster <- round(summary_table$VIF_within_cluster, 2)

# Export single CSV
write.csv(summary_table, file.path(output_dir, "variable_summary_pca_loadings.csv"), 
          row.names = FALSE, na = "")
cat("✓ Saved combined table:", nrow(summary_table), "variables\n")

# =============================================================================
# TEXT SUMMARY PAGES IN PDF
# =============================================================================

cat("\n=== Adding Text Summary to PDF ===\n")

# Create text summary pages
plot.new()
text(0.5, 0.95, "PCA ANALYSIS SUMMARY REPORT", cex = 1.8, font = 2)
text(0.5, 0.90, "Nordic Cloudberry (Rubus chamaemorus) Environmental Variables", cex = 1.2)
text(0.5, 0.85, paste("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), cex = 1)
text(0.5, 0.80, paste("Correlation Method:", cor_method), cex = 1)
text(0.5, 0.75, paste("Sample Size:", nrow(env_data), "cells"), cex = 1)
text(0.5, 0.70, paste("Number of Variables:", ncol(env_data)), cex = 1)

# Variance Explained
y_pos <- 0.60
text(0.5, y_pos, "VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS", cex = 1.3, font = 2)
y_pos <- y_pos - 0.05
var_exp_text <- capture.output(print(pca_result$eig[1:min(10, nrow(pca_result$eig)), ], digits = 3))
for (i in 1:min(12, length(var_exp_text))) {
  text(0.1, y_pos, var_exp_text[i], cex = 0.7, family = "mono", adj = 0)
  y_pos <- y_pos - 0.03
}

# Key findings
y_pos <- 0.15
text(0.5, y_pos, "KEY FINDINGS", cex = 1.3, font = 2)
y_pos <- y_pos - 0.04
text(0.5, y_pos, paste("First PC explains", round(pca_result$eig[1, 2], 2), "% of variance"), cex = 1)
y_pos <- y_pos - 0.04
text(0.5, y_pos, paste("First two PCs explain", round(sum(pca_result$eig[1:2, 2]), 2), "% of variance"), cex = 1)
y_pos <- y_pos - 0.04
text(0.5, y_pos, paste("Identified", n_clusters, "correlation-based variable clusters"), cex = 1)

# Page 2: Top Contributing Variables
plot.new()
text(0.5, 0.95, "TOP 20 CONTRIBUTING VARIABLES (PC1 & PC2)", cex = 1.5, font = 2)
y_pos <- 0.88
top20_text <- capture.output(print(head(loadings_df[, c("Variable", "Dim.1", "Dim.2", "Contribution_PC1_PC2")], 20), 
                                   row.names = FALSE, digits = 3))
for (i in 1:min(23, length(top20_text))) {
  text(0.05, y_pos, top20_text[i], cex = 0.65, family = "mono", adj = 0)
  y_pos <- y_pos - 0.03
}

# Page 3: Variable Clusters
plot.new()
text(0.5, 0.95, "VARIABLE CLUSTERS (based on correlation)", cex = 1.5, font = 2)
y_pos <- 0.88
for (i in 1:n_clusters) {
  cluster_vars <- cluster_df$Variable[cluster_df$Cluster == i]
  text(0.05, y_pos, paste("Cluster", i, "(", length(cluster_vars), "variables):"), 
       cex = 1, font = 2, adj = 0)
  y_pos <- y_pos - 0.03
  
  # Print variables in this cluster
  for (var in cluster_vars) {
    if (y_pos < 0.05) {
      plot.new()
      text(0.5, 0.95, "VARIABLE CLUSTERS (continued)", cex = 1.5, font = 2)
      y_pos <- 0.88
    }
    text(0.08, y_pos, paste("-", var), cex = 0.8, adj = 0)
    y_pos <- y_pos - 0.025
  }
  y_pos <- y_pos - 0.02
  
  if (y_pos < 0.1 && i < n_clusters) {
    plot.new()
    text(0.5, 0.95, "VARIABLE CLUSTERS (continued)", cex = 1.5, font = 2)
    y_pos <- 0.88
  }
}

# Page 4: High Correlation Pairs
plot.new()
text(0.5, 0.95, "HIGH CORRELATION PAIRS (|r| > 0.75)", cex = 1.5, font = 2)

high_cor_pairs <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(ncol(cor_matrix)-1)) {
  for (j in (i+1):ncol(cor_matrix)) {
    if (abs(cor_matrix[i, j]) > 0.75) {
      high_cor_pairs <- rbind(high_cor_pairs, data.frame(
        Variable1 = rownames(cor_matrix)[i],
        Variable2 = colnames(cor_matrix)[j],
        Correlation = cor_matrix[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
}

high_cor_pairs <- high_cor_pairs[order(abs(high_cor_pairs$Correlation), decreasing = TRUE), ]

y_pos <- 0.88
if (nrow(high_cor_pairs) > 0) {
  high_cor_text <- capture.output(print(high_cor_pairs, row.names = FALSE, digits = 3))
  for (i in 1:min(30, length(high_cor_text))) {
    if (y_pos < 0.05) {
      plot.new()
      text(0.5, 0.95, "HIGH CORRELATION PAIRS (continued)", cex = 1.5, font = 2)
      y_pos <- 0.88
    }
    text(0.05, y_pos, high_cor_text[i], cex = 0.65, family = "mono", adj = 0)
    y_pos <- y_pos - 0.03
  }
  if (length(high_cor_text) > 30) {
    text(0.5, y_pos - 0.05, paste("... and", length(high_cor_text) - 30, "more pairs"), cex = 1, font = 3)
  }
} else {
  text(0.5, 0.5, "No high correlation pairs found (|r| > 0.75)", cex = 1.2)
}

# Close PDF device
dev.off()

cat("Completed combined PDF report\n")

cat("===============================================================================\n")
cat("PCA ANALYSIS SUMMARY REPORT\n")
cat("SDM Environmental Variables\n")
cat("===============================================================================\n\n")

cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Correlation Method:", cor_method, "\n")
cat("Sample Size:", nrow(env_data), "cells\n")
cat("Number of Variables:", ncol(env_data), "\n\n")

cat("===============================================================================\n")
cat("VARIANCE EXPLAINED\n")
cat("===============================================================================\n\n")
print(pca_result$eig[1:min(10, nrow(pca_result$eig)), ])

cat("\n===============================================================================\n")
cat("TOP 20 CONTRIBUTING VARIABLES (PC1 & PC2)\n")
cat("===============================================================================\n\n")
print(head(loadings_df[, c("Variable", "Dim.1", "Dim.2", "Contribution_PC1_PC2")], 20), 
      row.names = FALSE)

cat("\n===============================================================================\n")
cat("VARIABLE CLUSTERS (based on correlation)\n")
cat("===============================================================================\n\n")
for (i in 1:n_clusters) {
  cluster_vars <- cluster_df$Variable[cluster_df$Cluster == i]
  cat("Cluster", i, "(", length(cluster_vars), "variables):\n")
  cat(paste("  -", cluster_vars, collapse = "\n"), "\n\n")
}

cat("\n===============================================================================\n")
cat("HIGH CORRELATION PAIRS (|r| > 0.75)\n")
cat("===============================================================================\n\n")

high_cor_pairs <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(ncol(cor_matrix)-1)) {
  for (j in (i+1):ncol(cor_matrix)) {
    if (abs(cor_matrix[i, j]) > 0.75) {
      high_cor_pairs <- rbind(high_cor_pairs, data.frame(
        Variable1 = rownames(cor_matrix)[i],
        Variable2 = colnames(cor_matrix)[j],
        Correlation = cor_matrix[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
}

high_cor_pairs <- high_cor_pairs[order(abs(high_cor_pairs$Correlation), decreasing = TRUE), ]
print(high_cor_pairs, row.names = FALSE)

cat("\n===============================================================================\n")
cat("OUTPUT FILES CREATED\n")
cat("===============================================================================\n\n")

output_files <- list.files(output_dir, pattern = "\\.(pdf|csv|txt)$")
cat(paste(output_files, collapse = "\n"), "\n")

cat("\n===============================================================================\n")
cat("END OF REPORT\n")
cat("===============================================================================\n")