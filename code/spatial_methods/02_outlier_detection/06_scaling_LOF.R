library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(ggspavis)
library(rrcovHD)
library(mvoutlier)
library(BiocNeighbors)
library(escheR)
library(patchwork)
library(dbscan)

plot_dir = here('plots', 'spatial_methods', '06_outlier_detection')
processed_dir = here('processed-data', 'spatial_methods')

## Download the spot-level data
spe <- fetch_data(type = "spe")
spe




localOutliers <- function(spe, k = 36, feature='sum_umi', samples='sample_id', log2=TRUE, z_threshold = 3, output_z=FALSE) {
  
  # log2 transform the sum_umi and sum_gene features
  if (log2) {
    feature_log2 <- paste0(feature, '_log2')
    colData(spe)[feature_log2] <- log2(colData(spe)[[feature]])
    feature2use <- feature_log2
  } else {
    feature2use <- feature
  }
  
  # Get a list of unique sample IDs
  unique_sample_ids <- unique(colData(spe)[[samples]])
  
  # Initialize list variables to store the results
  mod_z <- vector("list", length(unique_sample_ids))
  
  # Initialize a list to store each spaQC dataframe
  spaQC_list <- vector("list", length(unique_sample_ids))
  
  # Loop through each unique sample ID
  for (sample_id in seq_along(unique_sample_ids)) {
    # Subset the data for the current sample
    sample <- unique_sample_ids[sample_id]
    spe_subset <- subset(spe, , sample_id == sample)
    
    # Create a list of spatial coordinates and qc features
    spaQC <- colData(spe_subset)
    spaQC$coords <- spatialCoords(spe_subset)
    
    # Find nearest neighbors
    dnn <- findKNN(spatialCoords(spe_subset), k = k)$index
    
    # Initialize variables for the current sample
    mod_z[[sample_id]] <- rep(NA, nrow(spaQC))
    
    # Loop through each row in the nearest neighbor index matrix
    for(i in 1:nrow(dnn)) {
      dnn.idx <- dnn[i,]
      mod_z[[sample_id]][i] <- modified_z(spaQC[c(i, dnn.idx[dnn.idx != 0]),][[feature_log2]])[1]
    }



# ============ Scaled up ==============
df <- data.frame(
  mito = spe.subset$z.mito,
  umi = spe.subset$z.umi,
  gene = spe.subset$z.gene
)

# Create a list of spatial coordinates and qc features
spaQC_total <- colData(spe)
spaQC_total$coords <- spatialCoords(spe)

# Get a list of unique sample IDs
unique_sample_ids <- unique(spe$sample_id)

# Initialize list variables to store the results
var.umi <- vector("list", length(unique_sample_ids))
z.umi <- vector("list", length(unique_sample_ids))

# Initialize a list to store each spaQC dataframe
spaQC_list <- vector("list", length(unique_sample_ids))


# Loop through each unique sample ID
for(sample_id in seq_along(unique_sample_ids)) {
  
  # Subset the data for the current sample
  sample <- unique_sample_ids[sample_id]
  spe_subset <- subset(spe, ,sample_id == sample)
  
  # ===================== Find Neighbors ==============================
  # Create a list of spatial coordinates and qc features
  spaQC <- colData(spe_subset)
  spaQC$coords <- spatialCoords(spe_subset)