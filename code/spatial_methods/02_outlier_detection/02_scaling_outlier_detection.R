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
library(tictoc)

plot_dir = here('plots', 'spatial_methods', '02_outlier_detection')
processed_dir = here('processed-data', 'spatial_methods')

load(here('raw-data', "dlPFC_raw.Rdata"))
spe <- spe_raw
spe

rm(spe_raw)

# =========== Outlier function ==============
# as of 10/17/23, the spatialEco function was removed from CRAN due to a bad build
# so I'm copying their outlier function here.

outliers <- function(x, s = 1.4826) {
  e <- (length(x) - 1) / sqrt(length(x)) 
  mad <- function (x, center=stats::median(x), constant=s,
                   low=FALSE, high=FALSE) {
    n <- length(x)
    constant * if ((low || high) && n%%2 == 0) {
      if (low && high) 
        stop("'low' and 'high' cannot be both TRUE")
      n2 <- n%/%2 + as.integer(high)
      sort(abs(x - center), partial = n2)[n2]
    }
    else stats::median(abs(x - center))
  }                         
  return( ( (0.6745 * (x - stats::median(x))) / mad(x) ) )
}




# ===================== Find Neighbors ==============================
# create a list of spatial coordinates and qc features
spaQC <- colData(spe.subset)
spaQC$coords <- spatialCoords(spe.subset)

# Find nearest neighbors
dnn <- findKNN(spatialCoords(spe.subset), k=15)$index



# =========== Outlier detection based on Mito percent ==========
var.mito<- rep(NA,nrow(spaQC))
z.mito <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.mito[i] <- var( spaQC[c(i, dnn.idx[dnn.idx != 0]),]$expr_chrM_ratio, na.rm=TRUE)
  z.mito[i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$expr_chrM_ratio)[1]
}
z.mito[!is.finite(z.mito)] <- 0 

spe.subset$var.mito <- var.mito
spe.subset$z.mito <- z.mito
spe.subset$z_mito_outlier <- ifelse(spe.subset$z.mito > 3, TRUE, FALSE)





# ============ Scaled up ==============

# Create a list of spatial coordinates and qc features
spaQC_total <- colData(spe)
spaQC_total$coords <- spatialCoords(spe)

# Get a list of unique sample IDs
unique_sample_ids <- unique(spe$sample_id)

# Initialize list variables to store the results
var.mito <- vector("list", length(unique_sample_ids))
z.mito <- vector("list", length(unique_sample_ids))

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
  
  # Find nearest neighbors
  dnn <- findKNN(spatialCoords(spe_subset), k=15)$index
  
  # =========== Outlier detection based on Mito percent ==========
  # Initialize variables for the current sample
  var.mito[[sample_id]] <- rep(NA, nrow(spaQC))
  z.mito[[sample_id]] <- rep(NA, nrow(spaQC))
  
  # Loop through each row in the nearest neighbor index matrix
  for(i in 1:nrow(dnn)) {
    dnn.idx <- dnn[i,]
    var.mito[[sample_id]][i] <- var(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$expr_chrM_ratio, na.rm=TRUE)
    z.mito[[sample_id]][i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$expr_chrM_ratio)[1]
  }
  
  # Handle non-finite values
  z.mito[[sample_id]][!is.finite(z.mito[[sample_id]])] <- 0
  
  spaQC$var.mito <- var.mito[[sample_id]]
  spaQC$z.mito <- z.mito[[sample_id]]
  spaQC$z_mito_outlier <- ifelse(spaQC$z.mito[[sample_id]] > 3 | spaQC$z.mito[[sample_id]] < -3, TRUE, FALSE)
  
  # Store the modified spaQC dataframe in the list
  spaQC_list[[sample_id]] <- spaQC
  
  
}

spaQC_aggregated <- do.call(rbind, spaQC_list)

colData(spe) <- spaQC_aggregated
colnames(colData(spe))
# [1] "sample_id"              "in_tissue"              "array_row"              "array_col"              "10x_graphclust"        
# [6] "10x_kmeans_10_clusters" "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters" 
# [11] "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters"  "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                   
# [16] "sum_umi"                "sum_gene"               "expr_chrM"              "expr_chrM_ratio"        "ManualAnnotation"      
# [21] "subject"                "region"                 "sex"                    "age"                    "diagnosis"             
# [26] "sample_id_complete"     "count"                  "coords"                 "var.mito"               "z.mito"                
# [31] "z_mito_outlier"     
