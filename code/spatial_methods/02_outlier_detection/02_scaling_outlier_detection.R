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

## Download the spot-level data
spe <- fetch_data(type = "spe")
spe


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



# ============ Scaled up ==============
# log2 transform the sum_umi and sum_gene features
colData(spe)$sum_umi_log2 <- log2(spe$sum_umi)
colData(spe)$sum_gene_log2 <- log2(spe$sum_gene)

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
  
  # Find nearest neighbors
  dnn <- findKNN(spatialCoords(spe_subset), k=15)$index
  
  # =========== Outlier detection based on Mito percent ==========
  # Initialize variables for the current sample
  var.umi[[sample_id]] <- rep(NA, nrow(spaQC))
  z.umi[[sample_id]] <- rep(NA, nrow(spaQC))
  
  # Loop through each row in the nearest neighbor index matrix
  for(i in 1:nrow(dnn)) {
    dnn.idx <- dnn[i,]
    var.umi[[sample_id]][i] <- var(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_umi_log2, na.rm=TRUE)
    z.umi[[sample_id]][i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_umi_log2)[1]
  }
  
  # Handle non-finite values
  z.umi[[sample_id]][!is.finite(z.umi[[sample_id]])] <- 0
  
  spaQC$var.umi <- var.umi[[sample_id]]
  spaQC$z.umi <- z.umi[[sample_id]]
  spaQC$z_umi_outlier <- ifelse(spaQC$z.umi > 3 | spaQC$z.umi < -3, TRUE, FALSE)
  
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




# ====================== Explore these outliers ======================

spe.subset <- subset(spe, ,sample_id == spe$sample_id[1])
# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p3 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") |>
  add_ground(var = "z_umi_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p4 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") |>
  add_ground(var = "z_umi_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

pdf(width = 12.5, height = 12.5, here(plot_dir, 'Many_SpotPlots_sumUmi_k15.pdf'))
(p1+p2)/(p3+p4)
dev.off()