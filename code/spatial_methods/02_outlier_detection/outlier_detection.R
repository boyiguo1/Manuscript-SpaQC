#library(sp)
#library(spdep)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(ggspavis)
library(rrcovHD)
library(mvoutlier)
library(BiocNeighbors)
#library(spatialEco)
library(escheR)
library(patchwork)



plot_dir = here('plots', 'spatial_methods', '02_outlier_detection')


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

# ============= Playing around with multivariate outlier detection ============
# load Iris data
# data(iris)

# virginica <- iris[101:150, 1:3]

# result <- mvn(data = virginica, mvnTest = "hz", 
              # multivariateOutlierMethod = "quan")


# ========================================================================
# rrcovHD vignettes
# Notes: Seems like I can maybe input "~all_variables" as the function,
#        with the Spot x Gene matrix/dataframe  as the data? 
#
# Update: Actually would prob make senses to run just on typical qc variables.
#         Compare 1 vs many. 

# data(hemophilia)
# obj <- OutlierSign2(gr~.,data=hemophilia)
# obj
# 
# getDistance(obj)            # returns an array of distances
# getClassLabels(obj, 1)      # returns an array of indices for a given class
# getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
# getFlag(obj)                # returns an 0/1 array of flags
# plot(obj, class=2)          # standard plot function







# =========== subset data to one sample and plot standard QC metrics ==============

# get spe
spe <- fetch_data(type = "spe")
spe

# subset spe to one sample for testing purposes
spe.subset <- subset(spe, ,sample_id == "151507")
head(spatialCoords(spe.subset))
#.                   pxl_col_in_fullres pxl_row_in_fullres
# AAACAACGAATAGTTC-1               3276               2514
# AAACAAGTATCTCCCA-1               9178               8520
# AAACAATCTACTAGCA-1               5133               2878
# AAACACCAATAACTGC-1               3462               9581
# AAACAGCTTTCAGAAG-1               2779               7663
# AAACAGGGTCTATATT-1               3053               8143


# Visualize UMI and Mito 
p1 <- vis_gene(
  spe = spe.subset,
  geneid = "expr_chrM_ratio",
  sampleid = unique(spe.subset$sample_id)[1]
)

p2 <- vis_gene(
  spe = spe.subset,
  geneid = "sum_umi",
  sampleid = unique(spe.subset$sample_id)[1]
)


p3 <- vis_gene(
  spe = spe.subset,
  geneid = "sum_gene",
  sampleid = unique(spe.subset$sample_id)[1]
)

pdf(height=10, width=20, here(plot_dir, 'SpotPlots_QCmetrics_raw.pdf'))
(p1+p2+p3)
dev.off()



# ===================== Spatially aware QC ==============================
# create a list of spatial coordinates and qc features
spaQC <- colData(spe.subset)
spaQC$coords <- spatialCoords(spe.subset)

# Check the first few entries of your main data objects
head(spaQC)
head(spatialCoords(spe.subset))



# Find nearest neighbors
dnn <- findKNN(spatialCoords(spe.subset), k=15)$index
head(dnn)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15]
# [1,] 3170 1086 1339 2809 4044 3715 2274 3064 2557  3389    83  2602  2079   248  1275
# [2,] 1204 4115  454 2434  551 2692  658 1024 3599   786  4182  1665  2339  3540  1438
# [3,] 3339 1708  788 2514 1856 2699 1717 1423 3607  2154   671  3920  3090   113  3689
# [4,] 2628 3856 2114  665 3515 2931 1133 2751 2590  1975  3196  3400   570  1435  3610
# [5,] 2066 2202  510 1476 1377 3312  346  438 3878  1734  2773  2141  4211  3897  1933
# [6,] 1598 1613 2301 1493 4070  309  653 2390 4130  2443  3601  3813   466  3062   764

# NOTE: findKNN does NOT include the center spot in this output. Other methods do (RANN).


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

# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "expr_chrM_ratio") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.mito") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p3 <- make_escheR(spe.subset) |> 
  add_fill(var = "expr_chrM_ratio") |>
  add_ground(var = "z_mito_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p4 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.mito") |>
  add_ground(var = "z_mito_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

pdf(width = 12.5, height = 12.5, here(plot_dir, 'SpotPlots_mitoPercent_k15.pdf'))
(p1+p2)/(p3+p4)
dev.off()



# =========== Outlier detection based on sum UMI ==========
var.umi<- rep(NA,nrow(spaQC))
z.umi <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.umi[i] <- var( spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_umi, na.rm=TRUE)
  z.umi[i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_umi)[1]
}
z.umi[!is.finite(z.umi)] <- 0 

spe.subset$var.umi <- var.umi
spe.subset$z.umi <- z.umi
spe.subset$z_umi_outlier <- ifelse(spe.subset$z.umi > 3 | spe.subset$z.umi < -3, TRUE, FALSE)

# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p3 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi") |>
  add_ground(var = "z_umi_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

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

pdf(width = 12.5, height = 12.5, here(plot_dir, 'SpotPlots_UMI_k15.pdf'))
(p1+p2)/(p3+p4)
dev.off()



# =========== Outlier detection based on sum gene ==========
var.gene<- rep(NA,nrow(spaQC))
z.gene <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.gene[i] <- var( spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_gene, na.rm=TRUE)
  z.gene[i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_gene)[1]
}
z.gene[!is.finite(z.gene)] <- 0 

spe.subset$var.gene <- var.gene
spe.subset$z.gene <- z.gene
spe.subset$z_gene_outlier <- ifelse(spe.subset$z.gene > 3 | spe.subset$z.gene < -3, TRUE, FALSE)

# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_gene") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.gene") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p3 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_gene") |>
  add_ground(var = "z_gene_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p4 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.gene") |>
  add_ground(var = "z_gene_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

pdf(width = 12.5, height = 12.5, here(plot_dir, 'SpotPlots_Gene_k15.pdf'))
(p1+p2)/(p3+p4)
dev.off()



# ================ Let's now look at sum/detected genes on log2 scale ==================


colData(spe)$sum_umi_log2 <- log2(spe$sum_umi)
colData(spe)$sum_gene_log2 <- log2(spe$sum_gene)

# Extract data from colData of the SPE object
layer_data <- colData(spe)$layer_guess_reordered
sum_gene <- colData(spe)$sum_gene
sum_umi <- colData(spe)$sum_umi
sum_umi_log2 <- colData(spe)$sum_umi_log2
sum_gene_log2 <- colData(spe)$sum_gene_log2



# Check that the lengths of the extracted data are the same
if(length(expr_data) != length(layer_data)) {
  stop("Mismatch in the lengths of expr_data and layer_data!")
}

# Combine the expression data with the layer data
combined_data <- data.frame(
  layer = layer_data,
  sum_umi = sum_umi,
  sum_gene = sum_gene,
  sum_umi_log2 = sum_umi_log2,
  sum_gene_log2 = sum_gene_log2
)

# Let's first plot the distributions of log2 scaled counts 



p1 <- ggplot(combined_data, aes(x=sum_gene, fill=layer)) +
  geom_histogram(alpha=0.5, position="identity", bins=30) +
  theme_minimal() +
  labs(title="Raw count of Sum Genes", x="Sum Genes", y="Count", fill="Layer")

p2 <- ggplot(combined_data, aes(x=sum_gene_log2, fill=layer)) +
  geom_histogram(alpha=0.5, position="identity", bins=30) +
  theme_minimal() +
  labs(title="Log2 transformation of Sum Genes", x="log2(Sum Genes)", y="Count", fill="Layer")

pdf(height = 5, width=10, here(plot_dir, 'Histogram_sumGenes_vs_log2.pdf'))
(p1+p2)
dev.off()


p1 <- ggplot(combined_data, aes(x=sum_umi, fill=layer)) +
  geom_histogram(alpha=0.5, position="identity", bins=30) +
  theme_minimal() +
  labs(title="Raw count of Sum UMI", x="Sum UMI", y="Count", fill="Layer")

p2 <- ggplot(combined_data, aes(x=sum_umi_log2, fill=layer)) +
  geom_histogram(alpha=0.5, position="identity", bins=30) +
  theme_minimal() +
  labs(title="Log2 transformation of Sum UMI", x="Sum UMI", y="Count", fill="Layer")

pdf(height = 5, width=10, here(plot_dir, 'Histogram_sumUMI_vs_log2.pdf'))
(p1+p2)
dev.off()