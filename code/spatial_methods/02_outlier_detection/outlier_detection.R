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
data(iris)

virginica <- iris[101:150, 1:3]

result <- mvn(data = virginica, mvnTest = "hz", 
              multivariateOutlierMethod = "quan")


# ========================================================================
# rrcovHD vignettes
# Notes: Seems like I can maybe input "~all_variables" as the function,
#        with the Spot x Gene matrix/dataframe  as the data? 
#
# Update: Actually would prob make senses to run just on typical qc variables.
#         Compare 1 vs many. 

data(hemophilia)
obj <- OutlierSign2(gr~.,data=hemophilia)
obj

getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                # returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function









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
pdf(here(plot_dir, 'SpotPlot_sum_umi.pdf'))
vis_gene(
  spe = spe.subset,
  geneid = "sum_umi",
  sampleid = unique(spe.subset$sample_id)[1]
)
dev.off()

pdf(here(plot_dir, 'SpotPlot_mito_percent.pdf'))
vis_gene(
  spe = spe.subset,
  geneid = "expr_chrM_ratio",
  sampleid = unique(spe.subset$sample_id)[1]
)
dev.off()

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


# =========== Outlier detection based on Mito percent ==========
var.mito<- rep(NA,nrow(spaQC))
z.mito <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.mito[i] <- var( spaQC[c(i, dnn.idx[dnn.idx != 0], i),]$expr_chrM_ratio, na.rm=TRUE)
  z.mito[i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$expr_chrM_ratio)[1]
}
z.mito[!is.finite(z.mito)] <- 0 

spe.subset$var.mito <- var.mito
spe.subset$z.mito <- z.mito
spe.subset$z_mito_4 <- ifelse(spe.subset$z.mito > 4, TRUE, FALSE)

# Visualize
pdf(here(plot_dir, 'SpotPlot_z_mito_k15.pdf'))
make_escheR(spe.subset) |> 
  add_fill(var = "z.mito") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")
dev.off()

pdf(here(plot_dir, 'SpotPlot_var_mito_k15.pdf'))
make_escheR(spe.subset) |> 
  add_fill(var = "var.mito") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")
dev.off()

pdf(here(plot_dir, 'SpotPlot_z_mito_k15_outliers.pdf'))
make_escheR(spe.subset) |> 
  add_fill(var = "z.mito") |>
  add_ground(var = "z_mito_4", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")
dev.off()

pdf(here(plot_dir, 'SpotPlot_mito_k15_outliers.pdf'))
make_escheR(spe.subset) |> 
  add_fill(var = "expr_chrM_ratio") |>
  add_ground(var = "z_mito_4", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")
dev.off()



# =========== Outlier detection based on UMI ==========
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
spe.subset$z_umi_4 <- ifelse(spe.subset$z.umi > 4, TRUE, FALSE)

# Visualize
pdf(here(plot_dir, 'SpotPlot_z_umi_k15.pdf'))
make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")
dev.off()

pdf(here(plot_dir, 'SpotPlot_var_umi_k15.pdf'))
make_escheR(spe.subset) |> 
  add_fill(var = "var.umi") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")
dev.off()

pdf(here(plot_dir, 'SpotPlot_z_umi_k15_outliers.pdf'))
make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") |>
  add_ground(var = "z_umi_4", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")
dev.off()

pdf(here(plot_dir, 'SpotPlot_sum_umi_k15_outliers.pdf'))
p <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi") |>
  add_ground(var = "z_umi_4", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")
print(p)
dev.off()