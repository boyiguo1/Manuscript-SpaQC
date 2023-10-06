#library(sp)
#library(spdep)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(ggspavis)
library(rrcovHD)

plot_dir = here('plots', 'spatial_methods', '02_outlier_detection')



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

# create a list of spatial coordinates and qc features
spaQC <- colData(spe.subset)
spaQC$coords <- spatialCoords(spe.subset)

# Find nearest neighbors
dnn <- RANN::nn2(spatialCoords(spe.subset), k=50, searchtype="standard")$nn.idx 
head(dnn)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    1 3170 1086 1339 2809 4044 3715 2274 3064  2557
# [2,]    2 1204 4115  454 2434 2692  551  658 1024  3599
# [3,]    3 3339 1708  788 2514 1856 2699 1717 1423  3607
# [4,]    4 2628 3856 2114  665 3515 2931 1133 2751  2590
# [5,]    5 2066 2202  510 1476 1377 3312  346  438  3878
# [6,]    6 1598 1613 2301 1493 4070  309  653 4130  2390


var.mito<- rep(NA,nrow(spaQC))
z.mito <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.mito[i] <- var( spaQC[dnn.idx[dnn.idx != 0],]$expr_chrM_ratio, na.rm=TRUE)
  z.mito[i] <- outliers(spaQC[dnn.idx[dnn.idx != 0],]$expr_chrM_ratio)[1]
}
z.mito[!is.finite(z.mito)] <- 0 

spe.subset$var.mito <- var.mito
spe.subset$z.mito <- z.mito

# Visualize
pdf(here(plot_dir, 'SpotPlot_z_mito_k50.pdf'))
vis_gene(
  spe = spe.subset,
  geneid = "z.mito",
  sampleid = unique(spe.subset$sample_id)[1]
)
dev.off()

pdf(here(plot_dir, 'SpotPlot_var_mito_k50.pdf'))
vis_gene(
  spe = spe.subset,
  geneid = "var.mito",
  sampleid = unique(spe.subset$sample_id)[1]
)
dev.off()






var.umi<- rep(NA,nrow(spaQC))
z.umi <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.umi[i] <- var( spaQC[dnn.idx[dnn.idx != 0],]$sum_umi, na.rm=TRUE)
  z.umi[i] <- outliers(spaQC[dnn.idx[dnn.idx != 0],]$sum_umi)[1]
}
z.umi[!is.finite(z.umi)] <- 0 

spe.subset$var.umi <- var.umi
spe.subset$z.umi <- z.umi

# Visualize
pdf(here(plot_dir, 'SpotPlot_z_umi_k50.pdf'))
vis_gene(
  spe = spe.subset,
  geneid = "z.umi",
  sampleid = unique(spe.subset$sample_id)[1]
)
dev.off()

pdf(here(plot_dir, 'SpotPlot_var_umi_k50.pdf'))
vis_gene(
  spe = spe.subset,
  geneid = "var.umi",
  sampleid = unique(spe.subset$sample_id)[1]
)
dev.off()