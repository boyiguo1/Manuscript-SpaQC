library(BiocNeighbors)
library(sp)
library(spdep)
library(RANN)
#library(spatialEco)
library(here)

plot_dir = here('plots', 'spatial_methods', '01_neighbor_detection')

# =========================================================================================================
# All of this was copied from https://gis.stackexchange.com/questions/219255/spatial-outliers-detection-in-r

data(meuse)
coordinates(meuse) <- ~x+y
head(meuse)

nobs <- 10000
ndim <- 20
data <- matrix(runif(nobs*ndim), ncol=ndim)
# [,1]      [,2]       [,3]       [,4]      [,5]      [,6]      [,7]       [,8]      [,9]     [,10]      [,11]
# [1,] 0.03999569 0.3063762 0.04116342 0.90822575 0.6086831 0.3656955 0.3173032 0.80650036 0.2671365 0.0116891 0.85672607
# [2,] 0.90794949 0.7141987 0.98547411 0.88268829 0.2444410 0.3773601 0.1092714 0.05001264 0.8936426 0.2324360 0.21100879
# [3,] 0.09326530 0.1424537 0.97923706 0.42764098 0.2792876 0.1315959 0.2422186 0.41572319 0.1431045 0.8976944 0.87982672
# [4,] 0.52113904 0.5303953 0.43748758 0.16208741 0.6244782 0.8989258 0.6169904 0.86513675 0.9974696 0.8283996 0.95791028
# [5,] 0.82561443 0.6865004 0.51275845 0.38547003 0.2863589 0.3365394 0.1517484 0.90974300 0.3991544 0.7048711 0.08341231
# [6,] 0.35335989 0.6913951 0.29567498 0.02712041 0.6092858 0.8854710 0.2327442 0.20069371 0.2960986 0.6355422 0.35491418
# [,12]     [,13]      [,14]      [,15]     [,16]     [,17]      [,18]      [,19]     [,20]
# [1,] 0.72864376 0.5927613 0.27305200 0.88200284 0.2789034 0.6135097 0.51889792 0.78569385 0.3769640
# [2,] 0.45017472 0.4273462 0.06570014 0.35449949 0.3298815 0.6017861 0.01328981 0.66957644 0.2987101
# [3,] 0.35463226 0.5266667 0.09262628 0.19684934 0.2436558 0.5580423 0.19427395 0.43330679 0.1218896
# [4,] 0.69220408 0.2516254 0.69512095 0.97185191 0.4154661 0.4248056 0.50566537 0.20146782 0.2913667
# [5,] 0.61259099 0.9487847 0.96369114 0.03496125 0.7831489 0.4160534 0.76786320 0.43995124 0.2831979
# [6,] 0.01407575 0.7379854 0.55780016 0.76624037 0.4417086 0.0140884 0.70392668 0.05055324 0.7651475

fout <- findKNN(coordinates(meuse), k=36)
head(fout$index)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,] 5003 4801 1918 9828 7666 6026 3138 9599 7145  1958
# [2,] 7074 2970 2760 4899 7113 4927 8443 8199 6960  1033
# [3,] 3723 8168 6276 4218 2547 6875 7901  701 7434   522
# [4,] 1515 9582 5704 3946 5651 7263 6007 8149 8060  9487
# [5,] 2598 8635 1504 2081 5621 7075 6723 7899 8657  4175
# [6,] 3156 3296 1159 6042  681 4995 9047 5000 3658   652


# ========================= Integrating this with SPE ============================
# Let's see how this works with SPE objects. 

library(SpatialExperiment)
library(spatialLIBD)
library(ggspavis)

spe <- fetch_data(type = "spe")
spe

spe.subset <- subset(spe, ,sample_id == "151507")

head(spatialCoords(spe.subset))
#.                   pxl_col_in_fullres pxl_row_in_fullres
# AAACAACGAATAGTTC-1               3276               2514
# AAACAAGTATCTCCCA-1               9178               8520
# AAACAATCTACTAGCA-1               5133               2878
# AAACACCAATAACTGC-1               3462               9581
# AAACAGCTTTCAGAAG-1               2779               7663
# AAACAGGGTCTATATT-1               3053               8143

dnn <- findKNN(spatialCoords(spe.subset), k=36)$index
head(dnn)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15]
# [1,] 3170 1086 1339 2809 4044 3715 2274 3064 2557  3389    83  2602  2079   248  1275
# [2,] 1204 4115  454 2434  551 2692  658 1024 3599   786  4182  1665  2339  3540  1438
# [3,] 3339 1708  788 2514 2699 1856 1717 1423 3607  2154   671  3920  3090   113  3689
# [4,] 2628 3856 2114  665 3515 2931 1133 2751 2590  1975  3400  3196   570  1435  3610
# [5,] 2066 2202  510 1476 3312 1377  346  438 3878  1734  2773  2141  4211  3897  1933
# [6,] 1598 1613 2301 1493 4070  309  653 2390 4130  2443  3601  3813   466  3062   764

# Add discrete values for random spot (50) and it's neighbors
spe.subset$nn_standard <- 0 # add zeros
spe.subset$nn_standard[50] <- 1 # spot = 1

neighbors_of_first_spot <- dnn[50,]
spe.subset$nn_standard[neighbors_of_first_spot] <- 2 # nn = 2
colnames(colData(spe.subset))
# [1] "sample_id"                   "Cluster"                     "sum_umi"                     "sum_gene"                   
# [5] "subject"                     "position"                    "replicate"                   "subject_position"           
# [9] "discard"                     "key"                         "cell_count"                  "SNN_k50_k4"                 
# [13] "SNN_k50_k5"                  "SNN_k50_k6"                  "SNN_k50_k7"                  "SNN_k50_k8"                 
# [17] "SNN_k50_k9"                  "SNN_k50_k10"                 "SNN_k50_k11"                 "SNN_k50_k12"                
# [21] "SNN_k50_k13"                 "SNN_k50_k14"                 "SNN_k50_k15"                 "SNN_k50_k16"                
# [25] "SNN_k50_k17"                 "SNN_k50_k18"                 "SNN_k50_k19"                 "SNN_k50_k20"                
# [29] "SNN_k50_k21"                 "SNN_k50_k22"                 "SNN_k50_k23"                 "SNN_k50_k24"                
# [33] "SNN_k50_k25"                 "SNN_k50_k26"                 "SNN_k50_k27"                 "SNN_k50_k28"                
# [37] "GraphBased"                  "Maynard"                     "Martinowich"                 "layer_guess"                
# [41] "layer_guess_reordered"       "layer_guess_reordered_short" "expr_chrM"                   "expr_chrM_ratio"            
# [45] "SpatialDE_PCA"               "SpatialDE_pool_PCA"          "HVG_PCA"                     "pseudobulk_PCA"             
# [49] "markers_PCA"                 "SpatialDE_UMAP"              "SpatialDE_pool_UMAP"         "HVG_UMAP"                   
# [53] "pseudobulk_UMAP"             "markers_UMAP"                "SpatialDE_PCA_spatial"       "SpatialDE_pool_PCA_spatial" 
# [57] "HVG_PCA_spatial"             "pseudobulk_PCA_spatial"      "markers_PCA_spatial"         "SpatialDE_UMAP_spatial"     
# [61] "SpatialDE_pool_UMAP_spatial" "HVG_UMAP_spatial"            "pseudobulk_UMAP_spatial"     "markers_UMAP_spatial"       
# [65] "spatialLIBD"                 "ManualAnnotation"            "in_tissue"                   "array_row"                  
# [69] "array_col"                   "nn_standard"  

head(spe.subset$nn_standard)
# [1] 0 0 0 0 0 0

unique(spe.subset$nn_standard)
# [1] 0 1 2


# Visualize
pdf(here(plot_dir, 'SpotPlot_BiocNeighbors_k36.pdf'))
vis_clus(
  spe = spe.subset,
  clustervar = "nn_standard",
  sampleid = unique(spe.subset$sample_id)[1]
)
dev.off()

