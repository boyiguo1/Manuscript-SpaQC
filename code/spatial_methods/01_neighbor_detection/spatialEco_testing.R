library(sp)
library(spdep)
library(RANN)
library(spatialEco)
library(here)


plot_dir = here('plots', 'spatial_methods', '01_neighbor_detection')

# =========================================================================================================
# All of this was copied from https://gis.stackexchange.com/questions/219255/spatial-outliers-detection-in-r

data(meuse)
coordinates(meuse) <- ~x+y
head(meuse)
# coordinates cadmium copper lead zinc  elev       dist   om ffreq soil lime landuse dist.m
# 1 (181072, 333611)    11.7     85  299 1022 7.909 0.00135803 13.6     1    1    1      Ah     50
# 2 (181025, 333558)     8.6     81  277 1141 6.983 0.01222430 14.0     1    1    1      Ah     30
# 3 (181165, 333537)     6.5     68  199  640 7.800 0.10302900 13.0     1    1    1      Ah    150
# 4 (181298, 333484)     2.6     81  116  257 7.655 0.19009400  8.0     1    2    0      Ga    270
# 5 (181307, 333330)     2.8     48  117  269 7.480 0.27709000  8.7     1    2    0      Ah    380
# 6 (181390, 333260)     3.0     61  137  281 7.791 0.36406700  7.8     1    2    0      Ga    470

# Here we can calculate global outliers using modified z-score. 
#  None of the z values in this variable exceed the common threshold (ie., z>9) but, 
#  by plotting the z values you can see the spatial distribution of z.

( meuse$Zscore <- spatialEco::outliers(meuse$cadmium) )  
spplot(meuse, "Zscore", col.regions=cm.colors(10))

# Here we calculate local (1000m radius) variation using a distance-based neighbor 
#  variance and the local z-score. This is simply a statistic of each point including 
#  all neighbors within the specified distance. You could apply a large variety of 
#  statistics using this approach. In this example, we do identify one significant 
#  local outlier using the local modified z-score.

# this RNN::nn2 function appears to do something very similar/same to BiocNeighors
dnn <- RANN::nn2(coordinates(meuse), searchtype="radius", 
                 radius = 1000)$nn.idx 
var.cadmium <- rep(NA,nrow(meuse))
z.cadmium <- rep(NA,nrow(meuse))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.cadmium[i] <- var( meuse[dnn.idx[dnn.idx != 0],]$cadmium, na.rm=TRUE)
  z.cadmium[i] <- outliers(meuse[dnn.idx[dnn.idx != 0],]$cadmium)[1]
}
z.cadmium[!is.finite(z.cadmium)] <- 0 

meuse$var.cadmium <- var.cadmium
spplot(meuse, "var.cadmium", col.regions=cm.colors(10))

meuse$z.cadmium <- z.cadmium
spplot(meuse, "z.cadmium", col.regions=cm.colors(10))


# Here is where we calculate the local autocorrelation using the Local Moran's-I 
# or LISA statistic. First, we need to calculate a minimum search distance to 
# ensure that we do not have any empty sets (null Wij neighbor matrix) and then
# use the min distance to build the Wij spatial weights matrix.


all.linked <- max(unlist(nbdists(knn2nb(knearneigh(coordinates(meuse))), 
                                 coordinates(meuse))))
nb <- dnearneigh(meuse, 0, all.linked)

# Now we can use the Wij matrix to calculate the Local Moran's-I and
# create a point dataset, corresponding to meuse, to explore results. 
# We also print and plot the number of observations that are statistically significant.

mI <- localmoran(meuse@data[,"lead"], nb2listw(nb, style="W"))
LocalI <- meuse
LocalI@data <- data.frame(ID=rownames(LocalI@data), as.data.frame(mI))
names(LocalI@data)[6] <- "Pr"
spplot(LocalI, "Z.Ii", xlab="Local Morans-I", col.regions=topo.colors(30))   

cat(nrow( LocalI[LocalI@data[,"Pr"] < 0.05 ,]), "obs of", 
    nrow(LocalI), "are significant at p=0.05","\n")

plot(LocalI, pch=19)
points(LocalI[which(LocalI$Pr <= 0.05),], pch=19,col="red")

# The high and low have no distinct I values so, one must create a vector to distinguish 
# significant and high hotspots. Red represents hot spots or spatial outliers.

LocalI@data <- data.frame(LocalI@data, HotSpots=ifelse( mI[,5] <= 0.05 & mI[,4] >= mean(mI[,4]), 1, 0) )
LocalI@data$HotSpots <- as.factor(LocalI@data$HotSpots)

spplot(LocalI, "HotSpots", xlab="Local Moran’s-I Hot Spots", col.regions=c("blue","red") )



# =========================================================================================================
# This definitely looks useful. RANN::nn2 appears similar to BiocNeighbors that Boyi
# shared. It gives (IMO) more intuitive neighbor definitions, such as neighbores based on
# radius; although radius may not be the best metric. 
#
# TODO: it would be worth seeing if RANN:nn2 is faster than BiocNeighbors. 
# TODO: Generate some visuals for calculated neighbors (similar to above "Hot Spots")


# Visualize neighbor detection using various methods for RANN::nn2


data(meuse)
coordinates(meuse) <- ~x+y
head(meuse)
# coordinates cadmium copper lead zinc  elev       dist   om ffreq soil lime landuse dist.m
# 1 (181072, 333611)    11.7     85  299 1022 7.909 0.00135803 13.6     1    1    1      Ah     50
# 2 (181025, 333558)     8.6     81  277 1141 6.983 0.01222430 14.0     1    1    1      Ah     30
# 3 (181165, 333537)     6.5     68  199  640 7.800 0.10302900 13.0     1    1    1      Ah    150
# 4 (181298, 333484)     2.6     81  116  257 7.655 0.19009400  8.0     1    2    0      Ga    270
# 5 (181307, 333330)     2.8     48  117  269 7.480 0.27709000  8.7     1    2    0      Ah    380
# 6 (181390, 333260)     3.0     61  137  281 7.791 0.36406700  7.8     1    2    0      Ga    470



# ======== NN by Radius ========
dnn <- RANN::nn2(coordinates(meuse), searchtype="radius", 
                 radius = 1000)$nn.idx 
head(dnn)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    1    2    3    8    7    4   13    5   14     9
# [2,]    2    1    3    8    7   13    4   14    9    84
# [3,]    3    1    2    4    7    8    5    9   84    14
# [4,]    4    3    5    7    6    1    2    8   10    84
# [5,]    5    6    7    4   10   84   11    3    9     8
# [6,]    6    5   10    4   11    7   84   30    9     3


# Add discrete values for random spot (50) and it's neighbors
meuse$nn_radius <- 0 # add zeros
meuse$nn_radius[50] <- 1 # spot = 1

neighbors_of_first_spot <- dnn[50, -1]
meuse$nn_radius[neighbors_of_first_spot] <- 2 # nn = 2
head(meuse)
# coordinates cadmium copper lead zinc  elev       dist   om ffreq soil lime landuse dist.m nn_radius
# 1 (181072, 333611)    11.7     85  299 1022 7.909 0.00135803 13.6     1    1    1      Ah     50         1
# 2 (181025, 333558)     8.6     81  277 1141 6.983 0.01222430 14.0     1    1    1      Ah     30         2
# 3 (181165, 333537)     6.5     68  199  640 7.800 0.10302900 13.0     1    1    1      Ah    150         2
# 4 (181298, 333484)     2.6     81  116  257 7.655 0.19009400  8.0     1    2    0      Ga    270         2
# 5 (181307, 333330)     2.8     48  117  269 7.480 0.27709000  8.7     1    2    0      Ah    380         2
# 6 (181390, 333260)     3.0     61  137  281 7.791 0.36406700  7.8     1    2    0      Ga    470         0

# visualize
pdf(here(plot_dir,"RANN_example_neighbors_radius.pdf"))
spplot(meuse, "nn_radius", xlab="Example Radius NN", col.regions=c("darkgrey", "red","blue") )
dev.off()



# ======== NN by Standard (exact neighbors) ========
dnn <- RANN::nn2(coordinates(meuse), searchtype="standard")$nn.idx 
head(dnn)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    1    2    3    8    7    4   13    5   14     9
# [2,]    2    1    3    8    7   13    4   14    9    84
# [3,]    3    1    2    4    7    8    5    9   84    14
# [4,]    4    3    5    7    6    1    2    8   10    84
# [5,]    5    6    7    4   10   84   11    3    9     8
# [6,]    6    5   10    4   11    7   84   30    9     3


# Add discrete values for random spot (50) and it's neighbors
meuse$nn_standard <- 0 # add zeros
meuse$nn_standard[50] <- 1 # spot = 1

neighbors_of_first_spot <- dnn[50, -1]
meuse$nn_standard[neighbors_of_first_spot] <- 2 # nn = 2
head(meuse)
# coordinates cadmium copper lead zinc  elev       dist   om ffreq soil lime landuse dist.m nn_radius
# 1 (181072, 333611)    11.7     85  299 1022 7.909 0.00135803 13.6     1    1    1      Ah     50         1
# 2 (181025, 333558)     8.6     81  277 1141 6.983 0.01222430 14.0     1    1    1      Ah     30         2
# 3 (181165, 333537)     6.5     68  199  640 7.800 0.10302900 13.0     1    1    1      Ah    150         2
# 4 (181298, 333484)     2.6     81  116  257 7.655 0.19009400  8.0     1    2    0      Ga    270         2
# 5 (181307, 333330)     2.8     48  117  269 7.480 0.27709000  8.7     1    2    0      Ah    380         2
# 6 (181390, 333260)     3.0     61  137  281 7.791 0.36406700  7.8     1    2    0      Ga    470         0

# visualize
pdf(here(plot_dir,"RANN_example_neighbors_standard.pdf"))
spplot(meuse, "nn_standard", xlab="Example Exact NN (k=10)", col.regions=c("darkgrey", "red","blue") )
dev.off()


# Notes: It seems like radius and standard (exact) give very similar results. I'm not sure really where this would vary,
#  but it's something to keep in mind for later. 


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

dnn <- RANN::nn2(spatialCoords(spe.subset), k=15, searchtype="standard")$nn.idx 
head(dnn)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    1 3170 1086 1339 2809 4044 3715 2274 3064  2557
# [2,]    2 1204 4115  454 2434 2692  551  658 1024  3599
# [3,]    3 3339 1708  788 2514 1856 2699 1717 1423  3607
# [4,]    4 2628 3856 2114  665 3515 2931 1133 2751  2590
# [5,]    5 2066 2202  510 1476 1377 3312  346  438  3878
# [6,]    6 1598 1613 2301 1493 4070  309  653 4130  2390

# Add discrete values for random spot (50) and it's neighbors
spe.subset$nn_standard <- 0 # add zeros
spe.subset$nn_standard[50] <- 1 # spot = 1

neighbors_of_first_spot <- dnn[50, -1]
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
pdf(here(plot_dir, 'SpotPlot_neighbors_standard_k15.pdf'))
vis_clus(
  spe = spe.subset,
  clustervar = "nn_standard",
  sampleid = unique(spe.subset$sample_id)[1]
)
dev.off()
