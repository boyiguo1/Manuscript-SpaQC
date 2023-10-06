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

# Find nearest neighbors
dnn <- RANN::nn2(spatialCoords(spe.subset), k=15, searchtype="standard")$nn.idx 
head(dnn)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# [1,]    1 3170 1086 1339 2809 4044 3715 2274 3064  2557
# [2,]    2 1204 4115  454 2434 2692  551  658 1024  3599
# [3,]    3 3339 1708  788 2514 1856 2699 1717 1423  3607
# [4,]    4 2628 3856 2114  665 3515 2931 1133 2751  2590
# [5,]    5 2066 2202  510 1476 1377 3312  346  438  3878
# [6,]    6 1598 1613 2301 1493 4070  309  653 4130  2390

