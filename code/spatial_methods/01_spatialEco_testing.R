# This was taken from https://gis.stackexchange.com/questions/219255/spatial-outliers-detection-in-r

library(sp)
library(spdep)
library(RANN)
library(spatialEco)

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



# ========= 
# This definitely looks useful. RANN::nn2 appears similar to BiocNeighbors that Boyi
# shared. It gives (IMO) more intuitive neighbor definitions, such as neighbores based on
# radius; although radius may not be the best metric. 
#
# TODO: it would be worth seeing if RANN:nn2 is faster than BiocNeighbors. 
# TODO: Generate some visuals for calculated neighbors (similar to above "Hot Spots")





