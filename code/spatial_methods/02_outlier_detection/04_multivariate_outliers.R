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


plot_dir = here('plots', 'spatial_methods', '02_outlier_detection',"multivariate_methods")

## Download the spot-level data
spe <- fetch_data(type = "spe")

## This is a SpatialExperiment object
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



# ========= Playing around with multivariate methods ==========
# first thought is to true Sum UMI, unique UMI, and Mito% as features

spe.subset <- subset(spe, ,sample_id == unique(spe$sample_id[1]))

# ===================== Spatially aware QC ==============================
colData(spe.subset)$sum_umi_log2 <- log2(spe.subset$sum_umi)
colData(spe.subset)$sum_gene_log2 <- log2(spe.subset$sum_gene)


# create a list of spatial coordinates and qc features
spaQC <- colData(spe.subset)
spaQC$coords <- spatialCoords(spe.subset)

# Check the first few entries of your main data objects
head(spaQC)
head(spatialCoords(spe.subset))



# Find nearest neighbors
dnn <- findKNN(spatialCoords(spe.subset), k=36)$index
head(dnn)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15]
# [1,] 3170 1086 1339 2809 4044 3715 2274 3064 2557  3389    83  2602  2079   248  1275
# [2,] 1204 4115  454 2434  551 2692  658 1024 3599   786  4182  1665  2339  3540  1438
# [3,] 3339 1708  788 2514 1856 2699 1717 1423 3607  2154   671  3920  3090   113  3689
# [4,] 2628 3856 2114  665 3515 2931 1133 2751 2590  1975  3196  3400   570  1435  3610
# [5,] 2066 2202  510 1476 1377 3312  346  438 3878  1734  2773  2141  4211  3897  1933
# [6,] 1598 1613 2301 1493 4070  309  653 2390 4130  2443  3601  3813   466  3062   764

# NOTE: findKNN does NOT include the center spot in this output. Other methods do (RANN).


# =========== Modified Z based on Mito percent ==========
var.mito<- rep(NA,nrow(spaQC))
z.mito <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.mito[i] <- var( spaQC[c(i, dnn.idx[dnn.idx != 0]),]$expr_chrM_ratio, na.rm=TRUE)
  z.mito[i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$expr_chrM_ratio)[1]
}
z.mito[!is.finite(z.mito)] <- 0 


# ===========  Modified Z based on Sum UMI percent ==========
var.umi<- rep(NA,nrow(spaQC))
z.umi <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.umi[i] <- var( spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_umi_log2, na.rm=TRUE)
  z.umi[i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_umi_log2)[1]
}
z.umi[!is.finite(z.umi)] <- 0 


# ===========  Modified Z based on Sum UMI percent ==========
var.gene<- rep(NA,nrow(spaQC))
z.gene <- rep(NA,nrow(spaQC))  
for(i in 1:nrow(dnn)){
  dnn.idx <- dnn[i,] 
  var.gene[i] <- var( spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_gene_log2, na.rm=TRUE)
  z.gene[i] <- outliers(spaQC[c(i, dnn.idx[dnn.idx != 0]),]$sum_gene_log2)[1]
}
z.gene[!is.finite(z.gene)] <- 0 


# Add back to spe
spe.subset$z.mito <- z.mito
spe.subset$z.umi <- z.umi
spe.subset$z.gene <- z.gene



plot(spe.subset$z.mito, spe.subset$z.umi)
plot(spe.subset$z.mito, spe.subset$z.gene)
plot(spe.subset$z.umi, spe.subset$z.gene)



# ====== mvoutlier package =====

df <- data.frame(
  mito = spe.subset$z.mito,
  umi = spe.subset$z.umi,
  gene = spe.subset$z.gene
)

distances <- dd.plot(df, quan=1/2, alpha=0.00001)

df <- data.frame(
  mito = spe.subset$z.mito,
  umi = spe.subset$z.umi
)

colorData <- color.plot(df, quan=1, alpha=0.5)

pdf(here(plot_dir, "mvoutliers_output.pdf"))
outlier_spots <- aq.plot(df, delta=qchisq(0.99, df = ncol(df)),
                    quan = 1, alpha = 0.00001)
dev.off()



spe.subset$outliers <- outlier_spots[[1]]

# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi") |>
  add_ground(var = "outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")




pdf(here(plot_dir, 'SpotPlots_mvoutliers.pdf'))
p1
dev.off()


ggplot(as.data.frame(distances), aes(x=md.rob, y=md.cla, color=outliers)) +
  geom_point()




# ========= Local outlier factor using DBCAN ==========

library(dbscan)
df <- data.frame(
  mito = spe.subset$z.mito,
  umi = spe.subset$z.umi,
  gene = spe.subset$z.gene
)

outs <- lof(df, minPts=20)
outs_scaled <- lof(scale(df), minPts=20)

df$lof_scaled <- outs_scaled
df$lof <- outs


summary(outs_scaled)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9089  0.9945  1.0210  1.0525  1.0671  4.3245 

summary(outs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9075  0.9944  1.0211  1.0523  1.0670  4.3158 

hist(outs_scaled, breaks = 100)
hist(outs, breaks = 100)

scatter.smooth(outs,outs_scaled)

df$outliers <- ifelse(df$lof_scaled > 3, TRUE, FALSE)

p1 <- ggplot(df, aes(x=z.mito, y=z.umi, color=outliers)) +
  geom_point()

p2 <- ggplot(df, aes(x=z.mito, y=z.gene, color=outliers)) +
  geom_point()


p3 <- ggplot(df, aes(x=z.umi, y=z.gene, color=outliers)) +
  geom_point()

pdf(width=20, height=6, here(plot_dir,"LOF_outlier_scattplots.pdf"))
p1+p2+p3
dev.off()

# ======== SpotPlots of MV outliers =========


spe.subset$mvlocal_outliers <- df$outliers 

# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p3 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") |>
  add_ground(var = "mvlocal_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p4 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") |>
  add_ground(var = "mvlocal_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

pdf(width = 12.5, height = 12.5, here(plot_dir, 'SpotPlots_LOF_k20_LOF3.pdf'))
(p1+p2)/(p3+p4)
dev.off()



