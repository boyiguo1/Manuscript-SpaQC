devtools::install("/Users/mtotty2/Documents/R/SpotSweeper")
library(SpotSweeper)

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
library(dplyr)
library(tidyr)
library(scran)
library(scater)

plot_dir = here('plots', 'spatial_methods', '03_mito_variance')
processed_dir = here('processed-data')

## Download the spot-level data
spe <- fetch_data(type = "spe")
spe


# load SCZ data
spe.scz <- readRDS(here("raw-data", "test_spe_after_spot_qc_63.rds"))
spe.scz  


# ======= Variance of Maynard et al data ========

spe <- localOutliers(spe, 
                     features=c("sum_umi","sum_gene", "expr_chrM_ratio"), 
                     n_neighbors=36, 
                     data_output=TRUE,
                     method="multivariate",
                     cutoff=3)

spe <- localVariance(spe,
                     n_neighbors=36,
                     n_cores=6)

# plot variance of Maynard et al data
vis_grid_gene(spe, 
              geneid="expr_chrM_ratio_var", 
              spatial=FALSE,
              pdf_file=here(plot_dir, "chrM_var_new.pdf"),
              point_size =2.5) 

# ======= Variance of SCZ PNN data ========

# plot mito ratio
vis_grid_gene(spe.scz, 
              geneid="expr_chrM_ratio", 
              spatial=FALSE,
              pdf_file=here(plot_dir, "SCZ_PNN_mito_ratio.pdf")
              ) 

spe.scz <- localOutliers(spe.scz, 
                     features=c("sum_umi","sum_gene", "expr_chrM_ratio"), 
                     n_neighbors=36, 
                     data_output=TRUE,
                     method="multivariate",
                     cutoff=3,
                     log2=FALSE)

vis_grid_gene(spe.scz, 
              geneid="expr_chrM_ratio_var", 
              spatial=FALSE,
              pdf_file=here(plot_dir, "SCZ_PNN_mito_ratio_var.pdf"),
              point_size =.85
              ) 

pdf(here(plot_dir, "SCZ_PNN_mito_ratio_var_hist.pdf"))
hist(spe.scz$expr_chrM_ratio_var, breaks=100)
dev.off()

# subset to sample V13M06−279_A1
spe.scz.subset <- subset(spe.scz, ,sample_id == unique(spe.scz$sample_id)[12])

# plot histogram
pdf(here(plot_dir, "SCZ_PNN_mito_ratio_var_hist_V13M06−279_A1.pdf"))
hist(spe.scz.subset$expr_chrM_ratio_var, breaks=100)
dev.off()


# subset to sample V13M06−279_A1
spe.scz.subset <- subset(spe.scz, ,sample_id == unique(spe.scz$sample_id)[1])

# plot histogram
pdf(here(plot_dir, "SCZ_PNN_mito_ratio_var_hist_V12M14−053_A1.pdf"))
hist(spe.scz.subset$expr_chrM_ratio_var, breaks=100)
dev.off()



# save data
saveRDS(spe.scz, here(processed_dir, "spe_scz_local_outliers.rds"))



# ====== Trying various K values ======

# subset to sample V13M06−279_A1
spe.scz.subset <- subset(spe.scz, ,sample_id == unique(spe.scz$sample_id)[12])


spe.scz.subset <- localVariance(spe.scz.subset,
                     n_neighbors=36,
                     n_cores=8,
                     name="mito_var_k36")


spe.scz.subset <- localVariance(spe.scz.subset,
                         n_neighbors=18,
                         n_cores=8,
                         name="mito_var_k18")


spe.scz.subset <- localVariance(spe.scz.subset,
                         n_neighbors=6,
                         n_cores=8,
                         name="mito_var_k6")

# get column data from spe.scz and make into dataframe for the various Ks
df <- data.frame(
  k36 = spe.scz.subset$mito_var_k36,
  k18 = spe.scz.subset$mito_var_k18,
  k6 = spe.scz.subset$mito_var_k6
)

# plot histogram using ggplot
p1 <- ggplot(df, aes(x=k36)) + 
  geom_histogram(bins=250) +
  labs(title="Histogram of mito variance for various K=36",
       x="mito variance",
       y="count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")


p2 <- ggplot(df, aes(x=k18)) +
  geom_histogram(bins=250) +
  labs(title="Histogram of mito variance for various K=18",
       x="mito variance",
       y="count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")


p3 <- ggplot(df, aes(x=k6)) +
  geom_histogram(bins=250) +
  labs(title="Histogram of mito variance for various K=6",
       x="mito variance",
       y="count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

pdf(here(plot_dir, "SCZ_PNN_mito_ratio_var_V13M06−279_A1_variousKs_hist.pdf"), height=5, width=5)
p1/p2/p3
dev.off()


p1 <- vis_gene(spe.scz.subset, 
         geneid="mito_var_k36", 
         spatial=FALSE,
         point_size=1
)


p2 <- vis_gene(spe.scz.subset, 
         geneid="mito_var_k18", 
         spatial=FALSE,
         point_size=1
)


p3 <- vis_gene(spe.scz.subset, 
         geneid="mito_var_k6", 
         spatial=FALSE,
         point_size=1
)

pdf(here(plot_dir, "SCZ_PNN_mito_ratio_var_V13M06−279_A1_variousKs_spotplots.pdf"), height=7, width=13)
(p1+p2+p3)
dev.off()


# ========= Let's try to label some of these "outliers" ======
pdf(here(plot_dir, "SCZ_PNN_mito_ratio_var_V13M06−279_A1_k18_outlier_threhsold.pdf"), height=7, width=13)
ggplot(df, aes(x=k18)) +
  geom_histogram(bins=250) +
  labs(title="Histogram of mito variance for various K=18",
       x="mito variance",
       y="count") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0.0002, linetype="dashed", color="red") 
dev.off()

# this looks good. Now add a new colData for TRUE is < .0002
spe.scz.subset$mito_var_k18_outlier <- spe.scz.subset$mito_var_k18 < .0002

# plot using escheR
p1 <- vis_gene(spe.scz.subset, 
         geneid="expr_chrM_ratio", 
         spatial=FALSE,
         point_size=1)

p2 <- vis_gene(spe.scz.subset, 
               geneid="mito_var_k18", 
               spatial=FALSE,
               point_size=1)

p3 <- vis_clus(spe.scz.subset, 
         clustervar="mito_var_k18_outlier", 
         spatial=FALSE,
         point_size=1)

pdf(here(plot_dir, "SCZ_PNN_mito_ratio_var_V13M06−279_A1_k18_outliers.pdf"), height=4, width=12)
p1+p2+p3
dev.off()


pdf(here(plot_dir, "SCZ_PNN_mito_ratio_var_V13M06−279_A1_k18_outliers_scatter.pdf"), height=5, width=5)
plot(spe.scz.subset$mito_var_k18, spe.scz.subset$expr_chrM_ratio, pch=20, col=ifelse(spe.scz.subset$mito_var_k18_outlier, "red", "black"))
dev.off()

# add expr_chrM_mito to tibble and plot density scatter plot
library(ggpointdensity)
library(viridis)

dat <- tibble(ratio = spe.scz.subset$expr_chrM_ratio,
       var = spe.scz.subset$mito_var_k18,
       group = "mito_percent")

p1 <- ggplot(dat, aes(x=var, y=ratio)) + 
  geom_bin2d(bins=80) +
  scale_fill_continuous(type = "viridis")
