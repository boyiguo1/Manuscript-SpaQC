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

# plot variance of Maynard et al data
vis_grid_gene(spe, 
              geneid="expr_chrM_ratio_var", 
              spatial=FALSE,
              pdf_file=here(plot_dir, "chrM_var.pdf"),
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
