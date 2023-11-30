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

plot_dir = here('plots', 'spatial_methods', '02_outlier_detection','amygdala')
processed_dir = here('processed-data', 'spatial_methods')

load(here("raw-data","AMY_raw.Rdata"))
spe.amy <- spe
spe.amy
# class: SpatialExperiment 
# dim: 28412 39936 
# metadata(0):
#   assays(1): counts
# rownames(28412): ENSG00000243485 ENSG00000238009 ... ENSG00000278817 ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(39936): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTGTATTACACG-1
# TTGTTTGTGTAAATTC-1
# colData names(26): sample_id slide ... ManualAnnotation overlaps_tissue
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


load(here("raw-data","HYP_raw.Rdata"))
spe.hyp <- spe
spe.hyp
# class: SpatialExperiment 
# dim: 30256 39921 
# metadata(0):
#   assays(1): counts
# rownames(30256): ENSG00000243485 ENSG00000238009 ... ENSG00000278817 ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(39921): AAACAAGTATCTCCCA-1 AAACAATCTACTAGCA-1 ... TTGTTTGTATTACACG-1
# TTGTTTGTGTAAATTC-1
# colData names(28): sample_id in_tissue ... Pmask_dark_blue overlaps_tissue
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


colnames(colData(spe.amy))
# [1] "sample_id"              "slide"                  "array"                 
# [4] "brnum"                  "species"                "replicate"             
# [7] "in_tissue"              "array_row"              "array_col"             
# [10] "10x_graphclust"         "10x_kmeans_10_clusters" "10x_kmeans_2_clusters" 
# [13] "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters" 
# [16] "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters"  "10x_kmeans_8_clusters" 
# [19] "10x_kmeans_9_clusters"  "key"                    "sum_umi"               
# [22] "sum_gene"               "expr_chrM"              "expr_chrM_ratio"       
# [25] "ManualAnnotation"       "overlaps_tissue"  

colnames(colData(spe.hyp))
# [1] "sample_id"              "in_tissue"              "array_row"             
# [4] "array_col"              "10x_graphclust"         "10x_kmeans_10_clusters"
# [7] "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters" 
# [10] "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters" 
# [13] "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                   
# [16] "sum_umi"                "sum_gene"               "expr_chrM"             
# [19] "expr_chrM_ratio"        "ManualAnnotation"       "slide"                 
# [22] "array"                  "brnum"                  "species"               
# [25] "replicate"              "Nmask_dark_blue"        "Pmask_dark_blue"       
# [28] "overlaps_tissue" 

# ================ Visualize data first ============
# get only in tissue spots
spe.amy <- spe.amy[, colData(spe.amy)$in_tissue == TRUE]
spe.hyp <- spe.hyp[, colData(spe.hyp)$in_tissue == TRUE]

plot_order <- c()

vis_grid_gene(
  spe.amy,
  geneid = "expr_chrM",
  pdf_file = here(plot_dir, 'SpotPlots_AMY_mito.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)

vis_grid_gene(
  spe.amy,
  geneid = "expr_chrM_ratio",
  pdf_file = here(plot_dir, 'SpotPlots_AMY_mito_percent.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)

vis_grid_gene(
  spe.amy,
  geneid = "sum_umi",
  pdf_file = here(plot_dir, 'SpotPlots_AMY_sum_umi.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)

vis_grid_gene(
  spe.amy,
  geneid = "sum_gene",
  pdf_file = here(plot_dir, 'SpotPlots_AMY_sum_gene.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)


# ========= Hypo =========


vis_grid_gene(
  spe.hyp,
  geneid = "expr_chrM",
  pdf_file = here(plot_dir, 'SpotPlots_HYP_mito.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)

vis_grid_gene(
  spe.hyp,
  geneid = "expr_chrM_ratio",
  pdf_file = here(plot_dir, 'SpotPlots_HYP_mito_percent.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)

vis_grid_gene(
  spe.hyp,
  geneid = "sum_umi",
  pdf_file = here(plot_dir, 'SpotPlots_HYP_sum_umi.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)

vis_grid_gene(
  spe.hyp,
  geneid = "sum_gene",
  pdf_file = here(plot_dir, 'SpotPlots_HYP_sum_gene.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)


# ========= Basic qc using thresholds =======

# ===== amygdala =====
qc_mito <- colData(spe.amy)$expr_chrM_ratio > 0.3 # this is arbitrary
table(qc_mito)
# FALSE  TRUE 
# 29736   149 

qc_umi <- colData(spe.amy)$sum_umi < 250 # this is arbitrary
table(qc_umi)
# FALSE  TRUE 
# 28090  1795

qc_gene <- colData(spe.amy)$sum_gene < 50 # this is arbitrary
table(qc_gene)
# FALSE  TRUE 
# 29786    99 

colData(spe.amy)$qc_mito <- qc_mito
colData(spe.amy)$qc_umi <- qc_umi
colData(spe.amy)$qc_gene <- qc_gene

discard_new <- qc_umi | qc_mito | qc_gene
spe.amy$discard <- as.factor(discard_new)



# ===== Hypothalamus =====
qc_mito <- colData(spe.hyp)$expr_chrM_ratio > 0.3 # this is arbitrary
table(qc_mito)
# FALSE  TRUE 
# 33795  3495 

qc_umi <- colData(spe.hyp)$sum_umi < 250 # this is arbitrary
table(qc_umi)
# FALSE  TRUE 
# 37044   254

qc_gene <- colData(spe.hyp)$sum_gene < 50 # this is arbitrary
table(qc_gene)
# FALSE  TRUE 
# 37221    77 

colData(spe.hyp)$qc_mito <- qc_mito
colData(spe.hyp)$qc_umi <- qc_umi
colData(spe.hyp)$qc_gene <- qc_gene

discard_new <- qc_umi | qc_mito | qc_gene
spe.hyp$discard <- as.factor(discard_new)


# ======== Visualize discarded spots =======

vis_grid_clus(
  spe.amy,
  clustervar="discard",
  pdf_file = here(plot_dir, 'SpotPlots_AMY_discarded.pdf'),
  height = 24,
  width = 36,
  image_id = "lowres",
  point_size = 2.5,
)

vis_grid_clus(
  spe.hyp,
  clustervar="discard",
  pdf_file = here(plot_dir, 'SpotPlots_HYP_discarded.pdf'),
  height = 24,
  width = 36,
  image_id = "lowres",
  point_size = 2.5,
)


# =========== Spatially-aware QC with SpotSweeper ==============
# Notes for this section: Currently localOutlier only allows one feature input.
# should change this to accept a list of features that defaults in sum_umi, sum_gene, and expr_chrM_ratio.
#
# BiocNeighbors also throws a lot of errors. May be better to go with a 
# distance/radius based metric instead of nearest K since there will also be 
# equidistant neighbors. 

# ======== Amygdala ========

colnames(colData(spe.amy))
# [1] "sample_id"              "slide"                  "array"                  "brnum"                 
# [5] "species"                "replicate"              "in_tissue"              "array_row"             
# [9] "array_col"              "10x_graphclust"         "10x_kmeans_10_clusters" "10x_kmeans_2_clusters" 
# [13] "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters" 
# [17] "10x_kmeans_7_clusters"  "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                   
# [21] "sum_umi"                "sum_gene"               "expr_chrM"              "expr_chrM_ratio"       
# [25] "ManualAnnotation"       "overlaps_tissue"        "qc_mito"                "qc_umi"                
# [29] "qc_gene"                "discard_new"            "discard"  

spe.amy <- localOutliers(spe.amy,
                         features=c("sum_umi","sum_gene", "expr_chrM_ratio"), 
                         n_neighbors=36, 
                         data_output=TRUE,
                         method="multivariate",
                         cutoff=3)

colnames(colData(spe.amy))
# [1] "sample_id"              "slide"                  "array"                  "brnum"                 
# [5] "species"                "replicate"              "in_tissue"              "array_row"             
# [9] "array_col"              "10x_graphclust"         "10x_kmeans_10_clusters" "10x_kmeans_2_clusters" 
# [13] "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters"  "10x_kmeans_6_clusters" 
# [17] "10x_kmeans_7_clusters"  "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                   
# [21] "sum_umi"                "sum_gene"               "expr_chrM"              "expr_chrM_ratio"       
# [25] "ManualAnnotation"       "overlaps_tissue"        "qc_mito"                "qc_umi"                
# [29] "qc_gene"                "discard_new"            "discard"                "sum_umi_log2"         
# [33] "coords"                 "local_outliers"         "sum_umi_z"    

spe.hyp <- localOutliers(spe.hyp,
                         features=c("sum_umi","sum_gene", "expr_chrM_ratio"), 
                         n_neighbors=36, 
                         data_output=TRUE,
                         method="multivariate",
                         cutoff=3)


# plot new outliers 
spe.amy$local_outliers <- as.factor(spe.amy$local_outliers)
vis_grid_clus(
  spe.amy,
  clustervar="local_outliers",
  pdf_file = here(plot_dir, 'SpotPlots_AMY_local_outliers_mv.pdf'),
  height = 24,
  width = 36,
  image_id = "lowres",
  point_size = 2.5,
)

#spe.amy$sum_umi_log2 <- log2(spe.amy$sum_umi)

vis_grid_gene(
  spe.amy,
  geneid = "sum_umi_log2",
  pdf_file = here(plot_dir, 'SpotPlots_AMY_sum_umi_log2.pdf'),
  assayname = "counts",
  minCount = 0,
  height = 24,
  width = 36,
  point_size = 2.5,
)

vis_grid_gene(
  spe.amy,
  geneid = "sum_umi_z",
  pdf_file = here(plot_dir, 'SpotPlots_AMY_sum_umi_z.pdf'),
  assayname = "counts",
  minCount = -10,
  height = 24,
  width = 36,
  point_size = 2.5,
)

# ======== Hypothalamus ========
spe.hyp$local_outliers <- as.factor(spe.hyp$local_outliers)
vis_grid_clus(
  spe.hyp,
  clustervar="local_outliers",
  pdf_file = here(plot_dir, 'SpotPlots_HYP_local_outliers_mv.pdf'),
  height = 24,
  width = 36,
  image_id = "lowres",
  point_size = 2.5,
)