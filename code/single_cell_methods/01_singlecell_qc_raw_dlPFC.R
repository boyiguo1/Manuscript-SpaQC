library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(scater)
library(ggspavis)
library(scran)
library(patchwork)
library(ggside)
library(ggpubr)
library(dplyr)

# Save directories
plot_dir = here("plots","single_cell_methods", "dlPFC_raw")
processed_dir = here("processed-data","single_cell_methods")
raw_dir = here("raw-data")

# ============== Let's plot some metrics using the final published data =========

load(here(raw_dir, "dlPFC_raw.Rdata"))
spe <- spe_raw
spe
# class: SpatialExperiment 
# dim: 36601 149757 
# metadata(0):
#   assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817 ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(149757): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
# colData names(27): sample_id in_tissue ... sample_id_complete count
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor
rm(spe_raw)

colnames(colData(spe))
# [1] "sample_id"              "in_tissue"              "array_row"              "array_col"              "10x_graphclust"        
# [6] "10x_kmeans_10_clusters" "10x_kmeans_2_clusters"  "10x_kmeans_3_clusters"  "10x_kmeans_4_clusters"  "10x_kmeans_5_clusters" 
# [11] "10x_kmeans_6_clusters"  "10x_kmeans_7_clusters"  "10x_kmeans_8_clusters"  "10x_kmeans_9_clusters"  "key"                   
# [16] "sum_umi"                "sum_gene"               "expr_chrM"              "expr_chrM_ratio"        "ManualAnnotation"      
# [21] "subject"                "region"                 "sex"                    "age"                    "diagnosis"             
# [26] "sample_id_complete"     "count"   



# =============== Plotting ============
# let's start by looking at the layer_guess_reordered annotations

p1 <- vis_clus(
  spe = spe,
  clustervar = "ManualAnnotation",
  sampleid = unique(spe$sample_id)[1]
)

p2 <- vis_gene(
  spe = spe,
  geneid = "expr_chrM",
  sampleid = unique(spe$sample_id)[1]
)

p3 <- vis_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",
  sampleid =unique(spe$sample_id)[1]
)

p4 <- vis_gene(
  spe = spe,
  geneid = "sum_umi",
  sampleid =unique(spe$sample_id)[1]
)

p5 <- vis_gene(
  spe = spe,
  geneid = "sum_gene",
  sampleid =unique(spe$sample_id)[1]
)

p6 <- vis_gene(
  spe = spe,
  geneid = "count",
  sampleid =unique(spe$sample_id)[1]
)

pdf(height = 10, width=20, here(plot_dir, 'SpotPlots_dlPFC_typical_QC_metrics.pdf'))
(p1+p2+p3)/(p4+p5+p6)
dev.off()