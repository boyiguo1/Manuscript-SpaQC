library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(scater)
library(ggspavis)
library(scran)
library(patchwork)

# Save directories
plot_dir = here("plots","single_cell_methods")
processed_dir = here("processed-data","single_cell_methods")

# ============== Let's plot some metrics using the final published data =========

spe <- fetch_data(type = "spe")
# class: SpatialExperiment 
# dim: 33538 47681 
# metadata(0):
#   assays(2): counts logcounts
# rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475 ENSG00000268674
# rowData names(9): source type ... gene_search is_top_hvg
# colnames(47681): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ... TTGTTTCCATACAACT-1 TTGTTTGTGTAAATTC-1
# colData names(69): sample_id Cluster ... array_row array_col
# reducedDimNames(6): PCA TSNE_perplexity50 ... TSNE_perplexity80 UMAP_neighbors15
# mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


colnames(colData(spe))
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
# [69] "array_col"   


# =============== Plotting ============
# let's start by looking at the layer_guess_reordered annotations

p1 <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
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
  geneid = "cell_count",
  sampleid =unique(spe$sample_id)[1]
)

pdf(height = 10, width=20, here(plot_dir, 'SpotPlots_dlPFC_typical_QC_metrics.pdf'))
(p1+p2+p3)/(p4+p5+p6)
dev.off()

# =========== plotting average single cell QC metrics per layer ===========

# Extract data from colData of the SPE object
expr_data <- colData(spe)$expr_chrM
expr_ratio <- colData(spe)$expr_chrM_ratio
layer_data <- colData(spe)$layer_guess_reordered
cell_count <- colData(spe)$cell_count
sum_umi <- colData(spe)$sum_umi

# Check that the lengths of the extracted data are the same
if(length(expr_data) != length(layer_data)) {
  stop("Mismatch in the lengths of expr_data and layer_data!")
}

# Combine the expression data with the layer data
combined_data <- data.frame(
  expr_chrM = expr_data,
  expr_ratio = expr_ratio,
  layer = layer_data,
  sum_umi = sum_umi,
  cell_count = cell_count
)

# Create the boxplot
p1 <- ggplot(combined_data, aes(x = layer, y = expr_chrM, fill = layer)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.shape = NA) + # Boxplot without outliers
  theme(legend.position = "none") + # Remove legend
  labs(title = "Expression of expr_chrM per layer",
       x = "Layer",
       y = "Expression level")

# Create the boxplot
p2 <- ggplot(combined_data, aes(x = layer, y = sum_umi, fill = layer)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.shape = NA) + # Boxplot without outliers
  theme(legend.position = "none") + # Remove legend
  labs(title = "Library size",
       x = "Layer",
       y = "Sum UMI")

# Create the boxplot
p3 <- ggplot(combined_data, aes(x = layer, y = cell_count, fill = layer)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.shape = NA) + # Boxplot without outliers
  theme(legend.position = "none") + # Remove legend
  labs(title = "Cell counts",
       x = "Layer",
       y = "Cells")

# Create the boxplot
p4 <- ggplot(combined_data, aes(x = layer, y = expr_ratio, fill = layer)) +
  geom_boxplot(alpha = 0.8, width = 0.7, outlier.shape = NA) + # Boxplot without outliers
  theme(legend.position = "none") + # Remove legend
  labs(title = "Percent Mito",
       x = "Layer",
       y = "Percent mito (%)")

pdf(height = 15, width=7.5, here(plot_dir, 'Boxplots_dlPFC_typical_QC_metrics.pdf'))
p1/p2/p3/p4
dev.off()


# ============= Let's quantify some more typical QC metircs ===========

# histogram of library sizes
pdf(here(plot_dir, "Histogram_dlPFC_LibrarySize_.pdf"))
hist(colData(spe)$sum_umi, breaks = 80, xlim=c(0,15000))
dev.off()

# histogram of mito %
pdf(here(plot_dir, "Histogram_dlPFC_MitoPercent_.pdf"))
hist(colData(spe)$expr_chrM_ratio, breaks = 20)
dev.off()

#### Both library size and mito percent are relatively uniform with very few spots >30% mito. But is this equal for all layers?


p1 <- ggplot(combined_data, aes(x=expr_ratio, fill=layer)) +
  geom_histogram(alpha=0.5, position="identity", bins=30) +
  theme_minimal() +
  labs(title="Histogram of mito percent by Layer", x="Mito Percent", y="Count", fill="Layer")

p2 <- ggplot(combined_data, aes(x=expr_ratio, fill=layer)) +
  geom_density(alpha=0.5, adjust=1) +
  theme_minimal() +
  labs(title="Density of mito percent by Layer", x="Mito Percent", y="Density", fill="Layer")

pdf(height = 5, width=10, here(plot_dir, 'Histogram_dlPFC_mito_by_layer.pdf'))
(p1+p2)
dev.off()


p1 <- ggplot(combined_data, aes(x=sum_umi, fill=layer)) +
  geom_histogram(alpha=0.5, position="identity", bins=30) +
  theme_minimal() +
  labs(title="Histogram of mito percent by Layer", x="Sum UMI", y="Count", fill="Layer")

p2 <- ggplot(combined_data, aes(x=sum_umi, fill=layer)) +
  geom_density(alpha=0.5, adjust=1) +
  theme_minimal() +
  labs(title="Density of library size by Layer", x="Sum UMI", y="Density", fill="Layer")

pdf(height = 5, width=10, here(plot_dir, 'Histogram_dlPFC_sum_by_layer.pdf'))
(p1+p2)
dev.off()