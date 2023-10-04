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
library(tidyr)

# Save directories
plot_dir = here("plots","single_cell_methods")
processed_dir = here("processed-data","single_cell_methods")

# ============== Let's plot some metrics using the final published data =========

spe <- fetch_data(type = "spe")
spe
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



# ================= Let's now look at how discarding spots based on standard single cell QC metrics effects layers =============



threshold_x = 500
threshold_y = 0.3
p <- ggplot(combined_data, aes(x = sum_umi, y = expr_ratio, color=layer)) + 
  geom_point(size = 0.5, alpha=0.2) + 
  ggtitle("QC metrics") + 
  theme_bw() + 
  scale_color_manual(name = "", values = spatialLIBD::libd_layer_colors) +
  geom_vline(xintercept = threshold_x, color = "black", linetype="dashed") +
  geom_hline(yintercept = threshold_y, color = "black", linetype="dashed") +
  geom_xsidedensity() + 
  geom_ysidedensity() +
  stat_ellipse()


pdf(here(plot_dir, "QCscatter_MitoPercent_vs_SumUmi_elipse.pdf"))
p
dev.off()

# ======== Let's discard based on these conservative threshold and see how it looks =======

qc_mito <- colData(spe)$expr_chrM_ratio > 0.3 # this is arbitrary
table(qc_mito)
# FALSE  TRUE 
# 47401   280 

qc_umi <- colData(spe)$sum_umi < 500 # this is arbitrary
table(qc_umi)
# FALSE  TRUE 
# 47372   309 

colData(spe)$qc_mito <- qc_mito
colData(spe)$qc_umi <- qc_umi

discard_new <- qc_umi | qc_mito
spe$discard_new <- discard_new


p1 <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = unique(spe$sample_id)[1]
)

p2 <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = unique(spe$sample_id)[2]
)

p3 <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = unique(spe$sample_id)[3]
)

p4 <- vis_clus(
  spe = spe,
  clustervar = "discard_new",
  sampleid = unique(spe$sample_id)[1]
)

p5 <- vis_clus(
  spe = spe,
  clustervar = "discard_new",
  sampleid = unique(spe$sample_id)[2]
)

p6 <- vis_clus(
  spe = spe,
  clustervar = "discard_new",
  sampleid = unique(spe$sample_id)[3]
)


pdf(width=20, height=10,here(plot_dir, "SpotPlot_discard_new.pdf"))
(p1+p2+p3)/(p4+p5+p6)
dev.off()


# =========== Now let's quantify the discarded spots based on layer ID =========

# ==== Total percentage =====
# Extract the required columns from the spe object
discard_new_col <- spe$discard_new
layer_guess_ordered_col <- spe$layer_guess_reordered

# Filter out the discarded spots
discarded_indices <- which(discard_new_col == TRUE)
discarded_domains <- layer_guess_ordered_col[discarded_indices]

# Count the occurrences of each spatial domain for the discarded spots
domain_counts <- table(discarded_domains)

# Convert to percentages
domain_percentages <- (domain_counts / sum(domain_counts)) * 100

# Convert the domain counts to a data frame
df <- as.data.frame(domain_percentages)
colnames(df) <- c("SpatialDomain", "Percentage")

# Plot using ggplot2
p <- ggplot(df, aes(x = SpatialDomain, y = Percentage, fill = SpatialDomain)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "Percentage of Discarded Spots by Spatial Domain",
       x = "Spatial Domain",
       y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(width=6, height=3,here(plot_dir, "Barplot_discarded_by_layer.pdf"))
p
dev.off()



# ===== Average percentage across samples =====
# Extract the required columns
discard_new_col <- spe$discard_new
layer_guess_ordered_col <- spe$layer_guess_reordered
sample_id_col <- spe$sample_id

# Create a data frame
df <- data.frame(discard_new = discard_new_col,
                 layer_guess_ordered = layer_guess_ordered_col,
                 sample_id = sample_id_col)

df <- df[!is.na(df$layer_guess_ordered), ]

# Filter out the discarded spots
discarded_df <- df[df$discard_new == TRUE, ]

# Group by sample_id and layer_guess_ordered, then calculate percentages
summary_df <- discarded_df %>%
  group_by(sample_id, layer_guess_ordered) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n / sum(n)) * 100)
  
# Create a dataframe with every possible combination of sample_id and layer_guess_ordered
all_combinations <- expand(df, sample_id, layer_guess_ordered)

# Left join summary_df onto all_combinations, filling in missing values with 0
final_df <- left_join(all_combinations, summary_df, by = c("sample_id", "layer_guess_ordered")) %>%
  replace_na(list(n = 0, percentage = 0))

# Plot using ggplot2
p <- ggplot(final_df, aes(x = layer_guess_ordered, y = percentage, fill = layer_guess_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "Percentage of Discarded Spots by Spatial Domain for each Sample ID",
       x = "Spatial Domain",
       y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(width=8, height=4,here(plot_dir, "Boxplot_discarded_by_layer.pdf"))
p
dev.off()

