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
library(ggridges)

# Save directories
plot_dir = here("plots","single_cell_methods", "MAD")
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
id <- "151671"

p1 <- vis_clus(
  spe = spe,
  clustervar = "spatialLIBD",
  sampleid = id
)

p2 <- vis_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",
  sampleid = id
)


p3 <- vis_gene(
  spe = spe,
  geneid = "sum_gene",
  sampleid = id
)


pdf(height = 10, width=20, here(plot_dir, 'SpotPlots_dlPFC_typical_QC_metrics.pdf'))
(p1+p2+p3)
dev.off()



# ========== outlier detection using MAD ==============

qc.umi_mad <- isOutlier(spe$sum_umi, log=TRUE, type="lower")
table(qc.umi_mad)
# FALSE  TRUE 
# 47230   451 

qc.gene_mad <- isOutlier(spe$sum_gene, log=TRUE, type="lower")
table(qc.gene_mad)
# FALSE  TRUE
# 47171   510 

qc.mito_mad <- isOutlier(spe$expr_chrM_ratio, type="higher")
table(qc.mito_mad)
# qc.mito_mad
# FALSE  TRUE 
# 47598    83 

colData(spe)$qc_mito_mad <- qc.mito_mad
colData(spe)$qc_umi_mad <- qc.umi_mad
colData(spe)$qc_gene_mad <- qc.gene_mad

mad_outliers <- qc.umi_mad | qc.mito_mad | qc.gene_mad
spe$discard_mad <- mad_outliers

# get subset of first sample
spe.subset <- spe[,spe$sample_id == unique(spe$sample_id)[1]]

# plot using escheR
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "layer_guess_reordered") +
  scale_fill_manual(values = libd_layer_colors)

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "discard_mad") +
  scale_fill_manual(name = "",values = c("TRUE" = "red", "FALSE" = "grey")) 

pdf(here(plot_dir, 'Layers_vs_discarded_mad.pdf'), height = 7, width=14)
p1+p2
dev.off()

p1 <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = unique(spe$sample_id)[1]
)

p2 <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = unique(spe$sample_id)[5]
)

p3 <- vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = unique(spe$sample_id)[9]
)

p4 <- vis_clus(
  spe = spe,
  clustervar = "discard_mad",
  sampleid = unique(spe$sample_id)[1]
)

p5 <- vis_clus(
  spe = spe,
  clustervar = "discard_mad",
  sampleid = unique(spe$sample_id)[5]
)

p6 <- vis_clus(
  spe = spe,
  clustervar = "discard_mad",
  sampleid = unique(spe$sample_id)[9]
)


pdf(width=25, height=12,here(plot_dir, "SpotPlot_discard_mad.pdf"))
(p1+p2+p3)/(p4+p5+p6)
dev.off()



# ===== Average percentage discarded across samples =====
# Extract the required columns
discard_new_col <- spe$discard_mad
layer_guess_ordered_col <- spe$layer_guess_reordered
sample_id_col <- spe$sample_id
subject_col <- spe$subject
position_col <- spe$subject_position

# rewrite position_col such that if CONTAINS "pos0", make it "0um", if else "300um"
new_position_col <- ifelse(grepl("pos0", position_col), "0um", "300um")
unique(position_col)

# Create a data frame
df <- data.frame(discard_new = discard_new_col,
                 layer_guess_ordered = layer_guess_ordered_col,
                 sample_id = sample_id_col,
                 subject = subject_col,
                 position = new_position_col)

df <- df[!is.na(df$layer_guess_ordered), ]


summary_df <- df %>%
  group_by(sample_id, layer_guess_ordered, subject, position) %>%
  summarise(
    total_spots = n(),  # Total number of spots in each domain for each sample
    discarded_spots = sum(discard_new, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    .groups = 'keep'  # Dropping the grouping after summarise
  ) %>%
  mutate(
    percentage_discarded = (discarded_spots / total_spots) * 100  # Calculating the percentage
  ) 
head(summary_df)
# # A tibble: 6 × 7
# # Groups:   sample_id, layer_guess_ordered, subject, position [6]
#   sample_id layer_guess_ordered subject position total_spots discarded_spots
# <chr>     <fct>               <chr>   <chr>          <int>           <int>
# 1 151507    Layer1              Br5292  0um              817              29
# 2 151507    Layer2              Br5292  0um              305               0
# 3 151507    Layer3              Br5292  0um             1215              13
# 4 151507    Layer4              Br5292  0um              369              10
# 5 151507    Layer5              Br5292  0um              675               5
# 6 151507    Layer6              Br5292  0um              486               1


# Plot using ggplot2
p <- ggplot(summary_df, aes(x = layer_guess_ordered, y = percentage_discarded, fill = layer_guess_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "Percentage of Discarded Spots within each Spatial Domain",
       x = "Spatial Domain",
       y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(width=8, height=4,here(plot_dir, "Boxplot_discarded_by_layer.pdf"))
p
dev.off()



# stacked bar plot of discarded spot per sample by layer
pdf(width=8, height=4,here(plot_dir, "Barplot_discarded_by_layer.pdf"))
ggplot(summary_df, aes(x = sample_id, y = discarded_spots, fill = layer_guess_ordered)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "Total Number of Discarded Spots per Sample",
       x = "Sample ID",
       y = "Number of Spots") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~subject, 
             scales = "free_x",
             switch = "x",
             nrow=1)
dev.off()
