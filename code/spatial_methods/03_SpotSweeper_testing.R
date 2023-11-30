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

plot_dir = here('plots', 'spatial_methods', '02_outlier_detection', 'SpotSweeper')
processed_dir = here('processed-data')

## Download the spot-level data
spe <- fetch_data(type = "spe")
spe


# ====================== Testing with my newly made function ======================

# I turned the above code into functions in the private SpotSweeper package
# https://github.com/MicTott/SpotSweeper

spe <- localOutliers(spe, k= 6, feature='sum_umi', samples='sample_id', log2=TRUE, z_threshold=3, output_z=TRUE)



# ====================== Explore these outliers ======================

spe.subset <- subset(spe, ,sample_id == spe$sample_id[1])
# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_z") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p3 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") |>
  add_ground(var = "local_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p4 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_z") |>
  add_ground(var = "local_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

pdf(width = 12.5, height = 12.5, here(plot_dir, 'Many_SpotPlots_sumUmi_k6_test.pdf'))
(p1+p2)/(p3+p4)
dev.off()

# ====================== Explore these outliers ======================

# ===== Average percentage across samples =====
# Extract the required columns
discard_new_col <- spe$z_umi_outlier
layer_guess_ordered_col <- spe$layer_guess_reordered
sample_id_col <- spe$sample_id

# Create a data frame
df <- data.frame(discard_new = discard_new_col,
                 layer_guess_ordered = layer_guess_ordered_col,
                 sample_id = sample_id_col)

df <- df[!is.na(df$layer_guess_ordered), ]


summary_df <- df %>%
  group_by(sample_id, layer_guess_ordered) %>%
  summarise(
    total_spots = n(),  # Total number of spots in each domain for each sample
    discarded_spots = sum(discard_new, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    .groups = 'keep'  # Dropping the grouping after summarise
  ) %>%
  mutate(
    percentage_discarded = (discarded_spots / total_spots) * 100  # Calculating the percentage
  ) 

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

pdf(width=8, height=4,here(plot_dir, "Boxplot_percent_discarded_k6.pdf"))
p
dev.off()



# ============ Testing different Ks ==========


spe <- localOutliers(spe, k= 18, threshold=3)


# ====================== Explore these outliers ======================

spe.subset <- subset(spe, ,sample_id == spe$sample_id[1])
# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p3 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") |>
  add_ground(var = "z_umi_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p4 <- make_escheR(spe.subset) |> 
  add_fill(var = "z.umi") |>
  add_ground(var = "z_umi_outlier", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

pdf(width = 12.5, height = 12.5, here(plot_dir, 'Many_SpotPlots_sumUmi_k18.pdf'))
(p1+p2)/(p3+p4)
dev.off()

# ====================== Explore these outliers ======================


# ==== Total percentage =====
# Extract the required columns
discard_new_col <- spe$z_umi_outlier
layer_guess_ordered_col <- spe$layer_guess_reordered
sample_id_col <- spe$sample_id

# Create a data frame
df <- data.frame(discard_new = discard_new_col,
                 layer_guess_ordered = layer_guess_ordered_col,
                 sample_id = sample_id_col)

df <- df[!is.na(df$layer_guess_ordered), ]


summary_df <- df %>%
  group_by(sample_id, layer_guess_ordered) %>%
  summarise(
    total_spots = n(),  # Total number of spots in each domain for each sample
    discarded_spots = sum(discard_new, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    .groups = 'keep'  # Dropping the grouping after summarise
  ) %>%
  mutate(
    percentage_discarded = (discarded_spots / total_spots) * 100  # Calculating the percentage
  ) 

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

pdf(width=8, height=4,here(plot_dir, "Boxplot_percent_discarded_k18.pdf"))
p
dev.off()




# ============ Testing different Ks ==========

spe <- localOutliers(spe, 
                     features=c("sum_umi","sum_gene", "expr_chrM_ratio"), 
                     n_neighbors=36, 
                     data_output=TRUE,
                     method="multivariate",
                     cutoff=3)


# ====================== Explore these outliers =====================

spe.subset <- subset(spe, ,sample_id == spe$sample_id[1])

# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p2 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_z") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

p3 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") |>
  add_ground(var = "local_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient(low ="white",high =  "darkgreen")

p4 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_z") |>
  add_ground(var = "local_outliers", stroke = 1) +
  scale_color_manual(
    name = "", # turn off legend name for ground_truth
    values = c(
      "TRUE" = "red",
      "FALSE" = "transparent")
  ) +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

pdf(width = 12.5, height = 12.5, here(plot_dir, 'Many_SpotPlots_sumUmi_sumGene_k36_mv.pdf'))
(p1+p2)/(p3+p4)
dev.off()

# ====================== Explore these outliers ======================

# ==== Total percentage =====
# Extract the required columns
discard_new_col <- as.matrix(spe$local_outliers)
layer_guess_ordered_col <- spe$layer_guess_reordered
sample_id_col <- spe$sample_id

# Create a data frame
df <- data.frame(discard_new = discard_new_col,
                 layer_guess_ordered = layer_guess_ordered_col,
                 sample_id = sample_id_col)

df <- df[!is.na(df$layer_guess_ordered), ]


summary_df <- df %>%
  group_by(sample_id, layer_guess_ordered) %>%
  summarise(
    total_spots = n(),  # Total number of spots in each domain for each sample
    discarded_spots = sum(discard_new, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    .groups = 'keep'  # Dropping the grouping after summarise
  ) %>%
  mutate(
    percentage_discarded = (discarded_spots / total_spots) * 100  # Calculating the percentage
  ) 

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

pdf(width=8, height=4,here(plot_dir, "Boxplot_percent_discarded_k36_retest_mv.pdf"))
p
dev.off()




# ================= UMAPS ===================

# =========== Threshold outliers =========

# detect outliers based on SC methods
qc_mito <- colData(spe)$expr_chrM_ratio > 0.2 # this is arbitrary
table(qc_mito)
# FALSE  TRUE 
# 47401   280 

qc_umi <- colData(spe)$sum_umi < 5000 # this is arbitrary
table(qc_umi)
# FALSE  TRUE 
# 47372   309 

qc_gene <- colData(spe)$sum_gene < 1000 # this is arbitrary
table(qc_gene)
# FALSE  TRUE 
# 46699   982 

colData(spe)$qc_mito <- qc_mito
colData(spe)$qc_umi <- qc_umi
colData(spe)$qc_gene <- qc_gene

sc_outliers <- qc_umi | qc_mito | qc_gene
spe$discard_new <- sc_outliers


reducedDimNames(spe)
# [1] "PCA"               "TSNE_perplexity50" "TSNE_perplexity5"  "TSNE_perplexity20" "TSNE_perplexity80"
# [6] "UMAP_neighbors15" 


# UMAP
point_size <- 0.3
# ===== row 1 - expr_chrM_ratio =====
p1 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="layer_guess_reordered", point_size = point_size) +
  ggtitle("Percent Mitochondrial Genes (> 20%)") +
  theme(plot.title = element_text(size = 20, face = "bold"))

p2 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="expr_chrM_ratio", point_size = point_size) 
#+ scale_color_gradient(low = "white", high = "red")

p3 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="qc_mito", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red"))

# ===== row 2 - sum_umi =====
p4 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="layer_guess_reordered", point_size = point_size)+
  ggtitle("Sum UMI (< 1000)") +
  theme(plot.title = element_text(size = 20, face = "bold"))

p5 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="sum_umi", point_size = point_size) 
#+ scale_color_gradient(low = "white", high = "red")

p6 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="qc_umi", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red"))


# ==== row 3 - sum_gene =====
p7 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="layer_guess_reordered", point_size = point_size) +
  ggtitle("Sum Gene (< 1000)") +
  theme(plot.title = element_text(size = 20, face = "bold"))

p8 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="sum_gene", point_size = point_size) 
#+ scale_color_gradient(low = "white", high = "red")

p9 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="qc_gene", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red"))


# ==== row 4 - local outliers
p10 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="layer_guess_reordered", point_size = point_size) +
  ggtitle("Local outliers (SpotSweeper; > 3 Zs) ") +
  theme(plot.title = element_text(size = 20, face = "bold"))

p11 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="LOF", point_size = point_size) 
#+ scale_color_gradient(low = "white", high = "red")

p12 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="local_outliers", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red"))

pdf(width=16, height=12,here(plot_dir, "UMAPs.pdf"))
(p1+p2+p3)/(p4+p5+p6)/(p7+p8+p9)/(p10+p11+p12)
dev.off()

png(width=6000, height=5000,here(plot_dir, "UMAPs.png"), res=300)
(p1+p2+p3)/(p4+p5+p6)/(p7+p8+p9)/(p10+p11+p12)
dev.off()



# ========= Replicating Fig 5 =========
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009290

spe$clusters <- spe$layer_guess_reordered

p1 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="clusters", point_size = point_size) +
  ggtitle("Manually annotated cluster assignment") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p2 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="sum_gene", point_size = point_size) +
  ggtitle("Unique genes per cell") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p3 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="qc_gene", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red")) +
  ggtitle("Cells excluded by <1000 threshold") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p4 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                      colour_by="local_outliers", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red")) +
  ggtitle("Cells excluded by detecting local outliers (SpotSweeper)") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p5 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="sample_id", point_size = point_size) + 
  ggtitle("UMAP by Sample ID") +
  theme(plot.title = element_text(size = 10, face = "bold"))


pdf(width=15, height=8,here(plot_dir, "UMAPs_Fig5_wSamples.pdf"))
(p1+p2+plot_spacer())/(p3+p4+p5)
dev.off()

png(width=3200, height=2400,here(plot_dir, "UMAPs_Fig5.png"), res=300)
(p1+p2)/(p3+p4)
dev.off()


# ======= Plot expression of marker genes to make sure the outleirs are not ENDO ==========
rownames(spe) <- rowData(spe)$gene_name
features <- c("SNAP25","GFAP",
              "SLC17A7","GAD2",
              "MBP", "MOBP",
              "PECAM1", "TIE1"
              )

pdf(width=8, height=12,here(plot_dir, "Expession_LocalOutliers_marker_genes.pdf"))
plotExpression(spe, features=features,colour_by="local_outliers", x="local_outliers")
dev.off()


# ===== Check for batch effects =====

# plot UMAPs by sample ID
plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="sample_id", point_size = point_size) + 
  ggtitle("Prior to batch correction") +
  theme(plot.title = element_text(size = 10, face = "bold"))



