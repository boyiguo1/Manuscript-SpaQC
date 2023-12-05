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


# ========= Spatial outlier detection =======
spe <- localOutliers(spe, 
                     features=c("sum_umi","sum_gene", "expr_chrM_ratio"), 
                     n_neighbors=36, 
                     data_output=TRUE,
                     method="multivariate",
                     cutoff=3)


# ========= Standard outlier detection (MAD) =======
library(scuttle)

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
colnames(colData(spe))


# # ===== Average percentage discarded across samples =====
# Extract the required columns
discard_mad_col <- spe$discard_mad
discard_local_col <- spe$local_outliers
layer_guess_ordered_col <- spe$layer_guess_reordered
sample_id_col <- spe$sample_id
subject_col <- spe$subject




# Create a data frame
df <- data.frame(discard_mad = discard_mad_col,
                 discard_local = as.logical(discard_local_col),
                 layer_guess_ordered = layer_guess_ordered_col,
                 sample_id = sample_id_col,
                 subject = subject_col)

df <- df[!is.na(df$layer_guess_ordered), ]

summary_df <- df %>%
  group_by(sample_id, layer_guess_ordered, subject) %>%
  summarise(
    total_spots = n(),  # Total number of spots in each domain for each sample
    discarded_mad = sum(discard_mad, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    discarded_local = sum(discard_local, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    .groups = 'keep'  # Dropping the grouping after summarise
  ) %>%
  mutate(
    percentage_mad = (discarded_mad / total_spots) * 100  # Calculating the percentage
  ) %>%
  mutate(
    percentage_local = (discarded_local/ total_spots) * 100  # Calculating the percentage
  )



head(summary_df)
# # A tibble: 6 × 8
# # Groups:   sample_id, layer_guess_ordered, subject [6]
#   sample_id layer_guess_ordered subject total_spots discarded_spots_mad discarded_spots_local
# <chr>     <fct>               <chr>         <int>               <int>                 <int>
# 1 151507    Layer1              Br5292          817                  14                     1
# 2 151507    Layer2              Br5292          305                   0                     0
# 3 151507    Layer3              Br5292         1215                  10                     4
# 4 151507    Layer4              Br5292          369                   9                     1
# 5 151507    Layer5              Br5292          675                   3                     1
# 6 151507    Layer6              Br5292          486                   0                     0
# # ℹ 2 more variables: percentage_discarded_mad <dbl>, percentage_discarded_local <dbl>



p1 <- ggplot(summary_df, aes(x = layer_guess_ordered, y = percentage_mad, fill = layer_guess_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "Percentage Discarded: Global",
       x = "Spatial Domain",
       y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 35) +
  guides(fill = FALSE, color = FALSE)

p2 <- ggplot(summary_df, aes(x = layer_guess_ordered, y = percentage_local, fill = layer_guess_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "Percentage Discarded: Local",
       x = "Spatial Domain",
       y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 35) +
  guides(fill = FALSE, color = FALSE)

pdf(here(plot_dir, "mad_vs_local_discarded_boxplots.pdf"), width = 10, height = 5)
gridExtra::grid.arrange(p1, p2, nrow = 1)
dev.off()

# I now want to compare the total number of spots retained for sample for MAD vs local outlier methods
# update the summary_df so that it has the total number of spots retained for both mad and local methods
summary_df <- summary_df %>%
  mutate(
    retained_mad = total_spots - discarded_mad,
    retained_local = total_spots - discarded_local
  )

# make it so it's total across layers. like sum the layer info
summary_df <- summary_df %>%
  group_by(sample_id, subject) %>%
  summarise(
    total_spots = sum(total_spots),
    discarded_mad = sum(discarded_mad),
    discarded_local = sum(discarded_local),
    retained_mad = sum(retained_mad),
    retained_local = sum(retained_local),
    .groups = 'keep'
  ) 

# pivot mad and local so they are in a single column for easy plotting
summary_df <- summary_df %>%
  pivot_longer(cols = c("discarded_mad", "discarded_local"), 
               names_to = "method", 
               values_to = "discarded")
head(summary_df)
# A tibble: 6 × 7
# Groups:   sample_id, subject [3]
# sample_id subject total_spots discarded_mad discarded_local method         retained
# <chr>     <chr>         <int>         <int>           <int> <chr>             <int>
#   1 151507    Br5292         4221            45               7 retained_mad       4176
# 2 151507    Br5292         4221            45               7 retained_local     4214
# 3 151508    Br5292         4381           111               4 retained_mad       4270
# 4 151508    Br5292         4381           111               4 retained_local     4377
# 5 151509    Br5292         4788            59              10 retained_mad       4729
# 6 151509    Br5292         4788            59              10 retained_local     4778

# plot the discarded spots across methods using bar plots for each sample
p1 <- ggplot(summary_df, aes(x = sample_id, y = discarded, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("red", "blue")) +
  labs(title = "Discarded Spots",
       x = "Sample ID",
       y = "Number of Spots") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(here(plot_dir,"mad_vs_local_discarded_barplots.pdf"), width = 10, height = 5)
p1
dev.off()


# ### add arbitrary threshold outliers and plot umaps
# ======== Let's discard based on these conservative threshold and see how it looks =======

qc_mito <- colData(spe)$expr_chrM_ratio > 0.3 # this is arbitrary
table(qc_mito)
# FALSE  TRUE 
# 47401   280 

qc_umi <- colData(spe)$sum_umi < 500 # this is arbitrary
table(qc_umi)
# FALSE  TRUE 
# 47372   309 

qc_gene <- colData(spe)$sum_gene < 500 # this is arbitrary
table(qc_gene)

colData(spe)$qc_mito_threshold <- qc_mito
colData(spe)$qc_umi_threshold <- qc_umi
colData(spe)$qc_gene_threshold <- qc_gene

sc_outliers <- qc_umi | qc_mito | qc_gene
spe$discard_threshold <- sc_outliers


# ======== UMAPS ==========

spe$clusters <- spe$layer_guess_reordered
point_size <- 0.3
p1 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="clusters", point_size = point_size) +
  ggtitle("Manually annotated cluster assignment") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p2 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="sample_id", point_size = point_size) + 
  ggtitle("UMAP by Sample") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p3 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="sum_umi_log2", point_size = point_size) +
  ggtitle("Total library size") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p4 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="discard_threshold", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red")) +
  ggtitle("Spots excluded by arbitrary threshold") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p5 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="discard_mad", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red")) +
  ggtitle("Spots excluded by 3 MAD") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p6 <- plotReducedDim(spe, dimred="UMAP_neighbors15",
                     colour_by="local_outliers", point_size = point_size) + 
  scale_color_manual(values = c("grey", "red")) +
  ggtitle("Spots excluded by local outliers (SpotSweeper LOF)") +
  theme(plot.title = element_text(size = 10, face = "bold"))

(p1+p2+p3)/(p4+p5+p6)

pdf(here(plot_dir,"UMAPs_Fig5_wSamples.pdf"), width = 15, height = 7.5)
(p1+p2+p3)/(p4+p5+p6)
dev.off()

# ==== instead of using plotReducedDim, plot the UMAP using ggplot2

# get the UMAP coordinates
umap_df <- reducedDim(spe)
  