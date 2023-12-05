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
library(ggridges)

plot_dir = here('plots', 'manuscript_figures',"figure_1")
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

# ======== Arbitrary threshold QC =======

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

threshold_outliers <- qc_umi | qc_mito | qc_gene
spe$discard_threshold <- threshold_outliers


# # ===== Average percentage discarded across samples =====
# Extract the required columns
discard_mad_col <- spe$discard_mad
discard_local_col <- spe$local_outliers
discard_threshold_col <- spe$discard_threshold
layer_guess_ordered_col <- spe$layer_guess_reordered
sample_id_col <- spe$sample_id
subject_col <- spe$subject

# Create a data frame
df <- data.frame(discard_mad = discard_mad_col,
                 discard_threshold = discard_threshold_col,
                 discard_local = as.logical(discard_local_col),
                 layer_guess_ordered = layer_guess_ordered_col,
                 sample_id = sample_id_col,
                 subject = subject_col)

df <- df[!is.na(df$layer_guess_ordered), ]

discarded_df <- df %>%
  group_by(sample_id, layer_guess_ordered, subject) %>%
  summarise(
    total_spots = n(),  # Total number of spots in each domain for each sample
    discarded_mad = sum(discard_mad, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    discarded_local = sum(discard_local, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    discarded_threshold = sum(discard_threshold, na.rm = TRUE),  # Counting TRUE as 1 and FALSE as 0
    .groups = 'keep'  # Dropping the grouping after summarise
  ) %>%
  mutate(
    percentage_mad = (discarded_mad / total_spots) * 100  # Calculating the percentage
  ) %>%
  mutate(
    percentage_local = (discarded_local/ total_spots) * 100  # Calculating the percentage
  )%>%
  mutate(
    percentage_threshold = (discarded_threshold/ total_spots) * 100  # Calculating the percentage
  )



head(discarded_df)
# # A tibble: 6 × 10
# # Groups:   sample_id, layer_guess_ordered, subject [6]
# sample_id layer_guess_ordered subject total_spots discarded_mad discarded_local
# <chr>     <fct>               <chr>         <int>         <int>           <int>
#   1 151507    Layer1              Br5292          817            14               1
# 2 151507    Layer2              Br5292          305             0               0
# 3 151507    Layer3              Br5292         1215            10               4
# 4 151507    Layer4              Br5292          369             9               1
# 5 151507    Layer5              Br5292          675             3               1
# 6 151507    Layer6              Br5292          486             0               0
# # ℹ 4 more variables: discarded_threshold <int>, percentage_mad <dbl>, percentage_local <dbl>,
# #   percentage_threshold <dbl>


# ===== Make df of QC metrics across layers =====

# Extract data 
expr_data <- colData(spe)$expr_chrM
expr_ratio <- colData(spe)$expr_chrM_ratio
layer_data <- colData(spe)$layer_guess_reordered
cell_count <- colData(spe)$cell_count
sum_umi <- colData(spe)$sum_umi
sum_gene <- colData(spe)$sum_gene

# Combine the expression data with the layer data
qc_df <- data.frame(
  expr_chrM = expr_data,
  expr_ratio = expr_ratio,
  layer = layer_data,
  sum_umi = sum_umi,
  cell_count = cell_count,
  sum_gene = sum_gene
)


# drop "NA" layer
qc_df <- qc_df[!is.na(qc_df$layer),]

# reverse the y order so that layer 1 is at the top
qc_df$layer <- factor(qc_df$layer, levels = rev(levels(qc_df$layer)))




# ================== From here we can start to generate plots ==================

# ===== Panel A =====
# visualizng the manual layer annotations and common QC metrics

point_size <- 1.9
p1 <- vis_clus(
  spe = spe,
  clustervar = "spatialLIBD",
  sampleid = unique(spe$sample_id)[1],
  colors = libd_layer_colors,
  spatial=FALSE,
  point_size = point_size) + 
  theme_void() + 
  ggtitle("Manual layer annotations") + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=14)) 

p2 <- vis_gene(
  spe = spe,
  geneid = "expr_chrM_ratio",
  sampleid = unique(spe$sample_id)[1],
  spatial=FALSE,
  point_size = point_size) + 
  theme_void() + 
  ggtitle("Mito percentage") + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=14)) 


p3 <- vis_gene(
  spe = spe,
  geneid = "sum_gene",
  sampleid = unique(spe$sample_id)[1],
  spatial=FALSE,
  point_size = point_size) + 
  theme_void() + 
  ggtitle("Unique genes") + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=14)) 

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_A_1.pdf'))
p1
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_A_2.pdf'))
p2
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_A_3.pdf'))
p3
dev.off()


# ===== Panel B =====

# ridge plot of sum_umi with a threshold of 500
p4 <- ggplot(qc_df, aes(x = sum_umi, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", size=1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Library size per layer",
       x = "Library Size",
       y = "Spatial Domain")

# ridge plot of expr_ratio with a threshold of 0.3
p5 <- ggplot(qc_df, aes(x = expr_ratio, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "red", size=1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Mito percent per layer",
       x = "Mito %",
       y = "Spatial Domain")

# ridge plot of sum_gene with a threshold of 500
p6 <- ggplot(qc_df, aes(x = sum_gene, y = layer, fill = layer)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", size=1) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20)) +
  scale_fill_manual(values = libd_layer_colors) +
  labs(title = "Unique genes per layer",
       x = "Number of genes",
       y = "Spatial Domain")

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_B_1.pdf'))
p4
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_B_2.pdf'))
p5
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_B_3.pdf'))
p6
dev.off()


# ===== Panel C =====


point_size <- 1.75
p7 <- vis_clus(
  spe = spe,
  clustervar = "spatialLIBD",
  sampleid = unique(spe$sample_id)[1],
  colors = libd_layer_colors,
  spatial=FALSE,
  point_size = point_size) + 
  theme_void() + 
  ggtitle("Manual layer annotations") + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=14)) 


p8 <- vis_clus(
  spe = spe,
  clustervar = "discard_threshold",
  sampleid = unique(spe$sample_id)[1],
  colors = c("grey", "red"),
  spatial=FALSE,
  point_size = point_size) + 
  theme_void() + 
  ggtitle("Arbitrary threshold") + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=14)) 


p9 <- vis_clus(
  spe = spe,
  clustervar = "discard_mad",
  sampleid = unique(spe$sample_id)[1],
  colors = c("grey", "red"),
  spatial=FALSE,
  point_size = point_size) + 
  theme_void() + 
  ggtitle("3 MAD") + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size=14)) 

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_C_1.pdf'))
p7
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_C_2.pdf'))
p8
dev.off()

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_C_3.pdf'))
p9
dev.off()

# ===== Panel D ======

p <- ggplot(discarded_df, aes(x = layer_guess_ordered, y = percentage_threshold, fill = layer_guess_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "Arbitrary threshold",
       x = "Spatial Domain",
       y = "Percent discarded (%)") +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        legend.position = "none") 

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_D.pdf'))
p
dev.off()


# ===== Panel E ======

p <- ggplot(discarded_df, aes(x = layer_guess_ordered, y = percentage_mad, fill = layer_guess_ordered)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = layer_guess_ordered), alpha = 0.7) +  # Display all data points with matching layer colors
  scale_fill_manual(values = spatialLIBD::libd_layer_colors) +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "3 MAD",
       x = "Spatial Domain",
       y = "Percent discarded (%)") +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        legend.position = "none")

pdf(height = 5, width=5, here(plot_dir, 'Figure1_Panel_E.pdf'))
p
dev.off()



# sing the discarded df, get just the WM and L1 spots. Then do a scatterplot of total spots vs discarded spots for each layer
# ===== Panel F ======

# Get just the WM and L1 spots
discarded_df_WM_L1 <- discarded_df %>% 
  filter(layer_guess_ordered %in% c("Layer1"))

# Scatterplot of total spots vs discarded spots for each layer
ggplot(discarded_df_WM_L1, aes(x = total_spots, y = discarded_mad, color = layer_guess_ordered)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = spatialLIBD::libd_layer_colors) +
  labs(title = "3 MAD",
       x = "Total spots",
       y = "Discarded spots") +
  theme_bw(base_size=20) +
  theme(legend.position = "none")
