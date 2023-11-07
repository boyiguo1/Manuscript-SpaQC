#devtools::install("/Users/mtotty2/Documents/R/SpotSweeper")
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

plot_dir = here('plots', 'spatial_methods', '02_outlier_detection', 'SpotSweeper')
processed_dir = here('processed-data')

## Download the spot-level data
spe <- fetch_data(type = "spe")
spe




# ====================== Testing with my newly made function ======================

# I turned the above code into functions in the private SpotSweeper package
# https://github.com/MicTott/SpotSweeper

spe <- localOutliers(spe, k= 6, threshold=3)


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

pdf(width = 12.5, height = 12.5, here(plot_dir, 'Many_SpotPlots_sumUmi_k6.pdf'))
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

spe <- localOutliers(spe, k=36 , threshold=3)


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

pdf(width = 12.5, height = 12.5, here(plot_dir, 'Many_SpotPlots_sumUmi_k36.pdf'))
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

pdf(width=8, height=4,here(plot_dir, "Boxplot_percent_discarded_k36.pdf"))
p
dev.off()



