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


plot_dir = here('plots', 'spatial_methods', '02_outlier_detection', 'SpotSweeper')
processed_dir = here('processed-data', 'spatial_methods')

load(here('raw-data', "dlPFC_raw.Rdata"))
spe <- spe_raw
spe

rm(spe_raw)


# ====================== Testing with my newly made function ======================

# I turned the above code into functions in the private SpotSweeper package
# https://github.com/MicTott/SpotSweeper

spe <- localOutliers(spe, k= 15, threshold=2)


# ====================== Explore these outliers ======================

spe.subset <- subset(spe, ,sample_id == "Br2743_mid")
# Visualize
p1 <- make_escheR(spe.subset) |> 
  add_fill(var = "sum_umi_log2") +
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

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
  scale_fill_gradient2(low ="purple" , mid = "white",high =  "darkgreen")

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

pdf(width = 12.5, height = 12.5, here(plot_dir, 'Many_SpotPlots_sumUmi_k15.pdf'))
(p1+p2)/(p3+p4)
dev.off()