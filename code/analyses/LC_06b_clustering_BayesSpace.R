##########################################################
# LC analyses: spatially aware clustering using BayesSpace
# Lukas Weber, Jun 2022
##########################################################

# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(BayesSpace)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "06b_clustering_BayesSpace")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrectedHarmony.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# -----------------------------------------------------------
# BayesSpace: offset spatial coordinates for multiple samples
# -----------------------------------------------------------

# modify spatial coordinates by adding offsets to row and column indices per 
# sample for clustering across multiple samples
# https://edward130603.github.io/BayesSpace/articles/joint_clustering.html

row <- colData(spe)$array_row
col <- colData(spe)$array_col

# 8 samples in total
length(sample_ids)

# add offsets to rows (+100) and/or columns (+150) to get nonoverlapping coordinates
# note: this needs to be adjusted manually for different numbers of samples

for (s in 2:4) {
  ix <- colData(spe)$sample_id == sample_ids[s]
  row[ix] <- row[ix] + (100 * (s - 1))
}
for (s in 6:8) {
  ix <- colData(spe)$sample_id == sample_ids[s]
  row[ix] <- row[ix] + (100 * (s - 5))
}
for (s in 5:8) {
  ix <- colData(spe)$sample_id == sample_ids[s]
  col[ix] <- col[ix] + 150
}


# store spatial coordinates in format expected by BayesSpace
colData(spe)$row <- row
colData(spe)$col <- col


pal <- unname(palette.colors(8, palette = "Okabe-Ito"))


# check offsets give non-overlapping samples
df <- cbind.data.frame(colData(spe), spatialCoords(spe))
ggplot(df, aes(x = row, y = col, color = sample_id)) + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = pal) + 
  coord_fixed() + 
  ggtitle("BayesSpace: multiple samples offsets check") + 
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme_bw()

fn <- file.path(dir_plots, "BayesSpace_multipleSamplesOffsetsCheck")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 4)
ggsave(paste0(fn, ".png"), width = 6.75, height = 4)


# --------------
# run BayesSpace
# --------------

# spatially-aware clustering using BayesSpace
# aiming to identify LC vs. WM regions in an unsupervised manner


# run once on all batch-corrected samples combined

# arguments:
# q = number of clusters
# d = number of input dimensions to use
# use.dimred: using Harmony batch integrated dimensions
# nrep: 10,000 iterations

# runtime: ~30 mins on laptop

set.seed(123)
spe <- spatialCluster(
  spe, 
  q = 4, 
  use.dimred = "HARM", 
  d = 15, 
  platform = "Visium", 
  init.method = "mclust", 
  nrep = 10000, 
  burn.in = 1000
)

# results
table(colData(spe)$spatial.cluster)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LC_BayesSpace")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

