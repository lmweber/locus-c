################################
# LC project
# Script for downstream analyses
# Lukas Weber, Apr 2022
################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "05_downstream", "individual_spots")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrected.rds")
spe <- readRDS(fn_spe)

dim(spe)


# ------------------------------------------
# select manually annotated individual spots
# ------------------------------------------

# select manually annotated individual spots only

spe <- spe[, colData(spe)$annot_spot]

dim(spe)


# ----------
# clustering
# ----------

# re-run clustering on selected spots

# clustering pipeline from OSTA using Harmony batch corrected PCs

set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "HARM")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe) <- factor(clus)


# -------------
# plot clusters
# -------------

df <- cbind.data.frame(
  colData(spe), spatialCoords(spe), 
  reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"), reducedDim(spe, "HARM")
)

pal <- c(unname(palette.colors(n = 8, palette = "Okabe-Ito")), 
         unname(palette.colors(n = 36, palette = "Polychrome 36")))


# Harmony embeddings
ggplot(df, aes(x = HARM1, y = HARM2, color = label)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Clustering: Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_clusters_HARM_individual")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# PCA
ggplot(df, aes(x = PC1, y = PC2, color = label)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Clustering: PCA") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_clusters_PCA_individual")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# UMAP
ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Clustering: UMAP") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_clusters_UMAP_individual")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# x-y space (by sample)
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = label)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.5) + 
  scale_color_manual(values = pal) + 
  # coord_fixed() +  ## use 'aspect.ratio = 1' instead with 'scales = "free"'
  scale_y_reverse() + 
  ggtitle("Clustering: x-y space") + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "LC_clusters_XYspace_individual")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 7)
ggsave(paste0(fn, ".png"), width = 6.75, height = 7)


# --------------------------
# marker gene identification
# --------------------------

# plot selected genes across clusters

rownames(spe) <- rowData(spe)$gene_name

genes <- c("SNAP25", "SYT1", "TH", "SLC6A2", "TPH2", "SLC6A4", "DBH", "DDC")

plotExpression(spe, features = genes, x = "label", colour_by = "label", ncol = 2) + 
  scale_color_manual(values = pal, name = "label") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

fn <- file.path(dir_plots, "selectedMarkers_byCluster_individual")
ggsave(paste0(fn, ".pdf"), width = 4, height = 7, bg = "white")
ggsave(paste0(fn, ".png"), width = 4, height = 7, bg = "white")

