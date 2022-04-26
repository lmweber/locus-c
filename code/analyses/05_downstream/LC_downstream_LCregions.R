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
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)


# directory to save plots
dir_plots <- here("plots", "05_downstream", "LC_regions")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrected.rds")
spe <- readRDS(fn_spe)

dim(spe)


# ------------------------------------
# select manually annotated LC regions
# ------------------------------------

# select spots from manually annotated LC regions

spe <- spe[, colData(spe)$annot_region]

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

fn <- file.path(dir_plots, "LC_clusters_HARM_LCregions")
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

fn <- file.path(dir_plots, "LC_clusters_PCA_LCregions")
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

fn <- file.path(dir_plots, "LC_clusters_UMAP_LCregions")
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

fn <- file.path(dir_plots, "LC_clusters_XYspace_LCregions")
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

fn <- file.path(dir_plots, "selectedMarkers_byCluster_LCregions")
ggsave(paste0(fn, ".pdf"), width = 6, height = 7, bg = "white")
ggsave(paste0(fn, ".png"), width = 6, height = 7, bg = "white")


# ------------------------------
# plot clustering vs. annotation
# ------------------------------

# compare cluster membership with manually annotated individual spots
tab_spots <- table(colData(spe)$annot_spot, colData(spe)$label)
tab_spots

# compare cluster membership with thresholding on TH expression
tab_TH <- table(colData(spe)$TH > 1, colData(spe)$label)
tab_TH

# convert to proportions
tab_spots <- apply(tab_spots, 2, function(col) col / sum(col))
tab_TH <- apply(tab_TH, 2, function(col) col / sum(col))

rownames(tab_spots) <- c("nonspots", "spots")
rownames(tab_TH) <- c("THneg", "THpos")

colnames(tab_spots) <- paste0("cluster", colnames(tab_spots))
colnames(tab_TH) <- paste0("cluster", colnames(tab_TH))

df <- 
  rbind(tab_spots, tab_TH) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "cluster") %>% 
  pivot_longer(., cols = -cluster, 
               names_to = "annotation", values_to = "proportion") %>% 
  mutate(type = ifelse(annotation %in% c("THneg", "THpos"), 
                       "TH_expression_threshold", "manually_annotated_spots")) %>% 
  as.data.frame()

pal <- c("white", "navy")

ggplot(df, aes(x = annotation, y = cluster, fill = proportion)) + 
  facet_wrap(~type, scales = "free") + 
  geom_tile() + 
  scale_fill_gradientn(colors = pal)

fn <- file.path(dir_plots, "clusterAnnotationComparison_LCregions")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

