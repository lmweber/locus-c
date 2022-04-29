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
dir_plots <- here("plots", "06_downstream", "combined")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrected.rds")
spe <- readRDS(fn_spe)

dim(spe)


# using all spots


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

fn <- file.path(dir_plots, "LC_clusters_HARM_combined")
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

fn <- file.path(dir_plots, "LC_clusters_PCA_combined")
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

fn <- file.path(dir_plots, "LC_clusters_UMAP_combined")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# x-y space (by sample)
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = label)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
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

fn <- file.path(dir_plots, "LC_clusters_XYspace_combined")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 7)
ggsave(paste0(fn, ".png"), width = 6.75, height = 7)


# --------------------------
# marker gene identification
# --------------------------

# code from OSCA

rownames(spe) <- rowData(spe)$gene_name

marker.info <- scoreMarkers(spe, colLabels(spe))
marker.info


# plot top markers for each cluster

for (i in names(marker.info)) {
  chosen <- marker.info[[i]]
  ordered <- chosen[order(chosen$mean.AUC, decreasing = TRUE), ]
  head(ordered[, 1:4])
  plotExpression(spe, features = head(rownames(ordered)), 
                 x = "label", colour_by = "label") + 
    scale_color_manual(values = pal, name = "label") + 
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
  fn <- file.path(dir_plots, "markers_combined", paste0("markers_cluster", i))
  ggsave(paste0(fn, ".pdf"), width = 10, height = 6, bg = "white")
  ggsave(paste0(fn, ".png"), width = 10, height = 6, bg = "white")
}


# plot selected genes across clusters

genes <- c("SNAP25", "SYT1", "TH", "SLC6A2", "TPH2", "SLC6A4", "DBH", "DDC")

plotExpression(spe, features = genes, x = "label", colour_by = "label", ncol = 2) + 
  scale_color_manual(values = pal, name = "label") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

fn <- file.path(dir_plots, "selectedMarkers_byCluster_combined")
ggsave(paste0(fn, ".pdf"), width = 9, height = 7, bg = "white")
ggsave(paste0(fn, ".png"), width = 9, height = 7, bg = "white")


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

rownames(tab_spots) <- c("not_spots", "spots")
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
  mutate(cluster = factor(cluster, levels = unique(cluster))) %>% 
  as.data.frame()

pal <- c("white", "navy")

ggplot(df, aes(x = annotation, y = cluster, fill = proportion)) + 
  facet_wrap(~type, scales = "free") + 
  geom_tile() + 
  scale_fill_gradientn(colors = pal) + 
  ggtitle("Clustering in combined region")

fn <- file.path(dir_plots, "clusterAnnotationComparison_combinedRegions")
ggsave(paste0(fn, ".pdf"), width = 6, height = 5)
ggsave(paste0(fn, ".png"), width = 6, height = 5)

