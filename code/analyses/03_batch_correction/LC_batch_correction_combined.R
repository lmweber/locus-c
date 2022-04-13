#############################
# LC project
# Script for batch correction
# Lukas Weber, Apr 2022
#############################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)
library(harmony)
library(ggplot2)
library(ggnewscale)


# directory to save plots
dir_plots <- here("plots", "02_batch_correction", "combined")


# ---------
# load data
# ---------

# load saved SPE object from previous script

# using combined data from both LC and WM regions

fn_spe <- here("processed_data", "SPE", "LC_qualityControlled.rds")
spe <- readRDS(fn_spe)

dim(spe)


# add donor id
colData(spe)$donor_id <- 
  gsub("_.*$", "", colData(spe)$sample_id) %>% 
  factor(., levels = c("Br6522", "Br8153", "Br5459", "Br2701", "Br8079"))

table(colData(spe)$donor_id)

colData(spe)$round_id <- 
  factor(colData(spe)$round_id, levels = c("round1", "round2", "round3"))

table(colData(spe)$round_id)


# to do: add donor_id and round_id (as factor) in previous scripts

# to do: add sample_part_id without NAs


# ------------------------
# without batch correction
# ------------------------

# run dimension reduction and clustering without batch correction
# demonstrates that batch correction is required in this dataset


# run standard clustering pipeline from OSTA

# normalization and logcounts
set.seed(123)
qclus <- quickCluster(spe)
table(qclus)
spe <- computeSumFactors(spe, cluster = qclus)
spe <- logNormCounts(spe)

# feature selection
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]
spe <- spe[!is_mito, ]
dec <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
head(top_hvgs)

# dimensionality reduction
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe) <- factor(clus)


# ---------
# PCA plots
# ---------

df <- cbind.data.frame(colData(spe), spatialCoords(spe), 
                       reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"))

pal <- c(unname(palette.colors(n = 8, palette = "Okabe-Ito")), 
         unname(palette.colors(n = 36, palette = "Polychrome 36")))

# donor ID
ggplot(df, aes(x = PC1, y = PC2, color = donor_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_PCA_donorID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# round ID
ggplot(df, aes(x = PC1, y = PC2, color = round_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_PCA_roundID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# sample ID
ggplot(df, aes(x = PC1, y = PC2, color = sample_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_PCA_sampleID")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 5)
ggsave(paste0(fn, ".png"), width = 6.75, height = 5)


# sample and part ID
ggplot(df, aes(x = PC1, y = PC2, color = sample_part_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_PCA_samplePartID")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5)
ggsave(paste0(fn, ".png"), width = 7, height = 5)


# clustering
ggplot(df, aes(x = PC1, y = PC2, color = label)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_PCA_clustering")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# ----------
# UMAP plots
# ----------

# donor ID
ggplot(df, aes(x = UMAP1, y = UMAP2, color = donor_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_UMAP_donorID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# round ID
ggplot(df, aes(x = UMAP1, y = UMAP2, color = round_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_UMAP_roundID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# sample ID
ggplot(df, aes(x = UMAP1, y = UMAP2, color = sample_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_UMAP_sampleID")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 5)
ggsave(paste0(fn, ".png"), width = 6.75, height = 5)


# sample and part ID
ggplot(df, aes(x = UMAP1, y = UMAP2, color = sample_part_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_UMAP_samplePartID")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5)
ggsave(paste0(fn, ".png"), width = 7, height = 5)


# clustering
ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_noBatchCorrection_UMAP_clustering")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# ----------------------------------
# run batch correction using Harmony
# ----------------------------------

# run Harmony on PCA dimensions to integrate sample IDs
# note: integrate on *sample IDs* based on plots above

pca_matrix <- reducedDim(spe, "PCA")
sample_ids <- colData(spe)$sample_id
stopifnot(nrow(pca_matrix) == length(sample_ids))

set.seed(123)
harmony_embeddings <- HarmonyMatrix(
  pca_matrix, 
  meta_data = sample_ids, 
  do_pca = FALSE
)

colnames(harmony_embeddings) <- paste0("HARM", seq_len(ncol(harmony_embeddings)))

dim(harmony_embeddings)
head(harmony_embeddings, 2)

# store in SPE object
reducedDims(spe) <- list(
  PCA = reducedDim(spe, "PCA"), 
  UMAP = reducedDim(spe, "UMAP"), 
  HARM = harmony_embeddings
)

reducedDims(spe)


# ----------------------------------------
# clustering with Harmony batch correction
# ----------------------------------------

# clustering using Harmony batch corrected PCs

# run standard clustering pipeline from OSTA

# clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "HARM")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe) <- factor(clus)


# ----------------------
# plots in Harmony space
# ----------------------

df <- cbind.data.frame(colData(spe), spatialCoords(spe), 
                       reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"), 
                       reducedDim(spe, "HARM"))

# donor ID
ggplot(df, aes(x = HARM1, y = HARM2, color = donor_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_Harmony_donorID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# round ID
ggplot(df, aes(x = HARM1, y = HARM2, color = round_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_Harmony_roundID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# sample ID
ggplot(df, aes(x = HARM1, y = HARM2, color = sample_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_Harmony_sampleID")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 5)
ggsave(paste0(fn, ".png"), width = 6.75, height = 5)


# sample and part ID
ggplot(df, aes(x = HARM1, y = HARM2, color = sample_part_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_Harmony_samplePartID")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5)
ggsave(paste0(fn, ".png"), width = 7, height = 5)


# clustering
ggplot(df, aes(x = HARM1, y = HARM2, color = label)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "LC_Harmony_clustering")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LC_batchCorrected")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

