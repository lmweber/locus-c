#############################
# LC project
# Script for batch correction
# Lukas Weber, Mar 2022
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
dir_plots <- here("plots", "02_batch_correction", "WMregions")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_preprocessed_QC.rds")
spe <- readRDS(fn_spe)

dim(spe)


# -----------------------------------
# clustering without batch correction
# -----------------------------------

# clustering on WM region spots without batch correction

# select WM spots
spe <- spe[, !colData(spe)$annot_region]
dim(spe)

# run standard clustering pipeline from OSTA

# normalization and logcounts
# note: use library size normalization
set.seed(123)
spe <- computeSumFactors(spe)
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


# plot sample IDs and clusters

df <- cbind.data.frame(colData(spe), spatialCoords(spe), 
                       reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"))

# sample IDs
ggplot(df, aes(x = UMAP1, y = UMAP2, color = sample_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  ggtitle("Sample IDs: no batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "WMregions_noBatchCorrection_sampleIDs")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5)
ggsave(paste0(fn, ".png"), width = 7, height = 5)


# clustering
ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  ggtitle("Clustering: no batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "WMregions_noBatchCorrection_clustering")
ggsave(paste0(fn, ".pdf"), width = 8, height = 5)
ggsave(paste0(fn, ".png"), width = 8, height = 5)


# ----------------------------------
# run batch correction using Harmony
# ----------------------------------

# run Harmony on PCA dimensions to integrate sample IDs

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


# plot sample IDs and clusters

df <- cbind.data.frame(colData(spe), spatialCoords(spe), 
                       reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"), 
                       reducedDim(spe, "HARM"))

# sample IDs in Harmony embedding space
ggplot(df, aes(x = HARM1, y = HARM2, color = sample_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  ggtitle("Harmony embeddings: sample IDs") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "WMregions_harmonyEmbeddings_sampleIDs")
ggsave(paste0(fn, ".pdf"), width = 6.5, height = 5)
ggsave(paste0(fn, ".png"), width = 6.5, height = 5)


# batch corrected clustering in UMAP space
ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  ggtitle("Batch corrected clustering") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "WMregions_harmonyBatchCorrected_clustering")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5)
ggsave(paste0(fn, ".png"), width = 7, height = 5)

