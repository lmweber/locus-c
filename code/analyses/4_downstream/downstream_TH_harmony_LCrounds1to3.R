##################################################################
# LC project
# Downstream analyses on TH+ spots: with Harmony batch integration
# Lukas Weber, Nov 2021
##################################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(ggplot2)
library(scater)
library(scran)
library(harmony)


# ---------
# load data
# ---------

# load SpatialExperiment object

fn_spe <- here("processed_data", "SPE", "LCrounds1to3_SPE_processed.rds")
spe <- readRDS(fn_spe)
dim(spe)

# subset to keep only TH+ spots

# define threshold as >n detected UMIs per spot
thresh_TH <- 0
spe <- spe[, colData(spe)$TH_sum > thresh_TH]

dim(spe)  ## 1786 spots across all samples

# sample IDs
sample_ids <- unique(colData(spe)$sample_id)


# ------------------------------------
# normalization and log-transformation
# ------------------------------------

# quick clustering for pool-based size factors
set.seed(123)
qclus <- quickCluster(spe)
table(qclus)

# calculate size factors
spe <- computeSumFactors(spe, cluster = qclus)
summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 20)

# calculate logcounts
spe <- logNormCounts(spe)
assayNames(spe)


# -----------------
# feature selection
# -----------------

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$symbol)
table(is_mito)
rowData(spe)$symbol[is_mito]

# remove mitochondrial genes
spe <- spe[!is_mito, ]
dim(spe)

# fit mean-variance relationship
dec <- modelGeneVar(spe)
# visualize mean-variance relationship
fit <- metadata(dec)
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)


# ------------------------
# dimensionality reduction
# ------------------------

# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)
reducedDimNames(spe)
dim(reducedDim(spe, "PCA"))

# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
reducedDimNames(spe)
dim(reducedDim(spe, "UMAP"))
# update column names for easier plotting
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)


# ------------------------------
# batch integration with Harmony
# ------------------------------

set.seed(100)

# check inputs for Harmony
dim(reducedDim(spe, "PCA"))
length(colData(spe)$sample_id)
stopifnot(nrow(reducedDim(spe, "PCA")) == length(colData(spe)$sample_id))

# skip PCA since already provided in input
harmony_embeddings <- HarmonyMatrix(
  data_mat = reducedDim(spe, "PCA"), 
  meta_data = colData(spe)$sample_id, 
  do_pca = FALSE, 
  verbose = TRUE
)

# Harmony outputs: embeddings can be used as inputs downstream (instead of PCs)
str(harmony_embeddings)
dim(harmony_embeddings)
head(harmony_embeddings, 3)

stopifnot(all(rownames(harmony_embeddings) == rownames(reducedDim(spe, "PCA"))))

# store Harmony embeddings in SPE object
reducedDims(spe)[["Harmony"]] <- harmony_embeddings
reducedDimNames(spe)
head(reducedDim(spe, "Harmony"), 3)


# ------------------------------------------
# clustering (TH+ spots; Harmony embeddings)
# ------------------------------------------

# graph-based clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "Harmony")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
# store cluster labels in column 'label' in colData
colLabels(spe) <- factor(clus)


# -------------
# plot clusters
# -------------

df <- cbind.data.frame(colData(spe), spatialData(spe), spatialCoords(spe), reducedDim(spe, "UMAP"))

# plot clustering: cluster labels
ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) + 
  geom_point(size = 0.3) + 
  ggtitle("UMAP: TH+ spots, Harmony, cluster labels") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos_harmony", "THpos_harmony_clustering_labels"), ".pdf"), width = 4.5, height = 4)
ggsave(paste0(here("plots", "THpos_harmony", "THpos_harmony_clustering_labels"), ".png"), width = 4.5, height = 4)

# plot clustering: sample IDs
ggplot(df, aes(x = UMAP1, y = UMAP2, color = sample_id)) + 
  geom_point(size = 0.3) + 
  ggtitle("UMAP: TH+ spots, Harmony, sample IDs") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos_harmony", "THpos_harmony_clustering_sampleIDs"), ".pdf"), width = 5.5, height = 4)
ggsave(paste0(here("plots", "THpos_harmony", "THpos_harmony_clustering_sampleIDs"), ".png"), width = 5.5, height = 4)

# plot clustering: labels by sample
ggplot(df, aes(x = x, y = y, color = label)) + 
  facet_wrap(~ sample_id, nrow = 2) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  ggtitle("TH+ spots, Harmony, cluster labels by sample") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos_harmony", "THpos_harmony_clustering_labels_by_sample"), ".pdf"), width = 12, height = 6.5)
ggsave(paste0(here("plots", "THpos_harmony", "THpos_harmony_clustering_labels_by_sample"), ".png"), width = 12, height = 6.5)


# --------------------------
# marker gene identification
# --------------------------

# code from OSCA

marker.info <- scoreMarkers(spe, colLabels(spe))
marker.info

# plot for each cluster

rownames(spe) <- rowData(spe)$symbol

for (i in names(marker.info)) {
  chosen <- marker.info[[i]]
  ordered <- chosen[order(chosen$mean.AUC, decreasing = TRUE), ]
  head(ordered[, 1:4])
  plotExpression(spe, features = head(rownames(ordered)), 
                 x = "label", colour_by = "label")
  ggsave(paste0(here("plots", "THpos_harmony", "THpos_harmony_clustering_markers_cluster"), i, ".pdf"), 
         width = 6, height = 6.5, bg = "white")
  ggsave(paste0(here("plots", "THpos_harmony", "THpos_harmony_clustering_markers_cluster"), i, ".png"), 
         width = 6, height = 6.5, bg = "white")
}

