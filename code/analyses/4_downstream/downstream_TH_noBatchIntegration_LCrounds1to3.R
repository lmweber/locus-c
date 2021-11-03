##################################
# LC project
# Downstream analyses on TH+ spots
# Lukas Weber, Nov 2021
##################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(ggplot2)
library(scater)
library(scran)


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


# ----------------------------------
# additional QC metrics on TH+ spots
# ----------------------------------

df <- cbind.data.frame(colData(spe), spatialData(spe), spatialCoords(spe))

# cell count (excluding round 2 samples)
ggplot(df[!grepl("round2", df$sample_id), ], aes(x = count)) + 
  geom_bar(fill = "navy") + 
  labs(x = "cell count (excluding round 2 samples)") + 
  ggtitle("QC metrics: TH+ spots") + 
  theme_bw()
ggsave(paste0(here("plots", "THpos", "THpos_cell_count"), ".pdf"), width = 5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_cell_count"), ".png"), width = 5, height = 4)

# sum
ggplot(df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "navy") + 
  labs(x = "sum UMIs per spot") + 
  ggtitle("QC metrics: TH+ spots") + 
  theme_bw()
ggsave(paste0(here("plots", "THpos", "THpos_sum"), ".pdf"), width = 5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_sum"), ".png"), width = 5, height = 4)

# detected
ggplot(df, aes(x = detected)) + 
  geom_histogram(color = "black", fill = "navy") + 
  labs(x = "detected genes per spot") + 
  ggtitle("QC metrics: TH+ spots") + 
  theme_bw()
ggsave(paste0(here("plots", "THpos", "THpos_detected"), ".pdf"), width = 5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_detected"), ".png"), width = 5, height = 4)


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


# ----------------------
# clustering (TH+ spots)
# ----------------------

# graph-based clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
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
  ggtitle("UMAP: TH+ spots, cluster labels") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos", "THpos_clustering_labels"), ".pdf"), width = 4.5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_clustering_labels"), ".png"), width = 4.5, height = 4)

# plot clustering: sample IDs
ggplot(df, aes(x = UMAP1, y = UMAP2, color = sample_id)) + 
  geom_point(size = 0.3) + 
  ggtitle("UMAP: TH+ spots, sample IDs") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos", "THpos_clustering_sampleIDs"), ".pdf"), width = 4.5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_clustering_sampleIDs"), ".png"), width = 4.5, height = 4)

# plot clustering: expression of TH
ggplot(df, aes(x = UMAP1, y = UMAP2, color = TH_sum)) + 
  geom_point(size = 0.3) + 
  ggtitle("UMAP: TH+ spots, expression of TH") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos", "THpos_clustering_THexpr"), ".pdf"), width = 4.5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_clustering_THexpr"), ".png"), width = 4.5, height = 4)

# plot clustering: sum UMIs
ggplot(df, aes(x = UMAP1, y = UMAP2, color = sum)) + 
  geom_point(size = 0.3) + 
  ggtitle("UMAP: TH+ spots, sum UMIs") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos", "THpos_clustering_sum"), ".pdf"), width = 4.5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_clustering_sum"), ".png"), width = 4.5, height = 4)

# plot clustering: detected
ggplot(df, aes(x = UMAP1, y = UMAP2, color = detected)) + 
  geom_point(size = 0.3) + 
  ggtitle("UMAP: TH+ spots, detected") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos", "THpos_clustering_detected"), ".pdf"), width = 4.5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_clustering_detected"), ".png"), width = 4.5, height = 4)

# plot clustering: mito percent
ggplot(df, aes(x = UMAP1, y = UMAP2, color = subsets_mito_percent)) + 
  geom_point(size = 0.3) + 
  ggtitle("UMAP: TH+ spots, mito percent") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "THpos", "THpos_clustering_mito"), ".pdf"), width = 4.5, height = 4)
ggsave(paste0(here("plots", "THpos", "THpos_clustering_mito"), ".png"), width = 4.5, height = 4)

