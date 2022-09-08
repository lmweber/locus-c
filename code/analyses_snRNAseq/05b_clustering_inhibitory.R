###################################################################
# LC snRNA-seq analyses: secondary clustering of inhibitory neurons
# Lukas Weber, Sep 2022
###################################################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(bluster)
library(ggplot2)
library(ggVennDiagram)


dir_plots <- here("plots", "snRNAseq", "05b_clustering_inhibitory")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_clustering")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)


# -------------------------
# Subset inhibitory neurons
# -------------------------

sce_full <- sce

# select inhibitory neuron clusters (from marker expression in heatmap)
clus_select <- c(26, 17, 14, 1, 8, 7, 24, 18)

ix_select <- colLabels(sce) %in% clus_select
table(ix_select)

sce <- sce[, ix_select]

dim(sce)


# --------------------
# Secondary clustering
# --------------------

# secondary clustering of inhibitory neurons

# clustering algorithm and parameters from OSCA
# two-stage clustering algorithm using high-resolution k-means and graph-based clustering

set.seed(100)
clus <- clusterCells(
  sce, 
  use.dimred = "PCA", 
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 2000), 
    second = NNGraphParam(k = 10)
  )
)
colLabels(sce) <- clus


# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# expression of key markers for NE, 5-HT, and cholinergic neuron populations
ix <- c(
  TH = which(rowData(sce)$gene_name == "TH"), 
  SLC6A2 = which(rowData(sce)$gene_name == "SLC6A2"), 
  DBH = which(rowData(sce)$gene_name == "DBH"), 
  SLC6A4 = which(rowData(sce)$gene_name == "SLC6A4"), 
  TPH2 = which(rowData(sce)$gene_name == "TPH2"), 
  SLC5A7 = which(rowData(sce)$gene_name == "SLC5A7"), 
  CHAT = which(rowData(sce)$gene_name == "CHAT"), 
  ACHE = which(rowData(sce)$gene_name == "ACHE"), 
  BCHE = which(rowData(sce)$gene_name == "BCHE"), 
  SLC18A3 = which(rowData(sce)$gene_name == "SLC18A3"), 
  PRIMA1 = which(rowData(sce)$gene_name == "PRIMA1")
)

n_clus <- length(table(colLabels(sce)))

res_list <- list()
for (k in seq_len(n_clus)) {
  res_list[[k]] <- rowMeans(logcounts(sce)[ix, colLabels(sce) == k])
}
res_mat <- do.call("rbind", res_list)
rownames(res_mat) <- seq_len(n_clus)
colnames(res_mat) <- names(ix)

cbind(
  n = table(colLabels(sce)), 
  table(colLabels(sce), colData(sce)$Sample), 
  res_mat
)


# --------------------
# Store cluster labels
# --------------------

# store secondary clustering labels in full SCE object

# check unique barcode IDs
table(duplicated(colnames(sce)))
table(duplicated(colnames(sce_full)))
table(colnames(sce) %in% colnames(sce_full))

# match and store cluster labels
clus_secondary <- rep(NA, ncol(sce_full))
names(clus_secondary) <- colnames(sce_full)
clus_secondary[colnames(sce)] <- colData(sce)$label

colData(sce_full)$label_secondary <- clus_secondary

# check
table(colData(sce_full)$label)
table(colData(sce_full)$label_secondary)


# -----------
# Save object
# -----------

# note saving 'sce_full' object

fn_out <- here("processed_data", "SCE", "sce_clustering_secondary")
saveRDS(sce_full, paste0(fn_out, ".rds"))

