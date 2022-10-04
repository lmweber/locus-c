##############################################################################
# LC snRNA-seq analyses: normalization, feature selection, dimension reduction
# Lukas Weber, Sep 2022
##############################################################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_qualityControlled")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)


# ---------------------------
# Normalization and logcounts
# ---------------------------

# normalization by deconvolution

set.seed(123)
quickclus <- quickCluster(sce)
table(quickclus)

sce <- computeSumFactors(sce, cluster = quickclus)
sce <- logNormCounts(sce)

summary(sizeFactors(sce))


# -----------------
# Feature selection
# -----------------

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$gene_name)
table(is_mito)
rowData(sce)$gene_name[is_mito]

# fit mean-variance relationship
dec <- modelGeneVar(sce, subset.row = !is_mito)
# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)


# ------------------------
# Dimensionality reduction
# ------------------------

# note: not applying any batch integration due to our interest in rare populations (LC-NE neurons)

# calculate PCA (on top HVGs)
set.seed(123)
sce <- runPCA(sce, subset_row = top_hvgs)

# calculate UMAP (on top PCs)
set.seed(123)
sce <- runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) <- paste0("UMAP", seq_len(ncol(reducedDim(sce, "UMAP"))))

reducedDimNames(sce)


# -----------
# Save object
# -----------

fn_out <- here("processed_data", "SCE", "sce_logcounts")
saveRDS(sce, paste0(fn_out, ".rds"))

