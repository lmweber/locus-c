#####################################################
# LC snRNA-seq analyses: SCE object for ExperimentHub
# Lukas Weber, Oct 2022
#####################################################


library(here)
library(SingleCellExperiment)


# ----------------
# Load SCE objects
# ----------------

# load SCE objects from previous script


# SCE object without gene filtering
fn <- here("processed_data", "SCE", "sce_logcounts")
sce_logcounts <- readRDS(paste0(fn, ".rds"))

dim(sce_logcounts)


# SCE object containing cluster labels
fn <- here("processed_data", "SCE", "sce_clustering_merged")
sce_clustering_merged <- readRDS(paste0(fn, ".rds"))

dim(sce_clustering_merged)


# ------------
# Remove zeros
# ------------

# remove genes with zero expression

ix_zeros <- rowSums(counts(sce_logcounts)) == 0
table(ix_zeros)

dim(sce_logcounts)

sce_logcounts <- sce_logcounts[!ix_zeros, ]

dim(sce_logcounts)


# check no new zeros introduced
table(rowSums(counts(sce_logcounts)) == 0, useNA = "always")
table(colSums(counts(sce_logcounts)) == 0, useNA = "always")


# -----------------------------------
# Create SCE object for ExperimentHub
# -----------------------------------

# select columns of colData and store cluster labels in SCE object without gene filtering

stopifnot(ncol(sce_logcounts) == ncol(sce_clustering_merged))
stopifnot(all(colnames(sce_logcounts) == colnames(sce_clustering_merged)))
stopifnot(all(colData(sce_logcounts)$Sample == colData(sce_clustering_merged)$Sample))
stopifnot(all(colData(sce_logcounts)$Barcode == colData(sce_clustering_merged)$Barcode))


# select columns in colData to keep

cols_keep_logcounts <- c(1:2, 8:12, 14:15)
colnames(colData(sce_logcounts))[cols_keep_logcounts]

cols_keep_clustering_merged <- c(27, 16, 31, 30, 28, 29, 26)
colnames(colData(sce_clustering_merged))[cols_keep_clustering_merged]


# store combined colData

colData(sce_logcounts) <- cbind(
  colData(sce_logcounts)[, cols_keep_logcounts], 
  colData(sce_clustering_merged)[, cols_keep_clustering_merged])


# check
head(colData(sce_logcounts), 2)


# -----------
# Save object
# -----------

# save SCE object for ExperimentHub

fn_out <- here("processed_data", "SCE", "LC_SingleCellExperiment_EHub")
saveRDS(sce_logcounts, paste0(fn_out, ".rds"))

