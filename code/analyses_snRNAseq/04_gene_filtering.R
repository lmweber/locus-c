#######################################
# LC snRNA-seq analyses: gene filtering
# Lukas Weber, Oct 2022
#######################################


library(here)
library(SingleCellExperiment)


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_logcounts")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)


# --------------------------
# Filter low-expressed genes
# --------------------------

# gene filtering for downstream DE testing steps

# note: keep both protein-coding and non-protein-coding genes

n_umis <- 30
ix_remove <- rowSums(counts(sce)) < n_umis
table(ix_remove)

dim(sce)

sce <- sce[!ix_remove, ]

dim(sce)


# -----------
# Save object
# -----------

fn_out <- here("processed_data", "SCE", "sce_filtered")
saveRDS(sce, paste0(fn_out, ".rds"))

