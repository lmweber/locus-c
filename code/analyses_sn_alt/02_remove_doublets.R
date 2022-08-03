########################################
# LC snRNA-seq analyses: remove doublets
# Lukas Weber, July 2022
########################################

# run in interactive session on JHPCE

# screen -S LC
# qrsh -l mem_free=20G,h_vmem=22G,h_fsize=200G -now n
# cd /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn_alt
# module load conda_R/devel
# R


library(SingleCellExperiment)
library(scDblFinder)
library(here)


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE_alt", "sce_raw")
sce <- readRDS(paste0(fn, ".rds"))

table(colData(sce)$Sample)


# ---------------
# Remove doublets
# ---------------

# identify and remove doublets using scDblFinder

# note: no random seed required

sce <- scDblFinder(sce, samples = "Sample")

head(colData(sce), 3)


# number of doublets per sample
table(colData(sce)$Sample, colData(sce)$scDblFinder.class)


# remove doublets
ix_dbl <- colData(sce)$scDblFinder.class == "doublet"
table(ix_dbl)

sce <- sce[, !ix_dbl]

dim(sce)


# -----------
# Save object
# -----------

fn_out <- here("processed_data", "SCE_alt", "sce_doubletsRemoved")
saveRDS(sce, paste0(fn_out, ".rds"))

