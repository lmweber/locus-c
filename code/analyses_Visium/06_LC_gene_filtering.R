####################################
# LC Visium analyses: gene filtering
# Lukas Weber, Oct 2022
####################################

# module load conda_R/4.2
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(here)
library(SpatialExperiment)


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_logcounts.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)


# --------------------------
# filter low-expressed genes
# --------------------------

# filter out genes with extremely low expression
# using simple threshold on total UMI counts summed across all spots

n_umis <- 80
ix_low_genes <- rowSums(counts(spe)) < n_umis
table(ix_low_genes)

spe <- spe[!ix_low_genes, ]
dim(spe)


# ------------
# filter zeros
# ------------

# new zeros may have been created after spot-level QC and filtering low-expressed genes

# remove genes with zero expression
ix_zero_genes <- rowSums(counts(spe)) == 0
table(ix_zero_genes)

if (sum(ix_zero_genes) > 0) {
  spe <- spe[!ix_zero_genes, ]
}

dim(spe)

# remove spots with zero expression
ix_zero_spots <- colSums(counts(spe)) == 0
table(ix_zero_spots)

if (sum(ix_zero_spots) > 0) {
  spe <- spe[, !ix_zero_spots]
}

dim(spe)


# check no zeros or NAs remaining
table(colData(spe)$in_tissue, useNA = "always")
table(rowSums(counts(spe)) == 0, useNA = "always")
table(colSums(counts(spe)) == 0, useNA = "always")


# ---------------
# save SPE object
# ---------------

fn_out <- here("processed_data", "SPE", "LC_filtered")
saveRDS(spe, paste0(fn_out, ".rds"))

