###############################################
# LC Visium analyses: dataset for ExperimentHub
# Lukas Weber, Oct 2022
###############################################

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


# --------------------------------
# clean dataset for public release
# --------------------------------

# remove unused columns from colData

# not using 'cell_count' (number of cells per spot)
# not using 'high_subsets_mito_percent' (QC filtering on mitochondrial percentage per spot)

ix_cols_remove <- c(8, 21)

colnames(colData(spe))[ix_cols_remove]

# remove columns
colData(spe) <- colData(spe)[, -ix_cols_remove]

# check
head(colData(spe), 2)


# ---------------
# save SPE object
# ---------------

fn_out <- here("processed_data", "SPE", "LC_SpatialExperiment_EHub")
saveRDS(spe, paste0(fn_out, ".rds"))

