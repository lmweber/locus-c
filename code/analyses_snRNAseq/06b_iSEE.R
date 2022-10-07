############################################
# LC snRNA-seq analyses: SCE object for iSEE
# Lukas Weber, Oct 2022
############################################


library(here)
library(SingleCellExperiment)


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "LC_singleNucleus_SCE_EHub")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)


# ----------
# Add colors
# ----------

# add colors for iSEE app


# -----------
# Save object
# -----------

# save SCE object for iSEE

fn_out <- here("processed_data", "SCE", "LC_singleNucleus_SCE_iSEE")
saveRDS(sce_logcounts, paste0(fn_out, ".rds"))

