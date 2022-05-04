#################################################################
# LC project
# Script for downstream analyses: differential expression testing
# Lukas Weber, Apr 2022
#################################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(BayesSpace)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "06_downstream", "DEtesting")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrected.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# ----------------
# pseudobulk spots
# ----------------

# pseudobulk spots within LC regions vs. WM regions

