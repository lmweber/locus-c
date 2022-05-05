########################################################
# LC project
# Script for downstream analyses: BayesSpace evaluations
# Lukas Weber, May 2022
########################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "06_downstream", "BayesSpace")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_BayesSpace.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)


sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# ---------------------
# calculate evaluations
# ---------------------

# evaluate performance of BayesSpace for clustering LC regions
# note: in this dataset we are only interested in the LC regions; ignore other
# clusters in the WM regions, which we assume are mostly noise

# calculate adjusted Rand index (ARI)

