#################################################################
# LC project
# Script for downstream analyses: differential expression testing
# Lukas Weber, May 2022
#################################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)
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

# note: using sample_id, not sample_part_id

table(colData(spe)$sample_id, colData(spe)$annot_region)

annot_region_fctr <- factor(as.numeric(colData(spe)$annot_region), labels = c("WM", "LC"))
table(annot_region_fctr)

ids <- DataFrame(
  annot_region_pseudo = annot_region_fctr, 
  sample_id_pseudo = colData(spe)$sample_id
)

# pseudobulk
spe_pseudo <- aggregateAcrossCells(spe, ids)

# add levels to colData and column names
levs <- paste(
  colData(spe_pseudo)$annot_region_pseudo, 
  colData(spe_pseudo)$sample_id_pseudo, 
  sep = "."
)
colData(spe_pseudo)$level_id <- levs
colnames(spe_pseudo) <- levs

# check
head(colData(spe_pseudo), 3)
counts(spe_pseudo)[1:3, ]


# recalculate logcounts
# using default library size scale factors
spe_pseudo <- logNormCounts(spe_pseudo, size.factors = NULL)

