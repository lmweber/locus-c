########################################
# LC analyses: identify SVGs using nnSVG
# Lukas Weber, Jun 2022
########################################

# qrsh -pe local 6 -l mem_free=10G,h_vmem=12G,h_fsize=200G
# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(nnSVG)
library(dplyr)
library(tidyr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "07_nnSVG")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_logcounts.rds")
spe <- readRDS(fn_spe)

dim(spe)

# remove samples where NE neurons were not captured
samples_remove <- "Br5459_LC_round2"
spe <- spe[, !(colData(spe)$sample_id %in% samples_remove)]
colData(spe)$sample_id <- droplevels(colData(spe)$sample_id)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# ---------
# run nnSVG
# ---------

# run nnSVG once per sample within LC annotated regions
# and store lists of top SVGs

res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {
  
  # select sample and LC annotated region
  ix <- colData(spe)$sample_id == sample_ids[s] & colData(spe)$annot_region
  spe_sub <- spe[, ix]
  
  # run nnSVG filtering for mitochondrial gene and low-expressed genes
  spe_sub <- filter_genes(spe_sub)
  
  # run nnSVG
  set.seed(123)
  spe_sub <- nnSVG(spe_sub, n_threads = 6)
  
  # store results
  res_list[[s]] <- rowData(spe_sub)
}


# ---------------
# combine results
# ---------------

# sum gene ranks across samples to generate overall ranking

