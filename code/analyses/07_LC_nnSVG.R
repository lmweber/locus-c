########################################
# LC analyses: identify SVGs using nnSVG
# Lukas Weber, Jun 2022
########################################

# qrsh -pe local 10 -l mem_free=6G,h_vmem=7G,h_fsize=200G
# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)
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


# sample-part IDs for parts that contain LC annotated regions
sample_part_ids <- c(
  "Br6522_LC_1_round1_single", 
  "Br6522_LC_2_round1_single", 
  "Br8153_LC_round2_left", "Br8153_LC_round2_right", 
  "Br2701_LC_round2_bottom", "Br2701_LC_round2_top", 
  "Br6522_LC_round3_left", "Br6522_LC_round3_right", 
  "Br8079_LC_round3_left", "Br8079_LC_round3_right", 
  "Br2701_LC_round3_left", "Br2701_LC_round3_right", 
  "Br8153_LC_round3_left"
)
sample_part_ids


# ---------
# run nnSVG
# ---------

# run nnSVG once per sample-part within LC annotated regions
# and store lists of top SVGs

res_list <- as.list(rep(NA, length(sample_part_ids)))
names(res_list) <- sample_part_ids

for (s in seq_along(sample_part_ids)) {
  
  # select sample-part and LC annotated region
  ix <- colData(spe)$sample_part_id == sample_part_ids[s] & colData(spe)$annot_region
  spe_sub <- spe[, ix]
  
  # run nnSVG filtering for mitochondrial gene and low-expressed genes
  spe_sub <- filter_genes(
    spe_sub, 
    filter_genes_ncounts = 2, 
    filter_genes_pcspots = 1
  )
  
  # re-calculate logcounts after filtering
  spe_sub <- logNormCounts(spe_sub)
  
  # run nnSVG
  set.seed(123)
  spe_sub <- nnSVG(spe_sub, n_threads = 10)
  
  # store results
  res_list[[s]] <- rowData(spe_sub)
}


# ---------------
# combine results
# ---------------

# sum gene ranks across sample-parts to generate overall ranking

# match results from each sample-part and store in correct rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_part_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_part_ids

for (s in seq_along(sample_part_ids)) {
  stopifnot(colnames(res_ranks)[s] == sample_part_ids[s])
  stopifnot(colnames(res_ranks)[s] == names(res_list)[s])
  
  rownames_s <- rownames(res_list[[s]])
  res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}

# keep only genes that were not filtered out in all sample-parts
res_ranks <- na.omit(res_ranks)

# calculate average ranks
avg_ranks <- sort(rowMeans(res_ranks))

# summary table
df_summary <- data.frame(
  gene_id = names(avg_ranks), 
  gene_name = rowData(spe)[names(avg_ranks), "gene_name"], 
  gene_type = rowData(spe)[names(avg_ranks), "gene_type"], 
  avg_rank = unname(avg_ranks), 
  row.names = names(avg_ranks)
)

