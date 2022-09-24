###############################################
# LC Visium analyses: identify SVGs using nnSVG
# Lukas Weber, Sep 2022
###############################################

# qrsh -pe local 10 -l mem_free=5G,h_vmem=6G,h_fsize=200G -now n
# module load conda_R/4.2
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(here)
library(SpatialExperiment)
library(nnSVG)
library(scater)
library(scran)


# directory to save outputs
dir_outputs <- here("outputs", "08_nnSVG")


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


# sample-part IDs for parts that contain LC annotated regions in samples above
sample_part_ids <- c(
  "Br6522_LC_1_round1_single", 
  "Br6522_LC_2_round1_single", 
  "Br8153_LC_round2_left", "Br8153_LC_round2_right", 
  "Br2701_LC_round2_top", "Br2701_LC_round2_bottom", 
  "Br6522_LC_round3_leftbottom", "Br6522_LC_round3_right", 
  "Br8079_LC_round3_left", "Br8079_LC_round3_right", 
  "Br2701_LC_round3_left", "Br2701_LC_round3_right", 
  "Br8153_LC_round3_left"
)
sample_part_ids


# ---------
# run nnSVG
# ---------

# run nnSVG once per sample-part and store lists of top SVGs

res_list <- as.list(rep(NA, length(sample_part_ids)))
names(res_list) <- sample_part_ids

for (s in seq_along(sample_part_ids)) {
  
  # select sample-part
  ix <- colData(spe)$sample_part_id == sample_part_ids[s]
  spe_sub <- spe[, ix]
  
  dim(spe_sub)
  
  # run nnSVG filtering for mitochondrial gene and low-expressed genes
  spe_sub <- filter_genes(
    spe_sub, 
    filter_genes_ncounts = 3, 
    filter_genes_pcspots = 0.5, 
    filter_mito = FALSE
  )
  
  # remove any zeros introduced by filtering
  ix_zeros <- colSums(counts(spe_sub)) == 0
  if (sum(ix_zeros) > 0) {
    spe_sub <- spe_sub[, !ix_zeros]
  }
  
  dim(spe_sub)
  
  # re-calculate logcounts after filtering
  spe_sub <- computeLibraryFactors(spe_sub)
  spe_sub <- logNormCounts(spe_sub)
  
  # run nnSVG
  set.seed(123)
  spe_sub <- nnSVG(spe_sub, n_threads = 10)
  
  # store results
  res_list[[s]] <- rowData(spe_sub)
}


# ------------
# save results
# ------------

# save nnSVG results
fn_out <- file.path(dir_outputs, "LC_nnSVG_results")
saveRDS(res_list, paste0(fn_out, ".rds"))

