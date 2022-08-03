##########################################
# LC analyses: normalization and logcounts
# Lukas Weber, Jun 2022
##########################################

# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)

# directory to save plots
dir_plots <- here("plots", "02_quality_control")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_filtered.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)


# ---------------------------
# normalization and logcounts
# ---------------------------

# normalization and log-transformed normalized counts (logcounts)
# using code from OSCA

# compare deconvolution size factors and library size factors
# note this is a multi-sample dataset

# deconvolution size factors
set.seed(100)
qclus <- quickCluster(spe)
table(qclus)
sf_deconv <- calculateSumFactors(spe, cluster = qclus)

# library size factors
sf_lib <- librarySizeFactors(spe)


stopifnot(length(sf_deconv) == length(sf_lib))


# plot comparing size factors

fn <- file.path(dir_plots, "LC_sizeFactors")

pdf(paste0(fn, ".pdf"), width = 6, height = 6)
plot(x = sf_lib, y = sf_deconv, pch = 16, cex = 0.3, col = qclus, log = "xy", 
     xlab = "Library size factors", ylab = "Deconvolution factors")
dev.off()

png(paste0(fn, ".png"), width = 6 * 200, height = 6 * 200, res = 200)
plot(x = sf_lib, y = sf_deconv, pch = 16, cex = 0.3, col = qclus, log = "xy", 
     xlab = "Library size factors", ylab = "Deconvolution factors")
dev.off()


# use library size factors and calculate logcounts

spe <- computeLibraryFactors(spe)
stopifnot(all(sizeFactors(spe) == sf_lib))

spe <- logNormCounts(spe)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LC_logcounts")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

