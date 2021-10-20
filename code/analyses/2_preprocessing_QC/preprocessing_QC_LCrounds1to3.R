#################################
# LC project
# Script for preprocessing and QC
# Lukas Weber, Oct 2021
#################################

# scran QC and splitting Visium slides into individual samples


# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)


# ---------
# load data
# ---------

# load previously built SpatialExperiment object

fn_spe <- here("processed_data", "SPE", "LCrounds1to3_SPE_raw.rds")
spe <- readRDS(fn_spe)


# -------------------------------
# spot-level quality control (QC)
# -------------------------------

# apply spot-level QC across all samples

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$symbol)
table(is_mito)
rowData(spe)$symbol[is_mito]

# calculate QC metrics using scater package
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

# plot histograms of QC metrics
par(mfrow = c(1, 4))
hist(colData(spe)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe)$detected, xlab = "detected", main = "Genes per spot")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mitochondrial", main = "Percent mito UMIs")
hist(colData(spe)$count, xlab = "number of cells", main = "No. cells per spot")
par(mfrow = c(1, 1))

# keep all spots for now
colData(spe)$discard <- FALSE

colData(spe)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LCrounds1to3_SPE_preprocessed.rds")
saveRDS(spe, fn_out)

