########################################################################
# LC snRNA-seq analyses: cluster identification using known marker genes
# Lukas Weber, Sep 2022
# using code by Matthew N Tran
########################################################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(jaffelab)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)


dir_plots <- here("plots", "snRNAseq", "06_cluster_identification", "DE_testing")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_clustering")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)

# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)

