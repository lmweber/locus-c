###################################
# LC snRNA-seq analyses: DE testing
# Lukas Weber, Sep 2022
###################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)


dir_plots <- here("plots", "snRNAseq", "07_DEtesting")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_clustering_merged")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)

# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# ----------
# DE testing
# ----------

# DE testing between NE neurons and all other nuclei

# create labels for NE neurons vs. all other nuclei
label_NE <- 
  ifelse(colData(sce)$label_merged == "NE", "NE", "other") %>% 
  factor(., levels = c("NE", "other"))

table(label_NE)

colData(sce)$label_NE <- label_NE


# store gene names in rownames for summary tables
rownames(sce) <- rowData(sce)$gene_name


# calculate DE tests
# note: not blocking by sample since 1 out of 3 samples contains almost zero NE neurons

marker_info <- findMarkers(sce, groups = colData(sce)$label_NE, lfc = 1)

marker_info$NE

table(marker_info$NE$FDR < 0.05)

sig <- marker_info$NE$FDR < 0.05
rownames(marker_info$NE[sig, ])

