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
library(tibble)
library(ggplot2)
library(ggnewscale)
library(ggrepel)


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

# pairwise DE testing between neuronal clusters


# store total UMI counts per gene
rowData(sce)$sum_gene <- rowSums(counts(sce))


# select neuronal clusters only (excluding ambiguous)
clus_neurons <- c(
  29, ## excitatory
  26, 17, 14, 1, 8, 7, 24, 18,  ## inhibitory
  6,  ## NE
  16  ## 5-HT
)
ix_neurons <- colData(sce)$label %in% clus_neurons
table(ix_neurons)

sce <- sce[, ix_neurons]
dim(sce)

# remove empty levels
colData(sce)$label <- droplevels(colData(sce)$label)


# calculate DE tests
# note: not blocking by sample since 1 out of 3 samples contains almost zero NE neurons
marker_info <- scoreMarkers(
  sce, 
  groups = colData(sce)$label, 
  row.data = rowData(sce)[, c("gene_id", "gene_name", "sum_gene")]
)
# marker_info <- findMarkers(
#   sce,
#   groups = colData(sce)$label,
#   lfc = 1,
#   direction = "up",
#   row.data = rowData(sce)[, c("gene_id", "gene_name", "sum_gene")]
# )

marker_info

marker_info[["6"]]
#View(as.data.frame(marker_info[["6"]]))
marker_info[["16"]]
#View(as.data.frame(marker_info[["16"]]))


# check some known genes

ix_known <- which(marker_info[["6"]]$gene_name %in% c("TH", "SLC6A2", "DBH"))
marker_info[["6"]][ix_known, ]

ix_known <- which(marker_info[["16"]]$gene_name %in% c("TPH2", "SLC6A4"))
marker_info[["16"]][ix_known, ]


# significant DE genes
# table(marker_info$NE$FDR < 1)
# table(marker_info$NE$FDR < 0.05)
# table(marker_info$NE$FDR < 1e-6)
# table(marker_info$NE$FDR < 1e-100)


# using mean Cohen's d statistic (i.e. standardized log-fold-change) across
# pairwise comparisons as ranking metric (from OSCA chapter 6)

ordered <- marker_info[["6"]][order(marker_info[["6"]]$mean.logFC.cohen, decreasing = TRUE), ]
head(ordered)
#View(as.data.frame(ordered))

ordered <- marker_info[["16"]][order(marker_info[["16"]]$mean.logFC.cohen, decreasing = TRUE), ]
head(ordered)
#View(as.data.frame(ordered))

