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

# pairwise DE testing between merged clusters


# store total UMI counts per gene
rowData(sce)$sum_gene <- rowSums(counts(sce))


# remove ambiguous neurons cluster
ix_ambiguous <- colData(sce)$label_merged == "neurons_ambiguous"
table(ix_ambiguous)

sce <- sce[, !ix_ambiguous]
dim(sce)

# remove empty level
colData(sce)$label_merged <- droplevels(colData(sce)$label_merged)


# calculate DE tests
# note: not blocking by sample since 1 out of 3 samples contains almost zero NE neurons
marker_info <- scoreMarkers(
  sce, 
  groups = colData(sce)$label_merged, 
  row.data = rowData(sce)[, c("gene_id", "gene_name", "sum_gene")]
)
# marker_info <- findMarkers(
#   sce, 
#   groups = colData(sce)$label_merged, 
#   lfc = 1, 
#   direction = "up", 
#   row.data = rowData(sce)[, c("gene_id", "gene_name", "sum_gene")]
# )

marker_info$NE
marker_info$`5HT`


# check some known genes

ix_known <- which(marker_info$NE$gene_name %in% c("TH", "SLC6A2", "DBH"))
marker_info$NE[ix_known, ]

ix_known <- which(marker_info$NE$gene_name %in% c("TPH2", "SLC6A4"))
marker_info$`5HT`[ix_known, ]


# significant DE genes
# table(marker_info$NE$FDR < 1)
# table(marker_info$NE$FDR < 0.05)
# table(marker_info$NE$FDR < 1e-6)
# table(marker_info$NE$FDR < 1e-100)


# using mean Cohen's d statistic (i.e. standardized log-fold-change) across
# pairwise comparisons as ranking metric (from OSCA chapter 6)

ordered <- marker_info$NE[order(marker_info$NE$mean.logFC.cohen, decreasing = TRUE), ]
head(ordered)
#View(as.data.frame(ordered))

ordered <- marker_info$`5HT`[order(marker_info$`5HT`$mean.logFC.cohen, decreasing = TRUE), ]
head(ordered)
#View(as.data.frame(ordered))

