##########################################################
# LC snRNA-seq analyses: Mulvey et al. (2018) rodent genes
# Lukas Weber, Sep 2022
##########################################################


library(here)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)


# directory to save plots
dir_plots <- here("plots", "snRNAseq", "08_rodent_genes")


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


# --------------------
# load mouse gene list
# --------------------

# list of 45 mouse genes for LC from Mulvey et al. (2018), Figure 2A

fn <- here("inputs", "Mulvey_markers", "Mulvey2018_Fig2A_markers.txt")
mouse_genes <- read.table(fn)[, 1]

length(mouse_genes)
mouse_genes


# converted to human genes using biomaRt (code below) and by searching for genes
# individually at: https://www.ncbi.nlm.nih.gov/gene/

human_genes <- c(
  "AGTR1", "ASB4", "CALCR", "CALR3", "CHODL", "CHRNA6", "CILP", "CYB561", 
  "DBH", "DDC", "DLK1", "ADGRE1", "EYA2", "SHISAL2B", "FAM183A", "FIBCD1", 
  "GAL", "GCH1", "GLRA2", "GNG4", "GPX3", "GTF2A1L", "HCRTR1", "IGSF5", "MAOA", 
  "MRAP2", "MYOM2", "NEUROG2", "SLC9B2", "NXPH4", "OVGP1", "PCBD1", "PHOX2A", 
  "PHOX2B", "PLA2G4D", "PTGER2", "SLC18A2", "SLC31A1", "SLC6A2", "STBD1", 
  "SYT17", "TH", "TM4SF1", "TM4SF5", "TRAF3IP2")

human_genes <- sort(human_genes)
length(human_genes)

# 42 out of 45 genes present in SCE object
sum(human_genes %in% rowData(sce)$gene_name)

# keep genes that are present in SCE object
ix_keep <- which(human_genes %in% rowData(sce)$gene_name)
human_genes <- human_genes[ix_keep]

stopifnot(length(human_genes) == sum(human_genes %in% rowData(sce)$gene_name))


# ------------
# plot heatmap
# ------------

# function from rafalib package
splitit <- function(x) split(seq(along = x), x)


# calculate matrix of values: genes x clusters

cell_idx <- splitit(sce$label)

dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$gene_name

mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[human_genes, i]))))


# number of nuclei per cluster
n <- table(colLabels(sce))


# row annotation
row_ha <- rowAnnotation(
  n = anno_barplot(as.numeric(n), gp = gpar(fill = "navy"), border = FALSE), 
  show_annotation_name = FALSE)


hm <- Heatmap(
  mat, 
  name = "mean\nlogcounts", 
  column_title = "Mulvey et al. genes", 
  column_title_gp = gpar(fontface = "bold"), 
  col = brewer.pal(n = 7, "OrRd"), 
  right_annotation = row_ha, 
  #row_order = cluster_pops_order, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  #row_split = cluster_pops_rev, 
  row_title = NULL, 
  column_names_gp = gpar(fontface = "italic"), 
  rect_gp = gpar(col = "gray50", lwd = 0.5))

hm


# save heatmap
fn <- file.path(dir_plots, "Mulvey_genes_clusters")

pdf(paste0(fn, ".pdf"), width = 8, height = 6.5)
hm
dev.off()

png(paste0(fn, ".png"), width = 8 * 200, height = 6.5 * 200, res = 200)
hm
dev.off()

