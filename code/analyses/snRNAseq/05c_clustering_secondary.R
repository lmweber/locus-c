###################################################################
# LC snRNA-seq analyses: secondary clustering of inhibitory neurons
# Lukas Weber, Jun 2023
###################################################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(bluster)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


dir_plots <- here("plots", "singleNucleus", "05c_clustering_inhibitory")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_clustering_merged")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)


# -------------------------
# Subset inhibitory neurons
# -------------------------

sce_full <- sce

# select inhibitory neuron clusters
# (identified based on marker expression heatmap from previous script)
clus_select <- c(24, 25, 14, 4, 8, 20, 17)

ix_select <- colLabels(sce) %in% clus_select
table(ix_select)

sce <- sce[, ix_select]

dim(sce)


# ------------------------------
# Re-calculate logcounts and PCA
# ------------------------------

# re-calculate logcounts and PCA on subset to ensure PCs capture the most relevant variation

# normalization using library size factors (instead of by deconvolution) since
# there is no major population variation within this subset
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)


# feature selection (within subset)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$gene_name)
table(is_mito)
rowData(sce)$gene_name[is_mito]

# fit mean-variance relationship
dec <- modelGeneVar(sce, subset.row = !is_mito)
# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)


# dimensionality reduction (within subset)

# note: not applying any batch integration (for consistency with main clustering)

# calculate PCA (on top HVGs)
set.seed(123)
sce <- runPCA(sce, subset_row = top_hvgs)

# calculate UMAP (on top PCs)
set.seed(123)
sce <- runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) <- paste0("UMAP", seq_len(ncol(reducedDim(sce, "UMAP"))))

reducedDimNames(sce)


# --------------------
# Secondary clustering
# --------------------

# secondary clustering of inhibitory neurons

# clustering algorithm and parameters from OSCA
# two-stage clustering algorithm using high-resolution k-means and graph-based clustering

set.seed(123)
clus <- clusterCells(
  sce, 
  use.dimred = "PCA", 
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 1000), 
    second = NNGraphParam(k = 10)
  )
)
colLabels(sce) <- clus


# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# expression of key markers for NE, 5-HT, and cholinergic neuron populations
ix <- c(
  DBH = which(rowData(sce)$gene_name == "DBH"), 
  TH = which(rowData(sce)$gene_name == "TH"), 
  SLC6A2 = which(rowData(sce)$gene_name == "SLC6A2"), 
  TPH2 = which(rowData(sce)$gene_name == "TPH2"), 
  SLC6A4 = which(rowData(sce)$gene_name == "SLC6A4"), 
  SLC5A7 = which(rowData(sce)$gene_name == "SLC5A7"), 
  CHAT = which(rowData(sce)$gene_name == "CHAT"), 
  ACHE = which(rowData(sce)$gene_name == "ACHE"), 
  BCHE = which(rowData(sce)$gene_name == "BCHE"), 
  SLC18A3 = which(rowData(sce)$gene_name == "SLC18A3"), 
  PRIMA1 = which(rowData(sce)$gene_name == "PRIMA1")
)

n_clus <- length(table(colLabels(sce)))

res_list <- list()
for (k in seq_len(n_clus)) {
  res_list[[k]] <- rowMeans(logcounts(sce)[ix, colLabels(sce) == k])
}
res_mat <- do.call("rbind", res_list)
rownames(res_mat) <- seq_len(n_clus)
colnames(res_mat) <- names(ix)

cbind(
  n = table(colLabels(sce)), 
  table(colLabels(sce), colData(sce)$Sample), 
  res_mat
)


# ------------
# Marker genes
# ------------

# marker genes for inhibitory neuron subpopulations

markers_inhib <- c(
  # neuron markers
  "SNAP25", "SYT1", 
  # inhibitory (GABAergic) neuron markers
  "GAD1", "GAD2", 
  # inhibitory subpopulations (from Keri Martinowich 2022-07-22)
  "SST", "KIT", "CALB1", "CALB2", "TAC1", "CNR1", "PVALB", "CORT", "VIP", "NPY", "CRHBP", "CCK" ## "HTR3A" (not present in data)
)

# second set of GABAergic neuron marker genes from Luskin and Li et al. (2022)
# and GAD1, GAD2

genes_gabaergic_luskin <- c("GAD1", "GAD2", 
                            "AGRP", "CALCA", "CCK", "GAL", "CARTPT", "NMB", 
                            "ADCYAP1", "BDNF", "PCSK1", "PCSK2", "PCSK1N", 
                            "PDYN", "PENK", "PNOC", "SST", "POMC", "TAC1")


# ----------------------------------------------
# Marker expression heatmap: inhibitory subtypes
# ----------------------------------------------

# marker labels
marker_labels_inhib <- c(
  rep("neuron", 2), 
  rep("inhibitory", 14))

marker_labels_inhib <- 
  factor(marker_labels_inhib, levels = unique(marker_labels_inhib))


# colors
colors_markers_inhib <- list(marker = c(
  neuron = "black", 
  inhibitory = "#AEC7E8"))


# number of nuclei per cluster
n <- table(colLabels(sce))


# heatmap data

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$label)
dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$gene_name


hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_inhib, i]))))


# row annotation
row_ha <- rowAnnotation(
  n = anno_barplot(as.numeric(n), gp = gpar(fill = "navy"), border = FALSE),
  show_annotation_name = FALSE)

# column annotation
col_ha <- columnAnnotation(
  marker = marker_labels_inhib, 
  show_annotation_name = FALSE, 
  show_legend = FALSE, 
  col = colors_markers_inhib)

# check range of values to select consistent range across plots
range(hm_mat)

hm <- Heatmap(
  hm_mat, 
  name = "mean\nlogcounts", 
  column_title = "LC secondary clustering mean marker expression", 
  column_title_gp = gpar(fontface = "bold"), 
  col = colorRamp2(seq(0, 6.4, length.out = 7), brewer.pal(n = 7, "OrRd")), 
  right_annotation = row_ha, 
  bottom_annotation = col_ha, 
  cluster_rows = TRUE, 
  cluster_columns = FALSE, 
  row_title = NULL, 
  column_split = marker_labels_inhib, 
  column_names_gp = gpar(fontface = "italic"), 
  rect_gp = gpar(col = "gray50", lwd = 0.5))

hm


# save heatmap
fn <- file.path(dir_plots, "clustering_inhibitory_subtypes")

pdf(paste0(fn, ".pdf"), width = 5.5, height = 4)
hm
dev.off()

png(paste0(fn, ".png"), width = 5.5 * 200, height = 4 * 200, res = 200)
hm
dev.off()


# --------------------------------------------------------------
# Marker expression heatmap:
# GABAergic neuron marker genes from Luskin and Li et al. (2022)
# --------------------------------------------------------------

# marker labels
marker_labels_luskin <- as.factor(rep("inhibitory", length(genes_gabaergic_luskin)))

# colors
colors_markers_luskin <- list(marker = c(inhibitory = "#AEC7E8"))

# number of nuclei per cluster
n <- table(colLabels(sce))


# heatmap data

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$label)
dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$gene_name


hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[genes_gabaergic_luskin, i]))))


# row annotation
row_ha <- rowAnnotation(
  n = anno_barplot(as.numeric(n), gp = gpar(fill = "navy"), border = FALSE),
  show_annotation_name = FALSE)

# column annotation
col_ha <- columnAnnotation(
  marker = marker_labels_luskin, 
  show_annotation_name = FALSE, 
  show_legend = FALSE, 
  col = colors_markers_luskin)

# check range of values to select consistent range across plots
range(hm_mat)

hm <- Heatmap(
  hm_mat, 
  name = "mean\nlogcounts", 
  column_title = "LC secondary clustering mean marker expression", 
  column_title_gp = gpar(fontface = "bold"), 
  col = colorRamp2(seq(0, 6.4, length.out = 7), brewer.pal(n = 7, "OrRd")), 
  right_annotation = row_ha, 
  bottom_annotation = col_ha, 
  cluster_rows = TRUE, 
  cluster_columns = FALSE, 
  row_title = NULL, 
  column_split = marker_labels_luskin, 
  column_names_gp = gpar(fontface = "italic"), 
  rect_gp = gpar(col = "gray50", lwd = 0.5))

hm


# save heatmap
fn <- file.path(dir_plots, "clustering_inhibitory_subtypes_luskin_genes")

pdf(paste0(fn, ".pdf"), width = 5.75, height = 4)
hm
dev.off()

png(paste0(fn, ".png"), width = 5.75 * 200, height = 4 * 200, res = 200)
hm
dev.off()


# --------------------
# Store cluster labels
# --------------------

# store secondary clustering labels in full SCE object

# check unique barcode IDs
table(duplicated(colnames(sce)))
table(duplicated(colnames(sce_full)))
table(colnames(sce) %in% colnames(sce_full))

# match and store cluster labels
clus_inhibitory <- rep(NA, ncol(sce_full))
names(clus_inhibitory) <- colnames(sce_full)
clus_inhibitory[colnames(sce)] <- colData(sce)$label

colData(sce_full)$label_inhibitory <- clus_inhibitory

# check
table(colData(sce_full)$label)
table(colData(sce_full)$label_inhibitory)
table(colData(sce_full)$label_inhibitory, useNA = "always")


# -----------
# Save object
# -----------

# note saving 'sce_full' object

fn_out <- here("processed_data", "SCE", "sce_clustering_secondary")
saveRDS(sce_full, paste0(fn_out, ".rds"))

