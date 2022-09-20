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
library(ComplexHeatmap)


dir_plots <- here("plots", "snRNAseq", "07_DEtesting")
dir_outputs <- here("outputs", "snRNAseq", "07_DEtesting")


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


# select neuronal clusters (excluding ambiguous) to test against for NE neurons
# excluding 5-HT since these are too similar to NE
clus_neurons_NE <- c(
  29, 21, 20, 23,  ## excitatory
  26, 17, 14, 1, 8, 7, 24, 18,  ## inhibitory
  6  ## NE
)
ix_neurons_NE <- colData(sce)$label %in% clus_neurons_NE
table(ix_neurons_NE)

sce_NE <- sce[, ix_neurons_NE]
dim(sce_NE)
# remove empty levels
colData(sce_NE)$label <- droplevels(colData(sce_NE)$label)


# select neuronal clusters (excluding ambiguous) to test against for 5-HT neurons
# excluding NE since these are too similar to 5-HT
clus_neurons_5HT <- c(
  29, 21, 20, 23,  ## excitatory
  26, 17, 14, 1, 8, 7, 24, 18,  ## inhibitory
  16  ## 5-HT
)
ix_neurons_5HT <- colData(sce)$label %in% clus_neurons_5HT
table(ix_neurons_5HT)

sce_5HT <- sce[, ix_neurons_5HT]
dim(sce_5HT)
# remove empty levels
colData(sce_5HT)$label <- droplevels(colData(sce_5HT)$label)


# calculate DE tests
# note: not blocking by sample since 1 out of 3 samples contains almost zero NE neurons
# testing for genes with log-fold-changes significantly greater than 1 (lfc = 1, direction = "up")

marker_info_NE <- findMarkers(
  sce_NE, 
  groups = colData(sce_NE)$label, 
  lfc = 1, 
  direction = "up", 
  row.data = rowData(sce_NE)[, c("gene_id", "gene_name", "sum_gene")], 
  add.summary = TRUE
)

marker_info_5HT <- findMarkers(
  sce_5HT, 
  groups = colData(sce_5HT)$label, 
  lfc = 1, 
  direction = "up", 
  row.data = rowData(sce_5HT)[, c("gene_id", "gene_name", "sum_gene")], 
  add.summary = TRUE
)

marker_info_NE
marker_info_5HT


# NE neuron cluster
marker_info_NE[["6"]]

# 5-HT neuron cluster
marker_info_5HT[["16"]]


# check some known genes

ix_known <- which(marker_info_NE[["6"]]$gene_name %in% c("TH", "SLC6A2", "DBH"))
marker_info_NE[["6"]][ix_known, ]

ix_known <- which(marker_info_5HT[["16"]]$gene_name %in% c("TPH2", "SLC6A4"))
marker_info_5HT[["16"]][ix_known, ]


# significant DE genes

table(marker_info_NE[["6"]]$FDR < 1)
table(marker_info_NE[["6"]]$FDR < 1e-100)

table(marker_info_5HT[["16"]]$FDR < 1)
table(marker_info_5HT[["16"]]$FDR < 1e-100)


# -------------------------
# Volcano plots: NE neurons
# -------------------------

# NE neurons

fdr <- marker_info_NE[["6"]]$FDR
logfc <- marker_info_NE[["6"]]$summary.logFC

names(fdr) <- names(logfc) <- marker_info_NE[["6"]]$gene_name


# select significant genes

thresh_fdr <- 0.05
thresh_logfc <- log2(2)
sig <- (fdr < thresh_fdr) & (logfc > thresh_logfc)
table(sig)

thresh_fdr <- 1e-20
thresh_logfc <- log2(4)
highlysig <- (fdr < thresh_fdr) & (logfc > thresh_logfc)
table(highlysig)


df <- data.frame(
  gene = names(fdr), 
  FDR = fdr, 
  log2FC = logfc, 
  highlysig = highlysig
)

pal <- c("black", "red")


# volcano plot
set.seed(123)
ggplot(df, aes(x = log2FC, y = -log10(FDR), color = highlysig, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$highlysig, ], size = 0.5) + 
  #geom_text_repel(data = df[df$highlysig, ], 
  #                size = 1.5, nudge_y = 0.1, 
  #                force = 0.1, force_pull = 0.1, min.segment.length = 0.1, 
  #                max.overlaps = 20) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(thresh_fdr), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = thresh_logfc, lty = "dashed", color = "royalblue") + 
  ggtitle("NE neuron cluster vs. all neuronal clusters") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "DEtesting_NEneuronsVsOtherNeuronal")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)


# -------------------
# Heatmap: NE neurons
# -------------------

# NE neurons

hmat <- marker_info_NE[["6"]][, c("gene_name", "self.average", "other.average", "FDR", "summary.logFC")]

# select significant
sig <- with(hmat, FDR < 0.05 & summary.logFC > 1)
hmat <- hmat[sig, ]

# order by FDR and select top n for plot
hmat <- hmat[order(hmat$FDR), ]

# format gene names and FDRs in row names
nms <- with(hmat, paste0(gene_name, " (", format(signif(FDR, 2)), ")"))
rownames(hmat) <- nms

hmat <- as.matrix(hmat[, c("self.average", "other.average")])
colnames(hmat) <- c("NE", "other")

# remove mitochondrial genes from heatmap
ix_mito <- grepl("^MT-", rownames(hmat))
table(ix_mito)
hmat <- hmat[!ix_mito, ]
dim(hmat)

# select top n
hmat <- hmat[1:100, ]

# create heatmap
hm <- Heatmap(
  hmat, 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  column_names_rot = 0, column_names_gp = gpar(fontsize = 10), column_names_centered = TRUE, 
  row_names_gp = gpar(fontsize = 9, fontface = "italic"), 
  column_title = "NE vs. other\nneuronal clusters", 
  column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  name = "mean\nlogcounts"
)

hm

# save heatmap
fn <- file.path(dir_plots, "DEtesting_heatmap_NEvsOtherNeuronal")

pdf(paste0(fn, ".pdf"), width = 3.75, height = 12)
hm
dev.off()

png(paste0(fn, ".png"), width = 3.75 * 200, height = 12 * 200, res = 200)
hm
dev.off()


# -----------------------
# Spreadsheet: NE neurons
# -----------------------

# save spreadsheet

cols <- c("gene_id", "gene_name", "sum_gene", "p.value", "FDR", "summary.logFC")
#cols <- c("gene_id", "gene_name", "sum_gene", "self.average", "other.average", "p.value", "FDR", "summary.logFC")
df <- marker_info_NE[["6"]][, cols]
colnames(df)[c(4, 6)] <- c("p_value", "logFC")

# select significant
sig <- with(df, FDR < 0.05 & logFC > 1)
table(sig)
df <- df[sig, ]

# order by FDR
df <- df[order(df$FDR), ]

df <- as.data.frame(df)
rownames(df) <- NULL


# save .csv file
fn <- file.path(dir_outputs, "DEtesting_NEvsOtherNeuronal.csv")
write.csv(df, file = fn, row.names = FALSE)


# ---------------------------
# Volcano plots: 5-HT neurons
# ---------------------------

# 5-HT neurons

fdr <- marker_info_5HT[["16"]]$FDR
logfc <- marker_info_5HT[["16"]]$summary.logFC

names(fdr) <- names(logfc) <- marker_info_5HT[["16"]]$gene_name


# select significant genes

thresh_fdr <- 0.05
thresh_logfc <- log2(2)
sig <- (fdr < thresh_fdr) & (logfc > thresh_logfc)
table(sig)

thresh_fdr <- 1e-20
thresh_logfc <- log2(4)
highlysig <- (fdr < thresh_fdr) & (logfc > thresh_logfc)
table(highlysig)


df <- data.frame(
  gene = names(fdr), 
  FDR = fdr, 
  log2FC = logfc, 
  highlysig = highlysig
)

pal <- c("black", "red")


# volcano plot
set.seed(123)
ggplot(df, aes(x = log2FC, y = -log10(FDR), color = highlysig, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$highlysig, ], size = 0.5) + 
  #geom_text_repel(data = df[df$highlysig, ], 
  #                size = 1.5, nudge_y = 0.1, 
  #                force = 0.1, force_pull = 0.1, min.segment.length = 0.1, 
  #                max.overlaps = 20) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(thresh_fdr), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = thresh_logfc, lty = "dashed", color = "royalblue") + 
  ggtitle("5-HT neuron cluster vs. all neuronal clusters") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "DEtesting_5HTneuronsVsOtherNeuronal")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)


# ---------------------
# Heatmap: 5-HT neurons
# ---------------------

# 5-HT neurons

hmat <- marker_info_5HT[["16"]][, c("gene_name", "self.average", "other.average", "FDR", "summary.logFC")]

# select significant
sig <- with(hmat, FDR < 0.05 & summary.logFC > 1)
hmat <- hmat[sig, ]

# order by FDR and select top n for plot
hmat <- hmat[order(hmat$FDR), ]

# format gene names and FDRs in row names
nms <- with(hmat, paste0(gene_name, " (", format(signif(FDR, 2)), ")"))
rownames(hmat) <- nms

hmat <- as.matrix(hmat[, c("self.average", "other.average")])
colnames(hmat) <- c("5-HT", "other")

# remove mitochondrial genes from heatmap
ix_mito <- grepl("^MT-", rownames(hmat))
table(ix_mito)
hmat <- hmat[!ix_mito, ]
dim(hmat)

# select top n
hmat <- hmat[1:100, ]

# create heatmap
hm <- Heatmap(
  hmat, 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  column_names_rot = 0, column_names_gp = gpar(fontsize = 10), column_names_centered = TRUE, 
  row_names_gp = gpar(fontsize = 9, fontface = "italic"), 
  column_title = "5-HT vs. other\nneuronal clusters", 
  column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  name = "mean\nlogcounts"
)

hm

# save heatmap
fn <- file.path(dir_plots, "DEtesting_heatmap_5HTvsOtherNeuronal")

pdf(paste0(fn, ".pdf"), width = 3.75, height = 12)
hm
dev.off()

png(paste0(fn, ".png"), width = 3.75 * 200, height = 12 * 200, res = 200)
hm
dev.off()


# -------------------------
# Spreadsheet: 5-HT neurons
# -------------------------

# save spreadsheet

cols <- c("gene_id", "gene_name", "sum_gene", "p.value", "FDR", "summary.logFC")
#cols <- c("gene_id", "gene_name", "sum_gene", "self.average", "other.average", "p.value", "FDR", "summary.logFC")
df <- marker_info_5HT[["16"]][, cols]
colnames(df)[c(4, 6)] <- c("p_value", "logFC")

# select significant
sig <- with(df, FDR < 0.05 & logFC > 1)
table(sig)
df <- df[sig, ]

# order by FDR
df <- df[order(df$FDR), ]

df <- as.data.frame(df)
rownames(df) <- NULL


# save .csv file
fn <- file.path(dir_outputs, "DEtesting_5HTvsOtherNeuronal.csv")
write.csv(df, file = fn, row.names = FALSE)

