###############################################################
# LC snRNA-seq analyses: clustering and supervised thresholding
# Lukas Weber, Oct 2022
###############################################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(bluster)
library(ggplot2)
library(ggVennDiagram)


dir_plots <- here("plots", "snRNAseq", "05a_clustering_main")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_filtered")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)


# ----------
# Clustering
# ----------

# clustering algorithm and parameters from OSCA
# two-stage clustering algorithm using high-resolution k-means and graph-based clustering

set.seed(121)
clus <- clusterCells(
  sce, 
  use.dimred = "PCA", 
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 2000), 
    second = NNGraphParam(k = 10)
  )
)
colLabels(sce) <- clus


# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# expression of key markers for NE, 5-HT, and cholinergic neuron populations
ix <- c(
  TH = which(rowData(sce)$gene_name == "TH"), 
  SLC6A2 = which(rowData(sce)$gene_name == "SLC6A2"), 
  DBH = which(rowData(sce)$gene_name == "DBH"), 
  SLC6A4 = which(rowData(sce)$gene_name == "SLC6A4"), 
  TPH2 = which(rowData(sce)$gene_name == "TPH2"), 
  SLC5A7 = which(rowData(sce)$gene_name == "SLC5A7")#, 
  #CHAT = which(rowData(sce)$gene_name == "CHAT"), 
  #ACHE = which(rowData(sce)$gene_name == "ACHE"), 
  #BCHE = which(rowData(sce)$gene_name == "BCHE"), 
  #SLC18A3 = which(rowData(sce)$gene_name == "SLC18A3"), 
  #PRIMA1 = which(rowData(sce)$gene_name == "PRIMA1")
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


# -----------------------
# Supervised thresholding
# -----------------------

# identify NE neuron nuclei based on positive expression of TH, SLC6A2, DBH

ix_TH <- which(rowData(sce)$gene_name == "TH")
ix_SLC6A2 <- which(rowData(sce)$gene_name == "SLC6A2")
ix_DBH <- which(rowData(sce)$gene_name == "DBH")

# number of nuclei
rbind(
  THpos = table(counts(sce)[ix_TH, ] > 0), 
  SLC6A2pos = table(counts(sce)[ix_SLC6A2, ] > 0), 
  DBHpos = table(counts(sce)[ix_DBH, ] > 0), 
  all3pos = table(counts(sce)[ix_TH, ] > 0 & counts(sce)[ix_SLC6A2, ] > 0 & counts(sce)[ix_DBH, ] > 0)
)

# identify nuclei
colData(sce)$posTH <- counts(sce)[ix_TH, ] > 0
colData(sce)$posSLC6A2 <- counts(sce)[ix_SLC6A2, ] > 0
colData(sce)$posDBH <- counts(sce)[ix_DBH, ] > 0

# identify nuclei: more strict thresholds
colData(sce)$posStrictTH <- counts(sce)[ix_TH, ] > 1
colData(sce)$posStrictSLC6A2 <- counts(sce)[ix_SLC6A2, ] > 1
colData(sce)$posStrictDBH <- counts(sce)[ix_DBH, ] > 1


# intersection set
colData(sce)$supervisedNE3of3 = colData(sce)$posTH & colData(sce)$posSLC6A2 & colData(sce)$posDBH
table(colData(sce)$supervisedNE3of3)
table(colData(sce)$Sample, colData(sce)$supervisedNE3of3)
# union set (> 0)
colData(sce)$supervisedNEunion3of3 = colData(sce)$posTH | colData(sce)$posSLC6A2 | colData(sce)$posDBH
table(colData(sce)$supervisedNEunion3of3)
table(colData(sce)$Sample, colData(sce)$supervisedNEunion3of3)
# union set (> 1)
colData(sce)$supervisedNEunionStrict3of3 = colData(sce)$posStrictTH | colData(sce)$posStrictSLC6A2 | colData(sce)$posStrictDBH
table(colData(sce)$supervisedNEunionStrict3of3)
table(colData(sce)$Sample, colData(sce)$supervisedNEunionStrict3of3)

# "2 out of 3" set
colData(sce)$supervised_NE = colData(sce)$posTH & colData(sce)$posDBH
table(colData(sce)$supervised_NE)
table(colData(sce)$Sample, colData(sce)$supervised_NE)


# check expression of key markers
res <- rowMeans(logcounts(sce)[ix, colData(sce)$supervised_NE])
names(res) <- names(ix)
res

res <- rowMeans(logcounts(sce)[ix, colData(sce)$supervisedNE3of3])
names(res) <- names(ix)
res

res <- rowMeans(logcounts(sce)[ix, colData(sce)$supervisedNEunion3of3])
names(res) <- names(ix)
res

res <- rowMeans(logcounts(sce)[ix, colData(sce)$supervisedNEunionStrict3of3])
names(res) <- names(ix)
res


# identify cholinergic interneurons based on positive expression of marker genes

ix_SLC5A7 = which(rowData(sce)$gene_name == "SLC5A7")
ix_CHAT = which(rowData(sce)$gene_name == "CHAT")
ix_ACHE = which(rowData(sce)$gene_name == "ACHE")
ix_BCHE = which(rowData(sce)$gene_name == "BCHE")
ix_SLC18A3 = which(rowData(sce)$gene_name == "SLC18A3")
ix_PRIMA1 = which(rowData(sce)$gene_name == "PRIMA1")

table(
  counts(sce)[ix_SLC5A7, ] > 0 & 
  counts(sce)[ix_CHAT, ] > 0 & 
  counts(sce)[ix_ACHE, ] > 0
)

table(
  counts(sce)[ix_SLC5A7, ] > 0 & 
  counts(sce)[ix_CHAT, ] > 0 & 
  counts(sce)[ix_ACHE, ] > 0 & 
  counts(sce)[ix_BCHE, ] > 0 & 
  counts(sce)[ix_SLC18A3, ] > 0 & 
  counts(sce)[ix_PRIMA1, ] > 0
)


# -----------------
# summarize results
# -----------------

# number of nuclei per cluster and sample

# unsupervised clustering
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)

# NE neuron cluster and 5-HT neuron cluster identified from marker genes above
clus_NE <- 6
clus_5HT <- 21

sum(colLabels(sce) == clus_NE)
sum(colLabels(sce) == clus_5HT)

tbl <- rbind(
  NE = table(colLabels(sce) == clus_NE, colData(sce)$Sample)[2, ], 
  `5HT` = table(colLabels(sce) == clus_5HT, colData(sce)$Sample)[2, ]
)
tbl
rowSums(tbl)


# supervised thresholding
table(colData(sce)$supervised_NE)
table(colData(sce)$supervised_NE, colData(sce)$Sample)[2, ]


# comparison between unsupervised clustering and supervised thresholding
table(
  unsupervised = colLabels(sce) == clus_NE, 
  supervised = colData(sce)$supervised_NE
)


# mitochondrial percentages in NE neuron clusters

# unsupervised
summary(colData(sce)$subsets_Mito_percent[colLabels(sce) == clus_NE])
# supervised
summary(colData(sce)$subsets_Mito_percent[colData(sce)$supervised_NE])


# for plotting

sce_plot <- sce
rownames(sce_plot) <- rowData(sce_plot)$gene_name

# unsupervised
sce_clusNE <- sce_plot[, colLabels(sce_plot) == clus_NE]
sce_clus5HT <- sce_plot[, colLabels(sce_plot) == clus_5HT]

# supervised
sce_supNE <- sce_plot[, colData(sce_plot)$supervised_NE]


genes_NE <- c("TH", "SLC6A2", "DBH")
genes_5HT <- c("TPH2", "SLC6A4")


# --------------------------------
# Venn diagrams comparing overlaps
# --------------------------------

colData(sce)$Key <- paste(colData(sce)$Sample, colData(sce)$Barcode, sep = "_")


# NE neurons: unsupervised clustering vs. supervised thresholding

x <- list(
  `clustering NE` = colData(sce)$Key[colLabels(sce) == clus_NE], 
  `supervised NE` = colData(sce)$Key[colData(sce)$supervised_NE]
)

ggVennDiagram(x) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = c("black", "black")) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

fn <- file.path(dir_plots, "overlap_NEneurons_clusteringVsSupervised")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 3.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 3.5)


# NE neurons: supervised thresholding on each of TH, SLC6A2, DBH

x <- list(
  `TH+` = colData(sce)$Key[colData(sce)$posTH], 
  `SLC6A2+` = colData(sce)$Key[colData(sce)$posSLC6A2], 
  `DBH+` = colData(sce)$Key[colData(sce)$posDBH]
)

ggVennDiagram(x) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = c("black", "black", "black")) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

fn <- file.path(dir_plots, "overlap_NEneurons_byMarker")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)


# NE neurons: supervised thresholding on each of TH, SLC6A2, DBH: more strict thresholds

x <- list(
  `TH++` = colData(sce)$Key[colData(sce)$posStrictTH], 
  `SLC6A2++` = colData(sce)$Key[colData(sce)$posStrictSLC6A2], 
  `DBH++` = colData(sce)$Key[colData(sce)$posStrictDBH]
)

ggVennDiagram(x) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = c("black", "black", "black")) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

fn <- file.path(dir_plots, "overlap_NEneurons_byMarker_strict")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)


# -------------------------------
# plot expression of marker genes
# -------------------------------

# unsupervised clustering

# plot expression of NE neuron marker genes
p <- gridExtra::grid.arrange(
  plotExpression(sce_clusNE, genes_NE, colour_by = "sum") + ggtitle("NE neuron cluster"), 
  plotExpression(sce_clusNE, genes_NE, colour_by = "detected") + ggtitle("NE neuron cluster"), 
  plotExpression(sce_clusNE, genes_NE, colour_by = "subsets_Mito_percent") + ggtitle("NE neuron cluster"), 
  ncol = 3
)
p
fn <- file.path(dir_plots, "clusteringNEneurons_expression")
ggsave(paste0(fn, ".pdf"), plot = p, width = 10, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 10, height = 4)


# plot expression of 5-HT neuron marker genes
p <- gridExtra::grid.arrange(
  plotExpression(sce_clus5HT, genes_5HT, colour_by = "sum") + ggtitle("5-HT neuron cluster"), 
  plotExpression(sce_clus5HT, genes_5HT, colour_by = "detected") + ggtitle("5-HT neuron cluster"), 
  plotExpression(sce_clus5HT, genes_5HT, colour_by = "subsets_Mito_percent") + ggtitle("5-HT neuron cluster"), 
  ncol = 3
)
p
fn <- file.path(dir_plots, "clustering5HTneurons_expression")
ggsave(paste0(fn, ".pdf"), plot = p, width = 10, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 10, height = 4)


# supervised thresholding

# plot expression of NE neuron marker genes
p <- gridExtra::grid.arrange(
  plotExpression(sce_supNE, genes_NE, colour_by = "sum") + ggtitle("NE neurons (supervised)"), 
  plotExpression(sce_supNE, genes_NE, colour_by = "detected") + ggtitle("NE neurons (supervised)"), 
  plotExpression(sce_supNE, genes_NE, colour_by = "subsets_Mito_percent") + ggtitle("NE neurons (supervised)"), 
  ncol = 3
)
p
fn <- file.path(dir_plots, "supervisedNEneurons_expression")
ggsave(paste0(fn, ".pdf"), plot = p, width = 10, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 10, height = 4)


# -------------------------
# plot UMAP representations
# -------------------------

# identify populations of interest
colData(sce)$unsupervised_NE <- colLabels(sce) == clus_NE
colData(sce)$unsupervised_5HT <- colLabels(sce) == clus_5HT


# unsupervised clustering

# NE neurons
plotReducedDim(sce, dimred = "UMAP", colour_by = "unsupervised_NE") + 
  scale_color_manual(values = c("navy", "red"), name = "NE neurons") + 
  ggtitle("Unsupervised clustering")

fn <- file.path(dir_plots, "UMAP_clusteringNEneurons")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)


# 5-HT neurons
plotReducedDim(sce, dimred = "UMAP", colour_by = "unsupervised_5HT") + 
  scale_color_manual(values = c("navy", "red"), name = "5-HT neurons") + 
  ggtitle("Unsupervised clustering")

fn <- file.path(dir_plots, "UMAP_clustering5HTneurons")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)


# all clusters
pal <- unname(palette.colors(36, "Polychrome 36"))
plotReducedDim(sce, dimred = "UMAP", colour_by = "label") + 
  scale_color_manual(values = pal, name = "cluster") + 
  ggtitle("Unsupervised clustering")

fn <- file.path(dir_plots, "UMAP_clustering")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4.75)
ggsave(paste0(fn, ".png"), width = 6, height = 4.75)


# supervised thresholding

# NE neurons
plotReducedDim(sce, dimred = "UMAP", colour_by = "supervised_NE") + 
  scale_color_manual(values = c("navy", "red"), name = "NE neurons") + 
  ggtitle("Supervised thresholding")

fn <- file.path(dir_plots, "UMAP_supervisedNEneurons")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)


# sample IDs
# note: batch integration was not included due to our interest in rare populations (LC-NE neurons)

# sample IDs
plotReducedDim(sce, dimred = "UMAP", colour_by = "Sample") + 
  ggtitle("Sample IDs")

fn <- file.path(dir_plots, "UMAP_sampleIDs")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)


# -----------
# Save object
# -----------

fn_out <- here("processed_data", "SCE", "sce_clustering_main")
saveRDS(sce, paste0(fn_out, ".rds"))

