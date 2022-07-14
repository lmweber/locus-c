##################################################################################################
# LC snRNA-seq analysis alternative workflows: unsupervised clustering and supervised thresholding
# Lukas Weber, July 2022
##################################################################################################

# run in interactive session on JHPCE

# screen -S LC
# qrsh -l mem_free=20G,h_vmem=22G,h_fsize=200G
# cd /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn_alt
# module load conda_R/devel
# R


library(SingleCellExperiment)
library(here)
library(DropletUtils)
library(scater)
library(scran)
library(bluster)
library(ggplot2)


dir_plots <- here("plots", "snRNAseq_alt")


# -----------------
# Create SCE object
# -----------------

# sample metadata
sample_info <- read.table(here("fastq", "snRNA-seq", "sample_libs_info.tsv"))


# note: nuclei calling using Cell Ranger (instead of emptyDrops)
# i.e. using "filtered" instead of "raw" files

sample_info$path <- file.path(
  here("processed_data", "cellranger"), 
  sample_info$V1, 
  "outs", 
  "filtered_feature_bc_matrix"
)
stopifnot(all(file.exists(sample_info$path)))


# build SCE
sce <- read10xCounts(
  samples = sample_info$path, 
  sample.names = paste0(sample_info$V3, "_", sample_info$V2), 
  type = "sparse", 
  col.names = TRUE
)


## Read in the gene information from the annotation GTF file
## code from: https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/build_basic_sce.R
gtf <-
  rtracklayer::import(
    "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
  )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[ , c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SCE object
rowRanges(sce) <- gtf[match_genes]

## Inspect object
sce

# Number of nuclei per sample
table(sce$Sample)


# save object
# fn_out <- here("processed_data", "SCE_alt", "sce_raw")
# saveRDS(sce, paste0(fn_out, ".rds"))


# load object
# fn <- here("processed_data", "SCE_alt", "sce_raw")
# sce <- readRDS(paste0(fn, ".rds"))


# --------------------
# Quality control (QC)
# --------------------

# perform QC on sum UMI counts and number of detected genes
# note: not using mitochondrial percentage, since this is high for biological reasons in LC-NE neurons

# store QC metrics
sce <- addPerCellQC(sce, subsets = list(Mito = which(seqnames(sce) == "chrM")))

# check distributions
quantile(colData(sce)$sum, seq(0, 1, by = 0.1))
quantile(colData(sce)$detected, seq(0, 1, by = 0.1))

# check 3 median absolute deviations (MADs)
reasons <- perCellQCFilters(sce)
attr(reasons$low_lib_size, "thresholds")
attr(reasons$low_n_features, "thresholds")

# note: 3 MADs is below minimum values observed for both sum UMIs and detected genes, 
# and we do not want to filter out high values since LC-NE neurons are very large; 
# so we keep all cells

colData(sce)$discard <- FALSE


# note high mitochondrial percentages
quantile(colData(sce)$subsets_Mito_percent, seq(0, 1, by = 0.1))
quantile(colData(sce)$subsets_Mito_percent, seq(0.9, 1, by = 0.01))
mean(colData(sce)$subsets_Mito_percent > 10)
mean(colData(sce)$subsets_Mito_percent > 20)


# plot QC metrics
p <- gridExtra::grid.arrange(
  plotColData(sce, x = "Sample", y = "sum", colour_by = "subsets_Mito_percent") + 
    scale_y_log10() + guides(color = guide_legend(title = "mito")) + ggtitle("Total count"), 
  plotColData(sce, x = "Sample", y = "detected", colour_by = "subsets_Mito_percent") + 
    scale_y_log10() + guides(color = guide_legend(title = "mito")) + ggtitle("Detected genes"), 
  plotColData(sce, x = "Sample", y = "subsets_Mito_percent") + 
    ggtitle("Mito percent"), 
  ncol=3
)

p

fn <- file.path(dir_plots, "QC_metrics")
ggsave(paste0(fn, ".pdf"), plot = p, width = 10, height = 3)
ggsave(paste0(fn, ".png"), plot = p, width = 10, height = 3)


# investigate high-mitochondrial nuclei
nonzero_TH <- counts(sce)[which(rowData(sce)$gene_name == "TH"), ] > 0
summary(colData(sce)$subsets_Mito_percent[nonzero_TH])


# --------------------------
# Filter low-expressed genes
# --------------------------

# note: keep both protein-coding and non-protein-coding genes

n_umis <- 30
ix_remove <- rowSums(counts(sce)) < n_umis
table(ix_remove)

dim(sce)

sce <- sce[!ix_remove, ]

dim(sce)


# ---------------------------
# Normalization and logcounts
# ---------------------------

# normalization by deconvolution

set.seed(123)
quickclus <- quickCluster(sce)
table(quickclus)

sce <- computeSumFactors(sce, cluster = quickclus)
sce <- logNormCounts(sce)

summary(sizeFactors(sce))


# -----------------
# Feature selection
# -----------------

# note: keep mitochondrial genes, since these are biologically meaningful in this dataset
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$gene_name)
table(is_mito)
rowData(sce)$gene_name[is_mito]

# fit mean-variance relationship
dec <- modelGeneVar(sce)
# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)


# ------------------------
# Dimensionality reduction
# ------------------------

# note: not applying any batch integration due to our interest in rare populations

# calculate PCA (on top HVGs)
set.seed(123)
sce <- runPCA(sce, subset_row = top_hvgs)

# calculate UMAP (on top PCs)
set.seed(123)
sce <- runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) <- paste0("UMAP", seq_len(ncol(reducedDim(sce, "UMAP"))))

reducedDimNames(sce)


# ----------
# Clustering
# ----------

# clustering algorithm and parameters from OSCA
# two-stage clustering using high-resolution k-means and graph-based clustering

set.seed(123)
clus <- clusterCells(
  sce, 
  use.dimred = "PCA", 
  BLUSPARAM = TwoStepParam(
    first = KmeansParam(centers = 2000), 
    second = NNGraphParam(k = 5)
  )
)
colLabels(sce) <- clus


# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# check expression of key markers
ix <- c(
  TH = which(rowData(sce)$gene_name == "TH"), 
  SLC6A2 = which(rowData(sce)$gene_name == "SLC6A2"), 
  DBH = which(rowData(sce)$gene_name == "DBH"), 
  TPH2 = which(rowData(sce)$gene_name == "TPH2"), 
  SLC6A4 = which(rowData(sce)$gene_name == "SLC6A4"), 
  SLC5A7 = which(rowData(sce)$gene_name == "SLC5A7")
)

n_clus <- length(table(colLabels(sce)))

res_list <- list()
for (k in seq_len(n_clus)) {
  res_list[[k]] <- rowMeans(logcounts(sce)[ix, colLabels(sce) == k])
}
res_mat <- do.call("rbind", res_list)
rownames(res_mat) <- seq_len(n_clus)
colnames(res_mat) <- names(ix)

res_mat


# -----------------------
# Supervised thresholding
# -----------------------

# identify NE neurons based on nonzero expression of TH, SLC6A2, DBH

ix_TH <- which(rowData(sce)$gene_name == "TH")
ix_SLC6A2 <- which(rowData(sce)$gene_name == "SLC6A2")
ix_DBH <- which(rowData(sce)$gene_name == "DBH")

# number of nuclei with positive expression of 1 of or all 3 of TH, SLC6A2, DBH
rbind(
  THpos = table(counts(sce)[ix_TH, ] > 0), 
  SLC6A2pos = table(counts(sce)[ix_SLC6A2, ] > 0), 
  DBHpos = table(counts(sce)[ix_DBH, ] > 0), 
  all3pos = table(counts(sce)[ix_TH, ] > 0 & counts(sce)[ix_SLC6A2, ] > 0 & counts(sce)[ix_DBH, ] > 0)
)

# identify NE neurons: intersection set
ix_NEneurons <- counts(sce)[ix_TH, ] > 0 & counts(sce)[ix_SLC6A2, ] > 0 & counts(sce)[ix_DBH, ] > 0
table(ix_NEneurons)

# identify NE neurons: union set
ix_union <- counts(sce)[ix_TH, ] > 0 | counts(sce)[ix_SLC6A2, ] > 0 | counts(sce)[ix_DBH, ] > 0
table(ix_union)


stopifnot(length(ix_NEneurons) == ncol(sce))

colData(sce)$supervisedNE <- as.factor(as.numeric(ix_NEneurons))

head(colData(sce), 2)


# check expression of key markers
res <- rowMeans(logcounts(sce)[ix, colData(sce)$supervisedNE == 1])
names(res) <- names(ix)
res


# -----------------
# summarize results
# -----------------

# number of nuclei per cluster and sample

# unsupervised clustering
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)

# NE neuron cluster and 5-HT neuron cluster identified from marker genes above
clus_NE <- 11
clus_5HT <- 41

sum(colLabels(sce) == clus_NE)
sum(colLabels(sce) == clus_5HT)

rbind(
  NE = table(colLabels(sce) == clus_NE, colData(sce)$Sample)[2, ], 
  `5HT` = table(colLabels(sce) == clus_5HT, colData(sce)$Sample)[2, ]
)


# supervised thresholding
table(colData(sce)$supervisedNE)
table(colData(sce)$supervisedNE, colData(sce)$Sample)[2, ]


# comparison between unsupervised clustering and supervised thresholding
table(
  unsupervised = colLabels(sce) == clus_NE, 
  supervised = colData(sce)$supervisedNE == 1
)


# for plotting
sce_plot <- sce
rownames(sce_plot) <- rowData(sce_plot)$gene_name

# unsupervised clustering
sce_clusNE <- sce_plot[, colLabels(sce_plot) == clus_NE]
sce_clus5HT <- sce_plot[, colLabels(sce_plot) == clus_5HT]

# supervised thresholding
sce_supNE <- sce_plot[, colData(sce_plot)$supervisedNE == 1]


genes_NE <- c("TH", "SLC6A2", "DBH")
genes_5HT <- c("TPH2", "SLC6A4")


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
colData(sce)$unsupervisedNE <- colLabels(sce) == clus_NE
colData(sce)$unsupervised5HT <- colLabels(sce) == clus_5HT


# unsupervised clustering

# NE neurons
plotReducedDim(sce, dimred = "UMAP", colour_by = "unsupervisedNE") + 
  scale_color_manual(values = c("navy", "red"), name = "NE neurons") + 
  ggtitle("Unsupervised clustering")

fn <- file.path(dir_plots, "clusteringNEneurons_UMAP")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)


# 5-HT neurons
plotReducedDim(sce, dimred = "UMAP", colour_by = "unsupervised5HT") + 
  scale_color_manual(values = c("navy", "red"), name = "5-HT neurons") + 
  ggtitle("Unsupervised clustering")

fn <- file.path(dir_plots, "clustering5HTneurons_UMAP")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)


# supervised thresholding

# NE neurons
plotReducedDim(sce, dimred = "UMAP", colour_by = "supervisedNE") + 
  scale_color_manual(values = c("navy", "red"), name = "NE neurons") + 
  ggtitle("Supervised thresholding")

fn <- file.path(dir_plots, "supervisedNEneurons_UMAP")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5)

