#####################################################################
# LC snRNA-seq analysis alternative workflow: unsupervised clustering
# Lukas Weber, July 2022
#####################################################################

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


# ------------------------
# Dimensionality reduction
# ------------------------

# note: not applying any batch integration due to our interest in rare populations

# PCA
set.seed(123)
sce <- runPCA(sce)

# UMAP
set.seed(123)
sce <- runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) <- paste0("UMAP", seq_len(ncol(reducedDim(sce, "UMAP"))))

reducedDimNames(sce)


# ----------
# Clustering
# ----------

# clustering algorithm and parameters from OSCA
# two-stage clustering using high-resolution k-means and graph-based clustering

set.seed(100)
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


# --------------------------
# summarize and plot results
# --------------------------

# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)

# clusters identified based on marker genes above
# cluster 21: NE neurons
# cluster 7: 5-HT neurons


sum(colLabels(sce) == 21)
sum(colLabels(sce) == 7)

rbind(
  NE = table(colLabels(sce) == 21, colData(sce)$Sample)[2, ], 
  `5HT` = table(colLabels(sce) == 7, colData(sce)$Sample)[2, ]
)


# for plotting
sce_plot <- sce
rownames(sce_plot) <- rowData(sce_plot)$gene_name

sce_NE <- sce_plot[, colLabels(sce_plot) == 21]
sce_5HT <- sce_plot[, colLabels(sce_plot) == 7]

genes_NE <- c("TH", "SLC6A2")
genes_5HT <- c("TPH2", "SLC6A4")


# plot expression of NE neuron marker genes
p <- gridExtra::grid.arrange(
  plotExpression(sce_NE, genes_NE, colour_by = "sum") + ggtitle("NE neuron cluster"), 
  plotExpression(sce_NE, genes_NE, colour_by = "detected") + ggtitle("NE neuron cluster"), 
  plotExpression(sce_NE, genes_NE, colour_by = "subsets_Mito_percent") + ggtitle("NE neuron cluster"), 
  ncol = 3
)

p

fn <- file.path(dir_plots, "clusterNEneurons_metrics")
ggsave(paste0(fn, ".pdf"), plot = p, width = 10, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 10, height = 4)


# plot expression of 5-H% neuron marker genes
p <- gridExtra::grid.arrange(
  plotExpression(sce_5HT, genes_5HT, colour_by = "sum") + ggtitle("5-HT neuron cluster"), 
  plotExpression(sce_5HT, genes_5HT, colour_by = "detected") + ggtitle("5-HT neuron cluster"), 
  plotExpression(sce_5HT, genes_5HT, colour_by = "subsets_Mito_percent") + ggtitle("5-HT neuron cluster"), 
  ncol = 3
)

p

fn <- file.path(dir_plots, "cluster5HTneurons_metrics")
ggsave(paste0(fn, ".pdf"), plot = p, width = 10, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 10, height = 4)


