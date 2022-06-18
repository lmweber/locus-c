##################################################################
# LC snRNA-seq analysis pipeline: testing some alternative methods
# adapting code from Matt N Tran, June 2022
# Lukas Weber, June 2022
##################################################################

# run in interactive session on JHPCE

# screen -S LC
# qrsh -l mem_free=20G,h_vmem=22G,h_fsize=200G
# cd /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn
# module load conda_R/4.1.x
# R


library(SingleCellExperiment)
library(here)
library(DropletUtils)
library(scater)
library(scran)
library(batchelor)
library(bluster)
library(harmony)


# -----------------
# Create SCE object
# -----------------

# Basic sample info
sample.info <- read.table(here("fastq", "snRNA-seq","sample_libs_info.tsv"))


# note: using Cell Ranger nuclei calling here
# i.e. modified to use "filtered" instead of "raw" files
# then require additional QC later on sum UMIs and detected genes

sample.info$path <- file.path(
  here("processed_data", "cellranger"),
  sample.info$V1,
  "outs",
  "filtered_feature_bc_matrix"
)
stopifnot(all(file.exists(sample.info$path)))


## Build basic SCE
sce <- read10xCounts(
  samples = sample.info$path,
  sample.names = paste0(sample.info$V3,"_",sample.info$V2),
  type = "sparse",
  col.names = TRUE
)


## Read in the gene information from the annotation GTF file
# (following Leo's method in https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/build_basic_sce.R)
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

# Re-name to sce.lc
sce.lc <- sce


# ---------------
# Quality control
# ---------------

# QC code from OSCA
# performing QC on all metrics since Cell Ranger nuclei calling was used

sce.lc <- addPerCellQC(sce.lc, subsets = list(Mito = which(seqnames(sce.lc) == "chrM")))

qc <- quickPerCellQC(colData(sce.lc), batch=colData(sce.lc)$Sample, 
                     sub.fields="subsets_Mito_percent")

# number of discarded nuclei per sample
colSums(as.matrix(qc))
table(qc$discard, colData(sce.lc)$Sample)

sce.lc <- sce.lc[, !qc$discard]

table(colData(sce.lc)$Sample)


# plotting code from OSCA
# gridExtra::grid.arrange(
#   plotColData(sce.lc, x="Sample", y="sum", colour_by="discard") +
#     scale_y_log10() + ggtitle("Total count"),
#   plotColData(sce.lc, x="Sample", y="detected", colour_by="discard") +
#     scale_y_log10() + ggtitle("Detected features"),
#   plotColData(sce.lc, x="Sample", y="subsets_Mito_percent",
#               colour_by="discard") + ggtitle("Mito percent"),
#   ncol=2
# )


# --------------------------
# Filter low-expressed genes
# --------------------------

# remove non-protein-coding genes

table(rowData(sce.lc)$gene_type)

sce.lc <- sce.lc[rowData(sce.lc)$gene_type == "protein_coding", ]

dim(sce.lc)
table(rowData(sce.lc)$gene_type)


# remove low-expressed genes

n_umis <- 10
ix_remove <- rowSums(counts(sce.lc)) < n_umis
table(ix_remove)

sce.lc <- sce.lc[!ix_remove, ]

dim(sce.lc)


# ---------------------------
# Normalization and logcounts
# ---------------------------

# code from OSCA

set.seed(123)
quickclus <- quickCluster(sce.lc)
sce.lc <- computeSumFactors(sce.lc, cluster = quickclus)
sce.lc <- logNormCounts(sce.lc)

summary(sizeFactors(sce.lc))


# -----------------
# Batch integration
# -----------------

# using Harmony

data_mat <- logcounts(sce.lc)
sample_ids <- colData(sce.lc)$Sample

set.seed(123)
harmony_dims <- HarmonyMatrix(
  data_mat = data_mat, 
  meta_data = sample_ids, 
  do_pca = TRUE
)

colnames(harmony_dims) <- paste0("HARM", seq_len(ncol(harmony_dims)))

dim(harmony_dims)
head(harmony_dims, 2)

# store in SCE object
reducedDim(sce.lc, "HARM") <- harmony_dims

reducedDimNames(sce.lc)


# ------------------------
# Dimensionality reduction
# ------------------------

# note: memory error

#set.seed(123)
#sce.lc <- runUMAP(sce.lc, dimred = "HARM")


# ----------
# Clustering
# ----------

# graph-based clustering

# g <- buildSNNGraph(sce.lc, k=10, use.dimred = "HARM")
# clust <- igraph::cluster_walktrap(g)$membership
# colLabels(sce.lc) <- factor(clust)
# 
# table(colLabels(sce.lc))
# 
# plotUMAP(sce.lc, colour_by ="label")


# graph-based clustering generates too many clusters
# use two-stage k-means / graph-based clustering from OSCA instead
# note: for more resolution, can increase centers and/or use lower k during graph construction

set.seed(123)
colLabels(sce.lc) <- clusterRows(reducedDim(sce.lc, "HARM"),
                                 TwoStepParam(KmeansParam(centers = 1000), NNGraphParam(k = 5)))

table(colLabels(sce.lc))


# number of nuclei per cluster per sample

table(colLabels(sce.lc))
table(colLabels(sce.lc), colData(sce.lc)$Sample)


# check expression of key markers

ix <- c(
  TH = which(rowData(sce.lc)$gene_name == "TH"), 
  SLC6A2 = which(rowData(sce.lc)$gene_name == "SLC6A2"), 
  TPH2 = which(rowData(sce.lc)$gene_name == "TPH2"), 
  SLC6A4 = which(rowData(sce.lc)$gene_name == "SLC6A4"), 
  SLC5A7 = which(rowData(sce.lc)$gene_name == "SLC5A7")
)


quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "TH", ], seq(0.9, 1, by = 0.01))
quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "SLC6A2", ], seq(0.9, 1, by = 0.01))
quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "TPH2", ], seq(0.9, 1, by = 0.01))
quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "SLC6A4", ], seq(0.9, 1, by = 0.01))
quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "SLC5A7", ], seq(0.9, 1, by = 0.01))


table(colData(sce.lc)[logcounts(sce.lc)[rowData(sce.lc)$gene_name == "TH", ] > 3, ]$label)

rowMeans(logcounts(sce.lc)[ix, colData(sce.lc)$sum > 50000])


rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 1])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 2])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 3])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 4])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 5])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 6])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 7])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 8])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 9])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 10])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 11])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 12])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 13])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 14])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 15])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 16])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 17])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 18])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 19])

