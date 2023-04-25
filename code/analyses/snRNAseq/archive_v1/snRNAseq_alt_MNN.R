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

set.seed(1000)
quickclus <- quickCluster(sce.lc)
sce.lc <- computeSumFactors(sce.lc, cluster = quickclus)
sce.lc <- logNormCounts(sce.lc)

summary(sizeFactors(sce.lc))


# -----------------
# Variance modeling
# -----------------

# code from OSCA
# note: keeping large number of HVGs for batch integration later

set.seed(1000)
dec.lc <- modelGeneVarByPoisson(sce.lc, block = sce.lc$Sample)
top.lc <- getTopHVGs(dec.lc, n = 5000)


# -----------------
# Batch integration
# -----------------

# using MNN

set.seed(1000)
merged.lc <- fastMNN(sce.lc, batch = sce.lc$Sample, subset.row = top.lc)

reducedDim(sce.lc, "MNN") <- reducedDim(merged.lc, "corrected")

# percentage of variance lost as diagnostic measure
metadata(merged.lc)$merge.info$lost.var


# ------------------------
# Dimensionality reduction
# ------------------------

# note: memory error

#set.seed(1000)
#sce.lc <- runUMAP(sce.lc, dimred = "MNN")


# ----------
# Clustering
# ----------

# graph-based clustering

# g <- buildSNNGraph(sce.lc, k=10, use.dimred = "MNN")
# clust <- igraph::cluster_walktrap(g)$membership
# colLabels(sce.lc) <- factor(clust)
# 
# table(colLabels(sce.lc))
# 
# plotUMAP(sce.lc, colour_by ="label")


# graph-based clustering generates too many clusters
# use two-stage k-means / graph-based clustering from OSCA instead
# note: for more resolution, can increase centers and/or use lower k during graph construction

set.seed(1000)
colLabels(sce.lc) <- clusterRows(reducedDim(sce.lc, "MNN"), 
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

rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 17])


# cluster 17 is 5-HT neurons (210 nuclei)
# ENSG00000180176 ENSG00000103546 ENSG00000139287 ENSG00000108576 ENSG00000115665 
# 0.03907259      0.06290066      3.97980795      2.87732649      0.09577232 


quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "TH", ], seq(0.9, 1, by = 0.01))
quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "SLC6A2", ], seq(0.9, 1, by = 0.01))
quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "TPH2", ], seq(0.9, 1, by = 0.01))
quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "SLC6A4", ], seq(0.9, 1, by = 0.01))
quantile(logcounts(sce.lc)[rowData(sce.lc)$gene_name == "SLC5A7", ], seq(0.9, 1, by = 0.01))


table(colData(sce.lc)[logcounts(sce.lc)[rowData(sce.lc)$gene_name == "TH", ] > 4, ]$label)

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
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 20])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 21])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 22])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 23])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 24])
rowMeans(logcounts(sce.lc)[ix, colLabels(sce.lc) == 25])


# ---------------------------------------
# Supervised identification of NE neurons
# ---------------------------------------

ix2 <- c(
  SNAP25 = which(rowData(sce.lc)$gene_name == "SNAP25"), 
  SYT1 = which(rowData(sce.lc)$gene_name == "SYT1")
)

ix_supervised_NE <- 
  logcounts(sce.lc)[rowData(sce.lc)$gene_name == "TH", ] > 1 & 
  logcounts(sce.lc)[rowData(sce.lc)$gene_name == "SLC6A2", ] > 1 & 
  colData(sce.lc)$sum > 10000
table(ix_supervised_NE)

rowMeans(logcounts(sce.lc)[ix, ix_supervised_NE])
rowMeans(logcounts(sce.lc)[ix2, ix_supervised_NE])

colData(sce.lc[ix, ix_supervised_NE])

colSums(table(colnames(sce.lc)[ix_supervised_NE], colData(sce.lc)$Sample[ix_supervised_NE]))


# plot expression
markers <- c("TH", "SLC6A2", "DBH", "SLC18A2", "GCH1", "DDC", 
             "TPH2", "SLC6A4", 
             "SLC5A7", "SLC18A3", "CHAT", "ACHE", 
             "SNAP25", "SYT1")

sce.lc_NE <- sce.lc
rownames(sce.lc_NE) <- rowData(sce.lc_NE)$gene_name
sce.lc_NE <- sce.lc_NE[markers, ix_supervised_NE]


png("expression_supervised_NE.png", width = 800, height = 400)
plotExpression(sce.lc_NE, markers)
dev.off()


# compare to nuclei in NE neuron cluster from hierarchical merged clustering

nuclei_NE <- sort(rownames(table(colnames(sce.lc)[ix_supervised_NE], colData(sce.lc)$Sample[ix_supervised_NE])))

nuclei_ref <- c("2_ACGATCACAGCCTATA-1", "2_ACGGAAGCAACCGTGC-1", "2_AGGACTTAGTGCAGGT-1", 
                "2_CACGAATCATCCAACA-1", "2_CATCCCACAACAGCTT-1", "2_CATGCCTGTAGGCAAC-1", 
                "2_CCCTCAAGTATTTCCT-1", "2_CGTTCTGTCTAAGAAG-1", "2_GAAGGACCATCATGAC-1", 
                "2_GACTATGGTAGAATAC-1", "2_GGACGTCCAAAGGGTC-1", "2_GGGCGTTTCTAGTACG-1", 
                "2_GGGCTCATCGCCGAGT-1", "2_GTAGAAAAGACTTGTC-1", "2_GTATTTCGTCTCAGGC-1", 
                "2_GTCAGCGAGTACTGGG-1", "2_TACTTACGTCGTTGCG-1", "2_TCTGTCGCATGCAGGA-1", 
                "2_TCTTTGATCATCTGTT-1", "2_TGTTGGACATTGAGCT-1", "2_TTAGGCATCGCCAATA-1", 
                "2_TTCTCTCCAAGACAAT-1", "2_TTGCTGCAGGGTTAGC-1", "2_TTGTTCAAGTCGTCTA-1", 
                "3_AATCGTGCAGCGTGCT-1", "3_ACCAAACTCTTCCAGC-1", "3_CAGATCAAGGCGATAC-1", 
                "3_CATTGCCGTCAGATTC-1", "3_CCGGGTAAGAAGCCTG-1", "3_CTCGAGGGTGCCTAAT-1", 
                "3_CTGCATCCAACCGCCA-1", "3_GGTAACTGTACCTAAC-1", "3_TATTGCTAGTGGATAT-1", 
                "3_TCTACATCATTCGATG-1", "3_TTGTGTTCACTAGTAC-1", "3_TTGTGTTTCACTCACC-1")

sum(nuclei_NE %in% nuclei_ref)

