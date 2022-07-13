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


# note: using all droplets
# i.e. using "raw" instead of "filtered" files
sample.info$path <- file.path(
  here("processed_data", "cellranger"),
  sample.info$V1,
  "outs",
  "raw_feature_bc_matrix"
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


# -----------
# Ambient RNA
# -----------

# investigate ambient RNA in empty droplets

# assume droplets with n_umis <= 100 are definitely empty and contain only ambient RNA
# (see Lun et al. 2019, section "Testing for deviations from the ambient profile")

colData(sce.lc)$sum <- colSums(counts(sce.lc))

dim(sce.lc)
anyNA(colData(sce.lc)$sum)

# remove droplets with zero expression
ix_zeros <- colData(sce.lc)$sum == 0
table(ix_zeros)
sce.lc <- sce.lc[, !ix_zeros]
dim(sce.lc)

# identify droplets with very low UMIs
colData(sce.lc)$is_low <- colData(sce.lc)$sum <= 100
table(colData(sce.lc)$is_low)

# distribution of these low counts
quantile(colData(sce.lc)$sum[colData(sce.lc)$is_low], seq(0, 1, by = 0.1))


# check genes of interest

ix_genes <- c(
  TH = which(rowData(sce.lc)$gene_name == "TH"), 
  SLC6A2 = which(rowData(sce.lc)$gene_name == "SLC6A2"), 
  TPH2 = which(rowData(sce.lc)$gene_name == "TPH2"), 
  SLC6A4 = which(rowData(sce.lc)$gene_name == "SLC6A4"), 
  SLC5A7 = which(rowData(sce.lc)$gene_name == "SLC5A7"), 
  SNAP25 = which(rowData(sce.lc)$gene_name == "SNAP25"), 
  SYT1 = which(rowData(sce.lc)$gene_name == "SYT1")
)

counts(sce.lc)[ix_genes, colData(sce.lc)$sum == 100]

