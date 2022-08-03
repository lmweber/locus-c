######################################################
# LC analyses: script to save SPE object for Shiny app
# Lukas Weber, Jun 2022
######################################################

# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_preprocessing.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)


# ------------------------------------
# additional information for Shiny app
# ------------------------------------

# using code from Leonardo Collado-Torres available at:
# http://research.libd.org/spatialLIBD/articles/TenX_data_download.html

## Spot info

# using expected column names
colData(spe)$key <- colData(spe)$key_id

colData(spe)$sum_umi <- colSums(counts(spe))
colData(spe)$sum_gene <- colSums(counts(spe) > 0)


## Gene info
## note: already added in preprocessing script


## Additional info for spatialLIBD

## Add information used by spatialLIBD
rowData(spe)$gene_search <- paste(
  rowData(spe)$gene_name, rowData(spe)$gene_id, sep = "; "
)
## Compute chrM expression and chrM expression ratio
is_mito <- which(seqnames(spe) == "chrM")
colData(spe)$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
colData(spe)$expr_chrM_ratio <- colData(spe)$expr_chrM / colData(spe)$sum_umi


## Spots over tissue
## Keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue]
dim(spe)


## Filter zeros (genes and spots)
## Note: do this after subsetting spots over tissue to avoid creating new zeros

## Remove genes with zero counts
no_expr <- which(rowSums(counts(spe)) == 0)
## Number of genes with no counts
length(no_expr)
## Proportion of genes with no counts
length(no_expr) / nrow(spe)
## Remove from object
spe <- spe[-no_expr, , drop = FALSE]

## Remove spots with zero counts
spots_no_counts <- which(colData(spe)$sum_umi == 0)
## Number of spots with no counts
length(spots_no_counts)
## Proportion of spots with no counts
length(spots_no_counts) / ncol(spe)
## Remove from object
spe <- spe[, -spots_no_counts, drop = FALSE]

dim(spe)


## Default cluster labels
## Add a column of default cluster labels
colData(spe)$all <- "all"


## Manual annotations
## Add a variable for saving manual annotations
colData(spe)$ManualAnnotation <- "NA"


## Genes of interest
## UMIs per spot for some key genes of interest
colData(spe)$counts_TH <- counts(spe)[which(rowData(spe)$gene_name == "TH"), ]
colData(spe)$counts_SLC6A2 <- counts(spe)[which(rowData(spe)$gene_name == "SLC6A2"), ]


# --------------------
# quality control (QC)
# --------------------

# note: keep all spots for Shiny app


# ---------------------------
# normalization and logcounts
# ---------------------------

# using scater and scran packages

# quick clustering for pool-based size factors
# with blocks by sample
set.seed(123)
qclus <- quickCluster(spe, block = colData(spe)$sample_id)

table(qclus)
table(colData(spe)$sample_id, qclus)

# calculate size factors
spe <- computeSumFactors(spe, cluster = qclus)
summary(sizeFactors(spe))

# note: remove small number of spots with size factors == 0
table(sizeFactors(spe) == 0)
sum(is.na(sizeFactors(spe)))
dim(spe)

spe <- spe[, !(sizeFactors(spe) == 0)]
dim(spe)

# calculate logcounts
spe <- logNormCounts(spe)
assayNames(spe)


## Genes of interest
## logcounts per spot for some key genes of interest
colData(spe)$logcounts_TH <- logcounts(spe)[which(rowData(spe)$gene_name == "TH"), ]
colData(spe)$logcounts_SLC6A2 <- logcounts(spe)[which(rowData(spe)$gene_name == "SLC6A2"), ]


# -----------
# save object
# -----------

# save as .rds and .RData
fn_out <- here("processed_data", "SPE", "LC_Shiny")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

