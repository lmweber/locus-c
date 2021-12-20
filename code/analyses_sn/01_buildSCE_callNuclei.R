### LC snRNA-seq analysis
### Build SCE from raw count matrices; perform nuclei calling
### qrsh -l bluejay,mf=32G,h_vmem=34G -pe local 2
### Adopted from https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/build_basic_sce.R
### Initiated MNT 13Dec2021

library(SingleCellExperiment)
library(DropletUtils)
library(rtracklayer)
library(BiocParallel)
library(jaffelab)
library(here)
library(sessioninfo)

here()
    # [1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"


### Create raw SCE ====
# Basic sample info
sample.info <- read.table(here("fastq", "snRNA-seq","sample_libs_info.tsv"))
    # First column is how the Cell Ranger output is organized

sample.info$path <- file.path(
  here("processed_data", "cellranger"),
  sample.info$V1,
  "outs",
  "raw_feature_bc_matrix"
)
stopifnot(all(file.exists(sample.info$path)))

## Build basic SCE (will add more subject-level colData later, once obtained)
Sys.time()
    # [1] "2021-12-15 15:53:02 EST"
sce <- read10xCounts(
  samples = sample.info$path,
  sample.names = paste0(sample.info$V3,"_",sample.info$V2),
  type = "sparse",
  col.names = TRUE
)
Sys.time()
    # [1] "2021-12-15 15:56:00 EST"

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
    # class: SingleCellExperiment 
    # dim: 36601 4857931 
    # metadata(1): Samples
    # assays(1): counts
    # rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
    #   ENSG00000277196
    # rowData names(6): source type ... gene_name gene_type
    # colnames(4857931): 1_AAACCCAAGAAACCCA-1 1_AAACCCAAGAAACTGT-1 ...
    #   3_TTTGTTGTCTTTGGAG-1 3_TTTGTTGTCTTTGGCT-1
    # colData names(2): Sample Barcode
    # reducedDimNames(0):
    # mainExpName: NULL
    # altExpNames(0):

table(sce$Sample)
    # Br2701_LC Br6522_LC Br8079_LC 
    #   1824762   1377516   1655653

## Size in Gb
lobstr::obj_size(sce) / 1024^3
    # 2.08 GB

# Re-name to sce.lc
sce.lc <- sce

## Save to /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/processed_data/SCE/
save(sce.lc, file=here("processed_data","SCE", "sce_raw_LC.rda"))



### Perform `emptyDrops()` for nuclei calling, within each sample ========================
  #   - Custom iteration: first perform `barcodeRanks()` to re-define the `lower=` param,
  #     which is by default set =100.  We think this is potentially too low and should
  #     capture that 'plateau' above the 'second knee point'.

  #   - This 'second knee point' can be estimated by restricting the fit space to (10,1e3)
  
sample.idx <- splitit(sce.lc$Sample)
e.out.custom <- list()

Sys.time()
    # [1] "2021-12-20 11:10:03 EST"
for(i in names(sample.idx)){
  
  # Re-define the `lower=` param by identifying the 'second knee point'
  newknee <- barcodeRanks(counts(sce.lc[ ,sample.idx[[i]]]),
                          fit.bounds=c(10,1e3))
  cat(paste0("'Second knee point' for ", i," is: ",
             metadata(newknee)$knee,"\n"))
  
  # Set `lower = newknee + 100` (to capture the 'plateau' mode of low UMI totals)
  cat(paste0("Simulating empty drops for: ", i,"... \n"))
    set.seed(109)
  e.out.custom[[i]] <- emptyDrops(counts(sce.lc[ ,sample.idx[[i]]]), niters=20000,
                                      lower = metadata(newknee)$knee + 100,
                                      BPPARAM=BiocParallel::MulticoreParam(2))
  cat(paste0("\n\t...Simulations complete. \n\t", Sys.time(), "\n\n\n"))
}
#'Second knee point' for Br2701_LC is: 136
#'Second knee point' for Br6522_LC is: 178
#'Second knee point' for Br8079_LC is: 155

for(i in 1:length(e.out.custom)){
  print(names(e.out.custom)[[i]])
  print(table(Signif = e.out.custom[[i]]$FDR <= 0.001,
              Limited = e.out.custom[[i]]$Limited))
}
    #[1] "Br2701_LC"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 23533     0
    #   TRUE    593  4464   # sum = 5057
    
    # [1] "Br6522_LC"
    #       Limited
    # Signif  FALSE TRUE
    #   FALSE  4799    0
    #   TRUE     88 3018    # sum = 3106
    
    # [1] "Br8079_LC"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 17194     0
    #   TRUE   1176  7070   # sum = 8246


## Save this with the raw SCE for interactive downstream analyses ===
README.custom <- "Object 'e.out.custom' is the output of emptyDrops() with lower= set to the quantified
  'second knee point' (+100) to better model empty/debris-containing droplets"
save(sce.lc, e.out.lc,
     e.out.custom, README.custom,
     file=here("processed_data","SCE", "sce_raw_LC.rda"))


## Filter based on this revised nuclei-calling approach and save into another .rda ===
BCs.keep <- unlist(
  lapply(e.out.custom, function(x){ rownames(x)[which(x$FDR <= 0.001)] }
         ))

sce.lc <- sce.lc[ ,BCs.keep]
save(sce.lc, file=here("processed_data","SCE", "sce_working_LC.rda"))


## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    #[1] "2021-12-20 12:07:57 EST"
proc.time()
    #    user   system  elapsed 
    #4249.531   43.424 6373.430 
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2021-12-20
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocIO                 1.4.0    2021-10-26 [2] Bioconductor
# BiocParallel         * 1.28.3   2021-12-09 [2] Bioconductor
# Biostrings             2.62.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
# crayon                 1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# DropletUtils         * 1.14.1   2021-11-08 [2] Bioconductor
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.0   2021-10-26 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicAlignments      1.30.0   2021-10-26 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# glue                   1.6.0    2021-12-17 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.0   2021-10-26 [2] Bioconductor
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# pillar                 1.6.4    2021-10-18 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
# RCurl                  1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
# restfulr               0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
# rhdf5                  2.38.0   2021-10-26 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rjson                  0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
# rlang                  0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# Rsamtools              2.10.0   2021-10-26 [2] Bioconductor
# rtracklayer          * 1.54.0   2021-10-26 [2] Bioconductor
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# scuttle                1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# XML                    3.99-0.8 2021-09-17 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/ntranngu/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────