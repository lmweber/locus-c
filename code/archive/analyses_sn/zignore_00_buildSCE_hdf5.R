### LC snRNA-seq analysis
### Build SCE from raw (HDF5) matrices; perform nuclei calling
### qrsh -l bluejay (default mf=2G,h_vmem=2G,h_fsize=10G)
### Adopted from https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/build_basic_sce.R
### Initiated MNT 13Dec2021

library(SingleCellExperiment)
library(DropletUtils)
library(rtracklayer)
library(here)
library(sessioninfo)

here()
    # [1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"

### ===

# Basic sample info
sample.info <- read.table(here("fastq", "snRNA-seq","sample_libs_info.tsv"))
    # First column is how the Cell Ranger output is organized

sample.info$path <- file.path(
  here("processed_data", "cellranger"),
  sample.info$V1,
  "outs",
  "raw_feature_bc_matrix.h5"
)
stopifnot(all(file.exists(sample.info$path)))

## Build basic SCE (attempt with H5files; will add more colData later)
Sys.time()
    # [1] "2021-12-13 13:13:03 EST"
sce <- read10xCounts(
  samples = sample.info$path,
  sample.names = paste0(sample.info$V3,"_",sample.info$V2),
  type = "HDF5",
  col.names = TRUE
)
Sys.time()
    # [1] "2021-12-13 13:13:40 EST"
    # wow this is quick with the DelayedMatrix

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

assay(sce, "counts")
    # <36601 x 4857931> sparse matrix of class DelayedMatrix and type "integer":

## Size in Gb
lobstr::obj_size(sce) / 1024^3
    # 0.805 GB

# Re-name to sce.lc
sce.lc <- sce

## Save to /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/processed_data/SCE/
save(sce.lc, file=here("processed_data","SCE", "sce_raw_delayedArray.rda"))


    ## MNT comment: Leo made a good point that efforts to work on HDF5/DelayedArray representation
    #  15Dec2021    is less time-saving, since it still doesn't reduce computational load/time for
    #               many steps in analyses - just the I/O, as seen above. Will be more relevant
    #               if we ever get to the 1M+ nuclei stage, such as in Me'n data, etc.
    #            -> See 01_buildSCE_callNuclei.R for the more standard pipeline


## Reproducibility information
print('Reproducibility information:')
Sys.time()
    # [1] "2021-12-13 13:48:28 EST"
proc.time()
    #    user   system  elapsed 
    # 222.646    5.646 2515.850 
options(width = 120)
session_info()
    #─ Session info ──────────────────────────────────────────────────────────────────
    # setting  value
    # version  R version 4.1.2 Patched (2021-11-04 r81138)
    # os       CentOS Linux 7 (Core)
    # system   x86_64, linux-gnu
    # ui       X11
    # language (EN)
    # collate  en_US.UTF-8
    # ctype    en_US.UTF-8
    # tz       US/Eastern
    # date     2021-12-13
    # pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
    # 
    # ─ Packages ──────────────────────────────────────────────────────────────────────
    # package              * version  date (UTC) lib source
    # beachmat               2.10.0   2021-10-26 [2] Bioconductor
    # Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
    # BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
    # BiocIO                 1.4.0    2021-10-26 [2] Bioconductor
    # BiocParallel           1.28.3   2021-12-09 [2] Bioconductor
    # Biostrings             2.62.0   2021-10-26 [2] Bioconductor
    # bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
    # cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
    # crayon                 1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
    # DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
    # DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
    # dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
    # DropletUtils         * 1.14.1   2021-11-08 [2] Bioconductor
    # edgeR                  3.36.0   2021-10-26 [2] Bioconductor
    # GenomeInfoDb         * 1.30.0   2021-10-26 [2] Bioconductor
    # GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
    # GenomicAlignments      1.30.0   2021-10-26 [2] Bioconductor
    # GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
    # HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
    # here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
    # IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
    # lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
    # limma                  3.50.0   2021-10-26 [2] Bioconductor
    # lobstr                 1.1.1    2019-07-02 [2] CRAN (R 4.1.0)
    # locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
    # Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
    # MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
    # matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
    # R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
    # R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
    # R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
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
    # sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
    # SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
    # sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
    # SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
    # XML                    3.99-0.8 2021-09-17 [2] CRAN (R 4.1.2)
    # XVector                0.34.0   2021-10-26 [2] Bioconductor
    # yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
    # zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
    # 
    # [1] /users/ntranngu/R/4.1.x
    # [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
    # [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
    # 
    # ─────────────────────────────────────────────────────────────────────────────────
