### LC snRNA-seq analysis
### Mito rate QC & doublet assessment
### Initiated: MNT 20Dec2021

library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(rtracklayer)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)

here()
    #[1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"



### Mito rate QC ============

# Load nuclei-called SCE
load(here("processed_data","SCE", "sce_working_LC.rda"), verbose=T)
    # sce.lc

# Compute mito rates at the sample-level
sample.idx <- splitit(sce.lc$Sample)

stats <- list()
for(i in names(sample.idx)){
  stats[[i]] <- perCellQCMetrics(sce.lc[ ,sample.idx[[i]]],
                                 subsets=list(Mito = which(seqnames(sce.lc) == "chrM")))
}
names(stats) <- names(sample.idx)


### Trick: Add a pseudo-count==1 for a 'MT transcript' ===
# Note: This was implemented because we realized samples with mito rate distributions that
#       were 'clean' and tightly distributed about 0 would yield a 3x MAD = 0, thus over-penalizing
#       nuclei even if they had a single MT transcript (throwing out upwards of 50% of the sample)

# First check computation of mito percent:
table(stats[[1]]$subsets_Mito_percent == (stats[[1]]$subsets_Mito_sum/stats[[1]]$sum)*100)
    # All TRUE

test.stats <- stats

for(i in 1:length(test.stats)){
  test.stats[[i]]$pseudo_subsets_Mito_sum <- test.stats[[i]]$subsets_Mito_sum + 1
  test.stats[[i]]$pseudo_subsets_Mito_percent <- test.stats[[i]]$pseudo_subsets_Mito_sum / (test.stats[[i]]$sum+1) * 100
}

## Lapply: MAD approach for mito rate thresholding
pseudo.high.mito <- lapply(test.stats, function(x) isOutlier(x$pseudo_subsets_Mito_percent, nmads=3, type="higher"))
pseudo.high.mito.table <- lapply(pseudo.high.mito, table)
# Percent dropped
sapply(pseudo.high.mito.table, function(x) round(x[2]/sum(x), 3))
    #Br2701_LC.TRUE Br6522_LC.TRUE Br8079_LC.TRUE 
    #         0.025          0.132          0.028 

# Thresholds
sapply(pseudo.high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 4)})
    #Br2701_LC.higher Br6522_LC.higher Br8079_LC.higher 
    #          8.1356           2.8406           8.5348 

    ## In this case, wouldn't have affected the samples too much bc their  
     # distributions aren't that tight ~0 as with some much cleaner samples ====
    high.mito <- lapply(stats, function(x) isOutlier(x$subsets_Mito_percent, nmads=3, type="higher"))
    high.mito.table <- lapply(high.mito, table)
    # Percent dropped
    sapply(high.mito.table, function(x) round(x[2]/sum(x), 3))
        #Br2701_LC.TRUE Br6522_LC.TRUE Br8079_LC.TRUE 
        #         0.026          0.134          0.027 
    
    # Thresholds
    sapply(high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 4)})
        #Br2701_LC.higher Br6522_LC.higher Br8079_LC.higher 
        #          7.9820           2.7090           8.3752 
    
     # end side note ====

## Bind [true] stats to colData
 # (we'll just keep the 'pseudo' result since this was made/meant to work with a range of data)
table(rownames(do.call("rbind", stats)) == colnames(sce.lc))
    # all 16409 TRUE
    
colData(sce.lc) <- cbind(colData(sce.lc),
                         do.call("rbind", stats),
                         do.call("c", pseudo.high.mito))
colnames(colData(sce.lc))[9] <- "high.mito"

# $sum == $total
sce.lc$total <- NULL

# Store original for comparison/plotting
sce.lc.unfiltered <- sce.lc
sce.lc <- sce.lc[ ,!sce.lc$high.mito]


## Plot some metrics
mitoCutoffs <- sapply(pseudo.high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 3)})
names(mitoCutoffs) <- gsub(".higher","", names(mitoCutoffs))

pdf(here("plots","snRNA-seq","LC-n3_QCmetrics_high-mitoColored_MNT.pdf"), height=4)
for(i in names(sample.idx)){
  grid.arrange(
    plotColData(sce.lc.unfiltered[ ,sample.idx[[i]]], y="sum", colour_by="high.mito", point_alpha=0.4) +
      scale_y_log10() + ggtitle(paste0("Total count: ", i)),
    plotColData(sce.lc.unfiltered[ ,sample.idx[[i]]], y="detected", colour_by="high.mito", point_alpha=0.4) +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce.lc.unfiltered[ ,sample.idx[[i]]], y="subsets_Mito_percent",
                colour_by="high.mito", point_alpha=0.4) +
      ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs[i],")")),
    ncol=3
  )
  # Mito rate vs n detected features
  print(
    plotColData(sce.lc.unfiltered[ ,sample.idx[[i]]], x="detected", y="subsets_Mito_percent",
                colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
      ggtitle(paste0("Sample: ", i,
                     ";     pre-QC nNuclei: ", ncol(sce.lc.unfiltered[ ,sce.lc.unfiltered$Sample==i]),";      ",
                     "nNuclei kept: ", ncol(sce.lc[ ,sce.lc$Sample==i])," (",
                     round(ncol(sce.lc[ ,sce.lc$Sample==i]) /
                             ncol(sce.lc.unfiltered[ ,sce.lc.unfiltered$Sample==i]) * 100, 2), "%)"
      ))
  )
  # Detected features vs total count
  print(
    plotColData(sce.lc.unfiltered[ ,sample.idx[[i]]], x="sum", y="detected",
                colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
      ggtitle(paste0("Sample: ", i,
                     ";     pre-QC nNuclei: ", ncol(sce.lc.unfiltered[ ,sce.lc.unfiltered$Sample==i]),";      ",
                     "nNuclei kept: ", ncol(sce.lc[ ,sce.lc$Sample==i])," (",
                     round(ncol(sce.lc[ ,sce.lc$Sample==i]) /
                             ncol(sce.lc.unfiltered[ ,sce.lc.unfiltered$Sample==i]) * 100, 2), "%)"
      ))
  )
}
dev.off()


# Save
save(sce.lc, file=here("processed_data","SCE", "sce_working_LC.rda"))


### Doublet score computation (no filtering here) ============




## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
#[1] "2021-12-20 14:28:37 EST"
proc.time()
#    user   system  elapsed 
# 181.646   11.122 6249.428 
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
# beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocIO                 1.4.0    2021-10-26 [2] Bioconductor
# BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
# BiocParallel         * 1.28.3   2021-12-09 [2] Bioconductor
# BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
# Biostrings             2.62.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
# colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
# cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
# crayon                 1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# DropletUtils         * 1.14.1   2021-11-08 [2] Bioconductor
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.0   2021-10-26 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicAlignments      1.30.0   2021-10-26 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                   1.6.0    2021-12-17 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# gridExtra            * 2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.0   2021-10-26 [2] Bioconductor
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
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
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# rtracklayer          * 1.54.0   2021-10-26 [2] Bioconductor
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scater               * 1.22.0   2021-10-26 [2] Bioconductor
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.4.3    2021-11-30 [2] CRAN (R 4.1.2)
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
