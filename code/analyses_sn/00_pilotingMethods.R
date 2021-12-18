### LC snRNA-seq analysis
### Build SCE from raw count matrices; perform nuclei calling
### qrsh -l bluejay,mf=32G,h_vmem=34G -pe local 2
#     (which points to this script with `R CMD BATCH`)
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


### MNT test 17Dec2021 ==================
load(here("processed_data","SCE","sce_raw_LC.rda"), verbose=T)
# sce.lc, e.out.lc

sample.idx <- splitit(sce.lc$Sample)
i <- "Br6522_LC"


## Create barcode rank plots from `DropletUtils::barcodeRanks()`; get stats ===
 #   MNT: trick the fitting to identify the 'second knee', as this will be used as
 #        `lower=` for `emptyDrops()`
BRdf.list <- list()
for(i in names(sample.idx)){
  BRdf.list[[i]] <- barcodeRanks(counts(sce.lc[ ,sample.idx[[i]]]),
                                 fit.bounds=c(10,1e3))
}

sapply(BRdf.list, metadata)
    # (default params without defined `fit.bounds`:)
    #             Br2701_LC Br6522_LC Br8079_LC
    # knee        27658     193       196      
    # inflection  8999      106       102
        # So interestingly this would have detected that 'second knee' for two samples

    # With `fit.bounds=c(10,1e3)`
    #             Br2701_LC Br6522_LC Br8079_LC
    # knee        136       178       155      
    # inflection  8999      106       102   *"The derivative is computed directly from
                                            # all points on the curve with total counts
                                            # greater than lower" - hence the 8999



# Plotting (adapted from the manual under `barcodeRanks()`) ===
for(i in names(BRdf.list)){
#png(here("plots","snRNA-seq",paste0("barcodeRankPlot_",i,"_dropletUtils-defaults.png")))
#png(here("plots","snRNA-seq",paste0("barcodeRankPlot_",i,"_dropletUtils-w-fitbounds.png")))
png(here("plots","snRNA-seq",paste0("barcodeRankPlot_",i,"_dropletUtils-w-fitbounds_plus100UMIs.png")))
  
  br.out <- BRdf.list[[i]]
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Barcode Rank", ylab="Total UMIs",
       main=paste0("Barcode rank plot for: ",i,"\n( fit.bounds=c(10,1e3) )"),
       cex.axis=0.8, cex.lab=1.2, las=1)
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    # k_2 + 100 to capture the huge mode of seeming empty drops (the 'plateau region')
    abline(h=metadata(br.out)$knee + 100, col="darkblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "darkblue", "forestgreen"),
         legend=c("knee_2", "knee+100", "inflection"))
dev.off()
}

## Apply these
#customLowers <- sapply(BRdf.list, function(x){metadata(x)$knee})
customLowers.100 <- sapply(BRdf.list, function(x){metadata(x)$knee + 100})
#e.out.custom <- list()
e.out.custom.100 <- list()
Sys.time()
    # [1] "2021-12-17 17:47:24 EST"
for(i in names(sample.idx)){
  cat(paste0("Simulating empty drops for: ", i,"... \n"))
  
  ## Run emptyDrops with this data-defined `lower=` param 
  set.seed(109)
  #e.out.custom[[i]] <- emptyDrops(counts(sce.lc[ ,sample.idx[[i]]]), niters=20000,
                              #lower=customLowers[i],
  e.out.custom.100[[i]] <- emptyDrops(counts(sce.lc[ ,sample.idx[[i]]]), niters=20000,
                              lower=customLowers.100[i],
                              BPPARAM=BiocParallel::MulticoreParam(2))
  cat(paste0("\n\t...Simulations complete. \n\t", Sys.time(), "\n\n\n"))
  Sys.time()
  #        ...Simulations complete. 
  # 2021-12-17 17:57:50
  # 2021-12-17 18:06:32
  # 2021-12-17 18:17:40   # ~10min / sample
}

for(i in 1:length(e.out.custom)){
  print(names(e.out.custom)[[i]])
  print(table(Signif = e.out.custom[[i]]$FDR <= 0.001,
              Limited = e.out.custom[[i]]$Limited))
}
    # [1] "Br2701_LC"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 65682     0
    #   TRUE   2173  9706   # sum = 11879

    # [1] "Br6522_LC"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 66676  2032   - looks much better but now many p-vals need more iterations
    #   TRUE      0   823   # (sum could be up to 2855)

    # [1] "Br8079_LC"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 66956     0
    #   TRUE    816  8930   # sum = 9746

for(i in 1:length(e.out.custom.100)){
  print(names(e.out.custom.100)[[i]])
  print(table(Signif = e.out.custom.100[[i]]$FDR <= 0.001,
              Limited = e.out.custom.100[[i]]$Limited))
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


## Save for now
save(#e.out.custom,
     e.out.custom.100,
     file=here("processed_data","SCE", #"emptyDrops-by-sample_customLower_k2.rda"
               "emptyDrops-by-sample_customLower_k2plus100.rda"
               ))



## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
# [1] "2021-12-17 19:29:10 EST"
proc.time()
#     user    system   elapsed 
# 6299.080   107.997 26844.828 
options(width = 120)
session_info()
#─ Session info ───────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2021-12-17
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
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
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

