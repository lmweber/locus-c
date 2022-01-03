### LC snRNA-seq analysis
### Feature selection with deviance & clustering
###     qrsh -l bluejay,mf=
### Initiated: MNT 31Dec2021

library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(scry)
library(batchelor)
library(bluster)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)

here()
    # [1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"



### Feature selection with deviance residuals =====

# Load working SCE
load(here("processed_data","SCE", "sce_working_LC.rda"), verbose=T)

sce.lc
    #class: SingleCellExperiment 
    # dim: 36601 15642 
    # metadata(1): Samples
    # assays(1): counts
    # rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
    #   ENSG00000277196
    # rowData names(6): source type ... gene_name gene_type
    # colnames(15642): 3_AAACCCAAGTAGAATC-1 3_AAACCCAAGTGCACCC-1 ...
    #   2_TTTGTTGTCAGGGATG-1 2_TTTGTTGTCCCTCTCC-1
    # colData names(9): Sample Barcode ... high.mito doubletScore
    # reducedDimNames(0):
    # mainExpName: NULL
    # altExpNames(0):

Sys.time()
    #[1] "2021-12-31 14:36:36 EST"
sce.lc <- devianceFeatureSelection(sce.lc,
                                   assay="counts", fam="binomial", sorted=F,
                                      # these are default params btw
                                   batch=as.factor(sce.lc$Sample))
Sys.time()
    #[1] "2021-12-31 14:38:08 EST"

# Btw:
table(is.na(rowData(sce.lc)$binomial_deviance))
    # FALSE  TRUE 
    # 29556  7045

# Observe:
plot(sort(rowData(sce.lc)$binomial_deviance, decreasing=T),
     type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance: LC (n=3)")
abline(v=2000, lty=2, col="red")


## Approximate GLM-PCA with PCA on null deviance residuals
Sys.time()
    #[1] "2021-12-31 15:13:20 EST"
sce.lc <- nullResiduals(sce.lc, assay="counts", fam="binomial",  # default params
                        type="deviance")#, batch=as.factor(sce.lc$Sample))
    # Errored out trying to pass `batch=` arg:
        # Error in res[, idx] <- .null_residuals(m[, idx], fam = fam, type = type,  : 
        # number of items to replace is not a multiple of replacement length
    
    # When leaving that arg out - successfully runs
        # [In addition:] Warning message:
        # In sqrt(x@x) : NaNs produced
Sys.time()
    #[1] "2021-12-31 15:20:56 EST"

    # The above adds a new computed assay called 'binomial_deviance_residuals'


## Take the top 2000 highly-deviant genes ('HDG's) and run PCA in that space
hdgs.lc <- rownames(sce.lc)[order(rowData(sce.lc)$binomial_deviance, decreasing=T)][1:2000]
hdgs.symbols <- rowData(sce.lc)$gene_name[match(hdgs.lc, rowData(sce.lc)$gene_id)]

Sys.time()
    #[1] "2021-12-31 16:20:14 EST"
set.seed(109)
sce.lc <-  runPCA(sce.lc, exprs_values="binomial_deviance_residuals",
                  subset_row=hdgs.lc, ncomponents=100,
                  name="GLMPCA_approx",
                  BSPARAM=BiocSingular::IrlbaParam())
Sys.time()
    #[1] "2021-12-31 16:24:58 EST"

# Visualize
plotReducedDim(sce.lc, dimred="GLMPCA_approx", colour_by="Sample",
               ncomponents=4, point_alpha=0.3, point_size=1.5)


## See how this changes with batchelor::reducedMNN ===
Sys.time()
    #[1] "2021-12-31 17:00:16 EST"
glmpca.mnn <- reducedMNN(reducedDim(sce.lc, "GLMPCA_approx"),
                         batch=as.factor(sce.lc$Sample),
                         merge.order=c("Br8079_LC", "Br2701_LC", "Br6522_LC")
                         )
Sys.time()
    #[1] "2021-12-31 17:00:37 EST"

# Store this
reducedDim(sce.lc, "GLMPCA_MNN") <- glmpca.mnn$corrected

# Vizualize
plotReducedDim(sce.lc, dimred="GLMPCA_MNN", colour_by="Sample",
               ncomponents=4, point_alpha=0.3, point_size=1.5)
    # Doesn't look too amazingly better--at least in the top 4 dims.  Might have to
    # do some clustering in both PC spaces to assess


# Save for now
save(sce.lc, file=here("processed_data","SCE", "sce_working_LC.rda"))




### Clustering ====================
  # Perform graph-based clustering, as in Tran-Maynard, et al. Neuron 2021

  # Just perform on top 50 PCs, since for now we're assessing GLM PCA
  # [and then with MNN on reduced dims]


## Make some other reduced dims for visualization ===

    ## Testing phase ====
    sce.test <- sce.lc
    
    # Take top 50 PCs for both to quicken computation
    reducedDim(sce.test, "glmpca_approx_50") <- reducedDim(sce.test, "GLMPCA_approx")[ ,1:50]
    reducedDim(sce.test, "glmpca_mnn_50") <- reducedDim(sce.test, "GLMPCA_MNN")[ ,1:50]
        # oh didn't have to use this - can just use 'n_dimred'
    
    ## Without MNN on GLMPCA ===
    # UMAP
    set.seed(109)
    sce.test.approx <- runUMAP(sce.test, dimred="glmpca_approx_50",
                               name="UMAP")
    # t-SNE
    set.seed(109)
    sce.test.approx <- runTSNE(sce.test.approx, dimred="glmpca_approx_50",
                               name="TSNE")
    
    # k-means clustering with k=13
    #   (== average of Cell Ranger's sample-specific n graph-based clusters)
    set.seed(109)
    clust.kmeans <- clusterCells(sce.test.approx, use.dimred="glmpca_approx_50", 
                                 BLUSPARAM=KmeansParam(centers=13))
    table(clust.kmeans)
        #   1    2    3    4    5    6    7    8    9   10   11   12   13 
        # 196  843  635  173  656  261  194 1187  725 1234 4177 4442  919
    sce.test.approx$kmeans.13 <- clust.kmeans
    
    ## WITH MNN on GLMPCA ===
    # UMAP
    set.seed(109)
    sce.test.mnn <- runUMAP(sce.test, dimred="glmpca_mnn_50",
                            name="UMAP")
    # t-SNE
    set.seed(109)
    sce.test.mnn <- runTSNE(sce.test.mnn, dimred="glmpca_mnn_50",
                            name="TSNE")
    
    # k-means (with k=13)
    set.seed(109)
    clust.kmeans <- clusterCells(sce.test.approx, use.dimred="glmpca_mnn_50", 
                                 BLUSPARAM=KmeansParam(centers=13))
    table(clust.kmeans)
        #   1    2    3    4    5    6    7    8    9   10   11   12   13 
        #1064  259  654 1674  637  268  190  121  973 1245  179 8103  275 
    sce.test.mnn$kmeans.13 <- clust.kmeans

    ## Add Tran-Maynard, et al. method (just use HDGs from above) =====
        sce.test <- multiBatchNorm(sce.test, batch=sce.test$Sample)
        set.seed(109)
        mnn.hold <-  fastMNN(sce.test, batch=as.factor(sce.test$Sample),
                             merge.order=c("Br8079_LC", "Br2701_LC", "Br6522_LC"),
                             subset.row=hdgs.lc, d=100,
                             correct.all=TRUE, get.variance=TRUE,
                             BSPARAM=BiocSingular::IrlbaParam())

        table(colnames(mnn.hold) == colnames(sce.test))  # all TRUE
        table(mnn.hold$batch == sce.test$Sample) # all TRUE
        
        # Add them to the SCE
        reducedDim(sce.test, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
        metadata(sce.test) <- metadata(mnn.hold)
        
        # UMAP
        set.seed(109)
        sce.test <- runUMAP(sce.test, dimred="PCA_corrected",
                            n_dimred=50, name="UMAP")
        # t-SNE
        set.seed(109)
        sce.test <- runTSNE(sce.test, dimred="PCA_corrected",
                                   n_dimred=50, name="TSNE")
        
        # k-means (with k=13)
        reducedDim(sce.test, "PCA_corrected_50") <- reducedDim(sce.test, "PCA_corrected")[ ,1:50]
        set.seed(109)
        clust.kmeans <- clusterCells(sce.test, use.dimred="PCA_corrected_50",
                                     BLUSPARAM=KmeansParam(centers=13))
        table(clust.kmeans)
            #   1    2    3    4    5    6    7    8    9   10   11   12   13 
            #2713  246  647  580  600  598  522 2091 4305  701  935 1292  412
        sce.test$kmeans.13 <- clust.kmeans
                
        # end added Tran-Maynard et al. chunk ====

    
    # How do these look?
    pdf(here("plots","snRNA-seq","LC-n3_reducedDims_GLMPCA_reducedMNN-test.pdf"), height=5, width=5)
        ## TSNE
        plotReducedDim(sce.test.approx, dimred="TSNE", colour_by="Sample",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle("t-SNE on GLMPCA (top 50)")
        plotReducedDim(sce.test.approx, dimred="TSNE", colour_by="kmeans.13",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle("t-SNE on GLMPCA (top 50)")
        
        plotReducedDim(sce.test.mnn, dimred="TSNE", colour_by="Sample",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle("t-SNE on GLMPCA-MNN (top 50)")
        plotReducedDim(sce.test.mnn, dimred="TSNE", colour_by="kmeans.13",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle("t-SNE on GLMPCA-MNN (top 50)")
        
        # Tran-Maynard
        plotReducedDim(sce.test, dimred="TSNE", colour_by="Sample",
                       point_alpha=0.2, point_size=0.8) +
            ggtitle("t-SNE on fastMNN-corrected (top 50) PCs\nfrom log-norm. counts")
        plotReducedDim(sce.test, dimred="TSNE", colour_by="kmeans.13",
                       point_alpha=0.2, point_size=0.8) +
            ggtitle("t-SNE on fastMNN-corrected (top 50) PCs\nfrom log-norm. counts")
        
        
        ## UMAP
        plotReducedDim(sce.test.approx, dimred="UMAP", colour_by="Sample",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on GLMPCA (top 50)")
        plotReducedDim(sce.test.approx, dimred="UMAP", colour_by="kmeans.13",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on GLMPCA (top 50)")
        
        plotReducedDim(sce.test.mnn, dimred="UMAP", colour_by="Sample",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on GLMPCA-MNN (top 50)")
        plotReducedDim(sce.test.mnn, dimred="UMAP", colour_by="kmeans.13",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on GLMPCA-MNN (top 50)")
        
        # Tran-Maynard
        plotReducedDim(sce.test, dimred="UMAP", colour_by="Sample",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on fastMNN-corrected (top 50) PCs\nfrom log-norm. counts")
        plotReducedDim(sce.test, dimred="UMAP", colour_by="kmeans.13",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on fastMNN-corrected (top 50) PCs\nfrom log-norm. counts")
        
    dev.off()
    
    
    ## Save these various test iterations for further exploration
    save(sce.test, sce.test.approx, sce.test.mnn, hdgs.lc,
         file=here("processed_data","SCE", "sce_reducedDim-tests_LC.rda"))
    
    ### END TESTING ========









## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
#[1] "2022-01-03 12:28:58 EST"
proc.time()
#     user    system   elapsed 
# 2518.426    23.341 11535.588 
options(width = 120)
session_info()
# ─ Session info ────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-01-03
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# batchelor            * 1.10.0   2021-10-26 [1] Bioconductor
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
# BiocParallel         * 1.28.3   2021-12-09 [2] Bioconductor
# BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# bluster              * 1.4.0    2021-10-26 [2] Bioconductor
# cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.2)
# codetools              0.2-18   2020-11-04 [3] CRAN (R 4.1.2)
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
# igraph                 1.2.10   2021-12-15 [2] CRAN (R 4.1.2)
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
# metapod                1.2.0    2021-10-26 [2] Bioconductor
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
# RcppAnnoy              0.0.19   2021-07-30 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
# ResidualMatrix         1.4.0    2021-10-26 [1] Bioconductor
# rhdf5                  2.38.0   2021-10-26 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rlang                  0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# RSpectra               0.16-0   2019-12-01 [2] CRAN (R 4.1.0)
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# Rtsne                  0.15     2018-11-10 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scater               * 1.22.0   2021-10-26 [2] Bioconductor
# scran                * 1.22.1   2021-11-14 [2] Bioconductor
# scry                 * 1.6.0    2021-10-26 [2] Bioconductor
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# uwot                   0.1.11   2021-12-02 [2] CRAN (R 4.1.2)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.4.3    2021-11-30 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/ntranngu/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ───────────────────────────────────────────────────────────────────────────

