### LC snRNA-seq analysis
### Marker detection with t-test
###     qsub -l bluejay,mf=76G,h_vmem=80G
### Initiated: MNT 07Apr2022

library(SingleCellExperiment)
library(scater)
library(scran)
library(scry)
library(batchelor)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)


### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

here()
# [1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"

source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")

## ===

load(here("processed_data","SCE", "sce_updated_LC.rda"), verbose=T)
    # sce.lc, annotationTab.lc, medianNon0.lc, hdgs.lc, cell_colors.lc

sce.lc
    # class: SingleCellExperiment 
    # dim: 36601 15642 
    # metadata(1): Samples
    # assays(3): counts binomial_deviance_residuals logcounts
    # rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
    #   ENSG00000277196
    # rowData names(7): source type ... gene_type binomial_deviance
    # colnames(15642): 3_AAACCCAAGTAGAATC-1 3_AAACCCAAGTGCACCC-1 ...
    #   2_TTTGTTGTCAGGGATG-1 2_TTTGTTGTCCCTCTCC-1
    # colData names(17): Sample Barcode ... mergedCluster.HC cellType.merged
    # reducedDimNames(5): GLMPCA_approx GLMPCA_MNN glmpca_mnn_50 UMAP TSNE
    # mainExpName: NULL
    # altExpNames(0):

table(sce.lc$cellType.merged)
    #ambig.lowNTx_A ambig.lowNTx_B ambig.lowNTx_C ambig.lowNTx_D ambig.lowNTx_E 
# 8657           2128            320            512            109 
# ambig.lowNTx_F          Astro     Endo.Mural        Excit_A        Excit_B 
# 95            430             99            295            305 
# Excit_C        Excit_D        Excit_E        Excit_F        Inhib_A 
# 70            267            168             32            414 
# Inhib_B        Inhib_C        Inhib_D        Inhib_E        Inhib_F 
# 222             99            148            137             50 
# Micro     Neuron.5HT      Neuron.NE          Oligo            OPC 
# 205             47             36            561            236

## doubletScore & sum distributions / cluster?
cellClust.idx <- splitit(sce.lc$cellType.merged)
sapply(cellClust.idx, function(x){round(quantile(sce.lc$doubletScore[x]), 2)})
    #

sapply(cellClust.idx, function(x){quantile(sce.lc$sum[x])})
    #



# First remove 'ambig.lowNTx' ('low N transcripts')-driven or -associated clusters
#     (this is a lot of nuclei btw...)
sce.lc <- sce.lc[ ,-grep("ambig.lowNTx_", sce.lc$cellType.merged)]
sce.lc$cellType.merged <- droplevels(sce.lc$cellType.merged)

# Remove 0 genes across all nuclei
sce.lc <- sce.lc[!rowSums(assay(sce.lc, "counts"))==0, ]  #33353, but after dropping, a 32223

## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.lc

assay(sce.lc, "logcounts") <- NULL
sizeFactors(sce.lc) <- NULL
sce.lc <- logNormCounts(sce.lc)

    
### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellClust.idx <- splitit(sce.lc$cellType.merged)
medianNon0.lc.19 <- lapply(cellClust.idx, function(x){
  apply(as.matrix(assay(sce.lc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.lc.19, table)
    #      Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A
    # FALSE 31823      30140   25715   26406   25868   26019   28262   25969   30767
    # TRUE    400       2083    6508    5817    6355    6204    3961    6254    1456
    # Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Micro Neuron.5HT Neuron.NE Oligo
    # FALSE   23496   28211   27000   28149   24313 31969      23810     26307 30183
    # TRUE     8727    4012    5223    4074    7910   254       8413      5916  2040
    # OPC
    # FALSE 31106
    # TRUE   1117




## Traditional t-test, pairwise ===
mod <- with(colData(sce.lc), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.lc.t.pw <- findMarkers(sce.lc, groups=sce.lc$cellType.merged,
                               assay.type="logcounts", design=mod, test="t",
                               direction="up", pval.type="all", full.stats=T)

sapply(markers.lc.t.pw, function(x){table(x$FDR<0.05)})
    #      Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A
# FALSE 31899      30853   32158   32182   31894   32199   32204   32039   32209
# TRUE    324       1370      65      41     329      24      19     184      14
# Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Micro Neuron.5HT Neuron.NE Oligo
# FALSE   32037   31923   32210   32190   32058 32016      31920     31945 31461
# TRUE      186     300      13      33     165   207        303       278   762
# OPC
# FALSE 32068
# TRUE    155




# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.lc.t.pw)){
  markers.lc.t.pw[[i]] <- cbind(markers.lc.t.pw[[i]],
                                medianNon0.lc.19[[i]][match(rownames(markers.lc.t.pw[[i]]),
                                                         names(medianNon0.lc.19[[i]]))])
  colnames(markers.lc.t.pw[[i]])[22] <- "non0median"
}

sapply(markers.lc.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    #     Astro.TRUE Endo.Mural.TRUE    Excit_A.TRUE    Excit_B.TRUE    Excit_C.TRUE 
    # 69             378              25              19             104 
    # Excit_D.TRUE    Excit_E.TRUE    Excit_F.TRUE    Inhib_A.TRUE    Inhib_B.TRUE 
    # 9               7              51               5             145 
    # Inhib_C.TRUE    Inhib_D.TRUE    Inhib_E.TRUE    Inhib_F.TRUE      Micro.TRUE 
    # 168               6               5              22              58 
    # Neuron.5HT.TRUE  Neuron.NE.TRUE      Oligo.TRUE        OPC.TRUE 
    # 96             101             483              63 


## Save these
save(markers.lc.t.pw, medianNon0.lc.19,
     file=here("processed_data","SCE",
               "markers-stats_LC-n3_findMarkers_19cellTypes.rda"))


    # As needed
    load(here("processed_data","SCE",
              "markers-stats_LC-n3_findMarkers_19cellTypes.rda"), verbose=T)


# Print these to pngs
markerList.t.pw <- lapply(markers.lc.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
 }
)

# Change to gene symbols
markerList.t.pw <- lapply(markerList.t.pw, function(x){
  rowData(sce.lc)$gene_name[match(x, rowData(sce.lc)$gene_id)]
  })

genes.top40.t <- lapply(markerList.t.pw, function(x){head(x, n=40)})

smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20]
left.set <- setdiff(names(genes.top40.t), smaller.set)
# Now remove those with no significant markers
smaller.set <- setdiff(smaller.set,
                       names(genes.top40.t)[lengths(genes.top40.t) == 0])

# Because scale_color_manual() has started keeping and printing unused
#     factor levels as 'NA' in the legend:
colors2print <- cell_colors.lc[-grep("ambig.lowNTx_", names(cell_colors.lc))]

# Smaller graphical window
#dir.create(here("plots","snRNA-seq","markers"))
    # Re-named Apr & before iteration to here("plots","snRNA-seq","zobsolete_Apr2022_markers")
for(i in smaller.set){
  png(here("plots","snRNA-seq","markers",
           paste0("LC_t_pairwise_topMarkers-", i, "_vlnPlots.png")), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         exprs_values = "logcounts",
                         features = head(genes.top40.t[[i]],n=4), 
                         features_name = i,
                         anno_name = "cellType.merged",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values = colors2print) +  
      ggtitle(label=paste0("LC ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}

# 20-40 markers
for(i in left.set){
  png(here("plots","snRNA-seq","markers",
           paste0("LC_t_pairwise_topMarkers-", i, "_vlnPlots.png")), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         exprs_values = "logcounts",
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType.merged",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values = colors2print) +  
      ggtitle(label=paste0("LC ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



### Cluster-vs-all-others single-nucleus-level iteration ========

# (Pre-process as above, if needed)
mod <- with(colData(sce.lc), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.lc.t.1vAll <- list()
for(i in levels(sce.lc$cellType.merged)){
  # Make temporary contrast
  sce.lc$contrast <- ifelse(sce.lc$cellType.merged==i, 1, 0)
  # Test cluster vs. all others
  markers.lc.t.1vAll[[i]] <- findMarkers(sce.lc, groups=sce.lc$contrast,
                                         assay.type="logcounts", design=mod, test="t",
                                         std.lfc=TRUE,
                                         direction="up", pval.type="all", full.stats=T)
}


class(markers.lc.t.1vAll[[1]])
    # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
    # -> we want the second entry, named "1"
    #    (for other purposes, might be interesting to look into that "0" entry, which
    #     is basically what genes are depleted in the cell type of interest)


## Do some reorganizing
markers.lc.t.1vAll <- lapply(markers.lc.t.1vAll, function(x){
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y){ y[ ,4] }) 
})

# Re-name std.lfc column and the entries; add non-0-median info
for(i in names(markers.lc.t.1vAll)){
  colnames(markers.lc.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.lc.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.lc.t.1vAll[[i]][["0"]] <- cbind(markers.lc.t.1vAll[[i]][["0"]],
                                          medianNon0.lc.19[[i]][match(rownames(markers.lc.t.1vAll[[i]][["0"]]),
                                                                   names(medianNon0.lc.19[[i]]))])
  colnames(markers.lc.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.lc.t.1vAll[[i]][["1"]] <- cbind(markers.lc.t.1vAll[[i]][["1"]],
                                          medianNon0.lc.19[[i]][match(rownames(markers.lc.t.1vAll[[i]][["1"]]),
                                                                   names(medianNon0.lc.19[[i]]))])
  colnames(markers.lc.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.lc.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}


## Save these
save(markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.19,
     file=here("processed_data","SCE",
               "markers-stats_LC-n3_findMarkers_19cellTypes.rda"))



## Marker numbers with the non-0-median filter
sapply(markers.lc.t.1vAll, function(x){
  table(x[[2]]$log.FDR < log(.001) & x[[2]]$non0median == TRUE)
})
    #       Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A
    # FALSE 32004      31172   28195   28413   31025   28941   31150   31635   31870
    # TRUE    219       1051    4028    3810    1198    3282    1073     588     353

    #       Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Micro Neuron.5HT Neuron.NE Oligo
    # FALSE   26518   30909   30528   30975   30950 32070      31128     31459 30996
    # TRUE     5705    1314    1695    1248    1273   153       1095       764  1227

    #         OPC
    # FALSE 31864
    # TRUE    359


## Print these to pngs
markerList.t.1vAll <- lapply(markers.lc.t.1vAll, function(x){
  rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
 }
)

# Change to gene symbols
markerList.t.1vAll <- lapply(markerList.t.1vAll, function(x){
  rowData(sce.lc)$gene_name[match(x, rowData(sce.lc)$gene_id)]
})

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(here("plots","snRNA-seq","markers",
           paste0("LC_t_1vALL_topMarkers-",i,"_vlnPlots.png")), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         exprs_values = "logcounts",
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType.merged",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values = colors2print) +  
      ggtitle(label=paste0("LC ", i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



## Write these top 40 lists to a csv ===
names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# Many of the PW results don't have 40 markers:
extend.idx <- names(which(lengths(markerList.t.pw) < 40))
for(i in extend.idx){
  markerList.t.pw[[i]] <- c(markerList.t.pw[[i]], rep("", 40-length(markerList.t.pw[[i]])))
}

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]
write.csv(top40genes, file=here("code","analyses_sn","top40genesLists_LC-n3_19cellTypes.csv"),
          row.names=FALSE)




### Some follow-up ================================================







## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    # [1] "2022-05-20 19:18:14 EDT"
proc.time()
    #     user    system   elapsed 
    # 1753.244   21.581 4234.483 
options(width = 120)
session_info()
    # ─ Session info ─────────────────────────────────────────────────────────────────
    # setting  value
    # version  R version 4.1.2 Patched (2021-11-04 r81138)
    # os       CentOS Linux 7 (Core)
    # system   x86_64, linux-gnu
    # ui       X11
    # language (EN)
    # collate  en_US.UTF-8
    # ctype    en_US.UTF-8
    # tz       US/Eastern
    # date     2022-05-20
    # pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
    # 
    # ─ Packages ─────────────────────────────────────────────────────────────────────
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
    # bluster                1.4.0    2021-10-26 [2] Bioconductor
    # cli                    3.3.0    2022-04-25 [2] CRAN (R 4.1.2)
    # cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
    # colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
    # cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
    # crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
    # DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
    # DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
    # DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
    # digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
    # dplyr                  1.0.9    2022-04-28 [2] CRAN (R 4.1.2)
    # dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
    # edgeR                  3.36.0   2021-10-26 [2] Bioconductor
    # ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
    # fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
    # farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
    # fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
    # gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
    # generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
    # GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
    # GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
    # GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
    # ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
    # ggplot2              * 3.3.6    2022-05-03 [2] CRAN (R 4.1.2)
    # ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
    # glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
    # googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
    # gridExtra            * 2.3      2017-09-09 [2] CRAN (R 4.1.0)
    # gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
    # here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
    # igraph                 1.3.1    2022-04-20 [2] CRAN (R 4.1.2)
    # IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
    # irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
    # jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
    # labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
    # lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
    # lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
    # limma                  3.50.3   2022-04-07 [2] Bioconductor
    # locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
    # magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
    # Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
    # MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
    # matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.1.2)
    # metapod                1.2.0    2021-10-26 [2] Bioconductor
    # munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
    # pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
    # pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
    # purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
    # R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
    # rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
    # RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.1.2)
    # Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
    # RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
    # ResidualMatrix         1.4.0    2021-10-26 [1] Bioconductor
    # rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
    # rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.1.2)
    # rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
    # S4Vectors            * 0.32.4   2022-03-24 [2] Bioconductor
    # ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
    # scales                 1.2.0    2022-04-13 [2] CRAN (R 4.1.2)
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
    # tibble                 3.1.7    2022-05-03 [2] CRAN (R 4.1.2)
    # tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
    # utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
    # vctrs                  0.4.1    2022-04-13 [2] CRAN (R 4.1.2)
    # vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
    # viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
    # viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
    # withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
    # XVector                0.34.0   2021-10-26 [2] Bioconductor
    # zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
    # 
    # [1] /users/ntranngu/R/4.1.x
    # [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
    # [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
    # 
    # ────────────────────────────────────────────────────────────────────────────────


