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

table(sce.lc$cellType.merged)
    #    ambig.lowNTx            Astro       Endo.Mural          Excit_A 
    #            3519              430              173              538 
    #         Excit_B          Excit_C          Excit_D          Excit_E 
    #             305              295               70               34 
    #         Excit_F          Inhib_A          Inhib_B          Inhib_C 
    #              32              444              222               99 
    #         Inhib_D            Micro       Neuron.5HT Neuron.5HT_noDDC 
    #             255              235               47               64 
    #  Neuron.ambig_A   Neuron.ambig_B   Neuron.ambig_C   Neuron.ambig_D 
    #             183               95              299              129 
    #  Neuron.mixed_A   Neuron.mixed_B        Neuron.NE          Oligo_A 
    #            5244               73               36              561 
    #         Oligo_B              OPC 
    #            2024              236 

## doubletScore & sum distributions / cluster?
cellClust.idx <- splitit(sce.lc$cellType.merged)
sapply(cellClust.idx, function(x){round(quantile(sce.lc$doubletScore[x]), 2)})
    #      ambig.lowNTx Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E
    # 0%           0.00  0.00       0.02    0.08    0.64    0.62    0.51    0.28
    # 25%          0.02  0.05       0.26    0.97    1.16    0.94    1.29    1.08
    # 50%          0.11  0.13       0.48    1.26    1.41    1.06    1.45    1.53
    # 75%          0.28  0.31       0.80    1.59    1.62    1.25    1.54    2.64
    # 100%         5.60  3.64       2.45    4.41    3.29    4.68    2.27    3.75

    #      Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Micro Neuron.5HT Neuron.5HT_noDDC
    # 0%      0.95    0.22    0.54    0.80    0.00  0.01       0.62             0.20
    # 25%     1.13    1.00    0.90    1.30    0.15  0.05       0.91             0.64
    # 50%     1.32    1.30    1.14    1.52    0.35  0.21       1.47             0.81
    # 75%     1.62    1.56    1.35    1.83    0.81  0.41       1.50             1.41
    # 100%    2.00    4.63    2.53    1.88    4.59  5.10       2.97             3.73

    #      Neuron.ambig_A Neuron.ambig_B Neuron.ambig_C Neuron.ambig_D Neuron.mixed_A
    # 0%             0.02           0.14           0.03           0.05           0.00
    # 25%            0.18           0.56           0.79           0.35           0.69
    # 50%            0.45           1.03           1.07           0.58           1.13
    # 75%            1.22           1.33           1.52           0.87           1.60
    # 100%           6.27           2.15           3.25           2.02           6.28

    #      Neuron.mixed_B Neuron.NE Oligo_A Oligo_B  OPC
    # 0%             0.04      0.71    0.02    0.00 0.16
    # 25%            0.36      0.93    0.12    0.08 0.48
    # 50%            0.50      0.96    0.22    0.22 0.76
    # 75%            1.07      1.03    0.39    0.80 1.12
    # 100%           6.51      2.19    9.67   11.04 3.46

sapply(cellClust.idx, function(x){quantile(sce.lc$sum[x])})
    #      ambig.lowNTx    Astro Endo.Mural   Excit_A Excit_B Excit_C Excit_D
    # 0%          239.0   461.00        662    886.00    6596    4936    9426
    # 25%         475.5  1493.50       1449   3918.25   18892   19536   21411
    # 50%         754.0  2468.50       2436  10885.00   26572   30326   34727
    # 75%        1445.5  4307.75       6176  33163.75   40847   42630   52002
    # 100%      25442.0 27657.00      43741 133387.00  142023  106030  182780
  
    #        Excit_E Excit_F   Inhib_A   Inhib_B Inhib_C Inhib_D   Micro Neuron.5HT
    # 0%     1920.00  8771.0   4047.00  14333.00  4450.0   647.0   465.0    14002.0
    # 25%    2742.75 18266.5  11668.25  37479.00 11127.5  1515.0  1209.0    43467.5
    # 50%    5752.00 27073.5  19723.50  47500.50 16461.0  2253.0  2101.0    56183.0
    # 75%   43216.25 38486.0  30649.50  60213.75 22221.0  3670.5  3975.5    66372.0
    # 100% 150106.00 73021.0 167360.00 101398.00 74217.0 10652.0 37667.0   134968.0
  
    #      Neuron.5HT_noDDC Neuron.ambig_A Neuron.ambig_B Neuron.ambig_C
    # 0%            1375.00          578.0         1592.0         1023.0
    # 25%           2971.25         1422.0         2680.5         2837.5
    # 50%           4548.50         2220.0         4455.0         4544.0
    # 75%           8243.75         4688.5         7578.5         7001.5
    # 100%         38921.00        19594.0        20806.0        79135.0

    #      Neuron.ambig_D Neuron.mixed_A Neuron.mixed_B Neuron.NE Oligo_A  Oligo_B
    # 0%              823          276.0           2689  15101.00    3431   240.00
    # 25%            2031         2862.0           6744  24208.75    5510  1030.75
    # 50%            3087         8665.5          12224  32792.50    7768  1799.00
    # 75%            4735        28019.5          23154  53717.50   11893  3019.25
    # 100%          21886       146018.0          70563 161314.00   50182 97039.00

    #           OPC
    # 0%     627.00
    # 25%   2645.50
    # 50%   4549.00
    # 75%   9329.25
    # 100% 46062.00


## sizeFactor distribution
sapply(cellClust.idx, function(x){round(quantile(sce.lc$sizeFactor[x]), 2)})
    #     ambig.lowNTx Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E
    # 0%           0.03  0.05       0.08    0.10    0.78    0.57    1.09    0.22
    # 25%          0.06  0.17       0.17    0.46    2.22    2.29    2.48    0.32
    # 50%          0.09  0.28       0.28    1.26    3.15    3.51    4.03    0.67
    # 75%          0.17  0.50       0.72    3.86    4.79    4.90    6.03    5.12
    # 100%         2.93  3.21       5.18   15.81   16.47   12.20   21.66   17.41

    #      Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Micro Neuron.5HT Neuron.5HT_noDDC
    # 0%      1.02    0.48    1.65    0.52    0.08  0.06       1.62             0.16
    # 25%     2.16    1.37    4.31    1.28    0.18  0.14       5.02             0.34
    # 50%     3.17    2.32    5.46    1.90    0.27  0.25       6.46             0.53
    # 75%     4.53    3.56    6.93    2.57    0.43  0.46       7.67             0.96
    # 100%    8.65   19.83   11.66    8.61    1.26  4.46      15.53             4.48

    #      Neuron.ambig_A Neuron.ambig_B Neuron.ambig_C Neuron.ambig_D Neuron.mixed_A
    # 0%             0.07           0.19           0.12           0.10           0.03
    # 25%            0.17           0.32           0.34           0.24           0.34
    # 50%            0.26           0.53           0.53           0.37           1.02
    # 75%            0.55           0.90           0.82           0.56           3.26
    # 100%           2.32           2.47           9.18           2.52          17.30

    #      Neuron.mixed_B Neuron.NE Oligo_A Oligo_B  OPC
    # 0%             0.31      1.79    0.40    0.03 0.07
    # 25%            0.78      2.87    0.64    0.12 0.31
    # 50%            1.41      3.87    0.90    0.21 0.53
    # 75%            2.66      6.27    1.38    0.35 1.11
    # 100%           8.12     18.71    5.82   11.26 5.34



# Remove 0 genes across all nuclei
sce.lc <- sce.lc[!rowSums(assay(sce.lc, "counts"))==0, ]  #


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.lc

assay(sce.lc, "logcounts") <- NULL
sizeFactors(sce.lc) <- NULL
sce.lc <- logNormCounts(sce.lc)


### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellClust.idx <- splitit(sce.lc$cellType.merged)
medianNon0.lc.26 <- lapply(cellClust.idx, function(x){
  apply(as.matrix(assay(sce.lc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.lc.26, table)
    #      ambig.lowNTx Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E
    # FALSE        33293 32953      32856   30185   27536   26845   26998   30675
    # TRUE            60   400        497    3168    5817    6508    6355    2678
    #       Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Micro Neuron.5HT Neuron.5HT_noDDC
    # FALSE   27099   28694   24626   29341   32823 32971      24940            31911
    # TRUE     6254    4659    8727    4012     530   382       8413             1442
    #       Neuron.ambig_A Neuron.ambig_B Neuron.ambig_C Neuron.ambig_D
    # FALSE          32766          31986          32145          32593
    # TRUE             587           1367           1208            760
    #       Neuron.mixed_A Neuron.mixed_B Neuron.NE Oligo_A Oligo_B   OPC
    # FALSE          30683          29546     27437   31313   33076 32236
    # TRUE            2670           3807      5916    2040     277  1117




## Traditional t-test, pairwise ===
mod <- with(colData(sce.lc), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.lc.t.pw <- findMarkers(sce.lc, groups=sce.lc$cellType.merged,
                               assay.type="logcounts", design=mod, test="t",
                               direction="up", pval.type="all", full.stats=T)

sapply(markers.lc.t.pw, function(x){table(x$FDR<0.05)})
    #




# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.lc.t.pw)){
  markers.lc.t.pw[[i]] <- cbind(markers.lc.t.pw[[i]],
                                medianNon0.lc.26[[i]][match(rownames(markers.lc.t.pw[[i]]),
                                                         names(medianNon0.lc.26[[i]]))])
  colnames(markers.lc.t.pw[[i]])[29] <- "non0median"
}

sapply(markers.lc.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    #    ambig.lowNTx.NA          Astro.TRUE     Endo.Mural.TRUE        Excit_A.TRUE 
    #                 NA                  79                  62                   2 
    #       Excit_B.TRUE        Excit_C.TRUE        Excit_D.TRUE        Excit_E.TRUE 
    #                 15                  18                  87                  38 
    #       Excit_F.TRUE        Inhib_A.TRUE        Inhib_B.TRUE        Inhib_C.TRUE 
    #                 47                   1                   6                 142 
    #       Inhib_D.TRUE          Micro.TRUE     Neuron.5HT.TRUE Neuron.5HT_noDDC.NA 
    #                  3                 123                  65                  NA 
    #Neuron.ambig_A.TRUE Neuron.ambig_B.TRUE   Neuron.ambig_C.NA Neuron.ambig_D.TRUE 
    #                  1                   5                  NA                   1 
    #  Neuron.mixed_A.NA Neuron.mixed_B.TRUE      Neuron.NE.TRUE        Oligo_A.TRUE 
    #                 NA                1017                  95                 484 
    #         Oligo_B.NA            OPC.TRUE 
    #                 NA                  66 


## Save these
save(markers.lc.t.pw, medianNon0.lc.26,
     file=here("processed_data","SCE",
               "markers-stats_LC-n3_findMarkers_26cellTypes.rda"))


    # As needed
    load(here("processed_data","SCE",
              "markers-stats_LC-n3_findMarkers_26cellTypes.rda"), verbose=T)


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

# Smaller graphical window
#dir.create(here("plots","snRNA-seq","markers"))
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
      scale_color_manual(values = c(tableau20, tableau10medium)) +  
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
      scale_color_manual(values = c(tableau20, tableau10medium)) +  
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
                                          medianNon0.lc.26[[i]][match(rownames(markers.lc.t.1vAll[[i]][["0"]]),
                                                                   names(medianNon0.lc.26[[i]]))])
  colnames(markers.lc.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.lc.t.1vAll[[i]][["1"]] <- cbind(markers.lc.t.1vAll[[i]][["1"]],
                                          medianNon0.lc.26[[i]][match(rownames(markers.lc.t.1vAll[[i]][["1"]]),
                                                                   names(medianNon0.lc.26[[i]]))])
  colnames(markers.lc.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.lc.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}


## Save these
save(markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.26,
     file=here("processed_data","SCE",
               "markers-stats_LC-n3_findMarkers_26cellTypes.rda"))



## Marker numbers with the non-0-median filter
sapply(markers.lc.t.1vAll, function(x){
  table(x[[2]]$log.FDR < log(.001) & x[[2]]$non0median == TRUE)
})
    #       ambig.lowNTx Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E
    # FALSE        33340 33084      33084   30970   28579   29132   31971   33022
    # TRUE            13   269        269    2383    4774    4221    1382     331
    #       Excit_F Inhib_A Inhib_B Inhib_C Inhib_D Micro Neuron.5HT Neuron.5HT_noDDC
    # FALSE   32870   29283   29154   31883   33149 33072      32544            33153
    # TRUE      483    4070    4199    1470     204   281        809              200
    #       Neuron.ambig_A Neuron.ambig_B Neuron.ambig_C Neuron.ambig_D
    # FALSE          33136          32871          32703          33070
    # TRUE             217            482            650            283
    #       Neuron.mixed_A Neuron.mixed_B Neuron.NE Oligo_A Oligo_B   OPC
    # FALSE          30824          31173     32500   31738   33194 32805
    # TRUE            2529           2180       853    1615     159   548


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
      scale_color_manual(values = c(tableau20, tableau10medium)) +  
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
# 1vAllOthers
# Many of the PW results don't have 40 markers:
extend.idx <- names(which(lengths(markerList.t.1vAll) < 40))
for(i in extend.idx){
  markerList.t.1vAll[[i]] <- c(markerList.t.1vAll[[i]], rep("", 40-length(markerList.t.1vAll[[i]])))
}

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]
write.csv(top40genes, file=here("code","analyses_sn","top40genesLists_LC-n3_26cellTypes.csv"),
          row.names=FALSE)




### Some follow-up ================================================
load(here("processed_data","SCE",
          "markers-stats_LC-n3_findMarkers_26cellTypes.rda"), verbose=T)

## Oligo_B: there are no pw markers; sizeFactors IQR: [0.12, 0.35] ===
#     Observation:
#         - there are no pairwise markers, and all the '1vAll' markers are more expressed
#           in Oligo_A;  -> are there any depleted markers?

head(markers.lc.t.1vAll[["Oligo_B"]][["Oligo_B_depleted"]])
printThese <- head(rownames(markers.lc.t.1vAll[["Oligo_B"]][["Oligo_B_depleted"]]),n=12)
printThese <- rowData(sce.lc)$gene_name[match(printThese, rowData(sce.lc)$gene_id)]

# Just for interactive exploration of some of these
png(here("plots","snRNA-seq","markers",
         "temp_explore_vlnPlots_Oligo_B_depleted.png"), height=600, width=800)
print(
  plotExpressionCustom(sce = sce.lc,
                       exprs_values = "logcounts",
                       #
                       features = head(printThese,12),
                       features_name = "",
                       anno_name = "cellType.merged",
                       ncol=4, point_alpha=0.4, point_size=0.9,
                       scales="free_y", swap_rownames="gene_name") +
    scale_color_manual(values = c(tableau20, tableau10medium)) +  
    ggtitle(label=paste0("Custom-selected markers: Oligo_B depleted genes")) +
    theme(plot.title = element_text(size = 12),
          axis.text.x = element_text(size=7))
  )
dev.off()

    ## Observation: there doesn't seem to be any genes uniquely depleted in these
    #         - top depleted genes are co-depleted in _B AND _A (and often other glial pops)


### Let's actually print these for the others that have little/no pw markers...
othersToCheck <- c("Excit_A", "Inhib_A", "Inhib_D", "Neuron.5HT_noDDC",
            paste0("Neuron.ambig_",c("A","B","C","D")))

for(i in othersToCheck){
  
  printThese <- head(rownames(markers.lc.t.1vAll[[i]][[paste0(i,"_depleted")]]),n=12)
  printThese <- rowData(sce.lc)$gene_name[match(printThese, rowData(sce.lc)$gene_id)]
  
  # Just for interactive exploration of some of these
  png(here("plots","snRNA-seq","markers",
           paste0("temp_explore_vlnPlots_",i,"_depleted.png")), height=600, width=800)
    print(
      plotExpressionCustom(sce = sce.lc,
                           exprs_values = "logcounts",
                           #
                           features = head(printThese,12),
                           features_name = "",
                           anno_name = "cellType.merged",
                           ncol=4, point_alpha=0.4, point_size=0.9,
                           scales="free_y", swap_rownames="gene_name") +
        scale_color_manual(values = c(tableau20, tableau10medium)) +  
        ggtitle(label=paste0("Custom-selected markers: ",i," depleted genes")) +
        theme(plot.title = element_text(size = 12),
              axis.text.x = element_text(size=7))
      )
  dev.off()
  
}


### Looser marker test to characterize some of the 'Neuron.ambig's & 'noDDC' pop ======
  # Motivation:
  #   with the standard 'pval.type="all"' argument, it requires that a marker is
  #   stat. sig. compared to each and every other cluster.  Can loosen this, using
  #   'pval.type="some"', to try and identify more pairwise markers

# (follow set-up, above)

load(here("processed_data","SCE",
          "markers-stats_LC-n3_findMarkers_26cellTypes.rda"), verbose=T)
    # markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.26

## Traditional t-test, pairwise ===
mod <- with(colData(sce.lc), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.lc.t.pw.some <- findMarkers(sce.lc, groups=sce.lc$cellType.merged,
                               assay.type="logcounts", design=mod, test="t",
                               direction="up", pval.type="some", min.prop=0.5, # default for pval.type='some'
                               full.stats=T)

sapply(markers.lc.t.pw.some, function(x){table(x$FDR<0.05)})
#




# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.lc.t.pw.some)){
  markers.lc.t.pw.some[[i]] <- cbind(markers.lc.t.pw.some[[i]],
                                medianNon0.lc.26[[i]][match(rownames(markers.lc.t.pw.some[[i]]),
                                                            names(medianNon0.lc.26[[i]]))])
  colnames(markers.lc.t.pw.some[[i]])[29] <- "non0median"
}

sapply(markers.lc.t.pw.some, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    #   ambig.lowNTx.TRUE            Astro.TRUE       Endo.Mural.TRUE          Excit_A.TRUE 
    #                   6                   216                   214                   336 
    #        Excit_B.TRUE          Excit_C.TRUE          Excit_D.TRUE          Excit_E.TRUE 
    #                1925                   908                   796                   284 
    #        Excit_F.TRUE          Inhib_A.TRUE          Inhib_B.TRUE          Inhib_C.TRUE 
    #                 385                  1160                   888                   885 
    #        Inhib_D.TRUE            Micro.TRUE       Neuron.5HT.TRUE Neuron.5HT_noDDC.TRUE 
    #                  88                   229                   511                   106 
    # Neuron.ambig_A.TRUE   Neuron.ambig_B.TRUE   Neuron.ambig_C.TRUE   Neuron.ambig_D.TRUE 
    #                 115                   195                   148                   143 
    # Neuron.mixed_A.TRUE   Neuron.mixed_B.TRUE        Neuron.NE.TRUE          Oligo_A.TRUE 
    #                 139                  1945                   629                  1139 
    #        Oligo_B.TRUE              OPC.TRUE 
    #                 114                   334 


    ## Compared to, with 'pval.type="all"'
            sapply(markers.lc.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
            #    ambig.lowNTx.NA          Astro.TRUE     Endo.Mural.TRUE        Excit_A.TRUE 
            #                 NA                  79                  62                   2 
            #       Excit_B.TRUE        Excit_C.TRUE        Excit_D.TRUE        Excit_E.TRUE 
            #                 15                  18                  87                  38 
            #       Excit_F.TRUE        Inhib_A.TRUE        Inhib_B.TRUE        Inhib_C.TRUE 
            #                 47                   1                   6                 142 
            #       Inhib_D.TRUE          Micro.TRUE     Neuron.5HT.TRUE Neuron.5HT_noDDC.NA 
            #                  3                 123                  65                  NA 
            #Neuron.ambig_A.TRUE Neuron.ambig_B.TRUE   Neuron.ambig_C.NA Neuron.ambig_D.TRUE 
            #                  1                   5                  NA                   1 
            #  Neuron.mixed_A.NA Neuron.mixed_B.TRUE      Neuron.NE.TRUE        Oligo_A.TRUE 
            #                 NA                1017                  95                 484 
            #         Oligo_B.NA            OPC.TRUE 
            #                 NA                  66 

            
            
save(markers.lc.t.pw.some, medianNon0.lc.26,
     file=here("processed_data","SCE",
          "markers-stats_LC-n3_findMarkers-PW-some_26cellTypes.rda"))

printSomeMarkers <- c("Excit_A", "Inhib_A", "Inhib_D", "Neuron.5HT_noDDC",
                paste0("Neuron.ambig_",c("A","B","C","D")), "Oligo_B")
            
# Print these to pngs
markerList.t.pw <- lapply(markers.lc.t.pw.some, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
}
)
# Change to gene symbols
markerList.t.pw <- lapply(markerList.t.pw, function(x){
  rowData(sce.lc)$gene_name[match(x, rowData(sce.lc)$gene_id)]
})
head(markerList.t.pw[["Inhib_D"]], n=10)

    ## These actually look like the '1vAllOthers' markers...
    markerList.t.1vAll <- lapply(markers.lc.t.1vAll, function(x){
      rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
    }
    )
    
    # Change to gene symbols
    markerList.t.1vAll <- lapply(markerList.t.1vAll, function(x){
      rowData(sce.lc)$gene_name[match(x, rowData(sce.lc)$gene_id)]
    })
    
    sapply(printSomeMarkers, function(x){
      length(intersect(markerList.t.pw[[x]],
                       markerList.t.1vAll[[x]]))
    })
        #        Excit_A          Inhib_A          Inhib_D Neuron.5HT_noDDC   Neuron.ambig_A 
        #            335             1159               87              102              115 
        # Neuron.ambig_B   Neuron.ambig_C   Neuron.ambig_D          Oligo_B 
        #            194              147              143              113 

            # & looking above, this is probably ~98% of them overlapping...

# Just print 20 of these 'some' markers
for(i in printSomeMarkers){
  png(here("plots","snRNA-seq","markers",
           paste0("LC_t_pairwise-SOME_topMarkers-", i, "_vlnPlots.png")), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         exprs_values = "logcounts",
                         features = head(markerList.t.pw[[i]],n=20), 
                         features_name = i,
                         anno_name = "cellType.merged",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values = c(tableau20, tableau10medium)) +  
      ggtitle(label=paste0("LC ", i, " top markers: single-nucleus-level p.w.-SOME t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}






## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    # [1] "2022-05-05 21:43:37 EDT"
proc.time()
    #     user    system   elapsed 
    # 2072.283    87.553 16963.639 
options(width = 120)
session_info()
    #─ Session info ────────────────────────────────────────────────────────────────────────────────
    # setting  value
    # version  R version 4.1.2 Patched (2021-11-04 r81138)
    # os       CentOS Linux 7 (Core)
    # system   x86_64, linux-gnu
    # ui       X11
    # language (EN)
    # collate  en_US.UTF-8
    # ctype    en_US.UTF-8
    # tz       US/Eastern
    # date     2022-05-05
    # pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
    # 
    # ─ Packages ────────────────────────────────────────────────────────────────────────────────────
    # package              * version  date (UTC) lib source
    # assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
    # batchelor            * 1.10.0   2021-10-26 [1] Bioconductor
    # beachmat               2.10.0   2021-10-26 [2] Bioconductor
    # beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
    # Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
    # BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
    # BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
    # BiocParallel           1.28.3   2021-12-09 [2] Bioconductor
    # BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
    # bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
    # bluster              * 1.4.0    2021-10-26 [2] Bioconductor
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
    # DropletUtils         * 1.14.2   2022-01-09 [2] Bioconductor
    # dynamicTreeCut       * 1.63-1   2016-03-11 [1] CRAN (R 4.1.2)
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
    # HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
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
    # pheatmap             * 1.0.12   2019-01-04 [2] CRAN (R 4.1.0)
    # pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
    # pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
    # purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
    # R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
    # R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
    # R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
    # R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
    # rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
    # RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.1.2)
    # Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
    # RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
    # ResidualMatrix         1.4.0    2021-10-26 [1] Bioconductor
    # rhdf5                  2.38.1   2022-03-10 [2] Bioconductor
    # rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
    # Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
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
    # ───────────────────────────────────────────────────────────────────────────────────────────────


