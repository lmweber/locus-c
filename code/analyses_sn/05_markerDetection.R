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


# # Just for interactive exploration of some of these
# plotExpressionCustom(sce = sce.lc,
#                      exprs_values = "logcounts",
#                      #
#                      features = head(Oligo_A_markers,4),
#                      features_name = "custom-selected",
#                      anno_name = "cellType.merged",
#                      ncol=2, point_alpha=0.4, point_size=0.9,
#                      scales="free_y", swap_rownames="gene_name") +
#   ggtitle(label=paste0("Lionel's custom markers of interest")) +
#   theme(plot.title = element_text(size = 12),
#         axis.text.x = element_text(size=7))


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



## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
# 
proc.time()
#     user    system   elapsed 
#
options(width = 120)
session_info()



