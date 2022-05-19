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
    #ambig.lowNTx_A ambig.lowNTx_B ambig.lowNTx_C ambig.lowNTx_D ambig.lowNTx_E 
    #          8657           2128            320            512            109 
    #         Astro     Endo.Mural        Excit_A        Excit_B        Excit_C 
    #           430             99            295            305             70 
    #       Excit_D        Excit_E        Excit_F        Inhib_A        Inhib_B 
    #           267             32            168            844            222 
    #       Inhib_C          Micro     Neuron.5HT      Neuron.NE          Oligo 
    #            99            205             47             36            561 
    #           OPC 
    #           236 

## doubletScore & sum distributions / cluster?
cellClust.idx <- splitit(sce.lc$cellType.merged)
sapply(cellClust.idx, function(x){round(quantile(sce.lc$doubletScore[x]), 2)})
    #      ambig.lowNTx_A ambig.lowNTx_B ambig.lowNTx_C ambig.lowNTx_D ambig.lowNTx_E
    # 0%             0.00           0.00           0.02           0.01           0.05
    # 25%            0.18           0.08           0.30           0.13           0.46
    # 50%            0.71           0.23           0.83           0.24           0.68
    # 75%            1.31           0.79           1.33           0.53           1.11
    # 100%           6.28          11.04           6.27           6.51           3.75

    #      Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A
    # 0%    0.00       0.14    0.62    0.64    0.51    0.94    0.95    0.21    0.00
    # 25%   0.05       0.26    0.94    1.16    1.29    1.21    1.13    0.85    0.64
    # 50%   0.13       0.44    1.06    1.41    1.45    1.39    1.32    1.14    1.11
    # 75%   0.31       0.79    1.25    1.62    1.54    1.69    1.62    1.51    1.44
    # 100%  3.64       2.20    4.68    3.29    2.27    4.41    2.00    3.14    4.63

    #      Inhib_B Inhib_C Micro Neuron.5HT Neuron.NE Oligo  OPC
    # 0%      0.54    0.80  0.01       0.62      0.71  0.02 0.16
    # 25%     0.90    1.30  0.05       0.91      0.93  0.12 0.48
    # 50%     1.14    1.52  0.17       1.47      0.96  0.22 0.76
    # 75%     1.35    1.83  0.45       1.50      1.03  0.39 1.12
    # 100%    2.53    1.88  5.10       2.97      2.19  9.67 3.46

sapply(cellClust.idx, function(x){quantile(sce.lc$sum[x])})
    #      ambig.lowNTx_A ambig.lowNTx_B ambig.lowNTx_C ambig.lowNTx_D ambig.lowNTx_E
    # 0%              239         240.00         578.00         692.00            823
    # 25%             935        1034.75        1897.50        1648.75           2173
    # 50%            2684        1763.00        3918.00        2232.50           3423
    # 75%           12478        2971.50        6972.75        3540.75           7090
    # 100%         146018       97039.00       36528.00       70563.00         150106

    #         Astro Endo.Mural Excit_A Excit_B Excit_C  Excit_D Excit_E  Excit_F
    # 0%     461.00       2772    4936    6596    9426   4765.0  8771.0   1353.0
    # 25%   1493.50       5793   19536   18892   21411  22112.5 18266.5   4824.0
    # 50%   2468.50       9373   30326   26572   34727  33001.0 27073.5  17998.5
    # 75%   4307.75      12830   42630   40847   52002  50180.5 38486.0  30538.0
    # 100% 27657.00      43741  106030  142023  182780 133387.0 73021.0 120749.0

    #       Inhib_A   Inhib_B Inhib_C Micro Neuron.5HT Neuron.NE Oligo      OPC
    # 0%      647.0  14333.00  4450.0   465    14002.0  15101.00  3431   627.00
    # 25%    3279.5  37479.00 11127.5  1139    43467.5  24208.75  5510  2645.50
    # 50%   10125.0  47500.50 16461.0  1882    56183.0  32792.50  7768  4549.00
    # 75%   23957.5  60213.75 22221.0  3154    66372.0  53717.50 11893  9329.25
    # 100% 167360.0 101398.00 74217.0 37667   134968.0 161314.00 50182 46062.00


## sizeFactor distribution
sapply(cellClust.idx, function(x){round(quantile(sce.lc$sizeFactor[x]), 2)})
    #     ambig.lowNTx_A ambig.lowNTx_B ambig.lowNTx_C ambig.lowNTx_D ambig.lowNTx_E
    # 0%             0.03           0.03           0.07           0.08           0.10
    # 25%            0.11           0.12           0.22           0.19           0.25
    # 50%            0.32           0.21           0.46           0.26           0.41
    # 75%            1.47           0.35           0.81           0.41           0.84
    # 100%          17.30          11.26           4.33           8.12          17.41

    #      Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A
    # 0%    0.05       0.33    0.57    0.78    1.09    0.56    1.02    0.16    0.08
    # 25%   0.17       0.69    2.29    2.22    2.48    2.59    2.16    0.57    0.39
    # 50%   0.28       1.11    3.51    3.15    4.03    3.85    3.17    2.13    1.20
    # 75%   0.50       1.52    4.90    4.79    6.03    5.87    4.53    3.56    2.78
    # 100%  3.21       5.18   12.20   16.47   21.66   15.81    8.65   14.31   19.83

    #      Inhib_B Inhib_C Micro Neuron.5HT Neuron.NE Oligo  OPC
    # 0%      1.65    0.52  0.06       1.62      1.79  0.40 0.07
    # 25%     4.31    1.28  0.13       5.02      2.87  0.64 0.31
    # 50%     5.46    1.90  0.22       6.46      3.87  0.90 0.53
    # 75%     6.93    2.57  0.36       7.67      6.27  1.38 1.11
    # 100%   11.66    8.61  4.46      15.53     18.71  5.82 5.34


# First remove 'ambig.lowNTx' ('low N transcripts')-driven or -associated clusters
#     (this is a lot of nuclei btw...)
sce.lc <- sce.lc[ ,-grep("ambig.lowNTx_", sce.lc$cellType.merged)]
sce.lc$cellType.merged <- droplevels(sce.lc$cellType.merged)

# Remove 0 genes across all nuclei
sce.lc <- sce.lc[!rowSums(assay(sce.lc, "counts"))==0, ]  #33353, but after dropping, a 32241

## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.lc

assay(sce.lc, "logcounts") <- NULL
sizeFactors(sce.lc) <- NULL
sce.lc <- logNormCounts(sce.lc)


    ## What if first re-assign 'inhib_g' & 'neuron.ambig_h' to their own? say, 'Inhib_D'
     # Motivation: that Inhib_A has NO top PW markers... 
    sce.lc$cellType.merged <- as.character(sce.lc$cellType.merged)
    sce.lc$cellType.merged[sce.lc$cellType %in% c("inhib_g","neuron.ambig_h")] <- "Inhib_D"
    sce.lc$cellType.merged <- factor(sce.lc$cellType.merged)
        
        #   Astro.TRUE Endo.Mural.TRUE    Excit_A.TRUE    Excit_B.TRUE    Excit_C.TRUE 
        #           74             387              29              22             133 
        # Excit_D.TRUE    Excit_E.TRUE    Excit_F.TRUE    Inhib_A.TRUE    Inhib_B.TRUE 
        #           11              55               8               1             171 
        # Inhib_C.TRUE    Inhib_D.TRUE      Micro.TRUE Neuron.5HT.TRUE  Neuron.NE.TRUE 
        #          183               2              59             101             104 
        #   Oligo.TRUE        OPC.TRUE 
        #          491              68 
    
        # [a bunch of other iterations in between...]
    
        ## LAST ATTEMPT: separate both 'inhib_d' & 'inhib_f'; keep 'inhib_g' and drop 'neuron.ambig_h'
        sce.lc$cellType.merged <- as.character(sce.lc$cellType.merged)
        sce.lc$cellType.merged[sce.lc$cellType =="inhib_d"] <- "Inhib_D"
        sce.lc$cellType.merged[sce.lc$cellType =="inhib_f"] <- "Inhib_E"
            sce.lc$cellType.merged[sce.lc$cellType =="inhib_g"] <- "Inhib_F"
            sce.lc$cellType.merged[sce.lc$cellType =="neuron.ambig_h"] <- "drop.lowNTx.neu"
            sce.lc <- sce.lc[ ,-grep("drop.", sce.lc$cellType.merged)]
        sce.lc$cellType.merged <- factor(sce.lc$cellType.merged)
            #			Astro.TRUE Endo.Mural.TRUE    Excit_A.TRUE    Excit_B.TRUE    Excit_C.TRUE 
            #             69             378              25              19             105 
            #   Excit_D.TRUE    Excit_E.TRUE    Excit_F.TRUE    Inhib_A.TRUE    Inhib_B.TRUE 
            #              9              51               7               5             139 
            #   Inhib_C.TRUE    Inhib_D.TRUE    Inhib_E.TRUE    Inhib_F.TRUE      Micro.TRUE 
            #            168               6               5              22              58 
            #Neuron.5HT.TRUE  Neuron.NE.TRUE      Oligo.TRUE        OPC.TRUE 
            #             96             101             489              64 
    
            ##  THIS LOOKS THE BEST --> proceed to make those changes in '04_clusterAgglomeration.R',
             #                          then re-run through this script
        
    
### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellClust.idx <- splitit(sce.lc$cellType.merged)
medianNon0.lc.16 <- lapply(cellClust.idx, function(x){
  apply(as.matrix(assay(sce.lc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.lc.16, table)
    #       Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A
    # FALSE 31841      30158   25733   26424   25886   26037   25987   28280   29485
    # TRUE    400       2083    6508    5817    6355    6204    6254    3961    2756
    
    #       Inhib_B Inhib_C Micro Neuron.5HT Neuron.NE Oligo   OPC
    # FALSE   23514   28229 31987      23828     26325 30201 31124
    # TRUE     8727    4012   254       8413      5916  2040  1117




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
                                medianNon0.lc.16[[i]][match(rownames(markers.lc.t.pw[[i]]),
                                                         names(medianNon0.lc.16[[i]]))])
  colnames(markers.lc.t.pw[[i]])[21] <- "non0median"
}

sapply(markers.lc.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    #   Astro.TRUE Endo.Mural.TRUE    Excit_A.TRUE    Excit_B.TRUE    Excit_C.TRUE 
    #           74             389              29              23             138 
    # Excit_D.TRUE    Excit_E.TRUE    Excit_F.TRUE      Inhib_A.NA    Inhib_B.TRUE 
    #           16              55               8              NA             174 
    # Inhib_C.TRUE      Micro.TRUE Neuron.5HT.TRUE  Neuron.NE.TRUE      Oligo.TRUE 
    #          184              59             101             105             495 
    #     OPC.TRUE 
    #           68 


## Save these
save(markers.lc.t.pw, medianNon0.lc.16,
     file=here("processed_data","SCE",
               "markers-stats_LC-n3_findMarkers_16cellTypes.rda"))


    # As needed
    load(here("processed_data","SCE",
              "markers-stats_LC-n3_findMarkers_16cellTypes.rda"), verbose=T)


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
      scale_color_manual(values = tableau20) +  
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
      scale_color_manual(values = tableau20) +  
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
                                          medianNon0.lc.16[[i]][match(rownames(markers.lc.t.1vAll[[i]][["0"]]),
                                                                   names(medianNon0.lc.16[[i]]))])
  colnames(markers.lc.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.lc.t.1vAll[[i]][["1"]] <- cbind(markers.lc.t.1vAll[[i]][["1"]],
                                          medianNon0.lc.16[[i]][match(rownames(markers.lc.t.1vAll[[i]][["1"]]),
                                                                   names(medianNon0.lc.16[[i]]))])
  colnames(markers.lc.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.lc.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}


## Save these
save(markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.16,
     file=here("processed_data","SCE",
               "markers-stats_LC-n3_findMarkers_16cellTypes.rda"))



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
      scale_color_manual(values = tableau20) +  
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
write.csv(top40genes, file=here("code","analyses_sn","top40genesLists_LC-n3_16cellTypes.csv"),
          row.names=FALSE)




### Some follow-up ================================================
load(here("processed_data","SCE",
          "markers-stats_LC-n3_findMarkers_16cellTypes.rda"), verbose=T)

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
    scale_color_manual(values = tableau20) +  
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
        scale_color_manual(values = tableau20) +  
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
          "markers-stats_LC-n3_findMarkers_16cellTypes.rda"), verbose=T)
    # markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.16

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
                                medianNon0.lc.16[[i]][match(rownames(markers.lc.t.pw.some[[i]]),
                                                            names(medianNon0.lc.16[[i]]))])
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

            
            
save(markers.lc.t.pw.some, medianNon0.lc.16,
     file=here("processed_data","SCE",
          "markers-stats_LC-n3_findMarkers-PW-some_16cellTypes.rda"))

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
      scale_color_manual(values = tableau20) +  
      ggtitle(label=paste0("LC ", i, " top markers: single-nucleus-level p.w.-SOME t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}






## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    #
proc.time()
    #     user    system   elapsed 
    #
options(width = 120)
session_info()
    #


