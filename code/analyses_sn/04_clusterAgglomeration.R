### LC snRNA-seq analysis
### Cluster agglomeration
###    Motivation is that we want a less-resolved atlas for purposes of
###    spatial spot-level deconvolution; something akin to the 'step 2' from
###    Tran-Maynard et al. Neuron 2021 (hierarchical clustering-based)
###     qrsh -l bluejay,mf=92G,h_vmem=96G
### Initiated: MNT 01Apr2022

library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(scry)
library(batchelor)
library(bluster)
library(pheatmap)
library(dynamicTreeCut)
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


## Load SCE
load(here("processed_data","SCE", "sce_updated_LC.rda"), verbose=T)
    # sce.lc, medianNon0.lc, hdgs.lc, annotationTab.lc

## Load graph & community detection that identified the 60 clusters
load(here("processed_data","SCE","graph_communities_glmpcamnn_LC.rda"), verbose=T)
    # snn.gr.glmpcamnn, clusters.glmpcamnn


## Assess pairwise cluster modularity (a measure of cluster separation) for
 #    the 60 graph-based clusters
mod.ratio <- pairwiseModularity(graph = snn.gr.glmpcamnn,
                                clusters = sce.lc$clusters.glmpcamnn,
                                as.ratio=TRUE)
colnames(mod.ratio) <- annotationTab.lc$cellType[match(colnames(mod.ratio),
                                                       annotationTab.lc$cluster)]
rownames(mod.ratio) <- annotationTab.lc$cellType[match(rownames(mod.ratio),
                                                       annotationTab.lc$cluster)]

# Plot
pdf(here("plots","snRNA-seq","clusterModularityRatio_LC-n3-60snnClusters.pdf"))
pheatmap(log2(mod.ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100),
         main="Modularity ratio for 60 SNN-clusters in LC (n=3)",
         fontsize_row=7, fontsize_col=7, angle_col=90)
grid::grid.text(label="log2(ratio)",x=0.96,y=0.65, gp=grid::gpar(fontsize=7))
dev.off()





### Tran-Maynard et al. Neuron 2021 method ================================

## Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.lc$cellType)
cluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.lc)$counts[ ,ii])
}
)

# Btw - total N transcripts / cluster:
sapply(clusIndexes, function(x){quantile(sce.lc$sum[x])}) # collapse chunk ====
    #      ambig.lowNTx Astro_A Astro_B Astro_C Astro_D Astro_E Endo.Mural_A Endo.Mural_B   Excit_A
    # 0%         239.00   461.0  1699.0    3365    5274    2223         2772       662.00   4765.00
    # 25%        447.75   923.5  2279.5    4030    6898    4021         4730      1176.25  19503.00
    # 50%        676.00  1311.0  2720.5    4508    7580    5570         8336      1582.00  29254.00
    # 75%       1120.00  1716.5  3039.5    5010    9217    8722        15046      2178.25  46196.25
    # 100%     25442.00 10876.0  4367.0   10198   13260   27657        43741      4919.00 133387.00
    
    #       Excit_B Excit_C   Excit_D  Excit_E  Excit_F Excit_G Excit_H Excit_I  Excit_J Excit_K Excit_L
    # 0%    4936.00    6596   6460.00   1116.0  11638.0    8808   32203   12266  29569.0  8771.0  9426.0
    # 25%  14316.75   15758  17894.25  10081.0  21738.0   28633   42111   27552  48452.5 18266.5 18836.0
    # 50%  23391.50   21552  24639.00  22227.0  30554.0   38305   53206   40118  58725.0 27073.5 23634.0
    # 75%  37051.25   29348  37863.50  37291.5  42491.5   49491   67837   56779  81222.0 38486.0 32343.5
    # 100% 90617.00   99198 120749.00 128727.0 142023.0  106030  182780  117807 133646.0 73021.0 52988.0
    
    #       Excit_M   Excit_N Inhib_A  Inhib_B Inhib_C   Inhib_D Inhib_E Inhib_F   Inhib_G Inhib_H
    # 0%    1353.00   1920.00  4450.0  14333.0    4071   4869.00   647.0    4047  13076.00   29538
    # 25%   2053.75   2742.75 11127.5  40492.5   11582  14700.75  1515.0    9357  33978.75   37949
    # 50%   3051.00   5752.00 16461.0  52043.0   19354  23123.00  2253.0   15461  43282.00   44892
    # 75%   4242.00  43216.25 22221.0  62016.0   30800  34026.25  3670.5   25945  76209.25   58259
    # 100% 21886.00 150106.00 74217.0 101398.0  109878 167360.00 10652.0   94071 124649.00   84401
    
    #      Inhib_I Macrophage Micro_A Micro_B Micro_C Neuron.5HT Neuron.5HT_noDDC Neuron.ambig_A
    # 0%   16811.0    5533.00    1959   810.0  465.00    14002.0          1375.00          276.0
    # 25%  28339.5    8276.75    3031  1438.0  815.75    43467.5          2971.25         1915.5
    # 50%  38143.5    9767.00    3527  1815.0  991.00    56183.0          4548.50         3053.5
    # 75%  46867.5   12005.00    4713  2089.5 1606.75    66372.0          8243.75         5355.5
    # 100% 94607.0   18780.00   24246 37667.0 4255.00   134968.0         38921.00       107202.0
    
    #      Neuron.ambig_B Neuron.ambig_C Neuron.ambig_D Neuron.ambig_E Neuron.ambig_F Neuron.ambig_G
    # 0%              692            886          578.0        1023.00           1813          823.0
    # 25%            1616           2475         1422.0        2167.75           3920         2042.5
    # 50%            1978           3919         2220.0        3334.50           6219         3180.0
    # 75%            2824           6221         4688.5        5444.25           8291         4761.5
    # 100%          23414          55227        19594.0       79135.00          36528        16659.0
    
    #      Neuron.ambig_H Neuron.mixed_A Neuron.mixed_B Neuron.NE Oligo_A Oligo_B Oligo_C Oligo_D
    # 0%           1592.0            885           2689  15101.00   892.0    3431   15554    8173
    # 25%          2680.5          13711           6744  24208.75  1342.5    4438   21327    9679
    # 50%          4455.0          26142          12224  32792.50  1651.5    5192   24410   11918
    # 75%          7578.5          41504          23154  53717.50  2211.5    6372   31306   13646
    # 100%        20806.0         146018          70563 161314.00 30386.0   37604   50182   24920
    
    #       Oligo_E Oligo_F Oligo_G Oligo_H Oligo_I   OPC_A   OPC_B OPC_C
    # 0%     240.00    1883    4148    1967    7628   627.0  7125.0  3744
    # 25%    759.25    2498    5579    2800    8944  1612.5 10102.0  5012
    # 50%   1195.50    2926    6523    3371   10432  2651.0 13224.5  5828
    # 75%   2475.50    3327    7759    4165   12197  3559.5 17257.5  6858
    # 100% 49835.00   43087   13975   97039   17725 46062.0 28898.0 21933
# end collapse chunk ====

# Compute LSFs at this level
sizeFactors.PB  <- librarySizeFactors(cluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(cluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB + 1)}))


## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# Just for observation
myplclust(tree.clusCollapsed, main="LC (n=3) SNN-cluster relationships (60 clusters)",
          cex.main=1, cex.lab=0.8, cex=0.6)

clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=210)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
    ## Cutting at 210 looks good for the main neuronal branch, but a lot of glial
    #    prelim clusters are dropped off (0's)

    # Cut at 425 for broad glia branch (will manually merge remaining dropped off)    
    glia.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                                  minClusterSize=2, deepSplit=1, cutHeight=425)
    unname(glia.treeCut[order.dendrogram(dend)])
    
    # Take those and re-assign to the first assignments
    clust.treeCut[order.dendrogram(dend)][c(38:60)] <- ifelse(glia.treeCut[order.dendrogram(dend)][c(38:60)] == 0,
                                                              0, glia.treeCut[order.dendrogram(dend)][c(38:60)] + max(clust.treeCut))

    unname(clust.treeCut[order.dendrogram(dend)])

# Add new labels to those prelimClusters cut off
clustersStill0 <- which(clust.treeCut[order.dendrogram(dend)]==0)
clust.treeCut[order.dendrogram(dend)][clustersStill0] <- max(clust.treeCut) +
  #c(1:length(clustersStill0))
  c(1,1, 2:11, 12,12, 13,13)
  # 'repeating' bc some of these branches are just pairs/would be merged at a higher cut height

# 'Re-write', since there are missing numbers
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(c(tableau20, tableau10medium)[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
dendextend::labels_colors(dend) <- cluster_colors[as.character(clust.treeCut[order.dendrogram(dend)])]

# Print for future reference
pdf(here("plots","snRNA-seq","hierarchicalClustering_LC-n3_SNNcluster-relationships.pdf"))
par(cex=0.7, font=2, mar=c(6,4,4,2))
plot(dend, cex.main=2,
     main="LC (n=3) prelim-SNN-cluster relationships \nby hierarchical clustering")
dev.off()


# Make reference for new cluster assignment
clusterRefTab.lc <- data.frame(graphClust = labels(dend),
                               merged = clust.treeCut[order.dendrogram(dend)])
# Add this to the annotationTab to this
annotationTab.lc$merged <- clusterRefTab.lc$merged[match(annotationTab.lc$cellType,
                                                         clusterRefTab.lc$graphClust)]


# Assign as 'mergedCluster.HC'
sce.lc$mergedCluster.HC <- factor(annotationTab.lc$merged[match(sce.lc$cellType,
                                                             annotationTab.lc$cellType)])

# Save
save(sce.lc, annotationTab.lc, medianNon0.lc, hdgs.lc,
     file=here("processed_data","SCE", "sce_updated_LC.rda"))


## Cluster modularity?
## Assess cluster modularity (a measure of cluster separation) 
mod.ratio.merged.HC <- pairwiseModularity(graph = snn.gr.glmpcamnn,
                                          clusters = sce.lc$mergedCluster.HC,
                                          as.ratio=TRUE)

# Plot
pdf(here("plots","snRNA-seq","clusterModularityRatio_LC-n3_26collapsedClusters-by-HC.pdf"))
# Add UMAP to this
print(
  plotReducedDim(sce.lc, dimred="UMAP", colour_by="mergedCluster.HC", text_by="mergedCluster.HC") +
    ggtitle(paste0("UMAP with 26 HC-merged clusters in LC (n=3)")) +
    scale_color_manual(values = c(tableau20, tableau10medium)) + labs(colour="New Cluster")
)
# Heatmap
pheatmap(log2(mod.ratio.merged.HC+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100),
         main="Modularity ratio for 26 HC-merged clusters in LC (n=3)",
         fontsize_row=7, fontsize_col=7, angle_col=90,
         display_numbers=T, number_format="%.1f", na_col="darkgrey")
grid::grid.text(label="log2(ratio)",x=0.97,y=0.64, gp=grid::gpar(fontsize=7))
dev.off()
    # Observations: overall, honestly this isn't too bad.
    #   - 18 and 7 share edge weights - this makes sense based on their HC
    #   - 16 == Neuron.5HT  &  19 == Neuron.5HT_noDDC
    #   - 18 <~ 24 ~> 20      Neuron.ambig_D <~ Neuron.ambig_E/_F ~> Neuron.NE
    #     (not really showing shared 'lineages' in the HC)


table(sce.lc$cellType, sce.lc$mergedCluster.HC)

# By eye and briefly looking at markers, these mergings look pretty good
#     other than:
table(droplevels(sce.lc$cellType[which(sce.lc$mergedCluster.HC == 1)]))
    # Excit_D        Excit_E        Inhib_G Neuron.ambig_A Neuron.mixed_A 
    #     114            219             50           2560           2301
    #     (thus the self-modularity ratio is pretty low)


### Since the below method performs poorly, print broad markers for annotation =====
  #       of these 26 merged clusters
source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")

markers.mathys.tran = list(
  'neuron' = c('SYT1', 'SNAP25', 'GRIN1'),
  'excit_neuron' = c('CAMK2A', 'NRGN','SLC17A7', 'SLC17A6', 'SLC17A8'),
  'inhib_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
  # Norepinephrine & serotonergic markers
  'neuron.NE' = c("TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC"), #SLC6A3 - saw no DAT
  'neuron.5HT' = c("SLC6A4", "TPH1", "TPH2", "DDC"),
  # SERT, serotonin T (aka 5-HTT); 
  'monoamine.metab' = c("COMT", "MAOA", "MAOB"),
  # MSN markers
  'MSNs.pan' = c("PPP1R1B","BCL11B"),# "CTIP2")
  'MSNs.D1' = c("DRD1", "PDYN", "TAC1"),
  'MSNs.D2' = c("DRD2", "PENK"),
  ## Non-neuronal:
  'oligodendrocyte' = c('MBP', 'MOBP', 'PLP1'),
  'oligo_precursor' = c('PDGFRA', 'VCAN', 'CSPG4'),
  'microglia' = c('CD74', 'CSF1R', 'C3'),
  'astrocyte' = c('GFAP', 'TNC', 'AQP4', 'SLC1A2'),
  'endothelial' = c('CLDN5', 'FLT1', 'VTN'),
  # Post-hoc from Tran-Maynard, et al. Neuron 2021
  'differn_committed_OPC' = c("SOX4", "BCAN", "GPR17", "TNS3"),
  'Tcell' = c('SKAP1', 'ITK', 'CD247'),
  'Mural' = c('COL1A2', 'TBX18', 'RBPMS'),
  'Macro' = c('CD163', 'SIGLEC1', 'F13A1')
)


pdf(here("plots","snRNA-seq",paste0("LC-n3_expression-violin_markers_26-HC-merged-graphClusters.pdf")),
    height=6, width=12)
for(i in 1:length(markers.mathys.tran)){
  print(
    plotExpressionCustom(sce = sce.lc,
                         exprs_values = "logcounts",
                         features = markers.mathys.tran[[i]], 
                         features_name = names(markers.mathys.tran)[[i]], 
                         #anno_name = "mergedCluster.HC",
                         anno_name = "cellType.merged",
                         ncol=2, point_alpha=0.4, point_size=0.9,
                         scales="free_y", swap_rownames="gene_name") +  
      ggtitle(label=paste0("LC-n3 HC-merged (26) clusters: ",
                           names(markers.mathys.tran)[[i]], " markers")) +
      theme(plot.title = element_text(size = 12),
            axis.text.x = element_text(size=7)) +
      scale_color_manual(values = c(tableau20, tableau10medium))
  )
}
dev.off()


# Annotation reassignment ===
annotationTab.lc <- annotationTab.lc[order(annotationTab.lc$merged), ]


# Reference (smaller table)
annTab.lc <- data.frame(cluster=c(1:26))
annTab.lc$cellType.merged <- NA

annTab.lc$cellType.merged[c(2,3,6,7, 15,23)] <- paste0("Excit_", c("A","B","C","D", "E","F"))
annTab.lc$cellType.merged[c(4,5,16,24)] <- paste0("Inhib_", c("A","B","C","D"))
annTab.lc$cellType.merged[c(1,18)] <- paste0("Neuron.mixed_", c("A","B"))
annTab.lc$cellType.merged[c(19,22,25,26)] <- paste0("Neuron.ambig_", c("A","B","C","D"))
annTab.lc$cellType.merged[c(14)] <- "ambig.lowNTx"

annTab.lc$cellType.merged[c(21)] <- "Neuron.NE"
annTab.lc$cellType.merged[c(17)] <- "Neuron.5HT"
annTab.lc$cellType.merged[c(20)] <- "Neuron.5HT_noDDC"

annTab.lc$cellType.merged[c(8,10)] <- paste0("Oligo_", c("A","B"))
annTab.lc$cellType.merged[c(12)] <- "OPC"
annTab.lc$cellType.merged[c(11)] <- "Micro"
annTab.lc$cellType.merged[c(9)] <- "Astro"
annTab.lc$cellType.merged[c(13)] <- "Endo.Mural"

# Add to the bigger reference
annotationTab.lc$cellType.merged <- annTab.lc$cellType.merged[match(annotationTab.lc$merged,
                                                                    annTab.lc$cluster)]

# Add to SCE
sce.lc$cellType.merged <- annotationTab.lc$cellType.merged[match(sce.lc$mergedCluster.HC,
                                                                 annotationTab.lc$merged)]
sce.lc$cellType.merged <- factor(sce.lc$cellType.merged)


# To distinguish original 60 clusters' annotations from these ones, `tolower()` the former
annotationTab.lc$cellType <- tolower(annotationTab.lc$cellType)
sce.lc$cellType <- factor(tolower(sce.lc$cellType))

### Save
save(sce.lc, annotationTab.lc, medianNon0.lc, hdgs.lc,
     file=here("processed_data","SCE", "sce_updated_LC.rda"))

    # --> re-print the marker plots, above, with these annotations





### Cluster agglomeration with `cut_at()` ================================
  # - This function relies on a user-defined k (number of clusters desired)

sce.test <- sce.lc

pdf(here("plots","snRNA-seq","explore_iterative_LC-n3_collapsedClusters-with-cut_at.pdf"))
for(i in seq(10,30, by=5)){
  colData(sce.test)[paste0("merged.cut_at", i)] <- factor(igraph::cut_at(clusters.glmpcamnn, n=i))
  table(colData(sce.test)[paste0("merged.cut_at", i)])
  
  # Plot on UMAP
  print(
    plotReducedDim(sce.test, dimred="UMAP", colour_by=paste0("merged.cut_at", i),
                   text_by=paste0("merged.cut_at", i)) +
      scale_color_manual(values = c(tableau20, tableau10medium)) +
      ggtitle(paste0("UMAP with ", i, " cut_at-merged clusters in LC (n=3)")) +
      labs(colour="New Cluster")
  )
  
  ## Compute & plot pw modularity ratio
  mod.ratio.temp  <- pairwiseModularity(graph = snn.gr.glmpcamnn,
                                        clusters = colData(sce.test)[paste0("merged.cut_at", i)][ ,1],
                                        as.ratio=TRUE)
  
  # Plot heatmap
  pheatmap(log2(mod.ratio.temp+1), cluster_rows=FALSE, cluster_cols=FALSE,
           color=colorRampPalette(c("white", "blue"))(100),
           main=paste0("Modularity ratio for ", i, " cut_at-merged clusters in LC (n=3)"),
           fontsize_row=7, fontsize_col=7, angle_col=90,
           display_numbers=T, number_format="%.1f", na_col="darkgrey")
  grid::grid.text(label="log2(ratio)",x=0.97,y=0.64, gp=grid::gpar(fontsize=7))
}
dev.off()

    ## Observations ===
    #    - even cluster 3 (with cut_at(...,n=30) really swallows everything,
    #      pretty indiscriminately...
    levels(droplevels(sce.test$cellType[which(sce.test$merged.cut_at30 == 3)]))
        # [1] "ambig.lowNTx"     "Endo.Mural_B"     "Excit_E"          "Excit_M"         
        # [5] "Excit_N"          "Inhib_E"          "Inhib_F"          "Inhib_G"         
        # [9] "Inhib_I"          "Micro_C"          "Neuron.5HT_noDDC" "Neuron.ambig_A"  
        # [13] "Neuron.ambig_B"   "Neuron.ambig_C"   "Neuron.ambig_E"   "Neuron.ambig_G"  
        # [17] "Neuron.ambig_H"   "Neuron.mixed_A"   "Oligo_A"          "Oligo_E"
    
    #    - Unsurprisingly cluster 3 gets a really low self score, but otherwise
    #      it's kind of hard to figure out what looks 'good', b/tw this and the HC method,
    #      just based on this PW modularity ratio heatmap

    #     --> not saving into the colData bc this attempt is just no good...


    
    
### Simplification just from previous 60 clusters-level annotation ======
sce.lc$cellType.collapsed <- as.factor(ss(as.character(sce.lc$cellType), "_", 1))
# Re-organize
colData(sce.lc) <- colData(sce.lc)[ ,c(1:13, 15, 14)]

# And then just merge microglia & macrophages together, since they're so similar
sce.lc$cellType.collapsed <- ifelse(sce.lc$cellType.collapsed %in% c("Micro","Macrophage"),
                                    "Micro.Macro", as.character(sce.lc$cellType.collapsed))

# Oh shoot, we lost the 'noDDC' set of 5-HT-like neurons:
sce.lc$cellType.collapsed[sce.lc$cellType == "Neuron.5HT_noDDC"] <- "Neuron.5HT_noDDC"

sce.lc$cellType.collapsed <- factor(sce.lc$cellType.collapsed)


# Save
save(sce.lc, annotationTab.lc, medianNon0.lc, hdgs.lc,
     file=here("processed_data","SCE", "sce_updated_LC.rda"))


## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    # [1] "2022-04-07 22:57:06 EDT"
proc.time()
    #     user    system   elapsed 
    # 1524.301    38.689 11076.920 
options(width = 120)
session_info()
    # ─ Session info ────────────────────────────────────────────────────────────────────────
    # setting  value
    # version  R version 4.1.2 Patched (2021-11-04 r81138)
    # os       CentOS Linux 7 (Core)
    # system   x86_64, linux-gnu
    # ui       X11
    # language (EN)
    # collate  en_US.UTF-8
    # ctype    en_US.UTF-8
    # tz       US/Eastern
    # date     2022-04-07
    # pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
    # 
    # ─ Packages ────────────────────────────────────────────────────────────────────────────
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
    # cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
    # cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
    # colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
    # cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
    # crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
    # DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
    # DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
    # DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
    # dendextend             1.15.2   2021-10-28 [2] CRAN (R 4.1.2)
    # digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
    # dplyr                  1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
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
    # ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
    # ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
    # glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
    # googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
    # gridExtra            * 2.3      2017-09-09 [2] CRAN (R 4.1.0)
    # gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
    # HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
    # here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
    # igraph                 1.3.0    2022-04-01 [2] CRAN (R 4.1.2)
    # IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
    # irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
    # jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
    # labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
    # lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
    # lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
    # limma                  3.50.1   2022-02-17 [2] Bioconductor
    # locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
    # magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
    # Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
    # MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
    # matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
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
    # tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
    # utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
    # vctrs                  0.4.0    2022-03-30 [2] CRAN (R 4.1.2)
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
    # ───────────────────────────────────────────────────────────────────────────────────────


