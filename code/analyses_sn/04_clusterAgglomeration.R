### LC snRNA-seq analysis
### Cluster agglomeration
###    Motivation is that we want a less-resolved atlas for purposes of
###    spatial spot-level deconvolution; something akin to the 'step 2' from
###    Tran-Maynard et al. Neuron 2021 (hierarchical clustering-based)
###     qrsh -l bluejay,mf=76G,h_vmem=80G
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

here()
# [1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"



## Load SCE
load(here("processed_data","SCE", "sce_updated_LC.rda"), verbose=T)
    # sce.lc, medianNon0.lc, hdgs.lc, annotationTab.lc

## Load graph & community detection that identified the 60 clusters
load(here("processed_data","SCE","graph_communities_glmpcamnn_LC.rda"), verbose=T)
    # snn.gr.glmpcamnn, clusters.glmpcamnn


## Assess cluster modularity (a measure of cluster separation) 
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



### Cluster agglomeration with `cut_at()` ================================
  # - This function relies on a user-defined k (number of clusters desired)







### Tran-Maynard et al. Neuron 2021 method ================================







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



