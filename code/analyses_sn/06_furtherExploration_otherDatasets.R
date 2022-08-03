### LC snRNA-seq analysis
### Exploration with other datasets
###     qsub -l bluejay,mf=76G,h_vmem=80G
### Initiated: MNT 06May2022

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
library(pheatmap)


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

load(here("processed_data","SCE", "sce_updated_LC.rda"), verbose=T)
# sce.lc, annotationTab.lc, medianNon0.lc, hdgs.lc, cell_colors.lc


## Exploring mouse markers of interest ===========================
# Provided by Lukas
markers.Mulvey.human <- c(
  "AGTR1", "ASB4", "CALCR", "CALR3", "CHODL", "CHRNA6", "CILP", "CYB561", 
  "DBH", "DDC", "DLK1", "ADGRE1", "EYA2", "SHISAL2B", "FAM183A", "FIBCD1", 
  "GAL", "GCH1", "GLRA2", "GNG4", "GPX3", "GTF2A1L", "HCRTR1", "IGSF5", "MAOA", 
  "MRAP2", "MYOM2", "NEUROG2", "SLC9B2", "NXPH4", "OVGP1", "PCBD1", "PHOX2A", 
  "PHOX2B", "PLA2G4D", "PTGER2", "SLC18A2", "SLC31A1", "SLC6A2", "STBD1", 
  "SYT17", "TH", "TM4SF1", "TM4SF5", "TRAF3IP2")

markers.Grimm.human <- c(
  "CBR3", "DNAH5", "SERPINE1", "LAYN", "TPH2", "RPH3AL", "NGB", "CYB561", 
  "GNAS", "SLC31A1", "TCP1", "PPIC", "COLEC10", "RAB3B", "MAOA", "PCBP3", 
  "TSPAN12", "FBP1", "DBH", "SERPINF1", "TXK", "SEC16B", "TRAF1", "PTGES", 
  "GGT5", "MMP2", "MCAM", "TFAP2A", "ACSL1", "UPB1", "UCP3", "COL5A1", 
  "ALDH1A1", "FBN1")
markers.Grimm.human <- sort(markers.Grimm.human)

table(markers.Mulvey.human %in% rowData(sce.lc)$gene_name)
table(markers.Grimm.human %in% rowData(sce.lc)$gene_name)
intersect(markers.Mulvey.human, markers.Grimm.human)
# [1] "CYB561"  "DBH"     "MAOA"    "SLC31A1"

markers2print <- list(Mulvey.etal = markers.Mulvey.human,
                      Grimm.etal = markers.Grimm.human)


# First remove 'ambig.lowNTx' ('low N transcripts')-driven or -associated clusters
#     (this is a lot of nuclei btw...)
sce.lc <- sce.lc[ ,-grep("ambig.lowNTx_", sce.lc$cellType.merged)]
sce.lc$cellType.merged <- droplevels(sce.lc$cellType.merged)

# Load marker stats so can order by enrichment statistics
load(here("processed_data","SCE",
          "markers-stats_LC-n3_findMarkers_19cellTypes.rda"), verbose=T)
    # markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.19


## Heatmap ===
for(i in names(markers2print)){
  pdf(here("plots","snRNA-seq","exploration",
           paste0("heatmap_", i, "Markers_snRNA-seq_by19MergedClusters.pdf")), height=6, width=10)
  ## Means version:
  current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[markers2print[[i]], ii])))
  # Set neuronal pops first
  neuronPosition <- c(grep("Excit", colnames(current_dat)),
                      grep("Inhib", colnames(current_dat)),
                      grep("Neuron", colnames(current_dat)))
  reorderedCols <- c(neuronPosition, setdiff(1:19, neuronPosition))
  current_dat <- current_dat[ ,reorderedCols]
  # Put NE neurons before the 5-HT ones
  current_dat <- current_dat[ ,c(1:12, 14, 13, 15:19)]
  
  markers.NE.enriched <- markers.lc.t.1vAll[["Neuron.NE"]][["Neuron.NE_enriched"]]
  rownames(markers.NE.enriched) <- rowData(sce.lc)$gene_name[match(rownames(markers.NE.enriched), rowData(sce.lc)$gene_id)]
  
  # Because of the similarities/confounding b/tw NE & 5-HT neurons:
      markers.5HT.enriched <- markers.lc.t.1vAll[["Neuron.5HT"]][["Neuron.5HT_enriched"]]
      rownames(markers.5HT.enriched) <- rowData(sce.lc)$gene_name[match(rownames(markers.5HT.enriched), rowData(sce.lc)$gene_id)]
      # Re-order to same as above
      markers.5HT.enriched <- markers.5HT.enriched[rownames(markers.NE.enriched), ]
      # Average the logFC
      markers.NE.enriched$std.logFC.alt <- (markers.NE.enriched$std.logFC + markers.5HT.enriched$std.logFC) / 2
      
  genes.ordered <- order(markers.NE.enriched[markers2print[[i]], ]$std.logFC.alt, decreasing=T)
  current_dat <- current_dat[genes.ordered, ]
  
  italicnames <- lapply(
    rownames(current_dat),
    function(x) bquote(italic(.(x))))
  
  # Print
  pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE,
           breaks = seq(0.02, 4, length.out = 101),
           color = viridis::viridis(100),
           main=paste0("\t\t\t",i," orthologous markers, 19 LC cell classes (means)"),
           labels_col = as.expression(italicnames),
           fontsize=12, fontsize_row = 15, fontsize_col=14)
  grid::grid.text(label="log2-\nExprs", x=0.97, y=0.57, gp=grid::gpar(fontsize=10))
  
  ## or medians version:
  current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMedians(dat[markers2print[[i]], ii])))
  # For some reason rownames aren't kept:
  rownames(current_dat) <- markers2print[[i]]
  # Set neuronal pops first
  neuronPosition <- c(grep("Excit", colnames(current_dat)),
                      grep("Inhib", colnames(current_dat)),
                      grep("Neuron", colnames(current_dat)))
  reorderedCols <- c(neuronPosition, setdiff(1:19, neuronPosition))
  current_dat <- current_dat[ ,reorderedCols]
  # Put NE neurons before the 5-HT ones
  current_dat <- current_dat[ ,c(1:12, 14, 13, 15:19)]
  current_dat <- current_dat[genes.ordered, ]
  
  # Print
  pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE,
           breaks = seq(0.02, 4, length.out = 101),
           color = viridis::viridis(100),
           main=paste0("\t\t\t\t",i," orthologous markers, 19 LC cell classes (medians)"),
           labels_col = as.expression(italicnames),
           fontsize=12, fontsize_row = 15, fontsize_col=14)
  grid::grid.text(label="log2-\nExprs", x=0.97, y=0.57, gp=grid::gpar(fontsize=10))
  
  dev.off()
}









### Broad markers heatmap ===
genes <- c('SNAP25','SLC17A7','SLC17A6','GAD1','GAD2',
           # NE neuron markers
           "TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC",
           # serotonergic markers (includes DDC but repetitive)
           "SLC6A4", "TPH2",  # (TPH1 not expressed by these clusters)
           
           ## Non-neuronal:
           # Astro
           'AQP4','GFAP',
           # Endo, Mural (RBPMS)
           'CLDN5','FLT1','RBPMS',
           # Macrophage, Microglia
           'CD163','C3',
           # Oligo
           'MBP',
           # OPC
           'PDGFRA','VCAN')

cell.idx <- splitit(sce.lc$cellType.merged)
dat <- as.matrix(assay(sce.lc, "logcounts"))
rownames(dat) <- rowData(sce.lc)$gene_name



pdf(here("plots","snRNA-seq","heatmap_broadMarkers_by19MergedClusters.pdf"), height=7, width=7.5)
par(mar=c(5,8,4,2))
## Means version:
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[genes, ii])))
# Set neuronal pops first
neuronPosition <- c(grep("Excit", colnames(current_dat)),
                    grep("Inhib", colnames(current_dat)),
                    grep("Neuron", colnames(current_dat)))
reorderedCols <- c(neuronPosition, setdiff(1:19, neuronPosition))
current_dat <- current_dat[ ,reorderedCols]
# Put NE neurons before the 5-HT ones
current_dat <- current_dat[ ,c(1:12, 14, 13, 15:19)]

italicnames <- lapply(
  rownames(current_dat),
  function(x) bquote(italic(.(x))))

# Print
pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         main="\t\t\t19 LC cell class broad marker expression profiles (means)",
         labels_col = as.expression(italicnames),
         fontsize=12, fontsize_row = 15, fontsize_col=14)
grid::grid.text(label="log2-\nExprs", x=0.96, y=0.63, gp=grid::gpar(fontsize=10))

## or medians version:
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMedians(dat[genes, ii])))
# For some reason rownames aren't kept:
rownames(current_dat) <- genes
# Set neuronal pops first
neuronPosition <- c(grep("Excit", colnames(current_dat)),
                    grep("Inhib", colnames(current_dat)),
                    grep("Neuron", colnames(current_dat)))
reorderedCols <- c(neuronPosition, setdiff(1:19, neuronPosition))
current_dat <- current_dat[ ,reorderedCols]
# Put NE neurons before the 5-HT ones
current_dat <- current_dat[ ,c(1:12, 14, 13, 15:19)]
# Print
pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE,
         breaks = seq(0.02, 4, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         main="\t\t\t\t19 LC cell class broad marker expression profiles (medians)",
         labels_col = as.expression(italicnames),
         fontsize=12, fontsize_row = 15, fontsize_col=14)
grid::grid.text(label="log2-\nExprs", x=0.96, y=0.63, gp=grid::gpar(fontsize=10))
dev.off()






### Print Lukas' pseudo-bulked (annotated 'LC vs WM') DE markers ===============
  # (these are all already subsetted for FDR < 0.05)
visiumDE_LCvWM <- read.csv(here("outputs","06_downstream","pseudobulkDE","pseudobulk_LCvsWM",
                                "LC_pseudobulkDE_LCvsWM_topGenes.csv"), sep=",")
    # 51 DE genes

visiumDE_withinLC <- read.csv(here("outputs","06_downstream","pseudobulkDE","pseudobulk_withinLCregions",
                                "LC_pseudobulkDE_withinLCregions_topGenes.csv"), sep=",")
    # 55 DE genes

list.visiumDE <- list(LCvsWM = visiumDE_LCvWM$gene_name,
                      withinAnnotatedLC = visiumDE_withinLC$gene_name)


## Print heatmaps, as above ===
for(i in names(list.visiumDE)){
  pdf(here("plots","snRNA-seq","exploration",
           paste0("heatmap_DEgenes_visium-", i, "_snRNA-seq_by19MergedClusters.pdf")), height=6, width=12)
  ## Means version:
  current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[list.visiumDE[[i]], ii])))
  # Set neuronal pops first
  neuronPosition <- c(grep("Excit", colnames(current_dat)),
                      grep("Inhib", colnames(current_dat)),
                      grep("Neuron", colnames(current_dat)))
  reorderedCols <- c(neuronPosition, setdiff(1:19, neuronPosition))
  current_dat <- current_dat[ ,reorderedCols]
  # Put NE neurons before the 5-HT ones
  current_dat <- current_dat[ ,c(1:12, 14, 13, 15:19)]
  
  italicnames <- lapply(
    rownames(current_dat),
    function(x) bquote(italic(.(x))))
  
  # Print
  pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE,
           breaks = seq(0.02, 4, length.out = 101),
           color = viridis::viridis(100),
           main=paste0("\t\t\t",i,"-DE genes, 19 LC cell classes (means)"),
           labels_col = as.expression(italicnames),
           fontsize=12, fontsize_row = 15, fontsize_col=14)
  grid::grid.text(label="log2-\nExprs", x=0.97, y=0.57, gp=grid::gpar(fontsize=10))
  
  ## or medians version:
  current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMedians(dat[list.visiumDE[[i]], ii])))
  # For some reason rownames aren't kept:
  rownames(current_dat) <- list.visiumDE[[i]]
  # Set neuronal pops first
  neuronPosition <- c(grep("Excit", colnames(current_dat)),
                      grep("Inhib", colnames(current_dat)),
                      grep("Neuron", colnames(current_dat)))
  reorderedCols <- c(neuronPosition, setdiff(1:19, neuronPosition))
  current_dat <- current_dat[ ,reorderedCols]
  # Put NE neurons before the 5-HT ones
  current_dat <- current_dat[ ,c(1:12, 14, 13, 15:19)]
  # Print
  pheatmap(t(current_dat), cluster_rows = FALSE, cluster_cols = FALSE,
           breaks = seq(0.02, 4, length.out = 101),
           color = viridis::viridis(100),
           main=paste0("\t\t\t\t",i,"-DE genes, 19 LC cell classes (medians)"),
           labels_col = as.expression(italicnames),
           fontsize=12, fontsize_row = 15, fontsize_col=14)
  grid::grid.text(label="log2-\nExprs", x=0.97, y=0.57, gp=grid::gpar(fontsize=10))
  
  dev.off()
}




## Within-region pairwise correlation? ===============
#  - How does it compare to cluster modularity ratio?

## Using cluster-vs-all-others stats, gather Cohen's d (~ t-statistic as employed in previous work)
#    See https://github.com/MarioniLab/scran/issues/57 for discussion/recommendation

library(pheatmap)
library(RColorBrewer)


clusters.d <- lapply(markers.lc.t.1vAll, function(x){x[[2]]$std.logFC})
# Take unique of top 100 / cell class
topGenes.space <- unique(unlist(lapply(clusters.d, function(x){head(names(x), 100)})))
# of length 1599

# out of curiosity:
length(intersect(topGenes.space, hdgs.lc))  # only 615


# Subset and combine into matrix
clusters.d <- lapply(clusters.d, function(x){x[topGenes.space]})
clusters.d.mat <- do.call(cbind, clusters.d)

# Correlation matrix
cor.d.topMarkers <- cor(clusters.d.mat)

# Reproduce cluster modularity (from '04_clusterAgglomeration.R') so it's in the same pdf
library(bluster)

load(here("processed_data","SCE","graph_communities_glmpcamnn_LC.rda"), verbose=T)
    # snn.gr.glmpcamnn, clusters.glmpcamnn

mod.ratio.merged.HC <- pairwiseModularity(graph = snn.gr.glmpcamnn,
                                          clusters = sce.lc$cellType.merged,
                                          as.ratio=TRUE)

# Remove 'ambig.lowNTx_' annotations
mod.ratio.merged.HC <- mod.ratio.merged.HC[-grep("ambig.lowNTx_", rownames(mod.ratio.merged.HC)),
                                           -grep("ambig.lowNTx_", colnames(mod.ratio.merged.HC))]

theSeq.all = seq(-1, 1, by = 0.025)
my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)


pdf(here("plots", "snRNA-seq", "heatmaps_correlation_CohensD_vs_clusterMod_19cellClasses.pdf"))
pheatmap(cor.d.topMarkers,
         color=my.col.all,
         breaks=theSeq.all,
         display_numbers=T, fontsize_number=7,
         fontsize_row=9, fontsize_col=9,
         main="Correlation of within-region cluster Cohen's D \n (top 100 enriched genes/cluster space)")
grid::grid.text(label="Pearson's r",x=0.96,y=0.48, gp=grid::gpar(fontsize=7))

# Cluster modularity
pheatmap(log2(mod.ratio.merged.HC+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100),
         main="Modularity ratio for 19 HC-merged clusters in LC (n=3)",
         fontsize_row=7, fontsize_col=7, angle_col=90,
         display_numbers=T, number_format="%.1f", na_col="darkgrey")
grid::grid.text(label="log2(ratio)",x=0.97,y=0.64, gp=grid::gpar(fontsize=7))
dev.off()




## With five regions from Neuron paper (Tran, Maynard, et al. 2021) ==========
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/markers-stats_all-regions-combined_SN-LEVEL-1vAll_MNT2021.rda", verbose=T)
    # FMstats.list, sampleNumNuclei, readme.mnt, ref.sampleInfo
readme.mnt
    #[1] "These stats are from region-specific specificity modeling (cluster-vs-all-others) at the single-nucleus level with
    #     'scran::findMarkers()'. The t-statistic is computed by sqrt(N.nuclei) * std.logFC."

## Create matrix of D's with region:subcluster identifiers
ds.list <- lapply(FMstats.list, function(x){
  sapply(x, function(y){y$std.logFC})
}
)
# Add back in region suffix
for(i in names(ds.list)){
  colnames(ds.list[[i]]) <- paste0(colnames(ds.list[[i]]), "_", i)
}
# cbind
ds.fullMat <- do.call(cbind, ds.list)

# Now change rownames to Ensembl IDs

# Load one of that project's SCEs to change back to/match Ensembl IDs:
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_MNT2021.rda", verbose=T)
# sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, clusterRefTab.dlpfc, ref.sampleInfo, annotationTab.dlpfc, cell_colors

table(rownames(ds.fullMat) %in% rowData(sce.dlpfc)$Symbol.uniq)
rownames(ds.fullMat) <- rowData(sce.dlpfc)$gene_id[match(rownames(ds.fullMat),
                                                         rowData(sce.dlpfc)$Symbol.uniq)]

# Now subset for intersecting gene space (basically what's shared b/tw two versions of the reference GTF)
geneU <- intersect(rownames(ds.fullMat), rownames(markers.lc.t.1vAll[[1]][[2]]))
    # of length 26514

ds.fullMat <- ds.fullMat[geneU, ]

# As above for just the within-LC correlations
clusters.d.lc <- lapply(markers.lc.t.1vAll, function(x){x[[2]]$std.logFC})
clusters.d.lc <- lapply(clusters.d.lc, function(x){x[geneU]})
clusters.d.lc <- do.call(cbind, clusters.d.lc)

table(rownames(ds.fullMat) == rownames(clusters.d.lc))  # all TRUE

# now cbind them
colnames(clusters.d.lc) <- paste0(colnames(clusters.d.lc), "_LocusC")
ds.fullMat <- cbind(ds.fullMat, clusters.d.lc)

# Top-100 iteration (ultimately the better version; less noise-driven)
clus_specific_indices = mapply(function(t) {
  oo = order(t, decreasing = TRUE)[1:100]
},
as.data.frame(ds.fullMat)
)
clus_ind = unique(as.numeric(clus_specific_indices))
length(clus_ind)  # originally 3715 unique in the reported 107 classes
    # now 4091

ds.defined <- ds.fullMat[clus_ind, ]
cor_d_defined <- cor(ds.defined)

my.col.all <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq.all)-1)


## To have iterations split up by neuronal vs non-neuronal pops
ds.defined.neu <- ds.defined
for(i in c("As", "Micro", "Endo", "Mural","Oligo", "OPC", "Tcell", "Macro")){
  ds.defined.neu <- ds.defined.neu[ ,-grep(i, colnames(ds.defined.neu))]
}

cor_d_defined.neu <- cor(ds.defined.neu)


# Add some cluster info for add'l heatmap annotations
clusterInfo <- data.frame(region=ss(colnames(ds.defined.neu), "_",3))
rownames(clusterInfo) <- colnames(ds.defined.neu)
clusterInfo$region[c(82:83)] <- "LocusC"

# Make some region cols
clusterCols <- list(region=tableau10medium[1:6])
names(clusterCols[["region"]]) <- levels(as.factor(clusterInfo$region))


# Non-neuronal cell classes
ds.defined.non <- ds.defined
glia.idx <- NA
for(i in c("As", "Micro", "Endo", "Mural","Oligo", "OPC", "Tcell", "Macro")){
  glia.idx <- c(glia.idx, grep(i, colnames(ds.defined.non)))
}
# Rm the empty NA
glia.idx <- glia.idx[-1]
glia.idx <- unique(glia.idx) # One duplicate ('Endo.Mural' LC population)
ds.defined.non <- ds.defined.non[ ,glia.idx]
cor_d_defined.non <- cor(ds.defined.non)

# Add some cluster info for add'l heatmap annotations
clusterInfo.glia <- data.frame(region=ifelse(is.na(ss(colnames(ds.defined.non), "_",3)),
                                             ss(colnames(ds.defined.non), "_",2),
                                             ss(colnames(ds.defined.non), "_",3))
)

rownames(clusterInfo.glia) <- colnames(ds.defined.non)

cor_d_defined.non <- cor(ds.defined.non)


## Print all iterations ===
pdf(here("plots", "snRNA-seq", "heatmaps_correlation_CohensD_LC-19_vs_RewardCircuitry-107-Neuron2021.pdf"), height=9, width=9)
# FULL correlatino heatmap ('6' brain regions)
pheatmap(cor_d_defined,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=4.5, fontsize_col=4.5,
         main="Correlation of Cohen's D across LC and 107 reward circuitry populations \n (Tran, Maynard, et al. Neuron 2021)"
         )

# Neuronal subset
pheatmap(cor_d_defined.neu,
         annotation_col=clusterInfo,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=5.5, fontsize_col=5.5,
         main="Correlation of Cohen's D across LC and 107 reward circuitry populations \n (Tran, Maynard, et al. Neuron 2021)"
         )
# Neuronal subset with numbers
pheatmap(cor_d_defined.neu,
         annotation_col=clusterInfo,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=5.5, fontsize_col=5.5,
         display_numbers=TRUE, fontsize_number=2.4,
         main="Correlation of Cohen's D across LC and 107 reward circuitry populations \n (Tran, Maynard, et al. Neuron 2021)"
)

# Non-neuronal subset with numbers
pheatmap(cor_d_defined.non,
         annotation_col=clusterInfo.glia,
         annotation_colors=clusterCols,
         #show_colnames=FALSE,
         color=my.col.all,
         breaks=theSeq.all,
         fontsize_row=5.7, fontsize_col=5.7,
         display_numbers=TRUE, fontsize_number=3.5,
         main="Correlation of Cohen's D across LC and 107 reward circuitry populations \n (Tran, Maynard, et al. Neuron 2021)"
)
dev.off()


## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
# [1] "2022-06-20 14:27:38 EDT"
proc.time()
#     user    system   elapsed 
#  133.151     6.222 13231.422 
options(width = 120)
session_info()
#─ Session info ────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-06-20
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────
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
# cli                    3.3.0    2022-04-25 [2] CRAN (R 4.1.2)
# cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
# DBI                    1.1.3    2022-06-18 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
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
# igraph                 1.3.2    2022-06-13 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.3   2022-04-07 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
# MASS                   7.3-56   2022-03-23 [3] CRAN (R 4.1.2)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.1.2)
# metapod                1.2.0    2021-10-26 [2] Bioconductor
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# nlme                   3.1-157  2022-03-25 [3] CRAN (R 4.1.2)
# pheatmap             * 1.0.12   2019-01-04 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer         * 1.1-3    2022-04-03 [2] CRAN (R 4.1.2)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.7 2022-06-09 [2] CRAN (R 4.1.2)
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
# segmented              1.6-0    2022-05-31 [1] CRAN (R 4.1.2)
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
# ───────────────────────────────────────────────────────────────────────────────────


