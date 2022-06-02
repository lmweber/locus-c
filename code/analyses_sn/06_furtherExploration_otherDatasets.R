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

# Because scale_color_manual() has started keeping and printing unused
#     factor levels as 'NA' in the legend:
colors2print <- cell_colors.lc[-grep("ambig.lowNTx_", names(cell_colors.lc))]


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





## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
# 
proc.time()
#     user    system   elapsed 
# 
options(width = 120)
session_info()


