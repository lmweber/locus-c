########################################################################
# LC snRNA-seq analyses: cluster identification using known marker genes
# Lukas Weber, Sep 2022
# using code by Matthew N Tran
########################################################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(jaffelab)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)


dir_plots <- here("plots", "snRNAseq", "06_cluster_identification", "known_markers")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_clustering")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)

# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# ------------------------------
# Marker expression violin plots
# ------------------------------

# plotting function from Matthew N Tran
plotExpressionCustom <- function(sce, features, features_name, anno_name = "cellType", 
                                 point_alpha = 0.2, point_size = 0.7, ncol = 2, xlab = NULL, 
                                 exprs_values = "logcounts", scales = "free_y", swap_rownames = NULL) {
  scater::plotExpression(sce, 
                         exprs_values = exprs_values, 
                         features = features, 
                         x = anno_name, 
                         colour_by = anno_name, 
                         ncol = ncol, 
                         xlab = xlab, 
                         point_alpha = point_alpha, 
                         point_size = point_size, 
                         add_legend = FALSE, 
                         scales = scales, 
                         swap_rownames = swap_rownames) + 
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median, 
                 geom = "crossbar", 
                 width = 0.3) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          strip.text = element_text(face = "italic")) + 
    ggtitle(label = paste0(features_name, " markers"))
}


# marker gene list from Matthew N Tran
markers.mathys.tran = list(
  'neuron' = c('SYT1', 'SNAP25', 'GRIN1'), 
  'excit_neuron' = c('CAMK2A', 'NRGN','SLC17A7', 'SLC17A6', 'SLC17A8'), 
  'inhib_neuron' = c('GAD1', 'GAD2', 'SLC32A1'), 
  # Norepinephrine & serotonergic markers
  'neuron.NE' = c("TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC"), #SLC6A3 - saw no DAT
  'neuron.5HT' = c("SLC6A4", "TPH1", "TPH2", "DDC"), 
  # SERT, serotonin T (aka 5-HTT)
  'monoamine.metab' = c("COMT", "MAOA", "MAOB"), 
  # MSN markers
  'MSNs.pan' = c("PPP1R1B","BCL11B"), # "CTIP2")
  'MSNs.D1' = c("DRD1", "PDYN", "TAC1"), 
  'MSNs.D2' = c("DRD2", "PENK"), 
  ## Non-neuronal
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


# palettes from scater
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", 
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5", 
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F", 
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", 
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2", 
                    "#CDCC5D", "#6DCCDA")
colors <- unique(c(tableau20, tableau10medium))


# plotting code from Matthew N Tran
fn <- here(dir_plots, paste0("clustersMarkersExpression_byPopulation_violins.pdf"))
pdf(fn, height = 6, width = 12)
for(i in 1:length(markers.mathys.tran)) {
  print(
    plotExpressionCustom(sce = sce, 
                         exprs_values = "logcounts", 
                         features = markers.mathys.tran[[i]], 
                         features_name = names(markers.mathys.tran)[[i]], 
                         anno_name = "label", 
                         ncol = 2, point_alpha = 0.4, point_size = 0.9, 
                         scales = "free_y", swap_rownames = "gene_name") + 
      ggtitle(label = paste0("LC clusters: ", 
                             names(markers.mathys.tran)[[i]], " markers")) + 
      theme(plot.title = element_text(size = 12), 
            axis.text.x = element_text(size=7)) + 
      scale_color_manual(values = colors)
  )
}
dev.off()


# --------------------------------------------
# Marker expression heatmap: broad populations
# --------------------------------------------

# # marker gene list from Matthew N Tran
# markers_broad <- c(
#   'SNAP25','SLC17A7','SLC17A6','GAD1','GAD2', 
#   # NE neuron markers
#   "TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC", 
#   # serotonergic markers (includes DDC but repetitive)
#   "SLC6A4", "TPH2",  # (TPH1 not expressed by these clusters)
#   ## Non-neuronal:
#   # Astro
#   'AQP4','GFAP', 
#   # Endo, Mural (RBPMS)
#   'CLDN5','FLT1','RBPMS', 
#   # Macrophage, Microglia
#   'CD163','C3', 
#   # Oligo
#   'MBP', 
#   # OPC
#   'PDGFRA','VCAN'
# )

# marker gene list from Matthew N Tran - updated LW
markers_broad <- c(
  # neuron markers
  "SNAP25", "SYT1", 
  # excitatory (glutamatergic) neuron markers
  "SLC17A7", "SLC17A6", "SLC17A8", ## alternative names: VGLUT1, VGLUT2, VGLUT3
  # inhibitory (GABAergic) neuron markers
  "GAD1", "GAD2", 
  # inhibitory subpopulations (from Keri Martinowich 2022-07-22)
  "SST", "KIT", "CALB1", "CALB2", "TAC1", "CNR1", ## "PVALB", "CORT", "VIP", "NPY", "CRHBP", "CCK", "HTR3A" (not present in our data)
  # NE neuron markers
  "DBH", "TH", "SLC6A2", "DDC", ## "SLC18A2", "GCH1", 
  # 5-HT (serotonin) markers
  "TPH2", "SLC6A4", ## "TPH1", 
  # astrocytes
  "GFAP", "AQP4", 
  # endothelial / mural (RBPMS)
  "CLDN5", "FLT1", "RBPMS", 
  # macrophages / microglia
  "CD163", "C3", 
  # oligodendrocytes
  "MBP", 
  # OPCs
  "PDGFRA", "VCAN"
)

# annotation data frame for heatmap
types_broad = c("neuron", "excitatory", "inhibitory", "inhibitory_subtypes", 
                "NE", "5HT", 
                "astrocytes", "endothelial_mural", "macrophages_microglia", 
                "oligodendrocytes", "OPCs")
annotation_broad <- data.frame(
  cluster = factor(c(
    rep(types_broad[[1]], 2), 
    rep(types_broad[[2]], 3), 
    rep(types_broad[[3]], 2), 
    rep(types_broad[[4]], 6), 
    rep(types_broad[[5]], 4), 
    rep(types_broad[[6]], 2), 
    rep(types_broad[[7]], 2), 
    rep(types_broad[[8]], 3), 
    rep(types_broad[[9]], 2), 
    rep(types_broad[[10]], 1), 
    rep(types_broad[[11]], 2)), 
    levels = types_broad, 
    labels = types_broad))
rownames(annotation_broad) <- markers_broad


# plotting code from Matthew N Tran

cell.idx <- splitit(sce$label)
dat <- as.matrix(assay(sce, "logcounts"))
rownames(dat) <- rowData(sce)$gene_name


## Medians version
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMedians(dat[markers_broad, ii])))
# For some reason rownames aren't kept
rownames(current_dat) <- markers_broad
# # Set neuronal pops first
# neuronPosition <- c(grep("Excit", colnames(current_dat)), 
#                     grep("Inhib", colnames(current_dat)), 
#                     grep("Neuron", colnames(current_dat)))
# reorderedCols <- c(neuronPosition, setdiff(1:19, neuronPosition))
# current_dat <- current_dat[ ,reorderedCols]
# # Put NE neurons before the 5-HT ones
# current_dat <- current_dat[ ,c(1:12, 14, 13, 15:19)]

italicnames <- lapply(
  rownames(current_dat), 
  function(x) bquote(italic(.(x))))

# order clusters
order_broad <- c(29, 26, 17, 14, 1, 8, 7, 24, 18, 21, 20, 23, 19, 30, 13, 3, 5, 2, 6, 16, 22, 25, 15, 11, 9, 10, 27, 4, 28, 12)
#order_broad <- c(27, 2, 29, 28, 17, 19, 21, 11, 5, 16, 3, 4, 30, 6, 7, 8, 15, 18, 22, 26, 9, 25, 1, 24, 13, 10, 14, 20, 12, 23)

set.seed(1)
pal <- sample(tableau20)[seq_along(types_broad)]
names(pal) <- types_broad

# create heatmap
p <- pheatmap(t(current_dat[, order_broad]), 
              annotation = annotation_broad, annotation_colors = list(cluster = pal), 
              cluster_rows = FALSE, cluster_cols = FALSE, 
              #breaks = seq(0.02, 4, length.out = 101), 
              color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100), 
              main = "LC clusters marker expression (medians)", 
              labels_col = as.expression(italicnames), 
              angle_col = 90, 
              fontsize = 12, fontsize_row = 15, fontsize_col = 14)
#grid::grid.text(label = "log2-\nExprs", x = 0.96, y = 0.63, gp = grid::gpar(fontsize = 10))

fn <- here(dir_plots, paste0("clustersMarkersExpression_heatmap_medians.pdf"))
pdf(fn, width = 10, height = 8)
#par(mar = c(5,8,4,2))
p
dev.off()

fn <- here(dir_plots, paste0("clustersMarkersExpression_heatmap_medians.png"))
png(fn, width = 10 * 200, height = 8 * 200, res = 200)
p
dev.off()


## Means version
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[markers_broad, ii])))
rownames(current_dat) <- markers_broad

italicnames <- lapply(
  rownames(current_dat), 
  function(x) bquote(italic(.(x))))

# create heatmap
p <- pheatmap(t(current_dat[, order_broad]), 
              annotation = annotation_broad, annotation_colors = list(cluster = pal), 
              cluster_rows = FALSE, cluster_cols = FALSE, 
              color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100), 
              main = "LC clusters marker expression (means)", 
              labels_col = as.expression(italicnames), 
              angle_col = 90, 
              fontsize = 12, fontsize_row = 15, fontsize_col = 14)

fn <- here(dir_plots, paste0("clustersMarkersExpression_heatmap_means.pdf"))
pdf(fn, width = 10, height = 8)
p
dev.off()

fn <- here(dir_plots, paste0("clustersMarkersExpression_heatmap_means.png"))
png(fn, width = 10 * 200, height = 8 * 200, res = 200)
p
dev.off()


# ----------------------------------------
# Alternative heatmap using ComplexHeatmap
# ----------------------------------------

hm_mat <- t(do.call(cbind, lapply(cell.idx, function(i) rowMeans(dat[markers_broad, i]))))

# markers to show
markers <- c(
  "SNAP25", "SYT1", 
  "SLC17A7", "SLC17A6", 
  "GAD1", "GAD2", 
  "SST", "KIT", "CALB1", "CALB2", "TAC1", "CNR1", 
  "DBH", "TH", "SLC6A2", "DDC", 
  "TPH2", "SLC6A4", 
  "GFAP", "AQP4", 
  "CLDN5", "FLT1", "RBPMS", 
  "CD163", "C3", 
  "MBP", 
  "PDGFRA", "VCAN"
)

hm_mat <- hm_mat[, markers]


# marker labels
marker_labels <- c(
  rep("neuron", 2), 
  rep("excitatory", 2), 
  rep("inhibitory", 8), 
  rep("NE", 4), 
  rep("5HT", 2), 
  rep("astrocytes", 2), 
  rep("endothelial_mural", 3), 
  rep("macrophages_microglia", 2), 
  rep("oligodendrocytes", 1), 
  rep("OPCs", 2))

marker_labels <- factor(marker_labels, levels = unique(marker_labels))

# colors: from tableau20 and tableau10medium
colors_markers <- list(marker = c(
  neuron = "black", 
  excitatory = "#1F77B4", 
  inhibitory = "#AEC7E8", 
  NE = "#D62728", 
  `5HT` = "#9467BD", 
  astrocytes = "#FF7F0E", 
  endothelial_mural = "#98DF8A", 
  macrophages_microglia = "#8C564B", 
  oligodendrocytes = "#9EDAE5", 
  OPCs = "#17BECF"))


# cluster labels
cluster_pops <- list(
  excitatory = 29, 
  inhibitory = c(26, 17, 14, 1, 8, 7, 24, 18), 
  neurons_ambiguous = c(21, 20, 23, 19, 30, 13, 3, 5, 2), 
  NE = 6, 
  `5HT` = 16, 
  astrocytes = c(22, 25), 
  endothelial_mural = 15, 
  macrophages_microglia = 11, 
  oligodendrocytes = c(9, 10, 27, 4), 
  OPCs = c(28, 12))
# cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))
# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops), times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]

cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))


# second set of cluster labels
neuron_pops <- ifelse(
  cluster_pops_rev %in% c("excitatory", "inhibitory", "NE", "5HT", "neurons_ambiguous"), 
  "neuronal", "non-neuronal") %>% 
  factor(., levels = c("neuronal", "non-neuronal"))


# colors: from tableau20 and tableau10medium
colors_clusters <- list(population = c(
  excitatory = "#1F77B4", 
  inhibitory = "#AEC7E8", 
  neurons_ambiguous = "gray50", 
  NE = "#D62728", 
  `5HT` = "#9467BD", 
  astrocytes = "#FF7F0E", 
  endothelial_mural = "#98DF8A", 
  macrophages_microglia = "#8C564B", 
  oligodendrocytes = "#9EDAE5", 
  OPCs = "#17BECF"))

colors_neurons <- list(class = c(
  neuronal = "black", 
  `non-neuronal` = "gray90"
))


# number of nuclei per cluster
n <- table(colLabels(sce))


# row annotation
row_ha <- rowAnnotation(
  n = anno_barplot(as.numeric(n), gp = gpar(fill = "navy"), border = FALSE), 
  class = neuron_pops, 
  population = cluster_pops_rev, 
  show_annotation_name = FALSE, 
  col = c(colors_clusters, colors_neurons))

# column annotation
col_ha <- columnAnnotation(
  marker = marker_labels, 
  show_annotation_name = FALSE, 
  show_legend = FALSE, 
  col = colors_markers)


hm <- Heatmap(
  hm_mat, 
  name = "mean\nlogcounts", 
  column_title = "LC clusters mean marker expression", 
  column_title_gp = gpar(fontface = "bold"), 
  col = brewer.pal(n = 7, "OrRd"), 
  right_annotation = row_ha, 
  bottom_annotation = col_ha, 
  row_order = cluster_pops_order, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  row_split = cluster_pops_rev, 
  row_title = NULL, 
  column_split = marker_labels, 
  column_names_gp = gpar(fontface = "italic"), 
  rect_gp = gpar(col = "gray50", lwd = 0.5))


# save heatmap
fn <- file.path(dir_plots, "clustering_heatmap_complex")

pdf(paste0(fn, ".pdf"), width = 8, height = 6.5)
hm
dev.off()

png(paste0(fn, ".png"), width = 8 * 200, height = 6.5 * 200, res = 200)
hm
dev.off()


# ---------
# UMAP plot
# ---------

# UMAP of clustering

labels_merged <- fct_collapse(colData(sce)$label, 
  excitatory = "29", 
  inhibitory = c("26", "17", "14", "1", "8", "7", "24", "18"), 
  neurons_ambiguous = c("21", "20", "23", "19", "30", "13", "3", "5", "2"), 
  NE = "6", 
  `5HT` = "16", 
  astrocytes = c("22", "25"), 
  endothelial_mural = "15", 
  macrophages_microglia = "11", 
  oligodendrocytes = c("9", "10", "27", "4"), 
  OPCs = c("28", "12"))

labels_merged <- fct_relevel(labels_merged, 
  c("excitatory", "inhibitory", "neurons_ambiguous", "NE", "5HT", "astrocytes", 
    "endothelial_mural", "macrophages_microglia", "oligodendrocytes", "OPCs"))

colData(sce)$labels_merged <- labels_merged


plotReducedDim(sce, dimred = "UMAP", colour_by = "labels_merged") + 
  scale_color_manual(values = colors_clusters[[1]], name = "clusters (merged)") + 
  theme_classic() + 
  ggtitle("LC clustering")

fn <- file.path(dir_plots, "clustering_UMAP_merged")
ggsave(paste0(fn, ".pdf"), width = 6.5, height = 4.75)
ggsave(paste0(fn, ".png"), width = 6.5, height = 4.75)


# ----------------------------------------------
# Marker expression heatmap: inhibitory subtypes
# ----------------------------------------------

# # marker gene list from Matthew N Tran
# markers_broad <- c(
#   'SNAP25','SLC17A7','SLC17A6','GAD1','GAD2', 
#   # NE neuron markers
#   "TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC", 
#   # serotonergic markers (includes DDC but repetitive)
#   "SLC6A4", "TPH2",  # (TPH1 not expressed by these clusters)
#   ## Non-neuronal:
#   # Astro
#   'AQP4','GFAP', 
#   # Endo, Mural (RBPMS)
#   'CLDN5','FLT1','RBPMS', 
#   # Macrophage, Microglia
#   'CD163','C3', 
#   # Oligo
#   'MBP', 
#   # OPC
#   'PDGFRA','VCAN'
# )

# marker gene list from Matthew N Tran - updated LW
markers_inhib <- c(
  # neuron markers
  "SNAP25", "SYT1", 
  # excitatory (glutamatergic) neuron markers
  "SLC17A7", "SLC17A6", "SLC17A8", ## alternative names: VGLUT1, VGLUT2, VGLUT3
  # inhibitory (GABAergic) neuron markers
  "GAD1", "GAD2", 
  # inhibitory subpopulations (from Keri Martinowich 2022-07-22)
  "PVALB", "SST", "CORT", "KIT", "VIP", "NPY", "CRHBP", "CALB1", "TAC1", "CCK", "CNR1", "CALB2", ## "HTR3A" (not present in our data)
  # NE neuron markers
  "DBH", "TH", "SLC6A2", "DDC", ## "SLC18A2", "GCH1", 
  # 5-HT (serotonin) markers
  "TPH2", "SLC6A4", ## "TPH1", 
  # cholinergic neurons
  "SLC5A7", "CHAT", "ACHE", "BCHE", "SLC18A3", "PRIMA1", 
  # astrocytes
  "GFAP", "AQP4", 
  # endothelial / mural (RBPMS)
  "CLDN5", "FLT1", "RBPMS", 
  # macrophages / microglia
  "CD163", "C3", 
  # oligodendrocytes
  "MBP", 
  # OPCs
  "PDGFRA", "VCAN"
)

# annotation data frame for heatmap
types_inhib = c("neuron", "excitatory", "inhibitory", "inhibitory_subtypes", 
                "NE", "5HT", "cholinergic", 
                "astrocytes", "endothelial_mural", "macrophages_microglia", 
                "oligodendrocytes", "OPCs")
annotation_inhib <- data.frame(
  cluster = factor(c(
    rep(types_inhib[[1]], 2), 
    rep(types_inhib[[2]], 3), 
    rep(types_inhib[[3]], 2), 
    rep(types_inhib[[4]], 12), 
    rep(types_inhib[[5]], 4), 
    rep(types_inhib[[6]], 2), 
    rep(types_inhib[[7]], 6), 
    rep(types_inhib[[8]], 2), 
    rep(types_inhib[[9]], 3), 
    rep(types_inhib[[10]], 2), 
    rep(types_inhib[[11]], 1), 
    rep(types_inhib[[12]], 2)), 
    levels = types_inhib, 
    labels = types_inhib))
rownames(annotation_inhib) <- markers_inhib


# plotting code from Matthew N Tran

cell.idx <- splitit(sce$label)
dat <- as.matrix(assay(sce, "logcounts"))
rownames(dat) <- rowData(sce)$gene_name


## Medians version
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMedians(dat[markers_inhib, ii])))
# For some reason rownames aren't kept
rownames(current_dat) <- markers_inhib
# # Set neuronal pops first
# neuronPosition <- c(grep("Excit", colnames(current_dat)), 
#                     grep("Inhib", colnames(current_dat)), 
#                     grep("Neuron", colnames(current_dat)))
# reorderedCols <- c(neuronPosition, setdiff(1:19, neuronPosition))
# current_dat <- current_dat[ ,reorderedCols]
# # Put NE neurons before the 5-HT ones
# current_dat <- current_dat[ ,c(1:12, 14, 13, 15:19)]

italicnames <- lapply(
  rownames(current_dat), 
  function(x) bquote(italic(.(x))))

# create heatmap
p <- pheatmap(t(current_dat), annotation = annotation_inhib, 
              cluster_rows = FALSE, cluster_cols = FALSE, 
              #breaks = seq(0.02, 4, length.out = 101), 
              color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100), 
              main = "LC clusters marker expression (medians)", 
              labels_col = as.expression(italicnames), 
              angle_col = 90, 
              fontsize = 12, fontsize_row = 15, fontsize_col = 14)
#grid::grid.text(label = "log2-\nExprs", x = 0.96, y = 0.63, gp = grid::gpar(fontsize = 10))

fn <- here(dir_plots, paste0("clustersMarkersExpression_inhibSubtypes_heatmap_medians.pdf"))
pdf(fn, width = 12, height = 8)
#par(mar = c(5,8,4,2))
p
dev.off()

fn <- here(dir_plots, paste0("clustersMarkersExpression_inhibSubtypes_heatmap_medians.png"))
png(fn, width = 12 * 200, height = 8 * 200, res = 200)
p
dev.off()


## Means version
current_dat <- do.call(cbind, lapply(cell.idx, function(ii) rowMeans(dat[markers_inhib, ii])))
rownames(current_dat) <- markers_inhib

italicnames <- lapply(
  rownames(current_dat), 
  function(x) bquote(italic(.(x))))

# create heatmap
p <- pheatmap(t(current_dat), annotation = annotation_inhib, 
              cluster_rows = FALSE, cluster_cols = FALSE, 
              color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100), 
              main = "LC clusters marker expression (means)", 
              labels_col = as.expression(italicnames), 
              angle_col = 90, 
              fontsize = 12, fontsize_row = 15, fontsize_col = 14)

fn <- here(dir_plots, paste0("clustersMarkersExpression_inhibSubtypes_heatmap_means.pdf"))
pdf(fn, width = 12, height = 8)
p
dev.off()

fn <- here(dir_plots, paste0("clustersMarkersExpression_inhibSubtypes_heatmap_means.png"))
png(fn, width = 12 * 200, height = 8 * 200, res = 200)
p
dev.off()


# ----------------
# number of nuclei
# ----------------

# plot number of nuclei per cluster: all samples

df <- data.frame(
  cluster = names(table(colLabels(sce))), 
  n_nuclei = as.numeric(table(colLabels(sce))))
df$cluster <- factor(df$cluster, levels = order_broad)

ggplot(df, aes(x = cluster, y = n_nuclei)) + 
  geom_bar(stat = "identity", fill = "navy") + 
  ylab("number of nuclei") + 
  ggtitle("Number of nuclei per cluster: all samples") + 
  theme_bw()

fn <- here(dir_plots, paste0("numberNuclei_allSamples"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)


# plot number of nuclei per cluster: per sample

tbl <- table(colLabels(sce), colData(sce)$Sample)
df <- as.data.frame(cbind(
  cluster = as.numeric(rownames(tbl)), 
  as.matrix(tbl)
))
df <- df %>% 
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "n_nuclei") %>% 
  mutate(cluster = factor(cluster, levels = unique(cluster))) %>% 
  mutate(sample = factor(sample, levels = c("Br6522_LC", "Br2701_LC", "Br8079_LC")))

ggplot(df, aes(x = cluster, y = n_nuclei)) + 
  facet_wrap(~sample, nrow = 3, scales = "free_x") + 
  geom_bar(stat = "identity", fill = "navy") + 
  ylab("number of nuclei") + 
  ggtitle("Number of nuclei per cluster: per sample") + 
  theme_bw()

fn <- here(dir_plots, paste0("numberNuclei_perSample"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 9)
ggsave(paste0(fn, ".png"), width = 6, height = 9)

