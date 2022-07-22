########################################################################
# LC snRNA-seq analyses: cluster identification using known marker genes
# Lukas Weber, July 2022
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
library(RColorBrewer)
library(pheatmap)


dir_plots <- here("plots", "snRNAseq_alt", "06_cluster_identification", "known_markers")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE_alt", "sce_clustering")
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


# -------------------------
# Marker expression heatmap
# -------------------------

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
  "SLC17A7", "SLC17A6", ## alternative names: VGLUT1=SLC17A7, VGLUT2=SLC17A6
  # inhibitory (GABAergic) neuron markers
  "GAD1", "GAD2", 
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
types = c("neuron", "excitatory", "inhibitory", "NE", "5HT", "cholinergic", 
          "astrocytes", "endothelial_mural", "macrophages_microglia", 
          "oligodendrocytes", "OPCs")
annotation <- data.frame(
  type = factor(c(
    rep(types[[1]], 2), 
    rep(types[[2]], 2), 
    rep(types[[3]], 2), 
    rep(types[[4]], 4), 
    rep(types[[5]], 2), 
    rep(types[[6]], 6), 
    rep(types[[7]], 2), 
    rep(types[[8]], 3), 
    rep(types[[9]], 2), 
    rep(types[[10]], 1), 
    rep(types[[11]], 2)), 
    levels = types, 
    labels = types))
rownames(annotation) <- markers_broad


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

# create heatmap
p <- pheatmap(t(current_dat), 
              annotation = annotation, 
              cluster_rows = FALSE, cluster_cols = FALSE, 
              breaks = seq(0.02, 4, length.out = 101), 
              color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100), 
              main = "LC clusters marker expression (medians)", 
              labels_col = as.expression(italicnames), 
              angle_col = 90, 
              fontsize = 12, fontsize_row = 15, fontsize_col = 14)
#grid::grid.text(label = "log2-\nExprs", x = 0.96, y = 0.63, gp = grid::gpar(fontsize = 10))

fn <- here(dir_plots, paste0("clustersMarkersExpression_heatmap.pdf"))
pdf(fn, width = 11, height = 9)
#par(mar = c(5,8,4,2))
p
dev.off()

fn <- here(dir_plots, paste0("clustersMarkersExpression_heatmap.png"))
png(fn, width = 11 * 200, height = 9 * 200, res = 200)
p
dev.off()


# ----------------
# number of nuclei
# ----------------

# plot number of nuclei per cluster: all samples

df <- data.frame(
  cluster = names(table(colLabels(sce))), 
  n_nuclei = as.numeric(table(colLabels(sce))))
df$cluster <- factor(df$cluster, levels = df$cluster)

ggplot(df, aes(x = cluster, y = n_nuclei)) + 
  geom_bar(stat = "identity", fill = "navy") + 
  ylab("number of nuclei") + 
  ggtitle("Number of nuclei per cluster: all samples") + 
  theme_bw()

fn <- here(dir_plots, paste0("numberNuclei_allSamples"))
ggsave(paste0(fn, ".pdf"), width = 7, height = 4)
ggsave(paste0(fn, ".png"), width = 7, height = 4)


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
  facet_wrap(~sample) + 
  geom_bar(stat = "identity", fill = "navy") + 
  ylab("number of nuclei") + 
  ggtitle("Number of nuclei per cluster: per sample") + 
  theme_bw()

fn <- here(dir_plots, paste0("numberNuclei_perSample"))
ggsave(paste0(fn, ".pdf"), width = 16, height = 3.5)
ggsave(paste0(fn, ".png"), width = 16, height = 3.5)

