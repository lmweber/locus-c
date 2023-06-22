########################################################################
# LC snRNA-seq analyses: cluster identification using known marker genes
# Lukas Weber, Jun 2023
# including code from Matthew N Tran
########################################################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)


dir_plots <- here("plots", "singleNucleus", "05b_cluster_identification")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_clustering_main")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)

# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# --------------
# Color palettes
# --------------

# color palettes from scater package
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", 
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5", 
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F", 
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9", 
                    "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA")

colors <- unique(c(tableau20, tableau10medium))


# ------------
# Marker genes
# ------------

# marker gene lists from Matthew N Tran and Keri Martinowich

markers_all <- c(
  # neuron markers
  "SNAP25", "SYT1", 
  # excitatory (glutamatergic) neuron markers
  "SLC17A6", "SLC17A8", ## "SLC17A7"; alternative names: "VGLUT1", "VGLUT2", "VGLUT3"
  # inhibitory (GABAergic) neuron markers
  "GAD1", "GAD2", 
  # inhibitory subpopulations (from Keri Martinowich 2022-07-22)
  "SST", "KIT", "CALB1", "CALB2", "TAC1", "CNR1", "PVALB", "CORT", "VIP", "NPY", "CRHBP", "CCK", ## "HTR3A" (not present in data)
  # cholinergic neurons
  "SLC5A7", "CHAT", "ACHE", "BCHE", "SLC18A3", "PRIMA1", 
  # additional miscellaneous markers for comparison with Luskin et al. (2022) (from Keri Martinowich 2022-10-20)
  "CALCA", "CALCR", "CARTPT", "GAL", "PENK", "PNOC", "SLC6A5", 
  # NE neuron markers
  "DBH", "TH", "SLC6A2", "SLC18A2", ## "DDC", "GCH1"
  # 5-HT (serotonin) markers
  "TPH2", "SLC6A4", ## "TPH1"
  # 5-HT other marker genes
  "FEV", "HTR1A", "HTR1B", "HCRTR2", 
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


# ------------------------------------
# Marker expression heatmap by cluster
# ------------------------------------

# markers for heatmap
markers <- c(
  "SNAP25", "SYT1", 
  "SLC17A6", 
  "GAD1", "GAD2", 
  "SST", "KIT", "CALB1", "CALB2", "TAC1", "CNR1", 
  "DBH", "TH", "SLC6A2", "SLC18A2", 
  "TPH2", "SLC6A4", 
  "GFAP", "AQP4", 
  "CLDN5", "FLT1", "RBPMS", 
  "CD163", "C3", 
  "MBP", 
  "PDGFRA", "VCAN"
)


# marker labels
marker_labels <- c(
  rep("neuron", 2), 
  rep("excitatory", 1), 
  rep("inhibitory", 8), 
  rep("NE", 4), 
  rep("5HT", 2), 
  rep("astrocytes", 2), 
  rep("endothelial_mural", 3), 
  rep("macrophages_microglia", 2), 
  rep("oligodendrocytes", 1), 
  rep("OPCs", 2))

marker_labels <- 
  factor(marker_labels, levels = unique(marker_labels))


# colors: selected from tableau20 and tableau10medium
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
  excitatory = c(30, 23, 26, 7), 
  inhibitory = c(24, 25, 14, 4, 8, 20, 17), 
  neurons_ambiguous = c(18, 19, 3, 22, 13, 1, 2), 
  NE = 6, 
  `5HT` = 21, 
  astrocytes = c(5, 29), 
  endothelial_mural = 16, 
  macrophages_microglia = 15, 
  oligodendrocytes = c(10, 12, 28, 9), 
  OPCs = c(27, 11))
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


# colors: selected from tableau20 and tableau10medium
colors_clusters <- list(population = c(
  excitatory = "#1F77B4", 
  inhibitory = "#AEC7E8", 
  neurons_ambiguous = "gray60", 
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


# heatmap data

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$label)
dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$gene_name


hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers, i]))))


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

hm


# save heatmap
fn <- file.path(dir_plots, "clustering_heatmap")

pdf(paste0(fn, ".pdf"), width = 8, height = 6.5)
hm
dev.off()

png(paste0(fn, ".png"), width = 8 * 200, height = 6.5 * 200, res = 200)
hm
dev.off()


# -----------------------------------------------------
# Marker expression heatmap: inhibitory and cholinergic
# -----------------------------------------------------

# marker expression heatmap with additional inhibitory and cholinergic markers

# re-using some objects from previous heatmap above


# markers for heatmap
markers_extended <- c(
  "SNAP25", "SYT1", 
  "SLC17A6", 
  "GAD1", "GAD2", 
  "SST", "KIT", "CALB1", "CALB2", "TAC1", "CNR1", "PVALB", "CORT", "VIP", "NPY", "CRHBP", "CCK", 
  "CALCA", "CALCR", "CARTPT", "GAL", "PENK", "PNOC", "SLC6A5", 
  "SLC5A7", "CHAT", "ACHE", "BCHE", "SLC18A3", "PRIMA1", 
  "SLC6A3", "ALDH1A1", "SLC26A7", 
  "DBH", "TH", "SLC6A2", "SLC18A2", 
  "TPH2", "SLC6A4", 
  "FEV", "HTR1A", "HTR1B", "HCRTR2", 
  "GFAP", "AQP4", 
  "CLDN5", "FLT1", "RBPMS", 
  "CD163", "C3", 
  "MBP", 
  "PDGFRA", "VCAN"
)


# marker labels
marker_labels_extended <- c(
  rep("neuron", 2), 
  rep("excitatory", 1), 
  rep("inhibitory", 2), 
  rep("inhibitory_subtypes", 12), 
  rep("miscellaneous", 7), 
  rep("cholinergic", 6), 
  rep("dopaminergic", 3), 
  rep("NE", 4), 
  rep("5HT", 2), 
  rep("5HT_other", 4), 
  rep("astrocytes", 2), 
  rep("endothelial_mural", 3), 
  rep("macrophages_microglia", 2), 
  rep("oligodendrocytes", 1), 
  rep("OPCs", 2))

marker_labels_extended <- 
  factor(marker_labels_extended, levels = unique(marker_labels_extended))


# colors: selected from tableau20 and tableau10medium
colors_markers_extended <- list(marker = c(
  neuron = "black", 
  excitatory = "#1F77B4", 
  inhibitory = "#AEC7E8", 
  inhibitory_subtypes = "lightskyblue", 
  miscellaneous = "darkslateblue", 
  cholinergic = "gold", 
  dopaminergic = "deeppink", 
  NE = "#D62728", 
  `5HT` = "#9467BD", 
  `5HT_other` = "purple4", 
  astrocytes = "#FF7F0E", 
  endothelial_mural = "#98DF8A", 
  macrophages_microglia = "#8C564B", 
  oligodendrocytes = "#9EDAE5", 
  OPCs = "#17BECF"))


# heatmap data

hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_extended, i]))))


# column annotation
col_ha <- columnAnnotation(
  marker = marker_labels_extended, 
  show_annotation_name = FALSE, 
  show_legend = TRUE, 
  col = colors_markers_extended)


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
  column_split = marker_labels_extended, 
  column_names_gp = gpar(fontface = "italic"), 
  rect_gp = gpar(col = "gray50", lwd = 0.5))

hm


# save heatmap
fn <- file.path(dir_plots, "clustering_heatmap_extended")

pdf(paste0(fn, ".pdf"), width = 13.5, height = 6.5)
hm
dev.off()

png(paste0(fn, ".png"), width = 13.5 * 200, height = 6.5 * 200, res = 200)
hm
dev.off()


# ---------
# Plot UMAP
# ---------

# UMAP of clustering

label_merged <- fct_collapse(colData(sce)$label, 
  excitatory = as.character(cluster_pops[[1]]), 
  inhibitory = as.character(cluster_pops[[2]]), 
  neurons_ambiguous = as.character(cluster_pops[[3]]), 
  NE = as.character(cluster_pops[[4]]), 
  `5HT` = as.character(cluster_pops[[5]]), 
  astrocytes = as.character(cluster_pops[[6]]), 
  endothelial_mural = as.character(cluster_pops[[7]]), 
  macrophages_microglia = as.character(cluster_pops[[8]]), 
  oligodendrocytes = as.character(cluster_pops[[9]]), 
  OPCs = as.character(cluster_pops[[10]]))

label_merged <- fct_relevel(label_merged, 
  c("excitatory", "inhibitory", "neurons_ambiguous", "NE", "5HT", "astrocytes", 
    "endothelial_mural", "macrophages_microglia", "oligodendrocytes", "OPCs"))

colData(sce)$label_merged <- label_merged


plotReducedDim(sce, dimred = "UMAP", colour_by = "label_merged") + 
  scale_color_manual(values = colors_clusters[[1]], name = "clusters (merged)") + 
  theme_classic() + 
  ggtitle("LC clustering")

fn <- file.path(dir_plots, "clustering_UMAP_merged")
ggsave(paste0(fn, ".pdf"), width = 6.5, height = 4.75)
ggsave(paste0(fn, ".png"), width = 6.5, height = 4.75)


# ---------------------
# Plot number of nuclei
# ---------------------

# plot number of nuclei per cluster (with ambiguous nuclei)

tbl <- table(colLabels(sce))

df <- data.frame(
  cluster = factor(names(tbl), levels = cluster_pops_order), 
  n_nuclei = as.numeric(tbl), 
  prop_nuclei = as.numeric(tbl) / sum(as.numeric(tbl)), 
  population = unname(cluster_pops_rev))

ggplot(df, aes(x = cluster, y = n_nuclei, fill = population, 
               label = paste0(format(round(prop_nuclei * 100, 1), nsmall = 1), "%"))) + 
  geom_bar(stat = "identity") + 
  geom_text(vjust = -0.25, size = 2, fontface = "bold") + 
  scale_fill_manual(values = colors_clusters$population, name = "population") + 
  labs(y = "number of nuclei") + 
  ggtitle("Number of nuclei per cluster") + 
  theme_bw()

fn <- here(dir_plots, paste0("numberNuclei_perCluster"))
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# plot number of nuclei per cluster (without ambiguous nuclei)

tbl <- table(colLabels(sce))

df <- data.frame(
  cluster = factor(names(tbl), levels = cluster_pops_order), 
  n_nuclei = as.numeric(tbl), 
  population = unname(cluster_pops_rev))

df <- df[!(df$population == "neurons_ambiguous"), ]
df$population <- droplevels(df$population)
colors_noAmbiguous <- colors_clusters$population[!(names(colors_clusters$population) == "neurons_ambiguous")]

df$prop_nuclei <- df$n_nuclei / sum(df$n_nuclei)

ggplot(df, aes(x = cluster, y = n_nuclei, fill = population, 
               label = paste0(format(round(prop_nuclei * 100, 1), nsmall = 1), "%"))) + 
  geom_bar(stat = "identity") + 
  geom_text(vjust = -0.25, size = 2, fontface = "bold") + 
  scale_fill_manual(values = colors_noAmbiguous, name = "population") + 
  labs(y = "number of nuclei") + 
  ggtitle("Number of nuclei per cluster (excluding ambiguous)") + 
  theme_bw()

fn <- here(dir_plots, paste0("numberNuclei_perCluster_noAmbiguous"))
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 4)
ggsave(paste0(fn, ".png"), width = 6.75, height = 4)


# plot number of nuclei per cluster (with ambiguous nuclei): per sample

tbl <- table(colLabels(sce), colData(sce)$Sample)
df <- data.frame(
  cluster = factor(rownames(tbl), levels = cluster_pops_order), 
  population = unname(cluster_pops_rev), 
  as.data.frame(cbind(as.matrix(tbl)))
)
df <- df %>% 
  pivot_longer(cols = -c(cluster, population), names_to = "sample", values_to = "n_nuclei") %>% 
  mutate(sample = factor(sample, levels = c("Br6522_LC", "Br2701_LC", "Br8079_LC")))

ggplot(df, aes(x = cluster, y = n_nuclei, fill = population)) + 
  facet_wrap(~sample, nrow = 3, scales = "free_x") + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = colors_clusters$population, name = "population") + 
  ylab("number of nuclei") + 
  ggtitle("Number of nuclei per cluster: per sample") + 
  theme_bw()

fn <- here(dir_plots, paste0("numberNuclei_perCluster_perSample"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 9)
ggsave(paste0(fn, ".png"), width = 6, height = 9)


# ----------------------
# QC metrics per cluster
# ----------------------

df <- as.data.frame(colData(sce)) %>% 
  select(sum, detected, subsets_Mito_percent, Sample, label, label_merged) %>% 
  mutate(label = factor(label, levels = cluster_pops_order))


# sum UMI counts per cluster

ggplot(df, aes(x = label, y = sum, color = label_merged)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = colors_clusters$population, name = "population") + 
  labs(x = "cluster", 
       y = "sum UMI counts") + 
  ggtitle("Sum UMI counts per cluster") + 
  theme_bw()

fn <- here(dir_plots, paste0("QCmetrics_perCluster_sumUMI"))
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# detected genes per cluster

ggplot(df, aes(x = label, y = detected, color = label_merged)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = colors_clusters$population, name = "population") + 
  labs(x = "cluster", 
       y = "detected genes") + 
  ggtitle("Detected genes per cluster") + 
  theme_bw()

fn <- here(dir_plots, paste0("QCmetrics_perCluster_detectedGenes"))
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# mitochondrial proportion per cluster

ggplot(df, aes(x = label, y = subsets_Mito_percent, color = label_merged)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = colors_clusters$population, name = "population") + 
  labs(x = "cluster", 
       y = "mitochondrial percentage") + 
  ggtitle("Mitochondrial percentage per cluster") + 
  theme_bw()

fn <- here(dir_plots, paste0("QCmetrics_perCluster_mitochondrial"))
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# additional details: NE neuron cluster vs. other neurons (excluding ambiguous)

df_metrics_neur <- df %>% 
  filter(label_merged %in% c("excitatory", "inhibitory", "NE", "5HT")) %>% 
  mutate(label_NE = ifelse(label_merged == "NE", "NE", "other_neurons")) %>% 
  mutate(label_NE = factor(label_NE, levels = c("NE", "other_neurons"), 
                           labels = c("NE neurons", "other neurons\n(exc. ambiguous)"))) %>% 
  pivot_longer(cols = c("sum", "detected"), 
               names_to = "metric", values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("sum", "detected"), 
                         labels = c("sum UMIs", "detected genes")))

ggplot(df_metrics_neur, aes(x = label_NE, y = value, fill = label_NE)) + 
  facet_wrap(~metric, nrow = 1) + 
  geom_violin() + 
  scale_y_log10() + 
  scale_fill_manual(values = c("#D62728", "cornflowerblue"), name = "populations") + 
  ggtitle("NE neurons vs. other neurons") + 
  theme_bw() + 
  theme(axis.title = element_blank())

fn <- here(dir_plots, paste0("QCmetrics_NEvsOtherNeuron"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 4.25)
ggsave(paste0(fn, ".png"), width = 6, height = 4.25)


# additional details: NE neuron cluster vs. all other (excluding ambiguous neurons)

df_metrics_all <- df %>% 
  filter(label_merged != "neurons_ambiguous") %>% 
  mutate(label_NE = ifelse(label_merged == "NE", "NE", "other")) %>% 
  mutate(label_NE = factor(label_NE, levels = c("NE", "other"), 
                           labels = c("NE neurons", "all other\n(exc. ambiguous)"))) %>% 
  pivot_longer(cols = c("sum", "detected"), 
               names_to = "metric", values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("sum", "detected"), 
                         labels = c("sum UMIs", "detected genes")))

ggplot(df_metrics_all, aes(x = label_NE, y = value, fill = label_NE)) + 
  facet_wrap(~metric, nrow = 1) + 
  geom_violin() + 
  scale_y_log10() + 
  scale_fill_manual(values = c("#D62728", "black"), name = "populations") + 
  ggtitle("NE neurons vs. all other") + 
  theme_bw() + 
  theme(axis.title = element_blank())

fn <- here(dir_plots, paste0("QCmetrics_NEvsAllOther"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 4.25)
ggsave(paste0(fn, ".png"), width = 6, height = 4.25)


# -----------
# Save object
# -----------

# save object with merged cluster labels

fn_out <- here("processed_data", "SCE", "sce_clustering_merged")
saveRDS(sce, paste0(fn_out, ".rds"))

