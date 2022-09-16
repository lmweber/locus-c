###################################
# LC snRNA-seq analyses: DE testing
# Lukas Weber, Sep 2022
###################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(tidyr)
library(forcats)
library(tibble)
library(ggplot2)
library(ggnewscale)
library(ggrepel)


dir_plots <- here("plots", "snRNAseq", "07_DEtesting")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_clustering_merged")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)

# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# ----------
# DE testing
# ----------

# pairwise DE testing between neuronal clusters


# store total UMI counts per gene
rowData(sce)$sum_gene <- rowSums(counts(sce))


# select neuronal clusters (excluding ambiguous)
clus_neurons <- c(
  29, ## excitatory
  26, 17, 14, 1, 8, 7, 24, 18,  ## inhibitory
  6,  ## NE
  16  ## 5-HT
)
ix_neurons <- colData(sce)$label %in% clus_neurons
table(ix_neurons)

sce <- sce[, ix_neurons]
dim(sce)

# remove empty levels
colData(sce)$label <- droplevels(colData(sce)$label)


# calculate DE tests
# note: not blocking by sample since 1 out of 3 samples contains almost zero NE neurons
# testing for genes with log-fold-changes significantly greater than 1 (lfc = 1, direction = "up")
marker_info <- findMarkers(
  sce, 
  groups = colData(sce)$label, 
  lfc = 1, 
  direction = "up", 
  row.data = rowData(sce)[, c("gene_id", "gene_name", "sum_gene")]
)
# marker_info <- scoreMarkers(
#   sce, 
#   groups = colData(sce)$label, 
#   row.data = rowData(sce)[, c("gene_id", "gene_name", "sum_gene")]
# )

marker_info

# NE neurons
marker_info[["6"]]
#View(as.data.frame(marker_info[["6"]]))

# 5-HT neurons
marker_info[["16"]]
#View(as.data.frame(marker_info[["16"]]))


# check some known genes

ix_known <- which(marker_info[["6"]]$gene_name %in% c("TH", "SLC6A2", "DBH"))
marker_info[["6"]][ix_known, ]

ix_known <- which(marker_info[["16"]]$gene_name %in% c("TPH2", "SLC6A4"))
marker_info[["16"]][ix_known, ]


# significant DE genes
table(marker_info[["6"]]$FDR < 1)
table(marker_info[["6"]]$FDR < 0.05)
table(marker_info[["6"]]$FDR < 1e-6)
table(marker_info[["6"]]$FDR < 1e-100)


# using mean Cohen's d statistic (i.e. standardized log-fold-change) across
# pairwise comparisons as ranking metric (from OSCA chapter 6)

# ordered <- marker_info[["6"]][order(marker_info[["6"]]$mean.logFC.cohen, decreasing = TRUE), ]
# head(ordered)
# #View(as.data.frame(ordered))
# 
# ordered <- marker_info[["16"]][order(marker_info[["16"]]$mean.logFC.cohen, decreasing = TRUE), ]
# head(ordered)
# #View(as.data.frame(ordered))


# -------------
# Volcano plots
# -------------

# NE neurons

fdr <- marker_info[["6"]]$FDR
logfc <- marker_info[["6"]]$summary.logFC

names(fdr) <- names(logfc) <- marker_info[["6"]]$gene_name


# identify significant genes (low FDR and high logFC)
thresh_fdr <- 1e-6  ## absolute scale
thresh_logfc <- log2(2)  ## log2 scale
sig <- (fdr < thresh_fdr) & (logfc > thresh_logfc)

# highly significant thresholds
highlysig <- (fdr < 1e-20) & (logfc > log2(4))

# number of significant genes
table(sig)
table(highlysig)


df <- data.frame(
  gene = names(fdr), 
  FDR = fdr, 
  log2FC = logfc, 
  sig = sig, 
  highlysig = highlysig
)

pal <- c("black", "red")


# volcano plot with labels
set.seed(123)
ggplot(df, aes(x = log2FC, y = -log10(FDR), color = highlysig, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$highlysig, ], size = 0.5) + 
  geom_text_repel(data = df[df$highlysig, ], 
                  size = 1.5, nudge_y = 0.1, 
                  force = 0.1, force_pull = 0.1, min.segment.length = 0.1, 
                  max.overlaps = 20) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(1e-20), lty = "dashed", color = "royalblue") + 
  #geom_vline(xintercept = -thresh_logfc, lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(4), lty = "dashed", color = "royalblue") + 
  ggtitle("NE neuron cluster vs. all neuronal clusters") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "DEtesting_NEneuronsVsAllOtherNeuronalClusters")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)

