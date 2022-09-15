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

# DE testing between NE neuron cluster and all other clusters


# store total UMI counts per gene
rowData(sce)$sum_gene <- rowSums(counts(sce))


# remove ambiguous neurons cluster
ix_ambiguous <- colData(sce)$label_merged == "neurons_ambiguous"
table(ix_ambiguous)

sce <- sce[, !ix_ambiguous]
dim(sce)

# remove empty level
colData(sce)$label_merged <- droplevels(colData(sce)$label_merged)


# calculate DE tests
# note: not blocking by sample since 1 out of 3 samples contains almost zero NE neurons
marker_info <- findMarkers(
  sce, 
  groups = colData(sce)$label_merged, 
  lfc = 1, 
  direction = "up", 
  row.data = rowData(sce)[, c("gene_id", "gene_name", "sum_gene")]
)

marker_info$NE
#View(as.data.frame(marker_info$NE))

# check some known genes
ix_known <- which(marker_info$NE$gene_name %in% c("TH", "SLC6A2", "DBH"))
marker_info$NE[ix_known, ]


# significant DE genes
table(marker_info$NE$FDR < 1)
table(marker_info$NE$FDR < 0.05)
table(marker_info$NE$FDR < 1e-6)
table(marker_info$NE$FDR < 1e-100)


# volcano plot

fdr <- marker_info$NE$FDR
logfc <- marker_info$NE$summary.logFC

names(fdr) <- names(logfc) <- marker_info$NE$gene_name


# identify significant genes (low FDR and high logFC)
thresh_fdr <- 1e-100  ## absolute scale
thresh_logfc <- log2(2)  ## log2 scale
sig <- (fdr < thresh_fdr) & (abs(logfc) > thresh_logfc)

# number of significant genes
table(sig)


df <- data.frame(
  gene = names(fdr), 
  FDR = fdr, 
  logFC = logfc, 
  sig = sig
)

pal <- c("black", "red")


# volcano plot with labels
set.seed(123)
ggplot(df, aes(x = logFC, y = -log10(FDR), color = sig, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$sig, ], size = 0.5) + 
  #geom_text_repel(data = df[df$sig, ], 
  #                size = 1.5, nudge_y = 0.1, 
  #                force = 0.1, force_pull = 0.1, min.segment.length = 0.1, 
  #                max.overlaps = 20) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(thresh_fdr), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -thresh_logfc, lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = thresh_logfc, lty = "dashed", color = "royalblue") + 
  ggtitle("NE neuron cluster vs. all other nuclei") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "DEtesting_NEneuronsVsAllOther")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)

