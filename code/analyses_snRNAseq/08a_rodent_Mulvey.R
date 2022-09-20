##########################################################
# LC snRNA-seq analyses: Mulvey et al. (2018) rodent genes
# Lukas Weber, Sep 2022
##########################################################


library(here)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)


# directory to save plots
dir_plots <- here("plots", "snRNAseq", "08_rodent_genes")


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


# --------------------
# load mouse gene list
# --------------------

# list of 45 mouse genes for LC from Mulvey et al. (2018), Figure 2A

fn <- here("inputs", "Mulvey_markers", "Mulvey2018_Fig2A_markers.txt")
mouse_genes <- read.table(fn)[, 1]

length(mouse_genes)
mouse_genes


# converted to human genes using biomaRt (code below) and by searching for genes
# individually at: https://www.ncbi.nlm.nih.gov/gene/

human_genes <- c(
  "AGTR1", "ASB4", "CALCR", "CALR3", "CHODL", "CHRNA6", "CILP", "CYB561", 
  "DBH", "DDC", "DLK1", "ADGRE1", "EYA2", "SHISAL2B", "FAM183A", "FIBCD1", 
  "GAL", "GCH1", "GLRA2", "GNG4", "GPX3", "GTF2A1L", "HCRTR1", "IGSF5", "MAOA", 
  "MRAP2", "MYOM2", "NEUROG2", "SLC9B2", "NXPH4", "OVGP1", "PCBD1", "PHOX2A", 
  "PHOX2B", "PLA2G4D", "PTGER2", "SLC18A2", "SLC31A1", "SLC6A2", "STBD1", 
  "SYT17", "TH", "TM4SF1", "TM4SF5", "TRAF3IP2")

human_genes <- sort(human_genes)
length(human_genes)

# 42 out of 45 genes present in SCE object
sum(human_genes %in% rowData(sce)$gene_name)

# keep genes that are present in SCE object
ix_keep <- which(human_genes %in% rowData(sce)$gene_name)
human_genes <- human_genes[ix_keep]

stopifnot(length(human_genes) == sum(human_genes %in% rowData(sce)$gene_name))


# --------------------
# calculate enrichment
# --------------------

# enrichment in NE cluster vs. other neuronal clusters

ix_genes <- rowData(sce)$gene_name %in% human_genes
table(ix_genes)

dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$gene_name

dat <- dat[ix_genes, ]


# subset neuronal clusters
clus_NE <- c(
  6  ## NE
)
clus_neuronalExcNE <- c(
  29, 21, 20, 23,  ## excitatory
  26, 17, 14, 1, 8, 7, 24, 18,  ## inhibitory
  16  ## 5-HT
)

ix_NE <- colData(sce)$label %in% clus_NE
ix_other <- colData(sce)$label %in% clus_neuronalExcNE

table(ix_NE)
table(ix_other)

stopifnot(length(ix_NE) == ncol(dat))
stopifnot(length(ix_other) == ncol(dat))

clusters <- rep("nonneuron", ncol(dat))
clusters[ix_NE] <- "NE"
clusters[ix_other] <- "other"

# order genes by expression within NE cluster
ord <- order(rowMeans(dat[, ix_NE]))
genes_ord <- rownames(dat)[ord]

# data frame for boxplots
df <- data.frame(
  t(dat[ord, ]), 
  clusters = clusters) %>% 
  filter(clusters %in% c("NE", "other")) %>% 
  pivot_longer(cols = -clusters, names_to = "gene", values_to = "logcounts") %>% 
  mutate(clusters = factor(clusters, levels = c("other", "NE"))) %>%  ## reversed for plot
  mutate(gene = factor(gene, levels = genes_ord))


# -------------
# plot boxplots
# -------------

pal <- c("darkmagenta", "gray60")
pal_rev <- rev(pal)

# note: not showing outliers
ggplot(df, aes(x = logcounts, y = gene, color = clusters, fill = clusters)) + 
  geom_boxplot(alpha = 0.75, outlier.shape = NA) + 
  scale_color_manual(values = pal_rev, 
                     guide = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values = pal_rev, 
                    guide = guide_legend(reverse = TRUE)) + 
  ggtitle("Mulvey et al. genes") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 9, face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "enrichment_Mulvey_NEcluster")
ggsave(paste0(fn, ".pdf"), width = 4.75, height = 6.5)
ggsave(paste0(fn, ".png"), width = 4.75, height = 6.5)

