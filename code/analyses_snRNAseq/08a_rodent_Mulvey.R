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

# function from rafalib package
splitit <- function(x) split(seq(along = x), x)


# calculate matrix of values: genes x clusters

cell_idx <- splitit(sce$label)

dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$gene_name

mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[human_genes, i]))))


# subset neuronal clusters
clus_NE <- c(
  6  ## NE
)
clus_neuronalExcNE <- c(
  29, 21, 20, 23,  ## excitatory
  26, 17, 14, 1, 8, 7, 24, 18,  ## inhibitory
  16  ## 5-HT
)

mat_NE <- mat[clus_NE, ]
mat_neuronalExcNE <- mat[clus_neuronalExcNE, ]

# grand mean across neuronal clusters excluding NE
meanNeuronalExcNE <- colMeans(mat_neuronalExcNE)


stopifnot(all(names(mat_NE) == colnames(mat_neuronalExcNE)))


ord <- order(mat_NE)

df <- data.frame(
  gene = names(mat_NE), 
  NE = mat_NE, 
  other = meanNeuronalExcNE) %>% 
  pivot_longer(
    cols = c("NE", "other"), 
    names_to = "clusters", 
    values_to = "mean"
  ) %>% 
  mutate(gene = factor(gene, levels = names(mat_NE)[ord]))


# -------------
# plot boxplots
# -------------

ggplot(df, aes(x = mean, y = gene, color = clusters, fill = clusters)) + 
  geom_boxplot(alpha = 0.75, outlier.size = 0.5) + 
  #scale_color_manual(values = pal_rev, name = "annotation", 
  #                   guide = guide_legend(reverse = TRUE)) + 
  #scale_fill_manual(values = pal_rev, name = "annotation", 
  #                  guide = guide_legend(reverse = TRUE)) + 
  labs(x = "mean logcounts") + 
  ggtitle("Mulvey et al. genes") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 9, face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "enrichment_Mulvey_NEcluster")
ggsave(paste0(fn, ".pdf"), width = 4.75, height = 6.5)
ggsave(paste0(fn, ".png"), width = 4.75, height = 6.5)

