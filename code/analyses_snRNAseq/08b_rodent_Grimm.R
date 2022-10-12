#########################################################
# LC snRNA-seq analyses: Grimm et al. (2004) rodent genes
# Lukas Weber, Oct 2022
#########################################################


library(here)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)


# directory to save plots
dir_plots <- here("plots", "singleNucleus", "08_rodent_genes")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_clustering_secondary")
sce <- readRDS(paste0(fn, ".rds"))

dim(sce)

table(colData(sce)$Sample)

# number of nuclei per cluster and sample
table(colLabels(sce))
table(colLabels(sce), colData(sce)$Sample)


# ------------------
# load rat gene list
# ------------------

# list of rat genes for LC from Grimm et al. (2004), Figure 4

# converted to human genes by searching individually at https://www.ncbi.nlm.nih.gov/gene/

human_genes <- c(
  "CBR3", "DNAH5", "SERPINE1", "LAYN", "TPH2", "RPH3AL", "NGB", "CYB561", 
  "GNAS", "SLC31A1", "TCP1", "PPIC", "COLEC10", "RAB3B", "MAOA", "PCBP3", 
  "TSPAN12", "FBP1", "DBH", "SERPINF1", "TXK", "SEC16B", "TRAF1", "PTGES", 
  "GGT5", "MMP2", "MCAM", "TFAP2A", "ACSL1", "UPB1", "UCP3", "COL5A1", 
  "ALDH1A1", "FBN1")

length(human_genes)

# 34 out of 34 genes present in SCE object
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
  30, 23, 26, 7, ## excitatory
  24, 25, 14, 4, 8, 20, 17, ## inhibitory
  21  ## 5-HT
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

# order genes by expression within NE cluster (and within other clusters for ties)
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

pal <- c("red", "gray30")
pal_rev <- rev(pal)

# note: not showing outliers
ggplot(df, aes(x = logcounts, y = gene, color = clusters, fill = clusters)) + 
  geom_boxplot(alpha = 0.75, outlier.shape = NA) + 
  scale_color_manual(values = pal_rev, 
                     guide = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values = pal_rev, 
                    guide = guide_legend(reverse = TRUE)) + 
  ggtitle("Grimm et al. genes") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.y = element_text(size = 9, face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "enrichment_Grimm_NEcluster")
ggsave(paste0(fn, ".pdf"), width = 4.75, height = 5.75)
ggsave(paste0(fn, ".png"), width = 4.75, height = 5.75)

