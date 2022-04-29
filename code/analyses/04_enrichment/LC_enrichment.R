##############################################
# LC project
# Script to plot enrichment of LC marker genes
# Lukas Weber, Apr 2022
##############################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "04_enrichment")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_qualityControlled.rds")
spe <- readRDS(fn_spe)

dim(spe)


# ------------------------------------
# calculate enrichment of marker genes
# ------------------------------------

# calculate enrichment of marker genes in manually annotated regions

# select marker genes of interest
marker_genes <- c("TH", "SLC6A2")


# LC regions
table(colData(spe)$annot_region)
# WM regions
table(!colData(spe)$annot_region)
# individual spots
table(colData(spe)$annot_spot)
# individual non-spots
table(!colData(spe)$annot_spot)


sample_ids <- levels(colData(spe)$sample_id)
sample_ids

enrichment <- matrix(NA, nrow = length(sample_ids), ncol = length(marker_genes))
rownames(enrichment) <- sample_ids
colnames(enrichment) <- marker_genes

enrichment_LC <- enrichment_WM <- enrichment_spots <- enrichment_nonspots <- enrichment


# select manually annotated LC regions
spe_LC <- spe[, colData(spe)$annot_region]
dim(spe_LC)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(marker_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_LC)[rowData(spe_LC)$gene_name == marker_genes[j], 
                                      colData(spe_LC)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_LC[i, j] <- mean_ij
  }
}


# select manually annotated WM regions
spe_WM <- spe[, !colData(spe)$annot_region]
dim(spe_WM)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(marker_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_WM)[rowData(spe_WM)$gene_name == marker_genes[j], 
                                      colData(spe_WM)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_WM[i, j] <- mean_ij
  }
}


# select manually annotated individual spots
spe_spots <- spe[, colData(spe)$annot_spot]
dim(spe_spots)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(marker_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_spots)[rowData(spe_spots)$gene_name == marker_genes[j], 
                                         colData(spe_spots)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_spots[i, j] <- mean_ij
  }
}


# select manually annotated individual nonspots
spe_nonspots <- spe[, !colData(spe)$annot_spot]
dim(spe_nonspots)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(marker_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_nonspots)[rowData(spe_nonspots)$gene_name == marker_genes[j], 
                                            colData(spe_nonspots)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_nonspots[i, j] <- mean_ij
  }
}


# ---------------
# plot enrichment
# ---------------

df_enrichment_LC <- as.data.frame(enrichment_LC)
df_enrichment_LC$region <- "LC"
df_enrichment_LC$sample <- rownames(df_enrichment_LC)

df_enrichment_WM <- as.data.frame(enrichment_WM)
df_enrichment_WM$region <- "WM"
df_enrichment_WM$sample <- rownames(df_enrichment_WM)

df_enrichment_spots <- as.data.frame(enrichment_spots)
df_enrichment_spots$region <- "spots"
df_enrichment_spots$sample <- rownames(df_enrichment_spots)

df_enrichment_nonspots <- as.data.frame(enrichment_nonspots)
df_enrichment_nonspots$region <- "nonspots"
df_enrichment_nonspots$sample <- rownames(df_enrichment_nonspots)


df_enrichment_LC <- pivot_longer(df_enrichment_LC, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")
df_enrichment_WM <- pivot_longer(df_enrichment_WM, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")
df_enrichment_spots <- pivot_longer(df_enrichment_spots, cols = -c(region, sample), 
                                    names_to = "gene", values_to = "mean")
df_enrichment_nonspots <- pivot_longer(df_enrichment_nonspots, cols = -c(region, sample), 
                                       names_to = "gene", values_to = "mean")


pal1 <- c("darkorange", "purple4")
pal2 <- c("dodgerblue", "black")


# order genes
meds <- colMedians(enrichment_spots, useNames = TRUE)
genes_ordered <- names(sort(meds, decreasing = TRUE))
all(genes_ordered %in% marker_genes)
length(genes_ordered) == length(marker_genes)

df1 <- 
  full_join(df_enrichment_LC, df_enrichment_WM) %>% 
  mutate(region = factor(region, levels = c("LC", "WM"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids)) %>% 
  mutate(gene = factor(gene, levels = genes_ordered)) %>% 
  as.data.frame()

df2 <- 
  full_join(df_enrichment_spots, df_enrichment_nonspots) %>% 
  mutate(region = factor(region, levels = c("spots", "nonspots"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids)) %>% 
  mutate(gene = factor(gene, levels = genes_ordered)) %>% 
  as.data.frame()


# LC regions vs. WM regions
p <- ggplot(df1, aes(x = gene, y = mean, color = region)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = pal1) + 
  labs(y = "mean logcounts") + 
  ggtitle("Enrichment: LC domains") + 
  theme_bw()

fn <- here(dir_plots, "enrichment_selected_domains")
ggsave(paste0(fn, ".pdf"), plot = p, width = 5, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 5, height = 4)


# annotated spots vs. non-spots
p <- ggplot(df2, aes(x = gene, y = mean, color = region)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = pal2) + 
  labs(y = "mean logcounts") + 
  ggtitle("Enrichment: spots") + 
  theme_bw()

fn <- here(dir_plots, "enrichment_selected_spots")
ggsave(paste0(fn, ".pdf"), plot = p, width = 5, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 5, height = 4)

