#########################################
# LC analyses: enrichment of marker genes
# Lukas Weber, Oct 2022
#########################################

# module load conda_R/4.2
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "Visium", "04_enrichment")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_logcounts.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)


# ------------------------------------
# calculate enrichment of marker genes
# ------------------------------------

# enrichment defined in terms of mean logcounts per spot in LC vs. non-LC regions

# select marker genes of interest
# for NE neurons and 5-HT neurons
marker_genes <- c("TH", "SLC6A2", "TPH2", "SLC6A4")

# LC region annotations
table(colData(spe)$annot_region)


enrichment <- matrix(NA, nrow = length(sample_ids), ncol = length(marker_genes))
rownames(enrichment) <- sample_ids
colnames(enrichment) <- marker_genes

enrichment_LC <- enrichment_non <- enrichment


# select annotated LC regions
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

# select annotated non-LC regions
spe_non <- spe[, !colData(spe)$annot_region]
dim(spe_non)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(marker_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_non)[rowData(spe_non)$gene_name == marker_genes[j], 
                                       colData(spe_non)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_non[i, j] <- mean_ij
  }
}


# ---------------
# plot enrichment
# ---------------

df_enrichment_LC <- as.data.frame(enrichment_LC)
df_enrichment_LC$region <- "LC"
df_enrichment_LC$sample <- rownames(df_enrichment_LC)

df_enrichment_non <- as.data.frame(enrichment_non)
df_enrichment_non$region <- "non"
df_enrichment_non$sample <- rownames(df_enrichment_non)


df_enrichment_LC <- pivot_longer(df_enrichment_LC, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")
df_enrichment_non <- pivot_longer(df_enrichment_non, cols = -c(region, sample), 
                                  names_to = "gene", values_to = "mean")


df <- 
  full_join(df_enrichment_LC, df_enrichment_non) %>% 
  mutate(regions = factor(region, 
                          levels = c("LC", "non"), 
                          labels = c("LC regions", "non-LC regions"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids))


# subset marker genes for NE neurons and 5-HT neurons

df_NE <- df %>% 
  filter(gene %in% c("TH", "SLC6A2")) %>% 
  mutate(gene = factor(gene, levels = c("TH", "SLC6A2"))) %>% 
  as.data.frame()

df_5HT <- df %>% 
  filter(gene %in% c("TPH2", "SLC6A4")) %>% 
  mutate(gene = factor(gene, levels = c("TPH2", "SLC6A4"))) %>% 
  as.data.frame()


pal_NE <- c("red", "gray30")
pal_5HT <- c("darkmagenta", "gray30")


# ------------------------
# plots without sample IDs
# ------------------------

# plot marker genes for NE neurons
set.seed(123)
ggplot(df_NE) + 
  geom_boxplot(aes(x = gene, y = mean, color = regions, fill = regions), 
               alpha = 0.5, outlier.shape = NA) + 
  geom_point(aes(x = gene, y = mean, color = regions), 
             position = position_jitterdodge()) + 
  scale_color_manual(values = pal_NE, name = "annotation") + 
  scale_fill_manual(values = pal_NE, name = "annotation") + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Enrichment") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(face = "italic"))

fn <- here(dir_plots, "enrichment_annotatedRegions_NEmarkers")
ggsave(paste0(fn, ".pdf"), width = 4.25, height = 4)
ggsave(paste0(fn, ".png"), width = 4.25, height = 4)


# plot marker genes for 5-HT neurons
set.seed(123)
ggplot(df_5HT) + 
  geom_boxplot(aes(x = gene, y = mean, color = regions, fill = regions), 
               alpha = 0.5, outlier.shape = NA) + 
  geom_point(aes(x = gene, y = mean, color = regions), 
             position = position_jitterdodge()) + 
  scale_color_manual(values = pal_5HT, name = "annotation") + 
  scale_fill_manual(values = pal_5HT, name = "annotation") + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Enrichment") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(face = "italic"))

fn <- here(dir_plots, "enrichment_annotatedRegions_5HTmarkers")
ggsave(paste0(fn, ".pdf"), width = 4.25, height = 4)
ggsave(paste0(fn, ".png"), width = 4.25, height = 4)


# ---------------------
# plots with sample IDs
# ---------------------

# plots with shapes for sample IDs

# plot marker genes for NE neurons
set.seed(123)
ggplot(df_NE) + 
  geom_boxplot(aes(x = gene, y = mean, color = regions, fill = regions), 
               alpha = 0.5, outlier.shape = NA) + 
  geom_point(aes(x = gene, y = mean, color = regions, shape = sample), 
             stroke = 0.75, position = position_jitterdodge()) + 
  scale_color_manual(values = pal_NE, name = "annotation") + 
  scale_fill_manual(values = pal_NE, name = "annotation") + 
  scale_shape_manual(values = 1:12) + 
  guides(color = guide_legend(order = 1)) + 
  guides(fill = guide_legend(order = 1)) + 
  guides(shape = guide_legend(order = 2)) + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Enrichment") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(face = "italic"))

fn <- here(dir_plots, "enrichment_annotatedRegions_NEmarkers_withSampleIDs")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 4.5)


# plot marker genes for 5-HT neurons
set.seed(123)
ggplot(df_5HT) + 
  geom_boxplot(aes(x = gene, y = mean, color = regions, fill = regions), 
               alpha = 0.5, outlier.shape = NA) + 
  geom_point(aes(x = gene, y = mean, color = regions, shape = sample), 
             stroke = 0.75, position = position_jitterdodge()) + 
  scale_color_manual(values = pal_5HT, name = "annotation") + 
  scale_fill_manual(values = pal_5HT, name = "annotation") + 
  scale_shape_manual(values = 1:12) + 
  guides(color = guide_legend(order = 1)) + 
  guides(fill = guide_legend(order = 1)) + 
  guides(shape = guide_legend(order = 2)) + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Enrichment") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(face = "italic"))

fn <- here(dir_plots, "enrichment_annotatedRegions_5HTmarkers_withSampleIDs")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 4.5)

