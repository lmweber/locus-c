#########################################
# LC analyses: enrichment of marker genes
# Lukas Weber, Jun 2022
#########################################

# module load conda_R/devel
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

fn_spe <- here("processed_data", "SPE", "LC_logcounts.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)


# ------------------------------------
# calculate enrichment of marker genes
# ------------------------------------

# enrichment defined in terms of mean logcounts per spot in LC vs. WM regions

# select marker genes of interest
# for NE neurons and 5-HT neurons
marker_genes <- c("TH", "SLC6A2", "TPH2", "SLC6A4")

# LC region annotations
table(colData(spe)$annot_region)


enrichment <- matrix(NA, nrow = length(sample_ids), ncol = length(marker_genes))
rownames(enrichment) <- sample_ids
colnames(enrichment) <- marker_genes

enrichment_LC <- enrichment_WM <- enrichment


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

# select annotated WM regions
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


# ---------------
# plot enrichment
# ---------------

df_enrichment_LC <- as.data.frame(enrichment_LC)
df_enrichment_LC$region <- "LC"
df_enrichment_LC$sample <- rownames(df_enrichment_LC)

df_enrichment_WM <- as.data.frame(enrichment_WM)
df_enrichment_WM$region <- "WM"
df_enrichment_WM$sample <- rownames(df_enrichment_WM)


df_enrichment_LC <- pivot_longer(df_enrichment_LC, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")
df_enrichment_WM <- pivot_longer(df_enrichment_WM, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")


df <- 
  full_join(df_enrichment_LC, df_enrichment_WM) %>% 
  mutate(regions = factor(region, 
                          levels = c("LC", "WM"), 
                          labels = c("LC regions", "WM regions"))) %>% 
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


pal <- c("darkmagenta", "gray30")


# plot marker genes for NE neurons
set.seed(123)
ggplot(df_NE, aes(x = gene, y = mean, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge()) + 
  scale_color_manual(values = pal, name = "annotation") + 
  scale_fill_manual(values = pal, name = "annotation") + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Enrichment") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(face = "italic"))

fn <- here(dir_plots, "enrichment_annotatedRegions_NEmarkers")
ggsave(paste0(fn, ".pdf"), width = 4, height = 4)
ggsave(paste0(fn, ".png"), width = 4, height = 4)


# plot marker genes for 5-HT neurons
set.seed(123)
ggplot(df_5HT, aes(x = gene, y = mean, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge()) + 
  scale_color_manual(values = pal, name = "annotation") + 
  scale_fill_manual(values = pal, name = "annotation") + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Enrichment") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(face = "italic"))

fn <- here(dir_plots, "enrichment_annotatedRegions_5HTmarkers")
ggsave(paste0(fn, ".pdf"), width = 4, height = 4)
ggsave(paste0(fn, ".png"), width = 4, height = 4)


# --------------------------------------
# plot excluding sample Br5459_LC_round2
# --------------------------------------

# excluding samples where NE neurons were not captured (Br5459_LC_round2)

df_NE_sub <- df_NE[!(df_NE$sample %in% "Br5459_LC_round2"), ]

# plot marker genes for NE neurons
set.seed(123)
ggplot(df_NE_sub, aes(x = gene, y = mean, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge()) + 
  scale_color_manual(values = pal, name = "annotation") + 
  scale_fill_manual(values = pal, name = "annotation") + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Enrichment") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        axis.text.x = element_text(face = "italic"))

fn <- here(dir_plots, "enrichment_annotatedRegions_NEmarkers_excBr5459")
ggsave(paste0(fn, ".pdf"), width = 4, height = 4)
ggsave(paste0(fn, ".png"), width = 4, height = 4)

