##################################################################
# LC project
# Script for downstream analyses: BayesSpace evaluations and plots
# Lukas Weber, May 2022
##################################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(mclust)
library(dplyr)
library(tidyr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "06_downstream", "BayesSpace")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_BayesSpace.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)


sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# -------------
# plot clusters
# -------------

df <- cbind.data.frame(
  colData(spe), spatialCoords(spe), 
  reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"), reducedDim(spe, "HARM")
)

df$spatial.cluster <- as.factor(df$spatial.cluster)

pal <- unname(palette.colors(8, palette = "Okabe-Ito"))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = spatial.cluster)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = pal) + 
  scale_y_reverse() + 
  ggtitle("BayesSpace clustering") + 
  labs(color = "cluster") + 
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "BayesSpace_clustering")
ggsave(paste0(fn, ".pdf"), width = 7.25, height = 4)
ggsave(paste0(fn, ".png"), width = 7.25, height = 4)


# -----------
# evaluations
# -----------

# evaluate performance of BayesSpace for clustering LC regions
# note: in this dataset we are only interested in the LC regions; ignore other
# clusters in the WM regions, which we assume are mostly noise


# select cluster ID matching most closely to LC regions
selected <- 5

# calculate adjusted Rand index (ARI)
clus <- as.numeric(colData(spe)$spatial.cluster == selected)
ref <- as.numeric(colData(spe)$annot_region)

# overall comparison
table(clustering = clus, 
      reference = ref)

ari_overall <- adjustedRandIndex(clus, ref)
ari_overall


# calculate precision, recall, and ARI per sample
precision <- recall <- f1 <- ari <- rep(NA, length(sample_ids))
names(precision) <- names(recall) <- names(f1) <- names(ari) <- sample_ids

for (i in seq_along(sample_ids)) {
  # select clustering and reference labels for sample i
  ix <- colData(spe)$sample_id == sample_ids[i]
  clus_i <- clus[ix]
  ref_i <- ref[ix]
  stopifnot(length(clus_i) == sum(ix))
  stopifnot(length(ref_i) == sum(ix))
  # calculate TP, FP, FN, TN, n
  true_positives <- sum(clus_i == 1 & ref_i == 1)
  false_positives <- sum(clus_i == 1 & ref_i == 0)
  false_negatives <- sum(clus_i == 0 & ref_i == 1)
  true_negatives <- sum(clus_i == 0 & ref_i == 0)
  n <- length(ref_i)
  # calculate precision, recall, ARI
  precision_i <- true_positives / (true_positives + false_positives)
  recall_i <- true_positives / (true_positives + false_negatives)
  f1_i <- 2 * (precision_i * recall_i) / (precision_i + recall_i)
  ari_i <- adjustedRandIndex(clus_i, ref_i)
  # store results
  precision[i] <- precision_i
  recall[i] <- recall_i
  f1[i] <- f1_i
  ari[i] <- ari_i
}

# summary table
df_summary <- data.frame(
  precision = precision, 
  recall = recall, 
  F1 = f1, 
  ARI = ari
)

df_summary


# ------------
# plot summary
# ------------

# clustering performance by sample

df <- df_summary %>% 
  mutate(sample_id = factor(rownames(df_summary), levels = rownames(df_summary))) %>% 
  pivot_longer(cols = -sample_id, names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = c("precision", "recall", "F1", "ARI")))

pal <- unname(palette.colors(4, palette = "Classic Tableau"))

set.seed(1)
ggplot(df, aes(x = metric, y = value, color = metric)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = sample_id), width = 0.15, size = 1.25, stroke = 0.75) + 
  scale_color_manual(values = pal) + 
  scale_shape_manual(values = 0:7) + 
  guides(color = guide_legend(order = 1)) + 
  guides(shape = guide_legend(order = 2)) + 
  ylim(c(0, 1)) + 
  ggtitle("BayesSpace clustering performance") + 
  theme_bw() + 
  theme(axis.title.y = element_blank())

fn <- file.path(dir_plots, "BayesSpace_clustering_performance")
ggsave(paste0(fn, ".pdf"), width = 5.75, height = 4.5)
ggsave(paste0(fn, ".png"), width = 5.75, height = 4.5)

