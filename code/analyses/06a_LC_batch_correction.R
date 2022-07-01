#############################################
# LC analyses: batch correction using Harmony
# Lukas Weber, Jun 2022
#############################################

# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)
library(harmony)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "06a_batch_correction")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_logcounts.rds")
spe <- readRDS(fn_spe)

dim(spe)

# remove samples where NE neurons were not captured
samples_remove <- "Br5459_LC_round2"
spe <- spe[, !(colData(spe)$sample_id %in% samples_remove)]
colData(spe)$sample_id <- droplevels(colData(spe)$sample_id)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# -------------------------------------------------
# dimensionality reduction without batch correction
# -------------------------------------------------

# run dimensionality reduction without batch correction
# to demonstrate that batch correction is required in this dataset


# using code from OSCA

# feature selection
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]
spe <- spe[!is_mito, ]
dec <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
head(top_hvgs)

# dimensionality reduction
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)


# ---------
# PCA plots
# ---------

df <- cbind.data.frame(
  colData(spe), 
  spatialCoords(spe), 
  reducedDim(spe, "PCA"), reducedDim(spe, "UMAP")
)

pal <- c(
  unname(palette.colors(n = 8, palette = "Okabe-Ito")), 
  unname(palette.colors(n = 36, palette = "Polychrome 36"))
)


# donor ID
ggplot(df, aes(x = PC1, y = PC2, color = donor_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "noBatchCorrection_PCA_donorID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# round ID
ggplot(df, aes(x = PC1, y = PC2, color = round_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "noBatchCorrection_PCA_roundID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# sample ID
ggplot(df, aes(x = PC1, y = PC2, color = sample_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "noBatchCorrection_PCA_sampleID")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 5)
ggsave(paste0(fn, ".png"), width = 6.75, height = 5)


# sample and part ID
ggplot(df, aes(x = PC1, y = PC2, color = sample_part_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "noBatchCorrection_PCA_samplePartID")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5)
ggsave(paste0(fn, ".png"), width = 7, height = 5)


# ----------
# UMAP plots
# ----------

# donor ID
ggplot(df, aes(x = UMAP1, y = UMAP2, color = donor_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "noBatchCorrection_UMAP_donorID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# round ID
ggplot(df, aes(x = UMAP1, y = UMAP2, color = round_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "noBatchCorrection_UMAP_roundID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# sample ID
ggplot(df, aes(x = UMAP1, y = UMAP2, color = sample_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "noBatchCorrection_UMAP_sampleID")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 5)
ggsave(paste0(fn, ".png"), width = 6.75, height = 5)


# sample and part ID
ggplot(df, aes(x = UMAP1, y = UMAP2, color = sample_part_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("No batch correction") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "noBatchCorrection_UMAP_samplePartID")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5)
ggsave(paste0(fn, ".png"), width = 7, height = 5)


# ----------------------------------
# run batch correction using Harmony
# ----------------------------------

# run Harmony on PCA dimensions to integrate sample IDs
# note: integrating on sample IDs based on plots above

pca_matrix <- reducedDim(spe, "PCA")
sample_ids <- colData(spe)$sample_id
stopifnot(nrow(pca_matrix) == length(sample_ids))

set.seed(123)
harmony_embeddings <- HarmonyMatrix(
  pca_matrix, 
  meta_data = sample_ids, 
  do_pca = FALSE
)

colnames(harmony_embeddings) <- paste0("HARM", seq_len(ncol(harmony_embeddings)))

dim(harmony_embeddings)
head(harmony_embeddings, 2)

# store in SPE object
reducedDims(spe)[["HARM"]] <- harmony_embeddings

reducedDims(spe)


# -------------------------------
# PCA plots in Harmony dimensions
# -------------------------------

df <- cbind.data.frame(
  colData(spe), 
  spatialCoords(spe), 
  reducedDim(spe, "PCA"), 
  reducedDim(spe, "UMAP"), 
  reducedDim(spe, "HARM")
)


# donor ID
ggplot(df, aes(x = HARM1, y = HARM2, color = donor_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "batchCorrectedHarmony_donorID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# round ID
ggplot(df, aes(x = HARM1, y = HARM2, color = round_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "batchCorrectedHarmony_roundID")
ggsave(paste0(fn, ".pdf"), width = 6.25, height = 5)
ggsave(paste0(fn, ".png"), width = 6.25, height = 5)


# sample ID
ggplot(df, aes(x = HARM1, y = HARM2, color = sample_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "batchCorrectedHarmony_sampleID")
ggsave(paste0(fn, ".pdf"), width = 6.75, height = 5)
ggsave(paste0(fn, ".png"), width = 6.75, height = 5)


# sample and part ID
ggplot(df, aes(x = HARM1, y = HARM2, color = sample_part_id)) + 
  geom_point(size = 0.2, alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  ggtitle("Harmony embeddings") + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "batchCorrectedHarmony_samplePartID")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5)
ggsave(paste0(fn, ".png"), width = 7, height = 5)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LC_batchCorrectedHarmony")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

