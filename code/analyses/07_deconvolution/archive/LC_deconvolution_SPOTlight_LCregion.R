#####################################
# LC project
# Script for deconvolution of ST data
# Lukas Weber, Apr 2022
#####################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(SingleCellExperiment)
library(here)
library(SPOTlight)
library(scater)
library(scran)
library(ggplot2)

# directory to save plots
dir_plots <- here("plots", "07_deconvolution")


# ------------
# load ST data
# ------------

# load saved SPE object from previous script

fn_spe_LC <- here("processed_data", "SPE", "LC_batchCorrected_LCregions.rds")
spe_LC <- readRDS(fn_spe_LC)


# use LC region only

spe <- spe_LC

# convert sample IDs to factor
sample_ids <- c(
  "Br6522_LC_1_round1", "Br6522_LC_2_round1", 
  "Br8153_LC_round2", "Br5459_LC_round2", "Br2701_LC_round2", 
  "Br6522_LC_round3", "Br8079_LC_round3", "Br2701_LC_round3", "Br8153_LC_round3"
)
colData(spe)$sample_id <- factor(colData(spe)$sample_id, levels = sample_ids)


# select a single sample
spe <- spe[, colData(spe)$sample_id == "Br6522_LC_1_round1"]


# -------------------
# load snRNA-seq data
# -------------------

# SCE object containing merged clusters

fn_sce <- here("processed_data", "SCE", "sce_updated_LC.rda")
load(fn_sce)

dim(sce.lc)
table(colData(sce.lc)$cellType)
table(colData(sce.lc)$mergedCluster.HC)
table(colData(sce.lc)$cellType.collapsed)

# keep high-confidence clusters only, i.e. remove "ambig.lowNTx", "Neuron.ambig", "Neuron.mixed"
ix_keep <- !(colData(sce.lc)$cellType.collapsed %in% c("ambig.lowNTx", "Neuron.ambig", "Neuron.mixed"))
sce <- sce.lc[, ix_keep]

# merge 5-HT clusters ("Neuron.5HT", "Neuron.5HT_noDDC") into a single cluster
ix_5HT_all <- colData(sce)$cellType.collapsed %in% c("Neuron.5HT", "Neuron.5HT_noDDC")
colData(sce)$cellType.collapsed[ix_5HT_all] <- "Neuron.5HT"


# remove 5-HT cluster for deconvolution in LC region
sce <- sce[, !ix_5HT_all]


# remove empty factor levels
colData(sce)$cellType.collapsed <- droplevels(colData(sce)$cellType.collapsed)
table(colData(sce)$cellType.collapsed)

dim(sce)


# -------------
# run SPOTlight
# -------------

# using code from SPOTlight vignette
# https://marcelosua.github.io/SPOTlight/articles/SPOTlight_kidney.html


# feature selection: HVGs
dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# get the top 3000 genes
hvg <- getTopHVGs(dec, n = 3000)
head(hvg)

# add cluster labels to SCE object
colLabels(sce) <- colData(sce)$cellType.collapsed
# get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))

# compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # filter and keep relevant marker genes, those with AUC > 0.7
  x <- x[x$mean.AUC > 0.7, ]
  # sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
head(mgs_df, 2)

dim(sce)

# cell downsampling
# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$cellType.collapsed)
# downsample to at most 36 per identity & subset
n_cells <- 36
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

table(colData(sce)$cellType.collapsed)
dim(sce)

# deconvolution
res <- SPOTlight(
  x = sce,
  y = spe,
  groups = sce$cellType.collapsed,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

# extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]

# extract NMF model fit
mod <- res$NMF


# -------------
# visualization
# -------------

# plot scatterpies

ct <- colnames(mat)
mat[mat < 0.1] <- 0

# define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

# add cluster labels for plotting
lvs <- unique(mgs_df$cluster)

plotSpatialScatterpie(
  x = spe,
  y = mat,
  cell_types = colnames(mat),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.6) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal), 
    labels = lvs)

