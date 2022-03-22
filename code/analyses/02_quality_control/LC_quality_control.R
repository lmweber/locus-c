######################################################
# LC project
# Script for quality control (QC) and batch correction
# Lukas Weber, Mar 2022
######################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(harmony)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_Shiny.rds")
spe <- readRDS(fn_spe)

dim(spe)


# -------------------
# check preprocessing
# -------------------

# convert sample IDs to factor for easier plotting
sample_ids <- c(
  "Br6522_LC_1_round1", "Br6522_LC_2_round1", 
  "Br8153_LC_round2", "Br5459_LC_round2", "Br2701_LC_round2", 
  "Br6522_LC_round3", "Br8079_LC_round3", "Br2701_LC_round3", "Br8153_LC_round3"
)
colData(spe)$sample_id <- factor(colData(spe)$sample_id, levels = sample_ids)


# check SPE object contains only spots over tissue
table(colData(spe)$in_tissue)
all(colData(spe)$in_tissue)

dim(spe)


# --------------------------------------------------------
# split into 2 SPE objects for LC regions and white matter
# --------------------------------------------------------

# spots in manually annotated LC regions
table(colData(spe)$annot_region)

spe_LC <- spe[, colData(spe)$annot_region]
spe_WM <- spe[, !colData(spe)$annot_region]

dim(spe_LC)
dim(spe_WM)


# ------------------------
# spot-level QC: all spots
# ------------------------

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]

# calculate QC metrics using scater package
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# note duplicate columns from previously
all(colData(spe)$sum == colData(spe)$sum_umi)
all(colData(spe)$detected == colData(spe)$sum_gene)
all(colData(spe)$subsets_mito_sum == colData(spe)$expr_chrM)
all(colData(spe)$subsets_mito_percent == colData(spe)$expr_chrM_ratio * 100)

# note wide range of values in combined object
summary(colData(spe)$sum)
summary(colData(spe)$detected)


# plot histograms of QC metrics
fn <- file.path(dir_plots, "01_quality_control", "QC_histograms_allSpots.pdf")
pdf(fn, width = 8, height = 2.5)
par(mfrow = c(1, 4))
hist(colData(spe)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe)$detected, xlab = "detected", main = "Genes per spot")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito UMIs")
hist(colData(spe)$cell_count, xlab = "no. cells", main = "No. cells per spot")
par(mfrow = c(1, 1))
dev.off()


# -------------------------
# spot-level QC: LC regions
# -------------------------

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_LC)$gene_name)
table(is_mito)
rowData(spe_LC)$gene_name[is_mito]

# calculate QC metrics using scater package
spe_LC <- addPerCellQC(spe_LC, subsets = list(mito = is_mito))

summary(colData(spe_LC)$sum)
summary(colData(spe_LC)$detected)


# plot histograms of QC metrics
fn <- file.path(dir_plots, "01_quality_control", "QC_histograms_LCregions.pdf")
pdf(fn, width = 8, height = 2.5)
par(mfrow = c(1, 4))
hist(colData(spe_LC)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe_LC)$detected, xlab = "detected", main = "Genes per spot")
hist(colData(spe_LC)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito UMIs")
hist(colData(spe_LC)$cell_count, xlab = "no. cells", main = "No. cells per spot")
par(mfrow = c(1, 1))
dev.off()


# select adaptive thresholds using 3 MAD methodology from OSCA
reasons <- perCellQCFilters(perCellQCMetrics(spe_LC, subsets = list(mito = is_mito)), 
                            sub.fields = "subsets_mito_percent")
colSums(as.matrix(reasons))

# check thresholds
thresh_low_lib_size <- attr(reasons$low_lib_size, "thresholds")["lower"]
thresh_low_lib_size  ## 37.1
thresh_low_n_features <- attr(reasons$low_n_features, "thresholds")["lower"]
thresh_low_n_features ## 25.2
thresh_high_subsets_mito_percent <- attr(reasons$high_subsets_mito_percent, "thresholds")["higher"]
thresh_high_subsets_mito_percent  ## 63.4

# store in SPE object
stopifnot(nrow(reasons) == nrow(colData(spe_LC)))
colData(spe_LC) <- cbind(colData(spe_LC), reasons)


# -------------------------
# spot-level QC: WM regions
# -------------------------

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_WM)$gene_name)
table(is_mito)
rowData(spe_WM)$gene_name[is_mito]

# calculate QC metrics using scater package
spe_WM <- addPerCellQC(spe_WM, subsets = list(mito = is_mito))

summary(colData(spe_WM)$sum)
summary(colData(spe_WM)$detected)


# plot histograms of QC metrics
fn <- file.path(dir_plots, "01_quality_control", "QC_histograms_WMregions.pdf")
pdf(fn, width = 8, height = 2.5)
par(mfrow = c(1, 4))
hist(colData(spe_WM)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe_WM)$detected, xlab = "detected", main = "Genes per spot")
hist(colData(spe_WM)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito UMIs")
hist(colData(spe_WM)$cell_count, xlab = "no. cells", main = "No. cells per spot")
par(mfrow = c(1, 1))
dev.off()


# --------------
# generate plots
# --------------

# plot QC metrics using code from ggspavis
df <- cbind.data.frame(colData(spe), spatialCoords(spe))

# sum UMIs
ggplot(df, aes(x = x, y = y, color = sum < 100)) + 
  facet_wrap(~ sample_id, nrow = 2) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "red")) + 
  ggtitle("Spot-level QC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
#ggsave(paste0(here("plots", "QC", "QC_samples_sum"), ".pdf"), width = 12, height = 6.5)
ggsave(paste0(here("plots", "QC", "QC_samples_sum"), ".png"), width = 12, height = 6.5)

# detected genes
ggplot(df, aes(x = x, y = y, color = detected < 100)) + 
  facet_wrap(~ sample_id, nrow = 2) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "red")) + 
  ggtitle("Spot-level QC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
#ggsave(paste0(here("plots", "QC", "QC_samples_detected"), ".pdf"), width = 12, height = 6.5)
ggsave(paste0(here("plots", "QC", "QC_samples_detected"), ".png"), width = 12, height = 6.5)

# cell count
ggplot(df, aes(x = x, y = y, color = cell_count > 10)) + 
  facet_wrap(~ sample_id, nrow = 2) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "red")) + 
  ggtitle("Spot-level QC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
#ggsave(paste0(here("plots", "QC", "QC_samples_count"), ".pdf"), width = 12, height = 6.5)
ggsave(paste0(here("plots", "QC", "QC_samples_count"), ".png"), width = 12, height = 6.5)

# expression of TH
ix_TH <- which(rowData(spe)$symbol == "TH")
df$TH_expr <- counts(spe)[ix_TH, ]

ggplot(df, aes(x = x, y = y, color = TH_expr)) + 
  facet_wrap(~ sample_id, nrow = 2) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "blue") + 
  ggtitle("TH expression") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
#ggsave(paste0(here("plots", "QC", "QC_samples_TH_expr"), ".pdf"), width = 12, height = 6.5)
ggsave(paste0(here("plots", "QC", "QC_samples_TH_expr"), ".png"), width = 12, height = 6.5)


# thresholding on TH
ggplot(df, aes(x = x, y = y, color = TH_expr > 0)) + 
  facet_wrap(~ sample_id, nrow = 2) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "red")) + 
  ggtitle("TH thresholding") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
#ggsave(paste0(here("plots", "QC", "QC_samples_TH_thresh_0"), ".pdf"), width = 12, height = 6.5)
ggsave(paste0(here("plots", "QC", "QC_samples_TH_thresh_0"), ".png"), width = 12, height = 6.5)

ggplot(df, aes(x = x, y = y, color = TH_expr > 1)) + 
  facet_wrap(~ sample_id, nrow = 2) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "red")) + 
  ggtitle("TH thresholding") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
#ggsave(paste0(here("plots", "QC", "QC_samples_TH_thresh_1"), ".pdf"), width = 12, height = 6.5)
ggsave(paste0(here("plots", "QC", "QC_samples_TH_thresh_1"), ".png"), width = 12, height = 6.5)

ggplot(df, aes(x = x, y = y, color = TH_expr > 2)) + 
  facet_wrap(~ sample_id, nrow = 2) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "red")) + 
  ggtitle("TH thresholding") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
#ggsave(paste0(here("plots", "QC", "QC_samples_TH_thresh_2"), ".pdf"), width = 12, height = 6.5)
ggsave(paste0(here("plots", "QC", "QC_samples_TH_thresh_2"), ".png"), width = 12, height = 6.5)


# -----------
# save object
# -----------

# keep all spots
colData(spe)$discard <- FALSE
colData(spe)

# store expression of TH in colData
colData(spe)$TH_sum <- counts(spe)[ix_TH, ]
colData(spe)

# save object
fn_out <- here("processed_data", "SPE", "LCrounds1to3_SPE_processed.rds")
saveRDS(spe, fn_out)


# -------------
# summary plots
# -------------

sample_ids <- unique(colData(spe)$sample_id)
sample_ids

for (s in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s] & colData(spe)$tissue == 1]
  df <- as.data.frame(cbind(colData(spe_sub), spatialCoords(spe_sub)))
  
  ix_TH <- which(toupper(rowData(spe_sub)$symbol) == "TH")
  stopifnot(length(ix_TH) == 1)
  df$TH <- counts(spe_sub)[ix_TH, ]
  
  ix_SNAP25 <- which(toupper(rowData(spe_sub)$symbol) == "SNAP25")
  stopifnot(length(ix_SNAP25) == 1)
  df$SNAP25 <- counts(spe_sub)[ix_SNAP25, ]
  
  
  # plot total UMI counts
  p <- ggplot(df, aes(x = y, y = x, color = sum)) + 
    geom_point(size = 0.75) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_gradient(low = "gray90", high = "red") + 
    ggtitle(sample_ids[s]) + 
    labs(color = "UMI counts") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  if (!dir.exists(here("plots", "summary", sample_ids[s]))) {
    dir.create(here("plots", "summary", sample_ids[s]), recursive = TRUE)
  }
  
  fn <- here("plots", "summary", sample_ids[s], 
             paste0(sample_ids[s], "_sum"))
  #ggsave(paste0(fn, ".pdf"), plot = p, width = 5.25, height = 4)
  ggsave(paste0(fn, ".png"), plot = p, width = 5.25, height = 4)
  
  
  # plot detected genes
  p <- ggplot(df, aes(x = y, y = x, color = detected)) + 
    geom_point(size = 0.75) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_gradient(low = "gray90", high = "darkorchid4") + 
    ggtitle(sample_ids[s]) + 
    labs(color = "detected genes") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- here("plots", "summary", sample_ids[s], 
             paste0(sample_ids[s], "_detected"))
  #ggsave(paste0(fn, ".pdf"), plot = p, width = 5.25, height = 4)
  ggsave(paste0(fn, ".png"), plot = p, width = 5.25, height = 4)
  
  
  # plot expression of TH
  p <- ggplot(df, aes(x = y, y = x, color = TH)) + 
    geom_point(size = 0.75) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_gradient(low = "gray90", high = "darkgreen") + 
    ggtitle(sample_ids[s]) + 
    labs(color = "TH") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- here("plots", "summary", sample_ids[s], 
             paste0(sample_ids[s], "_TH"))
  #ggsave(paste0(fn, ".pdf"), plot = p, width = 5.25, height = 4)
  ggsave(paste0(fn, ".png"), plot = p, width = 5.25, height = 4)
  
  
  # plot expression of SNAP25
  p <- ggplot(df, aes(x = y, y = x, color = SNAP25)) + 
    geom_point(size = 0.75) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_gradient(low = "gray90", high = "darkgreen") + 
    ggtitle(sample_ids[s]) + 
    labs(color = "SNAP25") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- here("plots", "summary", sample_ids[s], 
             paste0(sample_ids[s], "_SNAP25"))
  #ggsave(paste0(fn, ".pdf"), plot = p, width = 5.25, height = 4)
  ggsave(paste0(fn, ".png"), plot = p, width = 5.25, height = 4)
  
}


# ----------------------------------------
# QC / summary plots for manual annotation
# ----------------------------------------

df <- cbind.data.frame(colData(spe), spatialCoords(spe))

# convert NA to FALSE
df$annot_region[is.na(df$annot_region)] <- FALSE
df$annot_spot[is.na(df$annot_spot)] <- FALSE


# lasso regions
ggplot(df, aes(x = x, y = y, color = annot_region)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "red")) + 
  ggtitle("Manual annotation: lasso regions") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- here("plots", "annot", "annot_region")
ggsave(paste0(fn, ".png"), width = 12, height = 5.25)


# individual spots
ggplot(df, aes(x = x, y = y, color = annot_spot)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "red")) + 
  ggtitle("Manual annotation: individual spots") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- here("plots", "annot", "annot_spot")
ggsave(paste0(fn, ".png"), width = 12, height = 5.25)


# number of cells per spot: all spots

ggplot(df, aes(x = cell_count)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_histogram(color = "black", fill = "navy") + 
  ggtitle("Number of cells per spot: all spots") + 
  theme_bw()

fn <- here("plots", "annot", "ncells_all")
ggsave(paste0(fn, ".png"), width = 11.5, height = 5.25)


# number of cells per spot: manually annotated regions

# skip round 2 since spot counting failed for this round
df_rounds13 <- df[df$round_id %in% c("round1", "round3"), ]

df_sub <- df_rounds13[df_rounds13$annot_region == TRUE, ]

ggplot(df_sub, aes(x = cell_count)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_histogram(color = "black", fill = "navy") + 
  ggtitle("Number of cells per spot: manually annotated regions") + 
  theme_bw()

fn <- here("plots", "annot", "ncells_annot_region")
ggsave(paste0(fn, ".png"), width = 7, height = 5.25)


# number of cells per spot: manually annotated spots

df_sub <- df_rounds13[df_rounds13$annot_spot == TRUE, ]

ggplot(df_sub, aes(x = cell_count)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_histogram(color = "black", fill = "navy") + 
  ggtitle("Number of cells per spot: manually annotated individual spots") + 
  theme_bw()

fn <- here("plots", "annot", "ncells_annot_spot")
ggsave(paste0(fn, ".png"), width = 7, height = 5.25)


# --------------------------------------------
# summary plots for NE neurons and 5HT neurons
# --------------------------------------------

# expression plots for some key markers identifying NE neurons and 5HT neurons

# using logcounts
fn_spe <- here("processed_data", "SPE", "LCrounds1to3_SPE_shiny.rds")
spe <- readRDS(fn_spe)

colData(spe)$sample_id <- factor(colData(spe)$sample_id, levels = sample_ids)

ix_TH <- which(toupper(rowData(spe)$gene_name) == "TH")
ix_SLC6A4 <- which(toupper(rowData(spe)$gene_name) == "SLC6A4")
ix_TPH2 <- which(toupper(rowData(spe)$gene_name) == "TPH2")

colData(spe)$logcounts_TH <- logcounts(spe)[ix_TH, ]
colData(spe)$logcounts_SLC6A4 <- logcounts(spe)[ix_SLC6A4, ]
colData(spe)$logcounts_TPH2 <- logcounts(spe)[ix_TPH2, ]

df <- cbind.data.frame(colData(spe), spatialCoords(spe))

# convert NA to FALSE
df$annot_region[is.na(df$annot_region)] <- FALSE
df$annot_spot[is.na(df$annot_spot)] <- FALSE


# expression plots

ggplot(df, aes(x = x, y = y, color = logcounts_TH)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "darkgreen") + 
  ggtitle("Expression of TH") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- here("plots", "annot", "logcounts_TH")
ggsave(paste0(fn, ".png"), width = 12, height = 5.25)


ggplot(df, aes(x = x, y = y, color = logcounts_TH)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "red") + 
  geom_point(data = df[df$annot_spot == TRUE, ], 
             pch = 21, size = 0.5, alpha = 0.2, color = "black") + 
  ggtitle("Expression of TH / manually annotated spots") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- here("plots", "annot", "logcounts_TH_overlap_spots")
ggsave(paste0(fn, ".png"), width = 12, height = 5.25)


ggplot(df, aes(x = x, y = y, color = logcounts_TH)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "red") + 
  geom_point(data = df[df$annot_region == TRUE, ], 
             pch = 21, size = 0.35, alpha = 0.15, color = "black") + 
  ggtitle("Expression of TH / manually annotated regions") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- here("plots", "annot", "logcounts_TH_overlap_regions")
ggsave(paste0(fn, ".png"), width = 12, height = 5.25)


ggplot(df, aes(x = x, y = y, color = logcounts_SLC6A4)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "red") + 
  ggtitle("Expression of SLC6A4") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- here("plots", "annot", "logcounts_SLC6A4")
ggsave(paste0(fn, ".png"), width = 12, height = 5.25)


ggplot(df, aes(x = x, y = y, color = logcounts_TPH2)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "red") + 
  ggtitle("Expression of TPH2") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- here("plots", "annot", "logcounts_TPH2")
ggsave(paste0(fn, ".png"), width = 12, height = 5.25)

