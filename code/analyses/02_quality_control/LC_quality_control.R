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
library(ggnewscale)


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


# ------------------------------------------------------
# split into 2 SPE objects for LC regions and WM regions
# ------------------------------------------------------

# spots in manually annotated LC regions
table(colData(spe)$annot_region)

spe_LC <- spe[, colData(spe)$annot_region]
spe_WM <- spe[, !colData(spe)$annot_region]

dim(spe_LC)
dim(spe_WM)


# -----------------------------------------------------
# spot-level QC: all spots (LC and WM regions together)
# -----------------------------------------------------

# not used for final QC

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

# check number of spots
dim(spe_LC)
table(colData(spe_LC)$sample_id)
table(colData(spe_LC)$sample_part_id)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_LC)$gene_name)
table(is_mito)
rowData(spe_LC)$gene_name[is_mito]

# calculate QC metrics using scater package
spe_LC <- addPerCellQC(spe_LC, subsets = list(mito = is_mito))


# check summaries and compare with combined object
summary(colData(spe_LC)$sum)
summary(colData(spe)$sum)

summary(colData(spe_LC)$detected)
summary(colData(spe)$detected)


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

# store in subsetted SPE object
stopifnot(nrow(reasons) == nrow(colData(spe_LC)))
colData(spe_LC) <- cbind(colData(spe_LC), reasons)

# store in main SPE object
colData(spe)$discard <- FALSE
colData(spe)[colData(spe_LC)$key_id, "discard"] <- colData(spe_LC)$discard

# check
table(colData(spe_LC)$discard)
table(colData(spe)$discard)


# -------------------------
# spot-level QC: WM regions
# -------------------------

# check number of spots
dim(spe_WM)
table(colData(spe_WM)$sample_id)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe_WM)$gene_name)
table(is_mito)
rowData(spe_WM)$gene_name[is_mito]

# calculate QC metrics using scater package
spe_WM <- addPerCellQC(spe_WM, subsets = list(mito = is_mito))


# check summaries and compare with combined object
summary(colData(spe_WM)$sum)
summary(colData(spe_LC)$sum)
summary(colData(spe)$sum)

summary(colData(spe_WM)$detected)
summary(colData(spe_LC)$detected)
summary(colData(spe)$detected)


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


# select adaptive thresholds using 3 MAD methodology from OSCA
reasons <- perCellQCFilters(perCellQCMetrics(spe_WM, subsets = list(mito = is_mito)), 
                            sub.fields = "subsets_mito_percent")
colSums(as.matrix(reasons))

# check thresholds
thresh_low_lib_size <- attr(reasons$low_lib_size, "thresholds")["lower"]
thresh_low_lib_size  ## 1.3
thresh_low_n_features <- attr(reasons$low_n_features, "thresholds")["lower"]
thresh_low_n_features ## 1.4
thresh_high_subsets_mito_percent <- attr(reasons$high_subsets_mito_percent, "thresholds")["higher"]
thresh_high_subsets_mito_percent  ## 69.6

# store in subsetted SPE object
stopifnot(nrow(reasons) == nrow(colData(spe_WM)))
colData(spe_WM) <- cbind(colData(spe_WM), reasons)

# store in main SPE object
colData(spe)[colData(spe_WM)$key_id, "discard"] <- colData(spe_WM)$discard

# check
table(colData(spe_WM)$discard)
table(colData(spe_LC)$discard)
table(colData(spe)$discard)


# --------------------
# plot discarded spots
# --------------------

# summary of discarded spots in LC and WM regions
table(colData(spe_LC)$sample_part_id, colData(spe_LC)$discard)
table(colData(spe_WM)$sample_id, colData(spe_WM)$discard)

# plot all spots
df <- cbind.data.frame(colData(spe), spatialCoords(spe))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85") + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70") + 
  new_scale_color() + 
  geom_point(data = df[df$discard, , drop = FALSE], 
             aes(color = discard), size = 0.1) + 
  scale_color_manual(values = "red") + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_discard_combinedSampleMetrics")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LC_preprocessed_QC")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

