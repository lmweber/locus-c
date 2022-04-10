#################################
# LC project
# Script for quality control (QC)
# Lukas Weber, Apr 2022
#################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
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


# -------------------------------
# calculate spot-level QC metrics
# -------------------------------

# combined object containing both LC and WM regions
# calculate QC metrics; then split object to select QC thresholds in LC and WM

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]

# calculate QC metrics using scater package
spe <- addPerCellQCMetrics(spe, subsets = list(mito = is_mito))

# note: remove duplicate columns from previously
all(colData(spe)$sum == colData(spe)$sum_umi)
all(colData(spe)$detected == colData(spe)$sum_gene)
all(colData(spe)$subsets_mito_sum == colData(spe)$expr_chrM)
all(colData(spe)$subsets_mito_percent == colData(spe)$expr_chrM_ratio * 100)
ix_remove <- which(names(colData(spe)) %in% c("sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio"))
ix_remove
dim(colData(spe))
colData(spe) <- colData(spe)[, -ix_remove]
dim(colData(spe))


# note extremely wide range of values in combined object containing both LC and WM
summary(colData(spe)$sum)
summary(colData(spe)$detected)


# plot histograms of QC metrics
fn <- file.path(dir_plots, "01_quality_control", "QC_histograms_allSpots.pdf")
pdf(fn, width = 8, height = 2.5)
par(mfrow = c(1, 4))
hist(reasons$sum, xlab = "sum", main = "UMIs per spot")
hist(reasons$detected, xlab = "detected", main = "Genes per spot")
hist(reasons$subsets_mito_percent, xlab = "percent mito", main = "Percent mito UMIs")
hist(colData(spe)$cell_count, xlab = "no. cells", main = "No. cells per spot")
par(mfrow = c(1, 1))
dev.off()


# --------------------------------------
# split into 2 SPE objects for LC and WM
# --------------------------------------

# select QC thresholds separately in LC and WM regions

# spots in manually annotated LC regions
table(colData(spe)$annot_region)

spe_LC <- spe[, colData(spe)$annot_region]
spe_WM <- spe[, !colData(spe)$annot_region]

dim(spe_LC)
dim(spe_WM)


# compare summary values

summary(colData(spe)$sum)
summary(colData(spe)$detected)

summary(colData(spe_LC)$sum)
summary(colData(spe_LC)$detected)

summary(colData(spe_WM)$sum)
summary(colData(spe_WM)$detected)


# to do: plots comparing these summary values


# ---------------------------
# spot-level QC in LC regions
# ---------------------------

# check
dim(spe_LC)
table(colData(spe_LC)$sample_id)
table(colData(spe_LC)$sample_part_id)


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
reasons_LC <- perCellQCFilters(perCellQCMetrics(spe_LC, subsets = list(mito = is_mito)), 
                            sub.fields = "subsets_mito_percent")
colSums(as.matrix(reasons_LC))

# check thresholds
thresh_low_lib_size <- attr(reasons_LC$low_lib_size, "thresholds")["lower"]
thresh_low_lib_size  ## 37.1
thresh_low_n_features <- attr(reasons_LC$low_n_features, "thresholds")["lower"]
thresh_low_n_features ## 25.2
thresh_high_subsets_mito_percent <- attr(reasons_LC$high_subsets_mito_percent, "thresholds")["higher"]
thresh_high_subsets_mito_percent  ## 63.4

# store in subsetted SPE object
stopifnot(nrow(reasons_LC) == nrow(colData(spe_LC)))
colData(spe_LC) <- cbind(colData(spe_LC), reasons_LC)


# store in main SPE object
colData(spe)$low_lib_size <- FALSE
colData(spe)$low_n_features <- FALSE
colData(spe)$high_subsets_mito_percent <- FALSE
colData(spe)$discard <- FALSE
colData(spe)[colData(spe_LC)$key_id, "low_lib_size"] <- colData(spe_LC)$low_lib_size
colData(spe)[colData(spe_LC)$key_id, "low_n_features"] <- colData(spe_LC)$low_n_features
colData(spe)[colData(spe_LC)$key_id, "high_subsets_mito_percent"] <- colData(spe_LC)$high_subsets_mito_percent
colData(spe)[colData(spe_LC)$key_id, "discard"] <- colData(spe_LC)$discard

# check
rbind(table(colData(spe_LC)$low_lib_size), table(colData(spe)$low_lib_size))
rbind(table(colData(spe_LC)$low_n_features), table(colData(spe)$low_n_features))
rbind(table(colData(spe_LC)$high_subsets_mito_percent), table(colData(spe)$high_subsets_mito_percent))
rbind(table(colData(spe_LC)$discard), table(colData(spe)$discard))


# ---------------------------
# spot-level QC in WM regions
# ---------------------------

# check
dim(spe_WM)
table(colData(spe_WM)$sample_id)


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
reasons_WM <- perCellQCFilters(perCellQCMetrics(spe_WM, subsets = list(mito = is_mito)), 
                            sub.fields = "subsets_mito_percent")
colSums(as.matrix(reasons_WM))

# check thresholds
thresh_low_lib_size <- attr(reasons_WM$low_lib_size, "thresholds")["lower"]
thresh_low_lib_size  ## 1.3
thresh_low_n_features <- attr(reasons_WM$low_n_features, "thresholds")["lower"]
thresh_low_n_features ## 1.4
thresh_high_subsets_mito_percent <- attr(reasons_WM$high_subsets_mito_percent, "thresholds")["higher"]
thresh_high_subsets_mito_percent  ## 69.6


# note; adaptive thresholds are too low in this dataset for lib_size and n_features
quantile(colData(spe_WM)$sum, seq(0, 1, by = 0.1))
quantile(colData(spe_WM)$detected, seq(0, 1, by = 0.1))

# cut at arbitrary low value instead
low_lib_size <- colData(spe_WM)$sum < 10
low_n_features <- colData(spe_WM)$detected < 10
table(low_lib_size)
table(low_n_features)

high_subsets_mito_percent <- reasons_WM$high_subsets_mito_percent
table(high_subsets_mito_percent)

discard <- low_lib_size | low_n_features | high_subsets_mito_percent
table(discard)


# store in subsetted SPE object
colData(spe_WM) <- cbind(
  colData(spe_WM), 
  low_lib_size = low_lib_size, 
  low_n_features = low_n_features, 
  high_subsets_mito_percent = high_subsets_mito_percent, 
  discard = discard
)

# store in main SPE object (in correct rows of colData)
colData(spe)[colData(spe_WM)$key_id, "low_lib_size"] <- colData(spe_WM)$low_lib_size
colData(spe)[colData(spe_WM)$key_id, "low_n_features"] <- colData(spe_WM)$low_n_features
colData(spe)[colData(spe_WM)$key_id, "high_subsets_mito_percent"] <- colData(spe_WM)$high_subsets_mito_percent
colData(spe)[colData(spe_WM)$key_id, "discard"] <- colData(spe_WM)$discard

# check
rbind(table(colData(spe_LC)$low_lib_size), 
      table(colData(spe_WM)$low_lib_size), 
      table(colData(spe)$low_lib_size))
rbind(table(colData(spe_LC)$low_n_features), 
      table(colData(spe_WM)$low_n_features), 
      table(colData(spe)$low_n_features))
rbind(table(colData(spe_LC)$high_subsets_mito_percent), 
      table(colData(spe_WM)$high_subsets_mito_percent), 
      table(colData(spe)$high_subsets_mito_percent))
rbind(table(colData(spe_LC)$discard), 
      table(colData(spe_WM)$discard), 
      table(colData(spe)$discard))


# --------------------
# plot discarded spots
# --------------------

# summary of discarded spots in LC and WM regions
table(colData(spe_LC)$sample_part_id, colData(spe_LC)$discard)
table(colData(spe_WM)$sample_id, colData(spe_WM)$discard)

# to do: summary plot of these values


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

fn_out <- here("processed_data", "SPE", "LC_qualityControlled")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

