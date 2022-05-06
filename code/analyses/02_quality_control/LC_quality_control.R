#################################
# LC project
# Script for quality control (QC)
# Lukas Weber, May 2022
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
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)


# directory to save plots
dir_plots <- here("plots", "02_quality_control")


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
fn <- file.path(dir_plots, "QC_histograms_allSpots.pdf")
pdf(fn, width = 8, height = 2.5)
par(mfrow = c(1, 4))
hist(colData(spe)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe)$detected, xlab = "detected", main = "Genes per spot")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito UMIs")
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


# table of summary values

df_summary <- data.frame(
  sumUMI_combined = round(as.numeric(summary(colData(spe)$sum))), 
  sumUMI_LC = round(as.numeric(summary(colData(spe_LC)$sum))), 
  sumUMI_WM = round(as.numeric(summary(colData(spe_WM)$sum))), 
  detectedGenes_combined = round(as.numeric(summary(colData(spe)$detected))), 
  detectedGenes_LC = round(as.numeric(summary(colData(spe_LC)$detected))), 
  detectedGenes_WM = round(as.numeric(summary(colData(spe_WM)$detected)))
)
rownames(df_summary) <- names(summary(colData(spe)$sum))

df_summary


# ---------------------------
# spot-level QC in LC regions
# ---------------------------

# to do: use sample-specific QC thresholds

# check
dim(spe_LC)
table(colData(spe_LC)$sample_id)
table(colData(spe_LC)$sample_part_id)


# plot histograms of QC metrics
fn <- file.path(dir_plots, "QC_histograms_LCregions.pdf")
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

# to do: use sample-specific QC thresholds

# check
dim(spe_WM)
table(colData(spe_WM)$sample_id)


# plot histograms of QC metrics
fn <- file.path(dir_plots, "QC_histograms_WMregions.pdf")
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


# note: adaptive thresholds are too low in this dataset for lib_size and n_features
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


# -------------------
# plot summary values
# -------------------

# summary of discarded spots in LC and WM regions by sample
table(colData(spe_LC)$sample_part_id, colData(spe_LC)$discard)
table(colData(spe_WM)$sample_id, colData(spe_WM)$discard)


# summary values of QC metrics in LC and WM regions by sample

cols_keep <- c("sample_id", "sum", "detected", "annot_region", "annot_spot", "discard")

df_summary_by_sample <- 
  colData(spe)[, cols_keep] %>% 
  as.data.frame() %>% 
  mutate(region = ifelse(annot_region, "LC", "WM")) %>% 
  mutate(region = factor(region, levels = c("LC", "WM"))) %>% 
  mutate(spots = ifelse(annot_spot, "spot", "nonspot")) %>% 
  mutate(spots = factor(spots, levels = c("spot", "nonspot"))) %>% 
  group_by(sample_id, region) %>% 
  summarize(medSumUMI = median(sum), 
            medSumGenes = median(detected), 
            sumDiscard = sum(discard)) %>% 
  pivot_longer(., cols = c("medSumUMI", "medSumGenes", "sumDiscard"), 
               names_to = "metric", values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("medSumUMI", "medSumGenes", "sumDiscard"))) %>% 
  as.data.frame()


#pal <- unname(palette.colors(3, palette = "Alphabet"))
pal <- c("darkorchid1", "dodgerblue", "darkorange")

set.seed(123)
ggplot(df_summary_by_sample, 
       aes(x = metric, y = value, color = metric)) + 
  facet_wrap(~region) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.15) + 
  scale_color_manual(values = pal) + 
  ggtitle("QC summary") + 
  theme_bw()

fn <- file.path(dir_plots, "QC_summary_boxplots")
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# --------------------
# plot discarded spots
# --------------------

# select all spots
df <- cbind.data.frame(colData(spe), spatialCoords(spe))

# plot discarded spots
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated\nregion") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$discard, , drop = FALSE], 
             aes(color = discard), size = 0.1) + 
  scale_color_manual(values = "red") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level quality control") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_discard")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# ----------------
# plot annotations
# ----------------

# plot annotated regions
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray80", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "red", name = "annotated\nregion") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  scale_y_reverse() + 
  ggtitle("Annotated regions") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "annotations_regions")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# plot annotated regions and spots
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray80", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "red", name = "annotated\nregion") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_spot, , drop = FALSE], 
             aes(color = annot_spot), size = 0.1) + 
  scale_color_manual(values = "black", name = "annotated\nspot") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Annotations") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "annotations_regions_spots")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# ---------------
# plot QC metrics
# ---------------

# plot sum UMI counts (values)
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = sum)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_viridis_c(trans = "sqrt") + 
  scale_y_reverse() + 
  labs(color = "counts") + 
  ggtitle("Sum UMI counts") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_sumUMIs_values")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# QC plot: sum UMIs
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated\nregion") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$low_lib_size, , drop = FALSE], 
             aes(color = low_lib_size), size = 0.1) + 
  scale_color_manual(values = "darkorange") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: sum UMIs") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_sumUMIs")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# QC plot: detected genes
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated\nregion") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$low_n_features, , drop = FALSE], 
             aes(color = low_n_features), size = 0.1) + 
  scale_color_manual(values = "darkorange") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: detected genes") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_detectedGenes")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# QC plot: mitochondrial proportion
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated\nregion") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$high_subsets_mito_percent, , drop = FALSE], 
             aes(color = high_subsets_mito_percent), size = 0.1) + 
  scale_color_manual(values = "darkorange", name = "high_mito") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: mitochondrial proportion") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_mitochondrial")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# cell count
df$high_cell_count <- df$cell_count > 10

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated\nregion") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$high_cell_count, , drop = FALSE], 
             aes(color = high_cell_count), size = 0.1) + 
  scale_color_manual(values = "darkorange") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: cell count") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_cellCount")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# ------------------
# plot TH expression
# ------------------

# plot TH expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = TH)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray85", high = "red", trans = "sqrt", 
                       name = "TH counts") + 
  scale_y_reverse() + 
  ggtitle("TH expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "expression_TH")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# plot TH thresholding
thresh <- 3

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = TH >= thresh)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("gray85", "red"), 
                     name = paste0("TH counts >= ", thresh)) + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  scale_y_reverse() + 
  ggtitle("TH expression threshold") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "threshold_TH")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# ----------------------
# plot SLC6A2 expression
# ----------------------

# plot SLC6A2 expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = SLC6A2)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray85", high = "red", trans = "sqrt", 
                       name = "SLC6A2 counts") + 
  scale_y_reverse() + 
  ggtitle("SLC6A2 expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "expression_SLC6A2")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# plot SLC6A2 thresholding
thresh <- 3

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = SLC6A2 >= thresh)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("gray85", "red"), 
                     name = paste0("SLC6A2 counts >= ", thresh)) + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  scale_y_reverse() + 
  ggtitle("SLC6A2 expression threshold") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "threshold_SLC6A2")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# --------------------
# plot TPH2 expression
# --------------------

# add to SPE object
# to do: move into SPE object in previous script
colData(spe)$TPH2 <- counts(spe)[which(rowData(spe)$gene_name == "TPH2"), ]
df <- cbind.data.frame(colData(spe), spatialCoords(spe))

# plot TPH2 expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = TPH2)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray85", high = "red", trans = "sqrt", 
                       name = "TPH2 counts") + 
  scale_y_reverse() + 
  ggtitle("TPH2 expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "expression_TPH2")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# plot TPH2 thresholding
thresh <- 1

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = TPH2 >= thresh)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("gray85", "red"), 
                     name = paste0("TPH2 counts >= ", thresh)) + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  scale_y_reverse() + 
  ggtitle("TPH2 expression threshold") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "threshold_TPH2")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# ----------------------
# plot SLC6A4 expression
# ----------------------

# plot SLC6A4 expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = SLC6A4)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray85", high = "red", trans = "sqrt", 
                       name = "SLC6A4 counts") + 
  scale_y_reverse() + 
  ggtitle("SLC6A4 expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "expression_SLC6A4")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# plot SLC6A4 thresholding
thresh <- 3

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = SLC6A4 >= thresh)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("gray85", "red"), 
                     name = paste0("SLC6A4 counts >= ", thresh)) + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  scale_y_reverse() + 
  ggtitle("SLC6A4 expression threshold") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "threshold_SLC6A4")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# ------------------------------
# plot additional selected genes
# ------------------------------

genes_nicotinic_acetylcholine <- c(
  "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7", 
  "CHRNA9", "CHRNA10", "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4")

genes_serotonin <- c("HTR1A", "HTR2A")

genes <- c(genes_nicotinic_acetylcholine, genes_serotonin)

length(genes)
sum(genes %in% rowData(spe)$gene_name)


for (g in seq_along(genes)) {
  df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  ix_gene <- which(rowData(spe)$gene_name == genes[g])
  df$gene <- counts(spe)[ix_gene, ]
  
  p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                      color = gene)) + 
    facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
    geom_point(size = 0.1) + 
    scale_color_gradient(low = "gray85", high = "red", trans = "sqrt", 
                         name = paste0(genes[g], " counts")) + 
    scale_y_reverse() + 
    ggtitle(genes[g], " expression") + 
    theme_bw() + 
    theme(aspect.ratio = 1, 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "selected_genes", paste0("expression_", genes[g]))
  ggsave(paste0(fn, ".pdf"), plot = p, width = 7, height = 6.75)
  ggsave(paste0(fn, ".png"), plot = p, width = 7, height = 6.75)
}


# -----------------------------
# heatmap comparing annotations
# -----------------------------

# heatmap comparing spot-level annotation and expression thresholding

# calculate overlaps
thresh <- 2
tbl <- table(threshold_TH = colData(spe)$TH >= thresh, 
             annot_spot = colData(spe)$annot_spot)
tbl

# convert to proportions
tbl_prop <- apply(tbl, 2, function(col) col / sum(col))
tbl_prop

rownames(tbl) <- rownames(tbl_prop) <- 
  c(paste0("TH counts < ", thresh), paste0("TH counts >= ", thresh))
colnames(tbl) <- colnames(tbl_prop) <- 
  c("not annotated spot", "annotated spot")

# convert to matrices
class(tbl) <- "numeric"
class(tbl_prop) <- "numeric"

tbl <- cbind(as.data.frame(tbl), type = "number")
tbl_prop <- cbind(as.data.frame(tbl_prop), type = "proportion")

tbl$TH_expression <- rownames(tbl)
tbl_prop$TH_expression <- rownames(tbl_prop)


df <- rbind(tbl, tbl_prop) %>% 
  pivot_longer(., cols = -c(TH_expression, type), 
               names_to = "annotation", values_to = "value") %>% 
  as.data.frame()


pal <- c("white", "deepskyblue")

ggplot() + 
  geom_tile(data = df[df$type == "proportion", ], 
            aes(x = annotation, y = TH_expression, fill = value)) + 
  geom_text(data = df[df$type == "number", ], 
            aes(x = annotation, y = TH_expression, label = value)) + 
  scale_fill_gradientn(colors = pal, name = "proportion", 
                       limits = c(0, 1), breaks = c(0, 0.5, 1)) + 
  ggtitle("TH vs. annotation") + 
  theme_bw() + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())

fn <- file.path(dir_plots, "heatmap_annotationVsThreshold")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)


# --------------------------------------
# remove discarded spots from SPE object
# --------------------------------------

dim(spe)

spe <- spe[, !colData(spe)$discard]
dim(spe)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LC_qualityControlled")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

