##################################################
# LC project
# Script for additional quality control (QC) plots
# Lukas Weber, Mar 2022
##################################################

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

fn_spe <- here("processed_data", "SPE", "LC_qualityControlled.rds")
spe <- readRDS(fn_spe)


# ---------------
# plot QC metrics
# ---------------

df <- cbind.data.frame(colData(spe), spatialCoords(spe))

# sum UMIs
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85") + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70") + 
  new_scale_color() + 
  geom_point(data = df[df$low_lib_size, , drop = FALSE], 
             aes(color = low_lib_size), size = 0.1) + 
  scale_color_manual(values = "darkorange") + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: sum UMIs") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_sumUMIs_combinedSampleMetrics")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# detected genes
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85") + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70") + 
  new_scale_color() + 
  geom_point(data = df[df$low_n_features, , drop = FALSE], 
             aes(color = low_n_features), size = 0.1) + 
  scale_color_manual(values = "darkorange") + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: detected genes") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_detectedGenes_combinedSampleMetrics")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# mitochondrial proportion
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85") + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70") + 
  new_scale_color() + 
  geom_point(data = df[df$high_subsets_mito_percent, , drop = FALSE], 
             aes(color = high_subsets_mito_percent), size = 0.1) + 
  scale_color_manual(values = "darkorange") + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: mitochondrial proportion") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_mito_combinedSampleMetrics")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# cell count
df$high_cell_count <- df$cell_count > 10

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85") + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70") + 
  new_scale_color() + 
  geom_point(data = df[df$high_cell_count, , drop = FALSE], 
             aes(color = high_cell_count), size = 0.1) + 
  scale_color_manual(values = "darkorange") + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: cell count") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_cellCount_combinedSampleMetrics")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# ---------------
# plot expression
# ---------------

# TH expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = TH)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray85", high = "darkgreen") + 
  scale_y_reverse() + 
  ggtitle("TH expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_exprTH")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# thresholding on TH
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = TH > 1)) + 
  facet_wrap(~ sample_id, nrow = 3, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("gray85", "darkgreen")) + 
  scale_y_reverse() + 
  ggtitle("TH thresholding") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_thresholdTH")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)

