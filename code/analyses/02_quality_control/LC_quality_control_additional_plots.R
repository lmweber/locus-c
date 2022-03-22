#################################
# LC project
# Script for quality control (QC)
# Lukas Weber, Mar 2022
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

fn_spe <- here("processed_data", "SPE", "LC_preprocessed_QC.rds")
spe <- readRDS(fn_spe)


# --------------
# generate plots
# --------------

df <- cbind.data.frame(colData(spe), spatialCoords(spe))

# sum UMIs
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = sum < 100)) + 
  facet_wrap(~ sample_id, nrow = 3) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "darkorange")) + 
  ggtitle("Spot-level QC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_sumUMI")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# detected genes
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = detected < 100)) + 
  facet_wrap(~ sample_id, nrow = 3) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "darkorange")) + 
  ggtitle("Spot-level QC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_detectedGenes")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# cell count
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = cell_count > 10)) + 
  facet_wrap(~ sample_id, nrow = 3) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "darkorange")) + 
  ggtitle("Spot-level QC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_cellCount")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# expression of TH
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = TH)) + 
  facet_wrap(~ sample_id, nrow = 3) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray85", high = "darkgreen") + 
  ggtitle("TH expression") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_exprTH")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)


# thresholding on TH
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = TH > 0)) + 
  facet_wrap(~ sample_id, nrow = 3) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("gray85", "darkgreen")) + 
  ggtitle("TH thresholding") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "01_quality_control", "QC_thresholdTH")
ggsave(paste0(fn, ".pdf"), width = 7, height = 6.75)
ggsave(paste0(fn, ".png"), width = 7, height = 6.75)

