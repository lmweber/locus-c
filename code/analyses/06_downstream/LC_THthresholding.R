####################################################
# LC project
# Script for downstream analyses: thresholding spots
# Lukas Weber, May 2022
####################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "06_downstream", "THthresholding")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrected.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)


# remove samples where NE neurons were not captured (see TH enrichment plots)
samples_remove <- "Br5459_LC_round2"
spe <- spe[, !(colData(spe)$sample_id %in% samples_remove)]

colData(spe)$sample_id <- droplevels(colData(spe)$sample_id)

table(colData(spe)$sample_id)


sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# ---------------
# TH thresholding
# ---------------

# select TH+ spots to compare with spot-level annotation

thresh <- 2

colData(spe)$THpos <- colData(spe)$TH >= thresh

# number of TH+ spots
table(colData(spe)$THpos)
# number of TH+ spots by sample
table(colData(spe)$sample_id, colData(spe)$THpos)
# TH+ spots vs. spot-level annotation
table(THpos = colData(spe)$THpos, 
      annot_spot = colData(spe)$annot_spot)


# --------------
# plot by sample
# --------------

df <- cbind.data.frame(colData(spe), spatialCoords(spe))

df$THposAnnot <- df$THpos & df$annot_spot
df$THposOnly <- df$THpos & !(df$annot_spot)
df$AnnotOnly <- df$annot_spot & !(df$THpos)

df$selected <- "none"
df$selected[df$THposAnnot] <- "THposAnnot"
df$selected[df$THposOnly] <- "THposOnly"
df$selected[df$AnnotOnly] <- "AnnotOnly"

df$selected <- factor(df$selected, levels = c("THposAnnot", "THposOnly", "AnnotOnly", "none"))
table(df$selected)


pal <- c("black", "blue", "red", "gray80")


# plot each sample separately

for (i in seq_along(sample_ids)) {
  
  # select sample
  ix_sample <- df$sample_id == sample_ids[i]
  df_sub <- df[ix_sample, ]
  
  ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                     color = selected)) + 
    geom_point(size = 0.7) + 
    scale_color_manual(values = pal) + 
    scale_y_reverse() + 
    ggtitle(paste0("TH+ vs. spot annotation: ", sample_ids[i])) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, paste0("THthresholding_", sample_ids[i]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3.75)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3.75)
}

