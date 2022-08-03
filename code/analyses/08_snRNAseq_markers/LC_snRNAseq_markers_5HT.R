########################################################
# LC project
# Script to plot top snRNA-seq markers in Visium samples
# Lukas Weber, June 2022
########################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "08_snRNAseq_markers")


# ---------
# load data
# ---------

# load saved SPE object
fn_spe <- here("processed_data", "SPE", "LC_qualityControlled.rds")
spe <- readRDS(fn_spe)
dim(spe)

# remove samples from donors that were not included in snRNA-seq data (Br8153, Br5459)
samples_remove <- c("Br8153_LC_round2", "Br5459_LC_round2", "Br8153_LC_round3")
spe <- spe[, !(colData(spe)$sample_id %in% samples_remove)]
colData(spe)$sample_id <- droplevels(colData(spe)$sample_id)
dim(spe)

# check sample IDs
table(colData(spe)$sample_id)


# --------------
# load gene list
# --------------

# top markers for 5-HT neurons (1 vs all tests) from snRNA-seq data
# https://github.com/lmweber/locus-c/blob/main/code/analyses_sn/top40genesLists_LC-n3_19cellTypes.csv

markers_5HT <- c(
  "SLC6A4", "TPH2", "DDC", "AC104984.4", "FEV", "LINC01242", "SLC18A2", "NPNT", 
  "LINC01997", "ITGBL1", "GATA3", "SLC10A4", "AC122707.1", "GCH1", "GPR149", 
  "NDNF", "TMPRSS9", "GPC3", "LINC02082", "AL356737.2")

length(markers_5HT)


# ---------------
# plot expression
# ---------------

# plot each gene across samples

for (q in seq_along(markers_5HT)) {
  
  ix <- which(rowData(spe)$gene_name == markers_5HT[q])
  
  if (length(ix) == 0) next
  
  # select marker
  df <- as.data.frame(cbind(
    colData(spe)[, "sample_id", drop = FALSE], 
    spatialCoords(spe), 
    marker = counts(spe)[ix, ]
  ))
  
  p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = marker)) + 
    facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
    geom_point(size = 0.3) + 
    scale_y_reverse() + 
    scale_color_gradient(low = "gray85", high = "red", 
                         trans = "sqrt", breaks = range(df$marker), 
                         name = "counts") + 
    ggtitle(markers_5HT[q]) + 
    theme_bw() + 
    theme(title = element_text(face = "italic"), 
          legend.title = element_text(face = "plain"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- here(dir_plots, "5HT", markers_5HT[q])
  ggsave(paste0(fn, ".pdf"), plot = p, width = 7, height = 4.75)
  ggsave(paste0(fn, ".png"), plot = p, width = 7, height = 4.75)
}

