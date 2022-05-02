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
dir_plots <- here("plots", "06_downstream", "thresholdSNAP25")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrected.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# -------------------
# SNAP25 thresholding
# -------------------

# select spots to compare with spot-level annotation

ix <- which(rowData(spe)$gene_name == "SNAP25")
colData(spe)$SNAP25 <- counts(spe)[ix, ]

thresh <- 2

colData(spe)$SNAP25pos <- colData(spe)$SNAP25 >= thresh

# number of thresholded spots
table(colData(spe)$SNAP25pos)
# number of thresholded spots by sample
table(colData(spe)$sample_id, colData(spe)$SNAP25pos)
# thresholded spots vs. spot-level annotation
table(annot = colData(spe)$annot_spot, thresh = colData(spe)$SNAP25pos)


# --------------
# plot by sample
# --------------

df <- cbind.data.frame(
  colData(spe), spatialCoords(spe), 
  reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"), reducedDim(spe, "HARM")
)

df$SNAP25posAnnot <- df$SNAP25pos & df$annot_spot
df$SNAP25posOnly <- df$SNAP25pos & !(df$annot_spot)
df$AnnotOnly <- df$annot_spot & !(df$SNAP25pos)

df$selected <- "none"
df$selected[df$SNAP25posAnnot] <- "SNAP25posAnnot"
df$selected[df$SNAP25posOnly] <- "SNAP25posOnly"
df$selected[df$AnnotOnly] <- "AnnotOnly"

df$selected <- factor(df$selected, levels = c("SNAP25posAnnot", "SNAP25posOnly", "AnnotOnly", "none"))
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
    ggtitle(paste0("SNAP25+ vs. spot annotation: ", sample_ids[i])) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, paste0("LC_thresholdSNAP25_", sample_ids[i]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3.75)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3.75)
}


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LC_thresholdSNAP25")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

