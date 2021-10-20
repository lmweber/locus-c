#################################
# LC project
# Script for preprocessing and QC
# Lukas Weber, Oct 2021
#################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)


# ---------
# load data
# ---------

# load SpatialExperiment object

fn_spe <- here("processed_data", "SPE", "LCrounds1to3_SPE_raw.rds")
spe <- readRDS(fn_spe)


# -------------------------------
# spot-level quality control (QC)
# -------------------------------

# apply spot-level QC across all samples

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$symbol)
table(is_mito)
rowData(spe)$symbol[is_mito]

# calculate QC metrics using scater package
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

# plot histograms of QC metrics
par(mfrow = c(1, 4))
hist(colData(spe)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe)$detected, xlab = "detected", main = "Genes per spot")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mitochondrial", main = "Percent mito UMIs")
hist(colData(spe)$count, xlab = "number of cells", main = "No. cells per spot")
par(mfrow = c(1, 1))

# keep all spots for now
colData(spe)$discard <- FALSE

colData(spe)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LCrounds1to3_SPE_processed.rds")
saveRDS(spe, fn_out)


# -------------
# summary plots
# -------------

sample_ids <- unique(colData(spe)$sample_id)
sample_ids

for (s in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s] & colData(spe)$tissue == 1]
  df <- as.data.frame(cbind(colData(spe_sub), spatialData(spe_sub), spatialCoords(spe_sub)))
  
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
    ggtitle(human_markers[g]) + 
    labs(color = "UMI counts") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  if (!dir.exists(here("plots", "summary", sample_ids[s]))) {
    dir.create(here("plots", "summary", sample_ids[s]))
  }
  
  fn <- here("plots", "summary", sample_ids[s], 
             paste0(sample_ids[s], "_sum.pdf"))
  ggsave(fn, plot = p, width = 5.25, height = 4)
  
  
  # plot detected genes
  p <- ggplot(df, aes(x = y, y = x, color = detected)) + 
    geom_point(size = 0.75) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_gradient(low = "gray90", high = "darkorchid4") + 
    ggtitle(human_markers[g]) + 
    labs(color = "detected genes") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- here("plots", "summary", sample_ids[s], 
             paste0(sample_ids[s], "_detected.pdf"))
  ggsave(fn, plot = p, width = 5.25, height = 4)
  
  
  # plot expression of TH
  p <- ggplot(df, aes(x = y, y = x, color = TH)) + 
    geom_point(size = 0.75) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_gradient(low = "gray90", high = "darkgreen") + 
    ggtitle(human_markers[g]) + 
    labs(color = "TH") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- here("plots", "summary", sample_ids[s], 
             paste0(sample_ids[s], "_TH.pdf"))
  ggsave(fn, plot = p, width = 5.25, height = 4)
  
  
  # plot expression of SNAP25
  p <- ggplot(df, aes(x = y, y = x, color = SNAP25)) + 
    geom_point(size = 0.75) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_gradient(low = "gray90", high = "darkgreen") + 
    ggtitle(human_markers[g]) + 
    labs(color = "SNAP25") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- here("plots", "summary", sample_ids[s], 
             paste0(sample_ids[s], "_SNAP25.pdf"))
  ggsave(fn, plot = p, width = 5.25, height = 4)
  
}

