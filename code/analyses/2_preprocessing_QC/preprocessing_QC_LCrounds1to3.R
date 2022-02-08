#################################
# LC project
# Script for preprocessing and QC
# Lukas Weber, Feb 2022
#################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(ggplot2)


# ---------
# load data
# ---------

# load SpatialExperiment object

fn_spe <- here("processed_data", "SPE", "LCrounds1to3_SPE_raw.rds")
spe <- readRDS(fn_spe)

dim(spe)

# convert sample IDs to factor
sample_ids <- c(
  "Br6522_LC_1_round1", "Br6522_LC_2_round1", 
  "Br8153_LC_round2", "Br5459_LC_round2", "Br2701_LC_round2", 
  "Br6522_LC_round3", "Br8079_LC_round3", "Br2701_LC_round3", "Br8153_LC_round3"
)
colData(spe)$sample_id <- factor(colData(spe)$sample_id, levels = sample_ids)

# add sample IDs with parts
colData(spe)$sample_part_ids <- paste(colData(spe)$sample_id, colData(spe)$part_id, sep = "_")

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == TRUE]
dim(spe)


# -------------------------------
# spot-level quality control (QC)
# -------------------------------

# spot-level QC across all samples

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
hist(colData(spe)$cell_count, xlab = "number of cells", main = "No. cells per spot")
par(mfrow = c(1, 1))

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

