#####################################
# LC Visium analyses: quality control
# Lukas Weber, Sep 2022
#####################################

# module load conda_R/4.2
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(here)
library(SpatialExperiment)
library(scater)
library(scran)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)


# directory to save plots
dir_plots <- here("plots", "03_quality_control")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_preprocessing.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)


# -----------------
# spots over tissue
# -----------------

# check SPE object contains only spots over tissue
table(colData(spe)$in_tissue)


# ------------
# filter zeros
# ------------

# remove genes with zero expression
ix_zero_genes <- rowSums(counts(spe)) == 0
table(ix_zero_genes)

spe <- spe[!ix_zero_genes, ]
dim(spe)

# remove spots with zero expression
ix_zero_spots <- colSums(counts(spe)) == 0
table(ix_zero_spots)

spe <- spe[, !ix_zero_spots]
dim(spe)

# check no zeros or NAs remaining
table(colData(spe)$in_tissue, useNA = "always")
table(rowSums(counts(spe)) == 0, useNA = "always")
table(colSums(counts(spe)) == 0, useNA = "always")


# -------------------------------
# calculate spot-level QC metrics
# -------------------------------

# combined object containing both LC and non-LC regions

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]

# calculate QC metrics using scater package
spe <- addPerCellQCMetrics(spe, subsets = list(mito = is_mito))

head(colData(spe), 2)


# note wide range of values due to combined LC and non-LC regions
summary(colData(spe)$sum)
summary(colData(spe)$detected)
summary(colData(spe)$subsets_mito_percent)

# note large differences by sample
df_qc <- 
  colData(spe) %>% 
  as.data.frame() %>% 
  select(c("sample_id", "sum", "detected", "subsets_mito_percent")) %>% 
  group_by(sample_id) %>% 
  summarize(median_sum = median(sum), 
            median_detected = median(detected), 
            median_mito = median(subsets_mito_percent))

df_qc


# plot histograms of QC metrics

fn <- file.path(dir_plots, "QC_histograms_allSpots.pdf")
pdf(fn, width = 6.5, height = 2.5)
par(mfrow = c(1, 3))
hist(log10(colData(spe)$sum), xlab = "log10 sum", main = "UMIs per spot")
hist(log10(colData(spe)$detected), xlab = "log10 detected", main = "Genes per spot")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mito", main = "Percent mito UMIs")
par(mfrow = c(1, 1))
dev.off()

# number of cells per spot (not used)
fn <- file.path(dir_plots, "QC_cellCount_allSpots.pdf")
pdf(fn, width = 4, height = 4)
hist(colData(spe)$cell_count, xlab = "no. cells", main = "No. cells per spot")
dev.off()


# ----------------------
# plot QC metric medians
# ----------------------

# plot QC metric medians by sample

df_qc <- df_qc %>% 
  pivot_longer(cols = c("median_sum", "median_detected", "median_mito"), 
               names_to = "metric", values_to = "median") %>% 
  mutate(metric = gsub("^median_", "", metric)) %>% 
  mutate(metric = factor(metric, levels = c("sum", "detected", "mito")))

pal <- unname(palette.colors(10, "Tableau 10"))

set.seed(123)
ggplot(as.data.frame(df_qc), 
       aes(x = metric, y = median)) + 
  facet_wrap(~ metric, scales = "free") + 
  geom_boxplot(aes(group = metric), outlier.shape = NA, width = 0.5) + 
  geom_jitter(aes(color = sample_id, shape = sample_id), 
              width = 0.2, size = 2, stroke = 1) + 
  scale_color_manual(values = pal, name = "sample ID") + 
  scale_shape_manual(values = 1:12, name = "sample ID") + 
  ggtitle("QC metrics by sample") + 
  theme_bw() + 
  theme(axis.title.x = element_blank())

fn <- file.path(dir_plots, "QC_metrics_bySample")
ggsave(paste0(fn, ".pdf"), width = 7, height = 4)
ggsave(paste0(fn, ".png"), width = 7, height = 4)


# ---------------------------------
# calculate QC thresholds by sample
# ---------------------------------

# using "3 median absolution deviations (MADs)" methodology from OSCA 
# calculated separately for each sample

# note: not filtering on mitochondrial percentage in this dataset

df <- colData(spe)[, c("sample_id", "sum", "detected", "subsets_mito_sum", 
                       "subsets_mito_detected", "subsets_mito_percent", "total")]

sample_ids <- levels(colData(spe)$sample_id)
df_list <- vector("list", length(sample_ids))
names(df_list) <- sample_ids

for (s in seq_along(sample_ids)) {
  df_sub <- df[df$sample_id == sample_ids[s], ]
  reasons <- perCellQCFilters(df_sub, sub.fields = "subsets_mito_percent")
  df_list[[s]] <- reasons
}

df_reasons <- do.call("rbind", df_list)
stopifnot(nrow(df) == nrow(df_reasons))
df <- cbind(df, df_reasons)


# number of low-quality spots by sample
table(df$sample_id, df$low_lib_size)
table(df$sample_id, df$low_n_features)
table(df$sample_id, df$high_subsets_mito_percent)
table(df$sample_id, df$discard)


# note: not filtering on mitochondrial percentage in this dataset
# re-calculate discarded spots without filtering on mitochondrial percentage
df$discard <- df$low_lib_size | df$low_n_features
table(df$discard)
table(df$sample_id, df$discard)


# store in SPE object
colData(spe) <- cbind(
  colData(spe), 
  df[, c("low_lib_size", "low_n_features", "high_subsets_mito_percent", "discard")]
)


# ----------------------
# plot low-quality spots
# ----------------------

# plot low-quality spots in x-y coordinates

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe)))


# plot low number of UMI counts
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$low_lib_size, , drop = FALSE], 
             aes(color = low_lib_size), size = 0.1) + 
  scale_color_manual(values = "red", name = "low UMIs") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: low UMI counts") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_samples_lowUMIs")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot low number of detected genes
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$low_n_features, , drop = FALSE], 
             aes(color = low_n_features), size = 0.1) + 
  scale_color_manual(values = "red", name = "low detected") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: low detected genes") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_samples_lowDetected")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot high mitochondrial percentage
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$high_subsets_mito_percent, , drop = FALSE], 
             aes(color = high_subsets_mito_percent), size = 0.1) + 
  scale_color_manual(values = "red", name = "high mito") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: high mitochondrial percentage") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_samples_highMito")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot discarded spots
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray85", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "gray70", name = "annotated") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$discard, , drop = FALSE], 
             aes(color = discard), size = 0.1) + 
  scale_color_manual(values = "red", name = "discard") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Spot-level QC: discarded spots") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "QC_samples_discard")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# ----------------
# plot annotations
# ----------------

# plot manual annotations before removing low-quality spots

# plot annotated regions
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray80", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "red", name = "annotated") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  scale_y_reverse() + 
  ggtitle("Annotations") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "annotations_regions")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot annotated regions and spots
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(aes(color = in_tissue), size = 0.1) + 
  scale_color_manual(values = "gray80", name = "tissue") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_region, , drop = FALSE], 
             aes(color = annot_region), size = 0.1) + 
  scale_color_manual(values = "red", name = "annotated\nregions") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
  new_scale_color() + 
  geom_point(data = df[df$annot_spot, , drop = FALSE], 
             aes(color = annot_spot), size = 0.1) + 
  scale_color_manual(values = "black", name = "annotated\nspots") + 
  guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
  scale_y_reverse() + 
  ggtitle("Annotations") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "annotations_regionsAndSpots")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# -----------------------------------
# plot annotations: individual panels
# -----------------------------------

# plot manual annotations before removing low-quality spots

# plot annotated regions
for (s in seq_along(sample_ids)) {
  df_sub <- df[df$sample_id == sample_ids[s], ]
  
  ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
    geom_point(aes(color = in_tissue), size = 0.3) + 
    scale_color_manual(values = "gray80", name = "tissue") + 
    guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
    new_scale_color() + 
    geom_point(data = df_sub[df_sub$annot_region, , drop = FALSE], 
               aes(color = annot_region), size = 0.3) + 
    scale_color_manual(values = "red", name = "annotated") + 
    guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
    coord_fixed() + 
    scale_y_reverse() + 
    ggtitle(sample_ids[s]) + 
    theme_bw() + 
    theme(plot.title = element_text(face = "bold"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "annotations_regions", 
                  paste0("annotations_regions_", sample_ids[s]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3)
}


# plot annotated regions and spots
for (s in seq_along(sample_ids)) {
  df_sub <- df[df$sample_id == sample_ids[s], ]
  
  ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) + 
    geom_point(aes(color = in_tissue), size = 0.3) + 
    scale_color_manual(values = "gray80", name = "tissue") + 
    guides(color = guide_legend(override.aes = list(size = 2), order = 1)) + 
    new_scale_color() + 
    geom_point(data = df_sub[df_sub$annot_region, , drop = FALSE], 
               aes(color = annot_region), size = 0.3) + 
    scale_color_manual(values = "red", name = "annotated\nregions") + 
    guides(color = guide_legend(override.aes = list(size = 2), order = 2)) + 
    new_scale_color() + 
    geom_point(data = df_sub[df_sub$annot_spot, , drop = FALSE], 
               aes(color = annot_spot), size = 0.3) + 
    scale_color_manual(values = "black", name = "annotated\nspots") + 
    guides(color = guide_legend(override.aes = list(size = 2), order = 3)) + 
    coord_fixed() + 
    scale_y_reverse() + 
    ggtitle(sample_ids[s]) + 
    theme_bw() + 
    theme(plot.title = element_text(face = "bold"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "annotations_regionsAndSpots", 
                  paste0("annotations_regionsAndSpots_", sample_ids[s]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3)
}


# --------------------------------------
# remove discarded spots from SPE object
# --------------------------------------

dim(spe)

table(colData(spe)$discard)

spe <- spe[, !colData(spe)$discard]

dim(spe)


# ---------------
# save SPE object
# ---------------

fn_out <- here("processed_data", "SPE", "LC_qualityControlled")
saveRDS(spe, paste0(fn_out, ".rds"))

