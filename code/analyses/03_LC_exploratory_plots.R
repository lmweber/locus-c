################################
# LC analyses: exploratory plots
# Lukas Weber, Jun 2022
################################

# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)


# directory to save plots
dir_plots <- here("plots", "03_exploratory_plots")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_logcounts.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)


# ------------------------
# select genes of interest
# ------------------------

ix_TH <- which(rowData(spe)$gene_name == "TH")
ix_SLC6A2 <- which(rowData(spe)$gene_name == "SLC6A2")
ix_TPH2 <- which(rowData(spe)$gene_name == "TPH2")
ix_SLC6A4 <- which(rowData(spe)$gene_name == "SLC6A4")
ix_MOBP <- which(rowData(spe)$gene_name == "MOBP")
ix_MBP <- which(rowData(spe)$gene_name == "MBP")

df <- as.data.frame(cbind(
  colData(spe), 
  spatialCoords(spe), 
  TH = counts(spe)[ix_TH, ], 
  SLC6A2 = counts(spe)[ix_SLC6A2, ], 
  TPH2 = counts(spe)[ix_TPH2, ], 
  SLC6A4 = counts(spe)[ix_SLC6A4, ], 
  MOBP = counts(spe)[ix_MOBP, ], 
  MBP = counts(spe)[ix_MBP, ]
))


# ---------------------------------------
# plot genes of interest: multiple panels
# ---------------------------------------

# plot TH expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = TH)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                       name = "counts") + 
  scale_y_reverse() + 
  ggtitle("TH expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        plot.title = element_text(face = "italic"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "counts_TH")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot SLC6A2 expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = SLC6A2)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                       name = "counts") + 
  scale_y_reverse() + 
  ggtitle("SLC6A2 expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        plot.title = element_text(face = "italic"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "counts_SLC6A2")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot TPH2 expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = TPH2)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                       name = "counts") + 
  scale_y_reverse() + 
  ggtitle("TPH2 expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        plot.title = element_text(face = "italic"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "counts_TPH2")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot SLC6A4 expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = SLC6A4)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                       name = "counts") + 
  scale_y_reverse() + 
  ggtitle("SLC6A4 expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        plot.title = element_text(face = "italic"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "counts_SLC6A4")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot MOBP expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = MOBP)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                       name = "counts") + 
  scale_y_reverse() + 
  ggtitle("MOBP expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        plot.title = element_text(face = "italic"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "counts_MOBP")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# plot MBP expression
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = MBP)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                       name = "counts") + 
  scale_y_reverse() + 
  ggtitle("MBP expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        plot.title = element_text(face = "italic"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "counts_MBP")
ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
ggsave(paste0(fn, ".png"), width = 9, height = 4)


# -----------------------------------------
# plot genes of interest: individual panels
# -----------------------------------------

# plot TH expression
for (s in seq_along(sample_ids)) {
  df_sub <- df[df$sample_id == sample_ids[s], ]
  
  ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                     color = TH)) + 
    geom_point(size = 0.3) + 
    scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                         name = "counts", breaks = range(df_sub$TH)) + 
    coord_fixed() + 
    scale_y_reverse() + 
    ggtitle("TH") + 
    theme_bw() + 
    theme(plot.title = element_text(face = "bold.italic"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "TH", paste0("counts_TH_", sample_ids[s]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3)
}


# plot SLC6A2 expression
for (s in seq_along(sample_ids)) {
  df_sub <- df[df$sample_id == sample_ids[s], ]
  
  ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                     color = SLC6A2)) + 
    geom_point(size = 0.3) + 
    scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                         name = "counts", breaks = range(df_sub$SLC6A2)) + 
    coord_fixed() + 
    scale_y_reverse() + 
    ggtitle("SLC6A2") + 
    theme_bw() + 
    theme(plot.title = element_text(face = "bold.italic"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "SLC6A2", paste0("counts_SLC6A2_", sample_ids[s]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3)
}


# plot TPH2 expression
for (s in seq_along(sample_ids)) {
  df_sub <- df[df$sample_id == sample_ids[s], ]
  
  ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                     color = TPH2)) + 
    geom_point(size = 0.3) + 
    scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                         name = "counts", breaks = range(df_sub$TPH2)) + 
    coord_fixed() + 
    scale_y_reverse() + 
    ggtitle("TPH2") + 
    theme_bw() + 
    theme(plot.title = element_text(face = "bold.italic"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "TPH2", paste0("counts_TPH2_", sample_ids[s]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3)
}


# plot SLC6A4 expression
for (s in seq_along(sample_ids)) {
  df_sub <- df[df$sample_id == sample_ids[s], ]
  
  ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                     color = SLC6A4)) + 
    geom_point(size = 0.3) + 
    scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                         name = "counts", breaks = range(df_sub$SLC6A4)) + 
    coord_fixed() + 
    scale_y_reverse() + 
    ggtitle("SLC6A4") + 
    theme_bw() + 
    theme(plot.title = element_text(face = "bold.italic"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "SLC6A4", paste0("counts_SLC6A4_", sample_ids[s]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3)
}


# ------------------------------
# plot additional selected genes
# ------------------------------

genes_nicotinic_acetylcholine <- c(
  "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7", 
  "CHRNA9", "CHRNA10", "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4")

genes_serotonin <- c("HTR1A", "HTR2A")


# several of these genes have been filtered out due to low expression
table(genes_nicotinic_acetylcholine %in% rowData(spe)$gene_name)
table(genes_serotonin %in% rowData(spe)$gene_name)


for (g in seq_along(genes_nicotinic_acetylcholine)) {
  df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  # skip if gene has been filtered out
  if (!(genes_nicotinic_acetylcholine[g] %in% rowData(spe)$gene_name)) next
  ix_gene <- which(genes_nicotinic_acetylcholine[g] == rowData(spe)$gene_name)
  df$gene <- counts(spe)[ix_gene, ]
  
  p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                      color = gene)) + 
    facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
    geom_point(size = 0.1) + 
    scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                         name = "counts", breaks = range(df$gene)) + 
    scale_y_reverse() + 
    ggtitle(paste0(genes_nicotinic_acetylcholine[g], " expression")) + 
    theme_bw() + 
    theme(aspect.ratio = 1, 
          plot.title = element_text(face = "bold.italic"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "nicotinic_acetylcholine", 
                  paste0("counts_", genes_nicotinic_acetylcholine[g]))
  ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
  ggsave(paste0(fn, ".png"), width = 9, height = 4)
}


for (g in seq_along(genes_serotonin)) {
  df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  # skip if gene has been filtered out
  if (!(genes_serotonin[g] %in% rowData(spe)$gene_name)) next
  ix_gene <- which(genes_serotonin[g] == rowData(spe)$gene_name)
  df$gene <- counts(spe)[ix_gene, ]
  
  p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                      color = gene)) + 
    facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
    geom_point(size = 0.1) + 
    scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                         name = "counts", breaks = range(df$gene)) + 
    scale_y_reverse() + 
    ggtitle(paste0(genes_serotonin[g], " expression")) + 
    theme_bw() + 
    theme(aspect.ratio = 1, 
          plot.title = element_text(face = "bold.italic"), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, "serotonin", 
                  paste0("counts_", genes_serotonin[g]))
  ggsave(paste0(fn, ".pdf"), width = 9, height = 4)
  ggsave(paste0(fn, ".png"), width = 9, height = 4)
}


# ----------------------------------------
# heatmap comparing spot-level annotations
# ----------------------------------------

# heatmap comparing spot-level annotations and TH expression thresholding

ix_TH <- which(rowData(spe)$gene_name == "TH")

df <- as.data.frame(cbind(
  colData(spe), 
  TH = counts(spe)[ix_TH, ]
))

# calculate overlaps
n_umis <- 2
tbl <- table(threshold_TH = df$TH >= n_umis, 
             annot_spot = df$annot_spot)
tbl

# convert to proportions
tbl_prop <- apply(tbl, 2, function(col) col / sum(col))
tbl_prop

rownames(tbl) <- rownames(tbl_prop) <- 
  c(paste0("TH counts < ", n_umis), paste0("TH counts >= ", n_umis))
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
  scale_fill_gradientn(colors = pal, name = "column\nproportion", 
                       limits = c(0, 1), breaks = c(0, 0.5, 1)) + 
  ggtitle("Spot-level annotation vs. TH expression") + 
  theme_bw() + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 12), 
        panel.grid = element_blank())

fn <- file.path(dir_plots, "heatmap_spotAnnotationVsTHexpression")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

