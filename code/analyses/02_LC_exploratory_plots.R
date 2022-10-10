################################
# LC analyses: exploratory plots
# Lukas Weber, Oct 2022
################################

# module load conda_R/4.2
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
library(scater)


# directory to save plots
dir_plots <- here("plots", "Visium", "02_exploratory_plots")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_preprocessing.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)


# ------------------------
# select genes of interest
# ------------------------

genes_main <- c("DBH", "TH", "SLC6A2", "TPH2", "SLC6A4", "SLC5A7")

genes_additional <- c("SST", "VIP", "PVALB", "MOBP", "MBP", "GAD1", "GAD2")

genes_nicotinic_acetylcholine <- c(
  "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7", 
  "CHRNA9", "CHRNA10", "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4")

genes_serotonin <- c("HTR1A", "HTR2A")


genes_all <- c(
  genes_main, 
  genes_additional, 
  genes_nicotinic_acetylcholine, 
  genes_serotonin)


# ---------------------------------------
# plot genes of interest: multiple panels
# ---------------------------------------

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe)))

sapply(file.path(dir_plots, "multiple_panels", 
                 c("main", "additional", "nicotinic_acetylcholine", "serotonin")), 
       dir.create, recursive = TRUE)


for (g in seq_along(genes_all)) {
  
  ix_gene <- which(rowData(spe)$gene_name == genes_all[g])
  df_plot <- cbind(df, gene = counts(spe)[ix_gene, ])
  
  ggplot(df_plot, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                      color = gene)) + 
    facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
    geom_point(size = 0.1) + 
    scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                         name = "counts", breaks = range(df_plot$gene)) + 
    scale_y_reverse() + 
    ggtitle(paste0(genes_all[g], " expression")) + 
    theme_bw() + 
    theme(aspect.ratio = 1, 
          panel.grid = element_blank(), 
          plot.title = element_text(face = "italic"), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  # save in subdirectories
  if (genes_all[g] %in% genes_main) subdir <- "main"
  if (genes_all[g] %in% genes_additional) subdir <- "additional"
  if (genes_all[g] %in% genes_nicotinic_acetylcholine) subdir <- "nicotinic_acetylcholine"
  if (genes_all[g] %in% genes_serotonin) subdir <- "serotonin"
  
  fn <- file.path(dir_plots, "multiple_panels", subdir, 
                  paste0("counts_", genes_all[g]))
  
  ggsave(paste0(fn, ".pdf"), width = 8.75, height = 4)
  ggsave(paste0(fn, ".png"), width = 8.75, height = 4)
}


# -----------------------------------------
# plot genes of interest: individual panels
# -----------------------------------------

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe)))

for (g in seq_along(genes_main)) {
  
  ix_gene <- which(rowData(spe)$gene_name == genes_main[g])
  df_plot <- cbind(df, gene = counts(spe)[ix_gene, ])
  
  dir.create(file.path(dir_plots, "individual_panels", "main", genes_main[g]), 
             recursive = TRUE)
  
  for (s in seq_along(sample_ids)) {
    
    df_sub <- df_plot[df_plot$sample_id == sample_ids[s], ]
    
    ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = gene)) + 
      geom_point(size = 0.3) + 
      scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                           name = "counts", breaks = range(df_sub$gene)) + 
      coord_fixed() + 
      scale_y_reverse() + 
      ggtitle(genes_main[g]) + 
      theme_bw() + 
      theme(plot.title = element_text(face = "bold.italic"), 
            panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    fn <- file.path(dir_plots, "individual_panels", "main", genes_main[g], 
                    paste0("counts_", genes_main[g], "_", sample_ids[s]))
    ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3)
    ggsave(paste0(fn, ".png"), width = 4.5, height = 3)
  }
}


# ----------------------------------------
# calculate spot-level QC metric summaries
# ----------------------------------------

# without any removal of samples


# filter zeros

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


# calculate QC metrics

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]

# calculate QC metrics using scater package
spe <- addPerCellQCMetrics(spe, subsets = list(mito = is_mito))

head(colData(spe), 2)


# calculate QC summaries by sample
df_qc_summary_by_sample <- 
  colData(spe) %>% 
  as.data.frame() %>% 
  select(c("sample_id", "sum", "detected", "subsets_mito_percent")) %>% 
  group_by(sample_id) %>% 
  summarize(median_sum = median(sum), 
            median_detected = median(detected), 
            median_mito = median(subsets_mito_percent))

df_qc_summary_by_sample


# calculate overall QC summaries
df_qc_summary_overall <- 
  colData(spe) %>% 
  as.data.frame() %>% 
  select(c("sum", "detected", "subsets_mito_percent")) %>% 
  summarize(median_sum = median(sum), 
            median_detected = median(detected), 
            median_mito = median(subsets_mito_percent))

df_qc_summary_overall

