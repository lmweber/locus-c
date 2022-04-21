##########################################################################
# LC project
# Script for downstream analyses: unsupervised clustering using BayesSpace
# Lukas Weber, Apr 2022
##########################################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(BayesSpace)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "05_downstream", "BayesSpace")


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


# --------------
# run BayesSpace
# --------------

# spatially-aware unsupervised clustering using BayesSpace
# aiming to identify LC vs. WM regions in an unsupervised manner

# add spatial coordinates in format expected by BayesSpace
colData(spe)$row <- colData(spe)$array_row
colData(spe)$col <- colData(spe)$array_col


# run separately for each sample

colData(spe)$cluster.init <- NA
colData(spe)$spatial.cluster <- NA

for (i in seq_along(sample_ids)) {
  
  # select sample
  ix_sample <- colData(spe)$sample_id == sample_ids[i]
  spe_sub <- spe[, ix_sample]
  
  # run BayesSpace
  
  # arguments:
  # q = number of clusters
  # using Harmony batch integrated dimensions
  # nrep = 50,000 for final results
  
  set.seed(123)
  spe_sub <- spatialCluster(
    spe_sub, 
    q = 4, 
    use.dimred = "HARM", 
    d = 15, 
    platform = "Visium", 
    init.method = "mclust", 
    nrep = 5000, 
    burn.in = 1000
  )
  
  # store results in main object
  colData(spe)[ix_sample, "cluster.init"] <- colData(spe_sub)[, "cluster.init"]
  colData(spe)[ix_sample, "spatial.cluster"] <- colData(spe_sub)[, "spatial.cluster"]
}

# check
table(colData(spe)$spatial.cluster)


# -------------
# plot clusters
# -------------

df <- cbind.data.frame(
  colData(spe), spatialCoords(spe), 
  reducedDim(spe, "PCA"), reducedDim(spe, "UMAP"), reducedDim(spe, "HARM")
)

pal <- c(unname(palette.colors(n = 8, palette = "Okabe-Ito")))

#pal <- c("gray80", "black", "navy", "red")


# plot each sample separately

for (i in seq_along(sample_ids)) {
  
  # select sample
  ix_sample <- df$sample_id == sample_ids[i]
  df_sub <- df[ix_sample, ]
  
  df_sub$spatial.cluster <- as.factor(df_sub$spatial.cluster)
  
  ggplot(df_sub, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                     color = spatial.cluster)) + 
    geom_point(size = 0.7) + 
    scale_color_manual(values = pal) + 
    scale_y_reverse() + 
    ggtitle(paste0("BayesSpace: ", sample_ids[i])) + 
    labs(color = "cluster") + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, paste0("LC_clustering_BayesSpace_", sample_ids[i]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3.75)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 3.75)
}

