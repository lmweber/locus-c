###################################################################
# LC analyses: spot-level deconvolution using cell2location - plots
# Lukas Weber, Oct 2022
###################################################################


library(SpatialExperiment)
library(here)
library(ggplot2)
library(viridis)


# directory to save plots
dir_plots <- here("plots", "Visium", "10_deconvolution_cell2location")


# ---------------
# load SPE object
# ---------------

# load SPE object containing cell2location results

fn_spe <- here("processed_data", "SPE", "LC_cell2location.rds")
spe <- readRDS(fn_spe)

# check
head(colData(spe), 2)

# samples
sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# select columns containing deconvolved clusters
# (note: excludes ambiguous neuron clusters)
colnames <- colnames(colData(spe))
cols <- colnames[grepl("^means", colnames)]
cols


# --------------
# generate plots
# --------------

# plot each cluster in each Visium sample

for (s in seq_along(sample_ids)) {
  
  # select sample
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s]]
  
  df <- as.data.frame(cbind(colData(spe_sub), spatialCoords(spe_sub)))
  
  for (q in seq_along(cols)) {
    p <- ggplot(df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", 
                               color = cols[q])) + 
      geom_point(size = 0.5) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_viridis(option = "magma", name = "abundance") + 
      labs(title = paste0("Cluster ", gsub("^.*sf_", "", cols[q])), 
           subtitle = sample_ids[s]) + 
      theme_bw() + 
      theme(panel.background = element_rect(fill = "gray80"), 
            panel.grid = element_line(color = "gray80"), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    if (!dir.exists(here(dir_plots, sample_ids[s]))) {
      dir.create(here(dir_plots, sample_ids[s]), recursive = TRUE)
    }
    fn <- here(dir_plots, sample_ids[s], 
               paste0("cell2location_", sample_ids[s], "_cluster", 
                      gsub("^.*sf_", "", cols[q])))
    ggsave(paste0(fn, ".pdf"), plot = p, width = 4, height = 3)
    ggsave(paste0(fn, ".png"), plot = p, width = 4, height = 3)
  }
}


# plot each cluster across all Visium samples

max_abundance <- max(as.matrix(colData(spe)[, cols]))

for (q in seq_along(cols)) {
  
  # select cluster
  df <- as.data.frame(cbind(
    colData(spe)[, c("sample_id", cols[q]), drop = FALSE], 
    spatialCoords(spe)
  ))
  
  p <- ggplot(df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", 
                             color = cols[q])) + 
    facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
    geom_point(size = 0.1) + 
    scale_y_reverse() + 
    scale_color_viridis(option = "magma", name = "abundance", trans = "sqrt", 
                        limits = c(0, max_abundance)) + 
    labs(title = paste0("Cluster ", gsub("^.*sf_", "", cols[q]))) + 
    theme_bw() + 
    theme(panel.background = element_rect(fill = "gray80"), 
          panel.grid = element_line(color = "gray80"), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  if (!dir.exists(here(dir_plots, "all_samples"))) {
    dir.create(here(dir_plots, "all_samples"), recursive = TRUE)
  }
  fn <- here(dir_plots, "all_samples", 
             paste0("cell2location_cluster", gsub("^.*sf_", "", cols[q])))
  ggsave(paste0(fn, ".pdf"), plot = p, width = 7.5, height = 4)
  ggsave(paste0(fn, ".png"), plot = p, width = 7.5, height = 4)
}

