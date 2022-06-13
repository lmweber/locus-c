###################################################################
# LC analyses: spot-level deconvolution using cell2location - plots
# using merged reference clusters
# Lukas Weber, Jun 2022
###################################################################


library(SpatialExperiment)
library(here)
library(ggplot2)
library(viridis)


# directory to save plots
dir_plots <- here("plots", "07_deconvolution", "cell2location", "merged")


# ---------------
# load SPE object
# ---------------

# load SPE object containing cell2location results

fn_spe <- here("processed_data", "SPE", "LC_cell2location_merged.rds")
spe <- readRDS(fn_spe)

# check
head(colData(spe), 2)

# samples
sample_ids <- levels(colData(spe)$sample_id)
sample_ids

# names of columns containing deconvolved cell types
cols <- c(
  "meanscell_abundance_w_sf_Astro", 
  "meanscell_abundance_w_sf_Endo.Mural", 
  paste0("meanscell_abundance_w_sf_Excit_", LETTERS[1:6]), 
  paste0("meanscell_abundance_w_sf_Inhib_", LETTERS[1:6]), 
  "meanscell_abundance_w_sf_Micro", 
  "meanscell_abundance_w_sf_Neuron.5HT", 
  "meanscell_abundance_w_sf_Neuron.NE", 
  "meanscell_abundance_w_sf_Oligo", 
  "meanscell_abundance_w_sf_OPC"
)


# --------------
# generate plots
# --------------

# plot each cell type in each Visium sample

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
      labs(title = gsub("^.*sf_", "", cols[q]), 
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
               paste0(sample_ids[s], "_", gsub("^.*sf_", "", cols[q])))
    ggsave(paste0(fn, ".pdf"), plot = p, width = 4, height = 3)
    ggsave(paste0(fn, ".png"), plot = p, width = 4, height = 3)
  }
}


# plot each cell type across samples

for (q in seq_along(cols)) {
  
  # select cell type
  df <- as.data.frame(cbind(
    colData(spe)[, c("sample_id", cols[q]), drop = FALSE], 
    spatialCoords(spe)
  ))
  
  p <- ggplot(df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", 
                             color = cols[q])) + 
    facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
    geom_point(size = 0.25) + 
    scale_y_reverse() + 
    scale_color_viridis(option = "magma", name = "abundance", trans = "sqrt") + 
    labs(title = gsub("^.*sf_", "", cols[q])) + 
    theme_bw() + 
    theme(panel.background = element_rect(fill = "gray80"), 
          panel.grid = element_line(color = "gray80"), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- here(dir_plots, paste0("cell2location_merged_", gsub("^.*sf_", "", cols[q])))
  ggsave(paste0(fn, ".pdf"), plot = p, width = 7, height = 4.75)
  ggsave(paste0(fn, ".png"), plot = p, width = 7, height = 4.75)
}

