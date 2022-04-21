#####################################################################
# LC analyses: spot-level deconvolution using cell2location: plotting
# Lukas Weber, Apr 2022
#####################################################################


library(SpatialExperiment)
library(here)
library(ggplot2)
library(viridis)


# directory to save plots
dir_plots <- here("plots", "06_deconvolution", "cell2location", "LCregions")


# ---------------
# load SPE object
# ---------------

# load SPE object containing cell2location results from previous script

fn_spe <- here("processed_data", "SPE", "LC_cell2location_LCregions.rds")
spe <- readRDS(fn_spe)

# check
head(colData(spe), 2)

# samples
sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# names of columns containing deconvolved cell types
cols <- c(
  "q05cell_abundance_w_sf_Astro", 
  "q05cell_abundance_w_sf_Endo.Mural", 
  "q05cell_abundance_w_sf_Excit", 
  "q05cell_abundance_w_sf_Inhib", 
  "q05cell_abundance_w_sf_Micro.Macro", 
  "q05cell_abundance_w_sf_Neuron.NE", 
  "q05cell_abundance_w_sf_Oligo", 
  "q05cell_abundance_w_sf_OPC"
)


# --------------
# generate plots
# --------------

# plot each deconvolved cell type in each sample

for (s in seq_along(sample_ids)) {
  
  # select sample
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s]]
  
  df <- as.data.frame(cbind(colData(spe_sub), spatialCoords(spe_sub)))
  
  for (q in seq_along(cols)) {
    p <- ggplot(df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", 
                          color = cols[q])) + 
      geom_point(size = 1) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_viridis(option = "magma", 
                          name = gsub("^.*_", "", cols[q])) + 
      ggtitle(paste0(sample_ids[s], ": ", gsub("^.*_", "", cols[q]))) + 
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
               paste0(sample_ids[s], "_LCregions_", gsub("^.*_", "", cols[q])))
    ggsave(paste0(fn, ".pdf"), plot = p, width = 4, height = 3)
    ggsave(paste0(fn, ".png"), plot = p, width = 4, height = 3)
  }
}

