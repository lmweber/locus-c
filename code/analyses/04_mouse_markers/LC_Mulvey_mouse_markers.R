###################################################
# LC project
# Script to plot Mulvey et al. (2018) mouse markers
# Lukas Weber, Apr 2022
###################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(biomaRt)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "04_mouse_markers", "Mulvey")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_qualityControlled.rds")
spe <- readRDS(fn_spe)

dim(spe)


# ---------------------
# load marker gene list
# ---------------------

# list of 45 known marker genes for LC in mouse
# from Mulvey et al. (2018), Figure 2A: https://pubmed.ncbi.nlm.nih.gov/29791834/

# load marker gene names saved in text file
file_markers <- here("inputs", "Mulvey_markers", "Mulvey2018_Fig2A_markers.txt")
mouse_markers <- read.table(file_markers)[, 1]

length(mouse_markers)
mouse_markers


# ----------------------
# convert to human genes
# ----------------------

# convert to homologous human gene names

# using code from biomaRt vignette and 
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genes <- getLDS(
  attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_markers, 
  mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = TRUE
)

genes

human_markers <- genes[, 2]
human_markers <- sort(human_markers)

# retrieves homologous human gene names for 41 out of the 45 mouse marker genes
length(human_markers)
human_markers

# alternatively:
# human_markers <- c(
#   "ASB4", "CALCR", "CHODL", "CYB561", "DBH", "DDC", "DLK1", "EYA2", "FAM183A", 
#   "FIBCD1", "GAL", "GCH1", "GLRA2", "GNG4", "GPX3", "HCRTR1", "MAOA", "MRAP2", 
#   "MYOM2", "NXPH4", "OVGP1", "PCBD1", "PHOX2A", "PHOX2B", "SLC18A2", "SLC31A1",
#   "SLC6A2", "STBD1", "SYT17", "TH", "TM4SF1", "TRAF3IP2")
# 
# length(human_markers)
# length(human_markers %in% rowData(spe)$gene_name)


# ---------------
# plot expression
# ---------------

# plot expression for each marker in each sample

sample_ids <- levels(colData(spe)$sample_id)
sample_ids

for (s in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s]]
  df <- as.data.frame(cbind(colData(spe_sub), spatialCoords(spe_sub)))
  
  for (g in seq_along(human_markers)) {
    ix_marker <- which(toupper(rowData(spe_sub)$gene_name) == toupper(human_markers[g]))
    stopifnot(length(ix_marker) == 1)
    
    # skip gene if total UMI counts below threshold
    thresh <- 10
    if (sum(counts(spe_sub)[ix_marker, ]) <= thresh) {
      next
    }
    df$marker <- counts(spe_sub)[ix_marker, ]
    
    # plot UMI counts
    p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = marker)) + 
      geom_point(size = 0.33) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_gradient(low = "gray85", high = "red", 
                           trans = "sqrt", breaks = range(df$marker)) + 
      ggtitle(human_markers[g]) + 
      labs(color = "counts") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    if (!dir.exists(here(dir_plots, sample_ids[s]))) {
      dir.create(here(dir_plots, sample_ids[s]), recursive = TRUE)
    }
    fn <- here(dir_plots, sample_ids[s], paste0(sample_ids[s], "_", human_markers[g]))
    ggsave(paste0(fn, ".pdf"), plot = p, width = 4, height = 3)
    ggsave(paste0(fn, ".png"), plot = p, width = 4, height = 3)
  }
}

