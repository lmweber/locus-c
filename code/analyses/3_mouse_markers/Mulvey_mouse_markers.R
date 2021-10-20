#########################################
# LC project
# Mulvey et al. (2018) mouse marker genes
# Lukas Weber, Oct 2021
#########################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(biomaRt)
library(ggplot2)


# ---------
# load data
# ---------

# load SpatialExperiment object

fn_spe <- here("processed_data", "SPE", "LCrounds1to3_SPE_processed.rds")
spe <- readRDS(fn_spe)


# ---------------------
# load marker gene list
# ---------------------

# list of 45 known LC marker genes for mouse
# from Mulvey et al. (2018), Figure 2A: https://pubmed.ncbi.nlm.nih.gov/29791834/

# load marker gene names saved in text file
file_markers <- here("inputs", "marker_genes", "Mulvey2018_Fig2A_marker_genes.txt")
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


# ---------------
# plot expression
# ---------------

# plot expression for each marker in each sample

sample_ids <- unique(colData(spe)$sample_id)
sample_ids

for (s in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s] & colData(spe)$tissue == 1]
  df <- as.data.frame(cbind(colData(spe_sub), spatialData(spe_sub), spatialCoords(spe_sub)))
  for (g in seq_along(human_markers)) {
    ix_marker <- which(toupper(rowData(spe_sub)$symbol) == toupper(human_markers[g]))
    stopifnot(length(ix_marker) == 1)
    
    # skip gene if total UMI counts below threshold
    thresh <- 5
    if (sum(counts(spe_sub)[ix_marker, ]) <= thresh) {
      next
    }
    df$marker <- counts(spe_sub)[ix_marker, ]
    
    # plot UMI counts
    p <- ggplot(df, aes(x = y, y = x, color = marker)) + 
      geom_point(size = 0.8) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_gradient(low = "gray90", high = "blue") + 
      ggtitle(human_markers[g]) + 
      labs(color = "UMI counts") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    # save plot
    if (!dir.exists(here("plots", sample_ids[s]))) {
      dir.create(here("plots", sample_ids[s]))
    }
    fn <- here("plots", sample_ids[s], paste0(sample_ids[s], "_", human_markers[g], ".pdf"))
    ggsave(fn, plot = p, width = 5.25, height = 4)
  }
}

