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
library(dplyr)
library(tidyr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "04_mouse_markers")


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


# alternatively: if biomaRt site is unresponsive, load human marker names saved previously
# human_markers <- c(
#   "AGTR1", "ASB4", "CALCR", "CALR3", "CHODL", "CHRNA6",
#   "CILP", "CYB561", "DBH", "DDC", "DLK1", "EYA2",
#   "FAM183A", "FIBCD1", "GAL", "GCH1", "GLRA2", "GNG4",
#   "GPX3", "GTF2A1L", "HCRTR1", "MAOA", "MRAP2", "MYOM2",
#   "NEUROG2", "NXPH4", "OVGP1", "PCBD1", "PHOX2A", "PHOX2B",
#   "PLA2G4D", "PTGER2", "SLC18A2", "SLC31A1", "SLC6A2", "STBD1",
#   "SYT17", "TH", "TM4SF1", "TM4SF5", "TRAF3IP2")
# 
# length(human_markers)


# keep human markers that are present in SPE object (38 markers)
human_markers <- human_markers[toupper(human_markers) %in% toupper(rowData(spe)$gene_name)]
length(human_markers)


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
    
    # skip gene if zero expression
    if (sum(counts(spe_sub)[ix_marker, ]) == 0) {
      next
    }
    
    df$marker <- counts(spe_sub)[ix_marker, ]
    
    # plot UMI counts
    p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = marker)) + 
      geom_point(size = 0.3) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_gradient(low = "gray85", high = "red", 
                           trans = "sqrt", breaks = range(df$marker)) + 
      ggtitle(human_markers[g]) + 
      labs(color = "counts") + 
      theme_bw() + 
      theme(title = element_text(face = "italic"), 
            legend.title = element_text(face = "plain"), 
            panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    if (!dir.exists(here(dir_plots, "Mulvey", sample_ids[s]))) {
      dir.create(here(dir_plots, "Mulvey", sample_ids[s]), recursive = TRUE)
    }
    fn <- here(dir_plots, "Mulvey", sample_ids[s], 
               paste0(sample_ids[s], "_", human_markers[g]))
    ggsave(paste0(fn, ".pdf"), plot = p, width = 3.5, height = 2.75)
    ggsave(paste0(fn, ".png"), plot = p, width = 3.5, height = 2.75)
  }
}


# --------------------
# calculate enrichment
# --------------------

# calculate enrichment of markers in manually annotated regions

# LC regions
table(colData(spe)$annot_region)
# individual spots
table(colData(spe)$annot_spot)
# WM regions
table(!colData(spe)$annot_region)


sample_ids <- levels(colData(spe)$sample_id)
sample_ids

enrichment <- matrix(NA, nrow = length(sample_ids), ncol = length(human_markers))
rownames(enrichment) <- sample_ids
colnames(enrichment) <- human_markers

enrichment_LC <- enrichment_spots <- enrichment_WM <- enrichment


# select manually annotated LC regions
spe_LC <- spe[, colData(spe)$annot_region]
dim(spe_LC)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(human_markers)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_LC)[rowData(spe_LC)$gene_name == human_markers[j], 
                                      colData(spe_LC)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_LC[i, j] <- mean_ij
  }
}


# select manually annotated individual spots
spe_spots <- spe[, colData(spe)$annot_spot]
dim(spe_spots)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(human_markers)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_spots)[rowData(spe_spots)$gene_name == human_markers[j], 
                                         colData(spe_spots)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_spots[i, j] <- mean_ij
  }
}


# select manually annotated WM regions
spe_WM <- spe[, !colData(spe)$annot_region]
dim(spe_WM)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(human_markers)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_WM)[rowData(spe_WM)$gene_name == human_markers[j], 
                                      colData(spe_WM)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_WM[i, j] <- mean_ij
  }
}


# -----------------------
# plot enrichment summary
# -----------------------

df_enrichment_LC <- as.data.frame(enrichment_LC)
df_enrichment_LC$region <- "LC"
df_enrichment_LC$sample <- rownames(df_enrichment_LC)

df_enrichment_spots <- as.data.frame(enrichment_spots)
df_enrichment_spots$region <- "spots"
df_enrichment_spots$sample <- rownames(df_enrichment_spots)

df_enrichment_WM <- as.data.frame(enrichment_WM)
df_enrichment_WM$region <- "WM"
df_enrichment_WM$sample <- rownames(df_enrichment_WM)


df_enrichment_LC <- pivot_longer(df_enrichment_LC, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")
df_enrichment_spots <- pivot_longer(df_enrichment_spots, cols = -c(region, sample), 
                                    names_to = "gene", values_to = "mean")
df_enrichment_WM <- pivot_longer(df_enrichment_WM, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")

df <- 
  full_join(df_enrichment_LC, df_enrichment_spots) %>% 
  full_join(., df_enrichment_WM) %>% 
  mutate(region = factor(region, levels = c("spots", "LC", "WM"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids)) %>% 
  mutate(gene = as.factor(gene)) %>% 
  as.data.frame()

pal <- c("darkorange", "purple4", "dodgerblue")


# 3 panels for individual spots, LC regions, and WM regions
p <- ggplot(df, aes(x = gene, y = mean, group = gene, color = region)) + 
  facet_wrap(~region) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = pal) + 
  labs(y = "mean logcounts") + 
  ggtitle("Mulvey marker genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1))

fn <- here(dir_plots, "enrichment", "enrichment_Mulvey_annot_3panels")
ggsave(paste0(fn, ".pdf"), plot = p, width = 12, height = 3)
ggsave(paste0(fn, ".png"), plot = p, width = 12, height = 3)


# side-by-side: individual spots, LC regions, and WM regions
p <- ggplot(df, aes(x = gene, y = mean, color = region)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = pal) + 
  labs(y = "mean logcounts") + 
  ggtitle("Mulvey marker genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1))

fn <- here(dir_plots, "enrichment", "enrichment_Mulvey_annot_1panel")
ggsave(paste0(fn, ".pdf"), plot = p, width = 9, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 9, height = 4)

