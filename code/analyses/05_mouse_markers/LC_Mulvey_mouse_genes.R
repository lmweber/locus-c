#################################################
# LC project
# Script to plot Mulvey et al. (2018) mouse genes
# Lukas Weber, June 2022
#################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "05_mouse_markers")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrected.rds")
spe <- readRDS(fn_spe)

dim(spe)


sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# --------------------
# load mouse gene list
# --------------------

# list of 45 mouse genes for LC from Mulvey et al. (2018), Figure 2A

fn <- here("inputs", "Mulvey_markers", "Mulvey2018_Fig2A_markers.txt")
mouse_genes <- read.table(fn)[, 1]

length(mouse_genes)
mouse_genes


# converted to human genes using biomaRt (code below) and by searching for genes
# individually at: https://www.ncbi.nlm.nih.gov/gene/

human_genes <- c(
  "AGTR1", "ASB4", "CALCR", "CALR3", "CHODL", "CHRNA6", "CILP", "CYB561", 
  "DBH", "DDC", "DLK1", "ADGRE1", "EYA2", "SHISAL2B", "FAM183A", "FIBCD1", 
  "GAL", "GCH1", "GLRA2", "GNG4", "GPX3", "GTF2A1L", "HCRTR1", "IGSF5", "MAOA", 
  "MRAP2", "MYOM2", "NEUROG2", "SLC9B2", "NXPH4", "OVGP1", "PCBD1", "PHOX2A", 
  "PHOX2B", "PLA2G4D", "PTGER2", "SLC18A2", "SLC31A1", "SLC6A2", "STBD1", 
  "SYT17", "TH", "TM4SF1", "TM4SF5", "TRAF3IP2")

human_genes <- sort(human_genes)
length(human_genes)

# 42 out of 45 genes present in SPE object
sum(human_genes %in% rowData(spe)$gene_name)


# code to convert to human gene names using biomaRt (not used)

# using code from biomaRt vignette and 
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

# library(biomaRt)
# 
# human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# 
# genes <- getLDS(
#   attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_genes, 
#   mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = TRUE
# )
# 
# human_genes <- sort(genes[, 2])


# ---------------
# plot expression
# ---------------

# plot expression for each gene in each sample

for (s in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s]]
  df <- as.data.frame(cbind(colData(spe_sub), spatialCoords(spe_sub)))
  
  for (g in seq_along(human_genes)) {
    ix_gene <- which(rowData(spe_sub)$gene_name == human_genes[g])
    
    # skip gene if zero expression
    if (sum(counts(spe_sub)[ix_gene, ]) == 0) next
    
    df$gene <- counts(spe_sub)[ix_gene, ]
    
    # plot UMI counts
    p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                        color = gene)) + 
      geom_point(size = 0.3) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_gradient(low = "gray85", high = "red", 
                           trans = "sqrt", breaks = range(df$gene)) + 
      ggtitle(human_genes[g], 
              subtitle = sample_ids[s]) + 
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
               paste0(sample_ids[s], "_", human_genes[g]))
    ggsave(paste0(fn, ".pdf"), plot = p, width = 3.5, height = 3)
    ggsave(paste0(fn, ".png"), plot = p, width = 3.5, height = 3)
  }
}


# --------------------
# calculate enrichment
# --------------------

# enrichment defined as difference in mean logcounts per spot between regions,
# e.g. manually annotated LC vs. WM regions

# to do: alternatively: calculate using pseudobulked logcounts per region


# LC regions
table(colData(spe)$annot_region)
# WM regions
table(!colData(spe)$annot_region)
# individual spots
table(colData(spe)$annot_spot)
# individual nonspots
table(!colData(spe)$annot_spot)


enrichment <- matrix(NA, nrow = length(sample_ids), ncol = length(human_genes))
rownames(enrichment) <- sample_ids
colnames(enrichment) <- human_genes

enrichment_LC <- enrichment_WM <- enrichment_spots <- enrichment_nonspots <- enrichment


# select manually annotated LC regions
spe_LC <- spe[, colData(spe)$annot_region]
dim(spe_LC)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(human_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_LC)[rowData(spe_LC)$gene_name == human_genes[j], 
                                      colData(spe_LC)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_LC[i, j] <- mean_ij
  }
}

# select manually annotated WM regions
spe_WM <- spe[, !colData(spe)$annot_region]
dim(spe_WM)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(human_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_WM)[rowData(spe_WM)$gene_name == human_genes[j], 
                                      colData(spe_WM)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_WM[i, j] <- mean_ij
  }
}

# select manually annotated individual spots
spe_spots <- spe[, colData(spe)$annot_spot]
dim(spe_spots)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(human_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_spots)[rowData(spe_spots)$gene_name == human_genes[j], 
                                         colData(spe_spots)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_spots[i, j] <- mean_ij
  }
}

# select manually annotated individual nonspots
spe_nonspots <- spe[, !colData(spe)$annot_spot]
dim(spe_nonspots)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(human_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_nonspots)[rowData(spe_nonspots)$gene_name == human_genes[j], 
                                            colData(spe_nonspots)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_nonspots[i, j] <- mean_ij
  }
}


# ---------------
# plot enrichment
# ---------------

df_enrichment_LC <- as.data.frame(enrichment_LC)
df_enrichment_LC$region <- "LC"
df_enrichment_LC$sample <- rownames(df_enrichment_LC)

df_enrichment_WM <- as.data.frame(enrichment_WM)
df_enrichment_WM$region <- "WM"
df_enrichment_WM$sample <- rownames(df_enrichment_WM)

df_enrichment_spots <- as.data.frame(enrichment_spots)
df_enrichment_spots$region <- "spots"
df_enrichment_spots$sample <- rownames(df_enrichment_spots)

df_enrichment_nonspots <- as.data.frame(enrichment_nonspots)
df_enrichment_nonspots$region <- "nonspots"
df_enrichment_nonspots$sample <- rownames(df_enrichment_nonspots)


df_enrichment_LC <- pivot_longer(df_enrichment_LC, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")
df_enrichment_WM <- pivot_longer(df_enrichment_WM, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")
df_enrichment_spots <- pivot_longer(df_enrichment_spots, cols = -c(region, sample), 
                                    names_to = "gene", values_to = "mean")
df_enrichment_nonspots <- pivot_longer(df_enrichment_nonspots, cols = -c(region, sample), 
                                       names_to = "gene", values_to = "mean")


# order genes
meds <- colMedians(enrichment_LC, useNames = TRUE)
genes_ordered <- names(sort(meds, decreasing = TRUE))
all(genes_ordered %in% human_genes)
length(genes_ordered) == length(human_genes)

df1 <- 
  full_join(df_enrichment_LC, df_enrichment_WM) %>% 
  mutate(regions = factor(region, 
                          levels = c("LC", "WM"), 
                          labels = c("LC regions", "WM regions"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids)) %>% 
  na.omit() %>% 
  mutate(gene = factor(gene, levels = genes_ordered)) %>% 
  as.data.frame()

df1_rev <- 
  df1 %>% 
  mutate(regions = fct_rev(regions)) %>% 
  mutate(gene = fct_rev(gene))

df2 <- 
  full_join(df_enrichment_spots, df_enrichment_nonspots) %>% 
  mutate(regions = factor(region, 
                          levels = c("spots", "nonspots"), 
                          labels = c("annotated spots", "not annotated spots"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids)) %>% 
  na.omit() %>% 
  mutate(gene = factor(gene, levels = genes_ordered)) %>% 
  as.data.frame()

df2_rev <- 
  df2 %>% 
  mutate(regions = fct_rev(regions)) %>% 
  mutate(gene = fct_rev(gene))


pal <- c("darkmagenta", "gray30")
pal_rev <- rev(pal)


# LC regions vs. WM regions

# horizontal format
ggplot(df1, aes(x = gene, y = mean, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) + 
  scale_color_manual(values = pal, name = "annotation") + 
  scale_fill_manual(values = pal, name = "annotation") + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Mulvey et al. (2018) genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, 
                                   face = "italic", hjust = 1))

fn <- here(dir_plots, "enrichment", "Mulvey_enrichment_annotatedRegions_horizontal")
ggsave(paste0(fn, ".pdf"), width = 7.5, height = 4)
ggsave(paste0(fn, ".png"), width = 7.5, height = 4)


# vertical format
ggplot(df1_rev, aes(x = mean, y = gene, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) + 
  scale_color_manual(values = pal_rev, name = "annotation", 
                     guide = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values = pal_rev, name = "annotation", 
                    guide = guide_legend(reverse = TRUE)) + 
  labs(x = "mean logcounts per spot") + 
  ggtitle("Mulvey et al. (2018) genes") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 9, face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "enrichment", "Mulvey_enrichment_annotatedRegions_vertical")
ggsave(paste0(fn, ".pdf"), width = 5, height = 7.5)
ggsave(paste0(fn, ".png"), width = 5, height = 7.5)


# annotated spots vs. not annotated spots

# horizontal format
ggplot(df2, aes(x = gene, y = mean, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) + 
  scale_color_manual(values = pal, name = "annotation") + 
  scale_fill_manual(values = pal, name = "annotation") + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Mulvey et al. (2018) genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, 
                                   face = "italic", hjust = 1))

fn <- here(dir_plots, "enrichment", "Mulvey_enrichment_annotatedSpots_horizontal")
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# vertical format
ggplot(df2_rev, aes(x = mean, y = gene, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) + 
  scale_color_manual(values = pal_rev, name = "annotation", 
                     guide = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values = pal_rev, name = "annotation", 
                    guide = guide_legend(reverse = TRUE)) + 
  labs(x = "mean logcounts per spot") + 
  ggtitle("Mulvey et al. (2018) genes") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 9, face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "enrichment", "Mulvey_enrichment_annotatedSpots_vertical")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 7.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 7.5)

