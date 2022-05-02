################################################
# LC project
# Script to plot Grimm et al. (2004) mouse genes
# Lukas Weber, Apr 2022
################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "05_mouse_markers")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_qualityControlled.rds")
spe <- readRDS(fn_spe)

dim(spe)


# --------------------
# load mouse gene list
# --------------------

# list of mouse genes for LC from Grimm et al. (2004), Figure 4

# converted to human genes by searching individually at https://www.ncbi.nlm.nih.gov/gene/

mouse_genes <- c(
  "CBR3", "DNAH5", "SERPINE1", "LAYN", "TPH2", "RPH3AL", "NGB", "CYB561", 
  "GNAS", "SLC31A1", "TCP1", "PPIC", "COLEC10", "RAB3B", "MAOA", "PCBP3", 
  "TSPAN12", "FBP1", "DBH", "SERPINF1", "TXK", "SEC16B", "TRAF1", "PTGES", 
  "GGT5", "MMP2", "MCAM", "TFAP2A", "ACSL1", "UPB1", "UCP3", "COL5A1", 
  "ALDH1A1", "FBN1")

length(mouse_genes)

# 34 out of 34 genes present in SPE object
sum(mouse_genes %in% rowData(spe)$gene_name)


# ---------------
# plot expression
# ---------------

# plot expression for each gene in each sample

sample_ids <- levels(colData(spe)$sample_id)
sample_ids

for (s in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s]]
  df <- as.data.frame(cbind(colData(spe_sub), spatialCoords(spe_sub)))
  
  for (g in seq_along(mouse_genes)) {
    ix_gene <- which(rowData(spe_sub)$gene_name == mouse_genes[g])
    stopifnot(length(ix_gene) == 1)
    
    # skip gene if zero expression
    if (sum(counts(spe_sub)[ix_gene, ]) == 0) {
      next
    }
    
    df$gene <- counts(spe_sub)[ix_gene, ]
    
    # plot UMI counts
    p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = gene)) + 
      geom_point(size = 0.3) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_gradient(low = "gray85", high = "red", 
                           trans = "sqrt", breaks = range(df$gene)) + 
      ggtitle(mouse_genes[g]) + 
      labs(color = "counts") + 
      theme_bw() + 
      theme(title = element_text(face = "italic"), 
            legend.title = element_text(face = "plain"), 
            panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    if (!dir.exists(here(dir_plots, "Grimm", sample_ids[s]))) {
      dir.create(here(dir_plots, "Grimm", sample_ids[s]), recursive = TRUE)
    }
    fn <- here(dir_plots, "Grimm", sample_ids[s], 
               paste0(sample_ids[s], "_", mouse_genes[g]))
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
# WM regions
table(!colData(spe)$annot_region)
# individual spots
table(colData(spe)$annot_spot)
# individual nonspots
table(!colData(spe)$annot_spot)


sample_ids <- levels(colData(spe)$sample_id)
sample_ids

enrichment <- matrix(NA, nrow = length(sample_ids), ncol = length(mouse_genes))
rownames(enrichment) <- sample_ids
colnames(enrichment) <- mouse_genes

enrichment_LC <- enrichment_WM <- enrichment_spots <- enrichment_nonspots <- enrichment


# select manually annotated LC regions
spe_LC <- spe[, colData(spe)$annot_region]
dim(spe_LC)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(mouse_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_LC)[rowData(spe_LC)$gene_name == mouse_genes[j], 
                                      colData(spe_LC)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_LC[i, j] <- mean_ij
  }
}

# select manually annotated WM regions
spe_WM <- spe[, !colData(spe)$annot_region]
dim(spe_WM)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(mouse_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_WM)[rowData(spe_WM)$gene_name == mouse_genes[j], 
                                      colData(spe_WM)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_WM[i, j] <- mean_ij
  }
}

# select manually annotated individual spots
spe_spots <- spe[, colData(spe)$annot_spot]
dim(spe_spots)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(mouse_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_spots)[rowData(spe_spots)$gene_name == mouse_genes[j], 
                                         colData(spe_spots)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_spots[i, j] <- mean_ij
  }
}

# select manually annotated individual nonspots
spe_nonspots <- spe[, !colData(spe)$annot_spot]
dim(spe_nonspots)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(mouse_genes)) {
    # calculate mean logcounts for gene j, sample i
    mean_ij <- mean(logcounts(spe_nonspots)[rowData(spe_nonspots)$gene_name == mouse_genes[j], 
                                            colData(spe_nonspots)$sample_id == sample_ids[i]])
    # store in matrix (transposed)
    enrichment_nonspots[i, j] <- mean_ij
  }
}


# -----------------------
# plot enrichment summary
# -----------------------

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
meds <- colMedians(enrichment_spots, useNames = TRUE)
genes_ordered <- names(sort(meds, decreasing = TRUE))
all(genes_ordered %in% mouse_genes)
length(genes_ordered) == length(mouse_genes)

df1 <- 
  full_join(df_enrichment_LC, df_enrichment_WM) %>% 
  mutate(region = factor(region, levels = c("LC", "WM"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids)) %>% 
  mutate(gene = factor(gene, levels = genes_ordered)) %>% 
  as.data.frame()

df2 <- 
  full_join(df_enrichment_spots, df_enrichment_nonspots) %>% 
  mutate(region = factor(region, levels = c("spots", "nonspots"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids)) %>% 
  mutate(gene = factor(gene, levels = genes_ordered)) %>% 
  as.data.frame()


pal1 <- c("darkorange", "purple4")
pal2 <- c("dodgerblue", "black")


# LC regions vs. WM regions
p <- ggplot(df1, aes(x = gene, y = mean, color = region)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = pal1) + 
  labs(y = "mean logcounts") + 
  ggtitle("Grimm genes: annotated LC regions") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1))

fn <- here(dir_plots, "enrichment", "Grimm_enrichment_regions")
ggsave(paste0(fn, ".pdf"), plot = p, width = 7.5, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 7.5, height = 4)


# annotated spots vs. nonspots
p <- ggplot(df2, aes(x = gene, y = mean, color = region)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_color_manual(values = pal2) + 
  labs(y = "mean logcounts") + 
  ggtitle("Grimm genes: annotated spots") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1))

fn <- here(dir_plots, "enrichment", "Grimm_enrichment_spots")
ggsave(paste0(fn, ".pdf"), plot = p, width = 7.5, height = 4)
ggsave(paste0(fn, ".png"), plot = p, width = 7.5, height = 4)

