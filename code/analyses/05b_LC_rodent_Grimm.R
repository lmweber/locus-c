################################################
# LC analyses: Grimm et al. (2004) rodent genes
# Lukas Weber, Jun 2022
################################################

# module load conda_R/devel
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
dir_plots <- here("plots", "05_rodent_genes")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_logcounts.rds")
spe <- readRDS(fn_spe)

dim(spe)

# remove samples where NE neurons were not captured
samples_remove <- "Br5459_LC_round2"
spe <- spe[, !(colData(spe)$sample_id %in% samples_remove)]
colData(spe)$sample_id <- droplevels(colData(spe)$sample_id)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# ------------------
# load rat gene list
# ------------------

# list of rat genes for LC from Grimm et al. (2004), Figure 4

# converted to human genes by searching individually at https://www.ncbi.nlm.nih.gov/gene/

human_genes <- c(
  "CBR3", "DNAH5", "SERPINE1", "LAYN", "TPH2", "RPH3AL", "NGB", "CYB561", 
  "GNAS", "SLC31A1", "TCP1", "PPIC", "COLEC10", "RAB3B", "MAOA", "PCBP3", 
  "TSPAN12", "FBP1", "DBH", "SERPINF1", "TXK", "SEC16B", "TRAF1", "PTGES", 
  "GGT5", "MMP2", "MCAM", "TFAP2A", "ACSL1", "UPB1", "UCP3", "COL5A1", 
  "ALDH1A1", "FBN1")

length(human_genes)

# 31 out of 34 genes present in filtered SPE object
sum(human_genes %in% rowData(spe)$gene_name)

# keep genes that are present in filtered SPE object
ix_keep <- which(human_genes %in% rowData(spe)$gene_name)
human_genes <- human_genes[ix_keep]

stopifnot(length(human_genes) == sum(human_genes %in% rowData(spe)$gene_name))


# ---------------
# plot expression
# ---------------

# plot expression for each gene in each sample

for (s in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[s]]
  df <- as.data.frame(cbind(colData(spe_sub), spatialCoords(spe_sub)))
  
  for (g in seq_along(human_genes)) {
    ix_gene <- which(rowData(spe_sub)$gene_name == human_genes[g])
    
    # skip gene if zero expression for this sample
    if (sum(counts(spe_sub)[ix_gene, ]) == 0) next
    
    df$gene <- counts(spe_sub)[ix_gene, ]
    
    # plot UMI counts
    p <- ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                        color = gene)) + 
      geom_point(size = 0.3) + 
      coord_fixed() + 
      scale_y_reverse() + 
      scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                           breaks = range(df$gene), name = "counts") + 
      ggtitle(human_genes[g], 
              subtitle = sample_ids[s]) + 
      theme_bw() + 
      theme(plot.title = element_text(face = "italic"), 
            panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    if (!dir.exists(here(dir_plots, "Grimm", sample_ids[s]))) {
      dir.create(here(dir_plots, "Grimm", sample_ids[s]), recursive = TRUE)
    }
    fn <- here(dir_plots, "Grimm", sample_ids[s], 
               paste0("counts_", human_genes[g], "_", sample_ids[s]))
    ggsave(paste0(fn, ".pdf"), width = 3.75, height = 3.25)
    ggsave(paste0(fn, ".png"), width = 3.75, height = 3.25)
  }
}


# ------------------------------------
# calculate enrichment of marker genes
# ------------------------------------

# enrichment defined in terms of mean logcounts per spot in LC vs. WM regions

# LC region annotations
table(colData(spe)$annot_region)


enrichment <- matrix(NA, nrow = length(sample_ids), ncol = length(human_genes))
rownames(enrichment) <- sample_ids
colnames(enrichment) <- human_genes

enrichment_LC <- enrichment_WM <- enrichment


# select annotated LC regions
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

# select annotated WM regions
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


# ---------------
# plot enrichment
# ---------------

df_enrichment_LC <- as.data.frame(enrichment_LC)
df_enrichment_LC$region <- "LC"
df_enrichment_LC$sample <- rownames(df_enrichment_LC)

df_enrichment_WM <- as.data.frame(enrichment_WM)
df_enrichment_WM$region <- "WM"
df_enrichment_WM$sample <- rownames(df_enrichment_WM)


df_enrichment_LC <- pivot_longer(df_enrichment_LC, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")
df_enrichment_WM <- pivot_longer(df_enrichment_WM, cols = -c(region, sample), 
                                 names_to = "gene", values_to = "mean")


# order genes
meds <- colMedians(enrichment_LC, useNames = TRUE)
genes_ordered <- names(sort(meds, decreasing = TRUE))
all(genes_ordered %in% human_genes)
length(genes_ordered) == length(human_genes)


df <- 
  full_join(df_enrichment_LC, df_enrichment_WM) %>% 
  mutate(regions = factor(region, 
                          levels = c("LC", "WM"), 
                          labels = c("LC regions", "WM regions"))) %>% 
  mutate(sample_id = factor(sample, levels = sample_ids)) %>% 
  na.omit() %>% 
  mutate(gene = factor(gene, levels = genes_ordered)) %>% 
  as.data.frame()

df_rev <- df %>% 
  mutate(regions = fct_rev(regions)) %>% 
  mutate(gene = fct_rev(gene))


pal <- c("darkmagenta", "gray30")
pal_rev <- rev(pal)


# plot enrichment: horizontal format
ggplot(df, aes(x = gene, y = mean, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.75, outlier.size = 0.5) + 
  scale_color_manual(values = pal, name = "annotation") + 
  scale_fill_manual(values = pal, name = "annotation") + 
  labs(y = "mean logcounts per spot") + 
  ggtitle("Mulvey et al. (2018) genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, 
                                   face = "italic", hjust = 1))

fn <- here(dir_plots, "enrichment_Grimm_annotatedRegions_horizontal")
ggsave(paste0(fn, ".pdf"), width = 7.5, height = 4)
ggsave(paste0(fn, ".png"), width = 7.5, height = 4)


# plot enrichment: vertical format
ggplot(df_rev, aes(x = mean, y = gene, color = regions, fill = regions)) + 
  geom_boxplot(alpha = 0.75, outlier.size = 0.5) + 
  scale_color_manual(values = pal_rev, name = "annotation", 
                     guide = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values = pal_rev, name = "annotation", 
                    guide = guide_legend(reverse = TRUE)) + 
  labs(x = "mean logcounts per spot") + 
  ggtitle("Mulvey et al. (2018) genes") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 9, face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "enrichment_Grimm_annotatedRegions_vertical")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 7.5)
ggsave(paste0(fn, ".png"), width = 4.5, height = 7.5)

