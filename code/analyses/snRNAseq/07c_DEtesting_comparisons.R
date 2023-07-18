###############################################
# LC snRNA-seq analyses: DE testing comparisons
# Lukas Weber, Jul 2023
###############################################


library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)
library(ggrepel)


dir_plots <- here("plots", "singleNucleus", "07_DEtesting")

dir_outputs_snRNAseq <- here("outputs", "singleNucleus", "07_DEtesting")
dir_outputs_Visium <- here("outputs", "Visium", "08_pseudobulkDE")


# ------------
# Load results
# ------------

# load results from previous scripts

# snRNA-seq
fn_snRNAseq_NEvsOtherNeuronal <- here(dir_outputs_snRNAseq, "DEtesting_NEvsOtherNeuronal.csv")
fn_snRNAseq_NEvsAllOther <- here(dir_outputs_snRNAseq, "DEtesting_NEvsAllOther.csv")

# Visium
fn_Visium_pseudobulk <- here(dir_outputs_Visium, "LC_pseudobulkDE_all.csv")


# snRNA-seq
res_snRNAseq_NEvsOtherNeuronal <- read_csv(fn_snRNAseq_NEvsOtherNeuronal)
res_snRNAseq_NEvsAllOther <- read_csv(fn_snRNAseq_NEvsAllOther)

# Visium
res_Visium_pseudobulk <- read_csv(fn_Visium_pseudobulk)


# check
table(res_snRNAseq_NEvsOtherNeuronal$significant)
table(res_snRNAseq_NEvsAllOther$significant)
table(res_Visium_pseudobulk$significant)


# -------------
# plot overlaps
# -------------

# Venn diagrams of overlaps in sets of significant DE genes

genes_sig_snRNAseq_NEvsOtherNeuronal <- 
  res_snRNAseq_NEvsOtherNeuronal |> 
  filter(significant) |> 
  select(gene_name) |> 
  unlist() |> 
  unname()

genes_sig_snRNAseq_NEvsAllOther <- 
  res_snRNAseq_NEvsAllOther |> 
  filter(significant) |> 
  select(gene_name) |> 
  unlist() |> 
  unname()

genes_sig_Visium_pseudobulk <- 
  res_Visium_pseudobulk |> 
  filter(significant) |> 
  select(gene_name) |> 
  unlist() |> 
  unname()


# overlap set: Visium (pseudobulk) vs. snRNA-seq (NE vs. all other)
genes_overlap <- genes_sig_Visium_pseudobulk[genes_sig_Visium_pseudobulk %in% genes_sig_snRNAseq_NEvsAllOther]
genes_overlap
length(genes_overlap)


x_snRNAseq <- list(
  `snRNA-seq: NE vs. other neurons` = genes_sig_snRNAseq_NEvsOtherNeuronal, 
  `snRNA-seq: NE vs. all other` = genes_sig_snRNAseq_NEvsAllOther
)

x_snRNAseq_vs_Visium <- list(
  `Visium: pseudobulk` = genes_sig_Visium_pseudobulk, 
  `snRNA-seq: NE vs. all other` = genes_sig_snRNAseq_NEvsAllOther
)


# snRNA-seq (NE vs. other neurons) vs. snRNA-seq (NE vs. all other)
ggVennDiagram(x_snRNAseq, set_size = 3) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = c("black", "black")) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

fn <- file.path(dir_plots, "vennDiagram_DEgenes_snRNAseq")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 3.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 3.5)


# Visium (pseudobulk) vs. snRNA-seq (NE vs. all other)
ggVennDiagram(x_snRNAseq_vs_Visium, set_size = 3) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = c("black", "black")) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

fn <- file.path(dir_plots, "vennDiagram_DEgenes_Visium_vs_snRNAseq")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 3.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 3.5)


# -------------------------
# plot correlations (ranks)
# -------------------------

# rank correlations for sets of significant DE genes

# calculate ranks per gene
res_sig_snRNAseq_NEvsOtherNeuronal <- res_snRNAseq_NEvsOtherNeuronal[res_snRNAseq_NEvsOtherNeuronal$significant, ]
res_sig_snRNAseq_NEvsAllOther <- res_snRNAseq_NEvsAllOther[res_snRNAseq_NEvsAllOther$significant, ]
res_sig_Visium_pseudobulk <- res_Visium_pseudobulk[res_Visium_pseudobulk$significant, ]

res_sig_snRNAseq_NEvsOtherNeuronal$rank <- rank(res_sig_snRNAseq_NEvsOtherNeuronal$p_value, ties.method = "first")
res_sig_snRNAseq_NEvsAllOther$rank <- rank(res_sig_snRNAseq_NEvsAllOther$p_value, ties.method = "first")
res_sig_Visium_pseudobulk$rank <- rank(res_sig_Visium_pseudobulk$p_value, ties.method = "first")

cols_keep <- c("gene_id", "gene_name", "p_value", "rank")

df_snRNAseq_NEvsOtherNeuronal <- 
  cbind(res_sig_snRNAseq_NEvsOtherNeuronal[, cols_keep], results = "snRNAseq_NEvsOtherNeuronal")
df_snRNAseq_NEvsAllOther <- 
  cbind(res_sig_snRNAseq_NEvsAllOther[, cols_keep], results = "snRNAseq_NEvsAllOther")
df_Visium_pseudobulk <- 
  cbind(res_sig_Visium_pseudobulk[, cols_keep], results = "Visium_pseudobulk")

df_plot <- 
  rbind(df_snRNAseq_NEvsOtherNeuronal, df_snRNAseq_NEvsAllOther, df_Visium_pseudobulk) |> 
  select(-c("gene_name", "p_value")) |> 
  pivot_wider(names_from = "results", values_from = "rank")


# snRNA-seq (NE vs. other neurons) vs. snRNA-seq (NE vs. all other)
ggplot(df_plot, aes(x = snRNAseq_NEvsAllOther, y = snRNAseq_NEvsOtherNeuronal)) + 
  geom_point(pch = 1, color = "navy", size = 1.25, stroke = 0.7) + 
  coord_fixed() + 
  xlim(c(0, 433)) + 
  ylim(c(0, 327)) + 
  labs(x = "DE gene rank (snRNA-seq: NE vs. all other)", 
       y = "DE gene rank (snRNA-seq: NE vs. other neurons)") + 
  ggtitle("Ranks of DE genes (snRNA-seq)") + 
  theme_bw()

fn <- file.path(dir_plots, "correlations_DEgenes_snRNAseq")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 4.5)


# Visium (pseudobulk) vs. snRNA-seq (NE vs. all other)
ggplot(df_plot, aes(x = snRNAseq_NEvsAllOther, y = Visium_pseudobulk)) + 
  geom_point(pch = 1, color = "navy", size = 1.25, stroke = 0.7) + 
  coord_fixed() + 
  xlim(c(0, 433)) + 
  ylim(c(0, 437)) + 
  labs(x = "DE gene rank (snRNA-seq: NE vs. all other)", 
       y = "DE gene rank (Visium: pseudobulk)") + 
  ggtitle("Ranks of DE genes (Visium vs. snRNA-seq)") + 
  theme_bw()

fn <- file.path(dir_plots, "correlations_DEgenes_Visium_vs_snRNAseq")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 5.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 5.5)


# ----------------------
# plot correlations (FC)
# ----------------------

# fold change (FC) correlations for sets of significant DE genes

# extract FCs per gene
res_sig_snRNAseq_NEvsAllOther <- res_snRNAseq_NEvsAllOther[res_snRNAseq_NEvsAllOther$significant, ]
res_sig_Visium_pseudobulk <- res_Visium_pseudobulk[res_Visium_pseudobulk$significant, ]

cols_keep_snRNAseq <- c("gene_id", "gene_name", "p_value", "summary_logFC")
cols_keep_Visium <- c("gene_id", "gene_name", "p_value", "log2FC")

df_snRNAseq_NEvsAllOther <- 
  cbind(res_sig_snRNAseq_NEvsAllOther[, cols_keep_snRNAseq], results = "snRNAseq_NEvsAllOther") |> 
  rename(logFC = summary_logFC)
df_Visium_pseudobulk <- 
  cbind(res_sig_Visium_pseudobulk[, cols_keep_Visium], results = "Visium_pseudobulk") |> 
  rename(logFC = log2FC)

df_overlap <- 
  rbind(df_snRNAseq_NEvsAllOther, df_Visium_pseudobulk) |> 
  select(-c("p_value")) |> 
  pivot_wider(names_from = "results", values_from = "logFC") |> 
  na.omit()


# Visium (pseudobulk) vs. snRNA-seq (NE vs. all other)
ggplot(df_overlap, aes(x = snRNAseq_NEvsAllOther, y = Visium_pseudobulk, 
                    label = gene_name)) + 
  geom_point(pch = 1, color = "navy", size = 1.25, stroke = 0.7) + 
  geom_text_repel(size = 2.5, max.overlaps = 20, color = "navy", 
                  fontface = "italic") + 
  coord_fixed() + 
  labs(x = "DE gene log2FC (snRNA-seq: NE vs. all other)", 
       y = "DE gene log2FC (Visium: pseudobulk)") + 
  ggtitle("logFC of DE genes") + 
  theme_bw()

fn <- file.path(dir_plots, "correlations_logFC_DEgenes_Visium_vs_snRNAseq")
ggsave(paste0(fn, ".pdf"), width = 6.5, height = 4.75)
ggsave(paste0(fn, ".png"), width = 6.5, height = 4.75)


# ------------------------
# spreadsheet: overlap set
# ------------------------

# gene IDs for overlap set (ordered by FC in Visium data)
gene_ids_overlap <- df_overlap$gene_id[order(df_overlap$Visium_pseudobulk, decreasing = TRUE)]

tbl_Visium_pseudobulk <- as.data.frame(res_sig_Visium_pseudobulk)
tbl_snRNAseq_NEvsAllOther <- as.data.frame(res_sig_snRNAseq_NEvsAllOther)

rownames(tbl_Visium_pseudobulk) <- tbl_Visium_pseudobulk$gene_id
rownames(tbl_snRNAseq_NEvsAllOther) <- tbl_snRNAseq_NEvsAllOther$gene_id

# select genes in overlap set
tbl_Visium_pseudobulk <- tbl_Visium_pseudobulk[gene_ids_overlap, ]
tbl_snRNAseq_NEvsAllOther <- tbl_snRNAseq_NEvsAllOther[gene_ids_overlap, ]

colnames(tbl_Visium_pseudobulk)[-c(1:6)] <- 
  paste0(colnames(tbl_Visium_pseudobulk)[-c(1:6)], ".Visium_pseudobulk")
colnames(tbl_snRNAseq_NEvsAllOther)[-c(1:2)] <- 
  paste0(colnames(tbl_snRNAseq_NEvsAllOther)[-c(1:2)], ".snRNAseq_NEvsAllOther")

# combine data frames
df_out <- full_join(tbl_Visium_pseudobulk, tbl_snRNAseq_NEvsAllOther, 
                    by = c("gene_id", "gene_name"))

# save .csv file
fn <- file.path(dir_outputs_snRNAseq, 
                "DEtesting_overlap_VisiumPseudobulk_snRNAseqNEvsAllOther.csv")
write.csv(df_out, file = fn, row.names = FALSE)

