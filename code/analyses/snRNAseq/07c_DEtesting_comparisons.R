###############################################
# LC snRNA-seq analyses: DE testing comparisons
# Lukas Weber, Apr 2023
###############################################


library(here)
library(dplyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)


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


x_snRNAseq <- list(
  `snRNA-seq: NE vs. other neurons` = genes_sig_snRNAseq_NEvsOtherNeuronal, 
  `snRNA-seq: NE vs. all other` = genes_sig_snRNAseq_NEvsAllOther
)

x_snRNAseq_vs_Visium <- list(
  `Visium: pseudobulk` = genes_sig_Visium_pseudobulk, 
  `snRNA-seq: NE vs. all other` = genes_sig_snRNAseq_NEvsAllOther
)


# snRNA-seq (NE vs. other neurons) vs. snRNA-seq (NE vs. all other)
ggVennDiagram(x_snRNAseq) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = c("black", "black")) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

fn <- file.path(dir_plots, "vennDiagram_DEgenes_snRNAseq")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 3.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 3.5)


# Visium (pseudobulk) vs. snRNA-seq (NE vs. all other)
ggVennDiagram(x_snRNAseq_vs_Visium) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = c("black", "black")) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

fn <- file.path(dir_plots, "vennDiagram_DEgenes_Visium_vs_snRNAseq")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 3.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 3.5)

