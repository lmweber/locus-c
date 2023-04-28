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
dir_outputs <- here("outputs", "singleNucleus", "07_DEtesting")


# ------------
# Load results
# ------------

# load results from previous scripts

fn_NEvsOtherNeuronal <- here(dir_outputs, "DEtesting_NEvsOtherNeuronal.csv")
fn_NEvsAllOther <- here(dir_outputs, "DEtesting_NEvsAllOther.csv")

res_NEvsOtherNeuronal <- read_csv(fn_NEvsOtherNeuronal)
res_NEvsAllOther <- read_csv(fn_NEvsAllOther)

# check
table(res_NEvsOtherNeuronal$significant)
table(res_NEvsAllOther$significant)


# -------------
# plot overlaps
# -------------

# Venn diagram of overlaps in sets of significant DE genes

gene_names_sig_NEvsOtherNeuronal <- 
  res_NEvsOtherNeuronal |> 
  filter(significant) |> 
  select(gene_name) |> 
  unlist() |> 
  unname()

gene_names_sig_NEvsAllOther <- 
  res_NEvsAllOther |> 
  filter(significant) |> 
  select(gene_name) |> 
  unlist() |> 
  unname()

x <- list(
  NEvsOtherNeuronal = gene_names_sig_NEvsOtherNeuronal, 
  NEvsAllOther = gene_names_sig_NEvsAllOther
)

ggVennDiagram(x) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_color_manual(values = c("black", "black")) + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

fn <- file.path(dir_plots, "vennDiagram_DEgenes_NE_Visium_snRNAseq")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 3.5)
ggsave(paste0(fn, ".png"), width = 5.5, height = 3.5)

