########################################
# LC snRNA-seq analyses: quality control
# Lukas Weber, Oct 2022
########################################


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)


dir_plots <- here("plots", "singleNucleus", "03_quality_control")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE", "sce_doubletsRemoved")
sce <- readRDS(paste0(fn, ".rds"))

table(colData(sce)$Sample)


# --------------------
# Quality control (QC)
# --------------------

# perform QC on sum UMI counts and number of detected genes
# note: not using mitochondrial proportion (due to biological reasons in LC-NE neurons)

# store QC metrics
sce <- addPerCellQC(sce, subsets = list(Mito = which(seqnames(sce) == "chrM")))

# check distributions
range(colData(sce)$sum)
range(colData(sce)$detected)
quantile(colData(sce)$sum, seq(0, 1, by = 0.1))
quantile(colData(sce)$detected, seq(0, 1, by = 0.1))

# check 3 median absolute deviations (MADs)
reasons <- perCellQCFilters(sce)
attr(reasons$low_lib_size, "thresholds")
attr(reasons$low_n_features, "thresholds")

# note: 3 MADs is outside range of values (minimum and maximum) for both sum sum
# UMIs and detected genes, so we keep all cells

colData(sce)$discard <- FALSE


# note high mitochondrial percentages
range(colData(sce)$subsets_Mito_percent)
quantile(colData(sce)$subsets_Mito_percent, seq(0, 1, by = 0.1))
quantile(colData(sce)$subsets_Mito_percent, seq(0.9, 1, by = 0.01))
mean(colData(sce)$subsets_Mito_percent > 10)
mean(colData(sce)$subsets_Mito_percent > 20)


# plot QC metrics
p <- gridExtra::grid.arrange(
  plotColData(sce, x = "Sample", y = "sum", colour_by = "subsets_Mito_percent") + 
    scale_y_log10() + guides(color = guide_legend(title = "mito")) + ggtitle("Total count"), 
  plotColData(sce, x = "Sample", y = "detected", colour_by = "subsets_Mito_percent") + 
    scale_y_log10() + guides(color = guide_legend(title = "mito")) + ggtitle("Detected genes"), 
  plotColData(sce, x = "Sample", y = "subsets_Mito_percent") + 
    ggtitle("Mito percent"), 
  ncol = 3
)

p

fn <- file.path(dir_plots, "QC_metrics")
ggsave(paste0(fn, ".pdf"), plot = p, width = 12, height = 3.5)
ggsave(paste0(fn, ".png"), plot = p, width = 12, height = 3.5)


# -------------------------
# High-mitochondrial nuclei
# -------------------------

# investigate nuclei with high mitochondrial proportion of reads

nonzero_TH <- counts(sce)[which(rowData(sce)$gene_name == "TH"), ] > 0

# proportion of nuclei with nonzero expression of TH
mean(nonzero_TH)

# mitochondrial proportion in nuclei with nonzero expression of TH
summary(colData(sce)$subsets_Mito_percent[nonzero_TH])
quantile(colData(sce)$subsets_Mito_percent[nonzero_TH], seq(0, 1, by = 0.1))


# -----------
# Save object
# -----------

fn_out <- here("processed_data", "SCE", "sce_qualityControlled")
saveRDS(sce, paste0(fn_out, ".rds"))

