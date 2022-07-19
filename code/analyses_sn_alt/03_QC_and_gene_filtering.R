###########################################################
# LC snRNA-seq analyses: quality control and gene filtering
# Lukas Weber, July 2022
###########################################################

# run in interactive session on JHPCE

# screen -S LC
# qrsh -l mem_free=20G,h_vmem=22G,h_fsize=200G -now n
# cd /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn_alt
# module load conda_R/devel
# R


library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)

dir_plots <- here("plots", "snRNAseq_alt")


# ---------------
# Load SCE object
# ---------------

# load SCE object from previous script

fn <- here("processed_data", "SCE_alt", "sce_doubletFiltered")
sce <- readRDS(paste0(fn, ".rds"))

table(colData(sce)$Sample)


# --------------------
# Quality control (QC)
# --------------------

# perform QC on sum UMI counts and number of detected genes
# note: not using mitochondrial percentage, since this is high for biological reasons in LC-NE neurons

# store QC metrics
sce <- addPerCellQC(sce, subsets = list(Mito = which(seqnames(sce) == "chrM")))

# check distributions
quantile(colData(sce)$sum, seq(0, 1, by = 0.1))
quantile(colData(sce)$detected, seq(0, 1, by = 0.1))

# check 3 median absolute deviations (MADs)
reasons <- perCellQCFilters(sce)
attr(reasons$low_lib_size, "thresholds")
attr(reasons$low_n_features, "thresholds")

# note: 3 MADs is below minimum values observed for both sum UMIs and detected genes, 
# and we do not want to filter out high values since LC-NE neurons are very large; 
# so we keep all cells

colData(sce)$discard <- FALSE


# note high mitochondrial percentages
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
  ncol=3
)

p

fn <- file.path(dir_plots, "QC_metrics")
ggsave(paste0(fn, ".pdf"), plot = p, width = 10, height = 3)
ggsave(paste0(fn, ".png"), plot = p, width = 10, height = 3)


# investigate high-mitochondrial nuclei
nonzero_TH <- counts(sce)[which(rowData(sce)$gene_name == "TH"), ] > 0
summary(colData(sce)$subsets_Mito_percent[nonzero_TH])


# --------------------------
# Filter low-expressed genes
# --------------------------

# note: keep both protein-coding and non-protein-coding genes

n_umis <- 30
ix_remove <- rowSums(counts(sce)) < n_umis
table(ix_remove)

dim(sce)

sce <- sce[!ix_remove, ]

dim(sce)


# -----------
# Save object
# -----------

fn_out <- here("processed_data", "SCE_alt", "sce_QCandGeneFiltered")
saveRDS(sce, paste0(fn_out, ".rds"))

