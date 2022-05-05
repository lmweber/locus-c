############################################################################
# LC project
# Script for downstream analyses: pseudobulk differential expression testing
# Lukas Weber, May 2022
############################################################################

# adapting code by Leonardo Collado-Torres, Abby Spangler, and Sowmya Parthiban from:
# https://github.com/LieberInstitute/spatialDLPFC


# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)
library(limma)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(ComplexHeatmap)


# directory to save plots
dir_plots <- here("plots", "06_downstream", "pseudobulkDE")

# directory to save outputs
dir_outputs <- here("outputs", "06_downstream", "pseudobulkDE")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_batchCorrected.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)


# remove samples where NE neurons were not captured (see TH enrichment plots)
samples_remove <- "Br5459_LC_round2"
spe <- spe[, !(colData(spe)$sample_id %in% samples_remove)]

colData(spe)$sample_id <- droplevels(colData(spe)$sample_id)

table(colData(spe)$sample_id)


sample_ids <- levels(colData(spe)$sample_id)
sample_ids


# ----------------
# pseudobulk spots
# ----------------

# pseudobulk spots within LC regions and WM regions

# note: pseudobulk by sample_id

table(colData(spe)$sample_id, colData(spe)$annot_region)

annot_region_fctr <- factor(as.numeric(colData(spe)$annot_region), 
                            labels = c("WM", "LC"))
table(annot_region_fctr)

ids <- DataFrame(
  annot_region_pseudo = annot_region_fctr, 
  sample_id_pseudo = colData(spe)$sample_id
)

# pseudobulking
spe_pseudo <- aggregateAcrossCells(spe, ids)

# add levels to colData and column names
levs <- paste(
  colData(spe_pseudo)$annot_region_pseudo, 
  colData(spe_pseudo)$sample_id_pseudo, 
  sep = "."
)
colData(spe_pseudo)$level_id <- levs
colnames(spe_pseudo) <- levs

# check
head(colData(spe_pseudo), 3)
counts(spe_pseudo)[1:3, ]


# calculate logcounts for pseudobulked data
# using default library size scale factors
spe_pseudo <- logNormCounts(spe_pseudo, size.factors = NULL)


# filter extremely low-count genes
# using threshold of sum UMI counts across all samples combined
# e.g. threshold 80 UMI counts = 8 samples * 10 UMI counts per sample
thresh <- 80
ix_remove <- rowSums(counts(spe_pseudo)) < thresh
table(ix_remove)

spe_pseudo <- spe_pseudo[!ix_remove, ]

dim(spe_pseudo)


# -----------------------
# pseudobulked DE testing
# -----------------------

# define model: LC vs. WM regions with blocks by sample
model_formula <- ~annot_region_pseudo
model_matrix <- model.matrix(model_formula, data = colData(spe_pseudo))

model_formula
model_matrix

# calculate intra-block correlation using limma
corfit <- duplicateCorrelation(
  logcounts(spe_pseudo), 
  model_matrix, 
  block = colData(spe_pseudo)$sample_id_pseudo
)
corfit$consensus.correlation


# calculate DE tests per gene using limma
# replace row names with gene names instead of gene IDs for convenience later
mat <- logcounts(spe_pseudo)
rownames(mat) <- rowData(spe_pseudo)$gene_name

res <- eBayes(
  lmFit(
    mat, 
    design = model_matrix, 
    block = colData(spe_pseudo)$sample_id_pseudo, 
    correlation = corfit$consensus.correlation
  )
)


# extract p-values
p_vals <- res$p.value[, "annot_region_pseudoLC"]

# calculate FDRs
fdrs <- p.adjust(p_vals, method = "fdr")


# -------
# results
# -------

# number of significant DE genes
table(fdrs <= 0.05)
table(fdrs <= 1e-2)
table(fdrs <= 1e-3)

# top genes
sort(fdrs[fdrs <= 1e-3])


# check: plot most significant DE gene

# most significant DE gene
which.min(fdrs)
most_sig <- names(which.min(fdrs))

# using original spot-level SPE object
df <- cbind.data.frame(colData(spe), spatialCoords(spe))
df[, most_sig] <- counts(spe)[rowData(spe)$gene_name == most_sig, ]

ggplot(df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", 
                      color = most_sig)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt") + 
  scale_y_reverse() + 
  ggtitle(most_sig) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())


# ------------
# volcano plot
# ------------

# summarize results with volcano plot

logFC <- res$coefficients[, "annot_region_pseudoLC"]

stopifnot(length(fdrs) == length(logFC))
stopifnot(all(names(fdrs) == names(logFC)))

# identify significant genes (low FDR and high logFC)
thresh_sig <- 0.005
thresh_high <- 1
sig_high <- (fdrs <= thresh_sig) & (abs(logFC) >= thresh_high)

table(sig_high)


df <- data.frame(
  gene = names(fdrs), 
  FDR = fdrs, 
  logFC = logFC, 
  sig_high = sig_high
)

pal <- c("black", "red")

set.seed(123)
ggplot(df, aes(x = logFC, y = -log10(FDR), 
               color = sig_high, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$sig_high, ], size = 0.5) + 
  geom_text_repel(data = df[df$sig_high, ], 
                  size = 1.5, nudge_y = 0.1, 
                  force = 0.1, force_pull = 0.1, min.segment.length = 0.1) + 
  scale_color_manual(values = pal, name = "selected") + 
  xlim(c(-3, 3.5)) + 
  ylim(c(0, 4.5)) + 
  ggtitle("Pseudobulk DE tests: LC vs. WM") + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "pseudobulkDE_volcano")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


# -------
# heatmap
# -------

# calculate mean logcounts in LC and WM regions
# where mean is calculated as unweighted mean across samples
mat_WM <- mat[, 1:8]
mat_LC <- mat[, 9:16]

mean_WM <- rowMeans(mat_WM)
mean_LC <- rowMeans(mat_LC)

hmat <- cbind(
  WM = mean_WM, 
  LC = mean_LC
)


# select top genes
# ordered by FDR within set of selected genes from volcano plot above
top <- sort(fdrs[sig_high])
top_names <- names(top)

# order matrix rows
hmat <- hmat[top_names, ]

# format p-values in row names
nms <- paste0(names(top), " (", signif(top, 2), ")")
rownames(hmat) <- nms


# create heatmap
hm <- Heatmap(
  hmat, 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  row_names_gp = gpar(fontsize = 7), 
  name = "mean\nlogcounts"
)

hm


# save heatmap
fn <- file.path(dir_plots, "pseudobulkDE_heatmap")

pdf(paste0(fn, ".pdf"), width = 3.5, height = 13)
hm
dev.off()

png(paste0(fn, ".png"), width = 3.5 * 200, height = 13 * 200, res = 200)
hm
dev.off()


# -------
# MA plot
# -------

df <- data.frame(
  gene = names(fdrs), 
  mean_WM = mean_WM, 
  mean_LC = mean_LC, 
  mean = (mean_WM + mean_LC) / 2, 
  logFC = logFC, 
  FDR = fdrs, 
  sig_high = sig_high
)


pal <- c("black", "red")

set.seed(123)
ggplot(df, aes(x = mean, y = logFC, 
               color = sig_high, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$sig_high, ], size = 0.5) + 
  geom_text_repel(data = df[df$sig_high, ], 
                  size = 1.5, nudge_y = 0.1, 
                  force = 0.1, force_pull = 0.1, min.segment.length = 0.1) + 
  scale_color_manual(values = pal, name = "selected") + 
  labs(x = "mean logcounts (pseudobulked LC and WM)", 
       y = "log fold change") + 
  ggtitle("Pseudobulk DE tests: LC vs. WM") + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "pseudobulkDE_MAplot")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


# ----------------
# export gene list
# ----------------

df <- data.frame(
  gene_name = names(fdrs), 
  mean_logcounts_WM = mean_WM, 
  mean_logcounts_LC = mean_LC, 
  mean_logcounts_LCWM = (mean_WM + mean_LC) / 2, 
  logFC = logFC, 
  pvals = p_vals, 
  FDR = fdrs, 
  sig_high = sig_high
)

# selected genes (highly significant and large logFC)
df <- df[df$gene_name %in% top_names, ]
rownames(df) <- df$gene_name

# order genes
df <- df[top_names, ]

# save .csv file
fn <- file.path(dir_outputs, "LC_pseudobulkDE_topGenes.csv")
write.csv(df, file = fn, row.names = FALSE)

