################################################################
# LC analyses: pseudobulked differential expression (DE) testing
# Lukas Weber, Jun 2022
################################################################

# adapting code by Leonardo Collado-Torres, Abby Spangler, and Sowmya Parthiban from:
# https://github.com/LieberInstitute/spatialDLPFC


# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)
library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(ComplexHeatmap)


# directory to save plots
dir_plots <- here("plots", "07_pseudobulkDE")

# directory to save outputs
dir_outputs <- here("outputs", "07_pseudobulkDE")


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


# ----------------
# pseudobulk spots
# ----------------

# pseudobulk spots within manually annotated LC regions and WM regions by sample

table(colData(spe)$sample_id, colData(spe)$annot_region)

annot_region_fctr <- factor(as.numeric(colData(spe)$annot_region), 
                            labels = c("WM", "LC"))
table(annot_region_fctr)

ids <- DataFrame(
  annot_region_pseudo = annot_region_fctr, 
  sample_id_pseudo = colData(spe)$sample_id
)

table(ids$sample_id_pseudo, ids$annot_region_pseudo)

# pseudobulking
spe_pseudo <- aggregateAcrossCells(spe, ids)
dim(spe_pseudo)

# add levels to colData and column names
levs <- paste(
  colData(spe_pseudo)$annot_region_pseudo, 
  colData(spe_pseudo)$sample_id_pseudo, 
  sep = "."
)
colData(spe_pseudo)$level_id <- levs
colnames(spe_pseudo) <- levs

# check
dim(spe_pseudo)
head(colData(spe_pseudo), 3)
counts(spe_pseudo)[1:3, ]


# ---------
# logcounts
# ---------

# calculate logcounts for pseudobulked data
# using default library size scale factors

spe_pseudo <- logNormCounts(spe_pseudo, size.factors = NULL)


# ---------
# filtering
# ---------

# additional filtering of low-expressed genes
# using threshold on total UMI counts summed across all samples
# e.g. total 80 = 8 samples * 10 UMI counts per sample

n_umis <- 80
ix_remove <- rowSums(counts(spe_pseudo)) < n_umis
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
mat <- logcounts(spe_pseudo)

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

# calculate adjusted p-values (FDRs)
fdrs <- p.adjust(p_vals, method = "fdr")


# store corresponding gene names for convenience later
stopifnot(length(fdrs) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs) == rowData(spe_pseudo)$gene_id))

fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

names(fdrs) <- fdrs_gene_names


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

# using spot-level (not pseudobulked) SPE object for plotting
df <- as.data.frame(cbind(colData(spe), spatialCoords(spe)))
ix_most_sig <- which(rowData(spe)$gene_name == most_sig)
df$most_sig <- counts(spe)[ix_most_sig, ]

ggplot(df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", 
                      color = "most_sig")) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                       name = paste0(most_sig, " counts")) + 
  scale_y_reverse() + 
  ggtitle(most_sig) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())


# ----------------------------------------------
# volcano plot: standard significance thresholds
# ----------------------------------------------

# summarize results with volcano plot

# standard significance thresholds

logfc <- res$coefficients[, "annot_region_pseudoLC"]

stopifnot(length(fdrs) == length(logfc))
stopifnot(all(fdrs_gene_ids == names(logfc)))

names(logfc) <- names(fdrs)

# identify significant genes (low FDR and high logFC)
thresh_fdr <- 0.05
thresh_logfc <- log2(2)
sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

# number of significant genes
table(sig)


df <- data.frame(
  gene = names(fdrs), 
  FDR = fdrs, 
  logFC = logfc, 
  sig = sig
)

pal <- c("black", "red")


# volcano plot without labels
ggplot(df, aes(x = logFC, y = -log10(FDR), color = sig)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$sig, ], size = 0.5) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(0.05), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(2), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(2), lty = "dashed", color = "royalblue") + 
  ggtitle("LC vs. WM regions") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "pseudobulkDE_volcano")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)


# ----------------------------------------------------
# volcano plot: more stringent significance thresholds
# ----------------------------------------------------

# summarize results with volcano plot

# more stringent significance thresholds

# identify significant genes (low FDR and high logFC)
thresh_fdr <- 1e-3
thresh_logfc <- log2(3)
stringent <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

# number of significant genes
table(stringent)


df <- data.frame(
  gene = names(fdrs), 
  FDR = fdrs, 
  logFC = logfc, 
  stringent = stringent
)

pal <- c("black", "red")


# volcano plot without labels
ggplot(df, aes(x = logFC, y = -log10(FDR), color = stringent)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$stringent, ], size = 0.5) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(1e-3), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(3), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(3), lty = "dashed", color = "royalblue") + 
  ggtitle("LC vs. WM regions") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "pseudobulkDE_stringent_volcano")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)


# volcano plot with labels
set.seed(123)
ggplot(df, aes(x = logFC, y = -log10(FDR), color = stringent, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$stringent, ], size = 0.5) + 
  geom_text_repel(data = df[df$stringent, ], 
                  size = 1.5, nudge_y = 0.1, 
                  force = 0.1, force_pull = 0.1, min.segment.length = 0.1, 
                  max.overlaps = 20) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(1e-3), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(3), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(3), lty = "dashed", color = "royalblue") + 
  ggtitle("LC vs. WM regions") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "pseudobulkDE_stringent_volcanoWithLabels")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)


# -------
# heatmap
# -------

stopifnot(nrow(mat) == nrow(spe_pseudo))
stopifnot(all(rownames(mat) == rowData(spe_pseudo)$gene_id))
rownames(mat) <- rowData(spe_pseudo)$gene_name


# calculate mean logcounts across samples in LC and WM regions
# where mean is calculated as unweighted mean across samples
mat_LC <- mat[, 9:16]
mat_WM <- mat[, 1:8]

mean_LC <- rowMeans(mat_LC)
mean_WM <- rowMeans(mat_WM)

hmat <- cbind(
  LC = mean_LC, 
  WM = mean_WM
)


# select top genes
# ordered by FDR within set of significant genes (stringent set)
stopifnot(length(fdrs) == length(stringent))
top <- sort(fdrs[stringent])
top_names <- names(top)

# order matrix rows
hmat <- hmat[top_names, ]

# format FDRs in row names
nms <- paste0(names(top), " (", format(signif(top, 2)), ")")
rownames(hmat) <- nms


# create heatmap (horizontal format)
hm <- Heatmap(
  t(hmat), 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  row_names_side = "left", row_names_gp = gpar(fontsize = 10), 
  column_names_gp = gpar(fontsize = 9, fontface = "italic"), 
  name = "mean\nlogcounts"
)

hm

# save heatmap (horizontal format)
fn <- file.path(dir_plots, "pseudobulkDE_heatmap_horizontal")

pdf(paste0(fn, ".pdf"), width = 6.5, height = 3.25)
hm
dev.off()

png(paste0(fn, ".png"), width = 6.5 * 200, height = 3.25 * 200, res = 200)
hm
dev.off()


# create heatmap (vertical format)
hm <- Heatmap(
  hmat, 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  column_names_rot = 0, column_names_gp = gpar(fontsize = 10), column_names_centered = TRUE, 
  row_names_gp = gpar(fontsize = 9, fontface = "italic"), 
  column_title = "LC vs. WM regions", 
  column_title_gp = gpar(fontface = "bold"), 
  name = "mean\nlogcounts"
)

hm

# save heatmap (vertical format)
fn <- file.path(dir_plots, "pseudobulkDE_heatmap_vertical")

pdf(paste0(fn, ".pdf"), width = 3.75, height = 6.5)
hm
dev.off()

png(paste0(fn, ".png"), width = 3.75 * 200, height = 6.5 * 200, res = 200)
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
  logFC = logfc, 
  FDR = fdrs, 
  sig = sig
)


pal <- c("black", "red")

ggplot(df, aes(x = mean, y = logFC, color = sig, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$sig, ], size = 0.5) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = log2(2), lty = "dashed", color = "royalblue") + 
  geom_hline(yintercept = -log2(2), lty = "dashed", color = "royalblue") + 
  labs(x = "mean logcounts (pseudobulked LC and WM)", 
       y = "log fold change") + 
  ggtitle("Pseudobulked DE tests: LC vs. WM") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "pseudobulkDE_MAplot")
ggsave(paste0(fn, ".pdf"), width = 4.25, height = 4)
ggsave(paste0(fn, ".png"), width = 4.25, height = 4)


# -----------------
# export gene lists
# -----------------

stopifnot(length(fdrs_gene_ids) == nrow(rowData(spe_pseudo)))
stopifnot(all(fdrs_gene_ids == rowData(spe_pseudo)$gene_id))

# all genes
df_all <- data.frame(
  gene_id = fdrs_gene_ids, 
  gene_name = fdrs_gene_names, 
  source = rowData(spe_pseudo)$source, 
  type = rowData(spe_pseudo)$type, 
  gene_version = rowData(spe_pseudo)$gene_version, 
  gene_type = rowData(spe_pseudo)$gene_type, 
  mean_logcounts_WM = mean_WM, 
  mean_logcounts_LC = mean_LC, 
  mean_logcounts_LCWM = (mean_WM + mean_LC) / 2, 
  logFC = logfc, 
  pval = p_vals, 
  FDR = fdrs, 
  significant = sig, 
  stringent = stringent
)

# subset significant and stringent genes
df_sig <- df_all[df_all$significant, ]
df_stringent <- df_all[df_all$stringent, ]

dim(df_all)
dim(df_sig)
dim(df_stringent)


# order in most user-friendly way for each set (by gene IDs or FDRs)
df_all <- df_all[order(df_all$gene_id), ]
df_sig <- df_sig[order(df_sig$FDR), ]
df_stringent <- df_stringent[order(df_stringent$FDR), ]


# save .csv files
fn_all <- file.path(dir_outputs, "LC_pseudobulkDE_all.csv")
write.csv(df_all, file = fn_all, row.names = FALSE)

fn_sig <- file.path(dir_outputs, "LC_pseudobulkDE_sigGenes.csv")
write.csv(df_sig, file = fn_sig, row.names = FALSE)

fn_stringent <- file.path(dir_outputs, "LC_pseudobulkDE_stringentGenes.csv")
write.csv(df_stringent, file = fn_stringent, row.names = FALSE)

