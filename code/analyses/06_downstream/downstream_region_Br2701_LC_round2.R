#############################################################
# LC project
# Downstream analyses on spots from manually annotated region
# Lukas Weber, Jan 2022
#############################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(ggplot2)
library(scater)
library(scran)
library(harmony)
library(tidyverse)


# ---------
# load data
# ---------

# load SpatialExperiment object

fn_spe <- here("processed_data", "SPE", "LCrounds1to3_SPE_processed.rds")
spe <- readRDS(fn_spe)
dim(spe)

# select sample
spe <- spe[, colData(spe)$sample_id == "Br2701_LC_round2"]
dim(spe)

# keep only spots over tissue
spe <- spe[, spatialData(spe)$in_tissue == 1]
dim(spe)


# identify spots from manually annotated regions (from shiny app)

# separate parts of sample on Visium slide
fn_parts <- here("inputs", "annotations", "annotations_Br2701_LC_round2_parts.csv")
csv_parts <- read.csv(fn_parts)
colnames(csv_parts) <- c("sample_id", "barcode", "part")
colData(spe) <- DataFrame(left_join(as.data.frame(colData(spe)), csv_parts, by = "barcode"))

colData(spe)$part[is.na(colData(spe)$part)] <- "NA"
table(colData(spe)$part)

# manually annotated LC region
fn_LC <- here("inputs", "annotations", "annotations_Br2701_LC_round2_LC.csv")
csv_LC <- read.csv(fn_LC)
colnames(csv_LC) <- c("sample_id", "barcode", "annotated_region")
colData(spe) <- DataFrame(left_join(as.data.frame(colData(spe)), csv_LC, by = "barcode"))

colData(spe)$annotated_region[is.na(colData(spe)$annotated_region)] <- "NA"
table(colData(spe)$annotated_region)

head(colData(spe), 3)


# subset one part and manually annotated LC region
spe <- spe[, colData(spe)$part == "part2" & colData(spe)$annotated_region == "LC"]

dim(spe)  ## 204 spots


# ------------------------------------
# normalization and log-transformation
# ------------------------------------

# quick clustering for pool-based size factors
set.seed(123)
qclus <- quickCluster(spe)
table(qclus)

# calculate size factors
spe <- computeSumFactors(spe, cluster = qclus)
summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 20)

# calculate logcounts
spe <- logNormCounts(spe)
assayNames(spe)


# -----------------
# feature selection
# -----------------

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$symbol)
table(is_mito)
rowData(spe)$symbol[is_mito]

# remove mitochondrial genes
spe <- spe[!is_mito, ]
dim(spe)

# fit mean-variance relationship
dec <- modelGeneVar(spe)
# visualize mean-variance relationship
fit <- metadata(dec)
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)


# ------------------------
# dimensionality reduction
# ------------------------

# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)
reducedDimNames(spe)
dim(reducedDim(spe, "PCA"))

# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
reducedDimNames(spe)
dim(reducedDim(spe, "UMAP"))
# update column names for easier plotting
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)


# ----------
# clustering
# ----------

# graph-based clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)

# alternatively: k-means clustering with fixed/higher number of clusters
set.seed(123)
k <- 10
clus <- kmeans(reducedDim(spe, "PCA"), centers = k)$cluster
table(clus)

# store cluster labels in column 'label' in colData
colLabels(spe) <- factor(clus)


# -------------
# plot clusters
# -------------

library(RColorBrewer)
colors <- brewer.pal(12, "Set3")

#palette.pals()
#colors <- palette("Polychrome 36")
#colors <- palette("Tableau 10")

df <- cbind.data.frame(colData(spe), spatialData(spe), spatialCoords(spe), reducedDim(spe, "UMAP"))

# plot clustering: cluster labels
ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = colors) + 
  ggtitle("UMAP: cluster labels") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "annotated_region", "Br2701_LC_round2_annotated_region_UMAP"), ".pdf"), 
       width = 4.5, height = 4)
ggsave(paste0(here("plots", "annotated_region", "Br2701_LC_round2_annotated_region_UMAP"), ".png"), 
       width = 4.5, height = 4)


# plot clustering: cluster labels
ggplot(df, aes(x = y, y = x, color = label)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = colors) + 
  coord_fixed() + 
  scale_x_reverse() + 
  ggtitle("x-y coords: cluster labels") + 
  labs(x = "x", y = "y") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(paste0(here("plots", "annotated_region", "Br2701_LC_round2_annotated_region_xy"), ".pdf"), 
       width = 4.5, height = 4)
ggsave(paste0(here("plots", "annotated_region", "Br2701_LC_round2_annotated_region_xy"), ".png"), 
       width = 4.5, height = 4)


# --------------------------
# marker gene identification
# --------------------------

# code from OSCA

rownames(spe) <- rowData(spe)$symbol

marker.info <- scoreMarkers(spe, colLabels(spe))
marker.info

# plot for each cluster

for (i in names(marker.info)) {
  chosen <- marker.info[[i]]
  ordered <- chosen[order(chosen$mean.AUC, decreasing = TRUE), ]
  head(ordered[, 1:4])
  plotExpression(spe, features = head(rownames(ordered)), 
                 x = "label", colour_by = "label")
  ggsave(paste0(here("plots", "annotated_region", "Br2701_LC_round2_cluster"), i, ".pdf"), 
         width = 6, height = 6.5, bg = "white")
  ggsave(paste0(here("plots", "annotated_region", "Br2701_LC_round2_cluster"), i, ".png"), 
         width = 6, height = 6.5, bg = "white")
}


# plot genes of interest for all clusters

markers <- c("TH", "DBH", "SLC6A2", "SLC18A2", "DDC", "GCH1", "MAOA")
plotExpression(spe, features = markers, x = "label", colour_by = "label")
ggsave(here("plots", "annotated_region", "Br2701_LC_round2_markers.pdf"), 
       width = 6, height = 6.5, bg = "white")
ggsave(here("plots", "annotated_region", "Br2701_LC_round2_markers.png"), 
       width = 6, height = 6.5, bg = "white")

