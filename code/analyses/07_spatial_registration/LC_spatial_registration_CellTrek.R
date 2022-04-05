###################################################
# LC project
# Script for spatial registration of snRNA-seq data
# Lukas Weber, Apr 2022
###################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(CellTrek)
library(Seurat)
library(ggplot2)
library(ggnewscale)


# directory to save plots
dir_plots <- here("plots", "07_spatial_registration")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe_LC <- here("processed_data", "SPE", "LC_batchCorrected_LCregions.rds")
spe_LC <- readRDS(fn_spe_LC)

fn_spe_WM <- here("processed_data", "SPE", "LC_batchCorrected_WMregions.rds")
spe_WM <- readRDS(fn_spe_WM)

dim(spe_LC)
dim(spe_WM)

# note: merging into combined SPE is difficult due to checks on sample_id; save 
# combined SPE in previous scripts instead


# use LC regions object

spe <- spe_LC

# convert sample IDs to factor
sample_ids <- c(
  "Br6522_LC_1_round1", "Br6522_LC_2_round1", 
  "Br8153_LC_round2", "Br5459_LC_round2", "Br2701_LC_round2", 
  "Br6522_LC_round3", "Br8079_LC_round3", "Br2701_LC_round3", "Br8153_LC_round3"
)
colData(spe)$sample_id <- factor(colData(spe)$sample_id, levels = sample_ids)


# ---------------
# load SCE object
# ---------------

fn_sce <- here("processed_data", "SCE", "sce_updated_LC.rda")
load(fn_sce)

# to do: merge clusters


# -------------------------
# convert to Seurat objects
# -------------------------

# CellTrek is integrated into Seurat framework and requires Seurat objects

# convert SCE object

# create Seurat object
seur_nuc <- CreateSeuratObject(
  counts = counts(sce.lc), 
  project = "LC_nuc", 
  assay = "snRNAseq"
)
# add logcounts assay manually
seur_nuc@assays$snRNAseq@data <- logcounts(sce.lc)
# add colData manually
all(rownames(seur_nuc@meta.data) == rownames(colData(sce.lc)))
seur_nuc@meta.data <- cbind(seur_nuc@meta.data, colData(sce.lc))
# check
str(seur_nuc)
seur_nuc@assays$snRNAseq@counts[1:6, 1:6]
seur_nuc@assays$snRNAseq@data[1:6, 1:6]
head(seur_nuc@meta.data)


# convert SPE object

# to do: build Seurat spatial objects from raw data instead

# select one sample
spe <- spe_LC[, colData(spe_LC)$sample_id == "Br6522_LC_1_round1"]
dim(spe)
# simplify barcode names (required by CellTrek)
colnames(spe) <- gsub("^.*_", "", colnames(spe))

# first convert SPE to SCE
spe_as_sce <- SingleCellExperiment(
  assays = list(
    counts = counts(spe), 
    logcounts = logcounts(spe)), 
  rowData = rowData(spe), 
  colData = cbind(colData(spe), spatialCoords(spe))
)
# create Seurat object
seur_st <- CreateSeuratObject(
  counts = counts(spe_as_sce), 
  project = "LC_st", 
  assay = "spatial"
)
# add logcounts assay manually
seur_st@assays$spatial@data <- logcounts(spe_as_sce)
# add colData manually
all(rownames(seur_st@meta.data) == rownames(colData(spe)))
seur_st@meta.data <- cbind(seur_st@meta.data, colData(spe))
# add images manually
scale_factors <- list(
  spot = NULL, 
  fiducial = NULL, 
  hires = imgData(spe)$scaleFactor[1], 
  lowres = imgData(spe)$scaleFactor[2])
class(scale_factors) <- "scalefactors"
img <- new(
  Class = "VisiumV1", 
  image = array(), 
  scale.factors = scale_factors, 
  coordinates = data.frame(
    tissue = colData(spe)$in_tissue, 
    row = spatialCoords(spe)[, "pxl_row_in_fullres"], 
    col = spatialCoords(spe)[, "pxl_col_in_fullres"], 
    imagerow = colData(spe)$array_row, 
    imagecol = colData(spe)$array_col), 
  spot.radius = numeric(), 
  assay = "spatial", 
  key = "mysample_"
)
seur_st@images <- list(mysample_ = img)
# check
str(seur_st)
seur_st@assays$spatial@counts[1:6, 1:6]
seur_st@assays$spatial@data[1:6, 1:6]
head(seur_st@meta.data)


# ------------
# run CellTrek
# ------------

out_traint <- CellTrek::traint(
  st_data = seur_st, 
  sc_data = seur_nuc, 
  st_assay = "spatial", 
  sc_assay = "snRNAseq", 
  cell_names = "clusters.glmpcamnn", 
  coord_xy = c("imagerow", "imagecol")
)

# to do: check UMAP plot
# DimPlot(out_traint, group.by = "type")

out_celltrek <- CellTrek::celltrek(
  st_sc_int = out_traint, 
  int_assay = "traint", 
  sc_data = seur_nuc, 
  sc_assay = "scRNAseq", 
  reduction = "pca")$celltrek

out_celltrek$cell_type <- factor(
  out_celltrek$cell_type, 
  levels = sort(unique(out_celltrek$cell_type)))

# to do: add image for visualization
CellTrek::celltrek_vis(
  out_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cell_type:id_new), 
  out_celltrek@images$mysample@image, 
  out_celltrek@images$mysample@scale.factors$lowres)

# to do: export spatial registration results to import back to SPE

