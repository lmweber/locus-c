###########################################################
# LC analyses: spot-level deconvolution using cell2location
# Lukas Weber, Apr 2022
###########################################################

# references:
# - cell2location tutorial: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
# - Scanpy AnnData conversion from R: https://theislab.github.io/scanpy-in-R/
# - reticulate: https://rstudio.github.io/reticulate/

# cell2location installation (note: from specific git commit):
# https://github.com/BayraktarLab/cell2location


# The following code is intended to be pasted interactively into a Python REPL 
# session running within an R session (see links above). This allows us to 
# access previously saved R objects and convert these to Python AnnData objects.


# start interactive session on JHPCE GPU node:
# qrsh -l caracol,mem_free=128G,h_vmem=128G

# module load conda_R/4.1.x
# R


# ----------------------
# load Visium SPE object
# ----------------------

# load Visium SPE object for use within Python session with reticulate

library(here)
library(SpatialExperiment)

fn_spe <- here("processed_data", "SPE", "LC_qualityControlled.rds")
spe <- readRDS(fn_spe)
dim(spe)

# remove samples where NE neurons were not captured
samples_remove <- "Br5459_LC_round2"
spe <- spe[, !(colData(spe)$sample_id %in% samples_remove)]

colData(spe)$sample_id <- droplevels(colData(spe)$sample_id)

dim(spe)

# check sample IDs
table(colData(spe)$sample_id)
table(colData(spe)$sample_part_id)


# -------------------------
# load snRNA-seq SCE object
# -------------------------

# load snRNA-seq SCE object for use within Python session with reticulate

fn_sce <- here("processed_data", "SCE", "sce_updated_LC.rda")
load(fn_sce)

dim(sce.lc)

# check sample IDs
table(colData(sce.lc)$Sample)
# check annotated cell type clusters
table(colData(sce.lc)$cellType.merged)

# remove ambiguous cell type clusters ("ambig.lowNTx_A" to "ambig.lowNTx_F")
ambig_names <- paste0("ambig.lowNTx_", LETTERS[1:6])
ix_remove <- colData(sce.lc)$cellType.merged %in% ambig_names
table(ix_remove)
sce.lc <- sce.lc[, !ix_remove]
colData(sce.lc)$cellType.merged <- droplevels(colData(sce.lc)$cellType.merged)
table(colData(sce.lc)$cellType.merged)

dim(sce.lc)


# ---------------------------------
# extract components from R objects
# ---------------------------------

# extract components from SPE/SCE objects to create new AnnData objects later
# see https://theislab.github.io/scanpy-in-R/, section 4.4.1: "Creating AnnData from SingleCellExperiment"
# note: use standard data.frames instead of DataFrames so reticulate can access them

# SPE object
spe_counts <- assay(spe, "counts")
spe_logcounts <- assay(spe, "logcounts")
spe_rowdata <- as.data.frame(rowData(spe))
spe_coldata <- as.data.frame(colData(spe))
spe_coldata_combined <- cbind(as.data.frame(colData(spe)), spatialCoords(spe))

# SCE object
sce_counts <- assay(sce.lc, "counts")
sce_binomial_deviance_residuals <- assay(sce.lc, "binomial_deviance_residuals")
sce_logcounts <- assay(sce.lc, "logcounts")
sce_rowdata <- as.data.frame(rowData(sce.lc))
sce_coldata <- as.data.frame(colData(sce.lc))
sce_GLMPCA_approx <- reducedDim(sce.lc, "GLMPCA_approx")
sce_GLMPCA_MNN <- reducedDim(sce.lc, "GLMPCA_MNN")
sce_glmpca_mnn_50 <- reducedDim(sce.lc, "glmpca_mnn_50")
sce_UMAP <- reducedDim(sce.lc, "UMAP")
sce_TSNE <- reducedDim(sce.lc, "TSNE")


# --------------------
# start Python session
# --------------------

# start interactive Python session using reticulate with correct Python installation

library(reticulate)
use_python("/users/lweber/miniconda3/bin/python", required = TRUE)
reticulate::repl_python()


# ----------------------------
# start cell2location workflow
# ----------------------------

# using Python code from cell2location tutorial
# note requires GPU for reasonable runtime

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns


# define results folders

results_folder = '/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/processed_data/deconvolution/'
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


# ----------------------
# create AnnData objects
# ----------------------

# create AnnData objects from previously loaded R objects using reticulate
# see https://theislab.github.io/scanpy-in-R/, section 4.4.1: "Creating AnnData from SingleCellExperiment"

# note: using raw counts for cell2location

# SPE object
adata_spe = sc.AnnData(
    X = r.spe_counts.T, 
    obs = r.spe_coldata_combined, 
    var = r.spe_rowdata
)

# SCE object
adata_sce = sc.AnnData(
    X = r.sce_counts.T, 
    obs = r.sce_coldata, 
    var = r.sce_rowdata
)
adata_sce.obsm['GLMPCA_approx'] = r.sce_GLMPCA_approx
adata_sce.obsm['GLMPCA_MNN'] = r.sce_GLMPCA_MNN
adata_sce.obsm['glmpca_mnn_50'] = r.sce_glmpca_mnn_50
adata_sce.obsm['UMAP'] = r.sce_UMAP
adata_sce.obsm['TSNE'] = r.sce_TSNE


# -------------------------------
# continue cell2location workflow
# -------------------------------

# continue cell2location workflow with code adapted for objects above


# -------------
# preprocessing
# -------------

# update AnnData Visium object to match structure in cell2location workflow

# rename object
adata_vis = adata_spe
# rename 'sample_id' column
adata_vis.obs['sample'] = adata_vis.obs['sample_id']
# rename 'gene_id' column
adata_vis.var['SYMBOL'] = adata_vis.var['gene_id']

adata_vis.var_names = adata_vis.var['gene_id']
adata_vis.var_names.name = None

# note: mitochondrial genes have already been filtered out
# (otherwise see cell2location tutorial code for how to filter them out here)


# update AnnData snRNA-seq object to match structure in cell2location workflow

# rename object
adata_ref = adata_sce
# use ENSEMBL gene IDs
adata_ref.var['SYMBOL'] = adata_ref.var['gene_id']


# recommended gene filtering from cell2location workflow

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()
adata_ref


# --------------------------------------------
# estimation of reference cell type signatures
# --------------------------------------------

# prepare AnnData for the regression model
scvi.data.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='cellType.merged'
                       )
scvi.data.view_anndata_setup(adata_ref)

# create and train the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# plot ELBO loss history during training, removing first 20 epochs from the plot
mod.plot_history(20)


# in this section, we export the estimated cell abundance (summary of the posterior distribution)
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# save model
mod.save(f"{ref_run_name}", overwrite=True)

# save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


# examine QC plots
mod.plot_QC()


# the model and output h5ad can be loaded later like this:
#mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
#adata_file = f"{ref_run_name}/sc.h5ad"
#adata_ref = sc.read_h5ad(adata_file)


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


# ------------------------------
# Cell2location: spatial mapping
# ------------------------------

# note from tutorial documentation: user needs to choose a value for 'N_cells_per_location'
# based on prior knowledge of tissue cell density in this dataset


# find shared genes and subset both AnnData and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare AnnData for Cell2location model
scvi.data.setup_anndata(adata=adata_vis, batch_key="sample")
scvi.data.view_anndata_setup(adata_vis)


# note: change N_cells_per_location=30 depending on tissue cell density

# create and train the model

# note: using 'detection_alpha=20'
# see https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=3,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=20
)

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);


# in this section, we export the estimated cell abundance (summary of the posterior distribution)
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file


# the model and output h5ad can be loaded later like this:
#mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
#adata_file = f"{run_name}/sp.h5ad"
#adata_vis = sc.read_h5ad(adata_file)


# examine reconstruction accuracy to assess if there are any issues with mapping
# the plot should be roughly diagonal, strong deviations will signal problems
mod.plot_QC()

mod.plot_spatial_QC_across_batches()

plt.savefig('plot1.png')


# --------------------------
# export AnnData object to R
# --------------------------

# export components of AnnData object back to R to add to SpatialExperiment object

obs = adata_vis.obs
var = adata_vis.var
uns = adata_vis.uns
obsm = adata_vis.obsm

# close Python REPL session and go back to R session
exit

# now can access objects using reticulate 'py$' syntax
str(py$obs)
str(py$var)
names(py$uns)
str(py$obsm)
# main results are stored in 'py$obsm' as follows
str(py$obsm['means_cell_abundance_w_sf'])
str(py$obsm['stds_cell_abundance_w_sf'])
str(py$obsm['q05_cell_abundance_w_sf'])
str(py$obsm['q95_cell_abundance_w_sf'])

# using 'q05_cell_abundance_w_sf' following cell2location tutorial

# check results
head(py$obsm['q05_cell_abundance_w_sf'])
dim(py$obsm['q05_cell_abundance_w_sf'])
colnames(py$obsm['q05_cell_abundance_w_sf'])
dim(spe)
all(rownames(py$obsm['q05_cell_abundance_w_sf']) == colnames(spe))
length(rownames(py$obsm['q05_cell_abundance_w_sf']) == colnames(spe))

# cell type abundances (number of cells per spot)
summary(py$obsm['q05_cell_abundance_w_sf'])
# deciles
sapply(py$obsm['q05_cell_abundance_w_sf'], quantile, seq(0, 1, by = 0.1))


# add cell2location results to SpatialExperiment object
colData(spe) <- cbind(colData(spe), py$obsm['q05_cell_abundance_w_sf'])


# save SPE object for further plotting
fn <- here("processed_data", "SPE", "LC_cell2location.rds")
saveRDS(spe, file = fn)


# ---------------------------------------------------------------
# alternatively: plotting in Python using cell2location functions
# ---------------------------------------------------------------

# alternatively: modify cell2location plotting code in next section to create plots in Python

# note: Python code in next section does not work, since 'adata_vis' object does not 
# contain image data in 'spatial' slot


# -------------------------------------------------
# visualizing cell abundance in spatial coordinates
# -------------------------------------------------

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select one slide
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'Br6522_LC_1_round1')  ## to do: update tutorial code here

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  # show first 8 cell types                               ## to do: update tutorial code here
                  color=['B_Cycling', 'B_GC_LZ', 'T_CD4+_TfH_GC', 'FDC',  ## to do: update tutorial code here
                         'B_naive', 'T_CD4+_naive', 'B_plasma', 'Endo'],  ## to do: update tutorial code here
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )

plt.savefig('plot2.png')


# now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters                          ## to do: update tutorial code here
clust_labels = ['T_CD4+_naive', 'B_naive', 'FDC']  ## to do: update tutorial code here
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')  ## to do: update tutorial code here

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )

plt.savefig('plot3.png')

