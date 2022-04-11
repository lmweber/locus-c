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
# qrsh -l caracol,mem_free=64G,h_vmem=64G

# module load conda_R/4.1.x
# R


# --------------
# load R objects
# --------------

# load SPE object for use within Python session with reticulate

library(here)
library(SpatialExperiment)

fn_spe <- here("processed_data", "SPE", "LC_qualityControlled.rds")
spe <- readRDS(fn_spe)
dim(spe)

# select LC regions only
spe <- spe[, colData(spe)$annot_region]
dim(spe)

table(colData(spe)$sample_id)
table(colData(spe)$sample_part_id)


# load SCE object (snRNA-seq)

fn_sce <- here("processed_data", "SCE", "sce_updated_LC.rda")
load(fn_sce)

dim(sce.lc)


# ---------------------------------
# extract components from R objects
# ---------------------------------

# extract components from SPE/SCE objects to create new AnnData objects later
# see https://theislab.github.io/scanpy-in-R/, section 4.4.1: "Creating AnnData from SingleCellExperiment"
# note: use standard data.frames instead of DataFrames so reticulate can access them

# SPE object
spe_assays <- assays(spe)
spe_counts <- assay(spe, "counts")
spe_rowdata <- as.data.frame(rowData(spe))
spe_coldata <- as.data.frame(colData(spe))
spe_coldata_combined <- cbind(as.data.frame(colData(spe)), spatialCoords(spe))

# SCE object
sce_assays <- assays(sce.lc)
sce_counts <- assay(sce.lc, "counts")
sce_rowdata <- as.data.frame(rowData(sce.lc))
sce_coldata <- as.data.frame(colData(sce.lc))
sce_reduceddims <- reducedDims(sce.lc)


# --------------------
# start Python session
# --------------------

# start interactive Python session using reticulate with correct Python installation

library(reticulate)
use_python("/users/lweber/miniconda3/bin/python", required = TRUE)  ## skip on laptop
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
#adata_sce.obsm['umap'] = r.sce_reduceddims

# to do: add reducedDims components to adata.sce

# to do: check these objects are using correct data, e.g. counts not logcounts


# -------------------------------
# continue cell2location workflow
# -------------------------------

# to do: continue tutorial workflow below using objects loaded above

# to do: adapt following code below from tutorial


# Loading Visium and scRNA-seq reference data

adata_vis = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]

# rename genes to ENSEMBL
adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var_names = adata_vis.var['gene_ids']
adata_vis.var_names.name = None


# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


# Download data if not already here
# if not os.path.exists('/dcs04/hicks/data/lweber/cell2location/data/sc.h5ad'):
#     !cd /dcs04/hicks/data/lweber/cell2location/data/ && wget https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad


# Read data
adata_ref = sc.read(f'/dcs04/hicks/data/lweber/cell2location/data/sc.h5ad')

# Use ENSEMBL as gene IDs to make sure IDs are unique and correctly matched
adata_ref.var['SYMBOL'] = adata_ref.var.index
adata_ref.var.index = adata_ref.var['GeneID-2'].copy()
adata_ref.var_names = adata_ref.var['GeneID-2'].copy()
adata_ref.var.index.name = None
adata_ref.raw.var['SYMBOL'] = adata_ref.raw.var.index
adata_ref.raw.var.index = adata_ref.raw.var['GeneID-2'].copy()
adata_ref.raw.var.index.name = None


# before we estimate the reference cell type signature we recommend to perform very permissive genes selection
# in this 2D histogram orange rectangle lays over excluded genes.
# In this case, the downloaded dataset was already filtered using this method,
# hence no density under the orange rectangle
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()


# Estimation of reference cell type signatures (NB regression)

# prepare anndata for the regression model
scvi.data.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='Subset',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['Method']
                       )
scvi.data.view_anndata_setup(adata_ref)


# create and train the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# plot ELBO loss history during training, removing first 20 epochs from the plot
mod.plot_history(20)


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


# Examine QC plots

mod.plot_QC()


# The model and output h5ad can be loaded later like this:

mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


# Cell2location: spatial mapping

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
scvi.data.setup_anndata(adata=adata_vis, batch_key="sample")
scvi.data.view_anndata_setup(adata_vis)


# Note: change N_cells_per_location=30 depending on tissue cell density

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
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


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file


# The model and output h5ad can be loaded later like this:

mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)


# Examine reconstruction accuracy to assess if there are any issues with mapping
# the plot should be roughly diagonal, strong deviations will signal problems
mod.plot_QC()

mod.plot_spatial_QC_across_batches()

plt.savefig('plot1.png')


# Visualising cell abundance in spatial coordinates


# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select one slide
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma',
                  # show first 8 cell types
                  color=['B_Cycling', 'B_GC_LZ', 'T_CD4+_TfH_GC', 'FDC',
                         'B_naive', 'T_CD4+_naive', 'B_plasma', 'Endo'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )

plt.savefig('plot2.png')


# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['T_CD4+_naive', 'B_naive', 'FDC']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'V1_Human_Lymph_Node')

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

