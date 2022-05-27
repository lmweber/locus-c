library("SingleCellExperiment")
library("lobstr")
library("here")
library("scuttle")
library("sessioninfo")

load(here("processed_data", "SCE", "sce_updated_LC.rda"), verbose = TRUE)

lobstr::obj_size(sce.lc) / 1024 ^ 3
# 5.413584 B

assays(sce.lc)$counts <- NULL
lobstr::obj_size(sce.lc) / 1024 ^ 3
# 4.861437 B
assays(sce.lc)$binomial_deviance_residuals <- NULL
lobstr::obj_size(sce.lc) / 1024 ^ 3
# 0.5954937 B

## Make unique gene names
rownames(sce.lc) <-
  uniquifyFeatureNames(rowData(sce.lc)$gene_id, rowData(sce.lc)$gene_name)

save(
  sce.lc,
  cell_colors.lc,
  file = here("code", "iSEE_apps", "snRNA-seq_app", "sce_lc_small.Rdata")
)



## Code for selecting top genes to show:
x <-
  read.csv(
    here(
      "outputs",
      "06_downstream",
      "pseudobulkDE",
      "pseudobulk_LCvsWM",
      "LC_pseudobulkDE_LCvsWM_topGenes.csv"
    )
  )
paste(x$gene_name[seq_len(10)], collapse = "\n")
