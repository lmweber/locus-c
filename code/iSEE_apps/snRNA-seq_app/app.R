
library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")
library("SpatialExperiment")

sce <- readRDS("LC_singleNucleus_SCE_EHub.rds")

## Make unique gene names
rownames(sce) <-
  uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)

source("initial.R", print.eval = TRUE)

## Fix donor
sce$Donor <- gsub("_LC", "", sce$Sample)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
colData(sce) <- cbind(
  colData(sce)[, !colnames(colData(sce)) %in% c("Donor", "label_merged")],
  colData(sce)[, c("label_merged", "Donor")]
)

sce$Donor <- as.factor(sce$Donor)
sce$label_merged <- as.factor(sce$label_merged)

sce <- registerAppOptions(sce, color.maxlevels = length(levels(sce$label_merged)))

pal_clusters <- c(
  excitatory = "#1F77B4", 
  inhibitory = "#AEC7E8", 
  neurons_ambiguous = "gray60", 
  NE = "#D62728", 
  `5HT` = "#9467BD", 
  astrocytes = "#FF7F0E", 
  endothelial_mural = "#98DF8A", 
  macrophages_microglia = "#8C564B", 
  oligodendrocytes = "#9EDAE5", 
  OPCs = "#17BECF")

iSEE(sce,
     appTitle = "snRNA-seq LC 2022",
     initial = initial,
     colormap = ExperimentColorMap(colData = list(
       Donor = function(n) {
         cols <- paletteer::paletteer_d(
                palette = "RColorBrewer::Dark2",
                n = length(unique(sce$Donor))
            )
         cols <- as.vector(cols)
         names(cols) <- levels(sce$Donor)
         return(cols)
       },
       label_merged = function(n) {
         return(pal_clusters)
       }
     ))
)

