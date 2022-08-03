
library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")

load("sce_lc_small.Rdata", verbose = TRUE)

source("initial.R", print.eval = TRUE)

## Fix donor
sce.lc$Donor <- gsub("_LC", "", sce.lc$Sample)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
colData(sce.lc) <- cbind(
  colData(sce.lc)[, !colnames(colData(sce.lc)) %in% c("Donor", "cellType.merged")],
  colData(sce.lc)[, c("cellType.merged", "Donor")]
)

sce.lc$Donor <- as.factor(sce.lc$Donor)

sce.lc <- registerAppOptions(sce.lc, color.maxlevels = length(cell_colors.lc))

iSEE(sce.lc,
     appTitle = "snRNA-seq LC 2022",
     initial = initial,
     colormap = ExperimentColorMap(colData = list(
       Donor = function(n) {
         cols <- paletteer::paletteer_d(
                palette = "RColorBrewer::Dark2",
                n = length(unique(sce.lc$Donor))
            )
         cols <- as.vector(cols)
         names(cols) <- levels(sce.lc$Donor)
         return(cols)
       },
       cellType.merged = function(n) {
         return(cell_colors.lc)
       }
     ))
)