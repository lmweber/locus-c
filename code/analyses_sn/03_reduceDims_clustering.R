### LC snRNA-seq analysis
### Feature selection with deviance & clustering
###     qrsh -l bluejay,mf=
### Initiated: MNT 31Dec2021

library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(scry)
library(batchelor)
library(bluster)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)

here()
    # [1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"

source("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/plotExpressionCustom.R")


### Feature selection with deviance residuals =====

# Load working SCE
load(here("processed_data","SCE", "sce_working_LC.rda"), verbose=T)

sce.lc
    #class: SingleCellExperiment 
    # dim: 36601 15642 
    # metadata(1): Samples
    # assays(1): counts
    # rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
    #   ENSG00000277196
    # rowData names(6): source type ... gene_name gene_type
    # colnames(15642): 3_AAACCCAAGTAGAATC-1 3_AAACCCAAGTGCACCC-1 ...
    #   2_TTTGTTGTCAGGGATG-1 2_TTTGTTGTCCCTCTCC-1
    # colData names(9): Sample Barcode ... high.mito doubletScore
    # reducedDimNames(0):
    # mainExpName: NULL
    # altExpNames(0):

Sys.time()
    #[1] "2021-12-31 14:36:36 EST"
sce.lc <- devianceFeatureSelection(sce.lc,
                                   assay="counts", fam="binomial", sorted=F,
                                      # these are default params btw
                                   batch=as.factor(sce.lc$Sample))
Sys.time()
    #[1] "2021-12-31 14:38:08 EST"

# Btw:
table(is.na(rowData(sce.lc)$binomial_deviance))
    # FALSE  TRUE 
    # 29556  7045

# Observe:
plot(sort(rowData(sce.lc)$binomial_deviance, decreasing=T),
     type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance: LC (n=3)")
abline(v=2000, lty=2, col="red")


## Approximate GLM-PCA with PCA on null deviance residuals
Sys.time()
    #[1] "2021-12-31 15:13:20 EST"
sce.lc <- nullResiduals(sce.lc, assay="counts", fam="binomial",  # default params
                        type="deviance")#, batch=as.factor(sce.lc$Sample))
    # Errored out trying to pass `batch=` arg:
        # Error in res[, idx] <- .null_residuals(m[, idx], fam = fam, type = type,  : 
        # number of items to replace is not a multiple of replacement length
    
    # When leaving that arg out - successfully runs
        # [In addition:] Warning message:
        # In sqrt(x@x) : NaNs produced
Sys.time()
    #[1] "2021-12-31 15:20:56 EST"

    # The above adds a new computed assay called 'binomial_deviance_residuals'


## Take the top 2000 highly-deviant genes ('HDG's) and run PCA in that space
hdgs.lc <- rownames(sce.lc)[order(rowData(sce.lc)$binomial_deviance, decreasing=T)][1:2000]
hdgs.symbols <- rowData(sce.lc)$gene_name[match(hdgs.lc, rowData(sce.lc)$gene_id)]

Sys.time()
    #[1] "2021-12-31 16:20:14 EST"
set.seed(109)
sce.lc <-  runPCA(sce.lc, exprs_values="binomial_deviance_residuals",
                  subset_row=hdgs.lc, ncomponents=100,
                  name="GLMPCA_approx",
                  BSPARAM=BiocSingular::IrlbaParam())
Sys.time()
    #[1] "2021-12-31 16:24:58 EST"

# Visualize
plotReducedDim(sce.lc, dimred="GLMPCA_approx", colour_by="Sample",
               ncomponents=4, point_alpha=0.3, point_size=1.5)


## See how this changes with batchelor::reducedMNN ===
Sys.time()
    #[1] "2021-12-31 17:00:16 EST"
glmpca.mnn <- reducedMNN(reducedDim(sce.lc, "GLMPCA_approx"),
                         batch=as.factor(sce.lc$Sample),
                         merge.order=c("Br8079_LC", "Br2701_LC", "Br6522_LC")
                         )
Sys.time()
    #[1] "2021-12-31 17:00:37 EST"

# Store this
reducedDim(sce.lc, "GLMPCA_MNN") <- glmpca.mnn$corrected

# Vizualize
plotReducedDim(sce.lc, dimred="GLMPCA_MNN", colour_by="Sample",
               ncomponents=4, point_alpha=0.3, point_size=1.5)
    # Doesn't look too amazingly better--at least in the top 4 dims.  Might have to
    # do some clustering in both PC spaces to assess


# Save for now
save(sce.lc, file=here("processed_data","SCE", "sce_working_LC.rda"))




### Clustering ====================
  # Perform graph-based clustering, as in Tran-Maynard, et al. Neuron 2021

  # Just perform on top 50 PCs, since for now we're assessing GLM PCA
  # [and then with MNN on reduced dims]
  # MNT update 27Jan2022: Proceed with the GLM-PCA > MNN approach, based off of 1/13 project mtg
  #                       -> Clean up approximate-only method & fastMNN comparisons


## Make some other reduced dims for visualization ===
sce.test <- sce.lc
    
# Take top 50 PCs to quicken computation
reducedDim(sce.test, "glmpca_mnn_50") <- reducedDim(sce.test, "GLMPCA_MNN")[ ,1:50]
    # oh didn't have to use this - can just use 'n_dimred'


# UMAP
set.seed(109)
sce.test.mnn <- runUMAP(sce.test, dimred="glmpca_mnn_50",
                        name="UMAP")
# t-SNE
set.seed(109)
sce.test.mnn <- runTSNE(sce.test.mnn, dimred="glmpca_mnn_50",
                        name="TSNE")

# For logcounts
sce.test <- multiBatchNorm(sce.test, batch=sce.test$Sample)


## Save these various test iterations for further exploration
save(sce.test, sce.test.approx, sce.test.mnn, hdgs.lc,
     file=here("processed_data","SCE", "sce_reducedDim-tests_LC.rda"))
    


### Clustering (step 1?): Perform graph-based clustering in this optimal PC space ===
  #                    - take k=20 NN to build graph

snn.gr.glmpcamnn <- buildSNNGraph(sce.test.mnn, k=20, use.dimred="glmpca_mnn_50")
clusters.glmpcamnn <- igraph::cluster_walktrap(snn.gr.glmpcamnn)
table(clusters.glmpcamnn$membership)
    ## 60 prelim clusters:

    #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
    #2301   99  178  474 2560  179  189  182   73  151   69  119  114  439 3080   65  271  219 
    #  19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36 
    # 183  123  109 1140  127  237   66   59  159   80  108  137   47   64  104  173  148  162 
    #  37   38   39   40   41   42   43   44   45   46   47   48   49   50   51   52   53   54 
    #  61   57  255   36  113   33   31   53   53   30  137   64   89   31   32   75   95  137 
    #  55   56   57   58   59   60 
    #  39   50   49   46   54   34 



## doubletScore distributions?
cellClust.idx <- splitit(sce.test.mnn$clusters.glmcpcamnn)
sapply(cellClust.idx, function(x){quantile(sce.test.mnn$doubletScore[x])})

## For sake of filtering, compute some median expression info across graph-based clusters
medianNon0.lc <- lapply(cellClust.idx, function(x){
  apply(as.matrix(assay(sce.test, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

head(medianNon0.lc[[1]])
    # Good it kept the gene names

# Save for future reference
medianNon0.glmpcamnn <- medianNon0.lc
save(medianNon0.glmpcamnn, file=here("processed_data","SCE", "medianNon0_Booleans_glmpcamnn-graphClusters.rda"))


# With so many clusters, can we predict neuron vs non-neuronal by marker genes?
    markers.neurons <- c("SYT1", "SLC17A7", "SLC17A6","GAD2")
    markers.neurons <- rowData(sce.test.mnn)$gene_id[match(markers.neurons, rowData(sce.test.mnn)$gene_name)]
    names(markers.neurons) <- c("SYT1", "SLC17A7", "SLC17A6","GAD2")
    
    sapply(medianNon0.glmpcamnn, function(x){
      x[markers.neurons] == TRUE
    })
        # ~half (+) for pan-neuronal marker SYT1 (this usually is captured better than SNAP25)
    neuron.clusters <- which(sapply(medianNon0.glmpcamnn, function(x){ x[markers.neurons["SYT1"]] == TRUE }))
        # So actually 41/60 clusters
    
    table(sce.test.mnn$clusters.glmcpcamnn %in% neuron.clusters)["TRUE"] / ncol(sce.test.mnn) # 0.8554533
         # Looks like ~86% of these nuclei are putatively neuronal
    
    sce.test.mnn$SYT1 <- ifelse(sce.test.mnn$clusters.glmcpcamnn %in% neuron.clusters, "neuronal", "nonNeuronal")
    
    

    # Cluster x donor distribution:
    table(sce.neuron.mnn$clusters.glmcpcamnn, sce.neuron.mnn$Sample)
    # ====
    #    Br2701_LC Br6522_LC Br8079_LC
    # 1        791       685       825
    # 2         16        45        38
    # 3         77         1       100
    # 5        556       120      1884
    # 8         26        66        90
    # 9          3        57        13
    # 10         4         3       144
    # 12        13        27        79
    # 13         9         0       105
    # 14       127       257        55
    # 15       546       460      2074
    # 17       178         0        93
    # 18        86        31       102
    # 19        50         0       133
    # 20        27         2        94
    # 22       874        85       181
    # 23         0       127         0
    # 26         1        32        26
    # 27        43         6       110
    # 28         0        17        63
    # 31        14        28         5  **
    # 35        85         9        54
    # 36        38         1       123  *
    # 39        13         0       242
    # 40        12         0        24  ***
    # 41         2       111         0
    # 43        23         0         8
    # 47        52         0        85
    # 48        56         3         5
    # 49        65         0        24
    # 50        19         0        12
    # 51        15         3        14
    # 52         0         0        75
    # 53         0         0        95
    # 54        59         1        77
    # 55        22         0        17
    # 56        18         1        31
    # 57         0        49         0
    # 58         0        46         0
    # 59         0         1        53
    # 60        21         1        12
    
    # ====
    
    
    ## Broad markers of interest:
    markers.mathys.tran = list(
      'neuron' = c('SYT1', 'SNAP25', 'GRIN1'),
      'excit_neuron' = c('CAMK2A', 'NRGN','SLC17A7', 'SLC17A6', 'SLC17A8'),
      'inhib_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
      # Norepinephrine & serotonergic markers
      'neuron.NE' = c("TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC"), #SLC6A3 - saw no DAT
      'neuron.5HT' = c("SLC6A4", "TPH1", "TPH2", "DDC"),
                        # SERT, serotonin T (aka 5-HTT); 
      'monoamine.metab' = c("COMT", "MAOA", "MAOB"),
      # MSN markers
      'MSNs.pan' = c("PPP1R1B","BCL11B"),# "CTIP2")
      'MSNs.D1' = c("DRD1", "PDYN", "TAC1"),
      'MSNs.D2' = c("DRD2", "PENK"),
      ## Non-neuronal:
      'oligodendrocyte' = c('MBP', 'MOBP', 'PLP1'),
      'oligo_precursor' = c('PDGFRA', 'VCAN', 'CSPG4'),
      'microglia' = c('CD74', 'CSF1R', 'C3'),
      'astrocyte' = c('GFAP', 'TNC', 'AQP4', 'SLC1A2'),
      'endothelial' = c('CLDN5', 'FLT1', 'VTN'),
      # Post-hoc from Tran-Maynard, et al. Neuron 2021
      'differn_committed_OPC' = c("SOX4", "BCAN", "GPR17", "TNS3"),
      'Tcell' = c('SKAP1', 'ITK', 'CD247'),
      'Mural' = c('COL1A2', 'TBX18', 'RBPMS'),
      'Macro' = c('CD163', 'SIGLEC1', 'F13A1')
    )

    ## Since we're happy with the  approx. GLMPCA>MNN approach, just use the multiBatchNorm 
     #    logcounts from the fastMNN input:
    assay(sce.test.mnn, "logcounts") <- assay(sce.test, "logcounts")
    
    
    ## And expression violin plots - subset those neuronal
    sce.nonNeuron.mnn <- sce.test.mnn[ ,sce.test.mnn$SYT1=="nonNeuronal"]
    sce.nonNeuron.mnn$clusters.glmcpcamnn <- droplevels(sce.nonNeuron.mnn$clusters.glmcpcamnn)
    
    # Change names of features so they're more interpretable
    rownames(sce.nonNeuron.mnn) <- rowData(sce.nonNeuron.mnn)$gene_name
    
    pdf(here("plots","snRNA-seq",paste0("LC-n3_expression-violin_markers-nonNeuronal_GLMPCA-MNN-graphClusters.pdf")),
        height=6, width=12)
    for(i in 1:length(markers.mathys.tran)){
      print(
        #plotExpressionCustom(sce = sce.nonNeuron.mnn,
        plotExpressionCustom(sce = sce.nonNeuron,
                             exprs_values = "logcounts",
                             features = markers.mathys.tran[[i]], 
                             features_name = names(markers.mathys.tran)[[i]], 
                             #anno_name = "clusters.glmcpcamnn",
                             anno_name = "cellType",
                             ncol=2, point_alpha=0.4, point_size=0.9,
                             scales="free_y") +  
          ggtitle(label=paste0("LC-n3 24 non-neuronal (/60) clusters: ",
                               names(markers.mathys.tran)[[i]], " markers")) +
          theme(plot.title = element_text(size = 12),
                axis.text.x = element_text(size=7))
      )
    }
    dev.off()


# MNT 27Jan2022 add: once happy, clean up:
sce.lc <- sce.test.mnn
sce.lc$kmeans.13 <- NULL
sce.lc$clusters.fastmnn <- NULL
# Reduced dims
reducedDim(sce.lc, "glmpca_approx_50") <- NULL
# Not sure why had to re-name to something longer - change back:
medianNon0.lc <- medianNon0.glmpcamnn

# Then save
save(sce.lc, medianNon0.lc, hdgs.lc,
     file=here("processed_data","SCE", "sce_updated_LC.rda"))




### Annotate clusters =========
## Re-print more neuronal markers for the 41 putative neuronal clusters
sce.neuron <- sce.lc

sce.neuron <- sce.lc[ ,sce.lc$SYT1=="neuronal"]
sce.neuron$clusters.glmcpcamnn <- droplevels(sce.neuron$clusters.glmcpcamnn)

# Change names of features so they're more interpretable
rownames(sce.neuron) <- rowData(sce.neuron)$gene_name

pdf(here("plots","snRNA-seq",paste0("LC-n3_expression_broad-markers_Neuronal-GLMPCA-MNN-clusters.pdf")),
    height=6, width=14)
for(i in 1:length(markers.mathys.tran)){
  print(
    plotExpressionCustom(sce = sce.neuron,
                         exprs_values = "logcounts",
                         features = markers.mathys.tran[[i]], 
                         features_name = names(markers.mathys.tran)[[i]], 
                         #anno_name = "clusters.glmcpcamnn",
                         anno_name = "cellType",
                         ncol=2, point_alpha=0.4, scales="free_y") +  
      ggtitle(label=paste0("LC-n3 36 neuronal (/60) clusters: ",
                           names(markers.mathys.tran)[[i]], " markers")) +
      theme(plot.title = element_text(size = 12),
            axis.text.x = element_text(size=7))
  )
}
dev.off()


## Some definitely aren't actually neuronal; 15 is elusive (very large - 3000+)
cellClust.idx <- splitit(sce.lc$clusters.glmcpcamnn)
sapply(cellClust.idx, function(x){quantile(sce.lc$sum[x])})
    # Oh it seems to be driven by low n transcripts - probably mostly debris.
    # (maybe can see evidence of this in the emptyDrops scores...)

sapply(cellClust.idx, function(x){round(quantile(sce.lc$doubletScore[x]),2)})


# Re-assign and re-print the above neuronal & non-neuronal violin plots
sce.lc$neuron <- sce.lc$SYT1
sce.lc$neuron[sce.lc$clusters.glmcpcamnn %in% c(12,15, 22,26,28)] <- "nonNeuronal"
    # Fixed proportion neuronal: 0.5691727
    # (this is much closer to what was sorted)

sce.neuron <- sce.lc[ ,sce.lc$neuron=="neuronal"]
sce.neuron$clusters.glmcpcamnn <- droplevels(sce.neuron$clusters.glmcpcamnn)

sce.nonNeuron.mnn <- sce.lc[ ,sce.lc$neuron=="nonNeuronal"]
sce.nonNeuron.mnn$clusters.glmcpcamnn <- droplevels(sce.nonNeuron.mnn$clusters.glmcpcamnn)

# Change names of features so they're more interpretable
rownames(sce.nonNeuron.mnn) <- rowData(sce.nonNeuron.mnn)$gene_name
rownames(sce.neuron) <- rowData(sce.neuron)$gene_name



annotationTab.lc <- data.frame(cluster=c(1:60))
annotationTab.lc$cellType <- NA
annotationTab.lc$cellType[c(1,9)] <- paste0("Neuron.mixed_", c("A","B"))
annotationTab.lc$cellType[c(5,14,17,19, 36,47,52, 53)] <-
  paste0("Neuron.ambig_", c("A","B","C","D", "E","F","G","H"))

annotationTab.lc$cellType[c(3,8,10,13,18,20,41,43,49:51, 55,59,60)] <-
  paste0("Excit_", c("A","B","C","D", "E","F","G","H", "I","J","K","L","M","N"))
annotationTab.lc$cellType[c(2,23,27,35, 39, 54, 56:58)] <-
  paste0("Inhib_", c("A","B","C","D", "E","F","G","H", "I"))

annotationTab.lc$cellType[c(40)] <- "Neuron.NE"
annotationTab.lc$cellType[c(31)] <- "Neuron.5HT"
annotationTab.lc$cellType[c(48)] <- "Neuron.5HT_noDDC"

# Non-neuronal classes
annotationTab.lc$cellType[c(4,7,16,21,22,24, 30,34,37)] <-
  paste0("Oligo_", c("A","B","C","D", "E","F","G","H", "I"))
annotationTab.lc$cellType[c(12,32,45)] <- paste0("OPC_", c("A","B","C"))
annotationTab.lc$cellType[c(25,26,28)] <- paste0("Micro_", c("A","B","C"))
annotationTab.lc$cellType[46] <- c("Macrophage")
annotationTab.lc$cellType[c(6,29, 38,42,44)] <-
  paste0("Astro_", c("A","B","C","D","E"))
annotationTab.lc$cellType[c(11,33)] <- paste0("Endo.Mural_", c("A","B"))
annotationTab.lc$cellType[15] <- c("ambig.lowNTx")



sce.lc$cellType <- annotationTab.lc$cellType[match(sce.lc$clusters.glmpcamnn,
                                                   annotationTab.lc$cluster)]
sce.lc$cellType <- factor(sce.lc$cellType)

options(width=100)
table(sce.lc$cellType)
    #    ambig.lowNTx          Astro_A          Astro_B          Astro_C          Astro_D 
    #            3080              179              108               57               33 
    #         Astro_E     Endo.Mural_A     Endo.Mural_B          Excit_A          Excit_B 
    #              53               69              104              178              182 
    #         Excit_C          Excit_D          Excit_E          Excit_F          Excit_G 
    #             151              114              219              123              113 
    #         Excit_H          Excit_I          Excit_J          Excit_K          Excit_L 
    #              31               89               31               32               39 
    #         Excit_M          Excit_N          Inhib_A          Inhib_B          Inhib_C 
    #              54               34               99              127              159 
    #         Inhib_D          Inhib_E          Inhib_F          Inhib_G          Inhib_H 
    #             148              255              137               50               49 
    #         Inhib_I       Macrophage          Micro_A          Micro_B          Micro_C 
    #              46               30               66               59               80 
    #      Neuron.5HT Neuron.5HT_noDDC   Neuron.ambig_A   Neuron.ambig_B   Neuron.ambig_C 
    #              47               64             2560              439              271 
    #  Neuron.ambig_D   Neuron.ambig_E   Neuron.ambig_F   Neuron.ambig_G   Neuron.ambig_H 
    #             183              162              137               75               95 
    #  Neuron.mixed_A   Neuron.mixed_B        Neuron.NE          Oligo_A          Oligo_B 
    #            2301               73               36              474              189 
    #         Oligo_C          Oligo_D          Oligo_E          Oligo_F          Oligo_G 
    #              65              109             1140              237              137 
    #         Oligo_H          Oligo_I            OPC_A            OPC_B            OPC_C 
    #             173               61              119               64               53 

# Fix typo
sce.lc$clusters.glmpcamnn <- sce.lc$clusters.glmcpcamnn
sce.lc$clusters.glmcpcamnn <- NULL
colData(sce.lc) <- colData(sce.lc)[ ,c(1:9, 13, 10:12)]

# Then save
save(sce.lc, medianNon0.lc, hdgs.lc, annotationTab.lc,
     file=here("processed_data","SCE", "sce_updated_LC.rda"))

# Also save graph
save(snn.gr.glmpcamnn, clusters.glmpcamnn,
     file=here("processed_data","SCE","graph_communities_glmpcamnn_LC.rda"))


## Re-print broad marker plots with the cluster annotations
sce.neuron <- sce.lc[ ,sce.lc$neuron == "neuronal"]
sce.neuron$cellType <- droplevels(sce.neuron$cellType)

sce.nonNeuron <- sce.lc[ ,sce.lc$neuron=="nonNeuronal"]
sce.nonNeuron$cellType <- droplevels(sce.nonNeuron$cellType)

# Change names of features so they're more interpretable
rownames(sce.nonNeuron) <- rowData(sce.nonNeuron)$gene_name
rownames(sce.neuron) <- rowData(sce.neuron)$gene_name

    # --> re-print from above




## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
#[1] "2022-04-01 14:55:46 EDT"
proc.time()
#     user    system   elapsed 
# 1222.625    67.289 23592.692 
options(width = 120)
session_info()
# ─ Session info ──────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-04-01
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# batchelor            * 1.10.0   2021-10-26 [1] Bioconductor
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
# BiocParallel         * 1.28.3   2021-12-09 [2] Bioconductor
# BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# bluster              * 1.4.0    2021-10-26 [2] Bioconductor
# cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
# cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# dplyr                  1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# DropletUtils         * 1.14.2   2022-01-09 [2] Bioconductor
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# gridExtra            * 2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# igraph                 1.3.0    2022-04-01 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.1   2022-02-17 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# metapod                1.2.0    2021-10-26 [2] Bioconductor
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# ResidualMatrix         1.4.0    2021-10-26 [1] Bioconductor
# rhdf5                  2.38.1   2022-03-10 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.4   2022-03-24 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scater               * 1.22.0   2021-10-26 [2] Bioconductor
# scran                * 1.22.1   2021-11-14 [2] Bioconductor
# scry                 * 1.6.0    2021-10-26 [2] Bioconductor
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.4.0    2022-03-30 [2] CRAN (R 4.1.2)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/ntranngu/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ─────────────────────────────────────────────────────────────────────────