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



## Some marker expression ===
ne.markers <- c("TH","DBH", "SLC6A2", "SLC18A2", "DDC", "GCH1", "MAOA", "MAOB", "COMT",
                "SLC6A3", "SLC6A4", "TPH1", "TPH2")
              # DAT, dopamine transporter; SERT, serotonin T (aka 5-HTT); 
                                    # TPH1/2 catalyze the first step in converting tryptophan > 5-HT
ne.markers <- rowData(sce.test.mnn)$gene_id[match(ne.markers, rowData(sce.test.mnn)$gene_name)]
names(ne.markers) <-  c("TH","DBH", "SLC6A2", "SLC18A2", "DDC", "GCH1", "MAOA", "MAOB", "COMT",
                        "SLC6A3", "SLC6A4", "TPH1", "TPH2")

## approx. GLMPCA>MNN - just use the multiBatchNorm logcounts from the fastMNN input:
assay(sce.test.mnn, "logcounts") <- assay(sce.test, "logcounts")

# Exploratory:
    # pdf(here("plots","snRNA-seq",paste0("LC-n3_expression_markers-NEneuron_GLMPCA-MNN_UMAP-graphClusters.pdf")), height=5, width=5)
    # for(i in 1:length(ne.markers)){
    #   print(
    #     plotReducedDim(sce.test.mnn, dimred="UMAP",
    #                    colour_by=ne.markers[i], by_exprs_values="logcounts",
    #                    point_alpha=0.2, point_size=1.2, theme_size=7) +
    #       ggtitle(paste0("Clustering on approx. GLM-PCA-MNN PCs \nLog counts for ", ne.markers[i],": ",names(ne.markers)[i]))
    #   )
    # }
    # dev.off()


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
    
    
    
    # Re-print with the above NE neuron markers
    markers.neurons <- c(markers.neurons, ne.markers)
    
    pdf(here("plots","snRNA-seq",paste0("LC-n3_expression_markers-NEneuron_GLMPCA-MNN_UMAP-graphClusters.pdf")), height=5, width=5)
    # Cluster annotations
    print(
      plotReducedDim(sce.test.mnn, dimred="UMAP", colour_by="clusters.glmcpcamnn",
                   point_alpha=0.2, point_size=1.2,
                   text_by="clusters.glmcpcamnn", text_size=2, theme_size=8) +
      ggtitle(paste0("UMAP on GLM-PCA-MNN (top 50 PCs)"))
      )
    # Putative neuronal clusters
    print(
      plotReducedDim(sce.test.mnn, dimred="UMAP", colour_by="SYT1",
                      point_alpha=0.2, point_size=1.2, theme_size=8)
      )
    # Neuronal markers
    for(i in 1:length(markers.neurons)){
      print(
        plotReducedDim(sce.test.mnn, dimred="UMAP",
                       colour_by=markers.neurons[i], by_exprs_values="logcounts",
                       point_alpha=0.2, point_size=1.2, theme_size=8) +
          ggtitle(paste0("Clustering on approx. GLM-PCA-MNN PCs \nLog counts for ",
                         markers.neurons[i],": ",names(markers.neurons)[i]))
      )
    }
    dev.off()
    
    
    ## And expression violin plots - subset those neuronal
    sce.neuron.mnn <- sce.test.mnn[ ,sce.test.mnn$SYT1=="neuronal"]
    sce.neuron.mnn$clusters.glmcpcamnn <- droplevels(sce.neuron.mnn$clusters.glmcpcamnn)

    # Change names of features so they're more interpretable
    rownames(sce.neuron.mnn) <- rowData(sce.neuron.mnn)$gene_name
    
    pdf(here("plots","snRNA-seq",paste0("LC-n3_expression-violin_markers-NEneuron_GLMPCA-MNN-graphClusters.pdf")),
        height=10, width=12)
    plotExpressionCustom(sce = sce.neuron.mnn,
                         exprs_values = "logcounts",
                         features = names(markers.neurons),
                         #features = c("SLC6A2", "SLC6A3", "SLC6A4", "TPH1", "TPH2"),
                         features_name = "",
                         anno_name = "clusters.glmcpcamnn",
                         ncol=2, point_alpha=0.4, scales="free_y") +  
      ggtitle(label=paste0("LC-n3 41 neuronal (/60) clusters: NE-neuron markers")) +
      theme(plot.title = element_text(size = 12),
            axis.text.x = element_text(size=7))
    dev.off()
    
        # Looks like [GLM-PCA>MNN] clusters 31 &/or 40 are the putative NE neurons (maybe some of 36)
        # Indeed these come from the two donors with more evidence of containing NE neurons:
        
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

    ## And expression violin plots - subset those neuronal
    sce.nonNeuron.mnn <- sce.test.mnn[ ,sce.test.mnn$SYT1=="nonNeuronal"]
    sce.nonNeuron.mnn$clusters.glmcpcamnn <- droplevels(sce.nonNeuron.mnn$clusters.glmcpcamnn)
    
    # Change names of features so they're more interpretable
    rownames(sce.nonNeuron.mnn) <- rowData(sce.nonNeuron.mnn)$gene_name
    
    pdf(here("plots","snRNA-seq",paste0("LC-n3_expression-violin_markers-nonNeuronal_GLMPCA-MNN-graphClusters.pdf")),
        height=6, width=12)
    for(i in 1:length(markers.mathys.tran)){
      print(
        plotExpressionCustom(sce = sce.nonNeuron.mnn,
                             exprs_values = "logcounts",
                             features = markers.mathys.tran[[i]], 
                             features_name = names(markers.mathys.tran)[[i]], 
                             anno_name = "clusters.glmcpcamnn",
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
                         anno_name = "clusters.glmcpcamnn",
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










## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
#
proc.time()
#     user    system   elapsed 
# 
options(width = 120)
session_info()
# 

