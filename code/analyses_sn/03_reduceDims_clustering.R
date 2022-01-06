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


## Make some other reduced dims for visualization ===

    ## Testing phase ====
    sce.test <- sce.lc
    
    # Take top 50 PCs for both to quicken computation
    reducedDim(sce.test, "glmpca_approx_50") <- reducedDim(sce.test, "GLMPCA_approx")[ ,1:50]
    reducedDim(sce.test, "glmpca_mnn_50") <- reducedDim(sce.test, "GLMPCA_MNN")[ ,1:50]
        # oh didn't have to use this - can just use 'n_dimred'
    
    ## Without MNN on GLMPCA ===
    # UMAP
    set.seed(109)
    sce.test.approx <- runUMAP(sce.test, dimred="glmpca_approx_50",
                               name="UMAP")
    # t-SNE
    set.seed(109)
    sce.test.approx <- runTSNE(sce.test.approx, dimred="glmpca_approx_50",
                               name="TSNE")
    
    # k-means clustering with k=13
    #   (== average of Cell Ranger's sample-specific n graph-based clusters)
    set.seed(109)
    clust.kmeans <- clusterCells(sce.test.approx, use.dimred="glmpca_approx_50", 
                                 BLUSPARAM=KmeansParam(centers=13))
    table(clust.kmeans)
        #   1    2    3    4    5    6    7    8    9   10   11   12   13 
        # 196  843  635  173  656  261  194 1187  725 1234 4177 4442  919
    sce.test.approx$kmeans.13 <- clust.kmeans
    
    ## WITH MNN on GLMPCA ===
    # UMAP
    set.seed(109)
    sce.test.mnn <- runUMAP(sce.test, dimred="glmpca_mnn_50",
                            name="UMAP")
    # t-SNE
    set.seed(109)
    sce.test.mnn <- runTSNE(sce.test.mnn, dimred="glmpca_mnn_50",
                            name="TSNE")
    
    # k-means (with k=13)
    set.seed(109)
    clust.kmeans <- clusterCells(sce.test.approx, use.dimred="glmpca_mnn_50", 
                                 BLUSPARAM=KmeansParam(centers=13))
    table(clust.kmeans)
        #   1    2    3    4    5    6    7    8    9   10   11   12   13 
        #1064  259  654 1674  637  268  190  121  973 1245  179 8103  275 
    sce.test.mnn$kmeans.13 <- clust.kmeans

    ## Add Tran-Maynard, et al. method (just use HDGs from above) =====
        sce.test <- multiBatchNorm(sce.test, batch=sce.test$Sample)
        set.seed(109)
        mnn.hold <-  fastMNN(sce.test, batch=as.factor(sce.test$Sample),
                             merge.order=c("Br8079_LC", "Br2701_LC", "Br6522_LC"),
                             subset.row=hdgs.lc, d=100,
                             correct.all=TRUE, get.variance=TRUE,
                             BSPARAM=BiocSingular::IrlbaParam())

        table(colnames(mnn.hold) == colnames(sce.test))  # all TRUE
        table(mnn.hold$batch == sce.test$Sample) # all TRUE
        
        # Add them to the SCE
        reducedDim(sce.test, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
        metadata(sce.test) <- metadata(mnn.hold)
        
        # UMAP
        set.seed(109)
        sce.test <- runUMAP(sce.test, dimred="PCA_corrected",
                            n_dimred=50, name="UMAP")
        # t-SNE
        set.seed(109)
        sce.test <- runTSNE(sce.test, dimred="PCA_corrected",
                                   n_dimred=50, name="TSNE")
        
        # k-means (with k=13)
        reducedDim(sce.test, "PCA_corrected_50") <- reducedDim(sce.test, "PCA_corrected")[ ,1:50]
        set.seed(109)
        clust.kmeans <- clusterCells(sce.test, use.dimred="PCA_corrected_50",
                                     BLUSPARAM=KmeansParam(centers=13))
        table(clust.kmeans)
            #   1    2    3    4    5    6    7    8    9   10   11   12   13 
            #2713  246  647  580  600  598  522 2091 4305  701  935 1292  412
        sce.test$kmeans.13 <- clust.kmeans
                
        # end added Tran-Maynard et al. chunk ====

    
    ## How do these look?
    pdf(here("plots","snRNA-seq","LC-n3_reducedDims_GLMPCA_reducedMNN-test.pdf"), height=5, width=5)
        ## TSNE
        plotReducedDim(sce.test.approx, dimred="TSNE", colour_by="Sample",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle("t-SNE on GLMPCA (top 50)")
        plotReducedDim(sce.test.mnn, dimred="TSNE", colour_by="Sample",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle("t-SNE on GLMPCA-MNN (top 50)")
        # Tran-Maynard
        plotReducedDim(sce.test, dimred="TSNE", colour_by="Sample",
                       point_alpha=0.2, point_size=0.8) +
            ggtitle("t-SNE on fastMNN-corrected (top 50) PCs\nfrom log-norm. counts")
        
        ## UMAP
        plotReducedDim(sce.test.approx, dimred="UMAP", colour_by="Sample",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on GLMPCA (top 50)")
        
        plotReducedDim(sce.test.mnn, dimred="UMAP", colour_by="Sample",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on GLMPCA-MNN (top 50)")
        # Tran-Maynard
        plotReducedDim(sce.test, dimred="UMAP", colour_by="Sample",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle("UMAP on fastMNN-corrected (top 50) PCs\nfrom log-norm. counts")
        
    dev.off()
    
    
    ## A more full version with clusters and other metrics
    list.sces <- list(GLMPCA_approx = sce.test.approx,
                      GLMPCA_MNN = sce.test.mnn,
                      fastMNN = sce.test)
    
    for (i in names(list.sces)){
            # For some reason this is creating broken pdfs - just do manually (e.g. `i <- names(list.sces)[1]`)
        pdf(here("plots","snRNA-seq",paste0("LC-n3_reducedDims_",i,"_kmeans13.pdf")), height=5, width=5)
        ## TSNE
        plotReducedDim(list.sces[[i]], dimred="TSNE", colour_by="Sample",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle(paste0("t-SNE on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="TSNE", colour_by="kmeans.13",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle(paste0("t-SNE on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="TSNE", colour_by="detected",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle(paste0("t-SNE on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="TSNE", colour_by="doubletScore",
                       point_alpha=0.3, point_size=0.8) +
            ggtitle(paste0("t-SNE on ", i, " (top 50 PCs)"))
        
        ## UMAP
        plotReducedDim(list.sces[[i]], dimred="UMAP", colour_by="Sample",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle(paste0("UMAP on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="UMAP", colour_by="kmeans.13",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle(paste0("UMAP on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="UMAP", colour_by="detected",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle(paste0("UMAP on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="UMAP", colour_by="doubletScore",
                       point_alpha=0.2, point_size=1.5) +
            ggtitle(paste0("UMAP on ", i, " (top 50 PCs)"))
        dev.off()
    }
    
    
    
    
    ## Save these various test iterations for further exploration
    save(sce.test, sce.test.approx, sce.test.mnn, hdgs.lc,
         file=here("processed_data","SCE", "sce_reducedDim-tests_LC.rda"))
    
    ### END TESTING ========


### Clustering step 1: Perform graph-based clustering in this optimal PC space ===
  #                  - take k=20 NN to build graph

## Testing phase ====
    snn.gr.glmpcamnn <- buildSNNGraph(sce.test, k=20, use.dimred="glmpca_mnn_50")
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
    
    snn.gr.fastmnn <- buildSNNGraph(sce.test, k=20, use.dimred="PCA_corrected_50")
    clusters.fastmnn <- igraph::cluster_walktrap(snn.gr.fastmnn)
    table(clusters.fastmnn$membership)
        # 55 prelim clusters:
        #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
        # 134  172 2383 2405  572 2304  107  366  177  123  102  489  591  144  143  127  134   87 
        #  19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36 
        # 220  511   97   46   80   86   47   99   85  211   70   83   57   96  102  152  103   69 
        #  37   38   39   40   41   42   43   44   45   46   47   48   49   50   51   52   53   54 
        # 178  148  167  201   68 1298   50   38   47   54  128  125   72   53   27   36   28   50 
        #  55 
        # 100 
    
    options(width=120)
    print(tail(table(fastMNN = clusters.fastmnn$membership, glmpca_MNN = clusters.glmpcamnn$membership)[ ,45:60],
               n=30), zero.print=".")
        # Decently similar/concordant
    pairwiseRand(sce.test$clusters.glmcpcamnn, sce.test$clusters.fastmnn, mode="index")
        #[1] 0.4555515
      

    # Store both ('sce.test' has the TSNE/UMAP from fastMNN approach)
    sce.test.mnn$clusters.glmcpcamnn <- factor(clusters.glmpcamnn$membership)
    sce.test.mnn$clusters.fastmnn <- factor(clusters.fastmnn$membership)
    
    sce.test$clusters.glmcpcamnn <- factor(clusters.glmpcamnn$membership)
    sce.test$clusters.fastmnn <- factor(clusters.fastmnn$membership)
    
    # Print
    list.sces <- list(clusters.glmcpcamnn = sce.test.mnn,
                      clusters.fastmnn = sce.test)
    
    for (i in names(list.sces)){
      # For some reason this is creating broken pdfs - just do manually (e.g. `i <- names(list.sces)[1]`)
      pdf(here("plots","snRNA-seq",paste0("LC-n3_reducedDims_",i,"_graphClusters.pdf")), height=5, width=5)
        ## UMAP
        plotReducedDim(list.sces[[i]], dimred="UMAP", colour_by="Sample",
                       point_alpha=0.2, point_size=1.2) +
          ggtitle(paste0("UMAP on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="UMAP", colour_by=i,
                       point_alpha=0.2, point_size=1.2,
                       text_by=i, text_size=2) +
          ggtitle(paste0("UMAP on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="UMAP", colour_by="detected",
                       point_alpha=0.2, point_size=1.2) +
          ggtitle(paste0("UMAP on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="UMAP", colour_by="doubletScore",
                       point_alpha=0.2, point_size=1.2) +
          ggtitle(paste0("UMAP on ", i, " (top 50 PCs)"))
        
        ## TSNE
        plotReducedDim(list.sces[[i]], dimred="TSNE", colour_by="Sample",
                       point_alpha=0.3, point_size=0.8) +
          ggtitle(paste0("t-SNE on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="TSNE", colour_by=i,
                       point_alpha=0.3, point_size=0.8,
                       text_by=i, text_size=2) +
          ggtitle(paste0("t-SNE on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="TSNE", colour_by="detected",
                       point_alpha=0.3, point_size=0.8) +
          ggtitle(paste0("t-SNE on ", i, " (top 50 PCs)"))
        plotReducedDim(list.sces[[i]], dimred="TSNE", colour_by="doubletScore",
                       point_alpha=0.3, point_size=0.8) +
          ggtitle(paste0("t-SNE on ", i, " (top 50 PCs)"))
      
      dev.off()
    }
    

    ## Save these various test iterations for further exploration
    save(sce.test, sce.test.approx, sce.test.mnn, hdgs.lc,
         file=here("processed_data","SCE", "sce_reducedDim-tests_LC.rda"))

    
    
    
    ## Some marker expression ===
    ne.markers <- c("TH","DBH", "SLC6A2", "SLC18A2", "DDC", "GCH1")
    ne.markers <- rowData(sce.test)$gene_id[match(ne.markers, rowData(sce.test)$gene_name)]
    names(ne.markers) <-  c("TH","DBH", "SLC6A2", "SLC18A2", "DDC", "GCH1")
    
    ## approx. GLMPCA>MNN - just use the multiBatchNorm logcounts from the fastMNN input:
    assay(sce.test.mnn, "logcounts") <- assay(sce.test, "logcounts")
    
    pdf(here("plots","snRNA-seq",paste0("LC-n3_expression_markers-NEneuron_GLMPCA-MNN_UMAP-graphClusters.pdf")), height=5, width=5)
    for(i in 1:length(ne.markers)){
      print(
        plotReducedDim(sce.test.mnn, dimred="UMAP",
                       colour_by=ne.markers[i], by_exprs_values="logcounts",
                       point_alpha=0.2, point_size=1.2, theme_size=7) +
          ggtitle(paste0("Clustering on approx. GLM-PCA-MNN PCs \nLog counts for ", ne.markers[i],": ",names(ne.markers)[i]))
      )
    }
    dev.off()
    
    ## fastMNN
    pdf(here("plots","snRNA-seq",paste0("LC-n3_expression_markers-NEneuron_fastMNN-UMAP-graphClusters.pdf")), height=5, width=5)
    for(i in 1:length(ne.markers)){
      print(
        plotReducedDim(sce.test, dimred="UMAP",
                       colour_by=ne.markers[i], by_exprs_values="logcounts",
                       point_alpha=0.2, point_size=1.2, theme_size=7) +
          ggtitle(paste0("Clustering on fastMNN PCs \nLog counts for ", ne.markers[i],": ",names(ne.markers)[i]))
            )
      }
    dev.off()
    
    ## doubletScore distributions?
    cellClust.idx <- splitit(sce.test$clusters.glmcpcamnn)
    sapply(cellClust.idx, function(x){quantile(sce.test$doubletScore[x])})
    
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
    
    markers.neurons <- c("SYT1", "SLC17A7", "SLC17A6","GAD2")
    markers.neurons <- rowData(sce.test)$gene_id[match(markers.neurons, rowData(sce.test)$gene_name)]
    names(markers.neurons) <- c("SYT1", "SLC17A7", "SLC17A6","GAD2")
    
    sapply(medianNon0.glmpcamnn, function(x){
      x[markers.neurons] == TRUE
    })
        # ~half (+) for pan-neuronal marker SYT1 (this usually is captured better than SNAP25)
    neuron.clusters <- which(sapply(medianNon0.glmpcamnn, function(x){ x[markers.neurons["SYT1"]] == TRUE }))
        # So actually 41/60 clusters
    
    table(sce.test$clusters.glmcpcamnn %in% neuron.clusters)["TRUE"] / ncol(sce.test) # 0.8554533
         # Looks like ~86% of these nuclei are putatively neuronal
    
    sce.test$SYT1 <- ifelse(sce.test$clusters.glmcpcamnn %in% neuron.clusters, "neuronal", "nonNeuronal")
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
        height=6, width=12)
    plotExpressionCustom(sce = sce.neuron.mnn,
                         exprs_values = "logcounts",
                         features = names(markers.neurons),
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










## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
#[1] "2022-01-05 21:59:27 EST"
proc.time()
#     user    system   elapsed 
# 2600.322   571.665 15194.725 
options(width = 120)
session_info()
# ─ Session info ─────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-01-05
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────
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
# cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.2)
# colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
# cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
# crayon                 1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# DropletUtils         * 1.14.1   2021-11-08 [2] Bioconductor
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.0   2021-10-26 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                   1.6.0    2021-12-17 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# gridExtra            * 2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# igraph                 1.2.11   2022-01-04 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.0   2021-10-26 [2] Bioconductor
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# metapod                1.2.0    2021-10-26 [2] Bioconductor
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.6.4    2021-10-18 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
# RCurl                  1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
# ResidualMatrix         1.4.0    2021-10-26 [1] Bioconductor
# rhdf5                  2.38.0   2021-10-26 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rlang                  0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
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
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.4.3    2021-11-30 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/ntranngu/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────

