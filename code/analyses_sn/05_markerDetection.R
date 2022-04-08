### LC snRNA-seq analysis
### Marker detection with t-test
###     qsub -l bluejay,mf=76G,h_vmem=80G
### Initiated: MNT 07Apr2022

library(SingleCellExperiment)
library(scater)
library(scran)
library(scry)
library(batchelor)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)


### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

here()
# [1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"

source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")

## ===

load(here("processed_data","SCE", "sce_updated_LC.rda"), verbose=T)

table(sce.lc$cellType.merged)
    #

## doubletScore & sum distributions / cluster?
cellClust.idx <- splitit(sce.lc$cellType.merged)
sapply(cellClust.idx, function(x){round(quantile(sce.lc$doubletScore[x]), 2)})


sapply(cellClust.idx, function(x){quantile(sce.lc$sum[x])})
    #



# Remove 0 genes across all nuclei
sce.lc <- sce.lc[!rowSums(assay(sce.lc, "counts"))==0, ]  #


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.lc

assay(sce.lc, "logcounts") <- NULL
sizeFactors(sce.lc) <- NULL
sce.lc <- logNormCounts(sce.lc)


### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellClust.idx <- splitit(sce.lc$cellType.merged)
# medianNon0.lc.26 <- lapply(cellClust.idx, function(x){
#   apply(as.matrix(assay(sce.lc, "logcounts")), 1, function(y){
#     median(y[x]) > 0
#   })
# })
# 
# sapply(medianNon0.lc.26, table)
#     #
# 
# 
# # # Just for interactive exploration of some of these
# # plotExpressionCustom(sce = sce.lc,
# #                      exprs_values = "logcounts",
# #                      # 
# #                      features = head(Oligo_A_markers,4),
# #                      features_name = "custom-selected",
# #                      anno_name = "cellType.merged",
# #                      ncol=2, point_alpha=0.4, point_size=0.9,
# #                      scales="free_y", swap_rownames="gene_name") +
# #   ggtitle(label=paste0("Lionel's custom markers of interest")) +
# #   theme(plot.title = element_text(size = 12),
# #         axis.text.x = element_text(size=7))
# 
# 
# ## Traditional t-test, pairwise ===
# mod <- with(colData(sce.lc), model.matrix(~ Sample))
# mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`
# 
# # Run pairwise t-tests
# markers.lc.t.pw <- findMarkers(sce.lc, groups=sce.lc$cellType.merged,
#                                assay.type="logcounts", design=mod, test="t",
#                                direction="up", pval.type="all", full.stats=T)
# 
# sapply(markers.lc.t.pw, function(x){table(x$FDR<0.05)})
#     #
# 
# 
# 
# 
# # Add respective 'non0median' column to the stats for each set of markers
# for(i in names(markers.lc.t.pw)){
#   markers.lc.t.pw[[i]] <- cbind(markers.lc.t.pw[[i]],
#                                 medianNon0.lc.26[[i]][match(rownames(markers.lc.t.pw[[i]]),
#                                                          names(medianNon0.lc.26[[i]]))])
#   colnames(markers.lc.t.pw[[i]])[29] <- "non0median"
# }
# 
# sapply(markers.lc.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
#     #
# 
# 
# ## Save these
# save(markers.lc.t.pw, medianNon0.lc.26,
#      file=here("processed_data","SCE",
#                "markers-stats_LC-n3_findMarkers_26cellTypes.rda"))


    # As needed
    load(here("processed_data","SCE",
              "markers-stats_LC-n3_findMarkers_26cellTypes.rda"), verbose=T)


# Print these to pngs
markerList.t.pw <- lapply(markers.lc.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
 }
)

# Change to gene symbols
markerList.t.pw <- lapply(markerList.t.pw, function(x){
  rowData(sce.lc)$gene_name[match(x, rowData(sce.lc)$gene_id)]
  })

genes.top40.t <- lapply(markerList.t.pw, function(x){head(x, n=40)})

smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20]
left.set <- setdiff(names(genes.top40.t), smaller.set)
# Now remove those with no significant markers
smaller.set <- setdiff(smaller.set,
                       names(genes.top40.t)[lengths(genes.top40.t) == 0])

# Smaller graphical window
#dir.create(here("plots","snRNA-seq","markers"))
for(i in smaller.set){
  png(here("plots","snRNA-seq","markers",
           paste0("LC_t_pairwise_topMarkers-", i, "_vlnPlots.png")), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         exprs_values = "logcounts",
                         features = head(genes.top40.t[[i]],n=4), 
                         features_name = i,
                         anno_name = "cellType.merged",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values = c(tableau20, tableau10medium)) +  
      ggtitle(label=paste0("LC ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}

# 20-40 markers
for(i in left.set){
  png(here("plots","snRNA-seq","markers",
           paste0("LC_t_pairwise_topMarkers-", i, "_vlnPlots.png")), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         exprs_values = "logcounts",
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType.merged",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values = c(tableau20, tableau10medium)) +  
      ggtitle(label=paste0("LC ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



### Cluster-vs-all-others single-nucleus-level iteration ========

# (Pre-process as above, if needed)
mod <- with(colData(sce.lc), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.lc.t.1vAll <- list()
for(i in levels(sce.lc$cellType.merged)){
  # Make temporary contrast
  sce.lc$contrast <- ifelse(sce.lc$cellType.merged==i, 1, 0)
  # Test cluster vs. all others
  markers.lc.t.1vAll[[i]] <- findMarkers(sce.lc, groups=sce.lc$contrast,
                                         assay.type="logcounts", design=mod, test="t",
                                         std.lfc=TRUE,
                                         direction="up", pval.type="all", full.stats=T)
}


class(markers.lc.t.1vAll[[1]])
    # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
    # -> we want the second entry, named "1"
    #    (for other purposes, might be interesting to look into that "0" entry, which
    #     is basically what genes are depleted in the cell type of interest)


## Do some reorganizing
markers.lc.t.1vAll <- lapply(markers.lc.t.1vAll, function(x){
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y){ y[ ,4] }) 
})

# Re-name std.lfc column and the entries; add non-0-median info
for(i in names(markers.lc.t.1vAll)){
  colnames(markers.lc.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.lc.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.lc.t.1vAll[[i]][["0"]] <- cbind(markers.lc.t.1vAll[[i]][["0"]],
                                          medianNon0.lc.26[[i]][match(rownames(markers.lc.t.1vAll[[i]][["0"]]),
                                                                   names(medianNon0.lc.26[[i]]))])
  colnames(markers.lc.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.lc.t.1vAll[[i]][["1"]] <- cbind(markers.lc.t.1vAll[[i]][["1"]],
                                          medianNon0.lc.26[[i]][match(rownames(markers.lc.t.1vAll[[i]][["1"]]),
                                                                   names(medianNon0.lc.26[[i]]))])
  colnames(markers.lc.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.lc.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}


## Save these
save(markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.26,
     file=here("processed_data","SCE",
               "markers-stats_LC-n3_findMarkers_26cellTypes.rda"))



## Marker numbers with the non-0-median filter
sapply(markers.lc.t.1vAll, function(x){
  table(x[[2]]$log.FDR < log(.001) & x[[2]]$non0median == TRUE)
})
#


## Print these to pngs
markerList.t.1vAll <- lapply(markers.lc.t.1vAll, function(x){
  rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
 }
)

# Change to gene symbols
markerList.t.1vAll <- lapply(markerList.t.1vAll, function(x){
  rowData(sce.lc)$gene_name[match(x, rowData(sce.lc)$gene_id)]
})

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(here("plots","snRNA-seq","markers",
           paste0("LC_t_1vALL_topMarkers-",i,"_vlnPlots.png")), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         exprs_values = "logcounts",
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType.merged",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values = c(tableau20, tableau10medium)) +  
      ggtitle(label=paste0("LC ", i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



## Write these top 40 lists to a csv ===
names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# Many of the PW results don't have 40 markers:
extend.idx <- names(which(lengths(markerList.t.pw) < 40))
for(i in extend.idx){
  markerList.t.pw[[i]] <- c(markerList.t.pw[[i]], rep("", 40-length(markerList.t.pw[[i]])))
}
# 1vAllOthers
# Many of the PW results don't have 40 markers:
extend.idx <- names(which(lengths(markerList.t.1vAll) < 40))
for(i in extend.idx){
  markerList.t.1vAll[[i]] <- c(markerList.t.1vAll[[i]], rep("", 40-length(markerList.t.1vAll[[i]])))
}

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]
write.csv(top40genes, file=here("code","analyses_sn","top40genesLists_LC-n3_26cellTypes.csv"),
          row.names=FALSE)



## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
# 
proc.time()
#     user    system   elapsed 
#
options(width = 120)
session_info()



