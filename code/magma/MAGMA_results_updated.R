### MAGMA results
### code by Matthew N Tran
### last updated by Lukas M Weber, Oct 2022

### MAGMA with locus coeruleus (LC) snRNA-seq NE neuron cluster
  #   - Compile gene set analyses (step3) results
  #   - Explore gene-level stats for provided summary stats


library(rtracklayer)
library(GenomicRanges)
library(jaffelab)
library(SingleCellExperiment)
library(fields)
library(RColorBrewer)
library(grid)
library(here)


### Compile GSA (step 3) results ====================================
  #   (into heatmap similar to in Tran, Maynard, et al. Neuron 2021)

magmaStats <- list()

magmaStats[["LC"]][["PD.Nalls.2019"]] <- read.table(here("code","magma","Results","lc_PD.gsa.out"), header = TRUE)
magmaStats[["LC"]][["ADHD.Demontis.2018"]] <- read.table(here("code","magma","Results","lc_ADHD.gsa.out"), header = TRUE)
magmaStats[["LC"]][["AD.Jansen.2019"]] <- read.table(here("code","magma","Results","lc_Alzheimers.gsa.out"), header = TRUE)


## Merge to assess significance thresholds ===

# magmaStats_list = lapply(magmaStats, function(m) {
#   z = sapply(m, "[[", "P")
#   rownames(z) = m[[1]]$VARIABLE
#   z
# })

magmaStats_list = lapply(magmaStats, function(m) {
  z = sapply(m, "[[", "P")
  mat = matrix(z, nrow = 1)
  colnames(mat) = names(z)
  rownames(mat) = m[[1]]$VARIABLE
  mat
})

magmaStats_wide = as.data.frame(do.call("rbind", magmaStats_list))
magmaStats_wide$Region = rep("LC", 
                             times = sapply(magmaStats_list, nrow))
magmaStats_wide$CellType = rownames(magmaStats_wide)
## reshape to long
magmaStats_long = reshape2::melt(magmaStats_wide)
colnames(magmaStats_long)[3:4] = c("GWAS", "P")
dim(magmaStats_long)
    # note: 1 cell class x 3 GWAS'

table(p.adjust(magmaStats_long$P, "fdr") < 0.1)

betacut.fdr <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "fdr") < 0.1])

magmaStats_long$P.adj.fdr <- p.adjust(magmaStats_long$P, "fdr")
magmaStats_long[which(magmaStats_long$P.adj.fdr < 0.1), ]

magmaStats_long


# alternatively: Bonferroni significance
table(p.adjust(magmaStats_long$P, "bonf") < 0.05)


# For betas ===

# magmaStats_list.beta = lapply(magmaStats, function(m) {
#   z = sapply(m, "[[", "BETA")
#   rownames(z) = m[[1]]$VARIABLE
#   z
# })

magmaStats_list.beta = lapply(magmaStats, function(m) {
  z = sapply(m, "[[", "BETA")
  mat = matrix(z, nrow = 1)
  colnames(mat) = names(z)
  rownames(mat) = m[[1]]$VARIABLE
  mat
})

magmaStats_wide.beta = as.data.frame(do.call("rbind", magmaStats_list.beta))
magmaStats_wide.beta$Region = rep("LC", 
                                  times = sapply(magmaStats_list.beta, nrow))
magmaStats_wide.beta$CellType = rownames(magmaStats_wide.beta)
## reshape to long
magmaStats_long.beta = reshape2::melt(magmaStats_wide.beta)
colnames(magmaStats_long.beta)[3:4] = c("GWAS", "Beta")

# Check before appending
table(paste0(magmaStats_long$CellType,":",magmaStats_long$GWAS) ==
        paste0(magmaStats_long.beta$CellType,":", magmaStats_long.beta$GWAS))

# cbind Beta
magmaStats_long$Beta <- magmaStats_long.beta$Beta
# Reorder
magmaStats_long <- magmaStats_long[ ,c("Region", "CellType", "GWAS", "Beta", "P", "P.adj.fdr")]


## Print to CSV ===
write.csv(magmaStats_long, file = here("code","magma","Results",
                                       "table_magma-GSA_19mergedClusters_3xGWAS.csv"),
          row.names = FALSE, quote = FALSE)



## Make heatmap ===
midpoint = function(x) x[-length(x)] + diff(x)/2

MAGMAplot = function(region, Pthresh, fdrThresh, ...) {
  ## Set up -log10(p's)
  wide_p = sapply(magmaStats[[region]], function(x){cbind(-log10(x$P))})
  rownames(wide_p) <- magmaStats[[region]][[1]]$VARIABLE
  wide_p[wide_p > Pthresh] = Pthresh
  wide_p <- round(wide_p[rev(sort(rownames(wide_p))), ], 3)
  
  
  ## Set up betas
  wide_beta <- sapply(magmaStats[[region]], function(x){cbind(x$BETA)})
  rownames(wide_beta) <- magmaStats[[region]][[1]]$VARIABLE
  wide_beta <- round(wide_beta[rev(sort(rownames(wide_beta))), ], 2)
  
  # # Use empirical cutoff (independentp=0.05) for printing betas 
  wide_beta[wide_p < -log10(0.05)] = ""
  # and Bonf. cutoff for bolding
  customFont <- ifelse(wide_p < -log10(fdrThresh), 1, 2)
  customCex <- ifelse(wide_p < -log10(fdrThresh), 0.9, 1.0)
  
  
  ## Plot
  clusterHeights <- seq(0,160,length.out=nrow(wide_p)+1)
  mypal = c("white", colorRampPalette(brewer.pal(9,"YlOrRd"))(60))[1:30]
  xlabs <- colnames(wide_p)
  
  # Heatmap of p's
  image.plot(x = seq(0,ncol(wide_p),by=1), y = clusterHeights, z = as.matrix(t(wide_p)),
             col = mypal,xaxt="n", yaxt="n",xlab = "", ylab="", ...)
  axis(2, rownames(wide_p), at=midpoint(clusterHeights), las=1)
  axis(1, rep("", ncol(wide_p)), at = seq(0.5,ncol(wide_p)-0.5))
  text(x = seq(0.5,ncol(wide_p)-0.5), y=-1*max(nchar(xlabs))/2, xlabs,
       xpd=TRUE, srt=45, cex=1.2, adj= 1)
  abline(h=clusterHeights,v=0:ncol(wide_p))
  
  # Print top decile of betas
  text(x = rep(seq(0.5,ncol(wide_p)-0.5),each = nrow(wide_p)), 
       y = rep(midpoint(clusterHeights), ncol(wide_p)),
       as.character(wide_beta),
       # If Bonf, a little bigger
       cex=customCex,
       # If Bonf, 2 (bold)
       font=customFont)
}

# Plot
pdf(here("plots","snRNA-seq","heatmap_LC_MAGMA-GSAresults_AD-ADHD-PD_GWAS.pdf"), w=6)
par(mar=c(8.5,7.5,6,1), cex.axis=1.0, cex.lab=0.5)
MAGMAplot(region="LC", Pthresh=12, fdrThresh=betacut.fdr)
abline(v=7,lwd=3)
text(x = 1.5, y=185, "MAGMA gene set analyses: LC cell classes", xpd=TRUE, cex=1.5, font=2)
text(x = 1.5, y=175, "(Betas for empirically significant associations, p < 0.05)", xpd=TRUE, cex=1, font=1)
grid::grid.text(label="-log10(p-value)", x=0.92, y=0.825, gp=gpar(fontsize=9))
dev.off()




### Explore gene-level MAGMA (step 2) stats ==========================
  #   (particularly for those with expected but not observed GSA hits)

load(here("processed_data","SCE","sce_updated_LC.rda"), verbose=T)
    # sce.lc, annotationTab.lc, medianNon0.lc, hdgs.lc, cell_colors.lc

# Also load marker stats (are any high risk-scoring genes are cell type markers?)
load(here("processed_data","SCE","markers-stats_LC-n3_findMarkers_19cellTypes.rda"), verbose=T)
    # markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.19


adhd.gene.out <- read.table(here("code","magma","SNP_Data","ADHD_Demontis2019_LC_snp-wise.genes.out"),
                            sep="", header=T)
    # Not a tsv, so use the defeault `sep=""` (parses thru spaces and creates a data.frame)

head(adhd.gene.out)
    #              GENE CHR  START   STOP NSNPS NPARAM     N    ZSTAT        P
    # 1 ENSG00000237491   1 679127 755445    14      1 18603  1.39080 0.082142
    # 2 ENSG00000177757   1 717751 765217    33      1 18428  1.26290 0.103310
    # 3 ENSG00000228794   1 725518 813582    62      2 17884  0.98915 0.161300


adhd.loci <- c("ST3GAL3", "PTPRF", "SPAG16", "PCDH7", "LINC00461", "MEF2C", "FOXP2",
               "LINC01288", "SORCS3", "DUSP6", "SEMA6D", "LINC01572")
adhd.loci %in% rowData(sce.lc)$gene_name  # all TRUE
names(adhd.loci) <- rowData(sce.lc)$gene_id[match(adhd.loci, rowData(sce.lc)$gene_name)]

rownames(adhd.gene.out) <- adhd.gene.out$GENE

adhd.gene.out[names(adhd.loci), ]
    #                            GENE CHR     START      STOP NSNPS NPARAM     N
    # ENSG00000126091 ENSG00000126091   1  44136495  44406837   539     20 22519
    # ENSG00000142949 ENSG00000142949   1  43955858  44099337   309     16 22610
    # ENSG00000144451 ENSG00000144451   2 214114103 215285225  3094     67 22525
    # ENSG00000169851 ENSG00000169851   4  30687037  31158427  1203     55 22452
    # ENSG00000245526 ENSG00000245526   5  87793363  88021874   407     35 22422
    # ENSG00000081189 ENSG00000081189   5  88002934  88235074   386     27 22445
    # ENSG00000128573 ENSG00000128573   7 113691382 114343827   890     44 22312
    # ENSG00000254344 ENSG00000254344   8  34606546  34732882   199     10 21790
    # ENSG00000156395 ENSG00000156395  10 106366048 107035000  1915     51 22626
    # ENSG00000139318 ENSG00000139318  12  89731012  89781278   135     16 22100
    # ENSG00000137872 ENSG00000137872  15  47441298  48076420  1609     59 22503
    # ENSG00000261008 ENSG00000261008  16  72260180  72733913   863     25 22519
    #                   ZSTAT          P
    # ENSG00000126091 6.73570 8.1549e-12
    # ENSG00000142949 6.30600 1.4320e-10
    # ENSG00000144451 2.18370 1.4491e-02
    # ENSG00000169851 3.44830 2.8209e-04
    # ENSG00000245526 4.83110 6.7885e-07
    # ENSG00000081189 5.04970 2.2126e-07
    # ENSG00000128573 4.46260 4.0481e-06
    # ENSG00000254344 0.34547 3.6487e-01
    # ENSG00000156395 5.64590 8.2183e-09
    # ENSG00000139318 6.02850 8.2769e-10
    # ENSG00000137872 5.59630 1.0949e-08
    # ENSG00000261008 2.56750 5.1221e-03

    #       - so the z-score recapitulates 'significance' in all but one: ENSG00000254344 (LINC01288)
    #           - discrepancy b/tw available summary stats vs. full data results?

quantile(adhd.gene.out$ZSTAT, probs=seq(0.9,1,by=0.01))
round(quantile(adhd.gene.out$ZSTAT, probs=seq(0.9,1,by=0.01)),3)
    #   90%   91%   92%   93%   94%   95%   96%   97%   98%   99%  100% 
    # 1.760 1.832 1.909 2.001 2.099 2.216 2.352 2.514 2.783 3.191 7.101 
        # So there's 1% that have a gene-level Z-score higher than 3.19

round(quantile(adhd.gene.out$ZSTAT, probs=seq(0.95,1,by=0.0025)),3)
    #   95% 95.25%  95.5% 95.75%    96% 96.25%  96.5% 96.75%    97% 97.25%  97.5% 
    # 2.216  2.245  2.280  2.313  2.352  2.383  2.428  2.469  2.514  2.559  2.624 

    #97.75%    98% 98.25%  98.5% 98.75%    99% 99.25%  99.5% 99.75%   100% 
    # 2.703  2.783  2.861  2.960  3.074  3.191  3.331  3.622  3.949  7.101 


## Print these in the 19 populations ===
sce.lc <- sce.lc[ ,-grep("ambig.lowNTx_", sce.lc$cellType.merged)]
sce.lc$cellType.merged <- droplevels(sce.lc$cellType.merged)

colors2print <- cell_colors.lc[-grep("ambig.lowNTx_", names(cell_colors.lc))]

pdf(here("plots","snRNA-seq","exploration","MAGMA_ADHD-reported-loci_vlnPlots_19clusters.pdf"), width=9)
plotExpressionCustom(sce = sce.lc,
                     exprs_values = "logcounts",
                     features = adhd.loci, 
                     features_name = "",
                     anno_name = "cellType.merged",
                     ncol=4, point_alpha=0.4,
                     scales="free_y", swap_rownames="gene_name") +
  scale_color_manual(values = colors2print) +  
  ggtitle(label=paste0("ADHD 12 reported loci in LC (n=3) snRNA-seq clusters")) +
  theme(plot.title = element_text(size = 12))
dev.off()



# What about intersection of cell class markers x MAGMA-computed high Z-scores?
markerList.t.pw <- lapply(markers.lc.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
}
)

# Try at some Z cutoffs, based on above distributions
table(adhd.gene.out$ZSTAT >= 4) # 70
table(adhd.gene.out$ZSTAT >= 3.6) # 155
topGenes.adhd <- adhd.gene.out$GENE[adhd.gene.out$ZSTAT > 3.6]

topGenes.markers.adhd <- intersect(unlist(markerList.t.pw), topGenes.adhd)

names(topGenes.markers.adhd) <- rowData(sce.lc)$gene_name[match(topGenes.markers.adhd,
                                                                rowData(sce.lc)$gene_id)]
    #[1] "RPLP2"     "CUBN"      "LINC02163" "CEND1"     "COL19A1"   "KIZ"      
    #[7] "MEF2C"     "LINC01524"
        # At Z-stats >= 3.6:
        # [1] "SNRK"      "RPLP2"     "XRN2"      "CUBN"      "LINC02163" "KCNIP1"   
        # [7] "CEND1"     "COL19A1"   "KIZ"       "GPM6A"     "MEF2C"     "LINC01524"
        # [13] "RGS9"  

pdf(here("plots","snRNA-seq","exploration","MAGMA_ADHD_highZs-vs-markers_vlnPlots_19clusters.pdf"), width=9)
plotExpressionCustom(sce = sce.lc,
                     exprs_values = "logcounts",
                     features = names(topGenes.markers.adhd), 
                     features_name = "",
                     anno_name = "cellType.merged",
                     ncol=4, point_alpha=0.4,
                     scales="free_y", swap_rownames="gene_name") +
  scale_color_manual(values = colors2print) +  
  ggtitle(label=paste0("ADHD: high MAGMA-computed Z's (>= 3.6) in LC (n=3) snRNA-seq clusters")) +
  theme(plot.title = element_text(size = 12))
dev.off()




## Parkinson's ===
pd.gene.out <- read.table(here("code","magma","SNP_Data","PD_Nalls2019_LC_snp-wise.genes.out"),
                            sep="", header=T)

head(pd.gene.out)

round(quantile(pd.gene.out$ZSTAT, probs=seq(0.9,1,by=0.01)),3)
round(quantile(pd.gene.out$ZSTAT, probs=seq(0.95,1,by=0.0025)),3)
    #   95% 95.25%  95.5% 95.75%    96% 96.25%  96.5% 96.75%    97% 97.25%  97.5% 
    # 1.868  1.898  1.934  1.974  2.016  2.063  2.121  2.169  2.219  2.274  2.337 

    #97.75%    98% 98.25%  98.5% 98.75%    99% 99.25%  99.5% 99.75%   100% 
    # 2.399  2.473  2.566  2.684  2.867  3.048  3.342  3.734  5.088  7.916 

table(pd.gene.out$ZSTAT >= 5.0) # 83
    # 83 (90 reported index SNPs)
topGenes.pd <- pd.gene.out$GENE[pd.gene.out$ZSTAT >= 3.7]
table(topGenes.pd %in% rowData(sce.lc)$gene_id) # all there
topGenes.pd <- rowData(sce.lc)$gene_name[match(topGenes.pd, rowData(sce.lc)$gene_id)]
topGenes.pd

# Take those from Figure 2 with significant variants at p < 5e-9 (red; olive label)
fig2.genes <- c("PMVK", "ITPKB", "SIPA1L2", "GBA", "NUCKS1", "KRTCAP2", "STK39", "RAB29", "TMEM163",
                "SATB1", "IP6K2", "GAK", "BST1", "MCCC1", "TMEM175", "SNCA", "HLA-DRB5", "ELOVL7",
                "FAM47E", "STBD1", "CAMK2D", "SCARB2", "GPNMB", "LOC100131289", "CTSB", "BIN3", "FGF20",
                "SH3GL2", "ITGA8", "IGSF9B", "BAG2", "DLG2", "INPP5F", "LRRK2", "GCH1", "CHD9",
                "GALC", "SYT17", "CASC16", "VPS13C", "SETD1A", "RETREG3", "WNT3", "CRHR1", "RIT2", "SPPL2B"
                )
    # 46
table(fig2.genes %in% rowData(sce.lc)$gene_name)  # 45/46; "LOC100131289" is missing

    # Btw--  SH3GL2, CRHR1, labeled twice; so is FGF20, but on different chromosomes...


intersect(fig2.genes, topGenes.pd)
    # [1] "NUCKS1"  "STK39"   "RAB29"   "TMEM163" "GAK"     "BST1"    "MCCC1"  
    # [8] "TMEM175" "SNCA"    "ELOVL7"  "FAM47E"  "STBD1"   "SCARB2"  "GPNMB"  
    # [15] "SH3GL2"  "IGSF9B"  "LRRK2"   "CASC16"  "SETD1A"  "WNT3"    "CRHR1"  
    # [22] "RIT2" 
    # At just Z >= 5:
        # [1] "NUCKS1"  "RAB29"   "TMEM163" "GAK"     "BST1"    "TMEM175" "SNCA"   
        # [8] "ELOVL7"  "FAM47E"  "STBD1"   "SCARB2"  "IGSF9B"  "LRRK2"   "SETD1A" 
        # [15] "WNT3"    "CRHR1" 

# Where is GCH1?  (rate-limiting enzyme for BH4 production, an essential
#                  cofactor for tyrosine & tryptophan hydroxylases)
#     - It's like in the top ~1.5% of z-scores, along with TH
pd.gene.out$gene.symbol <- rowData(sce.lc)$gene_name[match(pd.gene.out$GENE, rowData(sce.lc)$gene_id)]
pd.gene.out[pd.gene.out$gene.symbol %in% c("GCH1", "TH"), ]
    #                  GENE CHR    START     STOP NSNPS NPARAM     N  ZSTAT         P gene.symbol
    # 16778 ENSG00000180176  11  2175159  2228107   372     48 65561 2.7565 0.0029210          TH
    # 21328 ENSG00000131979  14 55298726 55404544   517     41 65179 2.7693 0.0028086        GCH1

# Distribution of z-scores for those reported loci
round(quantile(pd.gene.out[pd.gene.out$gene.symbol %in% fig2.genes, ]$ZSTAT,
               probs=seq(0.1,1,by=0.1)), 3)
    #   10%   20%   30%   40%   50%   60%   70%   80%   90%  100% 
    # 0.291 0.699 2.037 2.693 3.579 4.859 5.322 5.892 7.003 7.916 

    ## Keep in mind that they tagged genes within -/+ 1MB of index SNPs, whereas these
     #    ZSTATs were computed by -35kb, +10kb of the genes


## Print the intersection
print.intersecting <- intersect(fig2.genes, topGenes.pd)

pdf(here("plots","snRNA-seq","exploration","MAGMA_Zabove3.7_and_PD-reported-loci_vlnPlots_19clusters.pdf"), width=9, height=8)
plotExpressionCustom(sce = sce.lc,
                     exprs_values = "logcounts",
                     features = print.intersecting, 
                     features_name = "",
                     anno_name = "cellType.merged",
                     ncol=4, point_alpha=0.4,
                     scales="free_y", swap_rownames="gene_name") +
  scale_color_manual(values = colors2print) +  
  ggtitle(label=paste0("PD risk loci (p<5x10-9) vs. MAGMA (Z>=3.7) in LC (n=3) snRNA-seq clusters")) +
  theme(plot.title = element_text(size = 12))
dev.off()


## All of those genes in the reported set that are PW markers
markerList.t.pw <- lapply(markerList.t.pw, function(x){
  rowData(sce.lc)$gene_name[match(x, rowData(sce.lc)$gene_id)]
})

print.these <- intersect(unlist(markerList.t.pw), fig2.genes)

pdf(here("plots","snRNA-seq","exploration","MAGMA_PD-reported-loci_vs-markers_vlnPlots_19clusters.pdf"), width=8)
plotExpressionCustom(sce = sce.lc,
                     exprs_values = "logcounts",
                     features = c("MAP4K4", "ABHD6"), 
                     features_name = "",
                     anno_name = "cellType.merged",
                     ncol=3, point_alpha=0.4,
                     scales="free_y", swap_rownames="gene_name") +
  scale_color_manual(values = colors2print) +  
  ggtitle(label=paste0("PD risk loci (p<5x10-9) that are top PW markers in LC (n=3) snRNA-seq clusters")) +
  theme(plot.title = element_text(size = 12))
dev.off()


## And then ignoring the reported set; highest Z's (>=3.7) and also PW markers
setdiff(intersect(topGenes.pd, unlist(markerList.t.pw)), print.intersecting)
    # Just two (at this Z threshold): are most expressed in oligos, but not uniquely
    # [1] "MAP4K4" "ABHD6" 



## Alzheimer's ===
ad.gene.out <- read.table(here("code","magma","SNP_Data","AD_Jansen2019_LC_snp-wise.genes.out"),
                          sep="", header=T)

head(ad.gene.out)


