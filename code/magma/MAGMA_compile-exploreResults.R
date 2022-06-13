### MAGMA with locus coeruleus (LC) snRNA-seq populations
  #   - Compile gene set analyses (step3) results
  #   - Explore gene-level stats for provided summary stats
# MNT 09Jun2022

library(rtracklayer)
library(GenomicRanges)
library(jaffelab)
library(SingleCellExperiment)
library(sessioninfo)
library(here)

here()
#[1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"

# ==



### Compile GSA (step 3) results ====================================
  #   (into heatmap similar to in Tran, Maynard, et al. Neuron 2021)

magmaStats <- list()

magmaStats[["LC"]][["PD.meta2019"]] <- read.table(here("code","magma","Results","lc_PD.gsa.out"), header=T)
magmaStats[["LC"]][["ADHD.PGC.iPsych"]] <- read.table(here("code","magma","Results","lc_ADHD.gsa.out"), header=T)
magmaStats[["LC"]][["AD.metaPh3"]] <- read.table(here("code","magma","Results","lc_Alzheimers.gsa.out"), header=T)


## Merge to assess significance thresholds ===
magmaStats_list = lapply(magmaStats, function(m) {
  z = sapply(m, "[[", "P")
  rownames(z) = m[[1]]$VARIABLE
  z
})

magmaStats_wide = as.data.frame(do.call("rbind", magmaStats_list))
magmaStats_wide$Region = rep("LC", 
                             times = sapply(magmaStats_list, nrow))
magmaStats_wide$CellType = rownames(magmaStats_wide)
## reshape to long
magmaStats_long = reshape2::melt(magmaStats_wide)
colnames(magmaStats_long)[3:4] = c("GWAS", "P")
dim(magmaStats_long)
# [1] 57    4     - 57 bc 19 reported cell classes x 3 GWAS'

table(p.adjust(magmaStats_long$P, "fdr") < 0.1)
    #FALSE  TRUE 
    #   54     3    (none at < 0.05)
betacut.fdr <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "fdr") < 0.1])
    # [1] 0.003785  (allowing for FDR < 0.1)

magmaStats_long$P.adj.fdr <- p.adjust(magmaStats_long$P, "fdr")
magmaStats_long[which(magmaStats_long$P.adj.fdr < 0.1), ]
    #    Region CellType       GWAS          P P.adj.fdr
    # 39     LC    Astro AD.metaPh3 0.00178800  0.050958
    # 53     LC    Micro AD.metaPh3 0.00099219  0.050958
    # 57     LC    Oligo AD.metaPh3 0.00378500  0.071915


# Bonferroni significance?
table(p.adjust(magmaStats_long$P, "bonf") < 0.05) # none, but otherwise Microglia for AD at < 0.1
#betacut.bonf <- max(magmaStats_long$P[p.adjust(magmaStats_long$P, "bonf") < 0.05])



# For betas ===
magmaStats_list.beta = lapply(magmaStats, function(m) {
  z = sapply(m, "[[", "BETA")
  rownames(z) = m[[1]]$VARIABLE
  z
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
          row.names=F, quote=F)





















### Explore gene-level MAGMA (step 2) stats ==========================
  #   (particularly for those with expected but not observed GSA hits)

load(here("processed_data","SCE","sce_updated_LC.rda"), verbose=T)
    # sce.lc, annotationTab.lc, medianNon0.lc, hdgs.lc, cell_colors.lc

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





## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
# 
proc.time()
#    user   system  elapsed 
# 
options(width = 120)
session_info()

