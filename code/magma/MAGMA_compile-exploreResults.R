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























### Explore gene-level MAGMA (step 2) stats ==========================
  #   (particularly for those with expected but not observed GSA hits)

load(here("processed_data","SCE","sce_updated_LC.rda"), verbose=T)
    # sce.lc, annotationTab.lc, medianNon0.lc, hdgs.lc, cell_colors.lc

adhd.gene.out <- read.table(here("code","magma","SNP_Data","ADHD_Demontis2019_LC_snp-wise.genes.out"),
                            sep="\t", header=T)
    # Not a tsv...
    # Just grep those genes from the 12 risk loci (see Supp Table 7)
head(adhd.gene.out)
    #  GENE.............CHR......START.......STOP...NSNPS..NPARAM......N........ZSTAT............P
    # 1 ENSG00000237491    1     679127     755445      14       1  18603       1.3908     0.082142
    # 2 ENSG00000177757    1     717751     765217      33       1  18428       1.2629      0.10331
    # 3 ENSG00000228794    1     725518     813582      62       2  17884      0.98915       0.1613


adhd.loci <- c("ST3GAL3", "PTPRF", "SPAG16", "PCDH7", "LINC00461", "MEF2C", "FOXP2",
               "LINC01288", "SORCS3", "DUSP6", "SEMA6D", "LINC01572")
adhd.loci %in% rowData(sce.lc)$gene_name  # all TRUE
names(adhd.loci) <- rowData(sce.lc)$gene_id[match(adhd.loci, rowData(sce.lc)$gene_name)]

sapply(names(adhd.loci), function(x){
  adhd.gene.out[grep(x, adhd.gene.out[ ,1]), ]
})
    # "ENSG00000126091    1   44136495   44406837     539      20  22519       6.7357   8.1549e-12" 
    # "ENSG00000142949    1   43955858   44099337     309      16  22610        6.306    1.432e-10" 
    # "ENSG00000144451    2  214114103  215285225    3094      67  22525       2.1837     0.014491" 
    # "ENSG00000169851    4   30687037   31158427    1203      55  22452       3.4483   0.00028209" 
    # "ENSG00000245526    5   87793363   88021874     407      35  22422       4.8311   6.7885e-07" 
    # "ENSG00000081189    5   88002934   88235074     386      27  22445       5.0497   2.2126e-07" 
    # "ENSG00000128573    7  113691382  114343827     890      44  22312       4.4626   4.0481e-06" 
    # "ENSG00000254344    8   34606546   34732882     199      10  21790      0.34547      0.36487" 
    # "ENSG00000156395   10  106366048  107035000    1915      51  22626       5.6459   8.2183e-09" 
    # "ENSG00000139318   12   89731012   89781278     135      16  22100       6.0285   8.2769e-10" 
    # "ENSG00000137872   15   47441298   48076420    1609      59  22503       5.5963   1.0949e-08" 
    # "ENSG00000261008   16   72260180   72733913     863      25  22519       2.5675    0.0051221" 

    #       - so the z-score recapitulates significance in all but one: ENSG00000254344 (LINC01288)
    #           - discrepancy b/tw available summary stats vs. full data results?








## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
# 
proc.time()
#    user   system  elapsed 
# 
options(width = 120)
session_info()

