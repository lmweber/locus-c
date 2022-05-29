### for MAGMA with locus coeruleus snRNA-seq populations
  # MNT 27May2022

library(rtracklayer)
library(GenomicRanges)
#BiocManager::install("liftOver")
library(liftOver)
library(jaffelab)
library(sessioninfo)
library(here)

here()
    #[1] "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c"

# ==

## Create gene location map ==================================================

## get GTF
gtf = import("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
## of length 2565061
gtf = gtf[gtf$type == "gene"]
length(gtf)
    ## == nrow(sce.lc)
    
annoTab.full <- as.data.frame(ranges(gtf))
annoTab.full$names <- mcols(gtf)$gene_id
annoTab.full$symbol <- mcols(gtf)$gene_name
annoTab.full$seqlevel <- as.character(seqnames(gtf))
annoTab.full$strand <- as.character(strand(gtf))
annoTab.full <- annoTab.full[ ,c("names", "seqlevel", "start", "end", "strand", "symbol")]
table(annoTab.full$seqlevel)
    #       chr1      chr10      chr11      chr12      chr13      chr14      chr15 
    #       3410       1394       2065       1928        788       1474       1267 
    #      chr16      chr17      chr18      chr19       chr2      chr20      chr21 
    #       1649       1992        781       2027       2541        964        555 
    #      chr22       chr3       chr4       chr5       chr6       chr7       chr8 
    #        901       1893       1533       1811       1827       1686       1495 
    #       chr9       chrM       chrX       chrY GL000009.2 GL000194.1 GL000195.1 
    #       1319         13       1148        111          1          2          2 
    # GL000205.2 GL000213.1 GL000218.1 GL000219.1 KI270711.1 KI270713.1 KI270721.1 
    #          1          1          1          1          1          2          1 
    # KI270726.1 KI270727.1 KI270728.1 KI270731.1 KI270734.1 
    #          2          4          6          1          3 

#dir.create(here("code","magma","annotation"))
write.table(annoTab.full, file=here("code","magma","annotation",
                                    "GRCh38_gencode.v32_Ensembl98_GENES_all-36601.gene.loc"), sep="\t",
            row.names=F, col.names=F, quote=F)


## Define all expressed genes
load(here("processed_data","SCE",
          "markers-stats_LC-n3_findMarkers_19cellTypes.rda"), verbose=T)
    # markers.lc.t.pw, markers.lc.t.1vAll, medianNon0.lc.19

# Just take first entry (stats are computed for all genes/cell class)
expressedGenes <- rownames(markers.lc.t.1vAll[["Astro"]][["Astro_enriched"]])



### Lift coordinates from GRCh38 (hg38) > GRCh37 (hg19) ===
  # (because GWAS are still run w/ summary stats in hg19 coords)
  # - Also only use those for expressed genes - otherwise MTC penalizes stats


gene_df = read.delim(here("code","magma","annotation",
                          "GRCh38_gencode.v32_Ensembl98_GENES_all-36601.gene.loc"),
                     header=FALSE)
colnames(gene_df)= c("GeneID", "Chr", "Start", "End", "Strand", "Symbol")

# add 'chr' prefix
gene_df$Chr = paste0("chr", gene_df$Chr)
gr = makeGRangesFromDataFrame(gene_df, keep=TRUE)
names(gr) = gr$GeneID

## Lift to hg19
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
lifted_list = range(liftOver(gr, ch))
    #Discarding unchained sequences: GL000009.2, GL000194.1, GL000195.1, GL000205.2,
    #GL000213.1, GL000218.1, GL000219.1, KI270711.1, KI270713.1, KI270721.1, KI270726.1,
    #KI270727.1, KI270728.1, KI270731.1, KI270734.1

table(lengths(lifted_list))
#  0     1     2     3     4     5     7    11    14    16    21    31
# 76 33339    99     9     6     2     2     1     1     1     1     1
lifted_list = lifted_list[lengths(lifted_list) == 1]

lifted = unlist(lifted_list)
lifted_df = as.data.frame(lifted)
lifted_df$Symbol = gene_df$Symbol[match(rownames(lifted_df), gene_df$GeneID)]

# Check:
lifted_df[lifted_df$Symbol=="GABRQ", ]
    #                seqnames     start       end width strand Symbol
    #ENSG00000268089     chrX 151806349 151826007 19659      +  GABRQ
    #   ^ good - in hg19, has a different EnsemblID ("ENSG00000147402")

    # And in hg38 coords:
    gene_df[gene_df$Symbol=="GABRQ", ]
    #               GeneID  Chr     Start       End Strand Symbol
    #36340 ENSG00000268089 chrX 152637895 152657542      +  GABRQ

    
## Rearrange cols of lifted_df / match format ===
table(expressedGenes %in% rownames(lifted_df))
    # FALSE  TRUE
    #   157 32066

head(lifted_df, n=3)
#                 seqnames  start    end width strand      Symbol
# ENSG00000243485     chr1  29554  31109  1556      + MIR1302-2HG
# ENSG00000237613     chr1  34554  36081  1528      -     FAM138A
# ENSG00000186092     chr1  65419  71585  6167      +       OR4F5

head(gene_df, n=3)
#            GeneID  Chr  Start    End Strand      Symbol
# 1 ENSG00000243485 chr1  29554  31109      + MIR1302-2HG
# 2 ENSG00000237613 chr1  34554  36081      -     FAM138A
# 3 ENSG00000186092 chr1  65419  71585      +       OR4F5     <- match this format for MAGMA input

lifted_df$width <- NULL # don't need
lifted_df$GeneID <- rownames(lifted_df)
lifted_df <- lifted_df[ ,c(6,1:5)]

# Remove the 'chr' prefix
lifted_df$seqnames <- ss(as.character(lifted_df$seqnames),"chr",2)

# Finally subset for those expressed genes
lifted_df <- lifted_df[intersect(rownames(lifted_df), expressedGenes), ]

## Write out for MAGMA SNP annotation step
write.table(lifted_df, file=here("code","magma","annotation",
                                    "GRCh38_gencode.v32_Ensembl98_LIFTED-to-hg19_expressedGenes.gene.loc"), sep="\t",
            row.names=F, col.names=F, quote=F)




### GWAS SNP location files ====================================================
  # ** As per manual, the column order should be: SNP ID, chromosome, bp position

## GWAS for AD (PGC-ALZ, IGAP, ADSP, UKB meta-meta-analysis: Jansen, et al.2019) ===
sumStats.AD <- read.table(here("code","magma","GWAS_Results","AD_sumstats_Jansenetal_2019sept.txt"), header=T)
class(sumStats.AD)
dim(sumStats.AD)
    #[1] 13367299       14
head(sumStats.AD)
unique(sumStats.AD$CHR)  # 1:22 (no MT genes)
snploc.AD <- sumStats.AD[ ,c("SNP", "CHR", "BP")]
write.table(snploc.AD, file=here("code","magma","GWAS_Results","Alzheimers_PGC-IGAP-ADSP-UKB_2019.snploc"),
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".AD"))



## GWAS for ADHD (PGC x iPSYCH ===
sumStats.ADHD <- read.table(here("code","magma","GWAS_Results",
                                 "daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta"),
                            header=T)
class(sumStats.ADHD)
dim(sumStats.ADHD)
    #[1] 8094094      19
head(sumStats.ADHD)
unique(sumStats.ADHD$CHR)  # 1:22 as well
snploc.ADHD <- sumStats.ADHD[ ,c("SNP", "CHR", "BP")]
write.table(snploc.ADHD, file=here("code","magma","GWAS_Results","ADHD_PGC_2018.snploc"),
            sep="\t", col.names=T, row.names=F, quote=F)
rm(list=ls(pattern=".ADHD"))



## GWAS for Parkinson's ===
sumStats.PD <- read.table(here("code","magma","GWAS_Results",
                                 "nallsEtAl2019_excluding23andMe_allVariants.tab"),
                            header=T)
class(sumStats.PD)
dim(sumStats.PD)
    #[1] 17510617        9
head(sumStats.PD)
    #             SNP A1 A2   freq       b     se      p N_cases N_controls
    # 1 chr11:88249377  T  C 0.9931  0.1575 0.1783 0.3771    7161       5356
    # 2  chr1:60320992  A  G 0.9336  0.0605 0.0456 0.1846   26421     442271
    # 3  chr2:18069070  T  C 0.9988 -0.6774 1.3519 0.6163     582        905

    # *** strangely, no rsIDs are provided - instead, the chr:bp position...
    #     MAGMA requires these three columns in the SNP location file, so will need
    #     reach out to the authors if this can be provided

# --> Let's go ahead and try using the provided SNP ID, but will have to create
#     the CHR & SNP positions - just call it 'test_[...].snploc'

sumStats.PD$CHR <- ss(sumStats.PD$SNP,":",1)
sumStats.PD$CHR <- ss(sumStats.PD$CHR, "chr", 2)

sumStats.PD$BP <- ss(sumStats.PD$SNP,":",2)

unique(sumStats.PD$CHR)  # 1:22 as well
snploc.PD <- sumStats.PD[ ,c("SNP", "CHR", "BP")]
write.table(snploc.PD, file=here("code","magma","GWAS_Results","test_Parkinsons_NallsEtAl_Lancet2019.snploc"),
            sep="\t", col.names=T, row.names=F, quote=F)

## Create an 'Neff' using METAL's recommended computation for meta-GWAS (https://doi.org/10.1093/bioinformatics/btq340)
 #      instead of sum(N_cases, N_controls)
sumStats.PD$N_effective <- 4/(1/sumStats.PD$N_cases + 1/sumStats.PD$N_controls)

# Save
write.table(sumStats.PD, file=here("code","magma","GWAS_Results",
                                   "nallsEtAl2019_excluding23andMe_allVariants_Neff-ADDED-MNT.tab"),
            sep="\t", col.names=T, row.names=F, quote=F)


rm(list=ls(pattern=".PD"))





### Gene marker sets (analagous to layer marker setup) ==============
  # We'll use the 'enriched' stats--i.e. '1vAll' with the medianNon0
  # Remember from `scran::findMarkers()`:
  #           - the 'log.FDR' are on the natural log scale, NOT BASE 10
  
library(SingleCellExperiment)

## Load marker stats
load(here("processed_data","SCE",
          "markers-stats_LC-n3_findMarkers_19cellTypes.rda"), verbose=T)
    # markers.lc.t.1vAll, medianNon0.lc

# Already in Ensembl IDs?
head(rownames(markers.lc.t.1vAll[["Astro"]][["Astro_enriched"]]))

## These stats now have both an '_enriched' & '_depleted' result - take the '_enriched'
sapply(markers.lc.t.1vAll, names)
markers.lc.enriched <- lapply(markers.lc.t.1vAll, function(x){x[[2]]})
    
sapply(markers.lc.enriched, function(x){table(x$log.FDR < log(1e-6) & x$non0median==TRUE)})
    #       Astro Endo.Mural Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A
    # FALSE 32016      31388   29649   29578   31468   30078   31652   31974   31952
    # TRUE    207        835    2574    2645     755    2145     571     249     271

    #       Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Micro Neuron.5HT Neuron.NE Oligo
    # FALSE   28191   31347   31202   31528   31772 32085      31706     31861 31085
    # TRUE     4032     876    1021     695     451   138        517       362  1138

    #         OPC
    # FALSE 31906
    # TRUE    317


lc.markerSet <- data.frame()
for(i in names(markers.lc.enriched)){
  lc.i <- data.frame(Set=rep(i, sum(markers.lc.enriched[[i]]$log.FDR < log(1e-6) &
                                         markers.lc.enriched[[i]]$non0median==TRUE)),
             Gene=rownames(markers.lc.enriched[[i]])[markers.lc.enriched[[i]]$log.FDR < log(1e-6) &
                                                          markers.lc.enriched[[i]]$non0median==TRUE])
  lc.markerSet <- rbind(lc.markerSet, lc.i)
}

head(lc.markerSet)
table(lc.markerSet$Set)  # looks good

nrow(lc.markerSet)  # [1] 19799
    # and btw:
    length(unique(lc.markerSet$Gene)) #[1] 7919

# Write out
write.table(lc.markerSet, file=here("code","magma","/lcMarkerSets_fdr1e-6.txt"), sep="\t",
            row.names=F, col.names=T, quote=F)





## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    # [1] "2022-05-29 16:56:28 EDT"
proc.time()
    #    user   system  elapsed 
    #  46.493    1.946 1202.705 
options(width = 120)
session_info()
    #─ Session info ─────────────────────────────────────────────────────────────────────────────
    # setting  value
    # version  R version 4.1.2 Patched (2021-11-04 r81138)
    # os       CentOS Linux 7 (Core)
    # system   x86_64, linux-gnu
    # ui       X11
    # language (EN)
    # collate  en_US.UTF-8
    # ctype    en_US.UTF-8
    # tz       US/Eastern
    # date     2022-05-29
    # pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
    # 
    # ─ Packages ─────────────────────────────────────────────────────────────────────────────────
    # package                           * version  date (UTC) lib source
    # AnnotationDbi                     * 1.56.2   2021-11-09 [2] Bioconductor
    # assertthat                          0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
    # Biobase                           * 2.54.0   2021-10-26 [2] Bioconductor
    # BiocFileCache                       2.2.1    2022-01-23 [2] Bioconductor
    # BiocGenerics                      * 0.40.0   2021-10-26 [2] Bioconductor
    # BiocIO                              1.4.0    2021-10-26 [2] Bioconductor
    # BiocManager                         1.30.18  2022-05-18 [2] CRAN (R 4.1.2)
    # BiocParallel                        1.28.3   2021-12-09 [2] Bioconductor
    # biomaRt                             2.50.3   2022-02-03 [2] Bioconductor
    # Biostrings                          2.62.0   2021-10-26 [2] Bioconductor
    # bit                                 4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
    # bit64                               4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
    # bitops                              1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
    # blob                                1.2.3    2022-04-10 [2] CRAN (R 4.1.2)
    # BSgenome                            1.62.0   2021-10-26 [2] Bioconductor
    # cachem                              1.0.6    2021-08-19 [2] CRAN (R 4.1.2)
    # cli                                 3.3.0    2022-04-25 [2] CRAN (R 4.1.2)
    # crayon                              1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
    # curl                                4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
    # DBI                                 1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
    # dbplyr                              2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
    # DelayedArray                        0.20.0   2021-10-26 [2] Bioconductor
    # digest                              0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
    # dplyr                               1.0.9    2022-04-28 [2] CRAN (R 4.1.2)
    # ellipsis                            0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
    # fansi                               1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
    # fastmap                             1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
    # filelock                            1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
    # fs                                  1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
    # gargle                              1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
    # generics                            0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
    # GenomeInfoDb                      * 1.30.1   2022-01-30 [2] Bioconductor
    # GenomeInfoDbData                    1.2.7    2021-11-01 [2] Bioconductor
    # GenomicAlignments                   1.30.0   2021-10-26 [2] Bioconductor
    # GenomicFeatures                   * 1.46.5   2022-02-27 [2] Bioconductor
    # GenomicRanges                     * 1.46.1   2021-11-18 [2] Bioconductor
    # glue                                1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
    # GO.db                             * 3.14.0   2021-11-01 [2] Bioconductor
    # googledrive                         2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
    # graph                               1.72.0   2021-10-26 [2] Bioconductor
    # gwascat                           * 2.26.0   2021-10-26 [1] Bioconductor
    # here                              * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
    # hms                                 1.1.1    2021-09-26 [2] CRAN (R 4.1.2)
    # Homo.sapiens                      * 1.3.1    2021-11-04 [2] Bioconductor
    # httr                                1.4.3    2022-05-04 [2] CRAN (R 4.1.2)
    # IRanges                           * 2.28.0   2021-10-26 [2] Bioconductor
    # jaffelab                          * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
    # KEGGREST                            1.34.0   2021-10-26 [2] Bioconductor
    # lattice                             0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
    # lifecycle                           1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
    # liftOver                          * 1.18.0   2021-10-29 [1] Bioconductor
    # limma                               3.50.3   2022-04-07 [2] Bioconductor
    # magrittr                            2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
    # MASS                                7.3-56   2022-03-23 [3] CRAN (R 4.1.2)
    # Matrix                              1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
    # MatrixGenerics                      1.6.0    2021-10-26 [2] Bioconductor
    # matrixStats                         0.62.0   2022-04-19 [2] CRAN (R 4.1.2)
    # memoise                             2.0.1    2021-11-26 [2] CRAN (R 4.1.2)
    # org.Hs.eg.db                      * 3.14.0   2021-11-01 [2] Bioconductor
    # OrganismDbi                       * 1.36.0   2021-10-26 [2] Bioconductor
    # pillar                              1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
    # pkgconfig                           2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
    # png                                 0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
    # prettyunits                         1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
    # progress                            1.2.2    2019-05-16 [2] CRAN (R 4.1.0)
    # purrr                               0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
    # R6                                  2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
    # rafalib                           * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
    # rappdirs                            0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
    # RBGL                                1.70.0   2021-10-26 [2] Bioconductor
    # RColorBrewer                        1.1-3    2022-04-03 [2] CRAN (R 4.1.2)
    # Rcpp                                1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
    # RCurl                               1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
    # readr                               2.1.2    2022-01-30 [2] CRAN (R 4.1.2)
    # restfulr                            0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
    # rjson                               0.2.21   2022-01-09 [2] CRAN (R 4.1.2)
    # rlang                               1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
    # rprojroot                           2.0.3    2022-04-02 [2] CRAN (R 4.1.2)
    # Rsamtools                           2.10.0   2021-10-26 [2] Bioconductor
    # RSQLite                             2.2.14   2022-05-07 [2] CRAN (R 4.1.2)
    # rtracklayer                       * 1.54.0   2021-10-26 [2] Bioconductor
    # S4Vectors                         * 0.32.4   2022-03-24 [2] Bioconductor
    # segmented                           1.5-0    2022-04-11 [1] CRAN (R 4.1.2)
    # sessioninfo                       * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
    # snpStats                            1.44.0   2021-10-26 [2] Bioconductor
    # stringi                             1.7.6    2021-11-29 [2] CRAN (R 4.1.2)
    # stringr                             1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
    # SummarizedExperiment                1.24.0   2021-10-26 [2] Bioconductor
    # survival                            3.3-0    2022-03-01 [3] CRAN (R 4.1.2)
    # tibble                              3.1.7    2022-05-03 [2] CRAN (R 4.1.2)
    # tidyselect                          1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
    # TxDb.Hsapiens.UCSC.hg19.knownGene * 3.2.2    2021-05-21 [2] Bioconductor
    # tzdb                                0.3.0    2022-03-28 [2] CRAN (R 4.1.2)
    # utf8                                1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
    # VariantAnnotation                   1.40.0   2021-10-26 [2] Bioconductor
    # vctrs                               0.4.1    2022-04-13 [2] CRAN (R 4.1.2)
    # XML                                 3.99-0.9 2022-02-24 [2] CRAN (R 4.1.2)
    # xml2                                1.3.3    2021-11-30 [2] CRAN (R 4.1.2)
    # XVector                             0.34.0   2021-10-26 [2] Bioconductor
    # yaml                                2.3.5    2022-02-21 [2] CRAN (R 4.1.2)
    # zlibbioc                            1.40.0   2021-10-26 [2] Bioconductor
    # 
    # [1] /users/ntranngu/R/4.1.x
    # [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
    # [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
    # 
    # ────────────────────────────────────────────────────────────────────────────────────────────

