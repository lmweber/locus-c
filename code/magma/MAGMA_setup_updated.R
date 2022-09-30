###############################################
# Setup for MAGMA with LC snRNA-seq populations
# code from Matthew N Tran
# adapted by Lukas Weber, Sep 2022
###############################################

# module load conda_R/4.2


library(rtracklayer)
library(GenomicRanges)
library(liftOver)
library(jaffelab)
library(here)
library(readr)


## Create gene location map ==================================================

## get GTF
gtf = import("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
gtf = gtf[gtf$type == "gene"]
length(gtf)


annoTab.full <- as.data.frame(ranges(gtf))
annoTab.full$names <- mcols(gtf)$gene_id
annoTab.full$symbol <- mcols(gtf)$gene_name
annoTab.full$seqlevel <- as.character(seqnames(gtf))
annoTab.full$strand <- as.character(strand(gtf))
annoTab.full <- annoTab.full[, c("names", "seqlevel", "start", "end", "strand", "symbol")]

table(annoTab.full$seqlevel)

write.table(
  annoTab.full, 
  file = here("code", "magma", "annotation", 
              "GRCh38_gencode.v32_Ensembl98_GENES_all-36601.gene.loc"), 
  sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)


## Define all expressed genes
## note: only genes that pass filtering
sce <- readRDS(here("processed_data", "SCE", "sce_clustering_merged.rds"))
dim(sce)
expressedGenes <- rownames(sce)
str(expressedGenes)


### Lift coordinates from GRCh38 (hg38) > GRCh37 (hg19) ===
  # (because GWAS are still run w/ summary stats in hg19 coords)
  # - Also only use those for expressed genes - otherwise MTC penalizes stats

gene_df <- read.delim(here("code", "magma", "annotation", 
                           "GRCh38_gencode.v32_Ensembl98_GENES_all-36601.gene.loc"), 
                      header = FALSE)
colnames(gene_df) <- c("GeneID", "Chr", "Start", "End", "Strand", "Symbol")

gr <- makeGRangesFromDataFrame(gene_df, keep = TRUE)
names(gr) <- gr$GeneID

## Lift to hg19
path <- system.file(package = "liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)
lifted_list <- range(liftOver(gr, ch))
    #Discarding unchained sequences: GL000009.2, GL000194.1, GL000195.1, GL000205.2, 
    #GL000213.1, GL000218.1, GL000219.1, KI270711.1, KI270713.1, KI270721.1, KI270726.1, 
    #KI270727.1, KI270728.1, KI270731.1, KI270734.1

table(lengths(lifted_list))
lifted_list <- lifted_list[lengths(lifted_list) == 1]

lifted <- unlist(lifted_list)
lifted_df <- as.data.frame(lifted)
lifted_df$Symbol = gene_df$Symbol[match(rownames(lifted_df), gene_df$GeneID)]

# Check:
lifted_df[lifted_df$Symbol == "GABRQ", ]
    #                seqnames     start       end width strand Symbol
    #ENSG00000268089     chrX 151806349 151826007 19659      +  GABRQ
    #   ^ good - in hg19, has a different EnsemblID ("ENSG00000147402")

    # And in hg38 coords:
gene_df[gene_df$Symbol == "GABRQ", ]
    #               GeneID  Chr     Start       End Strand Symbol
    #36340 ENSG00000268089 chrX 152637895 152657542      +  GABRQ


## Rearrange cols of lifted_df / match format ===
table(expressedGenes %in% rownames(lifted_df))

head(lifted_df, 3)
head(gene_df, 3)

lifted_df$width <- NULL # don't need
lifted_df$GeneID <- rownames(lifted_df)
lifted_df <- lifted_df[ , c(6,1:5)]

# Remove the 'chr' prefix
lifted_df$seqnames <- ss(as.character(lifted_df$seqnames), "chr", 2)

# Finally subset for those expressed genes
lifted_df <- lifted_df[intersect(rownames(lifted_df), expressedGenes), ]

## Write out for MAGMA SNP annotation step
write.table(
  lifted_df, file = here("code", "magma", "annotation", 
                         "GRCh38_gencode.v32_Ensembl98_LIFTED-to-hg19_expressedGenes.gene.loc"), 
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)


### GWAS SNP location files ====================================================
  # ** As per manual, the column order should be: SNP ID, chromosome, bp position

## GWAS for AD (PGC-ALZ, IGAP, ADSP, UKB meta-meta-analysis: Jansen, et al.2019) ===
sumStats.AD <- read_table(
  here("code", "magma", "GWAS_Results", 
       "AD_sumstats_Jansenetal_2019sept.txt"), 
  col_names = TRUE
)
class(sumStats.AD)
dim(sumStats.AD)
head(sumStats.AD)
unique(sumStats.AD$CHR)  # 1:22 (no MT genes)
snploc.AD <- sumStats.AD[ , c("SNP", "CHR", "BP")]
write.table(
  snploc.AD, file = here("code", "magma", "GWAS_Results", 
                         "Alzheimers_PGC-IGAP-ADSP-UKB_2019.snploc"), 
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


## GWAS for ADHD (PGC x iPSYCH ===
sumStats.ADHD <- read_table(
  here("code", "magma", "GWAS_Results", 
       "daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta"), 
  col_names = TRUE
)
class(sumStats.ADHD)
dim(sumStats.ADHD)
head(sumStats.ADHD)
unique(sumStats.ADHD$CHR)  # 1:22 as well

table(sumStats.ADHD$P < 5e-8) # Reported 304 variants (12 loci) surpassing this threshold
    #   FALSE    TRUE 
    # 8093777     317   - this is pretty good; unfortunately the PD GWAS/supplement doesn't report
    #                     the number of total SNPs reaching genome-wide signif; just 90 [index] SNPs

snploc.ADHD <- sumStats.ADHD[, c("SNP", "CHR", "BP")]
write.table(
  snploc.ADHD, file = here("code", "magma", "GWAS_Results", 
                           "ADHD_PGC_2019.snploc"), 
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


## GWAS for Parkinson's ===
sumStats.PD <- read_table(
  here("code", "magma", "GWAS_Results", 
       "nallsEtAl2019_excluding23andMe_allVariants.tab"), 
  col_names = TRUE
)
class(sumStats.PD)
dim(sumStats.PD)
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

sumStats.PD$CHR <- ss(sumStats.PD$SNP, ":", 1)
sumStats.PD$CHR <- ss(sumStats.PD$CHR, "chr", 2)

sumStats.PD$BP <- ss(sumStats.PD$SNP, ":", 2)

unique(sumStats.PD$CHR)  # 1:22 as well

## **The provided SNP IDs DO have to be compatible with those in the 1000 genomes**

### With SNPlocs package ===
  #   See (https://support.bioconductor.org/p/133105/) for the similar query
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# Following the manual:
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
snpcount(snps)

# Try it with the smallest autosome:
chr22_snps <- snpsBySeqname(snps, "22")
class(chr22_snps)
head(chr22_snps)

chr22_snps <- as.data.frame(chr22_snps)
chr22_snps$chr.bp <- paste0("chr", chr22_snps$seqnames, ":", chr22_snps$pos)

# Inspect with chr22 SNPs in summary stats
table(sumStats.PD$CHR == "22")
sumStats.PD.chr22 <- sumStats.PD[sumStats.PD$CHR == "22", ]

table(sumStats.PD.chr22$SNP %in% chr22_snps$chr.bp)
    #FALSE   TRUE 
    # 6710 226625   - 97.12%   - this is probably ok
    #               - generally 70-75% SNPs map to one gene in MAGMA step 1 in any case

sumStats.PD.keep <- data.frame()

## What about across all snps?
rownames(sumStats.PD) <- sumStats.PD$SNP

for(i in seqnames(snps)){
  cat("Querying chr: ", i, "...\n")
  temp.snps <- as.data.frame(snpsBySeqname(snps, i))
  temp.snps$chr.bp <- paste0("chr", temp.snps$seqnames, ":", temp.snps$pos)
  
  # Subset sumStats to quantify % intersecting per chromosome
  sumStats.temp <- sumStats.PD[sumStats.PD$CHR == i, ]
  
  cat(paste0("\tPercent of summary statistics SNPs in chr:", i, " with rsIDs:\n"))
  print(table(sumStats.temp$SNP %in% temp.snps$chr.bp)["TRUE"] / nrow(sumStats.temp) * 100)
  cat("\n")
  
  temp.df <- sumStats.temp[intersect(sumStats.temp$SNP, temp.snps$chr.bp), ]
  temp.df$rsID <- temp.snps$RefSNP_id[match(temp.df$SNP, temp.snps$chr.bp)]
  
  # rbind
  sumStats.PD.keep <- rbind(sumStats.PD.keep, temp.df)
}

# Total SNPs in summary stats with rsIDs:
dim(sumStats.PD.keep)
    #[1] 16934405       12
# Overall % SNPs keeping
(nrow(sumStats.PD.keep) / nrow(sumStats.PD)) * 100
    #[1] 96.70936

snploc.PD <- sumStats.PD.keep[, c("rsID", "CHR", "BP")]
write.table(
  snploc.PD, file = here("code", "magma", "GWAS_Results", 
                         "Parkinsons_NallsEtAl_Lancet2019.snploc"), 
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## Create an 'Neff' using METAL's recommended computation for meta-GWAS (https://doi.org/10.1093/bioinformatics/btq340)
 #      instead of sum(N_cases, N_controls)
sumStats.PD.keep$N_effective <- 4/(1/sumStats.PD.keep$N_cases + 1/sumStats.PD.keep$N_controls)

# Save
write.table(
  sumStats.PD.keep, file = here("code", "magma", "GWAS_Results", 
                                "nallsEtAl2019_excluding23andMe_allVariants_Neff-and-rsID-ADDED-MNT.tab"), 
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



### Updated to here ###



### Gene marker sets (analogous to layer marker setup) ==============
  # We'll use the 'enriched' stats--i.e. '1vAll' with the medianNon0
  # Remember from `scran::findMarkers()`:
  #           - the 'log.FDR' are on the natural log scale, NOT BASE 10
  
library(SingleCellExperiment)

## Load marker stats
load(here("processed_data","SCE_MT",
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

