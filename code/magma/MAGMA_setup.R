### for MAGMA with locus coeruleus snRNA-seq populations
  # MNT 27May2022

library(rtracklayer)
library(GenomicRanges)
#BiocManager::install("liftOver")
library(liftOver)
library(jaffelab)
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
# TODO





## GWAS for Parkinson's ===
# TODO





### Gene marker sets (analagous to layer marker setup) ==============
  # TODO (pasted from former work)
  # We'll just use the 'enriched' stats--i.e. '1vAll'
  # Update May 2021: Apply a filter for markers that median expression in respective subcluster != 0
  #                  (already computed and added to a different .rda file)
  # UPDATE FOR REVISION: the 'log.FDR' are on the natural log scale, NOT BASE 10
  #                      -> Will switch to log.FDR <= log(1e-6), since this is ~= log10(1e-12)
  #                         (log(1e-12) == -27.63, way too strict)
  
library(SingleCellExperiment)

## DLPFC ===
load("rdas/revision/markers-stats_DLPFC-n3_findMarkers-SN-LEVEL_MNT_v2_2021.rda", verbose=T)
    # markers.dlpfc.t.1vAll, medianNon0.dlpfc

# For EnsemblIDs
load("rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_MNT2021.rda", verbose=T)
    
## These stats now have both an '_enriched' & '_depleted' result - take the '_enriched'
names(markers.dlpfc.t.1vAll[[1]])
markers.dlpfc.enriched <- lapply(markers.dlpfc.t.1vAll, function(x){x[[2]]})
    
sapply(markers.dlpfc.enriched, function(x){table(x$log.FDR < log(1e-6) & x$non0median==TRUE)})
    #       Astro Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Inhib_A Inhib_B
    # FALSE 28620   26037   26423   25833   27548   26947   26999   27592   27788
    # TRUE    690    3273    2887    3477    1762    2363    2311    1718    1522
    
    #      Inhib_C Inhib_D Inhib_E Inhib_F Macrophage Micro Mural Oligo   OPC Tcell
    # FALSE   27371   26679   29211   29236      29110 28734 29148 28444 28393 29194
    # TRUE     1939    2631      99      74        200   576   162   866   917   116


dlpfc.markerSet <- data.frame()
for(i in names(markers.dlpfc.enriched)){
  dlpfc.i <- data.frame(Set=rep(i, sum(markers.dlpfc.enriched[[i]]$log.FDR < log(1e-6) &
                                         markers.dlpfc.enriched[[i]]$non0median==TRUE)),
             Gene=rownames(markers.dlpfc.enriched[[i]])[markers.dlpfc.enriched[[i]]$log.FDR < log(1e-6) &
                                                          markers.dlpfc.enriched[[i]]$non0median==TRUE])
  dlpfc.markerSet <- rbind(dlpfc.markerSet, dlpfc.i)
}

head(dlpfc.markerSet)
table(dlpfc.markerSet$Set)  # looks good

dlpfc.markerSet$Gene <- rowData(sce.dlpfc)$gene_id[match(dlpfc.markerSet$Gene, rownames(sce.dlpfc))]
    # (btw:)
    length(unique(dlpfc.markerSet$Gene)) #[1] 6925

# Write out
write.table(dlpfc.markerSet, file="./MAGMA/dlpfcMarkerSets_fdr1e-6.txt", sep="\t",
            row.names=F, col.names=T, quote=F)





## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    #
proc.time()
    #
options(width = 120)
session_info()


