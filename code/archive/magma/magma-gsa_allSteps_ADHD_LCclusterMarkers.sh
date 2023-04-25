#!/bin/bash
#$ -cwd
#$ -N magma_steps1-3-gsa_ADHD
#$ -o ./logs/magma-gsa_steps1-3_ADHD_MNT06Jun2022.o
#$ -e ./logs/magma-gsa_steps1-3_ADHD_MNT06Jun2022.e
#$ -l bluejay,mem_free=16G,h_vmem=20G


echo "**** Job starts ****"
date

## Load MAGMA
module load magma/1.10

## Set some variables/paths
model="snp-wise"
ANNO=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/annotation/GRCh38_gencode.v32_Ensembl98_LIFTED-to-hg19_expressedGenes.gene.loc
BFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur

setcol=1
genecol=2
gs_lc=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/marker_set_NE.txt

SUMMSTATS=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/GWAS_Results/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta


## Step 1 - Annotation (SNP : gene mapping)
magma --annotate window=35,10 --snp-loc ./GWAS_Results/ADHD_PGC_2019.snploc --gene-loc $ANNO --out SNP_Data/ADHD_Demontis2019_LC

## Step 2 - Gene analysis (from SNP-level summary stats)
magma --bfile $BFILE --gene-annot SNP_Data/ADHD_Demontis2019_LC.genes.annot --pval $SUMMSTATS use=SNP,P ncol=Neff --gene-model ${model} --out SNP_Data/ADHD_Demontis2019_LC_${model}

## Step 3 - Gene set analyses (using gene-level output)
magma --gene-results SNP_Data/ADHD_Demontis2019_LC_snp-wise.genes.raw --set-annot $gs_lc gene-col=${genecol} set-col=${setcol} --out Results/lc_ADHD


echo "**** Job ends ****"
date

