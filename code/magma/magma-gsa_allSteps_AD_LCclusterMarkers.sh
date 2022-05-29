#!/bin/bash
#$ -cwd
#$ -N magma_steps1-2_AD
#$ -o ./logs/magma-gsa_steps1-2_AD_MNT29May2022.o
#$ -e ./logs/magma-gsa_steps1-2_AD_MNT29May2022.e
#$ -l bluejay,mem_free=16G,h_vmem=20G

echo "**** Job starts ****"
date

## List current modules for reproducibility
module list

## Load MAGMA
module load magma/1.10

## Set some variables/paths
model="snp-wise"
ANNO=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/annotation/GRCh38_gencode.v32_Ensembl98_LIFTED-to-hg19_expressedGenes.gene.loc
BFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur

setcol=1
genecol=2
# TODO gs_lc=

SUMMSTATS=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/magma/GWAS_Results/AD_sumstats_Jansenetal_2019sept.txt


echo "MNT comment 29May2022: just run steps 1-2 to confirm concordance b/tw GWAS SNP and reference SNP IDs"

## Step 1 - Annotation (SNP : gene mapping)
magma --annotate window=35,10 --snp-loc ./GWAS_Results/Alzheimers_PGC-IGAP-ADSP-UKB_2019.snploc --gene-loc $ANNO --out SNP_Data/AD_Jansen2019_LC

## Step 2 - Gene analysis (from SNP-level summary stats)
magma --bfile $BFILE --gene-annot SNP_Data/AD_Jansen2019_LC.genes.annot --pval $SUMMSTATS use=SNP,P ncol=Neff --gene-model ${model} --out SNP_Data/AD_Jansen2019_LC_${model}


## Step 3 - Gene set analyses (using gene-level output)
# magma --gene-results SNP_Data/PD_Nalls2019_LC_snp-wise.genes.raw --set-annot $gs_lc gene-col=${genecol} set-col=${setcol} --out Results/lc_PD


echo "**** Job ends ****"
date
