#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l mem_free=8G,h_vmem=10G,h_fsize=300G

# qsub code/spaceranger/filename.sh


############################
# Script to run Space Ranger
############################

# Round 3 LC samples: KMay_2021-07-09


# load spaceranger module (if required)
# module use /jhpce/shared/jhpce/modulefiles/libd
# module load spaceranger


# run in outputs directory
cwd=$(pwd)
mkdir -p processed_data/spaceranger/KMay_2021-07-09
cd processed_data/spaceranger/KMay_2021-07-09


# run spaceranger count for each sample

spaceranger count \
--id=Br6522_LC_round3 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/KMay_2021-07-09/Br6522_LC_round3 \
--image=../../../images/round3/2021-06-15_V10B01-003_LC_A1.tif \
--slide=V10B01-003 \
--area=A1 \
--loupe-alignment=../../../alignment/round3/V10B01-003-A1.json \
--nosecondary \
--jobmode=local \
--localcores=8 \
--localmem=64 \
--localvmem=80


spaceranger count \
--id=Br8079_LC_round3 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/KMay_2021-07-09/Br8079_LC_round3 \
--image=../../../images/round3/2021-06-15_V10B01-003_LC_B1.tif \
--slide=V10B01-003 \
--area=B1 \
--loupe-alignment=../../../alignment/round3/V10B01-003-B1.json \
--nosecondary \
--jobmode=local \
--localcores=8 \
--localmem=64 \
--localvmem=80


spaceranger count \
--id=Br2701_LC_round3 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/KMay_2021-07-09/Br2701_LC_round3 \
--image=../../../images/round3/2021-06-15_V10B01-003_LC_C1.tif \
--slide=V10B01-003 \
--area=C1 \
--loupe-alignment=../../../alignment/round3/V10B01-003-C1.json \
--nosecondary \
--jobmode=local \
--localcores=8 \
--localmem=64 \
--localvmem=80


spaceranger count \
--id=Br8153_LC_round3 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/KMay_2021-07-09/Br8153_LC_round3 \
--image=../../../images/round3/2021-06-15_V10B01-003_LC_D1.tif \
--slide=V10B01-003 \
--area=D1 \
--loupe-alignment=../../../alignment/round3/V10B01-003-D1.json \
--nosecondary \
--jobmode=local \
--localcores=8 \
--localmem=64 \
--localvmem=80


# restore working directory
cd $cwd

