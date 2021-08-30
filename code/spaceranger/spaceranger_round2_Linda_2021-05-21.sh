#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l mem_free=8G,h_vmem=10G,h_fsize=300G

# qsub code/spaceranger/filename.sh


############################
# Script to run Space Ranger
############################

# Round 2 LC samples: Linda_2021-05-21


# load spaceranger module (if required)
# module use /jhpce/shared/jhpce/modulefiles/libd
# module load spaceranger


# run in outputs directory
cwd=$(pwd)
mkdir -p processed_data/spaceranger/Linda_2021-05-21
cd processed_data/spaceranger/Linda_2021-05-21


# run spaceranger count for each sample

spaceranger count \
--id=Br8153_LC_round2 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/Linda_2021-05-21/Br8153_LC \
--image=../../../images/round2/V10U24-093_A1_Br8153_LC_1.tif \
--slide=V10U24-093 \
--area=A1 \
--loupe-alignment=../../../alignment/round2/V10U24-093-A1.json \
--nosecondary \
--jobmode=local \
--localcores=8 \
--localmem=64 \
--localvmem=80


spaceranger count \
--id=Br5459_LC_round2 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/Linda_2021-05-21/Br5459_LC \
--image=../../../images/round2/V10U24-093_B1_Br5459_LC_2.tif \
--slide=V10U24-093 \
--area=B1 \
--loupe-alignment=../../../alignment/round2/V10U24-093-B1.json \
--nosecondary \
--jobmode=local \
--localcores=8 \
--localmem=64 \
--localvmem=80


spaceranger count \
--id=Br2701_LC_round2 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/Linda_2021-05-21/Br2701_LC \
--image=../../../images/round2/V10U24-093_C1_Br2701_LC_3.tif \
--slide=V10U24-093 \
--area=C1 \
--loupe-alignment=../../../alignment/round2/V10U24-093-C1.json \
--nosecondary \
--jobmode=local \
--localcores=8 \
--localmem=64 \
--localvmem=80


# restore working directory
cd $cwd

