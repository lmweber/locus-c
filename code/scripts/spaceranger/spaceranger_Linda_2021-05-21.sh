#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub code/scripts/spaceranger/spaceranger_Linda_2021-05-21.sh

############################
# Script to run Space Ranger
############################

# locations of files:
# -------------------
# spaceranger reference: /dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A/
# fastq (first 2 samples MiSeq):    ./fastq/MiSeq_2020-08-12/
# fastq (first 2 samples NextSeq):  ./fastq/NextSeq_2020-09-25/
# fastq (new samples from Linda):   ./fastq/Linda_2021-05-21/
# images (split by sample):         ./images/split/
# image alignment files from Loupe: ./images_align/
# sample information spreadsheets:  ./sample_info/


# load spaceranger module (or use local installation)
#module use /jhpce/shared/jhpce/modulefiles/libd
#module load spaceranger

# run in output directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
mkdir -p processed_data/spaceranger/Linda_2021-05-21
cd processed_data/spaceranger/Linda_2021-05-21


# run spaceranger count for each sample

spaceranger count \
--id=Br8153_LC \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=fastq/Linda_2021-05-21/Br8153_LC \
--image=images/split/V10U24-093_A1_Br8153_LC_1.tif \
--slide=V10U24-093 \
--area=A1 \
--loupe-alignment=images_align/V10U24-093-A1.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=Br5459_LC \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=fastq/Linda_2021-05-21/Br5459_LC \
--image=images/split/V10U24-093_B1_Br5459_LC_2.tif \
--slide=V10U24-093 \
--area=B1 \
--loupe-alignment=images_align/V10U24-093-B1.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=Br2701_LC \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=fastq/Linda_2021-05-21/Br2701_LC \
--image=images/split/V10U24-093_C1_Br2701_LC_3.tif \
--slide=V10U24-093 \
--area=C1 \
--loupe-alignment=images_align/V10U24-093-C1.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd

