#!/bin/bash
#$ -cwd
#$ -pe local 6
#$ -l mem_free=10G,h_vmem=12G,h_fsize=300G

# qsub code/scripts/spaceranger/filename.sh


############################
# Script to run Space Ranger
############################

# Round 1 LC samples: NextSeqMiSeq


# load spaceranger module (if required)
# module use /jhpce/shared/jhpce/modulefiles/libd
# module load spaceranger


# run in outputs directory
cwd=$(pwd)
mkdir -p processed_data/spaceranger/NextSeqMiSeq
cd processed_data/spaceranger/NextSeqMiSeq


# run spaceranger count for each sample

spaceranger count \
--id=Br6522_LC_1_round1 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/NextSeq_2020-09-25/LC_1,../../../fastq/MiSeq_2020-08-12/LC_1 \
--image=../../../images/round1/Lieber-Institute_OTS-20-7043_1_2.tif \
--slide=V19B23-076 \
--area=B1 \
--loupe-alignment=../../../images_align/V19B23-076-B1.json \
--jobmode=local \
--localcores=6 \
--localmem=60


spaceranger count \
--id=Br6522_LC_2_round1 \
--transcriptome=/dcs04/hicks/data/lweber/data/refdata-gex-GRCh38-2020-A \
--fastqs=../../../fastq/NextSeq_2020-09-25/LC_2,../../../fastq/MiSeq_2020-08-12/LC_2 \
--image=../../../images/round1/Lieber-Institute_OTS-20-7043_1_3.tif \
--slide=V19B23-076 \
--area=C1 \
--loupe-alignment=../../../images_align/V19B23-076-C1.json \
--jobmode=local \
--localcores=6 \
--localmem=60


# restore working directory
cd $cwd

