#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_NextSeq.sh

############################
# Script to run Space Ranger
############################

# locations of files:
# -------------------
# spaceranger reference: /dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A
# fastq (MiSeq): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ
# fastq (NextSeq): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ_NextSeq/
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw
# manual alignment files from Loupe: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/images_raw_align_json

# summary spreadsheet: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Visium LC pilot_072120 Master.xlsx
# - contains sample ID, sample name, slide serial number, capture area ID


# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
mkdir -p ../outputs/NextSeq
cd ../outputs/NextSeq


# run spaceranger count for each sample
# note: only 2 samples for NextSeq

spaceranger count \
--id=LC_1 \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ_NextSeq/LC_1 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_2.tif \
--slide=V19B23-076 \
--area=B1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/images_raw_align_json/V19B23-076-B1.json \
--jobmode=local \
--localcores=8 \
--localmem=64

spaceranger count \
--id=LC_2 \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ_NextSeq/LC_2 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_3.tif \
--slide=V19B23-076 \
--area=C1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/images_raw_align_json/V19B23-076-C1.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd

