#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=12G,h_fsize=100G
#$ -V
#$ -cwd

############################
# Script to run Space Ranger
############################

# locations of files:
# -------------------
# spaceranger reference: /dcl02/leased/shicks/spaceranger/refdata-gex-GRCh38-2020-A
# fastq: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ
# images (screenshots): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Snapshots_Plus_Loupe
# - image files and Loupe manual alignment json files
# - add argument --loupe-alignment to use manual alignment files (see Andrew Jaffe's script 'pilot_align.sh')
# - see 'pilot_align.sh' for which screenshot matches to which sample
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw
# - no manual alignment files (use automatic alignment instead)
# summary spreadsheet: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Visium LC pilot_072120 Master.xlsx
# - contains sample ID, sample name, slide serial number, capture area ID

# run on JHPCE cluster
# qsub -V -cwd -pe local 8 -l mem_free=30G,h_vmem=40G,h_fsize=300G scripts/spaceranger.sh

# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
cd ../outputs

# run spaceranger count
spaceranger count \
--id=DLPFC \
--transcriptome=/dcl02/leased/shicks/spaceranger/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ/DLPFC \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_1.tif \
--slide=V19B23-076 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/json_after_manual_loupe/V19B23-076-A1.json \
--localcores=8 \
--localmem=64

spaceranger count \
--id=LC_1 \
--transcriptome=/dcl02/leased/shicks/spaceranger/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ/LC_1 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_2.tif \
--slide=V19B23-076 \
--area=B1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/json_after_manual_loupe/V19B23-076-B1.json \
--localcores=8 \
--localmem=64

spaceranger count \
--id=LC_2 \
--transcriptome=/dcl02/leased/shicks/spaceranger/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ/LC_2 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_3.tif \
--slide=V19B23-076 \
--area=C1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/json_after_manual_loupe/V19B23-076-C1.json \
--localcores=8 \
--localmem=64

spaceranger count \
--id=LC_3 \
--transcriptome=/dcl02/leased/shicks/spaceranger/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ/LC_3 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_4.tif \
--slide=V19B23-076 \
--area=D1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/json_after_manual_loupe/V19B23-076-D1.json \
--localcores=8 \
--localmem=64

# restore working directory
cd $cwd

