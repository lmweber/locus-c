#!/bin/bash

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
mkdir -p LC_1
cd LC_1
spaceranger count \
--id=LC_1 \
--transcriptome=/dcl02/leased/shicks/spaceranger/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ/LC_1 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_2.tif \
--slide=V19B23-076 \
--area=B1 \
--nosecondary \
--localcores=8 \
--localmem=200 \
--localvmem=400

# restore working directory
cd $cwd

