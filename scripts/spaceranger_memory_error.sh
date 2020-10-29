#!/bin/bash
#$ -l mem_free=11G,h_vmem=22G,h_fsize=100G
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger_memory_error.sh


# --------------------------------------------------------------------------
# Minimal script to reproduce Space Ranger memory error when using module on 
# JHPCE cluster
# 
# error OCCURS when using module
# error DOES NOT OCCUR when NOT using module (comment out lines 18 and 19)
# --------------------------------------------------------------------------

# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
cd ..
mkdir -p outputs_memory_error
cd outputs_memory_error


# run spaceranger count
spaceranger count \
--id=DLPFC \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ/DLPFC \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/Lieber-Institute_OTS-20-7043_1_1.tif \
--slide=V19B23-076 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/raw_manual_align_json/V19B23-076-A1.json \
--jobmode=local \
--localcores=1 \
--localmem=10

