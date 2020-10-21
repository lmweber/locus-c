#!/bin/bash

############################
# Script to run Space Ranger
############################

# example code to run Space Ranger
# adapted from Andrew Jaffe's script 'pilot_align.sh'

# run on JHPCE cluster
# qrsh -pe local 6 -l mem_free=7G,h_vmem=14G,h_fsize=100G -now n

# using local installation of spaceranger in home directory


# run spaceranger count
spaceranger count --id=LC_1 \
--transcriptome=../spaceranger/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ/LC_1 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Snapshots_Plus_Loupe/Lieber-Institute_OTS-20-7043_Pt1.jpg \
--slide=V19B23-076 --area=B1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Snapshots_Plus_Loupe/V19B23-076-B1.json \
--localcores=6 \
--localmem=40 \
--localvmem=80

