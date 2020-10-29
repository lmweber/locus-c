#!/bin/bash
#$ -pe local 8
#$ -l mem_free=9G,h_vmem=18G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger_testing.sh

################################################
# Minimal script to run Space Ranger for testing
################################################

# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
cd ..
mkdir -p outputs_testing
cd outputs_testing


# run spaceranger
spaceranger testrun --id=tiny

