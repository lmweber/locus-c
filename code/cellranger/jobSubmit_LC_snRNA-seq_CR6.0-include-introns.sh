#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -pe local 4
#$ -N snRNA-seq_LC_CR6.1
#$ -o logs/snRNA-seq_LC_CR6.1-includeIntrons.$TASK_ID.txt
#$ -e logs/snRNA-seq_LC_CR6.1-includeIntrons.$TASK_ID.txt
#$ -m e
#$ -t 1-3
#$ -tc 3

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## List current modules for reproducibility
module list

## load CellRanger
module load cellranger/6.1.1

## Locate file
FASTQDIR=/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/fastq/snRNA-seq
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' ${FASTQDIR}/sample_libs_info.tsv | awk "NR==${SGE_TASK_ID}")

# Fixed FASTQ file sample prefix ("-" instead of "_")
FIXED=$(awk 'BEGIN {FS="\t"} {print $4}' ${FASTQDIR}/sample_libs_info.tsv | awk "NR==${SGE_TASK_ID}")

echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQDIR}/${SAMPLE} \
    --sample=${FIXED} \
    --jobmode=local \
    --localcores=4 \
    --localmem=20 \
    --include-introns

## Move output
echo "Moving data to new location"
date
mv ${SAMPLE} /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/processed_data/cellranger/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
