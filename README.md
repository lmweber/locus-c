# Code repository for human locus coeruleus (LC) project analyses

This repository contains code to reproduce analyses and figures in our project on characterizing the gene expression landscape of the human locus coeruleus using single-nucleus RNA-sequencing (snRNA-seq) and spatially-resolved transcriptomics (SRT).


## Contents

Code scripts are in the [code/] directory in the following subdirectories:

- `spaceranger`: run Space Ranger to align sequencing reads (FASTQ files) for SRT data
- `cellranger`: run Cell Ranger to align sequencing reads (FASTQ files) for snRNA-seq data
- `analyses`: R scripts for analyses of SRT data
- `analyses_sn_alt`: R scripts for analyses of snRNA-seq data


# Internal info

## Location of files

Files are located on JHPCE at:

- `/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/`

Previously files were also at:

- `/dcs04/hicks/data/lweber/locus-c/`
- `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/locus-c/`


## FASTQ files

FASTQ files are in the following subdirectories:

- Round 1: `fastq/NextSeq_2020-09-25` and `fastq/MiSeq_2020-08-12`
- Round 2: `fastq/Linda_2021-05-21`
- Round 3: `fastq/KMay_2021-07-09`

Backup of FASTQ files is also at:

- `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/`


## Image files

Image files are in the following directories:

- Round 1: `/images/round1/`
- Round 2: `/images/round2/`
- Round 3: `/images/round3/`

Full set of image files for round 3 is also located at:

- Round 3: `/dcs04/lieber/marmaypag/visiumImages_LIBD001/raw-data/2021-06-15_V10B01-003/`


## Alignment files

Image alignment files from manual alignment in Loupe are in the following directories:

- Round 1: `/alignment/round1/`
- Round 2: `/alignment/round2/`
- Round 3: `/alignment/round3/`


## Other input files

- `inputs/`: other external input files


## Code scripts

- `code/spaceranger/`: scripts to run Space Ranger
- `code/other/`: various other shell scripts
- `code/analyses/`: main R analysis scripts
- `code/exploratory/`: RMarkdown scripts used for initial exploratory analyses


## Output files

- `processed_data/spaceranger/`: Space Ranger output files
- `processed_data/VistoSeg/`: VistoSeg output files
- `processed_data/SPE/`: SpatialExperiment object containing data from all 3 rounds
- `outputs/plots/`: saved plot files

