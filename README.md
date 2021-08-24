# Locus coeruleus (LC) 10x Visium project

Repository for locus coeruleus (LC) project


## Location of files

Two main copies of files are on JHPCE at:

- `/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/`
- `/dcs04/hicks/data/lweber/locus-c/`

Previously files were at:

-`/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/locus-c/`


## FASTQ files

FASTQ files are in the following directories:

- Round 1: `MiSeq_2020-08-12` and `NextSeq_2020-09-25`

- Round 2: `Linda_2021-05-21`

- Round 3: `KMay_2021-07-09`


## Code files

- `code/scripts/`: shell scripts to run Space Ranger and other shell scripts
- `code/analysis/`: R and RMarkdown scripts for analyses


## Data files

- `raw_data/fastq/`: FASTQ files for each sample
- `raw_data/images/`: TIF images for each sample
- `raw_data/images_align/`: manually aligned JSON files from Loupe for raw image files (for Space Ranger)
- `processed_data/spaceranger/`: output files from Space Ranger
- `inputs/`: other input files


## Output files

- `outputs/objects/`: intermediate output objects from R/RMarkdown scripts
- `outputs/plots/`: plot files from R/RMarkdown scripts

