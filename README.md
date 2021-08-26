# Locus coeruleus (LC) 10x Visium project

Repository for locus coeruleus (LC) project


## Location of files

Files are located at JHPCE at:

- `/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/`
- `/dcs04/hicks/data/lweber/locus-c/`

Backup of raw data files is also at:

- `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/`

Previously files were at:

-`/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/locus-c/`


## FASTQ files

FASTQ files are in the following subdirectories:

- Round 1: `MiSeq_2020-08-12` and `NextSeq_2020-09-25`
- Round 2: `Linda_2021-05-21`
- Round 3: `KMay_2021-07-09`


## Image files

Image files are in the following directories:

- Rounds 1 and 2: `/dcs04/hicks/data/lweber/locus-c/images/split/`
- Round 3: `/dcs04/lieber/marmaypag/visiumImages_LIBD001/raw-data/2021-06-15_V10B01-003/`


## Code scripts

- `code/scripts/`: Space Ranger and other shell scripts
- `code/analysis/`: R and RMarkdown analysis scripts


## Data files

- `raw_data/fastq/`: FASTQ files
- `raw_data/images/`: TIF image files
- `raw_data/images_align/`: fiducial alignment JSON files from Loupe (required for Space Ranger)
- `processed_data/spaceranger/`: output files from Space Ranger
- `inputs/`: other external input files


## Output files

- `outputs/objects/`: saved output objects from R/RMarkdown scripts
- `outputs/plots/`: saved plots

