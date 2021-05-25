# Locus coeruleus (LC) 10x Visium project

Repository for locus coeruleus (LC) project.


## Location of files

Main directory containing git repository and data files is located at: `/dcs04/hicks/data/lweber/locus-c/`

(Previously this was at: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/locus-c/` - moved due to space limitations.)


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


## Other miscellaneous files

- `sample_info/`: spreadsheets with sample information
- `Andrew_MiSeq/`: backup of files from initial runs on MiSeq data by Andrew Jaffe

