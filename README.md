# Locus coeruleus (LC) 10x Visium project

Repository for our code on locus coeruleus project


## Files

- location of our files: `/dcl02/leased/shicks/locus_c/`
- location of repository: `/dcl02/leased/shicks/locus_c/locus-c`


### Scripts

- scripts to run `spaceranger count` using raw images and manually aligned .json files (by Lukas Weber and Stephanie Hicks): [scripts/run_spaceranger/](scripts/run_spaceranger)

- additional scripts saved in other subdirectories in [scripts/](scripts/)


## LIBD files

- main directory: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/`


### Sample information

- summary spreadsheet: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Visium LC pilot_072120 Master.xlsx`

- contains sample ID, sample name, slide serial number, capture area ID


### FASTQ files

- MiSeq FASTQ files (4 samples: DLPFC, LC_1, LC_2, LC): `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ/`

- NextSeq FASTQ files (2 samples: LC_1, LC_2): `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/FASTQ_NextSeq/`


### Image files

- raw image files (`.tif` files, ~3 GB per sample): `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw/`

- manual alignment .json files from Loupe for raw image files (run by Lukas Weber): `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/images_raw_align_json/`

- smaller screenshots of image files (`.jpeg`) and manual alignment .json files from Loupe (run by Andrew Jaffe): `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Snapshots_Plus_Loupe/`


### Reference transcriptome

- reference transcriptome downloaded from 10x website: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A/`


### Scripts

- initial script to run `spaceranger count` using screenshots and manually aligned .json files (by Andrew Jaffe): `dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/pilot_align.sh`


## Authors

- Lukas Weber
- Stephanie Hicks

