# Locus coeruleus (LC) 10x Visium project

Repository for our code on locus coeruleus project


## Files

- location of our files: `/dcl02/leased/shicks/locus_c`
- location of Space Ranger reference: `/dcl02/leased/shicks/spaceranger`


## LIBD files

### Sample information 

Summary spreadsheet: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Visium LC pilot_072120 Master.xlsx`

- contains sample ID, sample name, slide serial number, capture area ID


### Image files 

location of LIBD files: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot`

- Raw image files (`.tif`): `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Raw`
- Screen shots of image files (`.jpeg`) and Loupe manual alignment json files (run by Andrew Jaffe): `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/Snapshots_Plus_Loupe`
- Manually aligned files (`.json`) by Lukas Weber: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/Images/json_after_manual_loupe`


### Scripts 

- Initial script by Andrew with code to run `spaceranger count`: `dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/MiSeq_Pilot/pilot_align.sh`
- Updated script by Lukas and Stephanie to run `spaceranger count` with manually aligned json files from images (this directory): `scripts/spaceranger.sh`


## Authors

- Lukas Weber
- Stephanie Hicks