# Code repository for locus coeruleus (LC) data analyses

This repository contains code scripts to reproduce analyses and figures in our manuscript:

- Weber and Divecha et al. (2023), *The gene expression landscape of the human locus coeruleus revealed by single-nucleus and spatially-resolved transcriptomics.* [eLife Reviewed Preprint](https://elifesciences.org/reviewed-preprints/84628).


## Overview

We used a combination of spatially-resolved transcriptomics (10x Genomics Visium platform) and single-nucleus RNA-sequencing to characterize the molecular landscape of the LC region and the transcriptomic profile of LC neurons in the human brain.

The dataset is freely accessible in both interactive web-based and downloadable formats.

Our manuscript describes the data, analyses, and provides links to access the data.


## Contents

Scripts to reproduce analyses and figures in manuscript:

- [code/analyses/Visium/](code/analyses/Visium/): scripts for analysis workflow for Visium data
- [code/analyses/snRNAseq/](code/analyses/snRNAseq/): scripts for analysis workflow for snRNA-seq data


Preprocessing scripts:

- [code/preprocessing/spaceranger/](code/preprocessing/spaceranger/): scripts to run Space Ranger for pre-processing Visium data
- [code/preprocessing/cellranger/](code/preprocessing/cellranger/): scripts to run Cell Ranger for pre-processing Visium data


## Additional information

Path to project directory on our compute cluster: `/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/`

Sample information: [sample_info/](sample_info/)

Web summaries from Space Ranger and Cell Ranger: [web_summaries/](web_summaries/)

