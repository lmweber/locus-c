# Code repository for human locus coeruleus (LC) analyses

This repository contains code scripts to reproduce analyses and figures in our manuscript:

- Weber and Divecha et al. (2023), *The gene expression landscape of the human locus coeruleus revealed by single-nucleus and spatially-resolved transcriptomics.* [eLife Reviewed Preprint](https://elifesciences.org/reviewed-preprints/84628).


## Overview

We applied spatially-resolved transcriptomics (10x Genomics Visium) and single-nucleus RNA-sequencing to generate transcriptome-scale, spatially-resolved gene expression data from the locus coeruleus (LC) in neurotypical adult human brain donors.

The dataset is freely accessible in both interactive web-based and downloadable formats.

Our manuscript describes the data, our analyses, and provides links to access the data.


## Contents

Scripts to reproduce analyses and figures:

- [code/analyses/Visium/](code/analyses/Visium/): scripts for analysis workflow for Visium SRT data
- [code/analyses/snRNAseq/](code/analyses/snRNAseq/): scripts for analysis workflow for snRNA-seq data


Preprocessing scripts:

- [code/preprocessing/spaceranger/](code/preprocessing/spaceranger/): scripts to run Space Ranger for pre-processing Visium SRT data
- [code/preprocessing/cellranger/](code/preprocessing/cellranger/): scripts to run Cell Ranger for pre-processing Visium SRT data


## Additional information

Path to project directory on our compute cluster: `/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/`

Sample information: [sample_info/](sample_info/)

Web summaries from Space Ranger and Cell Ranger: [web_summaries/](web_summaries/)

