# Code repository for human locus coeruleus (LC) analyses

This repository contains code scripts to reproduce analyses and figures in our manuscript:

- Weber and Divecha et al. (2022), *The gene expression landscape of the human locus coeruleus revealed by single-nucleus and spatially-resolved transcriptomics.* [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2022.10.28.514241v1).


## Overview

We applied spatially-resolved transcriptomics (10x Genomics Visium) and single-nucleus RNA-sequencing to generate transcriptome-wide, spatially-resolved gene expression data from the locus coeruleus (LC) in neurotypical adult human brain donors.

The dataset is freely accessible in both interactive web-based and downloadable formats.

Our manuscript describes the data, our analyses, and provides links to access the data.


## Contents

This repository contains code scripts to reproduce analyses and figures, as follows:

- [code/analyses/](code/analyses/): scripts for analysis workflow for Visium SRT data
- [code/analyses_snRNAseq/](code/analyses_snRNAseq/): scripts for analysis workflow for snRNA-seq data


Preprocessing scripts:

- [code/spaceranger/](code/spaceranger/): scripts to run Space Ranger for pre-processing Visium SRT data
- [code/cellranger/](code/cellranger/): scripts to run Cell Ranger for pre-processing Visium SRT data


## Additional information

Path to project directory containing files on our compute cluster: `/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/`

Sample information: [sample_info/](sample_info/)

Web summaries from Space Ranger and Cell Ranger: [web_summaries/](web_summaries/)

