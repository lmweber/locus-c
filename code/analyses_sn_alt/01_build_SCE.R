##################################################################
# LC snRNA-seq analyses: build SCE object from Cell Ranger outputs
# Lukas Weber, Sep 2022
# using code by Matthew N Tran and Leonardo Collado-Torres
##################################################################

# run in interactive session on JHPCE

# screen -S LC
# qrsh -l mem_free=20G,h_vmem=22G,h_fsize=200G -now n
# cd /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_snRNAseq
# module load conda_R/4.2
# R


library(SingleCellExperiment)
library(DropletUtils)
library(here)


# -----------------
# Create SCE object
# -----------------

# using Cell Ranger "filtered" outputs for nuclei calling


# sample metadata
# table contains the following info: 1c_m = Br6522, 2c_m = Br8079, 3c_m = Br2701
sample_info <- read.table(here("fastq", "snRNA-seq", "sample_libs_info.tsv"))

# Cell Ranger output files
sample_info$path <- file.path(
  here("processed_data", "cellranger"), 
  sample_info$V1, 
  "outs", 
  "filtered_feature_bc_matrix"
)
stopifnot(all(file.exists(sample_info$path)))


# build SCE object
sce <- read10xCounts(
  samples = sample_info$path, 
  sample.names = paste0(sample_info$V3, "_", sample_info$V2), 
  type = "sparse", 
  col.names = TRUE
)


## Read in the gene information from the annotation GTF file
## code from: https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/build_basic_sce.R
gtf <-
  rtracklayer::import(
    "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
  )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[ , c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SCE object
rowRanges(sce) <- gtf[match_genes]

## Inspect object
sce


# number of nuclei per sample
table(sce$Sample)


# -----------
# Save object
# -----------

fn_out <- here("processed_data", "SCE", "sce_raw")
saveRDS(sce, paste0(fn_out, ".rds"))

