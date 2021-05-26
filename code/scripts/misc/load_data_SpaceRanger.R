#############################################
# Code to load data
# old version from Space Ranger documentation
# not converting directly to VisiumExperiment
#############################################


# to run on cluster
# module load conda_R/4.0
# R


# ---------
# load data
# ---------

# Using "filtered" feature-barcode matrix from Space Ranger, which contains only spots under tissue.

# More details: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices

library(Matrix)

# select sample
sample_name <- "LC_1"

# path to Space Ranger output files
dir_outputs <- "~/data/MiSeq_Pilot/outputs/NextSeq"
#dir_outputs <- "../outputs/NextSeq"

# modified code from Space Ranger documentation
matrix_dir <- file.path(dir_outputs, sample_name, "outs", "filtered_feature_bc_matrix")
barcode_path <- file.path(matrix_dir, "barcodes.tsv.gz")
features_path <- file.path(matrix_dir, "features.tsv.gz")
matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix_path)

feature_names <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
barcode_names <- read.delim(barcode_path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- barcode_names$V1
rownames(mat) <- feature_names$V1

colnames(feature_names) <- c("gene_id", "gene_name", "feature_type")
colnames(barcode_names) <- "barcode_id"

head(feature_names)
head(barcode_names)

