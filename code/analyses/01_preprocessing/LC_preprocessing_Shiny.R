####################################################
# LC project
# Preprocessing (with additional info for Shiny app)
# Lukas Weber, Mar 2022
####################################################

# note: using SpatialExperiment version 1.5.4 (Bioconductor version 3.15)

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(scater)
library(scran)


# -------------------
# experiment metadata
# -------------------

# sample IDs, additional sample info, paths to input files, other metadata

# sample IDs
sample_ids <- c(
  "Br6522_LC_1_round1", "Br6522_LC_2_round1", 
  "Br8153_LC_round2", "Br5459_LC_round2", "Br2701_LC_round2", 
  "Br6522_LC_round3", "Br8079_LC_round3", "Br2701_LC_round3", "Br8153_LC_round3"
)

# number of samples per round
n_round1 <- 2
n_round2 <- 3
n_round3 <- 4

rounds <- c(
  rep("round1", n_round1), 
  rep("round2", n_round2), 
  rep("round3", n_round3)
)

# paths to Space Ranger input files
paths_spaceranger <- here("processed_data", "spaceranger", 
  c(rep("NextSeqMiSeq", n_round1), 
    rep("Linda_2021-05-21", n_round2), 
    rep("KMay_2021-07-09", n_round3))
)

# paths to VistoSeg input files (number of cells per spot)
paths_vistoseg <- here("inputs", "VistoSeg", 
  c(rep("round1", n_round1), 
    rep("round2", n_round2), 
    rep("round3", n_round3))
)

# paths and filenames for manual annotation files
# note: multiple files per sample for samples containing multiple parts
files_annot <- list(
  Br6522_LC_1_round1 = here("inputs", "annotations", "Br6522_LC_1_round1", "Br6522_LC_1_round1_lasso_spots.csv"), 
  Br6522_LC_2_round1 = here("inputs", "annotations", "Br6522_LC_2_round1", "Br6522_LC_2_round1_lasso_spots.csv"), 
  Br8153_LC_round2 = c(here("inputs", "annotations", "Br8153_LC_round2", "Br8153_LC_round2_left_lasso_spots.csv"), 
                       here("inputs", "annotations", "Br8153_LC_round2", "Br8153_LC_round2_right_lasso_spots.csv")), 
  Br5459_LC_round2 = c(here("inputs", "annotations", "Br5459_LC_round2", "Br5459_LC_round2_left_lasso_spots.csv"), 
                       here("inputs", "annotations", "Br5459_LC_round2", "Br5459_LC_round2_right_lasso_spots.csv")), 
  Br2701_LC_round2 = c(here("inputs", "annotations", "Br2701_LC_round2", "Br2701_LC_round2_top_lasso_spots.csv"), 
                       here("inputs", "annotations", "Br2701_LC_round2", "Br2701_LC_round2_bottom_lasso_spots.csv")), 
  Br6522_LC_round3 = c(here("inputs", "annotations", "Br6522_LC_round3", "Br6522_LC_round3_left_lasso_spots.csv"), 
                       here("inputs", "annotations", "Br6522_LC_round3", "Br6522_LC_round3_right_lasso_spots.csv")), 
  Br8079_LC_round3 = c(here("inputs", "annotations", "Br8079_LC_round3", "Br8079_LC_round3_left_lasso_spots.csv"), 
                       here("inputs", "annotations", "Br8079_LC_round3", "Br8079_LC_round3_right_lasso_spots.csv")), 
  Br2701_LC_round3 = c(here("inputs", "annotations", "Br2701_LC_round3", "Br2701_LC_round3_left_lasso_spots.csv"), 
                       here("inputs", "annotations", "Br2701_LC_round3", "Br2701_LC_round3_right_lasso_spots.csv")), 
  Br8153_LC_round3 = c(here("inputs", "annotations", "Br8153_LC_round3", "Br8153_LC_round3_left_lasso_spots.csv"))
)

# part IDs per sample
part_ids <- list(
  Br6522_LC_1_round1 = "single", 
  Br6522_LC_2_round1 = "single", 
  Br8153_LC_round2 = c("left", "right"), 
  Br5459_LC_round2 = c("left", "right"), 
  Br2701_LC_round2 = c("top", "bottom"), 
  Br6522_LC_round3 = c("left", "right"), 
  Br8079_LC_round3 = c("left", "right"), 
  Br2701_LC_round3 = c("left", "right"), 
  Br8153_LC_round3 = "left"
)

# number of parts per sample
n_parts <- unname(sapply(files_annot, length))


stopifnot(length(sample_ids) == length(rounds))
stopifnot(length(sample_ids) == length(paths_spaceranger))
stopifnot(length(sample_ids) == length(paths_vistoseg))
stopifnot(length(sample_ids) == length(files_annot))
stopifnot(length(sample_ids) == length(part_ids))
stopifnot(length(sample_ids) == length(n_parts))
stopifnot(all(sample_ids == names(files_annot)))
stopifnot(all(sample_ids == names(part_ids)))


# combined data frame (except lists)
df_samples <- data.frame(
  sample_id = sample_ids, 
  round_id = rounds, 
  n_parts = n_parts, 
  path_spaceranger = paths_spaceranger, 
  path_vistoseg = paths_vistoseg
)

df_samples


# ----------------------------------
# load data from spaceranger outputs
# ----------------------------------

spe <- read10xVisium(
  samples = here(df_samples$path_spaceranger, df_samples$sample_id, "outs"), 
  sample_id = df_samples$sample_id, 
  type = "sparse", 
  data = "raw", 
  images = c("hires", "lowres"), 
  load = TRUE
)

spe
dim(spe)
head(spatialCoords(spe))


# -------------------------------------
# add additional sample info in colData
# -------------------------------------

# round IDs

# get number of spots per sample (with samples in correct order)
n_spots <- table(colData(spe)$sample_id)[sample_ids]
n_spots
stopifnot(length(rounds) == length(n_spots))

# repeat round IDs for each spot
rep_rounds <- rep(rounds, times = n_spots)
stopifnot(length(rep_rounds) == ncol(spe))

colData(spe)$round_id <- rep_rounds


# key IDs (unique combination of sample IDs and barcode IDs)

key_ids <- paste(colData(spe)$sample_id, rownames(colData(spe)), sep = "_")

colData(spe)$key_id <- key_ids
rownames(colData(spe)) <- key_ids


# -------------------------------
# add VistoSeg outputs in colData
# -------------------------------

# load as list and then stack rows to ensure spots for each sample are in correct order
list_vs <- as.list(rep(NA, length(df_samples$sample_id)))

for (i in seq_along(df_samples$sample_id)) {
  fn <- file.path(df_samples$path_vistoseg[i], df_samples$sample_id[i], "tissue_spot_counts.csv")
  spot_counts <- read.csv(fn)
  barcodes_ord <- gsub("^.*_", "", rownames(colData(spe[, colData(spe)$sample_id == df_samples$sample_id[i]])))
  spot_counts_ord <- spot_counts[match(barcodes_ord, spot_counts$barcode), ]
  rownames(spot_counts_ord) <- NULL
  list_vs[[i]] <- spot_counts_ord
}

vistoseg_all <- do.call("rbind", list_vs)

stopifnot(ncol(spe) == nrow(vistoseg_all))


# check and ignore duplicated columns
all(gsub("^.*_", "", rownames(colData(spe))) == vistoseg_all$barcode)
all(colData(spe)$in_tissue == as.logical(vistoseg_all$tissue))
all(colData(spe)$array_row == vistoseg_all$row)
all(colData(spe)$array_col == vistoseg_all$col)

# some spatialCoords are off by one pixel, likely due to rounding
# since the difference is a maximum of one pixel this is okay
all(spatialCoords(spe)[, "pxl_row_in_fullres"] == vistoseg_all$imagerow)
table(spatialCoords(spe)[, "pxl_row_in_fullres"] == vistoseg_all$imagerow)
max(abs(spatialCoords(spe)[, "pxl_row_in_fullres"] - vistoseg_all$imagerow))

all(spatialCoords(spe)[, "pxl_col_in_fullres"] == vistoseg_all$imagecol)
table(spatialCoords(spe)[, "pxl_col_in_fullres"] == vistoseg_all$imagecol)
max(abs(spatialCoords(spe)[, "pxl_col_in_fullres"] - vistoseg_all$imagecol))


# store in SPE object
colData(spe)$cell_count <- vistoseg_all$count


# ---------------------------------
# add manual annotations in colData
# ---------------------------------

# load .csv files containing manual annotations and store as follows:
# - one column with region-level annotations
# - one column with spot-level annotations
# - one column identifying parts of sample

colData(spe)$annot_region <- as.logical(NA)
colData(spe)$part_id <- as.character(NA)
colData(spe)$annot_spot <- as.logical(NA)

# loop over samples
for (i in seq_along(files_annot)) {
  fns <- files_annot[[i]]
  # loop over parts per sample
  for (k in seq_along(fns)) {
    df <- read.csv(fns[k])
    # check that annotation names are as expected, i.e. sample_name_lasso and 
    # sample_name_lasso_spots (check this manually for each sample)
    print(table(df$ManualAnnotation))
    # drop NAs
    df <- df[!is.na(df$ManualAnnotation), ]
    print(dim(df))
    # store key IDs in rownames
    rownames(df) <- paste(df$sample_id, df$spot_name, sep = "_")
    # get key IDs and store annotations in colData columns
    # combined region(s)
    keys_region <- rownames(df)[grepl("lasso", df$ManualAnnotation)]
    print(length(keys_region))
    colData(spe)[keys_region, "annot_region"] <- TRUE
    # part IDs
    colData(spe)[keys_region, "part_id"] <- part_ids[[i]][k]
    # individual spots
    keys_spot <- rownames(df)[grepl("spot", df$ManualAnnotation)]
    print(length(keys_spot))
    colData(spe)[keys_spot, "annot_spot"] <- TRUE
  }
}

# replace NA values with FALSE
colData(spe)$annot_region[is.na(colData(spe)$annot_region)] <- FALSE
colData(spe)$annot_spot[is.na(colData(spe)$annot_spot)] <- FALSE

# combined sample IDs and part IDs
is_nas <- is.na(colData(spe)$part_id)
colData(spe)$sample_part_id <- as.character(NA)
colData(spe)[!is_nas, "sample_part_id"] <- 
  paste(colData(spe)$sample_id, colData(spe)$part_id, sep = "_")[!is_nas]


# -----------------------------
# additional info for Shiny app
# -----------------------------

# using code from Leonardo Collado-Torres available at:
# http://research.libd.org/spatialLIBD/articles/TenX_data_download.html

## Spot info

colData(spe)$key <- paste(colData(spe)$sample_id, rownames(colData(spe)), sep = "_")
colData(spe)$sum_umi <- colSums(counts(spe))
colData(spe)$sum_gene <- colSums(counts(spe) > 0)


## Gene info

## Read in the gene information from the annotation GTF file provided by 10x Genomics
gtf <- rtracklayer::import(
  here("data", "refdata-gex-GRCh38-2020-A/genes/genes.gtf")
)
## Subject to genes only
gtf <- gtf[gtf$type == "gene"]
## Set the names to be the gene IDs
names(gtf) <- gtf$gene_id
## Match the genes
match_genes <- match(rownames(spe), gtf$gene_id)
## They should all be present if you are using the correct GTF file from 10x Genomics
stopifnot(all(!is.na(match_genes)))
## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", 
                             "gene_name", "gene_type")]
## Add the gene info to our SPE object
rowRanges(spe) <- gtf[match_genes]
## Inspect the gene annotation data we added
rowRanges(spe)


## Additional info for spatialLIBD

## Add information used by spatialLIBD
rowData(spe)$gene_search <- paste(
  rowData(spe)$gene_name, rowData(spe)$gene_id, sep = "; "
)
## Compute chrM expression and chrM expression ratio
is_mito <- which(seqnames(spe) == "chrM")
colData(spe)$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
colData(spe)$expr_chrM_ratio <- colData(spe)$expr_chrM / colData(spe)$sum_umi


## Filter zeros (genes and spots)

dim(spe)

## Remove genes with no counts
no_expr <- which(rowSums(counts(spe)) == 0)
## Number of genes with no counts
length(no_expr)
## Proportion of genes with no counts
length(no_expr) / nrow(spe)
## Filter object
spe <- spe[-no_expr, , drop = FALSE]

## Remove spots with no counts
spots_no_counts <- which(colData(spe)$sum_umi == 0)
## Number of spots with no counts
length(spots_no_counts)
## Proportion of spots with no counts
length(spots_no_counts) / ncol(spe)
## Filter object
spe <- spe[, -spots_no_counts, drop = FALSE]

dim(spe)


## Spots over tissue
## Keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue]
dim(spe)


## Default cluster labels
## Add a column of default cluster labels
colData(spe)$all <- "all"


## Manual annotations
## Add a variable for saving manual annotations
colData(spe)$ManualAnnotation <- "NA"


## Genes of interest
## UMIs per spot for genes of interest
colData(spe)$TH <- counts(spe)[which(rowData(spe)$gene_name == "TH"), ]
colData(spe)$DBH <- counts(spe)[which(rowData(spe)$gene_name == "DBH"), ]
colData(spe)$SLC6A2 <- counts(spe)[which(rowData(spe)$gene_name == "SLC6A2"), ]
colData(spe)$SLC6A4 <- counts(spe)[which(rowData(spe)$gene_name == "SLC6A4"), ]
colData(spe)$SLC18A2 <- counts(spe)[which(rowData(spe)$gene_name == "SLC18A2"), ]
colData(spe)$DDC <- counts(spe)[which(rowData(spe)$gene_name == "DDC"), ]
colData(spe)$GCH1 <- counts(spe)[which(rowData(spe)$gene_name == "GCH1"), ]
colData(spe)$MAOA <- counts(spe)[which(rowData(spe)$gene_name == "MAOA"), ]


# --------------------
# quality control (QC)
# --------------------

# using scater package

# note: keep all spots for Shiny app
# QC for downstream analyses is performed in next script


# ---------------------------
# normalization and logcounts
# ---------------------------

# using scran package

# quick clustering for pool-based size factors
# with blocks by sample
set.seed(123)
qclus <- quickCluster(spe, block = colData(spe)$sample_id)
table(qclus)

table(colData(spe)$sample_id, qclus)

# calculate size factors
spe <- computeSumFactors(spe, cluster = qclus)
summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 40)

# note: remove small number of spots with size factors == 0
table(sizeFactors(spe) > 0)
dim(spe)
spe <- spe[, sizeFactors(spe) > 0]
dim(spe)

# calculate logcounts
spe <- logNormCounts(spe)
assayNames(spe)


## Genes of interest
## logcounts per spot for genes of interest
colData(spe)$TH_logcounts <- logcounts(spe)[which(rowData(spe)$gene_name == "TH"), ]
colData(spe)$DBH_logcounts <- logcounts(spe)[which(rowData(spe)$gene_name == "DBH"), ]
colData(spe)$SLC6A2_logcounts <- logcounts(spe)[which(rowData(spe)$gene_name == "SLC6A2"), ]
colData(spe)$SLC6A4_logcounts <- logcounts(spe)[which(rowData(spe)$gene_name == "SLC6A4"), ]
colData(spe)$SLC18A2_logcounts <- logcounts(spe)[which(rowData(spe)$gene_name == "SLC18A2"), ]
colData(spe)$DDC_logcounts <- logcounts(spe)[which(rowData(spe)$gene_name == "DDC"), ]
colData(spe)$GCH1_logcounts <- logcounts(spe)[which(rowData(spe)$gene_name == "GCH1"), ]
colData(spe)$MAOA_logcounts <- logcounts(spe)[which(rowData(spe)$gene_name == "MAOA"), ]


# -------------------------
# save object for Shiny app
# -------------------------

# save as .rds and .RData
fn_out <- here("processed_data", "SPE", "LC_Shiny")
saveRDS(spe, paste0(fn_out, ".rds"))
save(spe, file = paste0(fn_out, ".RData"))

