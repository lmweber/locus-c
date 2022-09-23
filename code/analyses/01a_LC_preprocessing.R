###################################
# LC Visium analyses: preprocessing
# Lukas Weber, Sep 2022
###################################

# note: using R 4.2, Bioconductor 3.15, SpatialExperiment 1.6.1

# module load conda_R/4.2
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(here)
library(SpatialExperiment)


# -------------------
# experiment metadata
# -------------------

# sample IDs, other metadata, paths to input files

# sample IDs
sample_ids <- c(
  "Br6522_LC_1_round1", "Br6522_LC_2_round1", 
  "Br8153_LC_round2", "Br5459_LC_round2", "Br2701_LC_round2", 
  "Br6522_LC_round3", "Br8079_LC_round3", "Br2701_LC_round3", "Br8153_LC_round3"
)

# donor IDs
donor_ids <- gsub("_.*$", "", sample_ids)

# round IDs
n_round1 <- 2
n_round2 <- 3
n_round3 <- 4

round_ids <- c(
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


stopifnot(length(sample_ids) == length(donor_ids))
stopifnot(length(sample_ids) == length(round_ids))
stopifnot(length(sample_ids) == length(paths_spaceranger))
stopifnot(length(sample_ids) == length(paths_vistoseg))
stopifnot(length(sample_ids) == length(files_annot))
stopifnot(all(sample_ids == names(files_annot)))


# combined data frame of sample-level information
df_samples <- data.frame(
  sample_id = sample_ids, 
  donor_id = donor_ids, 
  round_id = round_ids, 
  path_spaceranger = paths_spaceranger, 
  path_vistoseg = paths_vistoseg
)

df_samples


# -----------------------------------
# load data from Space Ranger outputs
# -----------------------------------

spe <- read10xVisium(
  samples = here(df_samples$path_spaceranger, df_samples$sample_id, "outs"), 
  sample_id = df_samples$sample_id, 
  type = "sparse", 
  data = "filtered", 
  images = c("hires", "lowres"), 
  load = TRUE
)

spe
dim(spe)
table(colData(spe)$sample_id)
head(spatialCoords(spe))


# -------------------------------------
# add additional information in colData
# -------------------------------------

# convert sample IDs to factor
colData(spe)$sample_id <- factor(colData(spe)$sample_id, levels = sample_ids)


# get number of spots per sample (with samples in correct order)
n_spots <- table(colData(spe)$sample_id)[sample_ids]
n_spots
stopifnot(length(round_ids) == length(n_spots))


# donor IDs for each spot
rep_donor_ids <- rep(donor_ids, times = n_spots)
stopifnot(length(rep_donor_ids) == ncol(spe))

colData(spe)$donor_id <- factor(rep_donor_ids, levels = unique(donor_ids))


# round IDs for each spot
rep_round_ids <- rep(round_ids, times = n_spots)
stopifnot(length(rep_round_ids) == ncol(spe))

colData(spe)$round_id <- factor(rep_round_ids, levels = unique(round_ids))


# key IDs (unique combination of sample IDs and barcode IDs)
key_ids <- paste(colData(spe)$sample_id, rownames(colData(spe)), sep = "_")

# store in column names, colData, and row names of colData
stopifnot(length(key_ids) == ncol(spe))
colnames(spe) <- key_ids
colData(spe)$key_id <- key_ids

stopifnot(all(rownames(colData(spe)) == key_ids))


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

colData(spe)$annot_region <- as.logical(NA)
colData(spe)$annot_spot <- as.logical(NA)

# loop over samples
for (i in seq_along(files_annot)) {
  fns <- files_annot[[i]]
  # loop over files per sample
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
    # individual spots
    keys_spot <- rownames(df)[grepl("spot", df$ManualAnnotation)]
    print(length(keys_spot))
    colData(spe)[keys_spot, "annot_spot"] <- TRUE
  }
}

# replace NA values for annotations with FALSE
colData(spe)$annot_region[is.na(colData(spe)$annot_region)] <- FALSE
colData(spe)$annot_spot[is.na(colData(spe)$annot_spot)] <- FALSE


# ---------------
# create part IDs
# ---------------

# create part IDs to identify tissue parts per sample

# manually segmenting areas based on spatial coordinates in interactive plots

part_ids_all <- rep(NA, ncol(spe))

for (s in seq_along(sample_ids)) {
  ix <- colData(spe)$sample_id == sample_ids[s]
  spe_sub <- spe[, ix]
  
  if (sample_ids[s] == "Br6522_LC_1_round1") {
    # 1 part
    part_ids_all[ix] <- "single"
  } else if (sample_ids[s] == "Br6522_LC_2_round1") {
    # 1 part
    part_ids_all[ix] <- "single"
  } else if (sample_ids[s] == "Br8153_LC_round2") {
    # 2 parts
    cnd <- with(colData(spe_sub), array_row < 50)
    part_ids_all[ix] <- ifelse(cnd, "left", "right")
  } else if (sample_ids[s] == "Br5459_LC_round2") {
    # 2 parts
    cnd <- with(colData(spe_sub), (array_row < 26 & array_col > 20) | 
                                  (array_row < 28 & array_col <= 20))
    part_ids_all[ix] <- ifelse(cnd, "left", "right")
  } else if (sample_ids[s] == "Br2701_LC_round2") {
    # 2 parts
    cnd <- with(colData(spe_sub), (array_row < 30 & array_col > 88) | 
                                  (array_row >= 30 & array_col > 82))
    part_ids_all[ix] <- ifelse(cnd, "top", "bottom")
  } else if (sample_ids[s] == "Br6522_LC_round3") {
    # 3 parts
    cnd1 <- with(colData(spe_sub), array_row > 50)
    cnd2 <- with(colData(spe_sub), array_col > 90)
    part_ids_all[ix] <- ifelse(cnd1, "right", ifelse(cnd2, "lefttop", "leftbottom"))
  } else if (sample_ids[s] == "Br8079_LC_round3") {
    # 2 parts
    cnd <- with(colData(spe_sub), array_row < 41)
    part_ids_all[ix] <- ifelse(cnd, "left", "right")
  } else if (sample_ids[s] == "Br2701_LC_round3") {
    # 2 parts
    cnd <- with(colData(spe_sub), array_row < 40)
    part_ids_all[ix] <- ifelse(cnd, "left", "right")
  } else if (sample_ids[s] == "Br8153_LC_round3") {
    # 2 parts
    cnd <- with(colData(spe_sub), array_row < 45)
    part_ids_all[ix] <- ifelse(cnd, "left", "right")
  }
}

stopifnot(length(part_ids_all) == ncol(spe))

# store part IDs
colData(spe)$part_id <- part_ids_all

# store combined sample IDs and part IDs
colData(spe)$sample_part_id <- paste(colData(spe)$sample_id, colData(spe)$part_id, sep = "_")


# ---------------------------
# additional gene information
# ---------------------------

# using code from Leonardo Collado-Torres available at:
# http://research.libd.org/spatialLIBD/articles/TenX_data_download.html


## Gene info

## Read in the gene information from the annotation GTF file provided by 10x Genomics
gtf <- rtracklayer::import(
  here("data", "refdata-gex-GRCh38-2020-A/genes/genes.gtf")
)
## Genes only
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
## Add the gene info to SPE object
rowRanges(spe) <- gtf[match_genes]
## Inspect the gene annotation data
rowRanges(spe)


# ---------------
# save SPE object
# ---------------

fn_out <- here("processed_data", "SPE", "LC_preprocessing")
saveRDS(spe, paste0(fn_out, ".rds"))

