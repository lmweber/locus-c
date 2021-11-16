####################################################################
# LC project
# Script to build SpatialExperiment object (with info for Shiny app)
# Lukas Weber, Nov 2021
####################################################################

# module load conda_R/4.1.x
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)


# ---------------
# set up metadata
# ---------------

n_round1 <- 2
n_round2 <- 3
n_round3 <- 4

sample_ids <- c(
  "Br6522_LC_1_round1", "Br6522_LC_2_round1", 
  "Br8153_LC_round2", "Br5459_LC_round2", "Br2701_LC_round2", 
  "Br6522_LC_round3", "Br8079_LC_round3", "Br2701_LC_round3", "Br8153_LC_round3"
)

rounds <- c(
  rep("round1", n_round1), 
  rep("round2", n_round2), 
  rep("round3", n_round3)
)

paths_spaceranger <- here(
  "processed_data", 
  "spaceranger", 
  c(rep("NextSeqMiSeq", n_round1), 
    rep("Linda_2021-05-21", n_round2), 
    rep("KMay_2021-07-09", n_round3))
)

paths_vistoseg <- here(
  "processed_data", 
  "VistoSeg", 
  c(rep("round1", n_round1), 
    rep("round2", n_round2), 
    rep("round3", n_round3))
)

stopifnot(length(sample_ids) == length(rounds))
stopifnot(length(sample_ids) == length(paths_spaceranger))
stopifnot(length(sample_ids) == length(paths_vistoseg))

df_samples <- data.frame(
  sample_id = sample_ids, 
  round_id = rounds, 
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

# update column names in spatialCoords
colnames(spatialCoords(spe)) <- c("x", "y")


# -----------------------------------------
# add additional sample metadata in colData
# -----------------------------------------

# number of spots per sample (with samples in correct order)
n_spots <- table(colData(spe)$sample_id)[sample_ids]
n_spots

stopifnot(length(rounds) == length(n_spots))

# round IDs for each spot
rep_rounds <- rep(rounds, times = n_spots)

stopifnot(length(rep_rounds) == ncol(spe))

colData(spe)$round_id <- rep_rounds


# -------------------------------
# add VistoSeg outputs in colData
# -------------------------------

# load as list and then stack rows to ensure spots for each sample are in correct order
list_vs <- as.list(rep(NA, length(df_samples$sample_id)))

for (i in seq_along(df_samples$sample_id)) {
  fn <- file.path(df_samples$path_vistoseg[i], df_samples$sample_id[i], "tissue_spot_counts.csv")
  spot_counts <- read.csv(fn)
  barcodes_ord <- rownames(colData(spe[, colData(spe)$sample_id == df_samples$sample_id[i]]))
  spot_counts_ord <- spot_counts[match(barcodes_ord, spot_counts$barcode), ]
  rownames(spot_counts_ord) <- NULL
  list_vs[[i]] <- spot_counts_ord
}

vistoseg_all <- do.call("rbind", list_vs)

stopifnot(ncol(spe) == nrow(vistoseg_all))
stopifnot(all(rownames(colData(spe)) == vistoseg_all$barcode))

# store in SpatialExperiment object
# note there is some redundancy with columns in spatialData but will leave this as a check
colData(spe) <- cbind(colData(spe), vistoseg_all)


# -----------
# save object
# -----------

fn_out <- here("processed_data", "SPE", "LCrounds1to3_SPE_raw.rds")
saveRDS(spe, fn_out)


# -----------------------------
# additional info for Shiny app
# -----------------------------

# using code from Leonardo Collado-Torres available at:
# http://research.libd.org/spatialLIBD/articles/TenX_data_download.html

## Spot info

colData(spe)$key <- paste(colData(spe)$sample_id, colData(spe)$barcode, sep = "_")
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
spe <- spe[, spatialData(spe)$in_tissue]

dim(spe)


## Manual annotations

## Add a variable for saving manual annotations
colData(spe)$ManualAnnotation <- "NA"


# -------------------------
# save object for Shiny app
# -------------------------

fn_out <- here("processed_data", "SPE", "LCrounds1to3_SPE_shiny.rds")
saveRDS(spe, fn_out)

