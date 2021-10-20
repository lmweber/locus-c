##########################################
# LC project
# Script to build SpatialExperiment object
# Lukas Weber, Sep 2021
##########################################

# Build a single SpatialExperiment object containing all samples


# module load conda_R/4.1.x
# Rscript filename.R


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

