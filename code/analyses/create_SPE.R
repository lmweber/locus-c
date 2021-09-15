###########################################
# LC project
# Script to create SpatialExperiment object
# Lukas Weber, Sep 2021
###########################################

# Create a single SpatialExperiment object containing all samples


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

paths <- here("processed_data", "spaceranger", 
              c(rep("NextSeqMiSeq", n_round1), 
                rep("Linda_2021-05-21", n_round2), 
                rep("KMay_2021-07-09", n_round3)))

stopifnot(length(sample_ids) == length(rounds))
stopifnot(length(sample_ids) == length(paths))

df_samples <- data.frame(
  sample_id = sample_ids, 
  round_id = rounds, 
  path = paths
)

df_samples


# ----------------------------------
# load data from spaceranger outputs
# ----------------------------------

spe <- read10xVisium(
  samples = here(df_samples$path, df_samples$sample_id, "outs"), 
  sample_id = df_samples$sample_id, 
  type = "sparse", 
  data = "raw", 
  images = c("hires", "lowres"), 
  load = TRUE
)

