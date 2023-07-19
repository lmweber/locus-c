#######################################
# LC analyses: number of cells per spot
# Lukas Weber, Jul 2023
#######################################

# module load conda_R/4.2
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "Visium", "02_exploratory_plots")


# ---------
# load data
# ---------

# load saved SPE object from previous script

fn_spe <- here("processed_data", "SPE", "LC_preprocessing.rds")
spe <- readRDS(fn_spe)

dim(spe)

table(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)


# ------------------------
# number of cells per spot
# ------------------------

files_cellsPerSpot <- list.files(here("cells_per_spot"), 
                                 pattern = ".csv", full.names = TRUE)

# re-order samples (6 samples from rounds 1 and 3)
files_cellsPerSpot <- files_cellsPerSpot[c(5:6, 1:4)]
samples_cellsPerSpot <- sample_ids[c(1:2, 6:9)]

# check order
files_cellsPerSpot
samples_cellsPerSpot


# store number of cells per spot in SPE object

colData(spe)$cell_count <- as.numeric(NA)

for (i in seq_along(files_cellsPerSpot)) {
  dat <- read_csv(files_cellsPerSpot[[i]])
  # keep spots over tissue
  dat <- dat %>% filter(tissue == 1)
  key_id <- paste0(samples_cellsPerSpot[i], "_", dat$barcode)
  # store cell counts from selected column in SPE object
  n <- dat %>% select(starts_with("N")) %>% unlist() %>% unname()
  colData(spe)[key_id, "cell_count"] <- n
}


# select spots within annotated LC regions only

spe_sub <- spe[, colData(spe)$annot_region]


# calculate median and IQR per sample
# medians: range 3-7
# 1st quartiles: range 2-4
# 3rd quartiles: range 5-11
for (i in seq_along(samples_cellsPerSpot)) {
  print(samples_cellsPerSpot[i])
  print(summary(colData(spe_sub)[colData(spe_sub)$sample_id == samples_cellsPerSpot[i], "cell_count"]))
}


df <- 
  colData(spe_sub) |> 
  as.data.frame() |> 
  filter(sample_id %in% samples_cellsPerSpot) |> 
  select(c("sample_id", "cell_count"))


# plot boxplots
ggplot(df, aes(x = sample_id, y = cell_count, color = sample_id)) + 
  geom_boxplot() + 
  scale_color_brewer(palette = "Set2") + 
  labs(y = "number of cells") + 
  ggtitle("Cells per spot within annotated LC regions") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 15, vjust = 0.75, hjust = 0.75))

fn <- file.path(dir_plots, "cells_per_spot", "n_cells_per_spot_boxplots")
ggsave(paste0(fn, ".pdf"), width = 7.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 7.5, height = 4.5)


# plot histograms
ggplot(df, aes(x = cell_count)) + 
  facet_wrap(~ sample_id, scales = "free") + 
  geom_histogram(binwidth = 1, color = "black", fill = "navy") + 
  labs(x = "number of cells") + 
  ggtitle("Cells per spot within annotated LC regions") + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "cells_per_spot", "n_cells_per_spot_histograms")
ggsave(paste0(fn, ".pdf"), width = 8, height = 5)
ggsave(paste0(fn, ".png"), width = 8, height = 5)

