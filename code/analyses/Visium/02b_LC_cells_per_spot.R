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

# check 
files_cellsPerSpot
samples_cellsPerSpot

# number of cells per spot in tissue regions
n_cellsPerSpot <- as.list(rep(NA, 6))
names(n_cellsPerSpot) <- samples_cellsPerSpot

for (i in seq_along(n_cellsPerSpot)) {
  dat <- read_csv(files_cellsPerSpot[[i]])
  n <- dat %>% filter(tissue == 1) %>% select(starts_with("N")) %>% unlist() %>% unname()
  n_cellsPerSpot[[i]] <- n
}

# calculate median and IQR per sample
# medians: range 3-5
# 1st quartiles: range 1-3
# 3rd quartiles: range 5-8
lapply(n_cellsPerSpot, summary)

df <- 
  data.frame(
    sample_id = rep(names(n_cellsPerSpot), times = sapply(n_cellsPerSpot, length)), 
    n = unname(unlist(n_cellsPerSpot))
  ) %>% 
  mutate(sample_id = factor(sample_id, levels = samples_cellsPerSpot))


# plot boxplots
ggplot(df, aes(x = sample_id, y = n, color = sample_id)) + 
  geom_boxplot() + 
  scale_y_sqrt() + 
  scale_color_brewer(palette = "Set2") + 
  labs(y = "number of nuclei (sqrt scale)") + 
  ggtitle("Number of nuclei per spot (in tissue)") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 15, vjust = 0.75, hjust = 0.75))

fn <- file.path(dir_plots, "cells_per_spot", "n_cells_per_spot_boxplots")
ggsave(paste0(fn, ".pdf"), width = 7.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 7.5, height = 4.5)


# trim at 98th percentile for plotting histograms
n_trimmed <- n_cellsPerSpot
for (i in seq_along(n_cellsPerSpot)) {
  n <- n_cellsPerSpot[[i]]
  n <- n[n < quantile(n, probs = 0.98)]
  n_trimmed[[i]] <- n
}

df_trimmed <- 
  data.frame(
    sample_id = rep(names(n_trimmed), times = sapply(n_trimmed, length)), 
    n = unname(unlist(n_trimmed))
  ) %>% 
  mutate(sample_id = factor(sample_id, levels = samples_cellsPerSpot))

# plot histograms
ggplot(df_trimmed, aes(x = n)) + 
  facet_wrap(~ sample_id, scales = "free") + 
  geom_histogram(binwidth = 1, color = "black", fill = "navy") + 
  scale_x_continuous(breaks = seq(0, 60, by = 5)) + 
  labs(x = "number of nuclei") + 
  ggtitle("Number of nuclei per spot (in tissue, trimmed)") + 
  theme_bw() + 
  theme(panel.grid = element_blank())

fn <- file.path(dir_plots, "cells_per_spot", "n_cells_per_spot_histograms")
ggsave(paste0(fn, ".pdf"), width = 8, height = 5)
ggsave(paste0(fn, ".png"), width = 8, height = 5)

