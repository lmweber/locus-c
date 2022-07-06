############################
# LC analyses: nnSVG results
# Lukas Weber, Jun 2022
############################

# module load conda_R/devel
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "08_nnSVG")

# directory to save outputs
dir_outputs <- here("outputs", "08_nnSVG")


# ---------
# load data
# ---------

# load results from previous script

fn_out <- here("outputs", "LC_nnSVG_results.rsd")
res_list <- readRDS(fn_out)

names(res_list)


# ---------------
# combine results
# ---------------

# sum gene ranks across sample-parts to generate overall ranking

# number of genes that passed filtering for each sample-part
sapply(res_list, nrow)

# note: exclude Br2701_LC_round2_top since too few genes passed filtering for this sample-part
ix_exc <- which(names(res_list) == "Br2701_LC_round2_top")
res_list <- res_list[-ix_exc]
sample_part_ids <- sample_part_ids[-ix_exc]


# match results from each sample-part and store in correct rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_part_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_part_ids

for (s in seq_along(sample_part_ids)) {
  stopifnot(colnames(res_ranks)[s] == sample_part_ids[s])
  stopifnot(colnames(res_ranks)[s] == names(res_list)[s])
  
  rownames_s <- rownames(res_list[[s]])
  res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}

# keep only genes that were not filtered out in all sample-parts
res_ranks <- na.omit(res_ranks)

# calculate average ranks
avg_ranks <- sort(rowMeans(res_ranks))


# summary table
df_summary <- data.frame(
  gene_id = names(avg_ranks), 
  gene_name = rowData(spe)[names(avg_ranks), "gene_name"], 
  gene_type = rowData(spe)[names(avg_ranks), "gene_type"], 
  avg_rank = unname(avg_ranks), 
  row.names = names(avg_ranks)
)


# -------------
# plot top SVGs
# -------------

ix_gene <- which(rowData(spe)$gene_name == "DBH")

df <- as.data.frame(cbind(
  colData(spe), 
  spatialCoords(spe), 
  gene = counts(spe)[ix_gene, ]
))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = gene)) + 
  facet_wrap(~ sample_id, nrow = 2, scales = "free") + 
  geom_point(size = 0.1) + 
  scale_color_gradient(low = "gray80", high = "red", trans = "sqrt", 
                       name = "counts", breaks = range(df$gene)) + 
  scale_y_reverse() + 
  ggtitle("SVG expression") + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        plot.title = element_text(face = "italic"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

