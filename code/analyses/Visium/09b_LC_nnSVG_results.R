############################
# LC analyses: nnSVG results
# Lukas Weber, Oct 2022
############################

# module load conda_R/4.2
# Rscript filename.R

# file location:
# /dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/


library(here)
library(SpatialExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)


# directory to save plots
dir_plots <- here("plots", "Visium", "09_nnSVG")

# directory to save outputs
dir_outputs <- here("outputs", "Visium", "09_nnSVG")


# ---------
# load data
# ---------

# load saved SPE object from previous scripts

fn_spe <- here("processed_data", "SPE", "LC_filtered.rds")
spe <- readRDS(fn_spe)

dim(spe)


# load saved nnSVG results from previous scripts

fn_out <- here(dir_outputs, "LC_nnSVG_results.rds")
res_list <- readRDS(fn_out)

names(res_list)


# sample-part IDs for parts that contain LC annotated regions
sample_part_ids <- c(
  "Br6522_LC_1_round1_single", 
  "Br6522_LC_2_round1_single", 
  "Br8153_LC_round2_left", "Br8153_LC_round2_right", 
  "Br2701_LC_round2_top", "Br2701_LC_round2_bottom", 
  "Br6522_LC_round3_leftbottom", "Br6522_LC_round3_right", 
  "Br8079_LC_round3_left", "Br8079_LC_round3_right", 
  "Br2701_LC_round3_left", "Br2701_LC_round3_right", 
  "Br8153_LC_round3_left"
)

stopifnot(all(names(res_list) == sample_part_ids))


# ------------------------------------
# combine results for multiple samples
# ------------------------------------

# combine results for multiple sample-part IDs by summing or binning gene ranks 
# to generate an overall ranking

# number of genes that passed filtering for each sample-part
sapply(res_list, nrow)


# match results from each sample-part and store in matching rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_part_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_part_ids

for (s in seq_along(sample_part_ids)) {
  stopifnot(colnames(res_ranks)[s] == sample_part_ids[s])
  stopifnot(colnames(res_ranks)[s] == names(res_list)[s])
  
  rownames_s <- rownames(res_list[[s]])
  res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}


# remove genes that were filtered out in all sample-parts
ix_allna <- apply(res_ranks, 1, function(r) all(is.na(r)))
res_ranks <- res_ranks[!ix_allna, ]

dim(res_ranks)


# calculate average ranks
# note missing values due to filtering for sample-parts
avg_ranks <- rowMeans(res_ranks, na.rm = TRUE)


# calculate number of sample-parts where each gene is within top 100 ranked SVGs
# for that sample-part
n_withinTop100 <- apply(res_ranks, 1, function(r) sum(r <= 100, na.rm = TRUE))


# summary table
df_summary <- data.frame(
  gene_id = names(avg_ranks), 
  gene_name = rowData(spe)[names(avg_ranks), "gene_name"], 
  gene_type = rowData(spe)[names(avg_ranks), "gene_type"], 
  overall_rank = rank(avg_ranks), 
  average_rank = unname(avg_ranks), 
  n_withinTop100 = unname(n_withinTop100), 
  row.names = names(avg_ranks)
)

# sort by average rank
df_summary <- df_summary[order(df_summary$average_rank), ]
head(df_summary)

# top n genes
top50genes <- df_summary$gene_name[1:50]
top100genes <- df_summary$gene_name[1:100]


# summary table of "replicated" SVGs (i.e. genes that are highly ranked in at
# least x sample-parts)
df_summaryReplicated <- df_summary[df_summary$n_withinTop100 >= 10, ]

# re-calculate rank within this set
df_summaryReplicated$overall_rank <- rank(df_summaryReplicated$average_rank)

dim(df_summaryReplicated)
head(df_summaryReplicated)

# top "replicated" SVGs
topSVGsReplicated <- df_summaryReplicated$gene_name


# combined table for supplementary table
df_summaryCombined <- df_summary

df_summaryCombined[rownames(df_summaryReplicated), "replicated_overall_rank"] <- 
  df_summaryReplicated$overall_rank

# re-order columns
df_summaryCombined <- df_summaryCombined[, c(1:3, 7, 4:6)]

dim(df_summaryCombined)
head(df_summaryCombined)


# ------------
# plot top SVG
# ------------

head(df_summary, 1)

ix_gene <- which(rowData(spe)$gene_name == "CAPS")

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
  ggtitle(paste0(rowData(spe)$gene_name[ix_gene], " expression")) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(), 
        plot.title = element_text(face = "italic"), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "topSVG_CAPS_expression")
ggsave(paste0(fn, ".pdf"), width = 7.5, height = 4)
ggsave(paste0(fn, ".png"), width = 7.5, height = 4)


# -----------------
# save spreadsheets
# -----------------

# save spreadsheet containing list of top SVGs

res_ranks_ord <- res_ranks[rownames(df_summary), ]

stopifnot(all(rownames(df_summary) == rownames(res_ranks_ord)))
stopifnot(nrow(df_summary) == nrow(res_ranks_ord))

df_all <- cbind(df_summary, res_ranks_ord)
rownames(df_all) <- NULL

# save .csv file
fn_all <- file.path(dir_outputs, "topSVGs_nnSVG_avgRanks.csv")
write.csv(df_all, file = fn_all, row.names = FALSE)


# save spreadsheet containing list of top "replicated" SVGs

res_ranksReplicated_ord <- res_ranks[rownames(df_summaryReplicated), ]

stopifnot(all(rownames(df_summaryReplicated) == rownames(res_ranksReplicated_ord)))
stopifnot(nrow(df_summaryReplicated) == nrow(res_ranksReplicated_ord))

df_replicated <- cbind(df_summaryReplicated, res_ranksReplicated_ord)
rownames(df_replicated) <- NULL

# save .csv file
fn_replicated <- file.path(dir_outputs, "topSVGs_nnSVG_avgRanks_replicated.csv")
write.csv(df_replicated, file = fn_replicated, row.names = FALSE)


# save combined spreadsheet for supplementary table

res_ranksCombined_ord <- res_ranks[rownames(df_summaryCombined), ]

stopifnot(all(rownames(df_summaryCombined) == rownames(res_ranksCombined_ord)))
stopifnot(nrow(df_summaryCombined) == nrow(res_ranksCombined_ord))

df_combined <- cbind(df_summaryCombined, res_ranksCombined_ord)
rownames(df_combined) <- NULL

# save .csv file
fn_combined <- file.path(dir_outputs, "topSVGs_nnSVG_avgRanks_combined.csv")
write.csv(df_combined, file = fn_combined, row.names = FALSE)


# -----------------------
# plot summaries of ranks
# -----------------------

# plot ranks per sample-part for all top SVGs from nnSVG

df_plot <- df_all %>% 
  filter(overall_rank <= 50) %>% 
  select(-c("gene_id", "gene_type")) %>% 
  pivot_longer(., cols = -c("gene_name", "overall_rank", "average_rank", "n_withinTop100"), 
               values_to = "rank_sample_part", names_to = "sample_part") %>% 
  mutate(sample_part = factor(sample_part, levels = sample_part_ids)) %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = factor(gene, levels = top100genes))


ggplot(df_plot) + 
  geom_boxplot(aes(x = rank_sample_part, y = gene, group = gene), 
               color = "navy") + 
  scale_y_discrete(limits = rev) + 
  labs(x = "rank (per tissue area)") + 
  ggtitle("Top SVGs (nnSVG)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "topSVGs_nnSVG_ranks_all")
ggsave(paste0(fn, ".pdf"), width = 6, height = 7)
ggsave(paste0(fn, ".png"), width = 6, height = 7)


# plot ranks per sample-part for "replicated" top SVGs from nnSVG

df_plot <- df_replicated %>% 
  select(-c("gene_id", "gene_type")) %>% 
  pivot_longer(., cols = -c("gene_name", "overall_rank", "average_rank", "n_withinTop100"), 
               values_to = "rank_sample_part", names_to = "sample_part") %>% 
  mutate(sample_part = factor(sample_part, levels = sample_part_ids)) %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = factor(gene, levels = top100genes))


ggplot(df_plot) + 
  geom_boxplot(aes(x = rank_sample_part, y = gene, group = gene), 
               color = "navy") + 
  scale_y_discrete(limits = rev) + 
  labs(x = "rank (per tissue area)") + 
  ggtitle("Top SVGs (nnSVG): replicated across tissue areas") + 
  theme_bw() + 
  theme(axis.text.y = element_text(face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "topSVGs_nnSVG_ranks_replicated")
ggsave(paste0(fn, ".pdf"), width = 6, height = 6)
ggsave(paste0(fn, ".png"), width = 6, height = 6)


# plot showing number of sample-parts where each top SVG is within top n SVGs
# for that sample-part

ranks_per_sample_part_Br8079 <- as.matrix(df_all[, 15:16])
ranks_per_sample_part_excBr8079 <- as.matrix(df_all[, c(7:14, 17:19)])
all_Br8079 <- apply(ranks_per_sample_part_Br8079, 1, function(r) !all(is.na(r))) & 
              apply(ranks_per_sample_part_excBr8079, 1, function(r) all(is.na(r)))

df_plot <- df_all %>% 
  mutate(all_Br8079 = all_Br8079) %>% 
  filter(overall_rank <= 50) %>% 
  select(c("gene_name", "n_withinTop100", "all_Br8079")) %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = factor(gene, levels = top100genes))


ggplot(df_plot, aes(x = n_withinTop100, y = gene, fill = all_Br8079)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("navy", "maroon"), name = "Br8079 only") + 
  scale_x_continuous(breaks = 0:13) + 
  scale_y_discrete(limits = rev) + 
  xlab("number of times within top 100 SVGs (13 tissue areas)") + 
  ggtitle("Top SVGs (nnSVG)") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "topSVGs_nnSVG_numberWithinTop100")
ggsave(paste0(fn, ".pdf"), width = 5, height = 7)
ggsave(paste0(fn, ".png"), width = 5, height = 7)

