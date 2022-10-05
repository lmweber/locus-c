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
dir_plots <- here("plots", "09_nnSVG")

# directory to save outputs
dir_outputs <- here("outputs", "09_nnSVG")


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
avg_ranks <- sort(rowMeans(res_ranks, na.rm = TRUE))


# summary table
df_summary <- data.frame(
  gene_id = names(avg_ranks), 
  gene_name = rowData(spe)[names(avg_ranks), "gene_name"], 
  gene_type = rowData(spe)[names(avg_ranks), "gene_type"], 
  avg_rank = unname(avg_ranks), 
  row.names = names(avg_ranks)
)

n <- 50
head(df_summary, n)


# top n genes
top_n_genes <- df_summary$gene_name[1:n]


# ------------
# plot top SVG
# ------------

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


# ----------------
# save spreadsheet
# ----------------

# save spreadsheet containing list of top SVGs

res_ranks_ord <- res_ranks[rownames(df_summary), ]

stopifnot(all(rownames(df_summary) == rownames(res_ranks_ord)))
stopifnot(nrow(df_summary) == nrow(res_ranks_ord))

df <- cbind(df_summary, res_ranks_ord)

rownames(df) <- NULL


# save .csv file
fn <- file.path(dir_outputs, "topSVGs_nnSVG_avgRanks.csv")
write.csv(df, file = fn, row.names = FALSE)


# -------------
# summary plots
# -------------

# plot ranks per sample-part per gene

df_plot <- df[1:n, ] %>% 
  select(-c("gene_id", "gene_type")) %>% 
  pivot_longer(., cols = -c("gene_name", "avg_rank"), values_to = "rank", names_to = "sample_part") %>% 
  mutate(sample_part = factor(sample_part, levels = sample_part_ids)) %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = factor(gene, levels = top_n_genes))


pal <- unname(c(palette.colors(8, "Okabe-Ito")[2:8], 
                palette.colors(36, "Polychrome 36")[c(3:4, 6:36)]))[1:13]

ggplot(df_plot) + 
  geom_boxplot(aes(x = rank, y = gene, group = gene), 
               alpha = 0, outlier.shape = NA) + 
  geom_point(aes(x = rank, y = gene, group = gene, color = sample_part), 
             shape = 1, stroke = 0.75) + 
  scale_y_discrete(limits = rev) + 
  scale_color_manual(values = pal) + 
  ggtitle("Top SVGs (nnSVG)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "topSVGs_nnSVG_ranks")
ggsave(paste0(fn, ".pdf"), width = 8, height = 7)
ggsave(paste0(fn, ".png"), width = 8, height = 7)


# plot number of sample-parts within top 50 SVGs per gene

length(sample_part_ids)
dim(df)

ranks_per_sample_part <- as.matrix(df[1:n, 5:17])
ranks_per_sample_part_notBr8079 <- as.matrix(df[1:n, c(5:12, 15:17)])

n_nonNA <- unname(apply(ranks_per_sample_part, 1, function(r) length(sample_part_ids) - sum(is.na(r))))
all_Br8079 <- apply(ranks_per_sample_part_notBr8079, 1, function(r) all(is.na(r)))

df_plot <- df[1:n, ] %>% 
  mutate(n_nonNA = n_nonNA) %>% 
  mutate(all_Br8079 = all_Br8079) %>% 
  select(c("gene_name", "n_nonNA")) %>% 
  rename(gene = gene_name) %>% 
  mutate(gene = factor(gene, levels = top_n_genes))

ggplot(df_plot, aes(x = n_nonNA, y = gene, fill = all_Br8079)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("navy", "maroon"), name = "Br8079 only") + 
  scale_x_continuous(breaks = 0:13) + 
  scale_y_discrete(limits = rev) + 
  xlab("number of sample-parts within top 50 SVGs") + 
  ggtitle("Top SVGs (nnSVG)") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(face = "italic"), 
        axis.title.y = element_blank())

fn <- here(dir_plots, "topSVGs_nnSVG_numberWithinTop")
ggsave(paste0(fn, ".pdf"), width = 5, height = 7)
ggsave(paste0(fn, ".png"), width = 5, height = 7)

