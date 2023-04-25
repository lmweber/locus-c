## For prettier violin plots with the aesthetics we typically go for
library(ggplot2)

plotExpressionCustom <- function(sce, features, features_name, anno_name = "cellType",
                                 point_alpha=0.2, point_size=0.7, ncol=2, xlab = NULL,
                                 exprs_values = "logcounts", scales = "free_y", swap_rownames=NULL){
  scater::plotExpression(sce, 
                         exprs_values = exprs_values, 
                         features = features,
                         x = anno_name, 
                         colour_by = anno_name,
                         ncol = ncol,
                         xlab = xlab,
                         point_alpha = point_alpha, 
                         point_size = point_size,
                         add_legend = F,
                         scales = scales,
                         swap_rownames = swap_rownames) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(face = "italic")) +  
  ggtitle(label=paste0(features_name, " markers"))
}

# MNT note: if want to add a custom title, just call this function and can 'overwrite' that with
#           another `+ ggtitle(...)`


