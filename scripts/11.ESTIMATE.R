#!/usr/bin/env Rscript 


################################
# ESTIMATE
################################

# skip to end of script and load function. 
# call from github if you remember to upload the function

library(tidyr)
library(ggpubr)
library(estimate)

load("/data/github/pca_network/results/TCGA_DFS/model.RData")
estimate = read.csv("/data/github/pca_network/data/ESTIMATE.txt", sep="\t", header = T, row.names = "ID")
rownames(estimate) = paste(rownames(estimate), "A", sep = "")
table(rownames(estimate) %in% rownames(cox))
estimate = merge(estimate, subset(cox, select=c(risk_category)), by=0)
estimate = tibble::column_to_rownames(estimate, "Row.names")

# Reshape the dataframe into a tidy format
tidy_df <- estimate %>%
  pivot_longer(cols = c(Stromal_score, Immune_score, ESTIMATE_score), names_to = "Variable", values_to = "Score")

# Create the violin plot
pdf("/data/github/pca_network/results/TCGA_DFS/ESTIMATE.pdf", width=8, height=6)
ggplot(tidy_df, aes(x = Variable, y = Score, fill = risk_category))+
  geom_split_violin(trim = FALSE, alpha = .4)+
  geom_boxplot(width = .2, alpha = .6,
               position = position_dodge(.25))+
  scale_fill_viridis_d(option = "E") +
  stat_summary(fun = "mean", geom = "point",
               position = position_dodge(width = 0.25)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .1,
               position = position_dodge(width = 0.25)) + theme_bw() +
  stat_compare_means(method = "t.test", label.y = 4250, aes(label = after_stat(p.signif)))+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 10, face = "bold"))
dev.off()

GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}