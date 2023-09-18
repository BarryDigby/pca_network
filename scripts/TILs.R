#!/usr/bin/env Rscript 

load("/data/github/pca_network/results/TCGA_DFS/model.RData")
library(msigdbr)
library(GSVA)
library(ggpubr)
library(ggplot2)
library(pheatmap)
library(data.table)
library(tidyr)
library(dplyr)

# Stage all gene sets for analysis
file_path = "/data/github/pca_network/data/TILs.txt"

data <- fread(file_path, sep = "\t", header = FALSE, fill = TRUE)
out_df <- data.frame(
  name = data$V1,
  genes = apply(data[, -1, with = FALSE], 1, function(x) paste(x, collapse = ", "))
)

long_df <- separate_rows(out_df, genes, sep = ",")

colnames(long_df) <- c("name", "gene")

long_df = long_df
long_df <- long_df %>% filter(gene != " ") %>% filter(gene != "NA")


TILs = split(
  trimws(long_df$gene),
  long_df$name
)

TIL_genes = trimws(unique(long_df$gene))

# read TCGA-PRAD
x = read.csv("/data/github/pca_network/results/TCGA_mrna_logcpm.txt", header=T, sep="\t")
tcga_df = x
colnames(tcga_df) = gsub("\\.", "-", colnames(tcga_df))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% TIL_genes),]
tcga_df = merge(tcga_df, mrna_attributes[,c(1,2)], by.x=0, by.y="ensembl_gene_id_version")
tcga_df = tibble::column_to_rownames(tcga_df, "external_gene_name")
tcga_df = tcga_df[,2:ncol(tcga_df)]
tcga_df = tcga_df[, rownames(cox)]
tcga_mat = as.matrix(tcga_df)

gsva = GSVA::gsva(tcga_mat, TILs, method = "gsva", kcdf = "Gaussian",  min.sz = 1, max.sz = 500, mx.diff = TRUE, verbose = FALSE)

order_cols = data.frame(id = rownames(cox),
                        risk = cox$risk_category)

order_cols$risk = sort(order_cols$risk, decreasing = TRUE)

ann_col = data.frame(row.names = order_cols$id,
                     Risk_strata = order_cols$risk)

col <- c("red2", "royalblue")
names(col) <- c("High risk", "Low risk")
ann_clr <- list(Risk_strata = col)

mat = t(scale(t(gsva), scale = T, center = T))

pheatmap(mat,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "none",
         annotation_col=ann_col,
         annotation_colors = ann_clr,
         color = hcl.colors(100, "RdBu",rev=T))

# heatmap is messy, make plots w t-test of pathways instead 
gsva = as.data.frame(t(gsva))
gsva$category = order_cols$risk

wide_df = gsva %>%
  pivot_longer(cols = colnames(gsva[1:ncol(gsva)-1]), names_to = "Variable", values_to = "Score")


pdf("/data/github/pca_network/results/TCGA_DFS/TILs.pdf", width=16, height=6)
ggplot(wide_df, aes(x = Variable, y = Score, fill = category))+
  geom_split_violin(trim = FALSE, alpha = .4)+
  geom_boxplot(width = .2, alpha = .6,
               position = position_dodge(.25))+
  scale_fill_viridis_d(option = "E") +
  stat_summary(fun = "mean", geom = "point",
               position = position_dodge(width = 0.25)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .1,
               position = position_dodge(width = 0.25)) + theme_bw() +
  stat_compare_means(method = "t.test", label.y = 1.5, aes(label = after_stat(p.signif)))+
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