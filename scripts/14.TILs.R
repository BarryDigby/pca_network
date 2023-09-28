#!/usr/bin/env Rscript 

library(msigdbr)
library(GSVA)
library(ggpubr)
library(ggplot2)
library(pheatmap)
library(data.table)
library(tidyr)
library(dplyr)

load("/data/github/pca_network/results/TCGA_DFS/model.RData")

source("https://raw.githubusercontent.com/BarryDigby/pca_network/main/data/geom_violin.R")
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
load("/data/github/pca_network/results/TCGA_DFS/model.RData")
mrna_attributes = read.csv("/data/github/pca_network/mrna/TCGA/ENS_BIO_HGNC.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$biotype == "protein_coding"),]
tcga_df = merge(mrna_df, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_version_id")
rownames(tcga_df) <- make.names(tcga_df$Gene, unique = TRUE)
#tcga_df = tibble::column_to_rownames(tcga_df, "Gene")
tcga_df = tcga_df[,2:ncol(tcga_df)]
tcga_df = tcga_df[, rownames(cox)]


tcga_meta = data.frame(row.names = colnames(tcga_df),
                       Risk_strata = cox$risk_category)

tcga_meta$Risk_strata = ifelse(tcga_meta$Risk_strata == "High risk", "High", "Low")
y <- edgeR::DGEList(tcga_df)
logcpm <- edgeR::cpm(y, normalized.lib.sizes = T, log=TRUE)
tcga_mat = logcpm

gsva = GSVA::gsva(tcga_mat, TILs, method = "ssgsea", kcdf = "Gaussian",  min.sz = 1, max.sz = 500, mx.diff = TRUE, verbose = FALSE)

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

# remove those that are not significant 
#wide_df = wide_df %>% filter(!Variable %in% c("Activated B cell", "Activated CD4 T cell", "Activated dendritic cell", "CD56bright natural killer cell",
                              #             "Immature  B cell", "Monocyte", "Neutrophil", "Plasmacytoid dendritic cell", "T follicular helper cell", "Type 17 T helper cell"))

wide_df$Variable = gsub("memeory", "memory", wide_df$Variable)

pdf("/data/github/pca_network/results/TCGA_DFS/TILs.pdf", width=16, height=6)
ggplot(wide_df, aes(x = Variable, y = Score, fill = category))+
  geom_split_violin(trim = FALSE, alpha = .4)+
  geom_boxplot(width = .2, alpha = .6,
               position = position_dodge(.25))+
  scale_fill_viridis_d(option = "D") +
  stat_summary(fun = "mean", geom = "point",
               position = position_dodge(width = 0.25)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .1,
               position = position_dodge(width = 0.25)) + theme_bw() +
  stat_compare_means(method = "t.test", label.y = 0.8, aes(label = after_stat(p.signif)))+
  theme(axis.text.x = element_text(angle = 55, vjust=1, hjust=1, size = 10, face = "bold"),
        axis.title.x = element_blank())
dev.off()
