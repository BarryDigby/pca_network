#!/usr/bin/env Rscript 


################################
# ESTIMATE
################################

library(tidyr)
library(ggpubr)
library(estimate)

load("/data/github/pca_network/results/TCGA_DFS/model.RData")
source("https://raw.githubusercontent.com/BarryDigby/pca_network/main/data/geom_violin.R")
estimate = read.csv("/data/github/pca_network/data/ESTIMATE.txt", sep="\t", header = T, row.names = "ID")
rownames(estimate) = paste(rownames(estimate), "A", sep = "")
table(rownames(estimate) %in% rownames(cox))
estimate = merge(estimate, subset(cox, select=c(risk_category)), by=0)
estimate = tibble::column_to_rownames(estimate, "Row.names")
colnames(estimate)[1:3] = c("Stromal", "Immune", "ESTIMATE")
# Reshape the dataframe into a tidy format
tidy_df <- estimate %>%
  pivot_longer(cols = c(Stromal, Immune, ESTIMATE), names_to = "Variable", values_to = "Score")

# Create the violin plot
pdf("/data/github/pca_network/results/TCGA_DFS/ESTIMATE.pdf", width=6, height=6)
ggplot(tidy_df, aes(x = Variable, y = Score, fill = risk_category))+
  geom_split_violin(trim = FALSE, alpha = .4)+
  geom_boxplot(width = .2, alpha = .6,
               position = position_dodge(.25))+
  scale_fill_viridis_d(option = "D") +
  stat_summary(fun = "mean", geom = "point",
               position = position_dodge(width = 0.25)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .1,
               position = position_dodge(width = 0.25)) + theme_bw() +
  stat_compare_means(method = "t.test", label.y = 4250, aes(label = after_stat(p.signif)))+
  theme(axis.text.x = element_text(angle = 55, vjust=1, hjust=1, size = 10, face = "bold"),
        axis.title.x = element_blank())
dev.off()