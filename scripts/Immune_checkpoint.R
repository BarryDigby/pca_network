#!/usr/bin/env Rscript

load("/data/github/pca_network/results/TCGA_DFS/model.RData")
source("https://raw.githubusercontent.com/BarryDigby/pca_network/main/data/geom_violin.R")
genes = c("IDO1", "LAG3", "CTLA4", "TNFRSF9", "ICOS", "CD80", "PDCD1LG2", "TIGIT", "CD70", "TNFSF9", "ICOSLG", "KIR3DL1", "CD86", "PDCD1", "LAIR1", "TNFRSF8", "TNFSF15", "TNFRSF14", "IDO2", "CD276", "CD40", "TNFRSF4", "TNFSF14", "HHLA2", "CD244", "CD274", "HAVCR2", "CD27", "BTLA", "LGALS9", "TMIGD2", "CD28", "CD48", "TNFRSF25", "CD40LG", "ADORA2A", "VTCN1", "CD160", "CD44", "TNFSF18", "TNFRSF18", "BTNL2", "CD200R1", "TNFSF4", "CD200", "NRP1", "C10orf54")
# read TCGA-PRAD
#x = read.csv("/data/github/pca_network/results/TCGA_mrna_logcpm.txt", header=T, sep="\t")
tcga_df = x
colnames(tcga_df) = gsub("\\.", "-", colnames(tcga_df))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
tcga_df = merge(tcga_df, mrna_attributes[,c(1,2)], by.x=0, by.y="ensembl_gene_id_version")
tcga_df = tibble::column_to_rownames(tcga_df, "external_gene_name")
tcga_df = tcga_df[,2:ncol(tcga_df)]
tcga_df = tcga_df[, rownames(cox)]
tcga_df = t(scale(t(tcga_df), scale = T, center=T))

tcga_df = as.data.frame(t(tcga_df))
tcga_df$category = cox$risk_category

# Reshape the dataframe into a tidy format
tidy_df <- tcga_df %>%
  pivot_longer(cols = colnames(tcga_df[,1:ncol(tcga_df)-1]), names_to = "Variable", values_to = "Score")

# remove ns for plot
tidy_df = tidy_df %>% filter(!Variable %in% c("ADORA2A", "CD160", "CD200", "CD274", "TNFSF15", "TNFSF4", "TNFSF9"))

# Create the violin plot
pdf("/data/github/pca_network/results/TCGA_DFS/Checkpoints.pdf", width=12, height=6)
ggplot(tidy_df, aes(x = Variable, y = Score, fill = category))+
  geom_split_violin(trim = FALSE, alpha = .4)+
  geom_boxplot(width = .2, alpha = .6,
               position = position_dodge(.25))+
  scale_fill_viridis_d(option = "D") +
  stat_summary(fun = "mean", geom = "point",
               position = position_dodge(width = 0.25)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .1,
               position = position_dodge(width = 0.25)) + theme_bw() +
  stat_compare_means(method = "t.test", label.y = 6, aes(label = after_stat(p.signif)))+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 10, face = "bold"))
dev.off()



library(immunedeconv)
tcga_df = x
colnames(tcga_df) = gsub("\\.", "-", colnames(tcga_df))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
tcga_df = merge(tcga_df, mrna_attributes[,c(1,2)], by.x=0, by.y="ensembl_gene_id_version")
tcga_df = tcga_df %>% filter(external_gene_name != "")
tcga_df = tibble::column_to_rownames(tcga_df, "external_gene_name")
tcga_df = tcga_df[,2:ncol(tcga_df)]
tcga_df = tcga_df[, rownames(cox)]

#tcga_df = 2^(tcga_df)
tcga_df = as.data.frame(t(scale(t(tcga_df), scale = T, center = T)))
res_mcp = as.data.frame(deconvolute(tcga_df, "mcp_counter", tumor = TRUE))
res_mcp = tibble::column_to_rownames(res_mcp, "cell_type")
res_mcp = as.data.frame(t(res_mcp))
res_mcp$category = cox$risk_category
#res_mcp = res_mcp[,which(colnames(res_mcp) != "Cancer associated fibroblast")]

tidy_df <- res_mcp %>%
  pivot_longer(cols = colnames(res_mcp[,1:ncol(res_mcp)-1]), names_to = "Variable", values_to = "Score")

tidy_df = tidy_df %>% filter(!Variable %in% c("Neutrophil", "T cell CD8+"))
# Create the violin plot
pdf("/data/github/pca_network/results/TCGA_DFS/MCPcounter.pdf", width=8, height=6)
ggplot(tidy_df, aes(x = Variable, y = Score, fill = category))+
  geom_split_violin(trim = FALSE, alpha = .4)+
  geom_boxplot(width = .2, alpha = .6,
               position = position_dodge(.25))+
  scale_fill_viridis_d(option = "D") +
  stat_summary(fun = "mean", geom = "point",
               position = position_dodge(width = 0.25)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .1,
               position = position_dodge(width = 0.25)) + theme_bw() +
  stat_compare_means(method = "t.test", label.y = 5, aes(label = after_stat(p.signif)))+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 10, face = "bold"))
dev.off()


