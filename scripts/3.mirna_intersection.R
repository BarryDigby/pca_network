#!/usr/bin/env Rscript
library(dplyr)
GSE21036 <- read.table("/data/github/pca_network/mirna/GSE21036/tumor_vs_normal.txt", header=T, sep="\t")
colnames(GSE21036)[1] = "miRNA"
GSE21036$experiment = "GSE21036"
GSE23022 <- read.table("/data/github/pca_network/mirna/GSE23022/tumor_vs_normal.txt", header=T, sep="\t")
GSE23022$miRNA = GSE23022$new_id
GSE23022$experiment = "GSE23022"
GSE23022 = subset(GSE23022, select=-c(new_id))
GSE36803 <- read.table("/data/github/pca_network/mirna/GSE36803/PCA_vs_Benign.txt", header=T, sep="\t")
GSE36803$miRNA = GSE36803$new_id
GSE36803 = subset(GSE36803, select=-c(new_id))
GSE36803$experiment = "GSE36803"
CLONE1 = read.csv("/data/github/pca_network/mirna/LNCaP/clone1_control.txt", header=T, sep="\t")
CLONE1$experiment = "CLONE1"
CLONE9 = read.csv("/data/github/pca_network/mirna/LNCaP/clone9_control.txt", header=T, sep="\t")
CLONE9$experiment = "CLONE9"
GSE45604 = read.csv("/data/github/pca_network/mirna/GSE45604/tumor_vs_normal.txt", header=T, sep="\t")
GSE45604$experiment = "GSE45604"
GSE46738 = read.csv("/data/github/pca_network/mirna/GSE46738/tumor_vs_normal.txt", sep="\t", header=T)
GSE46738$experiment = "GSE46738"
TCGA = read.csv("/data/github/pca_network/mirna/TCGA/tumor_v_normal.txt", sep="\t", header=T)
TCGA$experiment ="TCGA"

# intersect
intersection = rbind(GSE21036, GSE23022, GSE36803, CLONE1, CLONE9, GSE46738, GSE46738, TCGA)
intersection = intersection %>%
  group_by(miRNA) %>%
  mutate(
    ENZ = any(experiment %in% c("CLONE1", "CLONE9")),
    TVN = any(experiment %in% c("TCGA", "GSE46738", "GSE45604", "GSE36803", "GSE23022", "GSE21036"))
  ) %>%
  filter(ENZ & TVN)%>%
  ungroup() %>%
  dplyr::select(-ENZ, -TVN)
#intersection = intersection %>% group_by(miRNA) %>% filter(length(unique(experiment))>=2 & "CLONE1" %in% experiment) %>% ungroup()
intersection = intersection %>% group_by(miRNA) %>% filter(all(logFC>0) | all(logFC< 0)) %>% ungroup()
intersection = intersection %>% group_by(miRNA) %>% filter(all(adj.P.Val <= 0.05)) %>% ungroup()

write.table(intersection, "/data/github/pca_network/results/mirna_intersection.txt", sep="\t", row.names = FALSE, quote = FALSE)

library(ggvenn)
pdf("/data/github/pca_network/results/mirna_intersection.pdf", width=7, height=6)
ggvenn::ggvenn(list(GSE21036=GSE21036$miRNA, GSE23022=GSE23022$miRNA, GSE36803=GSE36803$miRNA, GSE45604=GSE45604$miRNA, GSE46738=GSE46738$miRNA, TCGA=TCGA$miRNA),
               show_percentage = F, set_name_size = 5)
dev.off()