#!/usr/bin/env Rscript 
library(dplyr)
LNCAP_clone1 = read.csv("/data/github/pca_network/mrna/LNCaP/clone1_control.txt", sep="\t", header=T)
LNCAP_clone9 = read.csv("/data/github/pca_network/mrna/LNCaP/clone9_control.txt", sep="\t", header=T)
# surgical castration only - do not use 
GSE88752_primary = read.csv("/data/github/pca_network/mrna/GSE88752/primary_vs_control.txt", sep="\t", header=T)
# surgical castration + enz treatment 
GSE88752_secondary = read.csv("/data/github/pca_network/mrna/GSE88752/secondary_vs_control.txt", sep="\t", header=T)
GSE143408_day7 = read.csv("/data/github/pca_network/mrna/GSE143408/7_vs_0.txt", sep="\t", header=T)
GSE143408_day14 = read.csv("/data/github/pca_network/mrna/GSE143408/14_vs_0.txt", sep="\t", header=T)
GSE143408_day21 = read.csv("/data/github/pca_network/mrna/GSE143408/21_vs_0.txt", sep="\t", header=T)
GSE78201_6MO = read.csv("/data/github/pca_network/mrna/GSE78201/6mo_vs_control.txt", sep="\t", header=T)
TCGA = read.csv("/data/github/pca_network/mrna/TCGA/tumor_vs_normal.txt", sep="\t", header=T)

# rearrange DFs for merge: gene + basic limma info.
LNCAP_clone1$experiment = "LNCaP_clone1"
LNCAP_clone9$experiment = "LNCaP_clone9"
#GSE88752_primary = GSE88752_primary[,c(9,2:7)]
#colnames(GSE88752_primary)[1] = "Gene"
#GSE88752_primary$experiment = "GSE88752_primary"
GSE88752_secondary = GSE88752_secondary[,c(9,2:7)]
colnames(GSE88752_secondary)[1] = "Gene"
GSE88752_secondary$experiment = "GSE88752_secondary"
GSE143408_day7 = GSE143408_day7[,c(8,2:7)]
colnames(GSE143408_day7)[1] = "Gene"
GSE143408_day7$experiment = "GSE143408_day7"
GSE143408_day14 = GSE143408_day14[,c(8,2:7)]
colnames(GSE143408_day14)[1] = "Gene"
GSE143408_day14$experiment = "GSE143408_day14"
GSE143408_day21 = GSE143408_day21[,c(8,2:7)]
colnames(GSE143408_day21)[1] = "Gene"
GSE143408_day21$experiment = "GSE143408_day21"
colnames(GSE78201_6MO)[1] = "Gene"
GSE78201_6MO$experiment = "GSE78201_6MO"
TCGA$experiment = "TCGA"

# intersect (n exp groups -1 ) for groupby()
# Remove Primary treatment - this is catration but no enzalutamide.
intersection = rbind(LNCAP_clone1, LNCAP_clone9, GSE88752_secondary, GSE143408_day7, GSE143408_day14, GSE143408_day21, GSE78201_6MO, TCGA)
# must be in TCGA, Clone1 and validated in one other enz resistance dataset.
#intersection = intersection %>% group_by(Gene) %>% filter(length(unique(experiment))>=3 & 'LNCaP_clone1' %in% experiment & 'TCGA' %in% experiment) %>% ungroup()
#intersection = intersection %>% group_by(Gene) %>% filter(length(unique(experiment))>=2 & 'LNCaP_clone1' %in% experiment) %>% ungroup()
intersection = intersection %>%
  group_by(Gene) %>%
  mutate(
    ENZ = any(experiment %in% c("LNCaP_clone1", "LNCaP_clone9", "GSE88752_secondary", "GSE143408_day7", "GSE143408_day14", "GSE143408_day21", "GSE78201_6MO")),
    TVN = any(experiment %in% c("TCGA"))
  ) %>%
  filter( ENZ & TVN ) %>%
  ungroup() %>%
  dplyr::select(-ENZ, -TVN)

intersection = intersection %>% group_by(Gene) %>% filter(all(logFC > 0) | all(logFC < 0)) %>% ungroup()
intersection = intersection %>% group_by(Gene) %>% filter(all(adj.P.Val <= 0.05)) %>% ungroup()
write.table(intersection, "/data/github/pca_network/results/mrna_intersection.txt", sep="\t", row.names = FALSE, quote = FALSE)
