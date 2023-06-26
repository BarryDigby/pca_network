#!/usr/bin/env Rscript
library(dplyr)
GSE21036 <- read.table("/data/github/pca_network/mirna/GSE21036/tumor_vs_normal.txt", header=T, sep="\t")
colnames(GSE21036)[1] = "miRNA"
GSE21036$experiment = "GSE21036"
GSE23022 <- read.table("/data/github/pca_network/mirna/GSE23022/tumor_vs_normal.txt", header=T, sep="\t")
GSE23022$experiment = "GSE23022"
GSE36803 <- read.table("/data/github/pca_network/mirna/GSE36803/PCA_vs_Benign.txt", header=T, sep="\t")
GSE36803$experiment = "GSE36803"

# remove _st string
GSE23022$miRNA = gsub("_st", "", GSE23022$miRNA)
GSE36803$miRNA = gsub("_st", "", GSE36803$miRNA)

# star not *
GSE21036$miRNA = gsub("\\*", "-star", GSE21036$miRNA)

# intersect
intersection = rbind(GSE21036, GSE23022, GSE36803)
intersection = intersection %>% group_by(miRNA) %>% filter(length(unique(experiment))>=2) %>% ungroup()
intersection = intersection %>% group_by(miRNA) %>% filter(all(logFC>0) | all(logFC<0)) %>% ungroup
write.table(intersection, "/data/github/pca_network/results/mirna_intersection.txt", sep="\t", row.names = FALSE, quote = FALSE)


circrna_targets = read.table("/data/github/pca_network/results/circrna_mirna.txt", sep="\t", header=T)

circrna_targets = circrna_targets[which(circrna_targets$miR_ID %in% unique(intersection$miRNA)), ]
