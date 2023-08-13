#!/usr/bin/env Rscript

de_mirs = read.csv("/data/github/pca_network/results/mirna_intersection.txt", sep="\t", header=T)

GSE21036_mat = read.csv("/data/github/pca_network/mirna/GSE21036/heatmap_counts.txt", header=T, sep="\t")
GSE21036_meta = read.csv("/data/github/pca_network/mirna/GSE21036/heatmap_meta.txt", header=T, sep="\t")
GSE21036_mat = GSE21036_mat[which(GSE21036_mat$SystematicName %in% de_mirs$miRNA),]

GSE23022_mat = read.csv("/data/github/pca_network/mirna/GSE23022/heatmap_counts.txt", header=T, sep="\t")
GSE23022_meta = read.csv("/data/github/pca_network/mirna/GSE23022/heatmap_meta.txt", header=T, sep="\t")
GSE23022_mat = GSE23022_mat[which(GSE23022_mat$SystematicName %in% de_mirs$miRNA),]

GSE36803_mat = read.csv("/data/github/pca_network/mirna/GSE36803/heatmap_counts.txt", header=T, sep="\t")
GSE36803_meta = read.csv("/data/github/pca_network/mirna/GSE36803/heatmap_meta.txt", header=T, sep="\t")
GSE36803_mat = GSE36803_mat[which(GSE36803_mat$SystematicName %in% de_mirs$miRNA),]

GSE45604_mat = read.csv("/data/github/pca_network/mirna/GSE45604/heatmap_counts.txt", header=T, sep="\t")
GSE45604_meta = read.csv("/data/github/pca_network/mirna/GSE45604/heatmap_meta.txt", header=T, sep="\t")
GSE45604_mat = GSE45604_mat[which(GSE45604_mat$SystematicName %in% de_mirs$miRNA),]

GSE46738_mat = read.csv("/data/github/pca_network/mirna/GSE46738/heatmap_counts.txt", header=T, sep="\t")
GSE46738_meta = read.csv("/data/github/pca_network/mirna/GSE46738/heatmap_meta.txt", header=T, sep="\t")
GSE46738_mat = GSE46738_mat[which(GSE46738_mat$SystematicName %in% de_mirs$miRNA),]

LNCaP_mat = read.csv("/data/github/pca_network/mirna/LNCaP/heatmap_counts.txt", header=T, sep="\t")
LNCaP_meta = read.csv("/data/github/pca_network/mirna/LNCaP/heatmap_meta.txt", header=T, sep="\t")
LNCaP_mat = LNCaP_mat[which(LNCaP_mat$SystematicName %in% de_mirs$miRNA),]

TCGA_mat = read.csv("/data/github/pca_network/mirna/TCGA/heatmap_counts.txt", header=T, sep="\t")
TCGA_meta = read.csv("/data/github/pca_network/mirna/TCGA/heatmap_meta.txt", header=T, sep="\t")
TCGA_mat = TCGA_mat[which(TCGA_mat$SystematicName %in% de_mirs$miRNA),]