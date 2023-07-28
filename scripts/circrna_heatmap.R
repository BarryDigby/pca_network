#!/usr/bin/env Rscript
library(GEOquery)
# GSE113153
gset <- getGEO("GSE113153", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21825", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
meta <- data.frame(Gleason = c(rep("High", 5), rep("Low", 5)),
                   Replicate = rep(c("1", "2", "3", "4", "5"), 2))
rownames(meta) = colnames(exprs(gset))
exp = as.data.frame(exprs(gset))
probes = read.table("/data/github/pca_network/data/arraystar_probes.csv", sep="\t", header=T)
probes = probes[,c(1,3)]
exp = merge(exp, probes, by.x=0, by.y="probeID")

# GSE118959
gset <- getGEO("GSE118959", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21825", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
meta1 <- data.frame(Treatment = c(rep("Control", 3), rep("Clone1", 3), rep("Clone9", 3)),
                   Replicate = rep(c("1", "2", "3"), 3))
rownames(meta) = colnames(exprs(gset))
exp1 = as.data.frame(exprs(gset))
exp1 = merge(exp1, probes, by.x=0, by.y="probeID", all.x = T)

# load overlaps 
intersection = read.table("/data/github/pca_network/results/circrna_intersection.txt", header=T, sep="\t")
library(tidyverse)

exp = exp[which(exp$Alias %in% intersection$circbaseID), ]
exp1 = exp1[which(exp1$Alias %in% intersection$circbaseID), ]
exp = exp[,c(2:ncol(exp))]
exp1 = exp1[,c(2:ncol(exp1))]
exp = exp %>% remove_rownames %>% column_to_rownames(var="Alias")
exp1 = exp1 %>% remove_rownames %>% column_to_rownames(var="Alias")

mat = cbind(exp, exp1)
mat = as.data.frame(scale(mat, scale = T, center = T))
library(pheatmap)

pheatmap(exp1, cluster_cols = T, cluster_rows = T, scale = "column", show_rownames = T)
