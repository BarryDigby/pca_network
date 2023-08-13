#!/usr/bin/env Rscript
library(GEOquery)
# GSE113153
gset <- getGEO("GSE113153", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21825", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
meta <- data.frame(Treatment = c(rep("Gleason > 8", 5), rep("Gleason < 8", 5)),
                   Replicate = rep(c("1", "2", "3", "4", "5"), 2),
                   Batch = rep("GSE113153", 10))
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
                   Replicate = rep(c("1", "2", "3"), 3),
                   Batch = rep("GSE118959", 9))
rownames(meta1) = colnames(exprs(gset))
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

# off by 2? short by 256 - 279 ? Clone 9 perhaps..
exp1 = exp1[which(rownames(exp1) %in% rownames(exp)),]

mat = cbind(exp, exp1)
master_meta = rbind(meta, meta1)
# limma removeBatch
mm <- model.matrix(~ 0 + Treatment + Batch, data=master_meta)
mat <- limma::removeBatchEffect(mat, batch=master_meta$Batch, design=mm)

library(pheatmap)
ann_col = data.frame(row.names = rownames(master_meta),
                     Sample = master_meta$Treatment)
col <- c("chartreuse", "skyblue", "royalblue", "purple", "orange")
col = rev(hcl.colors(5, "Dark 3"))
names(col) <- unique(ann_col$Sample)
ann_clr <- list(Sample = col) 

mat = t(scale(t(mat), center = T, scale = T))

pdf("/data/github/pca_network/results/circrna_heatmap.pdf", width=6, height = 6)
pheatmap(mat, cluster_cols = F, cluster_rows = T, scale = "none", show_rownames = F,
         annotation_col = ann_col, color = hcl.colors(100, "RdBu",rev=T),
         annotation_colors = ann_clr, show_colnames = F)
dev.off()
