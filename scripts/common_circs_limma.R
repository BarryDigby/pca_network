#!/usr/bin/env Rscript

GSE113153_res <- read.table("/data/github/pca_network/results/GSE113153_limma.txt", header=T)
GSE118959_Clone1_res <- read.table("/data/github/pca_network/results/GSE118959_Clone1_limma.txt", header = T)
GSE118959_Clone9_res <- read.table("/data/github/pca_network/results/GSE118959_Clone9_limma.txt", header=T)

probes <- read.table("/data/github/pca_network/data/arraystar_probes.csv", header=T, sep="\t")

common_circs <- Reduce(intersect, list(GSE113153_res$ID, GSE118959_Clone1_res$ID, GSE118959_Clone9_res$ID))

common_circs <- subset(probes, probes$probeID %in% common_circs)

regulation <- subset(GSE118959_Clone1_res, GSE118959_Clone1_res$ID %in% common_circs$probeID)
regulation <- ifelse(regulation$logFC > 0, "Up", "Down")
common_circs$Regulation <- regulation
common_circs <- subset(common_circs, select=c(probeID, circRNA, Alias, GeneSymbol, Regulation))
write.table(common_circs, "/data/github/pca_network/results/common_circs_limma.txt", quote=F, row.names=F, sep="\t")