#!/usr/bin/env Rscript

GSE113153 <- read.table("/data/github/pca_network/results/GSE113153.txt", header=T)
GSE118959_Clone1 <- read.table("/data/github/pca_network/results/GSE118959_Clone1.txt", header=T)
GSE118959_Clone9 <- read.table("/data/github/pca_network/results/GSE118959_Clone9.txt", header=T)


probes <- read.table("/data/github/pca_network/data/arraystar_probes.csv", header=T, sep="\t")
common_circs <- Reduce(intersect, list(GSE113153$Probe_ID, GSE118959_Clone1$Probe_ID, GSE118959_Clone9$Probe_ID))
common_circs <- subset(probes, probes$probeID %in% common_circs)


regulation <- subset(GSE118959_Clone1, GSE118959_Clone1$Probe_ID %in% common_circs$probeID)
regulation <- ifelse(regulation$Clone1_FC > 0, "Up", "Down")
common_circs$Regulation <- regulation
common_circs <- subset(common_circs, select=c(probeID, circRNA, Alias, GeneSymbol, Regulation))
write.table(common_circs, "/data/github/pca_network/results/common_circs.txt", quote=F, row.names=F, sep="\t")