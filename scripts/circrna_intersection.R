#!/usr/bin/env Rscript
library(dplyr)
GSE113153 = read.table("/data/github/pca_network/circrna/GSE113153/high_vs_low.txt", sep="\t", header=T)
GSE113153$experiment = "GSE113153"
CLONE9 = read.table("/data/github/pca_network/circrna/GSE118959/clone9_vs_control.txt", sep="\t", header=T)
CLONE9$experiment = "CLONE9"
CLONE1 = read.table("/data/github/pca_network/circrna/GSE118959/clone1_vs_control.txt", sep="\t", header=T)
CLONE1$experiment = "CLONE1"
master = rbind(GSE113153, CLONE1, CLONE9)

# probes annotated using 'arraystar_probes' but some Alias (i.e circBank IDs) are missing. Attempt to annotate manually with additional circBank file
arraystar = read.table("/data/github/pca_network/data/arraystar_probes.csv", sep="\t", header=T)
arraystar$locus = paste0(arraystar$chrom, sep=":", arraystar$txStart, sep="-", arraystar$txEnd, sep=":", arraystar$strand)
circbank = read.table("/data/github/pca_network/data/circBank_circrna_annotation.txt", header=T, sep="\t")
circbank$locus = paste0(circbank$position, sep=":", circbank$strand)

# Add arraystar locus using circRNA col, then locus as grouping key/.
subset_arraystar = subset(arraystar, select=c(circRNA, locus))
subset_circbank = subset(circbank, select=c(circbaseID, locus))
master = merge(master, subset_arraystar, by="circRNA")
master = merge(master, subset_circbank, by="locus", all.x=T)
# cannot recover the circrna below (458)
#sum(is.na(master$circbaseID))
#x = master[is.na(master$circbaseID), ]

# works like groupby() in python, use unique length to set condition 
# 2 or more datasets seems more robust before filtering by fold change direction.
intersection = master %>% group_by(circbaseID) %>% filter(length(unique(experiment))>=2) %>% ungroup()
intersection = na.omit(intersection)

# fold change consitency..
intersection = intersection %>% group_by(circbaseID) %>% filter(all(logFC>0) | all(logFC<0)) %>% ungroup

# tidy and save result. 
intersection = subset(intersection, select = c(circbaseID, logFC, AveExpr, t, P.Value, adj.P.Val, B, GeneSymbol, experiment))
write.table(intersection, "/data/github/pca_network/results/circrna_intersection.txt", quote = F, sep="\t", row.names = F)
