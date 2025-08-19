#!/usr/bin/env Rscript

# Script to intersect circRNAs from both studies

library(dplyr)
GSE113153 = read.table("/data/github/pca_network/circrna/GSE113153/high_vs_low.txt", sep="\t", header=T)
GSE113153$experiment = "GSE113153"
CLONE9 = read.table("/data/github/pca_network/circrna/GSE118959/clone9_vs_control.txt", sep="\t", header=T)
CLONE9$experiment = "CLONE9"
CLONE1 = read.table("/data/github/pca_network/circrna/GSE118959/clone1_vs_control.txt", sep="\t", header=T)
CLONE1$experiment = "CLONE1"
master = rbind(GSE113153, CLONE1, CLONE9)

# quick sanity check - performs the same as dplyr group_by()
#combined_circs = c(CLONE1$circbaseID, CLONE9$circbaseID, GSE113153$circbaseID)
#circ_counts = table(combined_circs)
#n_2 = names(circ_counts[circ_counts >= 2])
#master = master[which(master$circbaseID %in% n_2),]

# circ must show up in enz treat + t v n 
intersection = master %>% 
  group_by(circbaseID) %>% 
  filter(
    "GSE113153" %in% experiment & 
      (("CLONE1" %in% experiment | "CLONE9" %in% experiment))
    ) %>%
      ungroup()
#length(unique(experiment))>=2 & "GSE113153" %in% experiment) %>% ungroup()
intersection = na.omit(intersection)

# fold change consistency (direction) across experiments
intersection = intersection %>% group_by(circbaseID) %>% filter(all(logFC > 0) | all(logFC < 0)) %>% ungroup()
intersection = intersection %>% group_by(circbaseID) %>% filter(all(adj.P.Val <= 0.05)) %>% ungroup()

# tidy and save result. 
intersection = subset(intersection, select = c(circbaseID, logFC, AveExpr, t, P.Value, adj.P.Val, B, GeneSymbol, experiment))

write.table(intersection, "/data/github/pca_network/results/circrna_intersection.txt", quote = F, sep="\t", row.names = F)

library(ggvenn)
pdf("/data/github/pca_network/results/circrna_intersection.pdf", width=5, height=6)
ggvenn::ggvenn(list(Clone_1 = unique(CLONE1$circbaseID), Clone_9 = unique(CLONE9$circbaseID), GSE113153 = unique(GSE113153$circbaseID)), show_percentage = F)
dev.off()
