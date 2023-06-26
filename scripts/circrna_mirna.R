#!/usr/bin/env Rscript

circrna_intersection = read.table("/data/github/pca_network/results/circrna_intersection.txt", header=T, sep="\t")
circrna_mirna = read.table("/data/github/pca_network/data/circbank_miRNA_all_v1.txt", sep="\t", header=T)
circrna_mirna = circrna_mirna[,c(1,2)]
circrna_mirna = circrna_mirna[which(circrna_mirna$circbase_ID %in% unique(circrna_intersection$circbaseID)), ]
write.table(circrna_mirna, "/data/github/pca_network/results/circrna_mirna.txt", sep="\t", row.names = F, quote = F)
