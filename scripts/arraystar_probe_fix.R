#!/usr/bin/env Rscript

# Arraystar probes are missing a lot of circBase IDs. Reconcile them using circBank annotation file.

arraystar = read.csv("/data/github/pca_network/data/arraystar_probes.csv", sep="\t", header=T)
arraystar$circbank_position = paste0(arraystar$chrom, sep=":", arraystar$txStart, sep="-", arraystar$txEnd)

circbank = read.csv("/data/github/pca_network/data/circBank_circrna_annotation.txt", header=T, sep="\t")
circbank = circbank[,c(2,3)]

arraystar = merge(arraystar, circbank, by.x="circbank_position", by.y="position", all.x = T)
arraystar = arraystar[,c(1:4,14,5:13)]
write.table(arraystar, "/data/github/pca_network/data/arraystar_updated_circbaseID.csv", sep="\t", quote = F, row.names = F)
