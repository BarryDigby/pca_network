#!/usr/bin/env Rscript
library(dplyr)
circrna_intersection = read.table("/data/github/pca_network/results/circrna_intersection.txt", header=T, sep="\t")

## circbank miRNA predictions
circbank = read.table("/data/github/pca_network/data/circbank_miRNA_all_v1.txt", sep="\t", header=T)
circbank = circbank[,c(2,1)]
colnames(circbank) = c("circrna", "mirna")
# DE circs in circbank preds
circbank_in_de = circbank[which(circbank$circrna %in% circrna_intersection$circbaseID),]

# CDSC mirna preds - fromatted to match circbase IDs
cdsc_circbase = read.csv("/data/github/pca_network/data/cirbase_cdsc_merge-back-to-CDSC.txt", header=F, sep="\t")
cdsc_mirs = read.csv("/data/github/pca_network/data/CDSC_subset_cirbase_entries.txt", header=F, sep="\t")
cdsc = merge(cdsc_circbase, cdsc_mirs, by="V1")
cdsc = cdsc[,c(2,3)]
colnames(cdsc) = c("circrna", "mirna")
cdsc$mirna = gsub("^(\\d)", "miR-\\1", cdsc$mirna)
cdsc$mirna = paste0("hsa-", cdsc$mirna)
cdsc_in_de = cdsc[which(cdsc$circrna %in% circrna_intersection$circbaseID),]

master = rbind(cdsc_in_de, circbank_in_de)
master = master %>% unique()

write.table(master, "/data/github/pca_network/results/circrna_mirna_db.txt", sep="\t", row.names = F, quote = F)
