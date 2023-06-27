#!/usr/bin/env Rscript
library(dplyr)
# Overlap circrna - mirna from circbank with mirna from DE
de_mir = read.table("/data/github/pca_network/results/mirna_intersection.txt", header = T, sep="\t")
db_mir = read.table("/data/github/pca_network/results/circrna_mirna.txt", header=T, sep="\t")

# -1, -2, -5p. -3p are all valid naming conventions - do not try to simplify.
#table(de_mir$miRNA %in% db_mir$miR_ID)

# output dataframe with
# circna avg lfc mirna avg lfc
# miRNAs in both predicted circrna - miRNA and the miRNA that are DE
intersection = Reduce(intersect, list(db_mir$miR_ID, de_mir$miRNA))
# sanity check
#a = db_mir$miR_ID[db_mir$miR_ID %in% de_mir$miRNA]
#b = de_mir$miRNA[de_mir$miRNA %in% db_mir$miR_ID]
#intersect(a, b)

# load LFC data
circrna = read.table("/data/github/pca_network/results/circrna_intersection.txt", header=T, sep="\t")
avg_circ = circrna %>% group_by(circbaseID) %>% summarise(avg_circ_lfc = mean(logFC))
avg_mirna = de_mir %>% group_by(miRNA) %>% summarise(avg_mirna_lfc = mean(logFC))

# subset circrna-mirna to match common miRNA from DE,DB
intersection = db_mir[which(db_mir$miR_ID %in% intersection),]
intersection = merge(intersection, avg_circ, by.x="circbase_ID", by.y="circbaseID")
intersection = merge(intersection, avg_mirna, by.x="miR_ID", by.y="miRNA")

intersection = intersection[,c(2,3,1,4)]
intersection = intersection[order(intersection$circbase_ID),]
write.table(intersection, "/data/github/pca_network/results/circrna_mirna_network.txt", sep="\t", row.names = F)
