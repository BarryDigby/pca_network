#!/usr/bin/env Rscript
library(dplyr)
# Overlap circrna - mirna from circbank with mirna from DE
de_mir = read.table("/data/github/pca_network/results/mirna_intersection.txt", header = T, sep="\t")
db_mir = read.table("/data/github/pca_network/results/circrna_mirna_db.txt", header=T, sep="\t")


# common de miRs with circrna - mirna database preds
intersection = Reduce(intersect, list(db_mir$mirna, de_mir$miRNA))
# sanity check
#a = db_mir$miR_ID[db_mir$miR_ID %in% de_mir$miRNA]
#b = de_mir$miRNA[de_mir$miRNA %in% db_mir$miR_ID]
#intersect(a, b)

# load LFC data
circrna = read.table("/data/github/pca_network/results/circrna_intersection.txt", header=T, sep="\t")
avg_circ = circrna %>% group_by(circbaseID) %>% summarise(avg_circ_lfc = mean(logFC))
avg_mirna = de_mir %>% group_by(miRNA) %>% summarise(avg_mirna_lfc = mean(logFC))

# return circrna- mirnas, using de mirnas to subset
intersection = db_mir[which(db_mir$mirna %in% intersection),]

# append the mirna fold change info 
intersection = merge(intersection, avg_mirna, by.x="mirna", by.y="miRNA")
# append circrna fold change info 
intersection = merge(intersection, avg_circ, by.x="circrna", by.y="circbaseID")

intersection = intersection[,c(1,4,2,3)]
intersection = intersection[order(intersection$circrna),]
write.table(intersection, "/data/github/pca_network/results/circrna_mirna_network.txt", quote=F,  sep="\t", row.names = F)

library(ggvenn)
pdf("/data/github/pca_network/results/demirna_circtargets.pdf", width=5, height=6)
ggvenn::ggvenn(list(DEmiRNA=de_mir$miRNA, circRNA_targets=db_mir$mirna),
               show_percentage = F, set_name_size = 5, fill_color = c("royalblue", "red1"),
               fill_alpha = 0.5)
dev.off()