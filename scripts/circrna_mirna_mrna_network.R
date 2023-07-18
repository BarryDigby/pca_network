#!/usr/bin/env Rscript

# circrna - mirna connections with logFC. append to result for final network.
# validated by DE in silico statistical tests.
circ_mirna = read.csv("/data/github/pca_network/results/circrna_mirna_network.txt", header=T, sep="\t")

# mirna - gene connections (database predictions)
mir_db = read.csv("/data/github/pca_network/results/mirna_mrna_db.txt", header=T, sep="\t")

# mrna intersection enz resistance signatures
de_mrna = read.csv("/data/github/pca_network/results/mrna_intersection.txt", header=T, sep="\t")

# intersect DE genes with predicted DB genes
intersection = Reduce(intersect, list(de_mrna$Gene, mir_db$Target))

# subset limma results with 'confirmed' mrnas
de_mrna = de_mrna[which(de_mrna$Gene %in% intersection),]
#length(unique(de_mrna$Gene)) == length(intersection)

# average expression. 
de_mrna = de_mrna[,c(1,2,8)]
de_mrna = de_mrna %>% group_by(Gene) %>% summarise(avg_mrna_lfc = mean(logFC))

# link back to circrna - mirna avglfc network.
# add mirna to de_mrna to use as merge key. 
mir_db = mir_db[which(mir_db$Target %in% de_mrna$Gene),]
mir_db = mir_db %>% unique()
mir_db = merge(mir_db, de_mrna, by.x="Target", by.y="Gene")

circ_mirna = merge(circ_mirna, mir_db, by.x = "miR_ID", by.y = "ID")
circ_mirna = circ_mirna[,c(2,3,1,4,5,6)]
colnames(circ_mirna)[5] = "gene_ID"

write.table(circ_mirna, "/data/github/pca_network/results/circ_mirna_mrna_network.txt", sep="\t", row.names = F, quote=F)
