#!/usr/bin/env Rscript
library(dplyr)
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


circ_mirna = merge(circ_mirna, mir_db, by.x = "mirna", by.y = "ID")
circ_mirna = circ_mirna[,c(2,3,1,4,5,6)]
colnames(circ_mirna)[5] = "gene_ID"

# Filter rows where avg_circrna_lfc is positive, avg_mirna_lfc is negative, and avg_mrna_lfc is positive
condition_1 <- circ_mirna$avg_circ_lfc > 0 & circ_mirna$avg_mirna_lfc < 0 & circ_mirna$avg_mrna_lfc > 0
result_1 <- circ_mirna[condition_1, ]

# Filter rows where avg_circrna_lfc is negative, avg_mirna_lfc is positive, and avg_mrna_lfc is negative
condition_2 <- circ_mirna$avg_circ_lfc < 0 & circ_mirna$avg_mirna_lfc > 0 & circ_mirna$avg_mrna_lfc < 0
result_2 <- circ_mirna[condition_2, ]

master = rbind(result_1, result_2)
colnames(master)[5] = "mrna"

write.table(master, "/data/github/pca_network/results/circ_mirna_mrna_network.txt", sep="\t", row.names = F, quote=F)

# format for cytoscape

# fold change filtering 

#master <- master %>% filter(abs(avg_circ_lfc) > 1, abs(avg_mirna_lfc) > 1, abs(avg_mrna_lfc) > 1)

edgelist <- data.frame(source = character(),
                       target = character(),
                       stringsAsFactors = FALSE)

nodelist <- data.frame(id = character(),
                       type = character(),
                       logFC = numeric())

for (i in 1:nrow(master)) {
  row <- master[i,]
  circ <- row$circrna
  circ_fc <- row$avg_circ_lfc
  mir <- row$mirna
  mir_fc <- row$avg_mirna_lfc
  gene <- row$mrna
  gene_fc <- row$avg_mrna_lfc
  
  nodelist <- rbind(nodelist, data.frame(id = circ, type = "circRNA", logFC = circ_fc))
  nodelist <- rbind(nodelist, data.frame(id = mir, type = "miRNA", logFC = mir_fc))
  nodelist <- rbind(nodelist, data.frame(id = gene, type = "mRNA", logFC = gene_fc))
  edgelist <- rbind(edgelist, data.frame(source = circ, target = mir, stringsAsFactors = FALSE))
  edgelist <- rbind(edgelist, data.frame(source = mir, target = gene, stringsAsFactors = FALSE))
}

rownames(nodelist) <- NULL
rownames(edgelist) <- NULL

edgelist = edgelist %>% unique()
nodelist = nodelist %>% unique()

write.table(edgelist, "/data/github/pca_network/results/edgelist.csv", sep="\t", row.names = F, quote=F)
write.table(nodelist, "/data/github/pca_network/results/nodelist.csv", sep="\t", row.names = F, quote=F)
