#!/usr.bin/env Rscript

mirna_network = read.table("/data/github/pca_network/results/circrna_mirna_network.txt", header=T, sep="\t")
mirna_mrna = read.table("/data/github/pca_network/data/mirna_mrna.txt", header=F, sep="\t")
# reduce load on R. one is missing,  hsa-miR-375 no -3p or -5p on it.  
mirna_mrna = mirna_mrna[which(mirna_mrna$V1 %in% mirna_network$miR_ID),]

#hgnc symbols
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
info <- biomaRt::getBM(attributes=c("hgnc_symbol",
                           "refseq_mrna"),
                            filters = c("refseq_mrna"),
                            values = unique(mirna_mrna$V2),
                            mart = mart,
                            useCache=TRUE)
mirna_mrna = merge(mirna_mrna, info, by.x="V2", by.y="refseq_mrna", all.x = T)
mirna_mrna = mirna_mrna[,c(2,4)]
colnames(mirna_mrna) = c("miRNA", "Gene")
mirna_mrna = mirna_mrna[order(mirna_mrna$miRNA),]
write.table(mirna_mrna, "/data/github/pca_network/results/mirna_mrna.txt", quote = F, sep="\t", row.names = F)
