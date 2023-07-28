#!/usr.bin/env Rscript
library(dplyr)
# load circrna - mirna links (found by DE)
mirna_network = read.table("/data/github/pca_network/results/circrna_mirna_network.txt", header=T, sep="\t")
# load miRBase, miRNet predicted mirna - mrna genes. 
# subset this large dataset with your DE validated mirnas.
mirbase = read.table("/data/github/pca_network/data/mirbase.txt", header=F, sep="\t")
colnames(mirbase) = c("ID", "Target", "Score")
mirbase = mirbase[,c(1:2)]

# use miRNAs as input to website and download. 
mirnet = read.csv("/data/github/pca_network/data/mirnet2gene.csv", sep=",", header=T)
mirnet = mirnet[,c(1,3)]
mirnet$ID = gsub("mir", "miR", mirnet$ID)

mirtarbase = read.csv("/data/github/pca_network/data/miRTarBase.csv", header=T, sep=",")
mirtarbase = mirtarbase[which(mirtarbase$Species..miRNA.=="Homo sapiens"),]
mirtarbase = mirtarbase[,c(2,4)]
colnames(mirtarbase) = c("ID", "Target")

# reduce load on R. one is missing,  hsa-miR-375 no -3p or -5p on it? 
mirbase = mirbase[which(mirbase$ID %in% mirna_network$mirna),]
mirnet = mirnet[which(mirnet$ID %in% mirna_network$mirna),]
mirtarbase = mirtarbase[which(mirtarbase$ID %in% mirna_network$mirna),]

#hgnc symbols for miRBase
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
info <- biomaRt::getBM(attributes=c("hgnc_symbol",
                           "refseq_mrna"),
                            filters = c("refseq_mrna"),
                            values = unique(mirbase$Target),
                            mart = mart,
                            useCache=TRUE)
mirbase = merge(mirbase, info, by.x="Target", by.y="refseq_mrna", all.x = T)
mirbase = mirbase[,c(2,3)]
colnames(mirbase) = c("ID", "Target")

master = rbind(mirbase, mirnet, mirtarbase)
master = master %>% unique()

write.table(master, "/data/github/pca_network/results/mirna_mrna_db.txt", sep="\t", row.names = F, quote=F)
