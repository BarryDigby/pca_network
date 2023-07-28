#!/usr/bin/env Rscript


library(pathfindR)
library(dplyr)

network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
input = network[,c(5,6)]
input$adj.P.val = 0.001 # used for filtering >  0.05 in pathfindr
colnames(input) = c("Gene.Symbol", "logFC", "adj.P.Val")
input = input %>% unique()

keggpathway = pathfindR::run_pathfindR(input, gene_sets = "KEGG", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=2, iterations=10)
write.table(keggpathway, "/data/github/pca_network/results/kegg.txt", quote = F, row.names = F, sep="\t")


reactome = pathfindR::run_pathfindR(input, gene_sets = "Reactome", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=2, iterations=10)
write.table(reactome, "/data/github/pca_network/results/reactome.txt", quote = F, row.names = F, sep="\t")

go = pathfindR::run_pathfindR(input, gene_sets = "GO-All", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=2, iterations=5)
write.table(go, "/data/github/pca_network/results/go.txt", quote = F, row.names = F, sep="\t")
