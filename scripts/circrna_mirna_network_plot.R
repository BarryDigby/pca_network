#!/usr/bin/env Rscript

#first 10. . . 

library(igraph)
intersection = read.table("/data/github/pca_network/results/circrna_mirna_network.txt",header=T, sep="\t")
topology = intersection[1:10,c(1,3)]
fold_change = data.frame(node=c(intersection$circbase_ID[1:10],intersection$miR_ID[1:10]),
                         fc=c(intersection$avg_circ_lfc[1:10],intersection$avg_mirna_lfc[1:10]))
fold_change = fold_change[!duplicated(fold_change$node),]
g2 = graph.data.frame(topology, vertices=fold_change, directed = T)
g = simplify(g2)
palette = colorRampPalette(c("green", "white", "red"))(n=10000)
sf = max(abs(fold_change$fc))
node.colors = (fold_change$fc+sf)/(2*sf) * 10000
layout = layout_with_graphopt(g)
plot(g, vertex.color=palette[node.colors], vertex.label.dist=3.5, layout=layout, edge.arrow.size=0.5, edge.width=1)