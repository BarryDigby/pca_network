#!/usr/bin/env Rscript

# troubles with newer java installation. See function below. 
Sys.setenv(JAVA_HOME= "/usr/lib/jvm/java-17-openjdk-amd64")
#install.packages("pathfindR")
library(pathfindR)
library(dplyr)
library(ggplot2)
library(viridis)

network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
input = network[,c(5,6)]
input$adj.P.val = 0.001 # used for filtering >  0.05 in pathfindr
colnames(input) = c("Gene.Symbol", "logFC", "adj.P.Val")
input = input %>% unique()

# Run analysis and save. load for future sessions. 
#keggpathway = pathfindR::run_pathfindR(input, gene_sets = "KEGG", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=4, iterations=100)
#write.table(keggpathway, "/data/github/pca_network/results/kegg.txt", quote = F, row.names = F, sep="\t")

#reactome = pathfindR::run_pathfindR(input, gene_sets = "Reactome", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=4, iterations=100)
#write.table(reactome, "/data/github/pca_network/results/reactome.txt", quote = F, row.names = F, sep="\t")

#go = pathfindR::run_pathfindR(input, gene_sets = "GO-All", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=4, iterations=100)
#write.table(go, "/data/github/pca_network/results/go.txt", quote = F, row.names = F, sep="\t")

kegg     = read.csv("/data/github/pca_network/results/kegg.txt", header=T, sep="\t")
go       = read.csv("/data/github/pca_network/results/go.txt", header=T, sep="\t")
reactome = read.csv("/data/github/pca_network/results/reactome.txt", header=T, sep="\t")

melt = rbind(go,kegg,reactome)
melt = melt[which(melt$occurrence>50),]

# network plot for manuscript
pdf("/data/github/pca_network/results/pathfindR/network.pdf", width=22, height=16)
melt = cluster_enriched_terms(melt, plot_clusters_graph = T, use_description = T, use_active_snw_genes = T, vertex.label.cex=1, vertex.size.scaling=1.5)
dev.off()


# GO pathways
go = go[which(go$occurrence > 50),]
go_clst = cluster_enriched_terms(go, plot_dend = T, plot_clusters_graph = T, use_description = T, use_active_snw_genes = T)
save = c()
for(i in 1:nrow(go_clst)){
  row = as.data.frame(go_clst[i,])
  size = (length(strsplit(row$Up_regulated, ", ")[[1]]) + length(strsplit(row$Down_regulated, ", ")[[1]]))
  save = c(save, size)
}
go_clst$size = save

go_clst = go_clst[order(-go_clst$Fold_Enrichment),]

# plot top 10
go_clst = go_clst[1:10,]

ggplot(go_clst, aes(y=Term_Description,x=Fold_Enrichment)) +
  geom_point(aes(color=lowest_p,size=size)) +
  scale_color_gradientn(colours = rocket(100, begin = 0, end = 0.5, alpha = 0.8, direction = -1)) +
  labs( title = "GO Pathways",
        subtitle = "GO-BP, GO-MF, GO-CC",
        x='Fold Enrichment', y=NULL,
        color='P-value',size='Gene\nnumber'
  ) + theme_linedraw() + scale_size(range = c(3, 10)) + ggplot2::facet_grid(go_clst$Cluster ~ ., scales = "free_y", space = "free", drop = TRUE)  +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size=12),
    axis.text.y.left = element_text(size = 12),
    axis.title.x.bottom = element_text(size=10)
  )

tdat2 <- go_clst %>% 
  # 1. Remove grouping
  ungroup() %>% 
  # 2. Arrange by
  #   i.  facet group (tissue)
  #   ii. value (score)
  arrange(Fold_Enrichment, Term_Description) %>%
  # 3. Add order column of row numbers
  mutate(order = row_number())

# # tweak names for plot
# tdat2[1,]$Term_Description = "DNA-binding transcription activator activity,\n RNA polymerase II-specific"
# tdat2[4,]$Term_Description = "DNA-binding transcription factor activity,\n RNA polymerase II-specific"
# tdat2[9,]$Term_Description = "SCF-dependent proteasomal ubiquitin-dependent\n protein catabolic process"

pdf("/data/github/pca_network/results/pathfindR/go_top10.pdf", width=7, height=5)
ggplot(tdat2, aes(order, Fold_Enrichment)) +
  geom_point(aes(colour = lowest_p, size=size)) + 
  scale_x_continuous(
    breaks = tdat2$order,
    labels = tdat2$Term_Description,
    expand = c(0,0.5)) +
  # toggle 'expand' to pad axes, make room for points
  scale_y_continuous(expand = c(0, 0.35)) +
  #facet_grid(Cluster ~ ., scales = "free_y", space = "free", drop = TRUE) +
  coord_flip() +
  labs( title = "GO Pathways",
        subtitle = NULL,
        x=NULL, y="Fold Enrichment",
        color='P-value',size='Gene\nnumber'
  ) + 
  scale_color_gradientn(colours = rocket(100, begin = 0, end = 0.5, alpha = 0.8, direction = -1)) +
  theme_linedraw() + scale_size(range = c(2,7)) + theme(axis.text.y = element_text(face="bold"))
dev.off()


# KEGG 
kegg = kegg[which(kegg$occurrence>50),]
kegg_clst = cluster_enriched_terms(kegg, plot_dend = F, plot_clusters_graph = F, use_description = T, use_active_snw_genes = T)
save = c()
for(i in 1:nrow(kegg_clst)){
  row = as.data.frame(kegg_clst[i,])
  size = (length(strsplit(row$Up_regulated, ", ")[[1]]) + length(strsplit(row$Down_regulated, ", ")[[1]]))
  save = c(save, size)
}
kegg_clst$size = save

kegg_clst = kegg_clst[order(-kegg_clst$Fold_Enrichment),]
kegg_clst = kegg_clst[1:10,]
ggplot(kegg_clst, aes(y=Term_Description,x=Fold_Enrichment)) +
  geom_point(aes(color=lowest_p,size=size)) +
  scale_color_gradientn(colours = rocket(100, begin = 0, end = 0.5, alpha = 0.8, direction = -1)) +
  labs( title = "GO Pathways",
        subtitle = "GO-BP, GO-MF, GO-CC",
        x='Fold Enrichment', y=NULL,
        color='P-value',size='Gene\nnumber'
  ) + theme_linedraw() + scale_size(range = c(3, 10)) + ggplot2::facet_grid(kegg_clst$Cluster ~ ., scales = "free_y", space = "free", drop = TRUE)  +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size=12),
    axis.text.y.left = element_text(size = 12),
    axis.title.x.bottom = element_text(size=10)
  )
  
tdat2 <- kegg_clst %>% 
  # 1. Remove grouping
  ungroup() %>% 
  # 2. Arrange by
  #   i.  facet group (tissue)
  #   ii. value (score)
  arrange(Fold_Enrichment, Term_Description) %>%
  # 3. Add order column of row numbers
  mutate(order = row_number())

tdat2$Cluster = as.factor(tdat2$Cluster)

tdat2$order = rev(tdat2$order)

pdf("/data/github/pca_network/results/pathfindR/kegg_top10.pdf", width=8, height=4)
ggplot(tdat2, aes(order, Fold_Enrichment)) +
  geom_point(aes(colour = lowest_p, size=size)) + 
  scale_x_continuous(
    breaks = tdat2$order,
    labels = tdat2$Term_Description,
    expand = c(0,0.5)) +
  # toggle 'expand' to pad axes, make room for points
  scale_y_continuous(expand = c(0, 0.2)) +
  #facet_grid(Cluster ~ ., scales = "free_y", space = "free", drop = TRUE) +
  #coord_flip() +
  labs( title = "KEGG Pathways",
        subtitle = NULL,
        x=NULL, y="Fold Enrichment",
        color='P-value',size='Gene\nnumber'
  ) + 
  scale_color_gradientn(colours = rocket(100, begin = 0, end = 0.5, alpha = 0.8, direction = -1)) +
  theme_linedraw() + scale_size(range = c(4,6)) + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1, face="bold"))
dev.off()

# REACTOME
reactome = reactome[which(reactome$occurrence>50),]
reactome_clst = cluster_enriched_terms(reactome, plot_dend = F, plot_clusters_graph = F, use_description = T, use_active_snw_genes = T)
save = c()
for(i in 1:nrow(reactome_clst)){
  row = as.data.frame(reactome_clst[i,])
  size = (length(strsplit(row$Up_regulated, ", ")[[1]]) + length(strsplit(row$Down_regulated, ", ")[[1]]))
  save = c(save, size)
}
reactome_clst$size = save

reactome_clst = reactome_clst[order(-reactome_clst$Fold_Enrichment),]
reactome_clst = reactome_clst[1:10,]

ggplot(reactome_clst, aes(y=Term_Description,x=Fold_Enrichment)) +
  geom_point(aes(color=lowest_p,size=size)) +
  scale_color_gradientn(colours = rocket(100, begin = 0, end = 0.5, alpha = 0.8, direction = -1)) +
  labs( title = "GO Pathways",
        subtitle = "GO-BP, GO-MF, GO-CC",
        x='Fold Enrichment', y=NULL,
        color='P-value',size='Gene\nnumber'
  ) + theme_linedraw() + scale_size(range = c(3, 10)) + ggplot2::facet_grid(reactome_clst$Cluster ~ ., scales = "free_y", space = "free", drop = TRUE)  +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size=12),
    axis.text.y.left = element_text(size = 12),
    axis.title.x.bottom = element_text(size=10)
  )

tdat2 <- reactome_clst %>% 
  # 1. Remove grouping
  ungroup() %>% 
  # 2. Arrange by
  #   i.  facet group (tissue)
  #   ii. value (score)
  arrange(Fold_Enrichment, Term_Description) %>%
  # 3. Add order column of row numbers
  mutate(order = row_number())

tdat2$Cluster = as.factor(tdat2$Cluster)

pdf("/data/github/pca_network/results/pathfindR/reactome_top10.pdf", width=7, height=5)
ggplot(tdat2, aes(order, Fold_Enrichment)) +
  geom_point(aes(colour = lowest_p, size=size)) + 
  scale_x_continuous(
    breaks = tdat2$order,
    labels = tdat2$Term_Description,
    expand = c(0,0.8)) +
  # toggle 'expand' to pad axes, make room for points
  scale_y_continuous(expand = c(0, 0.5)) +
  #facet_grid(Cluster ~ ., scales = "free_y", space = "free", drop = TRUE) +
  coord_flip() +
  labs( title = "Reactome Pathways",
        subtitle = NULL,
        x=NULL, y="Fold Enrichment",
        color='P-value',size='Gene\nnumber'
  ) + 
  scale_color_gradientn(colours = rocket(100, begin = 0, end = 0.5, alpha = 0.8, direction = -1)) +
  theme_linedraw() + scale_size(range = c(3,6)) + theme(axis.text.y = element_text(face="bold"))
dev.off()




fetch_java_version <- function() {

  java_home <- Sys.getenv("JAVA_HOME", unset = NA)
  print(java_home)
  if (!is.na(java_home)) {
    java <- file.path(java_home, "bin", "java")
    if (identical(.Platform$OS.type, "windows")) {
      java <- paste0(java, ".exe")
    }
    if (!file.exists(java)) {
      print("executing file does not exist block")
      java <- ""
    }
  } else {
    java <- Sys.which("java")
  }
  print("printing java..")
  print(java)
  if (java == "") {
    stop(
      "Java version not detected. Please download and install Java from ",
      dQuote("https://www.java.com/en/")
    )
  }
  
  version <- system2(java, "-version", stderr = TRUE, stdout = TRUE)
  if (length(version) < 1) {
    stop(
      "Java version not detected. Please download and install Java from ",
      dQuote("https://www.java.com/en/")
    )
  }
  
  return(version)
}

fetch_java_version()
