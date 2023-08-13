#!/usr/bin/env Rscript

# troubles with newer java installation. See function below. 
Sys.setenv(JAVA_HOME= "/usr/lib/jvm/java-17-openjdk-amd64")
install.packages("pathfindR")
library(pathfindR)
library(dplyr)

network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
input = network[,c(5,6)]
input$adj.P.val = 0.001 # used for filtering >  0.05 in pathfindr
colnames(input) = c("Gene.Symbol", "logFC", "adj.P.Val")
input = input %>% unique()

keggpathway = pathfindR::run_pathfindR(input, gene_sets = "KEGG", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=4, iterations=1)
write.table(keggpathway, "/data/github/pca_network/results/kegg.txt", quote = F, row.names = F, sep="\t")


reactome = pathfindR::run_pathfindR(input, gene_sets = "Reactome", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=2, iterations=10)
write.table(reactome, "/data/github/pca_network/results/reactome.txt", quote = F, row.names = F, sep="\t")

go = pathfindR::run_pathfindR(input, gene_sets = "GO-All", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=2, iterations=5)
write.table(go, "/data/github/pca_network/results/go.txt", quote = F, row.names = F, sep="\t")

Sys.setenv(JAVA_HOME= "/usr/lib/jvm/java-17-openjdk-amd64")

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
