#!/usr/bin/env Rscript 

## Create 'expanded' network from prognostic signature.
## i.e for all circRNAs in the prognstic network, retrieve their parent genes. 
## unsure if we should expand at all sites (circ, miR, gene). . . 

prog_net = read.csv("/data/github/pca_network/results/prognostic_edgelist.csv", sep="\t", header=T)
prog_circ = prog_net$source[grep("hsa_circ", prog_net$source)]
prog_circ = unique(prog_circ)


## identify parent genes
limma_circs = read.csv("/data/github/pca_network/results/circrna_intersection.txt", sep="\t", header=T)
limma_circs = limma_circs[which(limma_circs$circbaseID %in% prog_circ),]
circ_parent = limma_circs[, c("circbaseID", "GeneSymbol")]

## Load EXPR data
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(limma)

mrna_query <- GDCquery(project = "TCGA-PRAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       experimental.strategy = "RNA-Seq")

mrna_df <- GDCprepare(mrna_query, directory = "~/Desktop/TCGA/mRNA/")
mrna_df <- assay(mrna_df)

## tidy matrix colnames 
delim_fn = function(x, n, i){
  do.call(c, lapply(x, function(X)
    paste(unlist(strsplit(X, "-"))[(n+1):(i)], collapse = "-")))
}

colnames(mrna_df) <- delim_fn(x = colnames(mrna_df), n = 0, i = 4)
mrna_df <- as.data.frame(mrna_df)

load("/data/github/pca_network/results/TCGA_DFS/model.RData")
mrna_attributes = read.csv("/data/github/pca_network/mrna/TCGA/ENS_BIO_HGNC.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$biotype == "protein_coding"),]
tcga_df = merge(mrna_df, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_version_id")
rownames(tcga_df) <- make.names(tcga_df$Gene, unique = TRUE)
#tcga_df = tibble::column_to_rownames(tcga_df, "Gene")
tcga_df = tcga_df[,2:ncol(tcga_df)]
tcga_df = tcga_df[, rownames(cox)]


tcga_meta = data.frame(row.names = colnames(tcga_df),
                       Risk_strata = cox$risk_category)

tcga_meta$Risk_strata = ifelse(tcga_meta$Risk_strata == "High risk", "High", "Low")
## Limma etc. 
design = model.matrix( ~ 0 + Risk_strata, data = tcga_meta )      

y <- edgeR::DGEList(tcga_df)
#keep <- edgeR::filterByExpr(y, design)
#y <- y[keep, ]
logcpm <- edgeR::cpm(y, normalized.lib.sizes = T, log=TRUE)
y <- edgeR::calcNormFactors(y)
v <- limma::voom(y, design, plot = F)



mrna_mat = logcpm[which(rownames(logcpm) %in% circ_parent$GeneSymbol), ]
mrna_mat = as.data.frame(t(mrna_mat))
mrna_mat = cbind(mrna_mat, cox$risk_category)


wide_df = mrna_mat %>%
  pivot_longer(cols = colnames(mrna_mat[1:ncol(mrna_mat)-1]), names_to = "Variable", values_to = "Score")

colnames(wide_df)[1] = "strata"

library(ggpubr)

ggboxplot(wide_df, x="strata", y="Score", facet.by = "Variable", color="black", fill="strata", add = "mean_se") + ylim(0, 15) + stat_compare_means(method="t.test", label.y = 14)




contrast <- limma::makeContrasts(
  high_v_low = Risk_strataHigh - Risk_strataLow,
  levels = colnames(design))

fit = limma::lmFit(v, design = design)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design = design), contrast))

h_v_l = limma::topTable(fit, number=Inf, p.value = 1, adjust.method = "BH", coef="high_v_low")
h_v_l$Gene = rownames(h_v_l)
table(h_v_l$Gene %in% circ_parent$GeneSymbol)

h_v_l = h_v_l[which(h_v_l$Gene %in% circ_parent$GeneSymbol),]

library(pathfindR)

input = h_v_l[,c("Gene", "logFC", "adj.P.Val")]
input$adj.P.Val = 0.001 # d not want pathfindr to filter and remove by P 
colnames(input) = c("Gene.Symbol", "logFC", "adj.P.Val")
input = input %>% unique()

# Run analysis and save. load for future sessions
keggpathway = pathfindR::run_pathfindR(input, gene_sets = "KEGG", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=4, iterations=100)
write.table(keggpathway, "/data/github/pca_network/results/parent_gene_kegg.txt", quote = F, row.names = F, sep="\t")

reactome = pathfindR::run_pathfindR(input, gene_sets = "Reactome", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=4, iterations=100)
write.table(reactome, "/data/github/pca_network/results/parent_gene_reactome.txt", quote = F, row.names = F, sep="\t")

go = pathfindR::run_pathfindR(input, gene_sets = "GO-All", plot_enrichment_chart = F, list_active_snw_genes = T, n_processes=4, iterations=100)
write.table(go, "/data/github/pca_network/results/parent_gene_go.txt", quote = F, row.names = F, sep="\t")

master_pathfindr = rbind(keggpathway, reactome, go)
clst_search = cluster_enriched_terms(master_pathfindr, plot_dend = F, plot_clusters_graph = F, use_description = T, use_active_snw_genes = T)
clst_search = clst_search[which(clst_search$occurrence>50), ]

sets <- msigdbr::msigdbr(species="Homo sapiens", category = "C2", subcategory = "CP:KEGG")
sets1 <-  msigdbr::msigdbr(species="Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
sets2 <- msigdbr::msigdbr(species="Homo sapiens", category = "C5", subcategory = "GO:BP")
sets3 <- msigdbr::msigdbr(species="Homo sapiens", category = "C5", subcategory = "GO:MF")
sets4 <- msigdbr::msigdbr(species="Homo sapiens", category = "C5", subcategory = "GO:CC")
master = rbind(sets, sets1, sets2, sets3, sets4)
# use all sets returned by both pathway analysis
my_sets = master[which(master$gs_exact_source %in% clst_search$ID),]
ids = unique(my_sets$gs_exact_source)
my_sets = split(
  my_sets$gene_symbol, # The genes we want split into pathways
  my_sets$gs_name # The pathways made as the higher levels of the list
)

var = GSVA::gsva(logcpm, my_sets, method = "gsva", kcdf = "Gaussian",  min.sz = 1, max.sz = 500, mx.diff = TRUE, verbose = FALSE)
var_df = as.data.frame(var)
out = c()
for(i in 1:nrow(var_df)){
  row = as.matrix(t(var_df[i,]))
  test = t.test(row ~ tcga_meta$Risk_strata, paired=F)
  high = test$estimate[1]
  low = test$estimate[2]
  pathway = colnames(row)
  p = test$p.value
  if( high > low & p < 0.05){
    out = c(out, pathway)
  }
}

var_df = var_df[out, ]
rownames(var_df) = sub("^[^_]*_", "", rownames(var_df))
ann_col = data.frame(row.names = colnames(var_df),
                     Risk = tcga_meta$Risk_strata)

col <- c("forestgreen", "purple")
names(col) <- c("Low", "High")
ann_clr <- list(Risk = col)

var_df = t(scale(t(var_df), center = T))
remove = "SENSORY_PROCESSING_OF_SOUND_BY_OUTER_HAIR_CELLS_OF_THE_COCHLEA"
var_df = var_df[which(!rownames(var_df) %in% remove),]
pheatmap(var_df,
         cluster_rows = T, 
         cluster_cols = F,
         show_colnames = FALSE,
         scale = "row",
         annotation_col=ann_col,
         annotation_colors = ann_clr,
         color = hcl.colors(100, "RdBu",rev=T))

# heatmap of parent genes
hm_mat = logcpm[input$Gene.Symbol, ]
tcga_meta$Risk_strata= sort(tcga_meta$Risk_strata, decreasing = TRUE)

ann_col = data.frame(row.names = rownames(tcga_meta),
                     Risk_strata = tcga_meta$Risk_strata)

col <- c("royalblue", "red2")
names(col) <- c("Low", "High")
ann_clr <- list(Risk_strata = col)

hm_mat = t(scale(t(hm_mat), scale = T, center = T))

pheatmap(hm_mat,
         cluster_rows = T, 
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "row",
         annotation_col=ann_col,
         annotation_colors = ann_clr,
         color = hcl.colors(100, "RdBu",rev=T))


# take a look at the prognostic value of the miRNAs in the network. 
mirs = unique(prog_net$source[grep("hsa-miR", prog_net$source)])

mirna_query <- GDCquery(project = "TCGA-PRAD",
                        data.category = "Transcriptome Profiling",
                        data.type = "miRNA Expression Quantification",
                        #workflow.type = "BCGSC miRNA Profiling",
                        experimental.strategy = "miRNA-Seq")

#GDCdownload(mirna_query, method = "api", files.per.chunk = 100,
#           directory = "~/Desktop/TCGA/miRNA/")

miR_df <- GDCprepare(mirna_query, directory = "~/Desktop/TCGA/miRNA/")

## remove columns we dont need, keep counts
rownames(miR_df) <- miR_df$miRNA_ID
miR_df <- miR_df[,-1]
number_cols <- ncol(miR_df)
subset <- seq(from = 1, to = number_cols, by = 3)
miR_df <- miR_df[, subset]

## Strip read_count, just want the 'cases' ID
colnames(miR_df) <- gsub(".*_","",colnames(miR_df))

## Limma etc. 
y <- edgeR::DGEList(miR_df)
y <- edgeR::calcNormFactors(y)
logcpm <- edgeR::cpm(y, normalized.lib.sizes = T, log=TRUE)
logcpm = as.data.frame(logcpm)

## match columns 
colnames(logcpm) =  substr(colnames(logcpm), 1, 16)

rownames(logcpm) = gsub("mir", "miR", rownames(logcpm))

logcpm$miRNA = rownames(logcpm)

alias = read.csv("/data/github/pca_network/data/mirbase_aliases_tabbed.txt", sep="\t", header=F)


## Use updated -3p/-5p notation instead of *
update_id = function(df) {
  out_vec = c()
  for (index in 1:length(df$miRNA)) {
    ## stage miRNA from tt
    mir = df$miRNA[index]
    ## if its end in * add \\* for str_detect REGEX
    if (grepl("\\*", mir)) {
      mir <- gsub("\\*", "\\\\*", mir)
    }
    # grab the corresponding row in Alias
    row = alias %>% filter_all(any_vars(str_detect(., paste0("\\b", gsub("\\*", "\\\\\\*", mir), "\\b"))))
    if (nrow(row) == 0) {
      # no match with alias, keep original label as it is valid.
      out_vec = c(out_vec, mir)
    } else{
      # Vectorise these.. so you can do n + 1
      if (nrow(row) == 1) {
        # match does not like escaping chars
        row_vec = unname(unlist(row[1, ]))
        idx = match(gsub("\\\\", "", mir), row_vec)
        idx = idx + 1
        out_st = row_vec[idx]
        # if it now contains a *, move to next 'latest version'. this is rare. 
        if(grepl("\\*", out_st)){out_st = row_vec[idx+1]}
        # sometimes the matching string is correct, but also matches alias in V3 (no adjacent text)
        if(out_st == "" || is.na(out_st)){out_st=gsub("\\\\", "", mir)}
        out_vec = c(out_vec, out_st)
      } else if (nrow(row)==2){
        row_vec =  c(unname(unlist(row[1, ])), unname(unlist(row[2, ])))
        idx = match(gsub("\\\\", "", mir), row_vec)
        idx = idx + 1
        out_st = row_vec[idx]
        # if it now contains a *, move to next 'latest version'. this is rare. 
        if(grepl("\\*", out_st)){out_st = row_vec[idx+1]}
        # sometimes the matching string is correct, but also matches alias in V3 (no adjacent text)
        if(out_st == "" || is.na(out_st)){out_st=gsub("\\\\", "", mir)}
        out_vec = c(out_vec, out_st)
      } else {
        # edeg cases that draw three rows
        row_vec =  c(unname(unlist(row[1, ])), unname(unlist(row[2, ])), unname(unlist(row[3, ])))
        idx = match(gsub("\\\\", "", mir), row_vec)
        idx = idx + 1
        out_st = row_vec[idx]
        # if it now contains a *, move to next 'latest version'. this is rare. 
        if(grepl("\\*", out_st)){out_st = row_vec[idx+1]}
        # sometimes the matching string is correct, but also matches alias in V3 (no adjacent text)
        if(out_st == "" || is.na(out_st)){out_st=gsub("\\\\", "", mir)}
        out_vec = c(out_vec, out_st)
      }
    }
  }
  return(out_vec)
}

library(dplyr)
library(stringr)
logcpm$new_id = update_id(logcpm)

## some manual mapping required here, e.g hsa-miR-26a-5p comes from hsa-miR-26a-1
mirs_update = c("hsa-miR-26a-1", "hsa-miR-30c-1", "hsa-miR-106a", "hsa-miR-1976", "hsa-miR-335-5p")
rownames(logcpm) = logcpm$new_id
sub = colnames(logcpm)[colnames(logcpm) %in% rownames(cox)]
logcpm = logcpm[mirs_update, sub]

mir_meta = data.frame(row.names = rownames(cox),
                      risk = cox$risk_category,
                      time = cox$days_to_follow_up,
                      bcr = cox$bcr_status)

mir_meta$risk = ifelse(mir_meta$risk == "High risk", "High", "Low")

logcpm = as.data.frame(t(logcpm))

logcpm = merge(logcpm, mir_meta, by=0)
logcpm = tibble::column_to_rownames(logcpm, "Row.names")
colnames(logcpm)[1:5] = gsub("\\-", "_", colnames(logcpm)[1:5])
library(RegParallel)
res <- RegParallel(
  data = logcpm,
  formula = 'Surv(time, bcr) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(logcpm)[1:5],
  blocksize = 1,
  p.adjust = "BH",
  cores = 1,
  nestedParallel = FALSE,
  conflevel = 95)

survminer::surv_cutpoint(logcpm, time = "time", event = "bcr", variables = colnames(logcpm)[1:5])


pheatmap(t(logcpm[,1:5]), scale = "row")
