#!/usr/bin/env Rscript

library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)

mirna_query <- GDCquery(project = "TCGA-PRAD",
                        data.category = "Transcriptome Profiling",
                        data.type = "miRNA Expression Quantification",
                        #workflow.type = "BCGSC miRNA Profiling",
                        experimental.strategy = "miRNA-Seq")

GDCdownload(mirna_query, method = "api", files.per.chunk = 100,
            directory = "~/Desktop/TCGA/miRNA/")

miR_df <- GDCprepare(mirna_query, directory = "~/Desktop/TCGA/miRNA/")

## remove columns we dont need, keep counts
rownames(miR_df) <- miR_df$miRNA_ID
miR_df <- miR_df[,-1]
number_cols <- ncol(miR_df)
subset <- seq(from = 1, to = number_cols, by = 3)
miR_df <- miR_df[, subset]

## Strip read_count, just want the 'cases' ID
colnames(miR_df) <- gsub(".*_","",colnames(miR_df))

## Match to metadata
miR_meta <- mirna_query[[1]][[1]]

miR_meta <- miR_meta[,c("cases", "sample_type")]
rownames(miR_meta) <- colnames(miR_df)
table(rownames(miR_meta) == miR_meta$cases)

## fix the levels that R thinks are there but are not
miR_meta$sample_type <- as.character(miR_meta$sample_type)
table(miR_meta$sample_type)

## Remove metastatic sample
metastatic_key <- miR_meta[which(miR_meta$sample_type == "Metastatic"),]

miR_meta <- miR_meta[!miR_meta$sample_type == metastatic_key$sample_type,]
miR_df <- miR_df[, -grep(paste0(metastatic_key$cases), colnames(miR_df))]

## Rename conditions
miR_meta$sample_type <- gsub("Primary solid Tumor", "Tumor", miR_meta$sample_type)
miR_meta$sample_type <- gsub("Solid Tissue Normal", "Normal", miR_meta$sample_type)
miR_meta$sample_type <- as.factor(miR_meta$sample_type)
levels(miR_meta$sample_type)
colnames(miR_meta) <- c("cases", "Condition")

## tidy vars
rm(mirna_query)
rm(subset)
rm(number_cols)
rm(metastatic_key)

## Limma etc. 
design = model.matrix( ~ 0 + Condition, data = miR_meta )                 

y <- edgeR::DGEList(miR_df)
keep <- edgeR::filterByExpr(y, design)
y <- y[keep, ]

# normalize and run voom transformation
y <- edgeR::calcNormFactors(y)
logcpm <- edgeR::cpm(y, normalized.lib.sizes = T, log=TRUE)
scaled <- t(scale(t(logcpm)))
v <- limma::voom(y, design, plot = F)



PCA <- prcomp(t(scaled), scale = F)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Status = miR_meta$Condition)

library(ggplot2)
pdf("/data/github/pca_network/mirna/TCGA/TCGA_mirna_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Status", palette = c("forestgreen", "orangered"),
                  title = "Log-transformed normalized expression data\n[TCGA miRNA]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

colnames(design) = c("normal","tumor")

contrast <- limma::makeContrasts(
  t_v_n = tumor - normal,
  levels = colnames(design))

fit = limma::lmFit(v, design = design)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design = design), contrast))

t_v_n = limma::topTable(fit, number=Inf, p.value = 0.05, adjust.method = "BH", coef="t_v_n")
t_v_n = tibble::rownames_to_column(t_v_n, "miRNA")
t_v_n$miRNA = gsub("mir", "miR", t_v_n$miRNA)

# mir IDs seem to be missing 5p , 3p etc etc. Annotated with mirbase aliases to make consistent with other datasets.

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
t_v_n$new_id = update_id(t_v_n)
t_v_n$miRNA = t_v_n$new_id
t_v_n = t_v_n[,c(1:7)]



write.table(t_v_n, "/data/github/pca_network/mirna/TCGA/tumor_v_normal.txt", sep="\t", row.names = F, quote = F)

