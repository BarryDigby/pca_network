#!/usr/bin/env Rscript

library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(limma)

mrna_query <- GDCquery(project = "TCGA-PRAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       experimental.strategy = "RNA-Seq")

GDCdownload(mrna_query, method = "api", files.per.chunk = 100,
            directory = "~/Desktop/TCGA/mRNA/")

mrna_df <- GDCprepare(mrna_query, directory = "~/Desktop/TCGA/mRNA/")

# metadata full
list_columns <- sapply(mrna_df@colData@listData, is.list)
# Step 2: Exclude list columns from the data
mrna_df@colData@listData <- mrna_df@colData@listData[!list_columns]
# Step 3: Convert the modified data to a dataframe
colData_df <- as.data.frame(mrna_df@colData@listData)

mrna_meta <- mrna_df$sample
mrna_meta <- cbind(mrna_meta, mrna_df$definition)
mrna_df <- assay(mrna_df)

## tidy matrix colnames 
delim_fn = function(x, n, i){
  do.call(c, lapply(x, function(X)
    paste(unlist(strsplit(X, "-"))[(n+1):(i)], collapse = "-")))
}

colnames(mrna_df) <- delim_fn(x = colnames(mrna_df), n = 0, i = 4)

mrna_meta <- as.data.frame(mrna_meta)
mrna_df <- as.data.frame(mrna_df)

## remove the metastatic sample from counts matrix and metadata 
metastatic_key <- mrna_meta[which(mrna_meta[,2] == "Metastatic"),]
mrna_meta <- mrna_meta[!mrna_meta[,2] == metastatic_key[,2],]
mrna_df <- mrna_df[, -grep(paste0(metastatic_key[,1]), colnames(mrna_df))]

## fix the levels that R thinks are there but are not
mrna_meta[,2] <- as.character(mrna_meta[,2])

## Rename conditions
mrna_meta[,2] <- gsub("Primary solid Tumor", "Tumor", mrna_meta[,2])
mrna_meta[,2] <- gsub("Solid Tissue Normal", "Normal", mrna_meta[,2])
mrna_meta[,2] <- as.factor(mrna_meta[,2])
levels(mrna_meta[,2])
colnames(mrna_meta) <- c("cases", "Condition")


## filter for protein coding genes in matrix (currently > 50,000 rows)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mrna_attributes <- getBM(attributes=c("external_gene_name",
                                      "ensembl_gene_id_version",
                                      "ensembl_gene_id",
                                      "gene_biotype"),
                         filters = c("ensembl_gene_id_version"),
                         values = rownames(mrna_df),
                         mart = mart,
                         useCache = T)

mrna_attributes <- mrna_attributes[which(mrna_attributes$gene_biotype == "protein_coding"),]

# annotate after DE.

## Limma etc. 
design = model.matrix( ~ 0 + Condition, data = mrna_meta )                 

y <- edgeR::DGEList(mrna_df)
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
                     Status = mrna_meta$Condition)

library(ggplot2)
pdf("/data/github/pca_network/mrna/TCGA/TCGA_mrna_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Status", palette = c("forestgreen", "orangered"),
                  title = "Log-transformed normalized expression data\n[TCGA mRNA]",
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
t_v_n = tibble::rownames_to_column(t_v_n, "Gene")
t_v_n = merge(t_v_n, mrna_attributes, by.x="Gene", by.y="ensembl_gene_id_version", all.x=T)
t_v_n$Gene = t_v_n$external_gene_name
t_v_n = t_v_n[which(t_v_n$Gene != ""),]
t_v_n = t_v_n[,c(1:7)]

write.table(t_v_n, "/data/github/pca_network/mrna/TCGA/tumor_vs_normal.txt", sep="\t", quote=F, row.names = F)