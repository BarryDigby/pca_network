#!/usr/bin/env Rscript
library(GEOquery)
library(dplyr)

gset <- getGEO("GSE116918", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL25318", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

belfast_probes = read.csv("/data/github/pca_network/data/belfast_signf_probes.txt", header=F, sep="\t")
colnames(belfast_probes) = c("probe", "ensembl_gene_id")
ex = ex[which(rownames(ex) %in% belfast_probes$probe),]
ex = as.data.frame(ex)
ex = merge(ex, belfast_probes, by.x=0, by.y="probe")
ex = ex[,2:ncol(ex)]
#ex = tibble::column_to_rownames(ex, "ensembl_gene_id")

# take the highest available expr value for each representative ENSG probe.
ex <- ex %>%
  group_by(ensembl_gene_id) %>%
  summarise_all(.funs = max) %>%
  ungroup()

# load Active.Genes
load("/data/github/pca_network/results/LASSO_cox.RData")
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% Active.Genes),]
mrna_attributes = subset(mrna_attributes, select=c(external_gene_name, ensembl_gene_id))

ex = merge(ex, mrna_attributes, by="ensembl_gene_id")
ex = tibble::column_to_rownames(ex, "external_gene_name")
ex = ex[,2:ncol(ex)]
mat = ex

belfast_rds = readRDS("/data/github/pca_network/data/Belfast_eSet.RDS")
mat2 = belfast_rds@assayData$exprs
meta = belfast_rds@phenoData@data
meta$days_to_metastasis = meta$time_to_metastasis * 30.44
rm(belfast_rds)

# Scale center log2 counts 
mat = t(scale(t(mat), scale = T, center = T))
mat = as.data.frame(t(mat))

# multiply coefficients
risk_score = rowSums( Active.Coefficients*mat )
mat$risk_score = risk_score
sub_meta = select(meta, select=c(days_to_metastasis, metastasis_status))
mat = cbind(mat, sub_meta)
colnames(mat)[11:12] <- c("days_to_metastasis", "metastasis_status")

mat$bcr_status <- meta$bcr_status

# COXPH

library(survival)
library(survminer)
mat$risk_score_cat = ifelse(mat$risk_score > mean(mat$risk_score), "high", "low")
surv_object <- Surv(mat$days_to_metastasis, mat$metastasis_status)
res2 = survfit(surv_object ~ risk_score_cat, data=mat)

ggsurvplot(res2,
           pval = TRUE, conf.int = T,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=mat,
           xlab="Time (days)")

surv_object <- Surv(mat$days_to_metastasis, mat$bcr_status)
res2 = survfit(surv_object ~ risk_score_cat, data=mat)

ggsurvplot(res2,
           pval = TRUE, conf.int = T,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=mat,
           xlab="Time (days)")



