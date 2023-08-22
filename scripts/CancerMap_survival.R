#!/usr/bin/env Rscript

cmap = readRDS("/data/github/pca_network/data/CancerMap_eSet.RDS")
mat = cmap@assayData$exprs
meta = cmap@phenoData@data

# remove normal samples from dataset
meta = meta[which(meta$sample_type!="Normal"),]
mat = mat[,which(colnames(mat) %in% rownames(meta))]

# data is normalised, on log2 scale
boxplot(mat[,1:20])

load("/data/github/pca_network/results/LASSO_cox.RData")
load("/data/github/pca_network/results/TCGA_survival.RData")

signf = signf_os[which(signf_os$gene %in% Active.Genes),]

mat = mat[which(rownames(mat) %in% signf$ensembl),]
signf = subset(signf, select=c(gene, ensembl))
mat = merge(mat, signf, by.x=0, by.y="ensembl")
mat = tibble::column_to_rownames(mat, "gene")
mat = mat[,2:ncol(mat)]

mat = as.data.frame(t(scale(t(mat), scale = T, center = T)))
boxplot(mat[,1:20])

mat = as.data.frame(t(mat))

risk = rowSums(Active.Coefficients*mat)
mat$risk_score = risk
sub_meta = subset(meta, select=c(time_to_bcr, bcr_status))
mat = cbind(mat, sub_meta)
mat$days_to_bcr = mat$time_to_bcr * 30.44
mat$risk_score_cat = ifelse(mat$risk_score > median(mat$risk_score), "high", "low")

surv_object <- Surv(mat$days_to_bcr, mat$bcr_status)
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
