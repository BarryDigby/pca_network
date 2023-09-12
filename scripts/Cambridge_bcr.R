#!/usr/bin/env Rscript 

library(survival)
library(survminer)

cambridge = readRDS("/data/github/pca_network/data/Cambridge_eSet.RDS")
mat = cambridge@assayData$exprs
meta = cambridge@phenoData@data
meta = meta[which(!is.na(meta$bcr_status)),]
mat = as.data.frame(mat[,rownames(meta)])
meta$days_to_follow_up = floor(meta$time_to_bcr*30.44)

#load("/data/github/pca_network/results/prognostic_model_os2.RData")

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% Active.Genes),]
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat = tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

mat = as.data.frame(scale(mat, scale = F, center = T))

risk = rowSums(mat*Active.Coefficients)
mat$risk_score = risk
mat$risk_category = ifelse(mat$risk_score > median(mat$risk_score), "high", "low")
mat = cbind(mat, subset(meta, select=c(days_to_follow_up, bcr_status)))

surv_object <- Surv(mat$days_to_follow_up, mat$bcr_status)
cox = coxph(surv_object ~ risk_category, data=mat)
res = survfit(surv_object ~ risk_category, data=mat)
logrank = survdiff(surv_object ~ risk_category, data=mat)

pdf("/data/github/pca_network/results/TCGA_DFS3/Cambridge_DFS.pdf", height=8,width=8)
ggsurvplot(res,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=mat,
           xlab="Time (days)")
dev.off()

ROC <- timeROC(T=mat$days_to_follow_up,
                     delta=mat$bcr_status,
                     marker=mat$risk_score,
                     cause=1,weighting="marginal",
                     times=c(365, floor(365*3), floor(365*5), floor(365*10)),
                     iid=TRUE)
ROC
confint(ROC, level = 0.95)
pdf("/data/github/pca_network/results/TCGA_DFS3/Cambridge_TimeROC.pdf")
plot(ROC, time=365)
plot(ROC, time=365*3)
plot(ROC, time=365*5)
dev.off()