#!/usr/bin/env Rscript
library(GEOquery)
library(dplyr)


belfast_rds = readRDS("/data/github/pca_network/data/Belfast_eSet.RDS")
mat = belfast_rds@assayData$exprs
meta = belfast_rds@phenoData@data
meta$days_to_follow_up = meta$time_to_bcr * 30.44
rm(belfast_rds)

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% Active.Genes),]
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat = tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

mat = as.data.frame(scale(mat, scale = F, center = T))

# multiply coefficients
risk_score = rowSums( Active.Coefficients*mat )
mat$risk_score = risk_score
sub_meta = select(meta, select=c(days_to_follow_up, bcr_status))
mat = cbind(mat, sub_meta)
colnames(mat)[4:5] <- c("days_to_follow_up", "bcr_status")


# COXPH

library(survival)
library(survminer)
mat$risk_category = ifelse(mat$risk_score > median(mat$risk_score), "high", "low")
surv_object <- Surv(mat$days_to_follow_up, mat$bcr_status)
res2 = survfit(surv_object ~ risk_category, data=mat)



ggsurvplot(res2,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=mat,
           xlab="Time (days)")



ROC <- timeROC(T=mat$days_to_follow_up,
                     delta=mat$bcr_status,
                     marker=mat$risk_score,
                     cause=1,weighting="marginal",
                     times=c(365, floor(365*3), floor(365*5), floor(365*10)),
                     iid=TRUE)
ROC
confint(ROC.train, level = 0.95)
pdf("/data/github/pca_network/results/TCGA_DFS3/Belfast_TimeROC.pdf")
plot(ROC.train, time=365)
plot(ROC.train, time=365*3)
plot(ROC.train, time=365*5)
dev.off()


