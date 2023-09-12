#!/usr/bin/env Rscript 

library(dplyr)
library(tidyr)
library(survminer)
library(survival)
library(survivalAnalysis)
library(pROC)
library(survivalROC)
library(ggpubr)
library(gtsummary) # format cox model 
library(RColorBrewer)

#########################################################################
# Use this dataframe to keep samples consistent
#########################################################################
atlas_mat = read.csv("/data/github/pca_network/data/prad_rna_cancer_sample.tsv", header=T, sep="\t")
atlas_mat = atlas_mat[,c(1,2,4)]
atlas_mat <- atlas_mat %>%
  pivot_wider(names_from = Sample, values_from = FPKM, values_fill = 0)
atlas_mat = tibble::column_to_rownames(atlas_mat, "Gene")

#########################################################################
# Load atlas metadata
# Remove missing OS status
#########################################################################
atlas_meta = read.csv("/data/github/pca_network/data/tcga_updated_meta.csv", header=T, sep=",")
rownames(atlas_meta) = atlas_meta$sample
atlas_meta = atlas_meta[which(atlas_meta$sample %in% colnames(atlas_mat)),]
rem = !(atlas_meta$vital=="Not Reported")
atlas_meta = atlas_meta[rem,]
atlas_mat = atlas_mat[,rem]
atlas_meta$status = ifelse(atlas_meta$vital=="Dead", 1, 0)

#########################################################################
# Convert ENS to symbols and subset matrix by Active.Genes,
# order correctly (FPKM)
# FPKM matrix now ready for expression plots
#########################################################################
load("/data/github/pca_network/results/prognostic_model_os2.RData")
ensv109 = read.csv("/data/github/pca_network/data/ensembl_v109_proteinatlas.csv", sep="\t", header = F)
colnames(ensv109) = c("ensembl_gene_id", "biotype", "hgnc_symbol")
ensv109 = ensv109[which(ensv109$hgnc_symbol %in% Active.Genes),]
atlas_mat = atlas_mat[which(rownames(atlas_mat) %in% ensv109$ensembl_gene_id),]
atlas_mat = merge(atlas_mat, ensv109[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
atlas_mat = tibble::column_to_rownames(atlas_mat, "hgnc_symbol")
atlas_mat = atlas_mat[,c(2:ncol(atlas_mat))]
atlas_mat = as.data.frame(t(atlas_mat))
atlas_mat = atlas_mat[rownames(atlas_meta),]

#########################################################################
# Use scaled centered logcpm STAR counts for CoxPH/Kaplan Meier
# coefficients came from this data - must be used on this again
#########################################################################
scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
colnames(scaled) = gsub("\\.", "-", colnames(scaled))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% Active.Genes),]
scaled = merge(scaled, mrna_attributes[,1:2], by.x=0, by.y="ensembl_gene_id_version")
scaled = tibble::column_to_rownames(scaled, "external_gene_name")
scaled = scaled[,c(2:ncol(scaled))]
scaled = scaled[,rownames(atlas_meta)]
scaled = as.data.frame(t(scaled))

#########################################################################
# Load TCGA RDS file and assess clinical meta vs atlas meta
#########################################################################
tcga_meta = readRDS("/data/github/pca_network/data/TCGA-PRAD_eSet.RDS")
tcga_meta = tcga_meta@phenoData@data

#########################################################################
# we are interested in Biochemical recurrence BCR
# match sample ids to atlas before merging..
# days to follow up is the same as time to BCR in TCGA meta
#########################################################################
rownames(tcga_meta) = paste(rownames(tcga_meta), "A", sep="")
tcga_meta = subset(tcga_meta, select=c(bcr_status))
atlas_meta = merge(atlas_meta, tcga_meta, by=0)
rownames(atlas_meta) = atlas_meta$Row.names
scaled = scaled[atlas_meta$Row.names,]

#########################################################################
# Multiply gene * coefficients = risk score
# append time to event, DFS status
#########################################################################
dfs_mat = scaled
risk = rowSums(Active.Coefficients*dfs_mat)
dfs_mat$risk_score = risk
dfs_mat$risk_category = ifelse(dfs_mat$risk_score > median(dfs_mat$risk_score), "high", "low")
dfs_mat = cbind(dfs_mat, subset(atlas_meta, select=c(days_to_follow_up, bcr_status)))

#########################################################################
# Kaplan-Meier DFS TCGA
#########################################################################
surv_object <- Surv(dfs_mat$days_to_follow_up, dfs_mat$bcr_status)
cox = coxph(surv_object ~ risk_category, data=dfs_mat)
res = survfit(surv_object ~ risk_category, data=dfs_mat)
logrank = survdiff(surv_object ~ risk_category, data=dfs_mat)

pdf("/data/github/pca_network/results/TCGA_DFS2/Risk_category_scaled_DFS.pdf", height=8,width=8)
ggsurvplot(res,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=dfs_mat,
           xlab="Time (days)")
dev.off()

#########################################################################
# Univariate Cox proportional hazards regression:
# Find LASSO Cox mRNAs ** with prognosis of PCa DFS
#########################################################################
options(scipen = 999)
result = data.frame(gene = character(),
                    best_cutoff = numeric(),
                    pvalue = numeric())
for(i in Active.Genes){
  mat = as.data.frame(atlas_mat[rownames(atlas_meta),which(colnames(atlas_mat)==i)])
  rownames(mat) = rownames(atlas_meta)
  colnames(mat) = i
  mat = merge(mat, atlas_meta, by=0)
  best_thresh = roc(mat$bcr_status, mat[,2])
  best_thresh = coords(best_thresh, x="best")
  mat$high_low = ifelse(mat[,2] > best_thresh[1,1], "High", "Low")
  res.cox = survfit(Surv(days_to_follow_up, bcr_status) ~ high_low, data=mat)
  res.cox.pval = surv_pvalue(res.cox)
  pval = res.cox.pval$pval
  if( pval < 0.05 ){
    row = data.frame(gene=i, best_cutoff=best_thresh[1,1], pvalue=pval)
    result = rbind(result,row)
  }
}

#########################################################################
# Scatterhist DFS FPKM - do all for supplementary
# consult result above for manuscript ** ones 
#########################################################################

signf_dfs = signf[which(signf$gene %in% Active.Genes),]
for(i in 1:nrow(signf_dfs)){
  row = as.data.frame(signf_dfs[i,])
  gene = row$gene
  
  mat = atlas_mat[rownames(atlas_meta),]
  mat = merge(mat, atlas_meta, by=0)
  mat$bcr_status = factor(mat$bcr_status)
  mat$years_to_follow_up = mat$days_to_follow_up/365.25
  best_thresh = roc(mat$bcr_status, mat[,gene])
  best_thresh = coords(best_thresh, x="best")
  mat$high_low = ifelse(mat[,gene] > best_thresh[1,1], "High", "Low")
  
  p <- mat %>% arrange(bcr_status) %>% ggscatterhist(mat,  x=paste0(gene), y="years_to_follow_up", palette = c("royalblue3","red1"),
                                                ylab = "Time after diagnosis (years)", xlab = "Expression level (FPKM)", fill = "bcr_status",
                                                color="bcr_status", shape="bcr_status", alpha = 0.9, ggtheme = theme_bw(), size = 2,
                                                margin.params = list(fill="bcr_status"), margin.plot = "density", legend = "top")
  
  p$sp <- p$sp + geom_vline(xintercept = best_thresh[1,1], linetype = "dashed", color = "black")
  
  pdf(paste0("/data/github/pca_network/results/TCGA_DFS2/",gene,"_fpkm_scatter.pdf"), height=5, width=8)
  print(p)
  dev.off()
}

#########################################################################
# Kaplan MeierDFS FPKM bcr 
#########################################################################

for(i in 1:nrow(signf_dfs)){
  
  row = as.data.frame(signf_dfs[i,])
  gene = row$gene
  
  mat = atlas_mat[rownames(atlas_meta),]
  mat = merge(mat, atlas_meta, by=0)
  best_thresh = roc(mat$bcr_status, mat[,gene])
  best_thresh = coords(best_thresh, x="best")
  mat$high_low = ifelse(mat[,gene] > best_thresh[1,1], "High", "Low")
  surv = survfit(Surv(days_to_follow_up, bcr_status) ~ high_low, data=mat)
  
  p <- ggsurvplot(surv, pval = TRUE, conf.int = F, risk.table = F, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "none", # Specify median survival
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("red1", "royalblue3"),
                  data=mat, xlab="Time (days)",ylab="Disease free survival probability")
  
  pdf(paste0("/data/github/pca_network/results/TCGA_DFS2/",gene,"_fpkm_bcr.pdf"), height=5, width=7)
  print(p)
  dev.off()
}

#########################################################################
# Kaplan MeierDFS scaled.. bcr 
#########################################################################
for(i in 1:nrow(signf_dfs)){
  
  row = as.data.frame(signf_dfs[i,])
  gene = row$gene
  
  mat = scaled[rownames(atlas_meta),]
  mat = merge(mat, atlas_meta, by=0)
  best_thresh = roc(mat$bcr_status, mat[,gene])
  best_thresh = coords(best_thresh, x="best")
  mat$high_low = ifelse(mat[,gene] > best_thresh[1,1], "High", "Low")
  surv = survfit(Surv(days_to_follow_up, bcr_status) ~ high_low, data=mat)
  
  p <- ggsurvplot(surv, pval = TRUE, conf.int = F, risk.table = F, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "none", # Specify median survival
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("red1", "royalblue3"),
                  data=mat, xlab="Time (days)",ylab="Disease free survival probability")
  
  pdf(paste0("/data/github/pca_network/results/TCGA_DFS2/",gene,"_scaled_bcr.pdf"), height=5, width=7)
  print(p)
  dev.off()
}
