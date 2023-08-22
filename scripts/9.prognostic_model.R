#!/usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(survival)
library(survminer)
library(pROC)
library(glmnet)

##########################################################################
# 1. Overall survival 
# 2. Disease-free survival
##########################################################################

##########################################################################
# load human protein atlas FPKM data (only tumor samples)
##########################################################################
atlas_mat = read.csv("/data/github/pca_network/data/prad_rna_cancer_sample.tsv", header=T, sep="\t")
atlas_mat = atlas_mat[,c(1,2,4)]
atlas_mat <- atlas_mat %>%
  pivot_wider(names_from = Sample, values_from = FPKM, values_fill = 0)
atlas_mat = tibble::column_to_rownames(atlas_mat, "Gene")

#########################################################################
# load TCGA metadata from Protein Atlas.
# Why protein atlas? There is no missing information for both 
# time to death (days_to_follow_up) and vital status (outcome survival)
# This is crucial for downstream LASSO cox modelling
# Two vital status are 'not reported' 
# encode vital as numeric
##########################################################################
atlas_meta = read.csv("/data/github/pca_network/data/tcga_updated_meta.csv", header=T, sep=",")
rownames(atlas_meta) = atlas_meta$sample
atlas_meta = atlas_meta[which(atlas_meta$sample %in% colnames(atlas_mat)),]
rem = !(atlas_meta$vital=="Not Reported")
atlas_meta = atlas_meta[rem,]
atlas_mat = atlas_mat[,rem]
atlas_meta$status = ifelse(atlas_meta$vital=="Dead", 1, 0)


#########################################################################
# stage mrnas in ceRNA network
# Atlas FPKM uses Ensembl v109 ensembl_gene_ids.
# Load parsed GTF file to make sure HGNC maps to correct ensembl id
#########################################################################
network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
genes = unique(network$mrna)
ensv109 = read.csv("/data/github/pca_network/data/ensembl_v109_proteinatlas.csv", sep="\t", header = F)
colnames(ensv109) = c("ensembl_gene_id", "biotype", "hgnc_symbol")
ensv109 = ensv109[which(ensv109$biotype=="protein_coding"),]
ensv109 = ensv109[which(ensv109$hgnc_symbol %in% genes),]
gene_ids = ensv109[,c(3,1)]


#########################################################################
# Univariate Cox proportional hazards regression:
# Find ceRNA mRNAs ** with prognosis of PCa (Overall survival)
#########################################################################
options(scipen = 999)
result = data.frame(gene = character(),
                    ensembl = character(),
                    best_cutoff = numeric(),
                    pvalue = numeric())
for(i in 1:nrow(gene_ids)){
  name = gene_ids[i,1]
  ens  = gene_ids[i,2]
  # 4 genes are not in atlas DF as ENS. Info lost. 
  if(ens %in% rownames(atlas_mat)){
    mat  = data.frame(t(atlas_mat[which(rownames(atlas_mat)==ens),]))
    mat  = merge(mat, atlas_meta, by.x=0, by.y="sample")
    colnames(mat)[2] = name
    best_thresh = roc(mat$status, mat[,2])
    best_thresh = coords(best_thresh, x="best")
    mat$high_low = ifelse(mat[,2] > best_thresh[1,1], "High", "Low")
    res.cox = survfit(Surv(days_to_follow_up, status) ~ high_low, data=mat)
    res.cox.pval = surv_pvalue(res.cox)
    pval = res.cox.pval$pval
    if( pval < 0.05 ){
      row = data.frame(gene=name, ensembl = ens, best_cutoff=best_thresh[1,1], pvalue=pval)
      result = rbind(result,row)
    }
  }
}

#########################################################################
# Protein atlas defines prognostic genes as being < 0.01.
# This helps the LASSO cox model by reducing input features.
#########################################################################

signf = result[which(result$pvalue < 0.01),]

##########################################################################
# Sanity check:
# check how many signf genes are prognostic according to protein atlas
##########################################################################
prog_pca = read.csv("/data/github/pca_network/results/prognostic_prostate.tsv", header=T, sep="\t")
table(prog_pca$Gene %in% signf$gene)


##########################################################################
# Sanity check:
# check how many signf genes are prognostic according to protein atlas
##########################################################################
prog_pca = read.csv("/data/github/pca_network/results/prognostic_prostate.tsv", header=T, sep="\t")
table(prog_pca$Gene %in% signf$gene)


##########################################################################
# Stage data for LASSO Cox:
# Use normalised scaled and centered logcpm STAR counts data for modelling
# need to convert ENSG version IDs to HGNC symbols
# match samples in metadata
#
# ! The final line makes sure the expression data matches the metadata.
# ! The results are not reproducible if this is altered.
##########################################################################
scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
colnames(scaled) = gsub("\\.", "-", colnames(scaled))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
version_to_id = merge(signf, mrna_attributes, by.x="ensembl", by.y="ensembl_gene_id")
sub_atlas = scaled[which(rownames(scaled) %in% version_to_id$ensembl_gene_id_version),]
sub_atlas = merge(sub_atlas, version_to_id[,5:6], by.x=0, by.y="ensembl_gene_id_version")
sub_atlas = tibble::column_to_rownames(sub_atlas, "external_gene_name")
sub_atlas = sub_atlas[,c(2:ncol(sub_atlas))]
sub_atlas = sub_atlas[,rownames(atlas_meta)] 


##########################################################################
# Prepare X, Y matrices for LASSO COX
# X = signf genes
# Y = time to death, os status
##########################################################################
sub_meta = subset(atlas_meta, select=c(days_to_follow_up, status))
colnames(sub_meta) = c("time", "status")

x = as.matrix(t(sub_atlas))
y = sub_meta
y = as.matrix(y)


##########################################################################
# LASSO Cox 
# run and save as RData object
##########################################################################
set.seed(123)
# cv.fit <- cv.glmnet(x, y, family="cox", alpha=1, maxit = 1000, lambda = NULL, type.measure = "deviance")
# fit = glmnet(x, y, family = "cox", alpha=1, maxit = 1000, lambda=NULL)
# Coefficients <- coef(fit, s = cv.fit$lambda.min)
# Active.Index <- which(Coefficients != 0)
# Active.Coefficients <- Coefficients[Active.Index]
# Active.Genes <- Coefficients@Dimnames[[1]][Active.Index]
# save(fit, cv.fit, Active.Coefficients, Active.Genes, signf, file="/data/github/pca_network/results/prognostic_model_os.RData")

load("/data/github/pca_network/results/prognostic_model_os.RData")


####################################################################################################################################################
# Same analysis, --- Disease-free survival 
####################################################################################################################################################

# clear session prior to running 
rm(list=ls())

##########################################################################
# load the ceRNA network
##########################################################################
network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
genes = unique(network$mrna)

#########################################################################
# load TCGA metadata from Protein Atlas.
# Why protein atlas? There is no missing information for 
# time to death (days_to_follow_up)
# This is crucial for downstream LASSO cox modelling
##########################################################################
atlas_meta = read.csv("/data/github/pca_network/data/tcga_updated_meta.csv", header=T, sep=",")
rownames(atlas_meta) = atlas_meta$sample

# merge BCR status using PCa DB dataset
pca_db = readRDS("/data/github/pca_network/data/TCGA-PRAD_eSet.RDS")
pca_db = pca_db@phenoData@data
pca_db = pca_db[,c(1,27)]
pca_db$sample_id = paste(pca_db$sample_id, "A", sep="")
table(atlas_meta$Row.names %in% pca_db$sample_id)
# 478 matching PCa DB <-> atlas meta. this is better than 464 valid TCGA eset time to event.
atlas_meta = merge(atlas_meta, pca_db, by.x=0, by.y="sample_id")
keep = atlas_meta$sample_type=="Primary Tumor"
atlas_meta = atlas_meta[keep,]

#########################################################################
# Use scaled centered logcpm STAR counts for CoxPH/Kaplan Meier
# coefficients came from this data - must be used on this again
#########################################################################
scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
colnames(scaled) = gsub("\\.", "-", colnames(scaled))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
scaled = merge(scaled, mrna_attributes[,1:2], by.x=0, by.y="ensembl_gene_id_version")
scaled = tibble::column_to_rownames(scaled, "external_gene_name")
scaled = scaled[,c(2:ncol(scaled))]
scaled = scaled[,atlas_meta$Row.names]
#scaled = as.data.frame(t(scaled))

##########################################################################
# load human protein atlas FPKM data (only tumor samples)
##########################################################################
# atlas_mat = read.csv("/data/github/pca_network/data/prad_rna_cancer_sample.tsv", header=T, sep="\t")
# atlas_mat = atlas_mat[,c(1,2,4)]
# atlas_mat <- atlas_mat %>%
#   pivot_wider(names_from = Sample, values_from = FPKM, values_fill = 0)
# atlas_mat = tibble::column_to_rownames(atlas_mat, "Gene")

#########################################################################
# load TCGA metadata from Protein Atlas.
# Why protein atlas? There is no missing information for 
# time to death (days_to_follow_up)
# This is crucial for downstream LASSO cox modelling
##########################################################################
# atlas_meta = read.csv("/data/github/pca_network/data/tcga_updated_meta.csv", header=T, sep=",")
# rownames(atlas_meta) = atlas_meta$sample
# atlas_meta = atlas_meta[which(atlas_meta$sample %in% colnames(atlas_mat)),]
# 
# # merge BCR status using PCa DB dataset
# pca_db = readRDS("/data/github/pca_network/data/TCGA-PRAD_eSet.RDS")
# pca_db = pca_db@phenoData@data
# pca_db = pca_db[,c(1,27)]
# pca_db$sample_id = paste(pca_db$sample_id, "A", sep="")
# table(atlas_meta$Row.names %in% pca_db$sample_id)
# # 478 matching PCa DB <-> atlas meta. this is better than 464 valid TCGA eset time to event.
# atlas_meta = merge(atlas_meta, pca_db, by.x=0, by.y="sample_id")
# atlas_mat = atlas_mat[,atlas_meta$Row.names]

##########################################################################
# ENG ID -> HGNC symbol
##########################################################################

# mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
# mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
# atlas_mat = merge(atlas_mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
# atlas_mat = tibble::column_to_rownames(atlas_mat, "external_gene_name")
# atlas_mat = atlas_mat[,c(2:ncol(atlas_mat))]
# table(atlas_meta$Row.names %in% colnames(atlas_mat))

#########################################################################
# Univariate Cox proportional hazards regression:
# Find ceRNA mRNAs ** with prognosis of PCa (Disease-free survival)
#########################################################################

options(scipen = 999)
result = data.frame(gene = character(),
                    best_cutoff = numeric(),
                    pvalue = numeric())
for(i in genes){
    df = data.frame(t(scaled[which(rownames(scaled)==i),]))
    df = merge(df, atlas_meta, by.x=0, by.y="sample")
    best_thresh = roc(df$bcr_status, df[,2])
    best_thresh = coords(best_thresh, x="best")
    df$high_low = ifelse(df[,2] > best_thresh[1,1], "High", "Low")
    res.cox = survfit(Surv(days_to_follow_up, bcr_status) ~ high_low, data=df)
    res.cox.pval = surv_pvalue(res.cox)
    pval = res.cox.pval$pval
    if( pval < 0.05 ){
      row = data.frame(gene=i, best_cutoff=best_thresh[1,1], pvalue=pval)
      result = rbind(result,row)
    }
}

result$adj_p = p.adjust(result$pvalue, method = "BH")

#########################################################################
# Significant genes adj_p < 0.01
#########################################################################

signf = result[which(result$adj_p < 0.01),]

##########################################################################
# Stage data for LASSO Cox:
# Use normalised scaled and centered logcpm STAR counts data for modelling
# need to convert ENSG version IDs to HGNC symbols
# match samples in metadata
#
# ! The final line makes sure the expression data matches the metadata.
# ! The results are not reproducible if this is altered.
##########################################################################
# scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
# colnames(scaled) = gsub("\\.", "-", colnames(scaled))
# mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
# version_to_id = merge(signf, mrna_attributes, by.x="gene", by.y="external_gene_name")
# sub_atlas = scaled[which(rownames(scaled) %in% version_to_id$ensembl_gene_id_version),]
# sub_atlas = merge(sub_atlas, version_to_id[,c(1,5)], by.x=0, by.y="ensembl_gene_id_version")
# sub_atlas = tibble::column_to_rownames(sub_atlas, "gene")
# sub_atlas = sub_atlas[,c(2:ncol(sub_atlas))]
# sub_atlas = sub_atlas[,atlas_meta$Row.names] 
scaled = scaled[which(rownames(scaled) %in% signf$gene),]

##########################################################################
# Prepare X, Y matrices for LASSO COX
# X = signf genes
# Y = time to DFS 
##########################################################################
sub_meta = subset(atlas_meta, select=c(days_to_follow_up, bcr_status))
colnames(sub_meta) = c("time", "status")

x = as.matrix(t(scaled))
y = sub_meta
y = as.matrix(y)


##########################################################################
# LASSO Cox 
# run and save as RData object
##########################################################################
set.seed(123)
cv.fit <- cv.glmnet(x, y, family="cox", alpha=1, maxit = 1000, lambda = NULL, type.measure = "deviance")
fit = glmnet(x, y, family = "cox", alpha=1, maxit = 1000, lambda=NULL)
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Genes <- Coefficients@Dimnames[[1]][Active.Index]
save(fit, cv.fit, Active.Coefficients, Active.Genes, signf, file="/data/github/pca_network/results/prognostic_model_dfs.RData")
