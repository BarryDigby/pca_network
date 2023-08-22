#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
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
# Use FPKM for plotting expression data (scatterhist)
#########################################################################
# atlas_mat = read.csv("/data/github/pca_network/data/prad_rna_cancer_sample.tsv", header=T, sep="\t")
# atlas_mat = atlas_mat[,c(1,2,4)]
# atlas_mat <- atlas_mat %>%
#   pivot_wider(names_from = Sample, values_from = FPKM, values_fill = 0)
# atlas_mat = tibble::column_to_rownames(atlas_mat, "Gene")

#########################################################################
# Load atlas metadata
# Remove missing OS status
#########################################################################
# atlas_meta = read.csv("/data/github/pca_network/data/tcga_updated_meta.csv", header=T, sep=",")
# rownames(atlas_meta) = atlas_meta$sample
# atlas_meta = atlas_meta[which(atlas_meta$sample %in% colnames(atlas_mat)),]
# rem = !(atlas_meta$vital=="Not Reported")
# atlas_meta = atlas_meta[rem,]
# atlas_mat = atlas_mat[,rem]
# atlas_meta$status = ifelse(atlas_meta$vital=="Dead", 1, 0)

#########################################################################
# Load cleaned tcga meta
#########################################################################
load("/data/github/pca_network/data/TCGA_meta_cleaned.RData")

#########################################################################
# Convert ENS to symbols and subset matrix by Active.Genes,
# order correctly (FPKM)
# FPKM matrix now ready for expression plots
#########################################################################
load("/data/github/pca_network/results/prognostic_model_os.RData")
# ensv109 = read.csv("/data/github/pca_network/data/ensembl_v109_proteinatlas.csv", sep="\t", header = F)
# colnames(ensv109) = c("ensembl_gene_id", "biotype", "hgnc_symbol")
# ensv109 = ensv109[which(ensv109$hgnc_symbol %in% Active.Genes),]
# atlas_mat = atlas_mat[which(rownames(atlas_mat) %in% ensv109$ensembl_gene_id),]
# atlas_mat = merge(atlas_mat, ensv109[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
# atlas_mat = tibble::column_to_rownames(atlas_mat, "hgnc_symbol")
# atlas_mat = atlas_mat[,c(2:ncol(atlas_mat))]
# atlas_mat = as.data.frame(t(atlas_mat))
# atlas_mat = atlas_mat[rownames(atlas_meta),]

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
scaled = scaled[,which(colnames(scaled) %in% master$Row.names)]
scaled = as.data.frame(t(scaled))

#########################################################################
# Load TCGA RDS file and assess clinical meta vs atlas meta
#########################################################################
# tcga_meta = readRDS("/data/github/pca_network/data/TCGA-PRAD_eSet.RDS")
# tcga_meta = tcga_meta@phenoData@data

#########################################################################
# Construct dataframe for analysis
#########################################################################

mat = scaled
risk = rowSums(Active.Coefficients*mat)
mat$risk_score = risk
mat$risk_category = ifelse(mat$risk_score > median(mat$risk_score), "high", "low")
mat = merge(mat, master, by.x=0, by.y="Row.names")

#########################################################################
# Perform univariate coxph, add signf to multivariate cox after. 
# Only 10 events, so need to bin patients e.g T2 vs T3 for path T
#########################################################################

mat$path_t2 = ifelse(grepl("T2", mat$ajcc_pathologic_t) == TRUE, 1,0)
mat$path_n = ifelse(grepl("N1", mat$ajcc_pathologic_n) == TRUE, 1,0)
mat$clin_m = ifelse(grepl("M0", mat$ajcc_clinical_m) == TRUE, 1,0)
mat$risk_cat_num = ifelse(mat$risk_category=="high", 1,0)

mat$risk_category = relevel(factor(mat$risk_category), ref="low")
vars_for_table = c("age","preop_psa", "gleason_score", "risk_category", "path_t2")

univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(days_to_follow_up, status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mat)})

forestmodel::forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)



forestmodel::forest_model(coxph(Surv(days_to_follow_up, status) ~ age + preop_psa + gleason_score + risk_score + path_t2, data=mat))
