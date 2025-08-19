#!/usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(survival)
library(RegParallel)
library(survminer)
library(timeROC)
library(ggplot2)
library(glmnet)
library(mfp)
library(CoxBoost)
library(boot)
library(parallel)
library(randomForestSRC)
library(rms)
library(mRMRe)
library(survival)
library(rsample)

## Load ceRNA network genes
network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
genes = unique(network$mrna)
rm(network)

## Stage datasets
## do not filter by network genes yet - external signatures use different genes.. 
## TCGA PRAD
tcga_rds <- readRDS("/data/github/pca_network/data/TCGA-PRAD_eSet.RDS")
mat <- tcga_rds@assayData$exprs
meta <- tcga_rds@phenoData@data
meta <- meta %>%
  filter(sample_type == "Primary") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

tcga_mat <- as.data.frame(scale(mat, scale=T, center=T))
tcga_mat <- cbind(tcga_mat, meta)
tcga_mat$study = "TCGA"
rm(mat, meta, tcga_rds, mrna_attributes)

## CIT
cit_rds <- readRDS("/data/github/pca_network/data/CIT_eSet.RDS")
mat = cit_rds@assayData$exprs
meta = cit_rds@phenoData@data
meta <- meta %>%
  filter(sample_type == "Primary") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

cit_mat = as.data.frame(scale(mat, scale = T, center = T))
cit_mat <- cbind(cit_mat, meta)
cit_mat$study <- "CIT"
rm(mat, meta, mrna_attributes, cit_rds)

## Cancer Map
cmap_rds <- readRDS("/data/github/pca_network/data/CancerMap_eSet.RDS")
mat = cmap_rds@assayData$exprs
meta = cmap_rds@phenoData@data
meta <- meta %>%
  filter(sample_type != "Normal") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

cmap_mat = as.data.frame(scale(mat, scale = T, center = T))
cmap_mat <- cbind(cmap_mat, meta)
cmap_mat$study <- "CMap"
rm(mat, meta, mrna_attributes, cmap_rds)

## Belfast
belfast_rds = readRDS("/data/github/pca_network/data/Belfast_eSet.RDS")
mat = belfast_rds@assayData$exprs
meta = belfast_rds@phenoData@data
meta <- meta %>%
  filter(sample_type == "Primary") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

belfast_mat = as.data.frame(scale(mat, scale = T, center = T))
belfast_mat <- cbind(belfast_mat, meta)
belfast_mat$study <- "Belfast"
rm(mat, meta, mrna_attributes)


## CPC
cpc = readRDS("/data/github/pca_network/data/CPC-Gene_eSet.RDS")
mat = cpc@assayData$exprs
meta = cpc@phenoData@data
meta <- meta %>%
  filter(sample_type == "Primary") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  mutate(gleason_score = NA) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

cpc_mat = as.data.frame(scale(mat, scale = T, center = T))
cpc_mat <- cbind(cpc_mat, meta)
cpc_mat$study <- "CPC"
rm(mat, meta, cpc, mrna_attributes)


## DKFZ
dkfz = readRDS("/data/github/pca_network/data/DKFZ_eSet.RDS")
mat = dkfz@assayData$exprs
meta = dkfz@phenoData@data
meta <- meta %>%
  filter(sample_type == "Primary") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

dkfz_mat = as.data.frame(scale(mat, scale = T, center = T))
dkfz_mat = cbind(dkfz_mat, meta)
dkfz_mat$study <- "DKFZ"
rm(mat, meta, mrna_attributes, dkfz)

## GSE54460
GSE54460 = readRDS("/data/github/pca_network/data/GSE54460_eSet.RDS")
mat = GSE54460@assayData$exprs
meta = GSE54460@phenoData@data
meta <- meta %>%
  filter(sample_type == "Primary") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

GSE54460_mat = as.data.frame(scale(mat, scale = T, center = T))
GSE54460_mat = cbind(GSE54460_mat, meta)
GSE54460_mat$study <- "GSE54460"
rm(mat, meta, GSE54460, mrna_attributes)


## Stockholm
stockholm = readRDS("/data/github/pca_network/data/Stockholm_eSet.RDS")
mat = stockholm@assayData$exprs
meta = stockholm@phenoData@data
meta <- meta %>%
  filter(sample_type == "Primary") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]


mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

stockholm_mat = as.data.frame(scale(mat, scale = T, center = T))
stockholm_mat = cbind(stockholm_mat, meta)
stockholm_mat$study <- "Stockholm"
rm(stockholm, mat, meta, mrna_attributes)

## Taylor
taylor = readRDS("/data/github/pca_network/data/Taylor_eSet.RDS")
mat = taylor@assayData$exprs
meta = taylor@phenoData@data
meta <- meta %>%
  filter(sample_type == "Primary") %>%
  filter(!(time_to_bcr == 0 | is.na(bcr_status) | is.na(time_to_bcr))) %>%
  dplyr::rename(days_to_follow_up = time_to_bcr) %>%
  mutate(days_to_follow_up = 30.44 * days_to_follow_up) %>%
  mutate(preop_psa = NA) %>%
  dplyr::select(preop_psa, gleason_score, days_to_follow_up, bcr_status)
mat <- mat[, which(colnames(mat) %in% rownames(meta))]

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var) 
mat <- mat %>% group_by(external_gene_name) %>% slice_max(order_by = variance, n = 1) %>% ungroup %>% dplyr::select(-variance)
mat <- tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

taylor_mat = as.data.frame(scale(mat, scale = T, center = T))
taylor_mat = cbind(taylor_mat, meta)
taylor_mat$study <- "Taylor"
rm(meta, mat, taylor, mrna_attributes)

## Ensure common genes across datasets
common_dataset_genes <- table(c(colnames(belfast_mat), colnames(cpc_mat), colnames(dkfz_mat), colnames(GSE54460_mat), colnames(stockholm_mat), colnames(taylor_mat), colnames(tcga_mat), colnames(cit_mat), colnames(cmap_mat)))
common_dataset_genes <- names(common_dataset_genes)[which(common_dataset_genes == 9)]

## Impuation of missing PSA and Gleason
vec <- common_dataset_genes
master_mat <- rbind(cit_mat[, vec], tcga_mat[, vec], cmap_mat[, vec], cpc_mat[, vec], GSE54460_mat[, vec], stockholm_mat[, vec], belfast_mat[, vec], taylor_mat[, vec], dkfz_mat[, vec])
set.seed(11384191)
master_mat$gleason_score <- factor(master_mat$gleason_score)
w <- transcan( ~ preop_psa + gleason_score, imputed=TRUE, data=master_mat, pl=FALSE, pr=FALSE)
attach(master_mat)
preop_psa_impute = impute(w, preop_psa, data=master_mat)
gleason_score_impute = impute(w, gleason_score, data=master_mat)
master_mat$preop_psa <- as.numeric(preop_psa_impute)
master_mat$gleason_score <- as.numeric(as.character(gleason_score_impute))

## find common network genes
network_genes <- intersect(colnames(master_mat), genes)
network_genes <- c(network_genes, c("gleason_score", "preop_psa", "bcr_status", "days_to_follow_up", "study"))

################################################################################
################################################################################
## Internal-external validation
################################################################################
################################################################################
##
## Strategy here is to run LOO for each study. Results in a lot of models, 
## which we can use to justify genes and methods selected for final model
## on all dataset. Internal-external validation: doi:10.1016/j.jclinepi.2015.04.005

mat <- master_mat[, network_genes]
save(mat, network_genes, file = "/data/github/pca_network/reviewer4/internal_external_start_point.RData")

# system call, run in parallel.. 


################################################################################
################################################################################
## Pre-filtering
################################################################################
################################################################################

# The transcriptome has been filtered using a circRNA-miRNA-mRNA knowledge-based
# filter. We apply additional filtering steps to the 226 genes commmon to all 
# datasets used in the derivation and validation dataset. 
# 1. No filtering 
# 2. Univariate filtering
# 3. Random Forest Variable Importance
# 4. Random Forest Minimum Depth
# 5. Random Forest Varible Hunting
# 6. mRMR 

# ! These steps take a long time to run, 

# 1. No filtering
gene_mask <- colnames(derivation_mat)[which(!(colnames(derivation_mat) %in% c("days_to_follow_up", "bcr_status")))]

# 2. Univariate Filtering
res <- c()
for(gene in gene_mask){
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", gene, " + preop_psa + gleason_score + study"))
  fit <- coxph(fx, derivation_mat)
  pv <- summary(fit)$coef[1,5]
  names(pv) <- gene
  res <- c(res, pv)
}
univariate_filter <- res[res < 0.01]
univariate_filter <- names(univariate_filter)
univariate_filter <- unique(c(univariate_filter, "gleason_score", "preop_psa", "study"))

# 3. Random Forest Variable Importance
# inclue study, important covariate here 
rf_mat <- derivation_mat #%>% dplyr::select(-"study")
rf_mat$study <- factor(rf_mat$study)
rf_imp <- rfsrc(Surv(days_to_follow_up, bcr_status) ~., rf_mat, ntree = 500, nsplit = 10, mtry = sqrt(223), nodesize = 3, importance = TRUE)
rf_imp_sample <- subsample(rf_imp, B = 100)
rf_imp_stats <- plot(rf_imp_sample)$stats
rf_importance <- colnames(rf_imp_stats)
rf_importance <- unique(c(rf_importance,"gleason_score", "preop_psa", "study"))

# 4. Random Forest Minimum Depth
rf_md <- var.select(Surv(days_to_follow_up, bcr_status) ~ ., rf_mat, method = "md", ntree = 500, nsplit = 10, nodesize = 3, splitrule = "logrank")
rf_minimal_depth <- rf_md$topvars
rf_minimal_depth <- unique(c(rf_minimal_depth, "gleason_score", "preop_psa", "study"))

# 5. Random Forest Variable Hunting
rf_vh <- var.select(Surv(days_to_follow_up, bcr_status) ~ ., rf_mat, method = "vh", ntree = 500, nrep = 25, nsplit = 10, K = 5, nodesize = 3, nstep =1, splitrule = "logrank")
rf_variable_hunter <- rf_vh$topvars
rf_variable_hunter <- unique(c(rf_variable_hunter,"gleason_score", "preop_psa", "study"))

# 6. mRMR
common <- c()

for( i in unique(derivation_mat$study)){
  
  mat <- derivation_mat %>% dplyr::filter(!(study %in% i))
  mr_mat <- mat[, gene_mask]
  mr_mat$study <- factor(mr_mat$study, ordered = TRUE)
  mr_mat$target <- Surv(mat$days_to_follow_up, mat$bcr_status)
  mr_data <- mRMR.data(data=mr_mat)
  mr_classic <- mRMR.classic(data=mr_data, target_indices=as.integer(ncol(mr_mat)), feature_count = 50)
  mr_genes <- mr_classic@feature_names[unlist(mRMRe::solutions(mr_classic))]
  common <- c(common, mr_genes)
  
}

mRMR_filter <- names(table(common)[table(common) >= 9])
mRMR_filter <- unique(c(mRMR_filter, c("gleason_score", "preop_psa", "study")))

# save filtering results
save(list = c("univariate_filter", "mRMR_filter", "rf_variable_hunter", "rf_importance", "rf_minimal_depth"), file = "/data/github/pca_network/rev4/cit/pre_filtering.RData")
#load("/data/github/pca_network/rev4/pre_filtering.RData")


################################################################################
################################################################################
## ML Methods for feature selection
################################################################################
################################################################################

# 1. Multivariate Coxph - 'baseline'
# 2. GLMnet Lasso
# 3. GLMnet Elastic Net
# 4. CoxBoost

# The inputs to each algorithm are as follows:
  
# 1. All (no filter)
# 2. Univariate Coxph
# 3. RF Variable Importance
# 4. RF Minimum Depth
# 5. RF Variable Hunter
# 6. mRMR Classic

# Resulting in 24 models for evaluation

# Create derviation subsets based on pre-filterin results
all_mat <- derivation_mat[, c(gene_mask, "days_to_follow_up", "bcr_status")]
univ_mat <- derivation_mat[, c(univariate_filter, "days_to_follow_up", "bcr_status")]
vimp_mat <- derivation_mat[, c(rf_importance, "days_to_follow_up", "bcr_status")]
vh_mat <- derivation_mat[, c(rf_variable_hunter, "days_to_follow_up", "bcr_status")]
md_mat <- derivation_mat[, c(rf_minimal_depth, "days_to_follow_up", "bcr_status")]
mrmr_mat <- derivation_mat[, c(mRMR_filter, "days_to_follow_up", "bcr_status")]

# GLMNet model matrices
# identify samples with time to follow up == 0
keep_patient <- rownames(derivation_mat)[which(derivation_mat$days_to_follow_up != 0)]
y <- as.matrix(derivation_mat[keep_patient, c("days_to_follow_up", "bcr_status")])
colnames(y) <- c("time", "status") 
all_x <- model.matrix( ~ 0 + ., derivation_mat[keep_patient, c(gene_mask)])
univ_x <- model.matrix( ~ 0 + ., derivation_mat[keep_patient, c(univariate_filter)])
vimp_x <- model.matrix( ~ 0 + ., derivation_mat[keep_patient, c(rf_importance)])
vh_x <- model.matrix(~ 0 + ., derivation_mat[keep_patient, c(rf_variable_hunter)])
md_x <- model.matrix(~0 + ., derivation_mat[keep_patient, c(rf_minimal_depth)])
mrmr_x <- model.matrix( ~ 0 + ., derivation_mat[keep_patient, c(mRMR_filter)])

## Multivariate Coxph

multivariate_cox_all <- coxph(Surv(days_to_follow_up, bcr_status) ~ ., all_mat, x = TRUE)
multivariate_cox_univ <- coxph(Surv(days_to_follow_up, bcr_status) ~ ., univ_mat, x = TRUE)
multivariate_cox_vimp <- coxph(Surv(days_to_follow_up, bcr_status) ~ ., vimp_mat, x = TRUE)
multivariate_cox_vh <- coxph(Surv(days_to_follow_up, bcr_status) ~ ., vh_mat, x = TRUE)
multivariate_cox_md <- coxph(Surv(days_to_follow_up, bcr_status) ~ ., md_mat, x = TRUE)
multivariate_cox_mrmr <- coxph(Surv(days_to_follow_up, bcr_status) ~ ., mrmr_mat, x = TRUE)
print("multivariate models done")

## GLMnet Lasso
## Folds are split according to each study in the derivation dataset
group_identifier <- as.numeric(as.factor(derivation_mat$study[which(rownames(derivation_mat) %in% keep_patient)]))
pf <- ifelse(grepl("gleason|preop", colnames(all_x)), 0, 1)
lasso_all <- glmnet::cv.glmnet(all_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop", colnames(univ_x)), 0, 1)
lasso_univ <- glmnet::cv.glmnet(univ_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop", colnames(vimp_x)), 0, 1)
lasso_vimp <- glmnet::cv.glmnet(vimp_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop", colnames(vh_x)), 0, 1)
lasso_vh <- glmnet::cv.glmnet(vh_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop", colnames(md_x)), 0, 1)
lasso_md <- glmnet::cv.glmnet(md_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop", colnames(mrmr_x)), 0, 1)
lasso_mrmr <- glmnet::cv.glmnet(mrmr_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
print("Lasso models done")

## GLMnet Elastic Net
# Folds split similarly to LASSO
elastic_net_cv <- function(x) {
  
  set.seed(123)
  alpha_seq = seq(0.1, 1, length = 10)
  fitEN <- list()
  pf <- ifelse(grepl("gleason|preop", colnames(x)), 0, 1)
  
  for (i in 1:length(alpha_seq)) {
    fitEN[[i]] <- cv.glmnet(x, y, family = "cox", alpha = alpha_seq[i], foldid = group_identifier, penalty.factor = pf)
  }
  
  idx <- which.min(sapply(fitEN, function(xx) {
    xx$cvm[xx$lambda == xx$lambda.min]
  }))
  
  return(fitEN[[idx]])
}

elastic_all <- elastic_net_cv(all_x)
elastic_univ <- elastic_net_cv(univ_x)
elastic_vimp <- elastic_net_cv(vimp_x)
elastic_vh <- elastic_net_cv(vh_x)
elastic_md <- elastic_net_cv(md_x)
elastic_mrmr <- elastic_net_cv(mrmr_x)
print("Elastic Net models done")

## Cox Boost 

cox_boost <- function(x, y){
  cv.res <- cv.CoxBoost(time=y[,1], status=y[,2], x=x, maxstepno=500,K=5,type="verweij",penalty=100)
  cbfit <- CoxBoost(time=y[,1], status=y[,2], x=x, stepno=cv.res$optimal.step,penalty=100)
  return(cbfit)
}

boost_all <- cox_boost(all_x, y)
boost_univ <- cox_boost(univ_x, y)
boost_vimp <- cox_boost(vimp_x, y)
boost_vh <- cox_boost(vh_x, y)
boost_md <- cox_boost(md_x, y)
boost_mrmr <- cox_boost(mrmr_x, y)
print("Coxboost models done")

# save models
save_models <- ls(.GlobalEnv)[grep("^(boost_|elastic_|lasso_|multivariate_)", ls(.GlobalEnv))]
save(list = save_models, file="/data/github/pca_network/rev4/models.RData")

# save processed datasets
save_datasets <- ls(.GlobalEnv)[grep("*_mat", ls(.GlobalEnv))]
save(list = save_datasets, file = "/data/github/pca_network/rev4/processed_data.RData")

