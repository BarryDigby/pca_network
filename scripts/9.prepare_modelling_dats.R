#!/usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(rms)

## Stage datasets
## duplicate gene IDs are filtered - highest variance selected. 
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
info <- biomaRt::getBM(attributes=c("hgnc_symbol",
                                    "ensembl_gene_id"),
                       mart = mart,
                       useCache=TRUE)


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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

tcga_mat <- as.data.frame(scale(mat, scale=T, center=T))
tcga_mat <- cbind(tcga_mat, meta)
tcga_mat$study = "TCGA"
rm(mat, meta, tcga_rds)

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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

cit_mat = as.data.frame(scale(mat, scale = T, center = T))
cit_mat <- cbind(cit_mat, meta)
cit_mat$study <- "CIT"
rm(mat, meta, cit_rds)

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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

cmap_mat = as.data.frame(scale(mat, scale = T, center = T))
cmap_mat <- cbind(cmap_mat, meta)
cmap_mat$study <- "CMap"
rm(mat, meta, cmap_rds)

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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

belfast_mat = as.data.frame(scale(mat, scale = T, center = T))
belfast_mat <- cbind(belfast_mat, meta)
belfast_mat$study <- "Belfast"
rm(mat, belfast_rds, meta)


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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

cpc_mat = as.data.frame(scale(mat, scale = T, center = T))
cpc_mat <- cbind(cpc_mat, meta)
cpc_mat$study <- "CPC"
rm(mat, meta, cpc)


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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

GSE54460_mat = as.data.frame(scale(mat, scale = T, center = T))
GSE54460_mat = cbind(GSE54460_mat, meta)
GSE54460_mat$study <- "GSE54460"
rm(mat, meta, GSE54460)


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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

stockholm_mat = as.data.frame(scale(mat, scale = T, center = T))
stockholm_mat = cbind(stockholm_mat, meta)
stockholm_mat$study <- "Stockholm"
rm(stockholm, mat, meta)

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

mat = merge(mat, info, by.x=0, by.y="ensembl_gene_id")
mat$variance <- apply(mat[ , 2:(ncol(mat) - 1)], 1, var)
mat <- mat %>% group_by(hgnc_symbol) %>% slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% ungroup %>% dplyr::select(-variance)
mat <- mat %>% filter(!(hgnc_symbol == ""))
mat <- tibble::column_to_rownames(mat, "hgnc_symbol")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

taylor_mat = as.data.frame(scale(mat, scale = T, center = T))
taylor_mat = cbind(taylor_mat, meta)
taylor_mat$study <- "Taylor"
rm(meta, mat, taylor)

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

## Finally, remove "-" characters from gene names, these will break model formulas
if(length(grepl("-",colnames(master_mat))) > 0){
  colnames(master_mat) <- gsub("-", ".", colnames(master_mat))
}

save(master_mat, file = "/data/github/pca_network/results/processed_datasets.RData")
