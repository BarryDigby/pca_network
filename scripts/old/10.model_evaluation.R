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


# loading models from RData can cause problems? Run 9 and come here
#load("/data/github/pca_network/rev4/models.RData")
load("/data/github/pca_network/rev4/pre_filtering.RData")
load("/data/github/pca_network/rev4/processed_data.RDS")

## Strategy ---
## Extract (non-zero) coefficients from models and place the genes
## in a multivariate coxph model.
## We cannot use lasso/elastic coxboost models for validation, the shrinkage penalties applied
## produce linear predictors with slope != 1 in the derivation dataset.

# Model Evaluation

# Lasso/Enet extraction
extract_lasso <- function(model){
  genes <- coef(model$glmnet.fit, s=model$lambda.min)[,1]
  genes <- names(genes)[genes > 0]
  genes <- c(genes[!grepl("^study", genes)], "study")
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(genes, collapse = " + ")))
  fit <- coxph(form, derivation_mat, x=TRUE)
  return(fit)
}

multivariate_lasso_all <- extract_lasso(lasso_all)
multivariate_lasso_md <- extract_lasso(lasso_md)
multivariate_lasso_mrmr <- extract_lasso(lasso_mrmr)
multivariate_lasso_univ <- extract_lasso(lasso_univ)
multivariate_lasso_vh <- extract_lasso(lasso_vh)
multivariate_lasso_vimp <- extract_lasso(lasso_vimp)

# Elastic
multivariate_elastic_all <- extract_lasso(elastic_all)
multivariate_elastic_md <- extract_lasso(elastic_md)
multivariate_elastic_mrmr <- extract_lasso(elastic_mrmr)
multivariate_elastic_univ <- extract_lasso(elastic_univ)
multivariate_elastic_vh <- extract_lasso(elastic_vh)
multivariate_elastic_vimp <- extract_lasso(elastic_vimp)

# Coxboost
extract_coxboost <- function(model){
  genes <- coef(model)
  genes <- names(genes)[genes > 0]
  genes <- c(genes[!grepl("^study", genes)], "study")
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(genes, collapse = " + ")))
  fit <- coxph(form, derivation_mat, x=TRUE)
  return(fit)
}

multivariate_coxboost_all <- extract_coxboost(boost_all)
multivariate_coxboost_md <- extract_coxboost(boost_md)
multivariate_coxboost_mrmr <- extract_coxboost(boost_mrmr)
multivariate_coxboost_univ <- extract_coxboost(boost_univ)
multivariate_coxboost_vh <- extract_coxboost(boost_vh)
multivariate_coxboost_vimp <- extract_coxboost(boost_vimp)

## Multivariate Models

results <- data.frame()

models_used <- data.frame(ls(.GlobalEnv)[grep("multivariate", ls(.GlobalEnv))])
gene_filter <- rep(c("gene_mask", "rf_minimal_depth", "mRMR_filter", "univariate_filter", "rf_variable_hunter", "rf_importance"), 4)
models_used$set <- gene_filter
colnames(models_used)[1] <- "model"

## Loop over each model, generate stats for derivation and validation datasets

for(row in 1:nrow(models_used)){
  
  model_name <- models_used[row, 1]
  model <- get(models_used[row, 1])
  genes <- names(coef(model))
  genes <- c(genes[!grepl("^study", genes)], "study")
  
  df_dev <- derivation_mat[, c(genes, "days_to_follow_up", "bcr_status")]
  df_val <- validation_mat[, c(genes, "days_to_follow_up", "bcr_status")]
  
  dev_lp <- predict(model, df_dev, type="lp")
  val_lp <- predict(model, df_val, type="lp")
  
  df_dev$lp <- dev_lp
  df_val$lp <- val_lp
  
  # model slope
  fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, df_dev)
  fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, df_val)
  
  out_dev <- c(round(as.numeric(coef(fit_dev)), 3), paste(round(confint(fit_dev)[,1], 2), round(confint(fit_dev)[,2], 2), sep = "-"), summary(fit_dev)$coef[,5]) 
  out_val <- c(round(as.numeric(coef(fit_val)), 3), paste(round(confint(fit_val)[,1], 2), round(confint(fit_val)[,2], 2), sep = "-"), summary(fit_val)$coef[,5]) 
  
  # model fit/misspecification
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(genes, collapse = "+"), "+offset(lp)"))
  fit_dev <- coxph(formula = form, df_dev)
  fit_val <- coxph(formula = form, df_val)
  
  out_dev <- c(out_dev, c(round(2*(diff(fit_dev$loglik)), 2), round(1-pchisq(2*(diff(fit_dev$loglik)), 8), 3)))
  out_val <- c(out_val, c(round(2*(diff(fit_val$loglik)), 2), round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3)))
  
  # Harrells Cindex
  out_dev <- c(out_dev, c(round(as.numeric(rcorr.cens(-1*dev_lp, Surv(df_dev$days_to_follow_up, df_dev$bcr_status))[1]),2)))
  out_val <- c(out_val, c(round(as.numeric(rcorr.cens(-1*val_lp, Surv(df_val$days_to_follow_up, df_val$bcr_status))[1]),2)))
  
  # Kaplan Meier 1 indicates high risk group
  df_dev$group <- as.numeric(dev_lp > median(dev_lp))
  df_val$group <- as.numeric(val_lp > median(val_lp))
  
  fit_dev <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = df_dev)
  fit_val <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = df_val)
  
  # Hazard Ratio
  df_dev$group <- as.factor(df_dev$group)
  df_val$group <- as.factor(df_val$group)
  df_dev$group <- relevel(df_dev$group, ref = "0")
  df_val$group <- relevel(df_val$group, ref = "0")
  
  fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = df_dev)
  fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = df_val)
  
  out_dev <- c(out_dev, c(round(as.numeric(coef(fit_dev)), 3), paste(round(confint(fit_dev)[,1], 2), round(confint(fit_dev)[,2], 2), sep = "-"), summary(fit_dev)$coef[,5]))
  out_val <- c(out_val, c(round(as.numeric(coef(fit_val)), 3), paste(round(confint(fit_val)[,1], 2), round(confint(fit_val)[,2], 2), sep = "-"), summary(fit_val)$coef[,5]))
  
  # AUC DEV
  ROC.1 <- timeROC(T=df_dev$days_to_follow_up,
                   delta=df_dev$bcr_status,
                   marker=df_dev$lp,
                   cause=1,weighting="marginal",
                   times=c(365),
                   iid=TRUE)
  
  ROC.2  <- timeROC(T=df_dev$days_to_follow_up,
                    delta=df_dev$bcr_status,
                    marker=df_dev$lp,
                    cause=1,weighting="marginal",
                    times=c(floor(365*3)),
                    iid=TRUE)
  
  ROC.3 <- timeROC(T=df_dev$days_to_follow_up,
                   delta=df_dev$bcr_status,
                   marker=df_dev$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*5)),
                   iid=TRUE)
  
  ROC.4 <- timeROC(T=df_dev$days_to_follow_up,
                   delta=df_dev$bcr_status,
                   marker=df_dev$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*8)),
                   iid=TRUE)
  
  out_dev <- c(out_dev, 
               c(round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.4$AUC[2]), 2), paste(confint(ROC.4)$CI_AUC[,1], confint(ROC.4)$CI_AUC[,2], sep = "-")))
  
  ROC.1 <- timeROC(T=df_val$days_to_follow_up,
                   delta=df_val$bcr_status,
                   marker=df_val$lp,
                   cause=1,weighting="marginal",
                   times=c(365),
                   iid=TRUE)
  
  ROC.2  <- timeROC(T=df_val$days_to_follow_up,
                    delta=df_val$bcr_status,
                    marker=df_val$lp,
                    cause=1,weighting="marginal",
                    times=c(floor(365*3)),
                    iid=TRUE)
  
  ROC.3 <- timeROC(T=df_val$days_to_follow_up,
                   delta=df_val$bcr_status,
                   marker=df_val$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*5)),
                   iid=TRUE)
  
  ROC.4 <- timeROC(T=df_val$days_to_follow_up,
                   delta=df_val$bcr_status,
                   marker=df_val$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*8)),
                   iid=TRUE)
  
  out_val <- c(out_val, 
               c(round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.4$AUC[2]), 2), paste(confint(ROC.4)$CI_AUC[,1], confint(ROC.4)$CI_AUC[,2], sep = "-")))
  
  out_dev <- c(out_dev, paste(genes, collapse = ", "))
  out_val <- c(out_val, paste(genes, collapse = ", "))
  
  results <- rbind(results, t(as.data.frame(out_dev)))
  rownames(results)[nrow(results)] <- paste(model_name, "derivation", sep = "_")
  results <- rbind(results, t(as.data.frame(out_val)))
  rownames(results)[nrow(results)] <- paste(model_name, "validation", sep = "_")
  
}

colnames(results) <- c(
  "Slope Coef",
  "Slope CI",
  "Slope Pval",
  "ChiSq",
  "ChiSq Pval",
  "Harrells C",
  "Hazard Ratio",
  "HR CI",
  "HR Pval",
  "AUC 1yr",
  "CI 1yr",
  "AUC 3yr",
  "CI 3yr",
  "AUC 5yr",
  "CI 5yr",
  "AUC 8yr",
  "CI 8yr",
  "genes"
)

################################################################################
################################################################################
## External PCA signatures
################################################################################
################################################################################

sigs <- readxl::read_xlsx("/data/github/pca_network/data/pcadb_signatures/PCaDB_Prognostic_Signatures_Gene_List.xlsx")

# sigs %>% group_by(Signature) %>% count() %>% arrange(n) %>%print(n=30)

ext_mat <- master_mat[rownames(validation_mat), ]

ext_results <- data.frame()

for (author in unique(sigs$Signature)){
  model_name <- author
  genes <- unique(sigs$`HGNC Symbol`[which(sigs$Signature == author)])
  genes <- intersect(c(genes, "days_to_follow_up", "bcr_status"), colnames(ext_mat))
  mat = ext_mat[, genes]
  
  genes <- colnames(mat)
  if(length(grepl("-",genes)) > 0){
    genes <- gsub("-", ".", genes)
    colnames(mat) <- gsub("-", ".", colnames(mat))
  }
  
  genes <- setdiff(genes, c("days_to_follow_up", "bcr_status"))
  
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(genes, collapse = "+")))
  model <- coxph(form, mat, x=TRUE)
  
  df_val <- mat
  val_lp <- predict(model, df_val, type="lp")
  df_val$lp <- val_lp
  
  fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, df_val)
  out_val <- c(round(as.numeric(coef(fit_val)), 3), paste(round(confint(fit_val)[,1], 2), round(confint(fit_val)[,2], 2), sep = "-"), summary(fit_val)$coef[,5]) 
  
  # model fit/misspecification
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(genes, collapse = "+"), "+offset(lp)"))
  fit_val <- coxph(formula = form, df_val)
  out_val <- c(out_val, c(round(2*(diff(fit_val$loglik)), 2), round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3)))
  
  # Harrells Cindex
  out_val <- c(out_val, c(round(as.numeric(rcorr.cens(-1*val_lp, Surv(df_val$days_to_follow_up, df_val$bcr_status))[1]),2)))
  
  
  # Hazard Ratio
  df_val$group <- as.numeric(val_lp > median(val_lp))
  df_val$group <- as.factor(df_val$group)
  df_val$group <- relevel(df_val$group, ref = "0")
  
  fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = df_val)
  
  out_val <- c(out_val, c(round(as.numeric(coef(fit_val)), 3), paste(round(confint(fit_val)[,1], 2), round(confint(fit_val)[,2], 2), sep = "-"), summary(fit_val)$coef[,5]))
  
  
  ROC.1 <- timeROC(T=df_val$days_to_follow_up,
                   delta=df_val$bcr_status,
                   marker=df_val$lp,
                   cause=1,weighting="marginal",
                   times=c(365),
                   iid=TRUE)
  
  ROC.2  <- timeROC(T=df_val$days_to_follow_up,
                    delta=df_val$bcr_status,
                    marker=df_val$lp,
                    cause=1,weighting="marginal",
                    times=c(floor(365*3)),
                    iid=TRUE)
  
  ROC.3 <- timeROC(T=df_val$days_to_follow_up,
                   delta=df_val$bcr_status,
                   marker=df_val$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*5)),
                   iid=TRUE)
  
  ROC.4 <- timeROC(T=df_val$days_to_follow_up,
                   delta=df_val$bcr_status,
                   marker=df_val$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*8)),
                   iid=TRUE)
  
  out_val <- c(out_val, 
               c(round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.4$AUC[2]), 2), paste(confint(ROC.4)$CI_AUC[,1], confint(ROC.4)$CI_AUC[,2], sep = "-")))
  
  
  ext_results <- rbind(ext_results, t(as.data.frame(out_val)))
  rownames(ext_results)[nrow(ext_results)] <- model_name
  
}


colnames(ext_results) <- c(
  "Slope Coef",
  "Slope CI",
  "Slope Pval",
  "ChiSq",
  "ChiSq Pval",
  "Harrells C",
  "Hazard Ratio",
  "HR CI",
  "HR Pval",
  "AUC 1yr",
  "CI 1yr",
  "AUC 3yr",
  "CI 3yr",
  "AUC 5yr",
  "CI 5yr",
  "AUC 8yr",
  "CI 8yr"
)

write.table(ext_results, "/data/github/pca_network/rev4/external_signatures.txt", quote=F, sep="\t")
write.table(results, "/data/github/pca_network/rev4/deriv_valid.txt", quote=F, sep="\t")
