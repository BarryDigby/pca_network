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
library(rms)
library(snowfall)
library(survAUC)
library(meta)
library(metafor)
library(survcomp)
library(randomForestSRC)
library(mRMRe)

## set the outdir for results - files + plots + script used to generate.
# useful if people want to tweak params 
outdir <- "/data/github/pca_network/results/"
outdir <- paste(outdir, Sys.Date(), sep="")
if(!dir.exists(outdir)){
  dir.create(outdir)
}

## General Rule of thumb - 15 events per predictor.
## 398/15 = 26.5.
## Aim for models with fewer predictors than the above.
set.seed(123789)
network_genes <- read.table("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
network_genes <- unique(network_genes$mrna)
load("/data/github/pca_network/results/processed_datasets.RData")
network_genes <- network_genes[which(network_genes %in% colnames(master_mat))]
mat <- master_mat[, c(network_genes, "preop_psa", "gleason_score", "days_to_follow_up", "bcr_status", "study")]
rm(master_mat)



# STEP 1



# LASSO model

## construct x and y matrix
keep_patient <- rownames(mat)[which(mat$days_to_follow_up != 0)]
y <- as.matrix(mat[keep_patient, c("days_to_follow_up", "bcr_status")])
colnames(y) <- c("time", "status") 
x <- model.matrix( ~ 0 + ., mat[keep_patient, c(network_genes, "gleason_score", "preop_psa")])

## use "study" as cross validation folds
group_identifier <- as.numeric(as.factor(mat$study[which(rownames(mat) %in% keep_patient)]))

## LASSO implementation
lambda_max <- max(glmnet(x, y, family = "cox", alpha = 1, relax = F)$lambda)
lambda_seq <- exp(seq(log(lambda_max), log(lambda_max/200), length = 10)) # smaller lambdas = sparser model
lasso_fit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_seq, foldid = group_identifier, type.measure = "C", relax = F)

## extract named vector of coefs
lasso_preds <- coef(lasso_fit$glmnet.fit, s=lasso_fit$lambda.min)[,1]
lasso_preds <- lasso_preds[lasso_preds != 0]

## function to remove offending predictors (schoenfeld residuals.)
remove_offending_predictors <- function(mat, lasso_preds) {
  repeat {
    # Create the formula
    fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(names(lasso_preds), collapse = " + ")))
    
    # Fit the Cox proportional hazards model with lasso coefficients
    fit_zph <- cox.zph(coxph(fx, data = mat, init = as.numeric(lasso_preds), iter.max = 0))
    
    # Identify predictors to remove based on Schoenfeld residuals
    remove <- rownames(fit_zph$table)[which(fit_zph$table[, 3] < 0.05)]
    
    # Remove "GLOBAL" if present
    if ("GLOBAL" %in% remove) {
      remove <- remove[remove != "GLOBAL"]
    }
    
    # Remove the offending predictors from lasso_preds
    lasso_preds <- lasso_preds[!(names(lasso_preds) %in% remove)]
    
    # Break the loop if no more offending predictors are found
    if (length(remove) == 0) {
      break
    }
  }
  
  # Final Cox proportional hazards model using LASSO coefficients
  fit <- coxph(fx, data = mat, init = as.numeric(lasso_preds), iter.max = 0)
  
  return(list(fit = fit, lasso_preds = lasso_preds))
}

remove_offending_predictors2 <- function(mat, lasso_preds, verbose = TRUE) {
  max_iter <- 10
  iter <- 1
  
  repeat {
    if (verbose) message("Schoenfeld check iteration: ", iter)
    
    fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(names(lasso_preds), collapse = " + ")))
    
    # Match order for init values
    X <- model.matrix(fx, data = mat)[, -1, drop = FALSE]
    init_vals <- lasso_preds[colnames(X)]
    
    # Catch failure in fit
    fit_try <- try(coxph(fx, data = mat, init = as.numeric(init_vals), iter.max = 0), silent = TRUE)
    if (inherits(fit_try, "try-error")) {
      warning("CoxPH model failed during Schoenfeld check. Returning current predictors.")
      break
    }
    
    # Schoenfeld residuals
    fit_zph <- try(cox.zph(fit_try), silent = TRUE)
    if (inherits(fit_zph, "try-error")) {
      warning("cox.zph failed — possibly singular or perfect separation. Returning current predictors.")
      break
    }
    
    remove <- rownames(fit_zph$table)[which(fit_zph$table[, 3] < 0.05)]
    remove <- setdiff(remove, "GLOBAL")
    
    if (verbose && length(remove) > 0) {
      message("Removing violating predictors: ", paste(remove, collapse = ", "))
    }
    
    lasso_preds <- lasso_preds[!(names(lasso_preds) %in% remove)]
    
    if (length(remove) == 0 || length(lasso_preds) <= 2 || iter >= max_iter) break
    
    iter <- iter + 1
  }
  
  # Final model
  final_fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(names(lasso_preds), collapse = " + ")))
  final_X <- model.matrix(final_fx, data = mat)[, -1, drop = FALSE]
  final_init <- lasso_preds[colnames(final_X)]
  
  final_fit <- try(coxph(final_fx, data = mat, init = as.numeric(final_init), iter.max = 0), silent = TRUE)
  
  if (inherits(final_fit, "try-error")) {
    warning("Final CoxPH model failed.")
    final_fit <- NULL
  }
  
  return(list(fit = final_fit, lasso_preds = lasso_preds))
}

schoenfeld_result <- remove_offending_predictors(mat, lasso_preds)

lasso_preds <- schoenfeld_result$lasso_preds
fit <- schoenfeld_result$fit

# STEP 2

lasso_master_result <- data.frame()
dev_mat <- mat

# Assess the models perfomance in the master dataset

# Linear predictors and probabilities
Xb_dev <- model.matrix(fit) %*% coef(fit)
Xbavg_dev <- sum(coef(fit)*fit$means) 
dev_lp <- Xb_dev - Xbavg_dev
dev_mat$lp <- dev_lp
dev_mat$prob <- predict(fit, dev_mat, "risk")

# model calibration slope
fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, dev_mat)
out_dev <- c(round(as.numeric(coef(fit_dev)), 4), round(confint(fit_dev)[,1], 4), round(confint(fit_dev)[,2], 4), summary(fit_dev)$coef[,5])

# model calibration slope in the large - want to be as close to zero. 
# code borrowed from predRupdate
z_val <- stats::qnorm(1 - (1-0.95)/2)
dev_expected_events <- mean(dev_mat$prob) * length(dev_mat$prob)

KM_dev <- summary(survival::survfit(Surv(dev_mat$days_to_follow_up, dev_mat$bcr_status) ~ 1), times = max(dev_mat$days_to_follow_up))
dev_log_OE_ratio <- log(KM_dev$n.event) - log(dev_expected_events)
dev_log_OE_ratio_se <- sqrt(KM_dev$surv / KM_dev$n.event)
dev_OE_ratio <- exp(dev_log_OE_ratio)
dev_OE_ratio_se <- sqrt(KM_dev$n.event * (1 - KM_dev$n.event / dev_expected_events) / dev_expected_events)
dev_OE_ratio_lower <- exp(dev_log_OE_ratio - (z_val * dev_log_OE_ratio_se))
dev_OE_ratio_upper <- exp(dev_log_OE_ratio + (z_val * dev_log_OE_ratio_se))
out_dev <- c(out_dev, c(round(dev_OE_ratio, 4), round(dev_OE_ratio_lower, 4), round(dev_OE_ratio_upper, 4)))


# model fit/misspecification - offset is same as coef = 1. Append to coef init vector 
form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(lasso_preds), collapse = "+"), " + offset(lp)"))
fit_dev <- coxph(formula = form, data = dev_mat, init = c(as.numeric(lasso_preds, 1)), iter.max=0)
out_dev <- c(out_dev, c(round(2*(diff(fit_dev$loglik)), 2), round(1-pchisq(2*(diff(fit_dev$loglik)), 8), 3)))


# Harrells Cindex
dev_c <- survcomp::concordance.index(x = dev_mat$lp, surv.time = dev_mat$days_to_follow_up, surv.event = dev_mat$bcr_status)
out_dev <- c(out_dev, round(dev_c$c.index, 4), round(dev_c$lower, 4), round(dev_c$upper, 4))

# Kaplan Meier 1 indicates high risk group
dev_mat$group <- as.numeric(dev_mat$lp > median(dev_mat$lp))
fit_dev <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = dev_mat)

# Hazard Ratio
dev_mat$group <- as.factor(dev_mat$group)
dev_mat$group <- relevel(dev_mat$group, ref = "0")
fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = dev_mat)
out_dev <- c(out_dev, c(round(as.numeric(coef(fit_dev)), 4), round(as.numeric(summary(fit_dev)$coef[,3]), 4),  paste(round(confint(fit_dev)[,1], 2), round(confint(fit_dev)[,2], 2), sep = "-"), summary(fit_dev)$coef[,5]))

# AUC
ROC.1 <- timeROC(T=dev_mat$days_to_follow_up,
                 delta=dev_mat$bcr_status,
                 marker=dev_mat$lp,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)

ROC.2  <- timeROC(T=dev_mat$days_to_follow_up,
                  delta=dev_mat$bcr_status,
                  marker=dev_mat$lp,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=dev_mat$days_to_follow_up,
                 delta=dev_mat$bcr_status,
                 marker=dev_mat$lp,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

out_dev <- c(out_dev, 
             c(round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
               round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
               round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-")))


out_dev <- c(out_dev, paste(names(coef(fit)), collapse = ", "))


lasso_master_result <- rbind(lasso_master_result, t(as.data.frame(out_dev)))

colnames(lasso_master_result) <- c(
  "Slope Coef", # possibly not relevant when carrying lasso coefs across. 
  "Slope LCI",
  "Slope UCI",
  "Slope Pval",
  "Slope ITL Coef",
  "Slope ITL LCI",
  "Slope ITL UCI",
  "ChiSq",
  "ChiSq Pval",
  "Harrells C",
  "Harrells C LCI",
  "Harrells C UCI",
  "Hazard Ratio",
  "HR SE",
  "HR CI",
  "HR Pval",
  "AUC 1yr",
  "CI 1yr",
  "AUC 3yr",
  "CI 3yr",
  "AUC 5yr",
  "CI 5yr",
  "genes"
)

lasso_master_result <- lasso_master_result%>%
  mutate(Dataset = "master_dataset")%>%
  dplyr::select(Dataset, everything())

write.table(lasso_master_result, file = paste(outdir, "/lasso_master.txt", sep=""), row.names = F, quote=F, sep="\t")

## STEP 3

# Internal-external cross validation of each study fold.

gen_lp <- function(dev_model, val_model, output){
  
  Xb_dev <- model.matrix(dev_model) %*% coef(dev_model)
  Xbavg_dev <- sum(coef(dev_model)*dev_model$means) 
  dev_lp <- Xb_dev - Xbavg_dev
  
  Xb_val <- model.matrix(val_model) %*% coef(dev_model) # determine Xb in validation data using B dev
  val_lp <- Xb_val - Xbavg_dev # center PI by using mean of PI of development data
  
  if(output == "dev"){
    return(dev_lp)
  }else if(output == "val"){
    return(val_lp)
  }
}


lasso_iecv_result <- data.frame()

for(i in unique(mat$study)){
  print(i)
  dev_mat <- mat %>% dplyr::filter(!(study %in% i))
  val_mat <- mat %>% dplyr::filter(study %in% i)
  dev_mat <- dev_mat[, c(names(lasso_preds), "days_to_follow_up", "bcr_status")]
  val_mat <- val_mat[, c(names(lasso_preds), "days_to_follow_up", "bcr_status")]
  
  # construct coxph model - use lasso coeficients
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(lasso_preds), collapse = " + ")))
  dev_fit <- coxph(fx, dev_mat, init = as.numeric(lasso_preds), iter.max=0)
  val_fit <- coxph(fx, val_mat, init = as.numeric(lasso_preds), iter.max=0)
  
  # Linear predictors and probabilities  
  dev_mat$lp <- gen_lp(dev_fit, val_fit, "dev")
  val_mat$lp <- gen_lp(dev_fit, val_fit, "val")
  dev_mat$prob <- predict(dev_fit, dev_mat, "risk")
  val_mat$prob <- predict(dev_fit, val_mat, "risk")
  
  # model calibration slope
  fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, dev_mat)
  fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, val_mat)
  
  out_dev <- c(round(as.numeric(coef(fit_dev)), 4), round(confint(fit_dev)[,1], 4), round(confint(fit_dev)[,2], 4), summary(fit_dev)$coef[,5])
  out_val <- c(round(as.numeric(coef(fit_val)), 4), round(confint(fit_val)[,1], 4), round(confint(fit_val)[,2], 4), summary(fit_val)$coef[,5])
  
  # model calibration slope in the large - want to be as close to zero. 
  # code borrowed from predRupdate
  z_val <- stats::qnorm(1 - (1-0.95)/2)
  dev_expected_events <- mean(dev_mat$prob) * length(dev_mat$prob)
  val_expected_events <- mean(val_mat$prob) * length(val_mat$prob)
  
  KM_dev <- summary(survival::survfit(Surv(dev_mat$days_to_follow_up, dev_mat$bcr_status) ~ 1), times = max(dev_mat$days_to_follow_up))
  dev_log_OE_ratio <- log(KM_dev$n.event) - log(dev_expected_events)
  dev_log_OE_ratio_se <- sqrt(KM_dev$surv / KM_dev$n.event)
  dev_OE_ratio <- exp(dev_log_OE_ratio)
  dev_OE_ratio_se <- sqrt(KM_dev$n.event * (1 - KM_dev$n.event / dev_expected_events) / dev_expected_events)
  dev_OE_ratio_lower <- exp(dev_log_OE_ratio - (z_val * dev_log_OE_ratio_se))
  dev_OE_ratio_upper <- exp(dev_log_OE_ratio + (z_val * dev_log_OE_ratio_se))
  
  KM_val <- summary(survival::survfit(Surv(val_mat$days_to_follow_up, val_mat$bcr_status) ~ 1), times = max(val_mat$days_to_follow_up))
  val_log_OE_ratio <- log(KM_val$n.event) - log(mean(val_mat$prob) * length(val_mat$prob))
  val_log_OE_ratio_se <- sqrt(KM_val$surv / KM_val$n.event)
  val_OE_ratio <- exp(val_log_OE_ratio)
  val_OE_ratio_se <- sqrt(KM_val$n.event * (1 - KM_val$n.event / val_expected_events) / val_expected_events)
  val_OE_ratio_lower <- exp(val_log_OE_ratio - (z_val * val_log_OE_ratio_se))
  val_OE_ratio_upper <- exp(val_log_OE_ratio + (z_val * val_log_OE_ratio_se))
  
  out_dev <- c(out_dev, c(round(dev_OE_ratio, 4), round(dev_OE_ratio_lower, 4), round(dev_OE_ratio_upper, 4)))
  out_val <- c(out_val, c(round(val_OE_ratio, 4), round(val_OE_ratio_lower, 4), round(val_OE_ratio_upper, 4)))
  
  # model fit/misspecification - offset is same as coef = 1. Append to coef init vector 
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(lasso_preds), collapse = "+"), " + offset(lp)"))
  fit_dev <- coxph(formula = form, data = dev_mat, init = c(as.numeric(lasso_preds, 1)), iter.max=0)
  fit_val <- coxph(formula = form, data = val_mat, init = c(as.numeric(lasso_preds, 1)), iter.max=0)
  
  out_dev <- c(out_dev, c(round(2*(diff(fit_dev$loglik)), 2), round(1-pchisq(2*(diff(fit_dev$loglik)), 8), 3)))
  out_val <- c(out_val, c(round(2*(diff(fit_val$loglik)), 2), round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3)))
  
  # Harrells Cindex
  dev_c <- survcomp::concordance.index(x = dev_mat$lp, surv.time = dev_mat$days_to_follow_up, surv.event = dev_mat$bcr_status)
  val_c <- survcomp::concordance.index(x = val_mat$lp, surv.time = val_mat$days_to_follow_up, surv.event = val_mat$bcr_status)
  out_dev <- c(out_dev, round(dev_c$c.index, 4), round(dev_c$lower, 4), round(dev_c$upper, 4))
  out_val <- c(out_val, round(val_c$c.index, 4), round(val_c$lower, 4), round(val_c$upper, 4))
  
  # Kaplan Meier 1 indicates high risk group
  dev_mat$group <- as.numeric(dev_mat$lp > median(dev_mat$lp))
  val_mat$group <- as.numeric(val_mat$lp > median(val_mat$lp))
  
  fit_dev <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = dev_mat)
  fit_val <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = val_mat)
  
  # Hazard Ratio
  dev_mat$group <- as.factor(dev_mat$group)
  val_mat$group <- as.factor(val_mat$group)
  dev_mat$group <- relevel(dev_mat$group, ref = "0")
  val_mat$group <- relevel(val_mat$group, ref = "0")
  
  fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = dev_mat)
  fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = val_mat)
  
  out_dev <- c(out_dev, c(round(as.numeric(coef(fit_dev)), 4), round(as.numeric(summary(fit_dev)$coef[,3]), 4),  paste(round(confint(fit_dev)[,1], 2), round(confint(fit_dev)[,2], 2), sep = "-"), summary(fit_dev)$coef[,5]))
  out_val <- c(out_val, c(round(as.numeric(coef(fit_val)), 4), round(as.numeric(summary(fit_val)$coef[,3]), 4),  paste(round(confint(fit_val)[,1], 2), round(confint(fit_val)[,2], 2), sep = "-"), summary(fit_val)$coef[,5]))
  
  # AUC
  ROC.1 <- timeROC(T=dev_mat$days_to_follow_up,
                   delta=dev_mat$bcr_status,
                   marker=dev_mat$lp,
                   cause=1,weighting="marginal",
                   times=c(365),
                   iid=TRUE)
  
  ROC.2  <- timeROC(T=dev_mat$days_to_follow_up,
                    delta=dev_mat$bcr_status,
                    marker=dev_mat$lp,
                    cause=1,weighting="marginal",
                    times=c(floor(365*3)),
                    iid=TRUE)
  
  ROC.3 <- timeROC(T=dev_mat$days_to_follow_up,
                   delta=dev_mat$bcr_status,
                   marker=dev_mat$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*5)),
                   iid=TRUE)
  
  out_dev <- c(out_dev, 
               c(round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-")))
  
  # validation 
  
  ROC.1 <- timeROC(T=val_mat$days_to_follow_up,
                   delta=val_mat$bcr_status,
                   marker=val_mat$lp,
                   cause=1,weighting="marginal",
                   times=c(365),
                   iid=TRUE)
  
  ROC.2  <- timeROC(T=val_mat$days_to_follow_up,
                    delta=val_mat$bcr_status,
                    marker=val_mat$lp,
                    cause=1,weighting="marginal",
                    times=c(floor(365*3)),
                    iid=TRUE)
  
  ROC.3 <- timeROC(T=val_mat$days_to_follow_up,
                   delta=val_mat$bcr_status,
                   marker=val_mat$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*5)),
                   iid=TRUE)
  
  out_val <- c(out_val, 
               c(round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
                 round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-")))
  
  out_dev <- c(out_dev, paste(names(coef(dev_fit)), collapse = ", "))
  out_val <- c(out_val, paste(names(coef(val_fit)), collapse = ", "))
  
  lasso_iecv_result <- rbind(lasso_iecv_result, t(as.data.frame(out_dev)))
  rownames(lasso_iecv_result)[nrow(lasso_iecv_result)] <- paste0(i, "_dev")
  lasso_iecv_result <- rbind(lasso_iecv_result, t(as.data.frame(out_val)))
  rownames(lasso_iecv_result)[nrow(lasso_iecv_result)] <- paste0(i, "_val")
}

colnames(lasso_iecv_result) <- c(
  "Slope Coef",
  "Slope LCI",
  "Slope UCI",
  "Slope Pval",
  "Slope ITL Coef",
  "Slope ITL LCI",
  "Slope ITL UCI",
  "ChiSq",
  "ChiSq Pval",
  "Harrells C",
  "Harrells C LCI",
  "Harrells C UCI",
  "Hazard Ratio",
  "HR SE",
  "HR CI",
  "HR Pval",
  "AUC 1yr",
  "CI 1yr",
  "AUC 3yr",
  "CI 3yr",
  "AUC 5yr",
  "CI 5yr",
  "genes"
)

lasso_iecv_result <- lasso_iecv_result%>%
  tibble::rownames_to_column(var = "Dataset")%>%
  dplyr::select(Dataset, everything())

write.table(lasso_iecv_result, file = paste(outdir, "/lasso_iecv.txt", sep=""), row.names = F, quote=F, sep="\t")

## Add information for plot
studies <- unique(mat$study)
res <- sapply(studies, function(study) {
  subset_mat <- mat[mat$study == study, ]
  n_events <- sum(subset_mat$bcr_status == 1, na.rm = TRUE)
  n_total <- nrow(subset_mat)
  sprintf("%s (%d/%d)", study, n_events, n_total)
}) 
lasso_iecv_result$info <- rep(res, each=2)
## random effects meta analysis
lasso_iecv_result <- lasso_iecv_result%>%
  arrange(Dataset)%>%
  tibble::column_to_rownames(var = "Dataset")

forest_plot <- function(model, effect, LCI, UCI, rlab){
  print(paste0("Filtering for ", model))
  df <- lasso_iecv_result %>% tibble::rownames_to_column(., "rownames") %>% dplyr::filter(grepl(paste0(model), rownames))
  df <- df %>% tibble::column_to_rownames(., "info")
  df <- df %>% dplyr::select(one_of(effect, LCI, UCI)) %>% mutate_all(as.numeric)
  met <- metagen(TE = df[, effect], lower = df[, LCI], upper = df[, UCI], studlab = rownames(df))
  ref_line <- ifelse(rlab == "Calibration Slope", 1, ifelse(rlab == "Harrell's C", NA, 0))
  pad = 0.1
  if (rlab == "Calibration Slope") {
    xlim_vals <- c(floor(min(met$lower, na.rm = TRUE)),
                   ceiling(max(met$upper, na.rm = TRUE)))
  } else if (rlab == "Calibration ITL") {
    xlim_vals <- c(-0.2,
                   max(met$upper, na.rm = TRUE) + pad)
  } else if (rlab == "Harrell's C") {
    xlim_vals <- c(max(0, min(met$lower, na.rm = TRUE) - pad), 1)
  }
  return(meta::forest(met, common = FALSE, rightlabs=c(rlab, "95% CI", "Weight"), leftcols=c("studlab"), leftlabs=c("Dataset (events/total)"),addrows.below.overall=3, ref=ref_line, xlim=xlim_vals))
}

pdf(file = paste(outdir, "/lasso_iecv_slope.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Slope Coef", "Slope LCI", "Slope UCI", "Calibration Slope")
dev.off()

pdf(file = paste(outdir, "/lasso_iecv_itl.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Slope ITL Coef", "Slope ITL LCI", "Slope ITL UCI", "Calibration ITL")
dev.off()

pdf(file = paste(outdir, "/lasso_iecv_harrells.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Harrells C", "Harrells C LCI", "Harrells C UCI", "Harrell's C")
dev.off()


## Step 4
rm(out_val)
lasso_ext_val <- data.frame()

## unseen, external dataset validation
## external val
library(GEOquery)
library(hgu133plus2.db)

GSE46602 <- getGEO("GSE46602")
meta <- GSE46602$GSE46602_series_matrix.txt.gz@phenoData@data
meta <- meta %>%
  dplyr::rename(
    `days_to_follow_up` = `bcr_free_time:ch1`, 
    `bcr_status` = `bcr:ch1`, 
    `gleason_score` = `gleason_grade:ch1`, 
    `preop_psa` = `preop_psa:ch1`) %>% 
  dplyr::select(days_to_follow_up, bcr_status, gleason_score, preop_psa)%>%
  mutate(across(everything(), ~ if_else(. == "NA", NA_character_, .)))%>%
  mutate(bcr_status = ifelse(bcr_status == "YES", 1, 0))%>%
  mutate(days_to_follow_up = 30.44 * as.numeric(days_to_follow_up))%>%
  mutate(across(everything(), as.numeric))

meta <- na.omit(meta)

ext_mat <- GSE46602$GSE46602_series_matrix.txt.gz@assayData$exprs
ext_mat <- as.data.frame(ext_mat[, rownames(meta)])

x <- mapIds(hgu133plus2.db, keys = rownames(ext_mat), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
x <- data.frame(probe = names(x), symbol = x, stringsAsFactors = FALSE)
#x <- x %>% filter(symbol %in% genes)

ext_mat <- ext_mat[x$probe, ]
ext_mat$Amean <- rowMeans(ext_mat)
ext_mat <- ext_mat %>%
  tibble::rownames_to_column("probe") %>%
  left_join(x, by = "probe")

ext_mat <- ext_mat %>%
  arrange(desc(Amean)) %>%
  distinct(symbol, .keep_all = TRUE)

ext_mat <- na.omit(ext_mat)
rownames(ext_mat) <- NULL

ext_mat <- ext_mat %>%
  dplyr::select(-Amean, -probe) %>%
  tibble::column_to_rownames("symbol")


ext_mat <- data.frame(t(ext_mat))
ext_mat <- data.frame(scale(ext_mat, scale = T, center = T))

ext_mat <- cbind(ext_mat, meta)

ext_mat <- ext_mat[, c(names(lasso_preds), "days_to_follow_up", "bcr_status")]

# External, unseen dataset results
fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(lasso_preds), collapse = " + ")))
ext_fit <- coxph(fx, data = ext_mat, init = as.numeric(lasso_preds), iter.max=0)

# fit is coxph on all 9 datasets above
ext_mat$lp <- gen_lp(fit, ext_fit, "val")
ext_mat$prob <- predict(ext_fit, ext_mat, "risk")


fit_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, ext_mat)
print("Calibration slope")
out_val <- c(round(as.numeric(coef(fit_ext)), 4), round(confint(fit_ext)[,1], 4), round(confint(fit_ext)[,2], 4), summary(fit_ext)$coef[,5])

z_val <- stats::qnorm(1 - (1-0.95)/2)
ext_expected_events <- mean(ext_mat$prob) * length(ext_mat$prob)

KM_ext <- summary(survival::survfit(Surv(ext_mat$days_to_follow_up, ext_mat$bcr_status) ~ 1), times = max(ext_mat$days_to_follow_up))
ext_log_OE_ratio <- log(KM_ext$n.event) - log(ext_expected_events)
ext_log_OE_ratio_se <- sqrt(KM_ext$surv / KM_ext$n.event)
ext_OE_ratio <- exp(ext_log_OE_ratio)
ext_OE_ratio_se <- sqrt(KM_ext$n.event * (1 - KM_ext$n.event / ext_expected_events) / ext_expected_events)
ext_OE_ratio_lower <- exp(ext_log_OE_ratio - (z_val * ext_log_OE_ratio_se))
ext_OE_ratio_upper <- exp(ext_log_OE_ratio + (z_val * ext_log_OE_ratio_se))

print("Calibration slope in the large")
out_val <- c(out_val, round(ext_OE_ratio, 4), round(ext_OE_ratio_lower, 4), round(ext_OE_ratio_upper, 4))

form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(lasso_preds), collapse = "+"), " + offset(lp)"))
fit_ext <- coxph(formula = form, data = ext_mat, init = c(as.numeric(lasso_preds, 1)), iter.max=0)

print("ChiSq model fit")
out_val <- c(out_val, round(2*(diff(fit_ext$loglik)), 2), round(1-pchisq(2*(diff(fit_ext$loglik)), 8), 3))

ext_c <- survcomp::concordance.index(x = ext_mat$lp, surv.time = ext_mat$days_to_follow_up, surv.event = ext_mat$bcr_status)
print("Harrells C")
out_val <- c(out_val, round(ext_c$c.index, 4), round(ext_c$lower, 4), round(ext_c$upper, 4))

ext_mat$group <- as.numeric(ext_mat$lp > median(ext_mat$lp))
fit_ext <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = ext_mat)
ext_mat$group <- as.factor(ext_mat$group)
ext_mat$group <- relevel(ext_mat$group, ref = "0")
fit_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = ext_mat)

print("Hazard Ratio")
out_val <- c(out_val, round(as.numeric(coef(fit_ext)), 4), round(as.numeric(summary(fit_ext)$coef[,3]), 4),  paste(round(confint(fit_ext)[,1], 2), round(confint(fit_ext)[,2], 2), sep = "-"), summary(fit_ext)$coef[,5])

ROC.1 <- timeROC(T=ext_mat$days_to_follow_up,
                 delta=ext_mat$bcr_status,
                 marker=ext_mat$lp,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)

ROC.2  <- timeROC(T=ext_mat$days_to_follow_up,
                  delta=ext_mat$bcr_status,
                  marker=ext_mat$lp,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=ext_mat$days_to_follow_up,
                 delta=ext_mat$bcr_status,
                 marker=ext_mat$lp,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

print("AUC values")
out_val <- c(out_val,
      round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
      round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
      round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-"))

out_val <- c(out_val, paste(names(coef(ext_fit)), collapse = ", "))

lasso_ext_val <- rbind(lasso_ext_val, t(as.data.frame(out_val)))

colnames(lasso_ext_val) <- c(
  "Slope Coef",
  "Slope LCI",
  "Slope UCI",
  "Slope Pval",
  "Slope ITL Coef",
  "Slope ITL LCI",
  "Slope ITL UCI",
  "ChiSq",
  "ChiSq Pval",
  "Harrells C",
  "Harrells C LCI",
  "Harrells C UCI",
  "Hazard Ratio",
  "HR SE",
  "HR CI",
  "HR Pval",
  "AUC 1yr",
  "CI 1yr",
  "AUC 3yr",
  "CI 3yr",
  "AUC 5yr",
  "CI 5yr",
  "genes"
)

lasso_ext_val <- lasso_ext_val%>%
  mutate(Dataset = "GSE46602_ext")%>%
  dplyr::select(Dataset, everything())

write.table(lasso_ext_val, file = paste(outdir, "/lasso_external.txt", sep=""), row.names = F, quote=F, sep="\t")

source_file <- "/data/github/pca_network/scripts/10-LASSO.R"
file.copy(source_file, outdir, overwrite = T)

