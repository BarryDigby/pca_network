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
# Load and prepare data
set.seed(123789)

# Read network genes and filter to those present in data
network_genes <- read.table("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header = TRUE, sep = "\t")
network_genes <- unique(network_genes$mrna)
load("/data/github/pca_network/results/processed_datasets.RData")
network_genes <- network_genes[which(network_genes %in% colnames(master_mat))]

mat <- master_mat[, c(network_genes, "preop_psa", "gleason_score", "days_to_follow_up", "bcr_status", "study")]
rm(master_mat)

# Keep only patients with follow-up time
keep_patient <- rownames(mat)[mat$days_to_follow_up != 0]
mat <- mat[keep_patient, ]
y <- as.matrix(mat[, c("days_to_follow_up", "bcr_status")])
colnames(y) <- c("time", "status")
x <- model.matrix(~ 0 + ., mat[, c(network_genes, "gleason_score", "preop_psa")])
n <- nrow(x)

snowfall::sfInit(parallel=TRUE, cpus=6)

# Apparent model (Step 1)
cat("Fitting apparent model\n")
cb.pen.apparent <- optimCoxBoostPenalty(time = y[, 1], status = y[, 2], x = x,
                                        minstepno = 5, maxstepno = 35,
                                        start.penalty = 9 * sum(y[, 2] == 1),
                                        iter.max = 10, upper.margin = 0.1,
                                        parallel=TRUE)

cb_fit_apparent <- CoxBoost(time = y[, 1], status = y[, 2], x = x,
                            stepno = cb.pen.apparent$cv.res$optimal.step,
                            penalty = cb.pen.apparent$penalty * 1.5)

coxboost_preds <- coef(cb_fit_apparent)
coxboost_preds <- coxboost_preds[coxboost_preds != 0]

# Schoenfeld filtering
remove_offending_predictors <- function(mat, lasso_preds) {
  repeat {
    fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(names(lasso_preds), collapse = " + ")))
    fit_attempt <- try({
      fit_zph <- cox.zph(coxph(fx, data = mat, init = as.numeric(lasso_preds), iter.max = 0))
    }, silent = TRUE)
    
    if (inherits(fit_attempt, "try-error")) {
      warning("cox.zph failed — possibly singular or perfect separation. Returning current predictors.")
      break
    }
    
    remove <- rownames(fit_zph$table)[which(fit_zph$table[, 3] < 0.05)]
    remove <- setdiff(remove, "GLOBAL")
    if (length(remove) == 0) break
    lasso_preds <- lasso_preds[!(names(lasso_preds) %in% remove)]
  }
  
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(names(lasso_preds), collapse = " + ")))
  fit <- coxph(fx, data = mat, init = as.numeric(lasso_preds), iter.max = 0)
  return(list(fit = fit, preds = lasso_preds))
}

# Apply Schoenfeld filtering to apparent model
sch_mat_app <- as.data.frame(cbind(x, y))
colnames(sch_mat_app)[(ncol(sch_mat_app) - 1):ncol(sch_mat_app)] <- c("days_to_follow_up", "bcr_status")
schoenfeld_result_app <- remove_offending_predictors(sch_mat_app, coxboost_preds)
fit <- schoenfeld_result_app$fit

# Apparent C-index
lp_app <- predict(fit, newdata = mat[, names(coef(fit)), drop = FALSE], type = "lp")
c_index_apparent <- survcomp::concordance.index(x = lp_app, surv.time = y[, "time"], surv.event = y[, "status"])$c.index



# STEP 2

coxboost_master_result <- data.frame()
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
form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coxboost_preds), collapse = "+"), " + offset(lp)"))
fit_dev <- coxph(formula = form, data = dev_mat, init = c(as.numeric(coxboost_preds, 1)), iter.max=0)
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


coxboost_master_result <- rbind(coxboost_master_result, t(as.data.frame(out_dev)))

colnames(coxboost_master_result) <- c(
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


coxboost_master_result <- coxboost_master_result%>%
  mutate(Dataset = "master_dataset")%>%
  dplyr::select(Dataset, everything())

write.table(coxboost_master_result, file = paste(outdir, "/coxboost_master.txt", sep=""), row.names = F, quote=F, sep="\t")

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


coxboost_iecv_result <- data.frame()

for(i in unique(mat$study)){
  print(i)
  dev_mat <- mat %>% dplyr::filter(!(study %in% i))
  val_mat <- mat %>% dplyr::filter(study %in% i)
  dev_mat <- dev_mat[, c(names(coxboost_preds), "days_to_follow_up", "bcr_status")]
  val_mat <- val_mat[, c(names(coxboost_preds), "days_to_follow_up", "bcr_status")]
  
  # construct coxph model - use lasso coeficients
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coxboost_preds), collapse = " + ")))
  dev_fit <- coxph(fx, dev_mat, init = as.numeric(coxboost_preds), iter.max=0)
  val_fit <- coxph(fx, val_mat, init = as.numeric(coxboost_preds), iter.max=0)
  
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
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coxboost_preds), collapse = "+"), " + offset(lp)"))
  fit_dev <- coxph(formula = form, data = dev_mat, init = c(as.numeric(coxboost_preds, 1)), iter.max=0)
  fit_val <- coxph(formula = form, data = val_mat, init = c(as.numeric(coxboost_preds, 1)), iter.max=0)
  
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
  
  coxboost_iecv_result <- rbind(coxboost_iecv_result, t(as.data.frame(out_dev)))
  rownames(coxboost_iecv_result)[nrow(coxboost_iecv_result)] <- paste0(i, "_dev")
  coxboost_iecv_result <- rbind(coxboost_iecv_result, t(as.data.frame(out_val)))
  rownames(coxboost_iecv_result)[nrow(coxboost_iecv_result)] <- paste0(i, "_val")
}

colnames(coxboost_iecv_result) <- c(
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

coxboost_iecv_result <- coxboost_iecv_result%>%
  tibble::rownames_to_column(var = "Dataset")%>%
  dplyr::select(Dataset, everything())

write.table(coxboost_iecv_result, file = paste(outdir, "/coxboost_iecv.txt", sep=""), row.names = F, quote=F, sep="\t")

## Add information for plot
studies <- unique(mat$study)
res <- sapply(studies, function(study) {
  subset_mat <- mat[mat$study == study, ]
  n_events <- sum(subset_mat$bcr_status == 1, na.rm = TRUE)
  n_total <- nrow(subset_mat)
  sprintf("%s (%d/%d)", study, n_events, n_total)
}) 
coxboost_iecv_result$info <- rep(res, each=2)
## random effects meta analysis
coxboost_iecv_result <- coxboost_iecv_result%>%
  arrange(Dataset) %>%
  tibble::column_to_rownames(var = "Dataset")

forest_plot <- function(model, effect, LCI, UCI, rlab){
  print(paste0("Filtering for ", model))
  df <- coxboost_iecv_result %>% tibble::rownames_to_column(., "rownames") %>% dplyr::filter(grepl(paste0(model), rownames))
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

pdf(file = paste(outdir, "/coxboost_iecv_slope.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Slope Coef", "Slope LCI", "Slope UCI", "Calibration Slope")
dev.off()

pdf(file = paste(outdir, "/coxboost_iecv_itl.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Slope ITL Coef", "Slope ITL LCI", "Slope ITL UCI", "Calibration ITL")
dev.off()

pdf(file = paste(outdir, "/coxboost_iecv_harrells.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Harrells C", "Harrells C LCI", "Harrells C UCI", "Harrell's C")
dev.off()



## Step 4
rm(out_val)
coxboost_ext_val <- data.frame()

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

ext_mat <- ext_mat[, c(names(coxboost_preds), "days_to_follow_up", "bcr_status")]

# External, unseen dataset results
fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coxboost_preds), collapse = " + ")))
ext_fit <- coxph(fx, data = ext_mat, init = as.numeric(coxboost_preds), iter.max=0)

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

form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coxboost_preds), collapse = "+"), " + offset(lp)"))
fit_ext <- coxph(formula = form, data = ext_mat, init = c(as.numeric(coxboost_preds, 1)), iter.max=0)

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

coxboost_ext_val <- rbind(coxboost_ext_val, t(as.data.frame(out_val)))

colnames(coxboost_ext_val) <- c(
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

coxboost_ext_val <- coxboost_ext_val%>%
  mutate(Dataset = "GSE46602_ext")%>%
  dplyr::select(Dataset, everything())

write.table(coxboost_ext_val, file = paste(outdir, "/coxboost_external.txt", sep=""), row.names = F, quote=F, sep="\t")

source_file <- "/data/github/pca_network/scripts/10-Coxboost.R"
file.copy(source_file, outdir, overwrite = T)

library(rms)

# datadist requirement
dd <- datadist(mat)
options(datadist = "dd")

# Refit model using rms
form <- formula(paste0("Surv(days_to_follow_up, bcr_status) ~ ", paste0(names(coef(fit)), collapse="+")))
fit_rms <- cph(form, init = platt_model$linear.predictors,
               iter.max = 0,
               data = mat,
               x = TRUE, y = TRUE, surv = TRUE)

shrink_fit <- shrink(fit_rms, type = "linear")

# Time points (in days, assuming your time variable is in days)
time_points <- c(1, 3, 5) * 365

# Loop over time points and store calibration objects
i = 0
for (time in time_points){
  i = i + 1
  cal_list[[i]] = rms::calibrate(fit_rms, u = time, cmethod = "KM", method = "boot", B = 500)
}


# Plot all on one page
par(mfrow = c(1, 3))
for (i in seq_along(time_points)) {
  plot(cal_list[[i]],
       xlab = "Predicted probability",
       ylab = "Observed probability",
       main = paste(time_points[i] / 365, "Years"))
}

