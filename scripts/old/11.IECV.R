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


load("/data/github/pca_network/results/processed_datasets.RData")
mat <- master_mat
rm(master_mat)
## Perform internal-external cross validation of model.
predictors <- readRDS("/data/github/pca_network/results/predictors.RDS")
predictors <- predictors$predictors

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

set.seed(123)

results <- data.frame()

for(i in unique(mat$study)){
  print(i)
  dev_mat <- mat %>% dplyr::filter(!(study %in% i))
  val_mat <- mat %>% dplyr::filter(study %in% i)
  dev_mat <- dev_mat[, c(predictors, "days_to_follow_up", "bcr_status")]
  val_mat <- val_mat[, c(predictors, "days_to_follow_up", "bcr_status")]
  
  # construct coxph model - use lasso coeficients
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(predictors, collapse = " + ")))
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
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(predictors, collapse = "+"), " + offset(lp)"))
  fit_dev <- coxph(formula = form, dev_mat)
  fit_val <- coxph(formula = form, val_mat)
  
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
  
  results <- rbind(results, t(as.data.frame(out_dev)))
  rownames(results)[nrow(results)] <- paste0(i, "_dev")
  results <- rbind(results, t(as.data.frame(out_val)))
  rownames(results)[nrow(results)] <- paste0(i, "_val")
}

colnames(results) <- c(
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



forest_plot <- function(model, effect, LCI, UCI){
  print(paste0("Filtering for ", model))
  df <- results %>% tibble::rownames_to_column(., "rownames") %>% dplyr::filter(grepl(paste0(model), rownames)) %>% tibble::column_to_rownames(., "rownames")
  df <- df %>% dplyr::select(one_of(effect, LCI, UCI)) %>% mutate_all(as.numeric)
  met <- metagen(TE = df[, effect], lower = df[, LCI], upper = df[, UCI], studlab = rownames(df))
  return(meta::forest(met, common = FALSE))
}

forest_plot("_val", "Slope Coef", "Slope LCI", "Slope UCI")
forest_plot("_val", "Slope ITL Coef", "Slope ITL LCI", "Slope ITL UCI")
forest_plot("_val", "Harrells C", "Harrells C LCI", "Harrells C UCI")

## performance on full dataset
dev_mat <- mat[, c(predictors, "days_to_follow_up", "bcr_status")]
dev_fit <- coxph(Surv(days_to_follow_up, bcr_status) ~ ., dev_mat)

dev_mat$lp <- gen_lp(dev_fit, val_fit, "dev")
dev_mat$prob <- predict(dev_fit, dev_mat, "risk")

fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, dev_mat)
print("Calibration slope")
cat(c(round(as.numeric(coef(fit_dev)), 4), round(confint(fit_dev)[,1], 4), round(confint(fit_dev)[,2], 4), summary(fit_dev)$coef[,5]))

z_val <- stats::qnorm(1 - (1-0.95)/2)
dev_expected_events <- mean(dev_mat$prob) * length(dev_mat$prob)

KM_dev <- summary(survival::survfit(Surv(dev_mat$days_to_follow_up, dev_mat$bcr_status) ~ 1), times = max(dev_mat$days_to_follow_up))
dev_log_OE_ratio <- log(KM_dev$n.event) - log(dev_expected_events)
dev_log_OE_ratio_se <- sqrt(KM_dev$surv / KM_dev$n.event)
dev_OE_ratio <- exp(dev_log_OE_ratio)
dev_OE_ratio_se <- sqrt(KM_dev$n.event * (1 - KM_dev$n.event / dev_expected_events) / dev_expected_events)
dev_OE_ratio_lower <- exp(dev_log_OE_ratio - (z_val * dev_log_OE_ratio_se))
dev_OE_ratio_upper <- exp(dev_log_OE_ratio + (z_val * dev_log_OE_ratio_se))

print("Calibration slope in the large")
cat(c(round(dev_OE_ratio, 4), round(dev_OE_ratio_lower, 4), round(dev_OE_ratio_upper, 4)))

form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coef(dev_fit)), collapse = "+"), " + offset(lp)"))
fit_dev <- coxph(formula = form, dev_mat)

print("ChiSq model fit")
cat(c(round(2*(diff(fit_dev$loglik)), 2), round(1-pchisq(2*(diff(fit_dev$loglik)), 8), 3)))

dev_c <- survcomp::concordance.index(x = dev_mat$lp, surv.time = dev_mat$days_to_follow_up, surv.event = dev_mat$bcr_status)
print("Harrells C")
cat(c(round(dev_c$c.index, 4), round(dev_c$lower, 4), round(dev_c$upper, 4)))

dev_mat$group <- as.numeric(dev_mat$lp > median(dev_mat$lp))
fit_dev <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = dev_mat)
dev_mat$group <- as.factor(dev_mat$group)
dev_mat$group <- relevel(dev_mat$group, ref = "0")
fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = dev_mat)

print("Hazard Ratio")
cat(c(round(as.numeric(coef(fit_dev)), 4), round(as.numeric(summary(fit_dev)$coef[,3]), 4),  paste(round(confint(fit_dev)[,1], 2), round(confint(fit_dev)[,2], 2), sep = "-"), summary(fit_dev)$coef[,5]))

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

print("AUC values")
cat(c(round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
      round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
      round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-")))



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

ext_mat <- ext_mat[, c(predictors, "days_to_follow_up", "bcr_status")]

# External, unseen dataset results
ext_fit <- coxph(Surv(days_to_follow_up, bcr_status) ~ ., ext_mat)

# dev fit is coxph on all 9 datasets above
ext_mat$lp <- gen_lp(dev_fit, ext_fit, "val")
ext_mat$prob <- predict(ext_fit, ext_mat, "risk")


fit_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, ext_mat)
print("Calibration slope")
cat(c(round(as.numeric(coef(fit_ext)), 4), round(confint(fit_ext)[,1], 4), round(confint(fit_ext)[,2], 4), summary(fit_ext)$coef[,5]))

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
cat(c(round(ext_OE_ratio, 4), round(ext_OE_ratio_lower, 4), round(ext_OE_ratio_upper, 4)))

form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coef(ext_fit)), collapse = "+"), " + offset(lp)"))
fit_ext <- coxph(formula = form, ext_mat)

print("ChiSq model fit")
cat(c(round(2*(diff(fit_ext$loglik)), 2), round(1-pchisq(2*(diff(fit_ext$loglik)), 8), 3)))

ext_c <- survcomp::concordance.index(x = ext_mat$lp, surv.time = ext_mat$days_to_follow_up, surv.event = ext_mat$bcr_status)
print("Harrells C")
cat(c(round(ext_c$c.index, 4), round(ext_c$lower, 4), round(ext_c$upper, 4)))

ext_mat$group <- as.numeric(ext_mat$lp > median(ext_mat$lp))
fit_ext <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = ext_mat)
ext_mat$group <- as.factor(ext_mat$group)
ext_mat$group <- relevel(ext_mat$group, ref = "0")
fit_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = ext_mat)

print("Hazard Ratio")
cat(c(round(as.numeric(coef(fit_ext)), 4), round(as.numeric(summary(fit_ext)$coef[,3]), 4),  paste(round(confint(fit_ext)[,1], 2), round(confint(fit_ext)[,2], 2), sep = "-"), summary(fit_ext)$coef[,5]))

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
cat(c(round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
      round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
      round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-")))

