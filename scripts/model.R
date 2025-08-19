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

## please use this function, not direct linear predictors for validation!!
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

elastic_net_cv <- function(x,y,group) {
  
  set.seed(123)
  alpha_seq = seq(0.1, 1, length = 10)
  
  lambda_max <- max(glmnet(x, y, family = "cox", alpha = 1)$lambda) 
  lambda_seq <- exp(seq(log(lambda_max), log(lambda_max / 100), length = 100))
  
  fitEN <- list()
  
  for (i in 1:length(alpha_seq)) {
    fitEN[[i]] <- cv.glmnet(x, y, family = "cox", alpha = alpha_seq[i], lambda = lambda_seq, foldid = group)
  }
  
  idx <- which.min(sapply(fitEN, function(xx) {
    xx$cvm[xx$lambda == xx$lambda.min]
  }))
  
  return(fitEN[[idx]])
}

lasso_cv <- function(x, y, group){
  lambda_max <- max(glmnet(x, y, family = "cox", alpha = 1)$lambda)
  lambda_seq <- exp(seq(log(lambda_max), log(lambda_max / 100), length = 100))
  
  lasso_fit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_seq,foldid = group)
  
  return(lasso_fit)
}

# dev mat here 
extract_lasso_elastic <- function(model, n){
  genes <- coef(model$glmnet.fit, s=model$lambda.min)[,1]
  genes <- names(genes)[abs(genes) > quantile(abs(genes))[n]]
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(genes, collapse = " + ")))
  fit <- coxph(form, dev_mat)
  return(fit)
}

cox_boost <- function(x, y){
  snowfall::sfInit(parallel=TRUE, cpus=6)
  cv.res <- cv.CoxBoost(time=y[,1], status=y[,2], x=x, maxstepno=100,parallel=TRUE,K=8,type="naive",penalty=200)
  cbfit <- CoxBoost(time=y[,1], status=y[,2], x=x, stepno=cv.res$optimal.step,penalty=200)
  return(cbfit)
}

# dev mat here 
extract_coxboost <- function(model, n){
  genes <- coef(model)
  genes <- names(genes)[abs(genes) > quantile(abs(genes))[n]]
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(genes, collapse = " + ")))
  fit <- coxph(form, dev_mat)
  return(fit)
}

set.seed(123)
load("/data/github/pca_network/reviewer4/internal_external_start_point.RData")

gene_mask <- colnames(mat)[which(!(colnames(mat) %in% c("preop_psa", "gleason_score", "days_to_follow_up", "bcr_status", "study")))]

results <- data.frame()

for(i in unique(mat$study)){
  print(i)
  dev_mat <- mat %>% dplyr::filter(!(study %in% i))
  val_mat <- mat %>% dplyr::filter(study %in% i)
  
  # Simple HR filter and loose Pval filter on the initial gene set.
  ## includes "study" as a random effect. 
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(gene_mask, collapse = " + "), " + preop_psa + gleason_score + strata(study)"))
  dev_fit <- coxph(fx, dev_mat)
  exp_filter <- names(exp(coef(dev_fit))[exp(coef(dev_fit)) > 1.1])
  exp_filter <- c(exp_filter, names(exp(coef(dev_fit))[exp(coef(dev_fit)) < 0.9]))
  # uncomment below if using coxme
  #coefficients <- fixef(dev_fit) 
  #se <- sqrt(diag(vcov(dev_fit))) 
  #z_scores <- coefficients / se
  #p_values <- 2 * pnorm(-abs(z_scores))
  p_values <- summary(dev_fit)$coef[,5]
  gene_pvalues <- p_values[p_values < 0.1]
  pass_filter <- intersect(names(gene_pvalues), exp_filter)
  
  # now construct filtered coxph model
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(pass_filter, collapse = " + ")))
  dev_coxph <- coxph(fx, dev_mat)
  
  # Export the models to global env for review
  assign(paste0(i, "_dev_coxph"), dev_coxph, envir=.GlobalEnv)
  
  # Pass the filtered genes to a LASSO & ENET model
  ## Stage lasso data 
  keep_patient <- rownames(dev_mat)[which(dev_mat$days_to_follow_up != 0)]
  y <- as.matrix(mat[keep_patient, c("days_to_follow_up", "bcr_status")])
  colnames(y) <- c("time", "status") 
  all_x <- model.matrix( ~ 0 + ., dev_mat[keep_patient, c(pass_filter)])
  group_identifier <- as.numeric(as.factor(dev_mat$study[which(rownames(dev_mat) %in% keep_patient)]))
  
  # Run LASSO and ENET
  lasso <- lasso_cv(all_x, y, group_identifier)
  elastic <- elastic_net_cv(all_x, y, group_identifier)
  coxboost <- cox_boost(all_x, y)
  
  # validation coxph in isolation
  fx <- as.formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coef(dev_coxph)), collapse = " + ")))
  val_coxph <- coxph(fx, val_mat)
  assign(paste0(i, "_val_coxph"), val_coxph, envir=.GlobalEnv)
  
  # now enter loop for quantile cutoffs for lasso elastic n coxboost
  for (coef_quantile in seq(1, 4)){
    dev_lasso <- extract_lasso_elastic(lasso, coef_quantile)
    dev_elastic <- extract_lasso_elastic(elastic, coef_quantile)
    dev_coxboost <- extract_coxboost(coxboost, coef_quantile)
  
    # Export the models to global env for review
    assign(paste0(i, "_dev_lasso_", coef_quantile), dev_lasso, envir=.GlobalEnv)
    assign(paste0(i, "_dev_elastic_", coef_quantile), dev_elastic, envir=.GlobalEnv)
    assign(paste0(i, "_dev_coxboost_", coef_quantile), dev_coxboost, envir=.GlobalEnv)
  
    # Now construct the validation models
    fx <- as.formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coef(dev_lasso)), collapse = " + ")))
    val_lasso <- coxph(fx, val_mat)
    fx <- as.formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coef(dev_elastic)), collapse = " + ")))
    val_elastic <- coxph(fx, val_mat)
    fx <- as.formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coef(dev_coxboost)), collapse = " + ")))
    val_coxboost <- coxph(fx, val_mat)
    
    # Export the models to global env for review
    assign(paste0(i, "_val_lasso_", coef_quantile), val_lasso, envir=.GlobalEnv)
    assign(paste0(i, "_val_elastic_", coef_quantile), val_elastic, envir=.GlobalEnv)
    assign(paste0(i, "_val_coxboost_", coef_quantile), val_coxboost, envir=.GlobalEnv)
  }
  
  # now enter loop to generate performance statistics:
  df <- data.frame(dev = grep(paste0(i, "_dev_"), ls(), value = TRUE), val = grep(paste0(i, "_val_"), ls(), value = TRUE))
  
  for (j in 1:nrow(df)){
  dev_name <- df[j, 1]
  val_name <- df[j, 2]
  dev_fit <- get(df[j, 1])
  val_fit <- get(df[j, 2])
  
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
  
  # Model Expected/Observed - same interpretation as calibration slope.
  # use model with genes + clinical covars, not lp. (dev_fit)
  #dev_expected <- sum(predict(dev_fit, dev_mat, type = "expected"))
  #dev_observed <- sum(dev_mat$bcr_status == 1)
  #dev_EO <- dev_expected / dev_observed
  ## SE
  #dev_var_cov <- vcov(dev_fit)
  #dev_dm <- model.matrix(dev_fit)
  #dev_var <- rowSums((dev_dm %*% dev_var_cov) * dev_dm)
  #dev_se <- sqrt(sum(dev_var))
  #
  #val_expected <- sum(predict(val_fit, val_mat, type = "expected"))
  #val_observed <- sum(val_mat$bcr_status == 1)
  #val_EO <- val_expected / val_observed
  ## SE
  #val_var_cov <- vcov(val_fit)
  #val_dm <- model.matrix(val_fit)
  #val_var <- rowSums((val_dm %*% val_var_cov) * val_dm)
  #val_se <- sqrt(sum(val_var))
  #
  #out_dev <- c(out_dev, c(round(dev_EO, 4), round(dev_se, 4)))
  #out_val <- c(out_val, c(round(val_EO, 4), round(val_se, 4)))
  
  # model fit/misspecification
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(names(coef(dev_fit)), collapse = "+"), " + offset(lp)"))
  fit_dev <- coxph(formula = form, dev_mat)
  fit_val <- coxph(formula = form, val_mat)

  out_dev <- c(out_dev, c(round(2*(diff(fit_dev$loglik)), 2), round(1-pchisq(2*(diff(fit_dev$loglik)), 8), 3)))
  out_val <- c(out_val, c(round(2*(diff(fit_val$loglik)), 2), round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3)))
  
  # Harrells Cindex - old way 
  #out_dev <- c(out_dev, c(round(as.numeric(rcorr.cens(-1*dev_mat$lp, Surv(dev_mat$days_to_follow_up, dev_mat$bcr_status))[1]),2)))
  #out_val <- c(out_val, c(round(as.numeric(rcorr.cens(-1*val_mat$lp, Surv(val_mat$days_to_follow_up, val_mat$bcr_status))[1]),2)))
  # newer way with variance, uncertain about variance produced here. 
  #fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, dev_mat)
  #fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, val_mat)
  #out_dev <- c(out_dev, round(concordance(fit_dev)$concordance, 4), round(concordance(fit_dev)$var, 4))
  #out_val <- c(out_val, round(concordance(fit_val)$concordance, 4), round(concordance(fit_val)$var, 4))
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
  
  # validation goes up to 6 years, rerun and use belfast for validation. 
  
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
  rownames(results)[nrow(results)] <- dev_name
  results <- rbind(results, t(as.data.frame(out_val)))
  rownames(results)[nrow(results)] <- val_name
  }
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
  df <- results %>% tibble::rownames_to_column(., "rownames") %>% dplyr::filter(grepl(paste0("_val_", model), rownames)) %>% tibble::column_to_rownames(., "rownames")
  df <- df %>% mutate_all(as.numeric)
  met <- metagen(TE = df[, effect], lower = df[, LCI], upper = df[, UCI], studlab = rownames(df))
  return(meta::forest(met, common = FALSE))
}


forest_plot("coxboost_4", "Slope Coef", "Slope LCI", "Slope UCI")

# elastic 4 has the best slope - 0.89. 
elastics <- grep("val_elastic_4", ls(), value = TRUE)
elastic_genes <- unlist(lapply(elastics, function(model_name) {
  model <- get(model_name)
  names(coef(model)) 
}))

ret_sig <- function(model){
  return(rownames(summary(model)$coef)[which(summary(model)$coef[,5] < 0.1)])
}

x <- c(ret_sig(TCGA_dev), ret_sig(Stockholm_dev), ret_sig(GSE54460_dev), ret_sig(CPC_dev), ret_sig(CMap_dev), ret_sig(Belfast_dev))

stable_sigs <- names(coef(DKFZ_val_coxboost_4))

fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(stable_sigs, collapse = " + ")))
fit <- coxph(fx, mat)
lp <- predict(fit, mat, type = "lp")

print(round(as.numeric(rcorr.cens(-1*lp, Surv(mat$days_to_follow_up, mat$bcr_status))[1]),2))


# calibration plot
dat <- mat[, c(stable_sigs, "days_to_follow_up", "bcr_status")]
dat$lp <- predict(fit, dat, type="lp")
d = rms::datadist(dat)
options(datadist="d")
units(d) = "Day"
f   <- cph(Surv(days_to_follow_up, bcr_status) ~ lp, data=dat, x=TRUE, y=TRUE, time.inc=2920, surv=TRUE)
cal <- calibrate(f, u=2920, B=300, conf.int = TRUE)
plot(cal, ylim=c(0, 1), xlim=c(0,1), subtitles=FALSE)
calkm <- calibrate(f, u=2920, cmethod='KM', B=3000)
plot(calkm, add = FALSE)




res <- c()
for(gene in gene_mask){
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", gene, " + preop_psa + gleason_score + (1|study)"))
  fit <- coxme::coxme(fx, mat)
  pv <- as.numeric(sub(".*\\s([0-9\\.]+)$", "\\1", capture.output(fit)[15]))
  names(pv) <- gene
  res <- c(res, pv)
}
univariate_filter <- res[res < 0.01]
univariate_filter <- names(univariate_filter)
univariate_filter <- unclass(na.omit(unique(c(univariate_filter, "study")))) # add these back for ML methods. 


### Random Forest Variable Importance

rf_mat <- mat #%>% dplyr::select(-"study")
rf_mat$study <- factor(rf_mat$study)
rf_imp <- rfsrc(Surv(days_to_follow_up, bcr_status) ~., rf_mat, ntree = 500, nsplit = 10, mtry = sqrt(223), nodesize = 3, importance = TRUE)
rf_imp_sample <- subsample(rf_imp, B = 100)
rf_imp_stats <- plot(rf_imp_sample)$stats
rf_importance <- colnames(rf_imp_stats)
rf_importance <- unique(c(rf_importance,"study"))

## Random Forest Minimal Depth

rf_md <- var.select(Surv(days_to_follow_up, bcr_status) ~ ., rf_mat, method = "md", ntree = 500, nsplit = 10, nodesize = 3, splitrule = "logrank")
rf_minimal_depth <- rf_md$topvars
rf_minimal_depth <- unique(c(rf_minimal_depth,"study"))



## Random Forest Variable Hunting

rf_vh <- var.select(Surv(days_to_follow_up, bcr_status) ~ ., rf_mat, method = "vh", ntree = 500, nrep = 25, nsplit = 10, K = 5, nodesize = 3, nstep =1, splitrule = "logrank")
rf_variable_hunter <- rf_vh$topvars
rf_variable_hunter <- unique(c(rf_variable_hunter, "study"))


## Minimum Redundancy Maximum Relevance Filter

common <- c()

for( i in unique(mat$study)){

  mr_mat <- mat %>% dplyr::filter(!(study %in% i))
  mr_mat2 <- mr_mat[, c(gene_mask, "study")]
  mr_mat2$study <- factor(mr_mat2$study, ordered = TRUE)
  mr_mat2$target <- Surv(mr_mat$days_to_follow_up, mr_mat$bcr_status)
  mr_data <- mRMR.data(data=mr_mat2)
  mr_classic <- mRMR.classic(data=mr_data, target_indices=as.integer(ncol(mr_mat2)), feature_count = 50)
  mr_genes <- mr_classic@feature_names[unlist(mRMRe::solutions(mr_classic))]
  common <- c(common, mr_genes)

}

mRMR_filter <- names(table(common)[table(common) >= 9])
mRMR_filter <- unique(c(mRMR_filter, c("study")))

save(list = c("univariate_filter", "mRMR_filter", "rf_variable_hunter", "rf_importance", "rf_minimal_depth"), file = "/data/github/pca_network/reviewer4/master/pre_filtering.RData")
print("saved filtered gene sets")

#load("/data/github/pca_network/rev4/belfast/pre_filtering.RData")
# Machine Learning Algorithms

# construct dataframe and X Y for glmnet
all_mat <- mat[, c(network_genes)]
univ_mat <- mat[, c(univariate_filter, "days_to_follow_up", "bcr_status")]
vimp_mat <- mat[, c(rf_importance, "days_to_follow_up", "bcr_status")]
vh_mat <- mat[, c(rf_variable_hunter, "days_to_follow_up", "bcr_status")]
md_mat <- mat[, c(rf_minimal_depth, "days_to_follow_up", "bcr_status")]
mrmr_mat <- mat[, c(mRMR_filter, "days_to_follow_up", "bcr_status")]

# identify samples with time to follow up == 0
keep_patient <- rownames(dev_mat)[which(dev_mat$days_to_follow_up != 0)]
y <- as.matrix(mat[keep_patient, c("days_to_follow_up", "bcr_status")])
colnames(y) <- c("time", "status") 

all_x <- model.matrix( ~ 0 + ., dev_mat[keep_patient, c(pass_filter)])
univ_x <- model.matrix( ~ 0 + ., mat[keep_patient, c(univariate_filter)])
vimp_x <- model.matrix( ~ 0 + ., mat[keep_patient, c(rf_importance)])
vh_x <- model.matrix(~ 0 + ., mat[keep_patient, c(rf_variable_hunter)])
md_x <- model.matrix(~0 + ., mat[keep_patient, c(rf_minimal_depth)])
mrmr_x <- model.matrix( ~ 0 + ., mat[keep_patient, c(mRMR_filter)])

## Multivariate Coxph
formula_builder <- function(x){
  vec <- gsub("study", "(1|study)", x)
  formula_str <- paste("Surv(days_to_follow_up, bcr_status) ~", paste(vec, collapse = " + "))
  formula <- as.formula(formula_str)
  return(formula)
}
multivariate_cox_all <- coxme::coxme(formula_builder(c(gene_mask, "gleason_score", "preop_psa", "study")), all_mat, x = TRUE)
multivariate_cox_univ <- coxme::coxme(formula_builder(univariate_filter), univ_mat, x = TRUE)
multivariate_cox_vimp <- coxme::coxme(formula_builder(rf_importance), vimp_mat, x = TRUE)
multivariate_cox_vh <-  coxme::coxme(formula_builder(rf_variable_hunter), vh_mat, x = TRUE)
multivariate_cox_md <- coxme::coxme(formula_builder(rf_minimal_depth), md_mat, x = TRUE)
multivariate_cox_mrmr <-  coxme::coxme(formula_builder(mRMR_filter), mrmr_mat, x = TRUE)
print("multivariate models done")

## GLMnet Lasso

group_identifier <- as.numeric(as.factor(dev_mat$study[which(rownames(dev_mat) %in% keep_patient)]))
pf <- ifelse(grepl("gleason|preop|study", colnames(all_x)), 0, 1)
lasso_all <- glmnet::cv.glmnet(all_x, y, family = "cox", nfolds = 8, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop|study", colnames(univ_x)), 0, 1)
lasso_univ <- glmnet::cv.glmnet(univ_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop|study", colnames(vimp_x)), 0, 1)
lasso_vimp <- glmnet::cv.glmnet(vimp_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop|study", colnames(vh_x)), 0, 1)
lasso_vh <- glmnet::cv.glmnet(vh_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop|study", colnames(md_x)), 0, 1)
lasso_md <- glmnet::cv.glmnet(md_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
pf <- ifelse(grepl("gleason|preop|study", colnames(mrmr_x)), 0, 1)
lasso_mrmr <- glmnet::cv.glmnet(mrmr_x, y, family = "cox", foldid = group_identifier, penalty.factor = pf)
print("Lasso models done")

## GLMnet Elastic Net

elastic_net_cv <- function(x) {

  set.seed(123)
  alpha_seq = seq(0.1, 1, length = 10)
  fitEN <- list()
  pf <- ifelse(grepl("gleason|preop|study", colnames(x)), 0, 1)
  
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
  snowfall::sfInit(parallel=TRUE, cpus=6)
  cv.res <- cv.CoxBoost(time=y[,1], status=y[,2], x=x, maxstepno=500,parallel=TRUE,K=8,type="verweij",penalty=100)
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

save_models <- ls(.GlobalEnv)[grep("^(boost_|elastic_|lasso_|multivariate_)", ls(.GlobalEnv))]
save(list = save_models, file="/data/github/pca_network/reviewer4/master/models.RData")

# Model Evaluation

# Lasso/Enet extraction
extract_lasso <- function(model){
  genes <- coef(model$glmnet.fit, s=model$lambda.min)[,1]
  genes <- names(genes)[genes != 0]
  #genes <- c(genes[!grepl("^study", genes)], "(1|study)")
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(genes, collapse = " + ")))
  fit <- coxph(form, dev_mat)
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
  genes <- c(genes[!grepl("^study", genes)], "(1|study)")
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~", paste(genes, collapse = " + ")))
  fit <- coxme::coxme(form, mat, x=TRUE)
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

for(row in 1:nrow(models_used)){

    model_name <- models_used[row, 1]
    model <- get(models_used[row, 1])
    genes <- names(coef(model))
    
    df_dev <- mat[, c(genes, "study", "days_to_follow_up", "bcr_status")]
    df_val <- mat[, c(genes, "days_to_follow_up", "bcr_status")]

    dev_lp <- ehahelper::predict_coxme(model, df_dev, type="lp")
    val_lp <- ehahelper::predict_coxme(model, df_val, type="lp")
    
    df_dev$lp <- dev_lp
    df_val$lp <- val_lp
    
    # model slope
    fit_dev <- coxme::coxme(Surv(days_to_follow_up, bcr_status) ~ lp + (1|study), df_dev)
    fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, df_val)
    
    out_dev <- c(round(as.numeric(coef(fit_dev)), 3), paste(round(confint(fit_dev)[,1], 2), round(confint(fit_dev)[,2], 2), sep = "-"), as.numeric(sub(".*\\s([0-9\\.]+)$", "\\1", capture.output(fit_dev)[15])))
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
    
    # AUC
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
    
    # validation goes up to 6 years, rerun and use belfast for validation. 
    
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

write.table(results, "/data/github/pca_network/reviewer4/master/results.txt", quote=F, sep="\t")
