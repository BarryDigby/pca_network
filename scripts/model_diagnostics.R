#!/usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(survival)
library(survminer)
library(pROC)
library(glmnet)
library(caret)
library(randomForestSRC)
library(CoxBoost)
library(MASS)
library(RegParallel)
library(timeROC)
library(car)
library(RegParallel)
library(pheatmap)


##########################################################################
# load the ceRNA network
##########################################################################
network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
genes = unique(network$mrna)

#########################################################################
# load TCGA metadata from Protein Atlas.
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
atlas_meta = merge(atlas_meta, pca_db, by.x="Row.names", by.y="sample_id")
keep = atlas_meta$sample_type=="Primary Tumor"
atlas_meta = atlas_meta[keep,]
# add Resectional status
# R0 no residual, R1 micro, R2 macro, RX uncertain.
resec = read.csv("~/Downloads/prad_tcga_clinical_data.tsv", header=T, sep="\t")
resec$Sample.ID = paste(resec$Sample.ID, "A", sep="")
table(resec$Sample.ID %in% atlas_meta$sample)
atlas_meta = merge(atlas_meta, subset(resec, select=c(Sample.ID, Surgical.Margin.Resection.Status)), by.x="sample", by.y="Sample.ID")
atlas_meta$years_to_follow_up <- as.character(floor(atlas_meta$years_to_follow_up))

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

#########################################################################
# Univariate Cox proportional hazards regression:
# Find ceRNA mRNAs ** with prognosis of PCa (Disease-free survival)
#########################################################################

# k blghe way 
univ_mat = as.data.frame(t(scaled))
univ_mat = cbind(univ_mat, atlas_meta[,c("days_to_follow_up", "bcr_status")])

res <- RegParallel(
  data = univ_mat,
  formula = 'Surv(days_to_follow_up, bcr_status) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(univ_mat)[1:256],
  blocksize = 50,
  p.adjust = "BH",
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)

res <- res[!is.na(res$P),]


# filter res DF
res = res[which(res$P.adjust <= 0.01),]

### check correlation of literature reviewed genes
lit_genes = c("TSPAN1", "REG4", "DDAH1", "CACNB3", "CTHRC1", "KIFC2", "INPP5E")
cor_mat = scaled[res$Variable,]
corr= cor(t(cor_mat))
library(viridis)
pheatmap(corr, display_numbers = F, cluster_rows = T, cluster_cols = T, scale = "none", show_colnames = F, col = hcl.colors(100, "Purple-Green",rev=T))

##########################################################################
# Stage data for LASSO Cox:
# Use normalised scaled and centered logcpm STAR counts data for modelling
# need to convert ENSG version IDs to HGNC symbols
# match samples in metadata
#
# ! The final line makes sure the expression data matches the metadata.
# ! The results are not reproducible if this is altered.
##########################################################################

scaled2 = scaled[which(rownames(scaled) %in% Active.Genes),]

mult = as.data.frame(t(scaled2))
mult = cbind(mult, atlas_meta$days_to_follow_up, atlas_meta$bcr_status)
colnames(mult)[17:18] <- c("days_to_follow_up", "bcr_status")

obj <- rfsrc(Surv(days_to_follow_up, bcr_status) ~ ., mult, splitrule = "logrank",
             ntree = 700, nodesize = 15, nsplit = 10, importance = TRUE)

o <- vimp(obj)
plot(o)


# cox boost
cph_boost = CoxBoost(mult$days_to_follow_up, mult$bcr_status, as.matrix(mult[1:9]))
cph_boost
cph_coefs = coefficients(cph_boost)
cph_coefs = coefficients(cph_boost)[cph_coefs]
cph_coefs[order(cph_coefs, decreasing = T)]

# overlapping coefs
overlap_genes = Reduce(intersect, list(names(cph_coefs), names(o$importance)[o$importance>=0.01]))

##########################################################################
# LASSO Cox
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
set.seed(11984191)
cv.fit <- cv.glmnet(x, y, family="cox", alpha=1, type.measure = "deviance", iter=1000)
fit = glmnet(x, y, family = "cox", alpha=1, iter=1000, lambda = NULL)
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Genes <- Coefficients@Dimnames[[1]][Active.Index]


cox = as.data.frame(t(scaled[res$Variable,]))

formula = Surv(atlas_meta$days_to_follow_up, atlas_meta$bcr_status) ~ .
model = StepReg::stepwiseCox(formula=formula, data=cox, selection = "bidirection", select="HQ", method = "efron")





cox = cbind(cox, atlas_meta[,c("days_to_follow_up", "bcr_status")])
cox_model = coxph(Surv(days_to_follow_up, bcr_status) ~ ., data=cox)



my.data$status1 <- ifelse(my.data$status==2,1,0)
data <- my.data
formula = Surv(time, status1) ~ . - status 

stepwiseCox(formula=formula,
            data=my.data,
            selection="bidirection",
            select="HQ",
            method="efron")

aic = stepAIC(cox_model, direction = "both")

summary(aic)

aic_coef = aic$coefficients[c(1:4,6:8)]
aic_coef

Active.Genes = names(aic_coef)
Active.Coefficients = as.numeric(aic_coef)

plot(cv.fit)
lbs_fun <- function(fit, ...) {
  L <- length(fit$lambda)
  x <- log(fit$lambda[L])
  y <- fit$beta[, L]
  labs <- names(y)
  text(x, y, labels=labs, ...)
}
# plot
plot(fit, xvar="lambda")
abline(v=c(log(cv.fit$lambda.1se), log(cv.fit$lambda.min)), lty=2)
# label
lbs_fun(fit)


cox = as.data.frame(scaled[Active.Genes,])
cox = as.data.frame(t(cox))
risk = rowSums(cox*Active.Coefficients)
cox$risk_score = risk
cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
cox = cbind(cox, atlas_meta[,c("days_to_follow_up", "bcr_status")])
res2 = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
res2_p = surv_pvalue(res2)
fuck = coxph(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
cox$bcr_status = as.numeric(cox$bcr_status)
time_intervals <- c(365, floor(365*3), floor(365*5))
ROC.train <- timeROC(T = cox$days_to_follow_up,
                     delta = cox$bcr_status,
                     marker = cox$risk_score,
                     cause = 1,
                     weighting = "marginal",
                     times = time_intervals,
                     iid = TRUE)
train_aucs = ROC.train$AUC


# use model coefs to initiate residuals.coxph

coxL <- coxph(Surv(atlas_meta$days_to_follow_up, atlas_meta$bcr_status) ~ . - SEC31B - SRRD, data=cox) 

test.ph = cox.zph(coxL)

ggcoxzph(test.ph, conf)
ggcoxzph2(test.ph)

plot(test.ph[1])

ggcoxdiagnostics(coxL, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())

ggcoxfunctional(Surv(atlas_meta$days_to_follow_up, atlas_meta$bcr_status) ~ REG4 + KIFC2 + CTHRC1 + SRRD + CACNB3 + TSPAN1 + INPP5E + DDAH1 + SEC31B, data = cox)
##########################################################################
# create a balanced partition over 1, 3 and 5 years
##########################################################################

colnames(atlas_meta)[duplicated(colnames(atlas_meta))] <- paste0(colnames(atlas_meta)[duplicated(colnames(atlas_meta))], "_duplicate")
# Create a new variable for follow-up period
clinical_data <- atlas_meta %>%
  mutate(follow_up_period = case_when(
    years_to_follow_up >= 0 & years_to_follow_up < 1 ~ "0-1",
    years_to_follow_up >= 1 & years_to_follow_up < 3 ~ "1-3",
    years_to_follow_up >= 3 & years_to_follow_up <= 5 ~ "3-5",
    TRUE ~ NA_character_
  ))

create_balanced_partition <- function(data) {
  balanced_indices <- createDataPartition(
    y = data$bcr_status,
    times = 1,
    p = 1,
    list = FALSE
  )
  return(data[balanced_indices, ])
}

balanced_data <- clinical_data %>%
  group_by(follow_up_period) %>%
  do(create_balanced_partition(.))

stop=FALSE
seed=0
for(i in 11060:500000){
if(stop){print("finished"); break}
seed=i
set.seed(seed)
partition_indices <- createDataPartition(
  y = rep(1:2, each = (nrow(balanced_data)/2)),
  times = 1,
  p = 0.5,
  list = FALSE
)

training_data <- balanced_data[partition_indices, ]
test_data <- balanced_data[-partition_indices, ]

scaled2 = scaled
train_mat = scaled2[, training_data$sample]
test_mat =  scaled2[, test_data$sample]

  sub_meta = subset(training_data, select=c(days_to_follow_up, bcr_status))
  colnames(sub_meta) = c("time", "status")
  x = as.matrix(t(train_mat[lit_genes,]))
  y = sub_meta
  y = as.matrix(y)
  # want it to bounce around a little insdie the training dataset
  for(j in 1:10){
  cv.fit <- cv.glmnet(x, y, family="cox", alpha=1, iter=1000, type.measure = "deviance", lambda=NULL, standardize=F) # data already scaled and centered
  fit = glmnet(x, y, family = "cox", iter=1000,alpha=1,lambda=NULL, standardize=F)
  Coefficients <- coef(fit, s = cv.fit$lambda.1se)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients <- Coefficients[Active.Index]
  Active.Genes <- Coefficients@Dimnames[[1]][Active.Index]
  if(length(Active.Genes) > 1){
    cox = as.data.frame(x[,Active.Genes])
    risk = rowSums(cox*Active.Coefficients)
    cox$risk_score = risk
    cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
    cox = cbind(cox, training_data[,c("days_to_follow_up", "bcr_status")])
    res2 = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
    res2_p = surv_pvalue(res2)
    cox$bcr_status = as.numeric(cox$bcr_status)
    time_intervals <- c(365, floor(365*3), floor(365*5))
    ROC.train <- timeROC(T = cox$days_to_follow_up,
                         delta = cox$bcr_status,
                         marker = cox$risk_score,
                         cause = 1,
                         weighting = "marginal",
                         times = time_intervals,
                         iid = TRUE)
    train_aucs = ROC.train$AUC
    # same but for coxboost
    # cox = as.data.frame(x[,cox_genes])
    # risk = rowSums(cox*cox_coef)
    # cox$risk_score = risk
    # cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
    # cox = cbind(cox, train_meta[,c("days_to_follow_up", "bcr_status")])
    # res21 = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
    # res21_p = surv_pvalue(res21)
    # cox$bcr_status = as.numeric(cox$bcr_status)
    # time_intervals <- c(365, floor(365*3), floor(365*5))
    # ROC.train1 <- timeROC(T = cox$days_to_follow_up,
    #                      delta = cox$bcr_status,
    #                      marker = cox$risk_score,
    #                      cause = 1,
    #                      weighting = "marginal",
    #                      times = time_intervals,
    #                      iid = TRUE)
    # train_aucs1 = ROC.train1$AUC
    # test
    x1 = as.matrix(t(test_mat[lit_genes,]))
    cox1 = as.data.frame(x1[,Active.Genes])
    risk = rowSums(cox1*Active.Coefficients)
    cox1$risk_score = risk
    cox1$risk_category = ifelse(cox1$risk_score > median(cox1$risk_score), "high", "low")
    cox1 = cbind(cox1, test_data[,c("days_to_follow_up", "bcr_status")])
    res3 = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox1)
    res3_p = surv_pvalue(res3)
    cox1$bcr_status = as.numeric(cox1$bcr_status)
    ROC.test <- timeROC(T = cox1$days_to_follow_up,
                         delta = cox1$bcr_status,
                         marker = cox1$risk_score,
                         cause = 1,
                         weighting = "marginal",
                         times = time_intervals,
                         iid = TRUE)
    test_aucs = ROC.test$AUC

    
    x2 = as.matrix(t(scaled[lit_genes,]))
    cox3 = as.data.frame(x2[,Active.Genes])
    risk = rowSums(cox3*Active.Coefficients)
    cox3$risk_score = risk
    cox3$risk_category = ifelse(cox3$risk_score > median(cox3$risk_score), "high", "low")
    cox3 = cbind(cox3, atlas_meta[,c("days_to_follow_up", "bcr_status")])
    res31 = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox3)
    res31_p = surv_pvalue(res31)
    cox3$bcr_status = as.numeric(cox3$bcr_status)
    ROC.full <- timeROC(T = cox3$days_to_follow_up,
                        delta = cox3$bcr_status,
                        marker = cox3$risk_score,
                        cause = 1,
                        weighting = "marginal",
                        times = time_intervals,
                        iid = TRUE)
    full_aucs = ROC.full$AUC
    
    if(res2_p$pval < 0.05 & res3_p$pval < 0.05 & res31_p$pval < 0.05){
    print(paste("Run number:", i))
    print(paste("Model number:", j))
    cat("LASSO Cox genes:", paste(Active.Genes), "\n")
    cat("LASSO-Cox Training logrank:", paste(res2_p$pval), "\n")
    cat("LASSO-Cox Testing logrank:", paste(res3_p$pval), "\n")
    cat("LASSO-Cox Full dataset logrank:", paste(res31_p$pval), "\n")
    cat("LASSO-Cox Training ROC:", paste(train_aucs), "\n")
    cat("LASSO-Cox Testing ROC:", paste(test_aucs), "\n")
    cat("LASSO-Cox Full dataset ROC:", paste(full_aucs), "\n")
    }
    # print("===================================================")
    # cat("CoxBoost genes:", paste(cox_genes), "\n")
    # cat("CoxBoost Training logrank:", paste(res21_p$pval), "\n")
    # cat("CoxBoost Testing logrank:", paste(res31_p$pval), "\n")
    # cat("CoxBoost Training ROC:", paste(train_aucs1), "\n")
    # cat("CoxBoost Testing ROC:", paste(test_aucs1), "\n")
    auc_cutoff = 0.65
    if(res2_p$pval < 0.05 & res3_p$pval < 0.05 & res31_p$pval < 0.05 &
       train_aucs[1] > auc_cutoff & train_aucs[2] > auc_cutoff & train_aucs[3] > auc_cutoff &
       test_aucs[1]> auc_cutoff & test_aucs[2] > auc_cutoff & test_aucs[3] > auc_cutoff &
       full_aucs[1] >auc_cutoff & full_aucs[2] > auc_cutoff & full_aucs[3] > auc_cutoff){
       # (res21_p$pval < 0.05 & res31_p$pval < 0.05 &
       #  train_aucs1[1] > 0.69 & train_aucs1[2] > 0.69 & train_aucs1[3] > 0.69 &
       #  test_aucs1[1]> 0.69 & test_aucs1[2] > 0.69 & test_aucs1[3] > 0.69)) {
      filename <- paste("~/Desktop/lessCACNB3_", paste(Active.Genes, collapse = "-"), "_", seed, ".RData", sep = "")
      save.image(file = filename)
      print("! File created !")
      stop=FALSE
      break
    }
    if(stop){break}
  }
  if(stop){break}
}
if(stop){break}
}







lit_genes_1se = read.csv("~/Desktop/LASSO_runs/lit_genes_lambda.1se.txt", header=T, sep="\t")
colnames(lit_genes_1se) = c("iteration", "model_no", "genes", "train", "test", "full")









load("~/Desktop/CTHRC1-CACNB3-KIFC2_11064.RData")

mat = cox # train
mat = cox1 # test
mat = cox3 # full

ggsurvplot(res31,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=mat,
           xlab="Time (days)")

ROC.1 <- timeROC(T=mat$days_to_follow_up,
                 delta=mat$bcr_status,
                 marker=mat$risk_score,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)
ROC.2  <- timeROC(T=mat$days_to_follow_up,
                  delta=mat$bcr_status,
                  marker=mat$risk_score,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=mat$days_to_follow_up,
                 delta=mat$bcr_status,
                 marker=mat$risk_score,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

plot(ROC.1, time=365, title = F, lwd=2)
plot(ROC.2, time=floor(365*3), col="blue", add = T, title=F, lwd=2)
plot(ROC.3, time=floor(365*5), col="forestgreen", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(ROC.1$AUC[2],2)),
               paste0("3 year AUC: ", round(ROC.2$AUC[2],2)),
               paste0("5 year AUC:", round(ROC.3$AUC[2],2)))
legend("bottomright", legend = my_legend,col=c("red","blue", "forestgreen"),lwd=2, cex=1)
