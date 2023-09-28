#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)
library(survminer)
library(survival)
library(rms)
library(survivalAnalysis)
library(pROC)
library(survivalROC)
library(ggpubr)
library(gtsummary)
library(RColorBrewer)

#########################################################################
# load model from previous script.
#########################################################################
load("/data/github/pca_network/results/TCGA_DFS/model.RData")

#########################################################################
# Load cleaned tcga meta, make sure resectional status is included
#########################################################################
load("/data/github/pca_network/data/TCGA_meta_cleaned.RData")
resec = read.csv("/data/github/pca_network/data/prad_tcga_clinical_data.tsv", header=T, sep="\t")
resec$Sample.ID = paste(resec$Sample.ID, "A", sep="")
table(resec$Sample.ID %in% master$Row.names)
master = merge(master, subset(resec, select=c(Sample.ID, Surgical.Margin.Resection.Status)), by.x="Row.names", by.y="Sample.ID")

#########################################################################
# Construct dataframe for analysis
# Tidy/create categoricals for analysis
#########################################################################

mat = master

age_cut = survminer::surv_cutpoint(mat, time="days_to_follow_up", event="bcr_status", variables="age")
mat$age = ifelse(mat$age >= age_cut$cutpoint$cutpoint, ">56", "<56")
mat = tibble::column_to_rownames(mat, "Row.names")
mat$psa <- ifelse(mat$preop_psa > 20, ">20ng/ml",
                  ifelse(mat$preop_psa >= 10 & mat$preop_psa <= 20, "10-20ng/ml", "<10ng/ml"))

mat$path_t = ifelse(mat$ajcc_pathologic_t == "T2a" | mat$ajcc_pathologic_t == "T2b" | mat$ajcc_pathologic_t == "T2c", "T2", mat$ajcc_pathologic_t)
mat$path_t = ifelse(mat$ajcc_pathologic_t == "T3a" | mat$ajcc_pathologic_t == "T3b", "T3", mat$path_t)
mat$clin_m = ifelse(mat$ajcc_clinical_m == "M0", "M0", "M1")
mat$path_n = ifelse(mat$ajcc_pathologic_n == "N1", "N1","N0")
mat$gleason = ifelse(mat$gleason_score == 6, "Gleason 6 or lower", mat$gleason_score)
mat$gleason = ifelse(mat$gleason_score == 7, "Gleason 7", mat$gleason)
mat$gleason = ifelse(mat$gleason_score > 7, "Gleason 8,9,10", mat$gleason)
mat$surgical_r = ifelse(mat$Surgical.Margin.Resection.Status == "R2" | mat$Surgical.Margin.Resection.Status == "R1",  "R1/R2", mat$Surgical.Margin.Resection.Status)


w <- transcan(~ age + gleason + psa + path_t + path_n + clin_m +
                surgical_r, imputed=TRUE, data=mat, pl=TRUE, pr=FALSE)


attach(mat)
mat2 = mat
psa = impute(w, psa, data=mat)
levels(factor(mat$psa))
mat2$psa = ifelse(psa == 1, "<10ng/ml",
                  ifelse(psa == 3, "10-20ng/ml",
                         ifelse(psa == 2, ">20ng/ml", psa)))

path_t = impute(w, path_t, data = mat)
levels(factor(mat$path_t))
mat2$path_t = ifelse(path_t == 1, "T2",
                     ifelse(path_t == 2, "T3", 
                            ifelse(path_t == 3, "T4", path_t)))


path_n = impute(w, path_n, data=mat)
levels(factor(mat$path_n))
mat2$path_n = ifelse(path_n == 1, "N0",
                     ifelse(path_n == 2, "N1", path_n))


clin_m = impute(w, clin_m, data=mat)
levels(factor(mat$clin_m))
mat2$clin_m = ifelse(clin_m == 1, "M0", clin_m)


surgical_r = impute(w, surgical_r, data=mat)
levels(factor(mat$surgical_r))
mat2$surgical_r = ifelse(surgical_r==1, "R0",
                         ifelse(surgical_r == 2, "R1/R2",
                                ifelse(surgical_r == 3 , "RX", surgical_r)))


mat2$risk_score = risk_scores

## univariate forest plot
vars_for_table = c("age","psa", "gleason", "path_t", "path_n", "clin_m", "surgical_r", "risk_score")
univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(days_to_follow_up, bcr_status)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mat2)})
forestmodel::forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)


#pdf("/data/github/pca_network/results/TCGA_DFS3/univariate_forest_model.pdf", width=8, height=5)
forestmodel::forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)

forestmodel::forest_model(coxph(Surv(days_to_follow_up, bcr_status) ~ age + psa + gleason +  path_t + path_n + clin_m + surgical_r + risk_score, data=mat2))

clin_mod = coxph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2)

# Harrell's C
rcorr.cens(-1*mat2$risk_score, Surv(mat2$days_to_follow_up, mat2$bcr_status))[1]


# nomogram doesnt save plot, save as pdf from window
pdf("/data/github/pca_network/results/TCGA_DFS/nomogram.pdf", width=8, height=6)
regplot::regplot(clin_mod, failtime = c(365,1095,1826, 3650), points = T, odds = FALSE, nsamp = 481, title = "")
dev.off()
# calibrate nomogram

d = rms::datadist(mat2)
options(datadist="d")
units(d) = "Day"
rms_cox = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356)
rms_cox2 = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*3)
rms_cox3 = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*5)
rms_cox4 = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)

cal1yr <- calibrate(rms_cox, B=1000, u=365*1, maxdim=3, conf.int=TRUE, cmethod = "KM")
cal3yr <- calibrate(rms_cox2, B=1000, u=365*3, maxdim=3, conf.int=TRUE, cmethod = "KM")
cal5yr <- calibrate(rms_cox3, B=1000, u=365*5, maxdim=3, conf.int=TRUE, cmethod = "KM")
cal8yr <- calibrate(rms_cox4, B=1000, u=365*8, maxdim=3, conf.int=TRUE, cmethod = "KM")

pdf("/data/github/pca_network/results/TCGA_DFS/nomogram_calibration.pdf", width=7, height=7)
plot(cal1yr, xlim=c(0,1), ylim=c(0,1), col="black", lwd=2)
plot(cal3yr,  xlim=c(0,1), ylim=c(0,1), col="green", add=T, lwd=2)
plot(cal5yr,  xlim=c(0,1), ylim=c(0,1), col="blue", add=T, lwd=2)
plot(cal8yr,  xlim=c(0,1), ylim=c(0,1), col="red", add=T, lwd=2)
my_legend =  c(paste0("1 year"),
               paste0("3 year"),
               paste0("5 year"),
               paste0("8 year"))
legend("bottomright", legend = my_legend,col=c("black","green", "blue", "red"),lwd=2, cex=1)
dev.off()


# nomogram ROC
time_points <- c(365, 3*365, 5*365, 8*365)

# Create ROC curves for each time point
roc_curves <- lapply(time_points, function(time) {
  # Subset data to include only observations within the specified time
  subset_data <- mat2[mat2$days_to_follow_up <= time, ]
  roc_obj <- roc(subset_data$bcr_status, predict(rms_cox4, newdata = subset_data))
  return(roc_obj)
})


pdf("/data/github/pca_network/results/TCGA_DFS/nomogram_roc_1358.pdf", width=7, height=6)
plot(roc_curves[[1]], time=365, col="black", title = F, lwd=2)
plot(roc_curves[[2]], time=floor(365*3), col="green", add = T, title=F, lwd=2)
plot(roc_curves[[3]], time=floor(365*5), col="blue", add = T, title=F, lwd=2)
plot(roc_curves[[4]], time=floor(365*8), col="red", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(roc_curves[[1]]$auc,2)),
               paste0("3 year AUC: ", round(roc_curves[[2]]$auc,2)),
               paste0("5 year AUC:", round(roc_curves[[3]]$auc,2)),
               paste0("8 year AUC:", round(roc_curves[[4]]$auc,2)))
legend("bottomright", legend = my_legend,col=c("black","green", "blue", "red"),lwd=2, cex=1)
dev.off()

# nom vs individual predictors 1 year
# changing time.inc does not matter, we specify time when it comes to plotting.
age = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
path_t = rms::cph(Surv(days_to_follow_up, bcr_status) ~  path_t, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
surgical_r = rms::cph(Surv(days_to_follow_up, bcr_status) ~ surgical_r, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
risk_score = rms::cph(Surv(days_to_follow_up, bcr_status) ~ risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
nom = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)

models <- c("age", "path_t", "surgical_r", "risk_score", "nom")

roc_curves <- lapply(models, function(model) {
  # Subset data to include only observations within the specified time
  subset_data <- mat2[mat2$days_to_follow_up <= 365, ]
  roc_obj <- roc(subset_data$bcr_status, predict(get(model), newdata = subset_data))
  return(roc_obj)
})

pdf("/data/github/pca_network/results/TCGA_DFS/nomogram_v_rest_1year.pdf", width=7, height=6)
plot(roc_curves[[1]], time=365, col="black", title = F, lwd=2)
plot(roc_curves[[2]], time=floor(365), col="yellow", add = T, title=F, lwd=2)
plot(roc_curves[[3]], time=floor(365), col="green", add = T, title=F, lwd=2)
plot(roc_curves[[4]], time=floor(365), col="blue", add = T, title=F, lwd=2)
plot(roc_curves[[5]], time=floor(365), col="red", add = T, title=F, lwd=2)
my_legend =  c(paste0("Age AUC: ", round(roc_curves[[1]]$auc,2)),
               paste0("Pathological T AUC: ", round(roc_curves[[2]]$auc,2)),
               paste0("Resection status AUC: ", round(roc_curves[[3]]$auc,2)),
               paste0("Prognostic index AUC: ", round(roc_curves[[4]]$auc,2)),
               paste0("Nomogram AUC: ", round(roc_curves[[5]]$auc,2)))
legend("bottomright", legend = my_legend,col=c("black","yellow", "green", "blue", "red"),lwd=2, cex=1)
dev.off()

# nom vs individual predictors 3 year
age = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
path_t = rms::cph(Surv(days_to_follow_up, bcr_status) ~  path_t, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
surgical_r = rms::cph(Surv(days_to_follow_up, bcr_status) ~ surgical_r, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
risk_score = rms::cph(Surv(days_to_follow_up, bcr_status) ~ risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
nom = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)

models <- c("age", "path_t", "surgical_r", "risk_score", "nom")

roc_curves <- lapply(models, function(model) {
  # Subset data to include only observations within the specified time
  subset_data <- mat2[mat2$days_to_follow_up <= 365*3, ]
  roc_obj <- roc(subset_data$bcr_status, predict(get(model), newdata = subset_data))
  return(roc_obj)
})

pdf("/data/github/pca_network/results/TCGA_DFS/nomogram_v_rest_3year.pdf", width=7, height=6)
plot(roc_curves[[1]], time=365, col="black", title = F, lwd=2)
plot(roc_curves[[2]], time=floor(365), col="yellow", add = T, title=F, lwd=2)
plot(roc_curves[[3]], time=floor(365), col="green", add = T, title=F, lwd=2)
plot(roc_curves[[4]], time=floor(365), col="blue", add = T, title=F, lwd=2)
plot(roc_curves[[5]], time=floor(365), col="red", add = T, title=F, lwd=2)
my_legend =  c(paste0("Age AUC: ", round(roc_curves[[1]]$auc,2)),
               paste0("Pathological T AUC: ", round(roc_curves[[2]]$auc,2)),
               paste0("Resection status AUC: ", round(roc_curves[[3]]$auc,2)),
               paste0("Prognostic index AUC: ", round(roc_curves[[4]]$auc,2)),
               paste0("Nomogram AUC: ", round(roc_curves[[5]]$auc,2)))
legend("bottomright", legend = my_legend,col=c("black","yellow", "green", "blue", "red"),lwd=2, cex=1)
dev.off()

# nom vs individual predictors 5 year
age = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
path_t = rms::cph(Surv(days_to_follow_up, bcr_status) ~  path_t, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
surgical_r = rms::cph(Surv(days_to_follow_up, bcr_status) ~ surgical_r, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
risk_score = rms::cph(Surv(days_to_follow_up, bcr_status) ~ risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
nom = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)

models <- c("age", "path_t", "surgical_r", "risk_score", "nom")

roc_curves <- lapply(models, function(model) {
  # Subset data to include only observations within the specified time
  subset_data <- mat2[mat2$days_to_follow_up <= 365*5, ]
  roc_obj <- roc(subset_data$bcr_status, predict(get(model), newdata = subset_data))
  return(roc_obj)
})

pdf("/data/github/pca_network/results/TCGA_DFS/nomogram_v_rest_5year.pdf", width=7, height=6)
plot(roc_curves[[1]], time=365, col="black", title = F, lwd=2)
plot(roc_curves[[2]], time=floor(365), col="yellow", add = T, title=F, lwd=2)
plot(roc_curves[[3]], time=floor(365), col="green", add = T, title=F, lwd=2)
plot(roc_curves[[4]], time=floor(365), col="blue", add = T, title=F, lwd=2)
plot(roc_curves[[5]], time=floor(365), col="red", add = T, title=F, lwd=2)
my_legend =  c(paste0("Age AUC: ", round(roc_curves[[1]]$auc,2)),
               paste0("Pathological T AUC: ", round(roc_curves[[2]]$auc,2)),
               paste0("Resection status AUC: ", round(roc_curves[[3]]$auc,2)),
               paste0("Prognostic index AUC: ", round(roc_curves[[4]]$auc,2)),
               paste0("Nomogram AUC: ", round(roc_curves[[5]]$auc,2)))
legend("bottomright", legend = my_legend,col=c("black","yellow", "green", "blue", "red"),lwd=2, cex=1)
dev.off()

# nom vs individual predictors 8 year
age = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
path_t = rms::cph(Surv(days_to_follow_up, bcr_status) ~  path_t, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
surgical_r = rms::cph(Surv(days_to_follow_up, bcr_status) ~ surgical_r, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
risk_score = rms::cph(Surv(days_to_follow_up, bcr_status) ~ risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)
nom = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)

models <- c("age", "path_t", "surgical_r", "risk_score", "nom")

roc_curves <- lapply(models, function(model) {
  # Subset data to include only observations within the specified time
  subset_data <- mat2[mat2$days_to_follow_up <= 365*8, ]
  roc_obj <- roc(subset_data$bcr_status, predict(get(model), newdata = subset_data))
  return(roc_obj)
})

pdf("/data/github/pca_network/results/TCGA_DFS/nomogram_v_rest_8year.pdf", width=7, height=6)
plot(roc_curves[[1]], time=365, col="black", title = F, lwd=2)
plot(roc_curves[[2]], time=floor(365), col="yellow", add = T, title=F, lwd=2)
plot(roc_curves[[3]], time=floor(365), col="green", add = T, title=F, lwd=2)
plot(roc_curves[[4]], time=floor(365), col="blue", add = T, title=F, lwd=2)
plot(roc_curves[[5]], time=floor(365), col="red", add = T, title=F, lwd=2)
my_legend =  c(paste0("Age AUC: ", round(roc_curves[[1]]$auc,2)),
               paste0("Pathological T AUC: ", round(roc_curves[[2]]$auc,2)),
               paste0("Resection status AUC: ", round(roc_curves[[3]]$auc,2)),
               paste0("Prognostic index AUC: ", round(roc_curves[[4]]$auc,2)),
               paste0("Nomogram AUC: ", round(roc_curves[[5]]$auc,2)))
legend("bottomright", legend = my_legend,col=c("black","yellow", "green", "blue", "red"),lwd=2, cex=1)
dev.off()


cindex <- validate(nom, mat2, event = mat2$bcr_status)




## DCA 
require(dcurves)
Nomogram <- coxph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2)
T_stage <- coxph(Surv(days_to_follow_up, bcr_status) ~ path_t, data=mat2)
Resection <- coxph(Surv(days_to_follow_up, bcr_status) ~ surgical_r, data=mat2)
Age <- coxph(Surv(days_to_follow_up, bcr_status) ~ age, data=mat2)
Risk <- coxph(Surv(days_to_follow_up, bcr_status) ~ risk_score, data=mat2)

# tbl_regression(Nomogram, exponentiate = TRUE)
# 
# mat2_u1 <- broom::augment( Nomogram, newdata = mat2 %>% mutate(days_to_follow_up = max(mat2$days_to_follow_up)), type.predict = "expected" ) %>% mutate( Nomogram = 1 - exp(-.fitted) ) 
# mat2_u2 <- broom::augment( T_stage, newdata = mat2 %>% mutate(days_to_follow_up = max(mat2$days_to_follow_up)), type.predict = "expected" ) %>% mutate( T_stage = 1 - exp(-.fitted) ) 
# mat2_u3 <- broom::augment( Resection, newdata = mat2 %>% mutate(days_to_follow_up = max(mat2$days_to_follow_up)), type.predict = "expected" ) %>% mutate( Resection = 1 - exp(-.fitted) ) 
# mat2_u4 <- broom::augment( Age, newdata = mat2 %>% mutate(days_to_follow_up = max(mat2$days_to_follow_up)), type.predict = "expected" ) %>% mutate( Age = 1 - exp(-.fitted) ) 
# mat2_u5 <- broom::augment( Risk, newdata = mat2 %>% mutate(days_to_follow_up = max(mat2$days_to_follow_up)), type.predict = "expected" ) %>% mutate( Risk = 1 - exp(-.fitted) ) 


mat2 <-
  mat2 %>%
  mutate(
    nomogram =
      1 - summary(survfit(Nomogram, newdata = mat2), times = 365*8)$surv[1, ],
    t_stage = 
      1 - summary(survfit(T_stage, newdata = mat2), times = 365*8)$surv[1, ],
    resection =
      1 - summary(survfit(Resection, newdata = mat2), times = 365*8)$surv[1, ],
    age =
      1 - summary(survfit(Age, newdata = mat2), times = 365*8)$surv[1, ],
    risk =
      1 - summary(survfit(Risk, newdata = mat2), times = 365*8)$surv[1, ],
  )

dca(Surv(days_to_follow_up, bcr_status) ~ nomogram + t_stage + resection + age + risk,
    data = mat2,
    time = 365*8,
    thresholds = seq(0, 1, 0.01)) %>%
  plot(smooth = FALSE)


df <- merge(x=mat2_u1,y=mat2_u2,by=".rownames", all.x = TRUE)
df <- merge(x=df,y=mat2_u3,by=".rownames", all.x = TRUE)
df <- merge(x=df,y=mat2_u4,by=".rownames", all.x = TRUE)
df <- merge(x=df,y=mat2_u5,by=".rownames", all.x = TRUE)

dca(Surv(days_to_follow_up, bcr_status) ~ Risk, 
    data = df,
    time = max(df$days_to_follow_up),
    thresholds = thresholds) %>%
  plot(smooth = TRUE, ylim=c(0.3, -0.3))


# set seed for random process
set.seed(112358)

# create a 10-fold cross validation set
rsample::vfold_cv(mat2, v = 4, repeats = 25) %>%
  # for each cut of the data, build logistic regression on the 90% (analysis set),
  # and perform DCA on the 10% (assessment set)
  rowwise() %>%
  mutate(
    # build regression model on analysis set
    cox_mod =
      coxph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data = rsample::analysis(splits)) %>%
      list(),
    # get predictions for assessment set
    df_assessment =
      broom::augment(
        cox_mod,
        newdata = rsample::assessment(splits),
        type.predict = "expected"
      ) %>%
     mutate( .fitted = 1 - exp(-.fitted) ) %>%
      list(),
    # calculate net benefit on assessment set
    dca_assessment =
      dca(Surv(days_to_follow_up, bcr_status) ~ .fitted,
          data = df_assessment,
          time = 365*5,
          thresholds = seq(0, 0.35, 0.01),
          label = list(.fitted = "Cross-validated Prediction Model")
      ) %>%
      as_tibble() %>%
      list()
  ) %>%
  # pool results from the 10-fold cross validation
  pull(dca_assessment) %>%
  bind_rows() %>%
  group_by(variable, label, threshold) %>%
  summarise(net_benefit = mean(net_benefit), .groups = "drop") %>%
  # plot cross validated net benefit values
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.014, 1)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  theme_bw()




# nomogram C-index
# reload mat2 etc to keep this clean
# load datdists..

surv.obj=with(mat2,Surv(days_to_follow_up,bcr_status))
cox.mod=rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + path_t + surgical_r + risk_score, data=mat2, surv = T, x=T, y=T, time.inc = 356*8)

###Create your survival estimates
estimates=survest(cox.mod,newdata=mat2,times=8*365)$surv


###Determine concordance
c_index = rcorr.cens(x=estimates,S=surv.obj)
c_index[1]









#########################################################################
# Perform univariate coxph, add signf to multivariate cox after. 
# Only 10 events, so need to bin patients e.g T2 vs T3 for path T
#########################################################################

univ = mat

# run code chunks below with roc to get cutoffs
mat$age_cat = ifelse(mat$age > age_thresh[1,1], "high", "low")
mat$gleason = ifelse(mat$gleason_score > glea_thresh[1,1], "Gleason high", "Gleason low")
mat$psa = ifelse(mat$preop_psa > psa_thresh[1,1], "PSA high", "PSA low")
mat$path_t = ifelse(mat$ajcc_pathologic_t == "T2a" | mat$ajcc_pathologic_t == "T2b" | mat$ajcc_pathologic_t == "T2c", "T2", "T3 + T4")
mat$path_n = ifelse(mat$ajcc_pathologic_n == "N0", "N0", "N1")
mat$clin_m = ifelse(mat$ajcc_clinical_m == "M0", "M0", "M1")
mat$risk_cat_num = ifelse(mat$risk_category=="high", 1,0)

mat$psa <- ifelse(mat$preop_psa > 20, ">20ng/ml",
                  ifelse(mat$preop_psa >= 10 & mat$preop_psa <= 20, "10-20ng/ml", "<10ng/ml"))

mat$path_t = ifelse(mat$ajcc_pathologic_t == "T2a" | mat$ajcc_pathologic_t == "T2b" | mat$ajcc_pathologic_t == "T2c", "T2", mat$ajcc_pathologic_t)
mat$path_t = ifelse(mat$ajcc_pathologic_t == "T3a" | mat$ajcc_pathologic_t == "T3b", "T3", mat$path_t)
mat$clin_m = ifelse(mat$ajcc_clinical_m == "M0", "M0", "M1")
mat$path_n = ifelse(mat$ajcc_pathologic_n == "N1", "N1","N0")
mat$gleason = ifelse(mat$gleason_score == 6, "Gleason 6 or lower", mat$gleason_score)
mat$gleason = ifelse(mat$gleason_score == 7, "Gleason 7", mat$gleason)
mat$gleason = ifelse(mat$gleason_score > 7, "Gleason 8,9,10", mat$gleason)
mat$surgical_r = ifelse(mat$Surgical.Margin.Resection.Status == "R2", NA, mat$Surgical.Margin.Resection.Status)

vars_for_table = c("age","psa", "gleason", "path_t", "path_n", "clin_m", "surgical_r", "risk_category")

univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(days_to_follow_up, bcr_status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mat)})

#pdf("/data/github/pca_network/results/TCGA_DFS3/univariate_forest_model.pdf", width=8, height=5)
forestmodel::forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)
#dev.off()
mult = mat[which(mat$Row.names %in% samples_for_multivariate),]

# TCGA-V1-A9OT-01A NA for path_n - which is a 1 bcr status for gleason 6
mat[which(mat$Row.names=="TCGA-V1-A9OT-01A"),]$path_n = "N0"

#pdf("/data/github/pca_network/results/TCGA_DFS3/multivariate_forest_model.pdf", height=5, width=8)
forestmodel::forest_model(coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + CTHRC1 + JAG2, data=univ_mat))
forestmodel::forest_model(coxph(Surv(days_to_follow_up, bcr_status) ~ age + psa + gleason + path_t + path_n + clin_m + surgical_r + risk_score, data=mat))
#dev.off()

mat$years_to_follow_up = floor(mat$days_to_follow_up/365.25)
#mat$path_t = ifelse(mat$ajcc_pathologic_t == "T2a" | mat$ajcc_pathologic_t == "T2b" | mat$ajcc_pathologic_t == "T2c", "T2", "T3 + T4")
mult_cox = coxph(Surv(days_to_follow_up, bcr_status) ~ path_t + surgical_r + risk_score, data=mat)
x = regplot::regplot(clin_mod, failtime = c(365,1095,1826), points = T)

d = rms::datadist(mat)
options(datadist="d")
rms_cox = rms::cph(Surv(days_to_follow_up, bcr_status) ~ path_t + risk_score, data=mat2, surv = T, x=T, y=T)
med = Quantile(rms_cox)
surv = Survival(rms_cox)

plot(nomogram(rms_cox, fun=function(x) med(lp=x), funlabel = "Median Survival Time"))
nom <- nomogram(rms_cox, fun=list(function(x) surv(365, x),
                                  function(x) surv(365*3, x),
                                  function(x) surv(365*5, x)),
                funlabel=c("1-year BCR Probability", 
                           "3-year BCR Probability",
                           "5-year BCR Probability"))
plot(nom)


# 1, 3, 5 yera? remake this qith a univariate model for each clinical var.. I think. 
# Load the pROC package
# Load the pROC package
# Load the pROC package
library(pROC)

risk_cox = rms::cph(Surv(days_to_follow_up, bcr_status) ~ risk_score, data=mat, surv = T, x=T, y=T)
patht_cox = rms::cph(Surv(days_to_follow_up, bcr_status) ~ path_t, data=mat, surv = T, x=T, y=T)
warnings()# Define the time points you're interested in (1, 3, and 5 years)
time_points <- c(365, 3*365, 5*365)  # Assuming days as the time unit, adjust as needed

# Create ROC curves for each time point
roc_curves <- lapply(time_points, function(time) {
  # Subset your data to include only observations within the specified time
  subset_data <- mat[mat$days_to_follow_up <= time, ]
  roc_obj <- roc(subset_data$bcr_status, predict(patht_cox, newdata = subset_data))
  return(roc_obj)
})


plot(roc_curves[[1]], time=365, col="red", title = F, lwd=2)
plot(roc_curves[[2]], time=floor(365*3), col="blue", add = T, title=F, lwd=2)
plot(roc_curves[[3]], time=floor(365*5), col="forestgreen", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(roc_curves[[1]]$auc,2)),
               paste0("3 year AUC: ", round(roc_curves[[2]]$auc,2)),
               paste0("5 year AUC:", round(roc_curves[[3]]$auc,2)))
legend("bottomright", legend = my_legend,col=c("red","blue", "forestgreen"),lwd=2, cex=1)



## GOO DRMS code below for generating calibration curves etc....
```{R}
mat = master
mat$psa <- ifelse(mat$preop_psa > 20, ">20ng/ml",
                  ifelse(mat$preop_psa >= 10 & mat$preop_psa <= 20, "10-20ng/ml", "<10ng/ml"))

mat$path_t = ifelse(mat$ajcc_pathologic_t == "T2a" | mat$ajcc_pathologic_t == "T2b" | mat$ajcc_pathologic_t == "T2c", "T2", mat$ajcc_pathologic_t)
mat$path_t = ifelse(mat$ajcc_pathologic_t == "T3a" | mat$ajcc_pathologic_t == "T3b", "T3", mat$path_t)
mat$clin_m = ifelse(mat$ajcc_clinical_m == "M0", "M0", "M1")
mat$path_n = ifelse(mat$ajcc_pathologic_n == "N1", "N1","N0")
mat$gleason = ifelse(mat$gleason_score == 6, "Gleason 6 or lower", mat$gleason_score)
mat$gleason = ifelse(mat$gleason_score == 7, "Gleason 7", mat$gleason)
mat$gleason = ifelse(mat$gleason_score > 7, "Gleason 8,9,10", mat$gleason)
mat$surgical_r = ifelse(mat$Surgical.Margin.Resection.Status == "R2" | mat$Surgical.Margin.Resection.Status == "R1",  "Residual tumor", mat$Surgical.Margin.Resection.Status)


w <- transcan(~ age + gleason + psa + path_t + path_n + clin_m +
                surgical_r,
              imputed=TRUE, data=mat, pl=FALSE, pr=FALSE)
```


```{r}

mat2 = mat

mat2$psa = impute(w, psa, data=mat)
mat2$psa = ifelse(mat2$psa == 1, "<10ng/ml",
                  ifelse(mat2$psa == 3, "10-20ng/ml", mat2$psa))


mat2$path_t = impute(w, path_t, data = mat)

mat2$path_t = ifelse(mat2$path_t == 1, "T2",
                     ifelse(mat2$path_t == 2, "T3", 
                            ifelse(mat2$path_t == 3, "T4", mat2$path_t)))
mat2$path_n = impute(w, path_n, data=mat)
mat2$path_n = ifelse(mat2$path_n == 1, "N0",
                     ifelse(mat2$path_n == 2, "N1", mat2$path_n))
mat2$clin_m = impute(w, clin_m, data=mat)
mat2$clin_m = ifelse(mat2$clin_m == 1, "M0", mat2$clin_m)
mat2$surgical_r = impute(w, surgical_r, data=mat)
mat2$surgical_r = ifelse(mat2$surgical_r==1, "R0",
                         ifelse(mat2$surgical_r == 2, "Residual tumor",
                                ifelse(mat2$surgical_r == 3 , "RX", mat2$surgical_r)))

mat2 = mat2[,c("days_to_follow_up", "bcr_status", "age", "gleason", "path_t", "path_n", "clin_m", "surgical_r")]
gene_mat = logcpm[keep,]
gene_mat =as.data.frame(scale(t(gene_mat), scale = T, center = T))
mat2 = cbind(mat2, gene_mat)


dd = datadist(mat2)
options(datadist = "dd")
rms_cox = rms::cph(Surv(days_to_follow_up, bcr_status) ~ age + gleason + surgical_r + path_n + clin_m + path_t + REG4 + CTHRC1 + INPP5E + DDAH1, data=mat2, x=TRUE, y=TRUE, surv=TRUE, time.inc = 365*5)
rms_cox
```

```{R}
heart <- hx + ekg %nin% c('normal','benign')
label(heart) <- 'Heart Disease Code'
map   <- (2*dbp + sbp)/3
label(map) <- 'Mean Arterial Pressure/10'
dd <- datadist(dd, heart, map)

f <- cph(S ~ rx + rcs(age,4) + rcs(wt,3) + pf.coded +
           heart + rcs(map,3) + rcs(hg,4) +
           rcs(sg,3) + rcs(sz,3) + rcs(log(ap),5) + bm,
         x=TRUE, y=TRUE, surv=TRUE, time.inc=5*12)
print(f, coefs=FALSE)
```

```{R}
anova(rms_cox, test = "LR")

```


```{R, figure.width=12, figure.height=18}
z <- predict(rms_cox, type='terms')
# required x=T above to store design matrix
S = Surv(days_to_follow_up, bcr_status)
f.short <- cph(S ~ z, x=TRUE, y=TRUE)
# store raw x, y so can get residuals
require(survival)   # or use survival::cox.zph(...)
phtest <- cox.zph(rms_cox, transform='identity')
phtest

ggcoxzph2(phtest[10])
```

```{R}
set.seed(1)  # so can reproduce results
v <- validate(rms_cox, B=300)
v

cal1yr <- calibrate(rms_cox, B=300, u=365*1, maxdim=3, conf.int=TRUE, cmethod = "KM")
cal3yr <- calibrate(rms_cox, B=300, u=365*5, maxdim=3, conf.int=TRUE, cmethod = "KM")
plot(cal1yr, xlim=c(0,1), ylim=c(0,1))
plot(cal3yr,  xlim=c(0,1), ylim=c(0,1))
```

```{R}

plot(summary(rms_cox))
#spar(ps=8)
surv  <- Survival(rms_cox)
surv1 <- function(x) surv(365,x)
surv3 <- function(x) surv(365*3,x)
surv5 <- function(x) surv(365*5,x)
quan  <- Quantile(rms_cox)
med   <- function(x) quan(lp=x)/365
ss    <- c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95)

nom <- nomogram(rms_cox,
                fun=list(surv1, surv3, surv5),
                funlabel=c('1-year Survival', '3-year Survival','5-year Survival'),
                fun.at=list(ss, ss, c(.5,1:6)))
plot(nom, xfrac=.65, lmgp=.35)
```









#######################
# nomogram 1, 3, 5 year
plot(roc(mat$bcr_status, predict(rms_cox)), lwd=2)
# Predict a patient: yhat gives the survival probability as if you were drawing a line tourself.
plot(Predict(rms_cox, type="predictions", time=365*5))

cal1 <- calibrate(rms_cox,X=T,Y=T, method='boot',B=228, u=c(365), cmethod = "KM")
cal2 <- calibrate(rms_cox, X=T, Y=T, method="boot", B=228, u=1024, cmethod="KM")
plot(cal1, xlim=c(0,1), ylim=c(0,1))
abline(coef = c(0,1), lty=2)


# Pathologic T: T1+T2
#pdf("/data/github/pca_network/results/TCGA_DFS3/univariate_clinical.pdf")
t2 = mat[which(mat$path_t==1),]

t2_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=t2)

ggsurvplot(t2_cox,
           
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=t2,
           xlab="Time (days) path T2")
t3_t4 = mat[which(mat$path_t!=1),]
t3_t4_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=t3_t4)
ggsurvplot(t3_t4_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=t3_t4,
           xlab="Time (days) path T3+T4")

N0 = mat[which(mat$path_n==1),]
N0_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=N0)
ggsurvplot(N0_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme

           palette = c("red1", "royalblue3"),
           data=N0,
           xlab="Time (days) path N0")

N1 = mat[which(mat$path_n!=1),]
N1_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=N1)
ggsurvplot(N1_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=N1,
           xlab="Time (days) path N1")

M0 = mat[which(mat$clin_m==1),]
M0_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=M0)
ggsurvplot(M0_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=M0,
           xlab="Time (days) clin M0")

M1 = mat[which(mat$clin_m!=1),]
M1_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=M1)
ggsurvplot(M1_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=M1,
           xlab="Time (days) clin M1")

age_thresh = roc(mat$risk_category, mat$age)
age_thresh = coords(age_thresh, "best")
low_age = mat[which(mat$age <= age_thresh[1,1]),]
low_age_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=low_age)
ggsurvplot(low_age_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=low_age,
           xlab="Time (days) low age")

high_age = mat[which(mat$age >= age_thresh[1,1]),]
high_age_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=high_age)
ggsurvplot(high_age_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=high_age,
           xlab="Time (days) high age 62.5")

glea_thresh = roc(mat$status, mat$gleason_score)
glea_thresh = coords(glea_thresh, "best")
low_gle = mat[which(mat$gleason_score <= glea_thresh[1,1]),]
low_gle_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=low_gle)
ggsurvplot(low_gle_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=low_gle,
           xlab="Time (days) low gleason")

high_gle = mat[which(mat$gleason_score >= glea_thresh[1,1]),]
high_gle_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=high_gle)
ggsurvplot(high_gle_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=high_gle,
           xlab="Time (days) high gleason 7.5")

psa_thresh = roc(mat$risk_category, mat$preop_psa)
psa_thresh = coords(psa_thresh, "best")
low_psa = mat[which(mat$preop_psa <= psa_thresh[1,1]),]
low_psa_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=low_psa)
ggsurvplot(low_psa_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=low_psa,
           xlab="Time (days) low psa 6.85")

high_psa = mat[which(mat$preop_psa >= psa_thresh[1,1]),]
high_psa_cox = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=high_psa)
ggsurvplot(high_psa_cox,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=high_psa,
           xlab="Time (days) high psa")
#dev.off()