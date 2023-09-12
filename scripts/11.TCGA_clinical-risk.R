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
library(gtsummary) # format cox model 
library(RColorBrewer)

#########################################################################
# Use FPKM for plotting expression data (scatterhist)
#########################################################################
atlas_mat = read.csv("/data/github/pca_network/data/prad_rna_cancer_sample.tsv", header=T, sep="\t")
atlas_mat = atlas_mat[,c(1,2,4)]
atlas_mat <- atlas_mat %>%
  pivot_wider(names_from = Sample, values_from = FPKM, values_fill = 0)
atlas_mat = tibble::column_to_rownames(atlas_mat, "Gene")

#########################################################################
# Load atlas metadata
# Remove missing OS status
#########################################################################
atlas_meta = read.csv("/data/github/pca_network/data/tcga_updated_meta.csv", header=T, sep=",")
rownames(atlas_meta) = atlas_meta$sample
atlas_meta = atlas_meta[which(atlas_meta$sample %in% colnames(atlas_mat)),]
rem = !(atlas_meta$vital=="Not Reported")
atlas_meta = atlas_meta[rem,]
atlas_mat = atlas_mat[,rem]
atlas_meta$status = ifelse(atlas_meta$vital=="Dead", 1, 0)

#########################################################################
# Load cleaned tcga meta
#########################################################################
load("/data/github/pca_network/data/TCGA_meta_cleaned.RData")
resec = read.csv("~/Downloads/prad_tcga_clinical_data.tsv", header=T, sep="\t")
resec$Sample.ID = paste(resec$Sample.ID, "A", sep="")
table(resec$Sample.ID %in% master$Row.names)
master = merge(master, subset(resec, select=c(Sample.ID, Surgical.Margin.Resection.Status)), by.x="Row.names", by.y="Sample.ID")
#master$years_to_follow_up <- as.character(floor(atlas_meta$years_to_follow_up))
load("~/Desktop/SEC31B-INPP5E-CACNB3-CTHRC1-KIFC2_462.RData")
#########################################################################
# Convert ENS to symbols and subset matrix by Active.Genes,
# order correctly (FPKM)
# FPKM matrix now ready for expression plots
#########################################################################
#load("/data/github/pca_network/results/prognostic_model_os2.RData")
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
# to generate risk scores
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

mat = as.data.frame(scaled[,Active.Genes]) ## change back to = scaled when done
risk = rowSums(Active.Coefficients*mat)
mat$risk_score = risk
mat$risk_category = ifelse(mat$risk_score > median(mat$risk_score), "high", "low")
mat = merge(mat, master, by.x=0, by.y="Row.names")

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

vars_for_table = c("age","psa", "gleason", "path_t", "path_n", "clin_m", "surgical_r", "risk_score")

univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(days_to_follow_up, bcr_status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mat)})

#pdf("/data/github/pca_network/results/TCGA_DFS3/univariate_forest_model.pdf", width=8, height=5)
forestmodel::forest_model(model_list = univ_models,covariates = vars_for_table,merge_models =T)
#dev.off()
mult = mat[which(mat$Row.names %in% samples_for_multivariate),]

# TCGA-V1-A9OT-01A NA for path_n - which is a 1 bcr status for gleason 6
mat[which(mat$Row.names=="TCGA-V1-A9OT-01A"),]$path_n = "N0"

#pdf("/data/github/pca_network/results/TCGA_DFS3/multivariate_forest_model.pdf", height=5, width=8)
forestmodel::forest_model(coxph(Surv(days_to_follow_up, bcr_status) ~ age + psa + gleason + path_t + path_n + clin_m + surgical_r + risk_score, data=mat))
#dev.off()

mat$years_to_follow_up = floor(mat$days_to_follow_up/365.25)
#mat$path_t = ifelse(mat$ajcc_pathologic_t == "T2a" | mat$ajcc_pathologic_t == "T2b" | mat$ajcc_pathologic_t == "T2c", "T2", "T3 + T4")
mult_cox = coxph(Surv(days_to_follow_up, bcr_status) ~ path_t + surgical_r + risk_score, data=mat)
x = regplot::regplot(mult_cox, failtime = c(365,1095,1826), points = T)

d = rms::datadist(mat)
options(datadist="d")
rms_cox = rms::cph(Surv(days_to_follow_up, bcr_status) ~ path_t + risk_score, data=mat, surv = T, x=T, y=T)
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