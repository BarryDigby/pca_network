#!/usr/bin/env Rscript 

library(survival)
library(sjPlot)
library(survminer)
library(timeROC)
library(ggpubr)
library(rms)
library(survAUC)
library(pec)
load("/data/github/pca_network/results/TCGA_DFS/model.RData")

counting <- function(df){
  cases <- nrow(df)
  bcr_overall <- table(df$bcr_status)
  bcr_events <- bcr_overall[[2]]
  year_one = df[which(df$days_to_follow_up < 365),]
  year_one =  nrow(year_one[which(year_one$bcr_status==1),])
  year_three = df[which(df$days_to_follow_up >= 365 & df$days_to_follow_up < 1095 ),]
  year_three =  nrow(year_three[which(year_three$bcr_status==1),])
  year_five = df[which(df$days_to_follow_up >= 1095 & df$days_to_follow_up < 1825 ),]
  year_five =  nrow(year_five[which(year_five$bcr_status==1),])
  year_eight = df[which(df$days_to_follow_up >= 1825 & df$days_to_follow_up < 2920 ),]
  year_eight=  nrow(year_eight[which(year_eight$bcr_status==1),])
  
  print(paste("Number of cases: ", cases))
  print(paste("Overall BCR: \n", bcr_overall))
  print(paste("BCR events: ", bcr_events))
  print(paste("Year one events: ", year_one))
  print(paste("Year three events: ", year_three))
  print(paste("Year five events: ", year_five))
  print(paste("Year eight events: ", year_eight))
}

# Belfast dataset
belfast_rds = readRDS("/data/github/pca_network/data/Belfast_eSet.RDS")
mat = belfast_rds@assayData$exprs
meta = belfast_rds@phenoData@data
meta$days_to_follow_up = meta$time_to_bcr * 30.44
rm(belfast_rds)

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat = tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

mat = as.data.frame(scale(mat, scale = T, center = T))

belfast_risk = predict(final_model, newdata = mat, type="lp")

sub_meta = subset(meta, select=c(days_to_follow_up, bcr_status))
mat = cbind(mat, sub_meta)

# Method 1: Regression on PI in validation data
# create new variables to match worked example:
fit_dev = final_model
df_dev = cox
Xb_dev <- model.matrix(fit_dev) %*% coef(fit_dev)
Xbavg_dev <- sum(coef(fit_dev)*fit_dev$means) 
PI_dev <- Xb_dev - Xbavg_dev # centered PI in development dataset (for discrimination later)
df_dev$PI_dev <- PI_dev

df_val <- mat
# fit model in validation data
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~  REG4 + SLC2A4 + CTHRC1 + JAG2, data=df_val) 
Xb_val <- model.matrix(fit_val) %*% coef(fit_dev) # determine Xb in validation data 
PI_val <- Xb_val - Xbavg_dev # center PI by using mean of PI of development data 
df_val$PI_val <- PI_val

fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val, data=df_val)
tab_model(fit_val, transform=NULL, show.r2 = FALSE)

fit_val_test <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val + offset(PI_val), data=df_val) 
tab_model(fit_val_test, transform=NULL, show.r2 = FALSE)

# Method 2: Check model misspecification/fit
# offset = set PI to 1 i,e ike fit_dev 
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + CTHRC1 + JAG2 + offset(PI_val), data=df_val)
round(2*(diff(fit_val$loglik)), 2) # Chi-squared value
round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3) # p-value


# Method 3: Measures of discrimination
# Derivation
rcorr.cens(-1*PI_dev, Surv(df_dev$days_to_follow_up, df_dev$bcr_status))[1] # Harrell's c-index
GHCI(PI_dev) # Gonen and Heller
# Validation
rcorr.cens(-1*PI_val, Surv(df_val$days_to_follow_up, df_val$bcr_status))[1] # Harrell's c-index
GHCI(df_val$PI_val) # Gonen and Heller

# Method 4: Kaplan-Meier curves for risk groups
# Derivation
df_dev$Risk_group_dev = ifelse(df_dev$PI_dev > median(df_dev$PI_dev), "High", "Low")
surv_object <- Surv(df_dev$days_to_follow_up, df_dev$bcr_status)
res2 = survfit(surv_object ~ Risk_group_dev, data=df_dev)
p1 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_dev,
                 xlab="Time (days)")
p1
# Validation
df_val$Risk_group_val = ifelse(df_val$PI_val > median(df_val$PI_val), "High", "Low")
surv_object <- Surv(df_val$days_to_follow_up, df_val$bcr_status)
res2 = survfit(surv_object ~ Risk_group_val, data=df_val)
p2 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_val,
                 xlab="Time (days)")
p2
library(dplyr)
df_val_KM <- df_val %>%
  group_split(Risk_group_val)
surv_val <- lapply(df_val_KM, function(x){
  obj_surv <- survfit(Surv(days_to_follow_up, bcr_status) ~ 1, data = x)
  data.frame(time=obj_surv$time, surv=obj_surv$surv)
})

p1$plot +
  geom_step(data=surv_val[[1]], aes(x = time, y = surv), linetype=2) +
  geom_step(data=surv_val[[2]], aes(x = time, y = surv), linetype=2)


# MEthod 5 -logrank - skipped - its in the plot above. 

# Method 6. Hazard Ratios between risk groups
# Derivation
df_dev$Risk_group_dev <- factor(df_dev$Risk_group_dev)
df_dev$Risk_group_dev <- relevel(df_dev$Risk_group_dev, ref = "Low")
fit_dev_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_dev), data=df_dev)
tab_model(fit_dev_hr, show.r2 = FALSE)
# validation
df_val$Risk_group_val <- factor(df_val$Risk_group_val)
df_val$Risk_group_val <- relevel(df_val$Risk_group_val, ref = "Low")
fit_val_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_val), data=df_val)
tab_model(fit_val_hr, show.r2 = FALSE)


########################################
# plots for manuscript
########################################

os_mat = df_val
os_mat = os_mat[order(os_mat$PI_val, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$status = factor(os_mat$bcr_status)

pdf("/data/github/pca_network/results/external/Belfast/risk_score_dist.pdf", width = 8, height = 3)
ggscatter(os_mat, y="PI_val", x="patients_inc_risk", color="Risk_group_val", fill="Risk_group_val", 
          ylab="Prognostic index", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = median(os_mat$PI_val), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/external/Belfast/scatter_dfs.pdf", width = 8, height = 3)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(status), x="patients_inc_risk", y="days_to_follow_up", shape="status", ylab = "Time (days)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$Risk_group_val)

col <- c("red", "royalblue3")
names(col) <- c("High", "Low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/external/Belfast/genes_heatmap.pdf", height = 3, width = 8)
pheatmap::pheatmap(t(os_mat[,c("REG4", "SLC2A4", "CTHRC1", "JAG2")]), labels_col = FALSE, color = col_palette, cluster_cols = F, cluster_rows = F,
                   scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
dev.off()

pdf("/data/github/pca_network/results/external/Belfast/high_vs_low.pdf", width=7, height=6)
p2
dev.off()


########################################
# plots for manuscript
########################################


ROC.1 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)
ROC.2  <- timeROC(T=df_val$days_to_follow_up,
                  delta=df_val$bcr_status,
                  marker=df_val$PI_val,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

ROC.4 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*8)),
                 iid=TRUE)
pdf("/data/github/pca_network/results/external/Belfast/ROC_1358.pdf", width=7, height=6)
plot(ROC.1, time=365, title = F, lwd=2)
plot(ROC.2, time=floor(365*3), col="blue", add = T, title=F, lwd=2)
plot(ROC.3, time=floor(365*5), col="forestgreen", add = T, title=F, lwd=2)
plot(ROC.4, time=floor(365*8), col="black", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(ROC.1$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.1)$CI_AUC), collapse = "-")),
               paste0("3 year AUC: ", round(ROC.2$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.2)$CI_AUC), collapse = "-")),
               paste0("5 year AUC:", round(ROC.3$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.3)$CI_AUC), collapse = "-")),
               paste0("8 year AUC:", round(ROC.4$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.4)$CI_AUC), collapse = "-")))
legend("bottomright", legend = my_legend,col=c("red","blue", "forestgreen", "black"),lwd=2, cex=1)
dev.off()


# CPC dataset

cpc = readRDS("/data/github/pca_network/data/CPC-Gene_eSet.RDS")
mat = cpc@assayData$exprs
meta = cpc@phenoData@data
meta = meta[which(!is.na(meta$bcr_status)),]
mat = as.data.frame(mat[,meta$sample_id])
meta$days_to_follow_up = floor(meta$time_to_bcr*30.44)


#load("/data/github/pca_network/results/prognostic_model_os2.RData")
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat = tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

mat = as.data.frame(scale(mat, scale = T, center = T))

cpc_risk = predict(final_model, newdata = mat, type="lp")

sub_meta = subset(meta, select=c(days_to_follow_up, bcr_status))
mat = cbind(mat, sub_meta)

# Method 1: Regression on PI in validation data
# create new variables to match worked example:
fit_dev = final_model
df_dev = cox
Xb_dev <- model.matrix(fit_dev) %*% coef(fit_dev)
Xbavg_dev <- sum(coef(fit_dev)*fit_dev$means) 
PI_dev <- Xb_dev - Xbavg_dev # centered PI in development dataset (for discrimination later)
df_dev$PI_dev <- PI_dev

df_val <- mat
# fit model in validation data
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~  REG4 + SLC2A4 + CTHRC1 + JAG2, data=df_val) 
Xb_val <- model.matrix(fit_val) %*% coef(fit_dev) # determine Xb in validation data 
PI_val <- Xb_val - Xbavg_dev # center PI by using mean of PI of development data 
df_val$PI_val <- PI_val

fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val, data=df_val)
tab_model(fit_val, transform=NULL, show.r2 = FALSE)

fit_val_test <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val + offset(PI_val), data=df_val) 
tab_model(fit_val_test, transform=NULL, show.r2 = FALSE)

# Method 2: Check model misspecification/fit
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + CTHRC1 + JAG2 + offset(PI_val), data=df_val)
round(2*(diff(fit_val$loglik)), 2) # Chi-squared value
round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3) # p-value


# Method 3: Measures of discrimination
# Derivation
rcorr.cens(-1*PI_dev, Surv(df_dev$days_to_follow_up, df_dev$bcr_status))[1] # Harrell's c-index
GHCI(PI_dev) # Gonen and Heller
# Validation
rcorr.cens(-1*PI_val, Surv(df_val$days_to_follow_up, df_val$bcr_status))[1] # Harrell's c-index
GHCI(df_val$PI_val) # Gonen and Heller

# Method 4: Kaplan-Meier curves for risk groups
# Derivation
df_dev$Risk_group_dev = ifelse(df_dev$PI_dev > median(df_dev$PI_dev), "High", "Low")
surv_object <- Surv(df_dev$days_to_follow_up, df_dev$bcr_status)
res2 = survfit(surv_object ~ Risk_group_dev, data=df_dev)
p1 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_dev,
                 xlab="Time (days)")
p1
# Validation
df_val$Risk_group_val = ifelse(df_val$PI_val > median(df_val$PI_val), "High", "Low")
surv_object <- Surv(df_val$days_to_follow_up, df_val$bcr_status)
res2 = survfit(surv_object ~ Risk_group_val, data=df_val)
p2 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_val,
                 xlab="Time (days)")
p2
library(dplyr)
df_val_KM <- df_val %>%
  group_split(Risk_group_val)
surv_val <- lapply(df_val_KM, function(x){
  obj_surv <- survfit(Surv(days_to_follow_up, bcr_status) ~ 1, data = x)
  data.frame(time=obj_surv$time, surv=obj_surv$surv)
})

p1$plot +
  geom_step(data=surv_val[[1]], aes(x = time, y = surv), linetype=2) +
  geom_step(data=surv_val[[2]], aes(x = time, y = surv), linetype=2)


# MEthod 5 -logrank - skipped

# Method 6. Hazard Ratios between risk groups
# Derivation
# Method 6. Hazard Ratios between risk groups
# Derivation
df_dev$Risk_group_dev <- factor(df_dev$Risk_group_dev)
df_dev$Risk_group_dev <- relevel(df_dev$Risk_group_dev, ref = "Low")
fit_dev_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_dev), data=df_dev)
tab_model(fit_dev_hr, show.r2 = FALSE)
# validation
df_val$Risk_group_val <- factor(df_val$Risk_group_val)
df_val$Risk_group_val <- relevel(df_val$Risk_group_val, ref = "Low")
fit_val_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_val), data=df_val)
tab_model(fit_val_hr, show.r2 = FALSE)



########################################
# plots for manuscript
########################################

os_mat = df_val
os_mat = os_mat[order(os_mat$PI_val, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$status = factor(os_mat$bcr_status)

pdf("/data/github/pca_network/results/external/CPC/risk_score_dist.pdf", width = 8, height = 3)
ggscatter(os_mat, y="PI_val", x="patients_inc_risk", color="Risk_group_val", fill="Risk_group_val", 
          ylab="Prognostic index", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = median(os_mat$PI_val), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/external/CPC/scatter_dfs.pdf", width = 8, height = 3)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(status), x="patients_inc_risk", y="days_to_follow_up", shape="status", ylab = "Time (days)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$Risk_group_val)

col <- c("red", "royalblue3")
names(col) <- c("High", "Low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/external/CPC/genes_heatmap.pdf", height = 3, width = 8)
pheatmap::pheatmap(t(os_mat[,c("REG4", "SLC2A4", "CTHRC1", "JAG2")]), labels_col = FALSE, color = col_palette, cluster_cols = F, cluster_rows = F,
                   scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
dev.off()

pdf("/data/github/pca_network/results/external/CPC/high_vs_low.pdf", width=7, height=6)
p2
dev.off()


########################################
# plots for manuscript
########################################


ROC.1 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)
ROC.2  <- timeROC(T=df_val$days_to_follow_up,
                  delta=df_val$bcr_status,
                  marker=df_val$PI_val,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

ROC.4 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*8)),
                 iid=TRUE)
pdf("/data/github/pca_network/results/external/CPC/ROC_1358.pdf", width=7, height=6)
plot(ROC.1, time=365, title = F, lwd=2)
plot(ROC.2, time=floor(365*3), col="blue", add = T, title=F, lwd=2)
plot(ROC.3, time=floor(365*5), col="forestgreen", add = T, title=F, lwd=2)
plot(ROC.4, time=floor(365*8), col="black", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(ROC.1$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.1)$CI_AUC), collapse = "-")),
               paste0("3 year AUC: ", round(ROC.2$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.2)$CI_AUC), collapse = "-")),
               paste0("5 year AUC:", round(ROC.3$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.3)$CI_AUC), collapse = "-")),
               paste0("8 year AUC:", round(ROC.4$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.4)$CI_AUC), collapse = "-")))
legend("bottomright", legend = my_legend,col=c("red","blue", "forestgreen", "black"),lwd=2, cex=1)
dev.off()


# GSE54460

GSE54460 = readRDS("/data/github/pca_network/data/GSE54460_eSet.RDS")
mat = GSE54460@assayData$exprs
meta = GSE54460@phenoData@data
meta = meta[which(!is.na(meta$bcr_status)),]
mat = as.data.frame(mat[,rownames(meta)])
meta$days_to_follow_up = floor(meta$time_to_bcr*30.44)

#load("/data/github/pca_network/results/prognostic_model_os.RData")

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat = tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

mat = as.data.frame(scale(mat, scale = T, center = T))

gse_risk = predict(final_model, newdata = mat, type="lp")

sub_meta = subset(meta, select=c(days_to_follow_up, bcr_status))
mat = cbind(mat, sub_meta)

# Method 1: Regression on PI in validation data
# create new variables to match worked example:
fit_dev = final_model
df_dev = cox
Xb_dev <- model.matrix(fit_dev) %*% coef(fit_dev)
Xbavg_dev <- sum(coef(fit_dev)*fit_dev$means) 
PI_dev <- Xb_dev - Xbavg_dev # centered PI in development dataset (for discrimination later)
df_dev$PI_dev <- PI_dev

df_val <- mat
# fit model in validation data
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~  REG4 + SLC2A4 + CTHRC1 + JAG2, data=df_val) 
Xb_val <- model.matrix(fit_val) %*% coef(fit_dev) # determine Xb in validation data 
PI_val <- Xb_val - Xbavg_dev # center PI by using mean of PI of development data 
df_val$PI_val <- PI_val

fit_dev = coxph(Surv(days_to_follow_up, bcr_status) ~ PI_dev, data=df_dev)

fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val, data=df_val)
tab_model(fit_val, transform=NULL, show.r2 = FALSE)

fit_val_test <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val + offset(PI_val), data=df_val) 
tab_model(fit_val_test, transform=NULL, show.r2 = FALSE)

# Method 2: Check model misspecification/fit
fit_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + CTHRC1 + JAG2 + offset(PI_dev), data=df_dev)
round(2*(diff(fit_dev$loglik)), 2) # Chi-squared value
round(1-pchisq(2*(diff(fit_dev$loglik)), 8), 3) # p-value

fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + CTHRC1 + JAG2 + offset(PI_val), data=df_val)
round(2*(diff(fit_val$loglik)), 2) # Chi-squared value
round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3) # p-value


# Method 3: Measures of discrimination
# Derivation
rcorr.cens(-1*PI_dev, Surv(df_dev$days_to_follow_up, df_dev$bcr_status))[1] # Harrell's c-index
GHCI(PI_dev) # Gonen and Heller
# Validation
rcorr.cens(-1*PI_val, Surv(df_val$days_to_follow_up, df_val$bcr_status))[1] # Harrell's c-index
GHCI(df_val$PI_val) # Gonen and Heller

# Method 4: Kaplan-Meier curves for risk groups
# Derivation
df_dev$Risk_group_dev = ifelse(df_dev$PI_dev > median(df_dev$PI_dev), "High", "Low")
surv_object <- Surv(df_dev$days_to_follow_up, df_dev$bcr_status)
res2 = survfit(surv_object ~ Risk_group_dev, data=df_dev)
p1 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_dev,
                 xlab="Time (days)")
p1
# Validation
df_val$Risk_group_val = ifelse(df_val$PI_val > median(df_val$PI_val), "High", "Low")
surv_object <- Surv(df_val$days_to_follow_up, df_val$bcr_status)
res2 = survfit(surv_object ~ Risk_group_val, data=df_val)
p2 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_val,
                 xlab="Time (days)")
p2
library(dplyr)
df_val_KM <- df_val %>%
  group_split(Risk_group_val)
surv_val <- lapply(df_val_KM, function(x){
  obj_surv <- survfit(Surv(days_to_follow_up, bcr_status) ~ 1, data = x)
  data.frame(time=obj_surv$time, surv=obj_surv$surv)
})

p1$plot +
  geom_step(data=surv_val[[1]], aes(x = time, y = surv), linetype=2) +
  geom_step(data=surv_val[[2]], aes(x = time, y = surv), linetype=2)


# MEthod 5 -logrank - skipped

# Method 6. Hazard Ratios between risk groups
# Derivation
df_dev$Risk_group_dev <- factor(df_dev$Risk_group_dev)
df_dev$Risk_group_dev <- relevel(df_dev$Risk_group_dev, ref = "Low")
fit_dev_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_dev), data=df_dev)
tab_model(fit_dev_hr, show.r2 = FALSE)
# validation
df_val$Risk_group_val <- factor(df_val$Risk_group_val)
df_val$Risk_group_val <- relevel(df_val$Risk_group_val, ref = "Low")
fit_val_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_val), data=df_val)
tab_model(fit_val_hr, show.r2 = FALSE)

########################################
# plots for manuscript
########################################

os_mat = df_val
os_mat = os_mat[order(os_mat$PI_val, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$status = factor(os_mat$bcr_status)

pdf("/data/github/pca_network/results/external/GSE54460/risk_score_dist.pdf", width = 8, height = 3)
ggscatter(os_mat, y="PI_val", x="patients_inc_risk", color="Risk_group_val", fill="Risk_group_val", 
          ylab="Prognostic index", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = median(os_mat$PI_val), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/external/GSE54460/scatter_dfs.pdf", width = 8, height = 3)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(status), x="patients_inc_risk", y="days_to_follow_up", shape="status", ylab = "Time (days)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$Risk_group_val)

col <- c("red", "royalblue3")
names(col) <- c("High", "Low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/external/GSE54460/genes_heatmap.pdf", height = 3, width = 8)
pheatmap::pheatmap(t(os_mat[,c("REG4", "SLC2A4", "CTHRC1", "JAG2")]), labels_col = FALSE, color = col_palette, cluster_cols = F, cluster_rows = F,
                   scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
dev.off()

pdf("/data/github/pca_network/results/external/GSE54460/high_vs_low.pdf", width=7, height=6)
p2
dev.off()


########################################
# plots for manuscript
########################################


ROC.1 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)
ROC.2  <- timeROC(T=df_val$days_to_follow_up,
                  delta=df_val$bcr_status,
                  marker=df_val$PI_val,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

ROC.4 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*8)),
                 iid=TRUE)
pdf("/data/github/pca_network/results/external/GSE54460/ROC_1358.pdf", width=7, height=6)
plot(ROC.1, time=365, title = F, lwd=2)
plot(ROC.2, time=floor(365*3), col="blue", add = T, title=F, lwd=2)
plot(ROC.3, time=floor(365*5), col="forestgreen", add = T, title=F, lwd=2)
plot(ROC.4, time=floor(365*8), col="black", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(ROC.1$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.1)$CI_AUC), collapse = "-")),
               paste0("3 year AUC: ", round(ROC.2$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.2)$CI_AUC), collapse = "-")),
               paste0("5 year AUC:", round(ROC.3$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.3)$CI_AUC), collapse = "-")),
               paste0("8 year AUC:", round(ROC.4$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.4)$CI_AUC), collapse = "-")))
legend("bottomright", legend = my_legend,col=c("red","blue", "forestgreen", "black"),lwd=2, cex=1)
dev.off()


# taylor

taylor = readRDS("/data/github/pca_network/data/Taylor_eSet.RDS")
mat = taylor@assayData$exprs
meta = taylor@phenoData@data
meta = meta[which(!is.na(meta$bcr_status)),]
meta = meta[which(meta$sample_type!="Metastasis"),]
mat = as.data.frame(mat[,rownames(meta)])
meta$days_to_follow_up = floor(meta$time_to_bcr*30.44)

#load("/data/github/pca_network/results/prognostic_model_os.RData")

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat = tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

mat = as.data.frame(scale(mat, scale = T, center = T))

taylor_risk = predict(final_model, newdata = mat, type="lp")

sub_meta = subset(meta, select=c(days_to_follow_up, bcr_status))
mat = cbind(mat, sub_meta)

# Method 1: Regression on PI in validation data
# create new variables to match worked example:
fit_dev = final_model
df_dev = cox
Xb_dev <- model.matrix(fit_dev) %*% coef(fit_dev)
Xbavg_dev <- sum(coef(fit_dev)*fit_dev$means) 
PI_dev <- Xb_dev - Xbavg_dev # centered PI in development dataset (for discrimination later)
df_dev$PI_dev <- PI_dev

df_val <- mat
# fit model in validation data
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~  REG4 + SLC2A4 + CTHRC1 + JAG2, data=df_val) 
Xb_val <- model.matrix(fit_val) %*% coef(fit_dev) # determine Xb in validation data 
PI_val <- Xb_val - Xbavg_dev # center PI by using mean of PI of development data 
df_val$PI_val <- PI_val

fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val, data=df_val)
tab_model(fit_val, transform=NULL, show.r2 = FALSE)

fit_val_test <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val + offset(PI_val), data=df_val) 
tab_model(fit_val_test, transform=NULL, show.r2 = FALSE)

# Method 2: Check model misspecification/fit
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + CTHRC1 + JAG2 + offset(PI_val), data=df_val)
round(2*(diff(fit_val$loglik)), 2) # Chi-squared value
round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3) # p-value


# Method 3: Measures of discrimination
# Derivation
rcorr.cens(-1*PI_dev, Surv(df_dev$days_to_follow_up, df_dev$bcr_status))[1] # Harrell's c-index
GHCI(PI_dev) # Gonen and Heller
# Validation
rcorr.cens(-1*PI_val, Surv(df_val$days_to_follow_up, df_val$bcr_status))[1] # Harrell's c-index
GHCI(df_val$PI_val) # Gonen and Heller

# Method 4: Kaplan-Meier curves for risk groups
# Derivation
df_dev$Risk_group_dev = ifelse(df_dev$PI_dev > median(df_dev$PI_dev), "High", "Low")
surv_object <- Surv(df_dev$days_to_follow_up, df_dev$bcr_status)
res2 = survfit(surv_object ~ Risk_group_dev, data=df_dev)
p1 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_dev,
                 xlab="Time (days)")
p1
# Validation
df_val$Risk_group_val = ifelse(df_val$PI_val > median(df_val$PI_val), "High", "Low")
surv_object <- Surv(df_val$days_to_follow_up, df_val$bcr_status)
res2 = survfit(surv_object ~ Risk_group_val, data=df_val)
p2 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_val,
                 xlab="Time (days)")
p2
library(dplyr)
df_val_KM <- df_val %>%
  group_split(Risk_group_val)
surv_val <- lapply(df_val_KM, function(x){
  obj_surv <- survfit(Surv(days_to_follow_up, bcr_status) ~ 1, data = x)
  data.frame(time=obj_surv$time, surv=obj_surv$surv)
})

p1$plot +
  geom_step(data=surv_val[[1]], aes(x = time, y = surv), linetype=2) +
  geom_step(data=surv_val[[2]], aes(x = time, y = surv), linetype=2)


# MEthod 5 -logrank - skipped

# Method 6. Hazard Ratios between risk groups
# Derivation
df_dev$Risk_group_dev <- factor(df_dev$Risk_group_dev)
df_dev$Risk_group_dev <- relevel(df_dev$Risk_group_dev, ref = "Low")
fit_dev_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_dev), data=df_dev)
tab_model(fit_dev_hr, show.r2 = FALSE)
# validation
df_val$Risk_group_val <- factor(df_val$Risk_group_val)
df_val$Risk_group_val <- relevel(df_val$Risk_group_val, ref = "Low")
fit_val_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_val), data=df_val)
tab_model(fit_val_hr, show.r2 = FALSE)

########################################
# plots for manuscript
########################################

os_mat = df_val
os_mat = os_mat[order(os_mat$PI_val, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$status = factor(os_mat$bcr_status)

pdf("/data/github/pca_network/results/external/Taylor/risk_score_dist.pdf", width = 8, height = 3)
ggscatter(os_mat, y="PI_val", x="patients_inc_risk", color="Risk_group_val", fill="Risk_group_val", 
          ylab="Prognostic index", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = median(os_mat$PI_val), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/external/Taylor/scatter_dfs.pdf", width = 8, height = 3)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(status), x="patients_inc_risk", y="days_to_follow_up", shape="status", ylab = "Time (days)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$Risk_group_val)

col <- c("red", "royalblue3")
names(col) <- c("High", "Low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/external/Taylor/genes_heatmap.pdf", height = 3, width = 8)
pheatmap::pheatmap(t(os_mat[,c("REG4", "SLC2A4", "CTHRC1", "JAG2")]), labels_col = FALSE, color = col_palette, cluster_cols = F, cluster_rows = F,
                   scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
dev.off()

pdf("/data/github/pca_network/results/external/Taylor/high_vs_low.pdf", width=7, height=6)
p2
dev.off()


########################################
# plots for manuscript
########################################


ROC.1 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)
ROC.2  <- timeROC(T=df_val$days_to_follow_up,
                  delta=df_val$bcr_status,
                  marker=df_val$PI_val,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

ROC.4 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*8)),
                 iid=TRUE)
pdf("/data/github/pca_network/results/external/Taylor/ROC_1358.pdf", width=7, height=6)
plot(ROC.1, time=365, title = F, lwd=2)
plot(ROC.2, time=floor(365*3), col="blue", add = T, title=F, lwd=2)
plot(ROC.3, time=floor(365*5), col="forestgreen", add = T, title=F, lwd=2)
plot(ROC.4, time=floor(365*8), col="black", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(ROC.1$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.1)$CI_AUC), collapse = "-")),
               paste0("3 year AUC: ", round(ROC.2$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.2)$CI_AUC), collapse = "-")),
               paste0("5 year AUC:", round(ROC.3$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.3)$CI_AUC), collapse = "-")),
               paste0("8 year AUC:", round(ROC.4$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.4)$CI_AUC), collapse = "-")))
legend("bottomright", legend = my_legend,col=c("red","blue", "forestgreen", "black"),lwd=2, cex=1)
dev.off()



# DKFZ 

dkfz = readRDS("/data/github/pca_network/data/DKFZ_eSet.RDS")
mat = dkfz@assayData$exprs
meta = dkfz@phenoData@data
meta = meta[which(!is.na(meta$bcr_status)),]
mat = as.data.frame(mat[,rownames(meta)])
meta$days_to_follow_up = floor(meta$time_to_bcr*30.44)

#load("/data/github/pca_network/results/prognostic_model_os.RData")

mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat = tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

mat = as.data.frame(scale(mat, scale = T, center = T))

mat = cbind(mat, subset(meta, select=c(days_to_follow_up, bcr_status)))

# Method 1: Regression on PI in validation data
# create new variables to match worked example:
fit_dev = final_model
df_dev = cox
Xb_dev <- model.matrix(fit_dev) %*% coef(fit_dev)
Xbavg_dev <- sum(coef(fit_dev)*fit_dev$means) 
PI_dev <- Xb_dev - Xbavg_dev # centered PI in development dataset (for discrimination later)
df_dev$PI_dev <- PI_dev

df_val <- mat
# fit model in validation data
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~  REG4 + SLC2A4 + CTHRC1 + JAG2, data=df_val) 
Xb_val <- model.matrix(fit_val) %*% coef(fit_dev) # determine Xb in validation data 
PI_val <- Xb_val - Xbavg_dev # center PI by using mean of PI of development data 
df_val$PI_val <- PI_val

fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val, data=df_val)
tab_model(fit_val, transform=NULL, show.r2 = FALSE)

fit_val_test <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val + offset(PI_val), data=df_val) 
tab_model(fit_val_test, transform=NULL, show.r2 = FALSE)

# Method 2: Check model misspecification/fit
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + CTHRC1 + JAG2 + offset(PI_val), data=df_val)
round(2*(diff(fit_val$loglik)), 2) # Chi-squared value
round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3) # p-value


# Method 3: Measures of discrimination
# Derivation
rcorr.cens(-1*PI_dev, Surv(df_dev$days_to_follow_up, df_dev$bcr_status))[1] # Harrell's c-index
GHCI(PI_dev) # Gonen and Heller
# Validation
rcorr.cens(-1*PI_val, Surv(df_val$days_to_follow_up, df_val$bcr_status))[1] # Harrell's c-index
GHCI(df_val$PI_val) # Gonen and Heller

# Method 4: Kaplan-Meier curves for risk groups
# Derivation
df_dev$Risk_group_dev = ifelse(df_dev$PI_dev > median(df_dev$PI_dev), "High", "Low")
surv_object <- Surv(df_dev$days_to_follow_up, df_dev$bcr_status)
res2 = survfit(surv_object ~ Risk_group_dev, data=df_dev)
p1 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_dev,
                 xlab="Time (days)")
p1
# Validation
df_val$Risk_group_val = ifelse(df_val$PI_val > median(df_val$PI_val), "High", "Low")
surv_object <- Surv(df_val$days_to_follow_up, df_val$bcr_status)
res2 = survfit(surv_object ~ Risk_group_val, data=df_val)
p2 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_val,
                 xlab="Time (days)")
p2
library(dplyr)
df_val_KM <- df_val %>%
  group_split(Risk_group_val)
surv_val <- lapply(df_val_KM, function(x){
  obj_surv <- survfit(Surv(days_to_follow_up, bcr_status) ~ 1, data = x)
  data.frame(time=obj_surv$time, surv=obj_surv$surv)
})

p1$plot +
  geom_step(data=surv_val[[1]], aes(x = time, y = surv), linetype=2) +
  geom_step(data=surv_val[[2]], aes(x = time, y = surv), linetype=2)


# MEthod 5 -logrank - skipped

# Method 6. Hazard Ratios between risk groups
# Derivation
df_dev$Risk_group_dev <- factor(df_dev$Risk_group_dev)
df_dev$Risk_group_dev <- relevel(df_dev$Risk_group_dev, ref = "Low")
fit_dev_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_dev), data=df_dev)
tab_model(fit_dev_hr, show.r2 = FALSE)
# validation
df_val$Risk_group_val <- factor(df_val$Risk_group_val)
df_val$Risk_group_val <- relevel(df_val$Risk_group_val, ref = "Low")
fit_val_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_val), data=df_val)
tab_model(fit_val_hr, show.r2 = FALSE)

########################################
# plots for manuscript
########################################

os_mat = df_val
os_mat = os_mat[order(os_mat$PI_val, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$status = factor(os_mat$bcr_status)

pdf("/data/github/pca_network/results/external/DKFZ/risk_score_dist.pdf", width = 8, height = 3)
ggscatter(os_mat, y="PI_val", x="patients_inc_risk", color="Risk_group_val", fill="Risk_group_val", 
          ylab="Prognostic index", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = median(os_mat$PI_val), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/external/DKFZ/scatter_dfs.pdf", width = 8, height = 3)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(status), x="patients_inc_risk", y="days_to_follow_up", shape="status", ylab = "Time (days)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$Risk_group_val)

col <- c("red", "royalblue3")
names(col) <- c("High", "Low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/external/DKFZ/genes_heatmap.pdf", height = 3, width = 8)
pheatmap::pheatmap(t(os_mat[,c("REG4", "SLC2A4", "CTHRC1", "JAG2")]), labels_col = FALSE, color = col_palette, cluster_cols = F, cluster_rows = F,
                   scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
dev.off()

pdf("/data/github/pca_network/results/external/DKFZ/high_vs_low.pdf", width=7, height=6)
p2
dev.off()


########################################
# plots for manuscript
########################################


ROC.1 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)
ROC.2  <- timeROC(T=df_val$days_to_follow_up,
                  delta=df_val$bcr_status,
                  marker=df_val$PI_val,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

pdf("/data/github/pca_network/results/external/DKFZ/ROC_1358.pdf", width=7, height=6)
plot(ROC.1, time=365, title = F, lwd=2)
plot(ROC.2, time=floor(365*3), col="blue", add = T, title=F, lwd=2)
plot(ROC.3, time=floor(365*5), col="forestgreen", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(ROC.1$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.1)$CI_AUC), collapse = "-")),
               paste0("3 year AUC: ", round(ROC.2$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.2)$CI_AUC), collapse = "-")),
               paste0("5 year AUC:", round(ROC.3$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.3)$CI_AUC), collapse = "-")))
legend("bottomright", legend = my_legend,col=c("red","blue", "forestgreen"),lwd=2, cex=1)
dev.off()


# stockholm 

stockholm = readRDS("/data/github/pca_network/data/Stockholm_eSet.RDS")
mat = stockholm@assayData$exprs
meta = stockholm@phenoData@data
meta = meta[which(!is.na(meta$bcr_status)),]
meta$days_to_follow_up = floor(meta$time_to_bcr*30.44)
# remove NA entry
meta = meta[which(!is.na(meta$days_to_follow_up)),]
mat = as.data.frame(mat[,rownames(meta)])



mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
mat = merge(mat, mrna_attributes[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
mat = tibble::column_to_rownames(mat, "external_gene_name")
mat = mat[,2:ncol(mat)]
mat = as.data.frame(t(mat))

mat = as.data.frame(scale(mat, scale = T, center = T))

mat = cbind(mat, subset(meta, select=c(days_to_follow_up, bcr_status)))

# Method 1: Regression on PI in validation data
# create new variables to match worked example:
fit_dev = final_model
df_dev = cox
Xb_dev <- model.matrix(fit_dev) %*% coef(fit_dev)
Xbavg_dev <- sum(coef(fit_dev)*fit_dev$means) 
PI_dev <- Xb_dev - Xbavg_dev # centered PI in development dataset (for discrimination later)
df_dev$PI_dev <- PI_dev

df_val <- mat
# fit model in validation data
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~  REG4 + SLC2A4 + CTHRC1 + JAG2, data=df_val) 
Xb_val <- model.matrix(fit_val) %*% coef(fit_dev) # determine Xb in validation data 
PI_val <- Xb_val - Xbavg_dev # center PI by using mean of PI of development data 
df_val$PI_val <- PI_val

fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val, data=df_val)
tab_model(fit_val, transform=NULL, show.r2 = FALSE)

fit_val_test <- coxph(Surv(days_to_follow_up, bcr_status) ~ PI_val + offset(PI_val), data=df_val) 
tab_model(fit_val_test, transform=NULL, show.r2 = FALSE)

# Method 2: Check model misspecification/fit
fit_val <- coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + CTHRC1 + JAG2 + offset(PI_val), data=df_val)
round(2*(diff(fit_val$loglik)), 2) # Chi-squared value
round(1-pchisq(2*(diff(fit_val$loglik)), 8), 3) # p-value


# Method 3: Measures of discrimination
# Derivation
rcorr.cens(-1*PI_dev, Surv(df_dev$days_to_follow_up, df_dev$bcr_status))[1] # Harrell's c-index
GHCI(PI_dev) # Gonen and Heller
# Validation
rcorr.cens(-1*PI_val, Surv(df_val$days_to_follow_up, df_val$bcr_status))[1] # Harrell's c-index
GHCI(df_val$PI_val) # Gonen and Heller

# Method 4: Kaplan-Meier curves for risk groups
# Derivation
df_dev$Risk_group_dev = ifelse(df_dev$PI_dev > median(df_dev$PI_dev), "High", "Low")
surv_object <- Surv(df_dev$days_to_follow_up, df_dev$bcr_status)
res2 = survfit(surv_object ~ Risk_group_dev, data=df_dev)
p1 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_dev,
                 xlab="Time (days)")
p1
# Validation
df_val$Risk_group_val = ifelse(df_val$PI_val > median(df_val$PI_val), "High", "Low")
surv_object <- Surv(df_val$days_to_follow_up, df_val$bcr_status)
res2 = survfit(surv_object ~ Risk_group_val, data=df_val)
p2 <- ggsurvplot(res2,
                 pval = TRUE, conf.int = F,
                 risk.table = T, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 surv.median.line = "none", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 palette = c("red1", "royalblue3"),
                 data=df_val,
                 xlab="Time (days)")
p2
library(dplyr)
df_val_KM <- df_val %>%
  group_split(Risk_group_val)
surv_val <- lapply(df_val_KM, function(x){
  obj_surv <- survfit(Surv(days_to_follow_up, bcr_status) ~ 1, data = x)
  data.frame(time=obj_surv$time, surv=obj_surv$surv)
})

p1$plot +
  geom_step(data=surv_val[[1]], aes(x = time, y = surv), linetype=2) +
  geom_step(data=surv_val[[2]], aes(x = time, y = surv), linetype=2)


# MEthod 5 -logrank - skipped

# Method 6. Hazard Ratios between risk groups
# Derivation
df_dev$Risk_group_dev <- factor(df_dev$Risk_group_dev)
df_dev$Risk_group_dev <- relevel(df_dev$Risk_group_dev, ref = "Low")
fit_dev_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_dev), data=df_dev)
tab_model(fit_dev_hr, show.r2 = FALSE)
# validation
df_val$Risk_group_val <- factor(df_val$Risk_group_val)
df_val$Risk_group_val <- relevel(df_val$Risk_group_val, ref = "Low")
fit_val_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ factor(Risk_group_val), data=df_val)
tab_model(fit_val_hr, show.r2 = FALSE)



########################################
# plots for manuscript
########################################

os_mat = df_val
os_mat = os_mat[order(os_mat$PI_val, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$status = factor(os_mat$bcr_status)

pdf("/data/github/pca_network/results/external/Stockholm/risk_score_dist.pdf", width = 8, height = 3)
ggscatter(os_mat, y="PI_val", x="patients_inc_risk", color="Risk_group_val", fill="Risk_group_val", 
          ylab="Prognostic index", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = median(os_mat$PI_val), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/external/Stockholm/scatter_dfs.pdf", width = 8, height = 3)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(status), x="patients_inc_risk", y="days_to_follow_up", shape="status", ylab = "Time (days)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = median(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$Risk_group_val)

col <- c("red", "royalblue3")
names(col) <- c("High", "Low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/external/Stockholm/genes_heatmap.pdf", height = 3, width = 8)
pheatmap::pheatmap(t(os_mat[,c("REG4", "SLC2A4", "CTHRC1", "JAG2")]), labels_col = FALSE, color = col_palette, cluster_cols = F, cluster_rows = F,
                   scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
dev.off()

pdf("/data/github/pca_network/results/external/Stockholm/high_vs_low.pdf", width=7, height=6)
p2
dev.off()


########################################
# plots for manuscript
########################################


ROC.1 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)
ROC.2  <- timeROC(T=df_val$days_to_follow_up,
                  delta=df_val$bcr_status,
                  marker=df_val$PI_val,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

ROC.4 <- timeROC(T=df_val$days_to_follow_up,
                 delta=df_val$bcr_status,
                 marker=df_val$PI_val,
                 cause=1,weighting="marginal",
                 times=c(floor(365*8)),
                 iid=TRUE)
pdf("/data/github/pca_network/results/external/Stockholm/ROC_1358.pdf", width=7, height=6)
plot(ROC.1, time=365, title = F, lwd=2)
plot(ROC.2, time=floor(365*3), col="blue", add = T, title=F, lwd=2)
plot(ROC.3, time=floor(365*5), col="forestgreen", add = T, title=F, lwd=2)
plot(ROC.4, time=floor(365*8), col="black", add = T, title=F, lwd=2)
my_legend =  c(paste0("1 year AUC: ", round(ROC.1$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.1)$CI_AUC), collapse = "-")),
               paste0("3 year AUC: ", round(ROC.2$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.2)$CI_AUC), collapse = "-")),
               paste0("5 year AUC:", round(ROC.3$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.3)$CI_AUC), collapse = "-")),
               paste0("8 year AUC:", round(ROC.4$AUC[2],2), "; 95% CI: ", paste(as.numeric(confint(ROC.4)$CI_AUC), collapse = "-")))
legend("bottomright", legend = my_legend,col=c("red","blue", "forestgreen", "black"),lwd=2, cex=1)
dev.off()
