#!/usr/bin/env Rscript

# ROC plots, scatterhist + survival for OS + DFS

load("/data/github/pca_network/results/TCGA_survival.RData")
signf= signf_os
# LASSO cox on signif genes 
library(glmnet)
library(dplyr)
library(ggpubr)
library(survival)
library(survminer)


# using gene counts from STAR (not FPKM protein ATLAS) to run LASSO 
# Run TCGA_mrna.R to produce 'scaled' and mrna_attributes. 
scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
colnames(scaled) = gsub("\\.", "-", colnames(scaled))
scaled = scaled[,which(colnames(scaled) %in% colnames(atlas_mat))]
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
version_to_id = merge(signf, mrna_attributes, by.x="ensembl", by.y="ensembl_gene_id")
sub_atlas = scaled[which(rownames(scaled) %in% version_to_id$ensembl_gene_id_version),]
sub_atlas = merge(sub_atlas, version_to_id[,5:6], by.x=0, by.y="ensembl_gene_id_version")
sub_atlas = tibble::column_to_rownames(sub_atlas, "external_gene_name")
sub_atlas = sub_atlas[,c(2:ncol(sub_atlas))]
#sub_atlas = sub_atlas[,rownames(atlas_meta)]

# Subset metadata for cox OS (time to follow, vital status) this can change if you edit metadata cols
sub_meta = atlas_meta[,c(9,8)]
sub_meta$vital = factor(ifelse(sub_meta$vital=="Dead", 1,0))
colnames(sub_meta) = c("time", "status")

x = as.matrix(t(sub_atlas))
y = sub_meta
y$status = as.numeric(as.character(y$status))

y = as.matrix(y)

# Run and save to RData
set.seed(123)
cv.fit <- cv.glmnet(x, y, family="cox", alpha=1, maxit = 1000, lambda = NULL, type.measure = "deviance")
fit = glmnet(x, y, family = "cox", alpha=1, maxit = 1000, lambda=NULL)
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Genes <- Coefficients@Dimnames[[1]][Active.Index]
# save(cv.fit, fit, Active.Genes, Active.Coefficients, file="/data/github/pca_network/results/LASSO_cox.RData")
# pdf("/data/github/pca_network/results/lasso_plot_fit.pdf", width = 8, height = 8)
# plot(cv.fit$glmnet.fit, "lambda", label=T)
# dev.off()
# pdf("/data/github/pca_network/results/lasso_log_lambda.pdf", width = 8, height = 8)
# plot(cv.fit)
# dev.off()
load("/data/github/pca_network/results/LASSO_cox.RData")


# Reload x if re-running code! 
x = as.data.frame(x)
x = x[,colnames(x) %in% Active.Genes]
risk_scores <- rowSums(Active.Coefficients * x)
x = cbind(x, sub_meta)
x$risk_score = risk_scores
x$risk_score_cat = ifelse(x$risk_score > mean(x$risk_score), "high", "low")

surv_object <- Surv(x$time, as.numeric(x$status))

res = coxph(Surv(x$time, as.numeric(x$status)) ~ risk_score_cat, data=x)
res2 = survfit(surv_object ~ risk_score_cat, data=x)
logrank = survdiff(surv_object ~ risk_score_cat, data=x)


# Kaplan Meier OS ~ risk score 
#pdf("/data/github/pca_network/results/tcga_os_surv_plots/risk_surv.pdf", height=8,width=8)
ggsurvplot(res2,
           pval = TRUE, conf.int = T,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=x,xlab="Time (days)")
#dev.off()


# ROC OS 1 , 3, 5, 10 year 
library(survivalROC)
nobs <- length(x$risk_score)
roc <- survivalROC(Stime=x$time,
                   status=x$status,
                   marker = x$risk_score,
                   method = 'NNE',
                   predict.time = 365,
                   span = 0.25*nobs^(-0.20))

roc2 <- survivalROC(Stime=x$time,
                    status=x$status,
                    marker = x$risk_score,
                    method = 'NNE',
                    predict.time = 365*3,
                    span = 0.25*nobs^(-0.20))

roc3 <- survivalROC(Stime=x$time,
                    status=x$status,
                    marker = x$risk_score,
                    method = 'NNE',
                    predict.time = 365*5,
                    span = 0.25*nobs^(-0.20))

roc4 <- survivalROC(Stime=x$time,
                    status=x$status,
                    marker = x$risk_score,
                    method = 'NNE',
                    predict.time = 365*10,
                    span = 0.25*nobs^(-0.20)) 

pdf("/data/github/pca_network/results/tcga_os_surv_plots/10_year_roc.pdf", width=8, height = 8)
plot(roc$FP, roc$TP, type="S", lwd=3, col=alpha("royalblue",0.7), xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(roc$AUC,3)),
     ylab="TP")
lines(roc2$FP, roc2$TP, type = "S", col=alpha("red",0.7), lwd=3)
lines(roc3$FP, roc3$TP, type = "S", col=alpha("grey9",0.7), lwd=3)
lines(roc4$FP, roc$TP, type="S", col=alpha("green",0.7), lwd=3)
abline(0,1, lwd=2)


# ROC OS overall 
risk_roc = roc(x$status, x$risk_score)
pdf("/data/github/pca_network/results/tcga_os_surv_plots/ROC_risk_score.pdf", width = 8, height = 8)
pROC::plot.roc(risk_roc, col="red", xlim = c(0.9,0.1), add = F)
legend("bottomright", legend = paste("AUC = ", round(risk_roc$auc, 3)))
dev.off()

plot_mt = x
plot_mt = plot_mt[order(plot_mt$risk_score, decreasing = F),]
plot_mt$patients_inc_risk = seq(1, nrow(plot_mt))
plot_mt$years = plot_mt$time/365
plot_mt$status = ifelse(plot_mt$status == 1, "Dead", "Alive")

# combine in inkscape.
pdf("/data/github/pca_network/results/tcga_os_surv_plots/scatter_surv_risk_ranked.pdf", width = 8, height = 2)
ggpubr::ggscatter(plot_mt, y="risk_score", x="patients_inc_risk", color="risk_score_cat", fill="risk_score_cat", 
                  ylab="Risk Score", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = mean(plot_mt$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = mean(plot_mt$risk_score), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/tcga_os_surv_plots/scatter_status_risk_ranked.pdf", width = 8, height = 4)
ggpubr::ggscatter(plot_mt %>% dplyr::arrange(status), x="patients_inc_risk", y="years", shape="status", ylab = "Time (years)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = mean(plot_mt$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

library(RColorBrewer)
col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(plot_mt),
                     Group = plot_mt$risk_score_cat)

col <- c("red", "royalblue3")
names(col) <- c("high", "low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/tcga_os_surv_plots/lasso_genes_heatmap.pdf", height = 4, width = 8)
pheatmap::pheatmap(t(plot_mt[,1:9]), labels_col = FALSE, color = col_palette, cluster_cols = F, 
                   scale = "column")
dev.off()

model_p <- data.frame(genes = Active.Genes, coef = Active.Coefficients)
model_p$coef <- round(model_p$coef, 5)
pdf("/data/github/pca_network/results/tcga_os_surv_plots/lasso_coefficients.pdf", width=8, height=8)
ggpubr::ggbarplot(model_p %>% dplyr::arrange((coef)), label = T, x="genes", y="coef", lab.vjust = 0.5, lab.hjust = 1.1, lab.col = "white",
                  fill = "royalblue", title = "Coefficients of LASSO Cox selected genes", ggtheme = theme_linedraw()) + rotate()
dev.off()




## multivariate COX ? missing M/N/T Stage for pathological or clinical - which to use. 

