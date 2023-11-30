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
library(boot)

#########################################################################
# load ceRNA network genes
#########################################################################
network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
genes = unique(network$mrna)

#########################################################################
# load TCGA metadata from Protein Atlas.
# Supplement additional metadata from other sources
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
resec = read.csv("/data/github/pca_network/data/prad_tcga_clinical_data.tsv", header=T, sep="\t")
resec$Sample.ID = paste(resec$Sample.ID, "A", sep="")
table(resec$Sample.ID %in% atlas_meta$sample)
atlas_meta = merge(atlas_meta, subset(resec, select=c(Sample.ID, Surgical.Margin.Resection.Status)), by.x="sample", by.y="Sample.ID")
atlas_meta$years_to_follow_up <- as.character(floor(atlas_meta$years_to_follow_up))

#########################################################################
# Use logcpm STAR counts for CoxPH/Kaplan Meier
#########################################################################
logcpm = read.csv("/data/github/pca_network/results/TCGA_mrna_logcpm.txt", header=T, row.names = 1, sep="\t")
colnames(logcpm) = gsub("\\.", "-", colnames(logcpm))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% genes),]
logcpm = merge(logcpm, mrna_attributes[,1:2], by.x=0, by.y="ensembl_gene_id_version")
logcpm = tibble::column_to_rownames(logcpm, "external_gene_name")
logcpm = logcpm[,c(2:ncol(logcpm))]
logcpm = logcpm[,atlas_meta$Row.names]

#########################################################################
# Univariate Cox proportional hazards regression:
# Find ceRNA mRNAs * with prognosis of PCa (disease free survival)
#########################################################################
univ_mat = as.data.frame(t(scale(t(logcpm), scale = T, center = T)))
univ_mat = as.data.frame(t(univ_mat))
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

res = res[which(res$P.adjust <= 0.05),]

genes = res$Variable

write.csv(res, "/data/github/pca_network/results/TCGA_DFS/ceRNA_genes_univariate_cox.csv", quote=F, row.names = FALSE)

##########################################################################
# Schoenfelds test - remove genes that violate Coxph assumptions
##########################################################################
univ_mat = univ_mat[,c(genes, "days_to_follow_up", "bcr_status")]
mult_cx = coxph(Surv(days_to_follow_up, bcr_status) ~  ., data=univ_mat)
test.cp = cox.zph(mult_cx, transform = "km")
test.cp

remove = rownames(test.cp$table)[which(test.cp$table[,"p"] <= 0.05)]
# save for plot
sec31b = ggcoxzph2(test.cp)[2]
slc22a3 = ggcoxzph2(test.cp)[12]

genes = setdiff(genes, remove)

univ_mat = univ_mat[,c(genes, "days_to_follow_up", "bcr_status")]
mult_cx = coxph(Surv(days_to_follow_up, bcr_status) ~  ., data=univ_mat)
test.cp = cox.zph(mult_cx)
test.cp

remove = rownames(test.cp$table)[which(test.cp$table[,"p"] <= 0.05)]

srrd = ggcoxzph2(test.cp)[2]
gad1 = ggcoxzph2(test.cp)[8]

genes = setdiff(genes, remove)

univ_mat = univ_mat[,c(genes, "days_to_follow_up", "bcr_status")]
mult_cx = coxph(Surv(days_to_follow_up, bcr_status) ~  ., data=univ_mat)
test.cp = cox.zph(mult_cx)
test.cp

# use altered coxzph plot for better y-axes, easier to see curves in plots.
source("https://raw.githubusercontent.com/BarryDigby/pca_network/main/data/modified_coxzph.R")

pdf("/data/github/pca_network/results/TCGA_DFS/schoenfeld_residuals.pdf", height=12,width=20)
gad1
sec31b
srrd
slc22a3
ggcoxzph2(test.cp)
dev.off()

# save schoenfeld residuals for supplementary
write.csv(test.cp$table, "/data/github/pca_network/results/TCGA_DFS/schoenfeld_residuals.csv", quote = F, row.names = T)

##########################################################################
# Plot hazard ratios
##########################################################################
univ_formulas <- sapply(genes, function(x) as.formula(paste('Surv(days_to_follow_up, bcr_status)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = univ_mat)})

#pdf("/data/github/pca_network/results/TCGA_DFS3/univariate_forest_model.pdf", width=8, height=5)
forestmodel::forest_model(model_list = univ_models,covariates = genes,merge_models =T)


##########################################################################
# Cox ph multiple regression
##########################################################################
set.seed(123)
mult_cox = coxph(Surv(days_to_follow_up, bcr_status) ~., data=univ_mat)
summary(mult_cox)

step = stepAIC(mult_cox, direction = "both")
summary(step)

step_df = summary(step)$coeff
step_df = cbind(step_df, confint(step))

write.csv(step_df, "/data/github/pca_network/results/TCGA_DFS/stepAIC_summary.csv", row.names = TRUE, quote = FALSE)

forestmodel::forest_model(coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + PAPSS1 + TRIM13 + INPP5E + JAG2 + CTHRC1, data=univ_mat))

final_model = coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + JAG2 + CTHRC1, data=univ_mat)
summary(final_model)
risk_scores = predict(final_model, newdata = univ_mat, type="lp")

model_df = summary(final_model)$coeff
genes = rownames(model_df)

##########################################################################
#  test linearity - read up on if these are useful for supp materials
##########################################################################

ggcoxdiagnostics(final_model, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

ggcoxdiagnostics(final_model, type = "schoenfeld",
                 linear.predictions = FALSE, ggtheme = theme_bw())

ggcoxdiagnostics(final_model, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

##########################################################################
# Stratify patients 
##########################################################################
cox = as.data.frame(univ_mat[,genes])
cox$risk_score = risk_scores
cox = cbind(cox, atlas_meta[,c("days_to_follow_up", "bcr_status")])
q = quantile(risk_scores, c(.25, .5, .75), type=1)
cox$risk_category = ifelse(cox$risk_score >= median(cox$risk_score), "High risk", "Low risk")
surv_object <- Surv(cox$days_to_follow_up, cox$bcr_status)
res2 = survfit(surv_object ~ risk_category, data=cox)

pdf("/data/github/pca_network/results/TCGA_DFS/high_vs_low.pdf", width=7, height=6)
ggsurvplot(res2,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "none", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=cox,
           xlab="Time (days)")
dev.off()

##########################################################################
# Heatmap, scatterplots 
##########################################################################
os_mat = cox
os_mat = os_mat[order(os_mat$risk_score, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$status = factor(os_mat$bcr_status)
os_mat$years = atlas_meta$years_to_follow_up

pdf("/data/github/pca_network/results/TCGA_DFS/risk_score_dist.pdf", width = 8, height = 3)
ggscatter(os_mat, y="risk_score", x="patients_inc_risk", color="risk_category", fill="risk_category", 
          ylab="Prognostic index", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = mean(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = mean(os_mat$risk_score), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/TCGA_DFS/scatter_dfs.pdf", width = 8, height = 3)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(status), x="patients_inc_risk", y="days_to_follow_up", shape="status", ylab = "Time (days)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = mean(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$risk_category)

col <- c("red", "royalblue3")
names(col) <- c("High risk", "Low risk")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/TCGA_DFS/genes_heatmap.pdf", height = 3, width = 8)
p = pheatmap::pheatmap(t(os_mat[,1:4]), labels_col = FALSE, color = col_palette, cluster_cols = F, 
                       scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
print(p)
dev.off()

##########################################################################
# Prediction accuracy 1, 3, 5 years
##########################################################################

ROC.1 <- timeROC(T=cox$days_to_follow_up,
                 delta=cox$bcr_status,
                 marker=cox$risk_score,
                 cause=1,weighting="marginal",
                 times=c(365),
                 iid=TRUE)
ROC.2  <- timeROC(T=cox$days_to_follow_up,
                  delta=cox$bcr_status,
                  marker=cox$risk_score,
                  cause=1,weighting="marginal",
                  times=c(floor(365*3)),
                  iid=TRUE)

ROC.3 <- timeROC(T=cox$days_to_follow_up,
                 delta=cox$bcr_status,
                 marker=cox$risk_score,
                 cause=1,weighting="marginal",
                 times=c(floor(365*5)),
                 iid=TRUE)

ROC.4 <- timeROC(T=cox$days_to_follow_up,
                 delta=cox$bcr_status,
                 marker=cox$risk_score,
                 cause=1,weighting="marginal",
                 times=c(floor(365*8)),
                 iid=TRUE)

pdf("/data/github/pca_network/results/TCGA_DFS/ROC_135.pdf", width=7, height=6)
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

save.image(file="/data/github/pca_network/results/TCGA_DFS/model.RData")
