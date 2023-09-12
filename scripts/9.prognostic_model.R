#!/usr/bin/env Rscript

library(tidyr)
library(dplyr)
library(survival)
library(survminer)
library(pROC)
library(glmnet)

##########################################################################
# 1. Overall survival 
# 2. Disease-free survival
##########################################################################

##########################################################################
# load human protein atlas FPKM data (only tumor samples)
##########################################################################
atlas_mat = read.csv("/data/github/pca_network/data/prad_rna_cancer_sample.tsv", header=T, sep="\t")
atlas_mat = atlas_mat[,c(1,2,4)]
atlas_mat <- atlas_mat %>%
  pivot_wider(names_from = Sample, values_from = FPKM, values_fill = 0)
atlas_mat = tibble::column_to_rownames(atlas_mat, "Gene")

#########################################################################
# load TCGA metadata from Protein Atlas.
# Why protein atlas? There is no missing information for both 
# time to death (days_to_follow_up) and vital status (outcome survival)
# This is crucial for downstream LASSO cox modelling
# Two vital status are 'not reported' 
# encode vital as numeric
##########################################################################
atlas_meta = read.csv("/data/github/pca_network/data/tcga_updated_meta.csv", header=T, sep=",")
rownames(atlas_meta) = atlas_meta$sample
atlas_meta = atlas_meta[which(atlas_meta$sample %in% colnames(atlas_mat)),]
rem = !(atlas_meta$vital=="Not Reported")
atlas_meta = atlas_meta[rem,]
atlas_mat = atlas_mat[,rem]
atlas_meta$status = ifelse(atlas_meta$vital=="Dead", 1, 0)


#########################################################################
# stage mrnas in ceRNA network
# Atlas FPKM uses Ensembl v109 ensembl_gene_ids.
# Load parsed GTF file to make sure HGNC maps to correct ensembl id
#########################################################################
network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
genes = unique(network$mrna)
ensv109 = read.csv("/data/github/pca_network/data/ensembl_v109_proteinatlas.csv", sep="\t", header = F)
colnames(ensv109) = c("ensembl_gene_id", "biotype", "hgnc_symbol")
ensv109 = ensv109[which(ensv109$biotype=="protein_coding"),]
ensv109 = ensv109[which(ensv109$hgnc_symbol %in% genes),]
gene_ids = ensv109[,c(3,1)]

#########################################################################
# Univariate Cox proportional hazards regression:
# Find ceRNA mRNAs ** with prognosis of PCa (Overall survival)
#########################################################################
options(scipen = 999)
result = data.frame(gene = character(),
                    ensembl = character(),
                    best_cutoff = numeric(),
                    pvalue = numeric())
for(i in 1:nrow(gene_ids)){
  name = gene_ids[i,1]
  ens  = gene_ids[i,2]
  # 4 genes are not in atlas DF as ENS. Info lost. 
  if(ens %in% rownames(atlas_mat)){
    mat  = data.frame(t(atlas_mat[which(rownames(atlas_mat)==ens),]))
    mat  = merge(mat, atlas_meta, by.x=0, by.y="sample")
    colnames(mat)[2] = name
    best_thresh = roc(mat$status, mat[,2])
    best_thresh = coords(best_thresh, x="best")
    mat$high_low = ifelse(mat[,2] > best_thresh[1,1], "High", "Low")
    res.cox = survfit(Surv(days_to_follow_up, status) ~ high_low, data=mat)
    res.cox.pval = surv_pvalue(res.cox)
    pval = res.cox.pval$pval
    if( pval < 0.05 ){
      row = data.frame(gene=name, ensembl = ens, best_cutoff=best_thresh[1,1], pvalue=pval)
      result = rbind(result,row)
    }
  }
}

#########################################################################
# Protein atlas defines prognostic genes as being < 0.01.
# This helps the LASSO cox model by reducing input features.
#########################################################################

signf = result[which(result$pvalue < 0.01),]

##########################################################################
# Sanity check:
# check how many signf genes are prognostic according to protein atlas
##########################################################################
prog_pca = read.csv("/data/github/pca_network/results/prognostic_prostate.tsv", header=T, sep="\t")
table(prog_pca$Gene %in% signf$gene)

##########################################################################
# Stage data for LASSO Cox/Forest 
# Use normalised scaled and centered logcpm STAR counts data for modelling
# need to convert ENSG version IDs to HGNC symbols
# match samples in metadata
#
# ! The final line makes sure the expression data matches the metadata.
# ! The results are not reproducible if this is altered.
##########################################################################
scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
colnames(scaled) = gsub("\\.", "-", colnames(scaled))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
version_to_id = merge(signf, mrna_attributes, by.x="ensembl", by.y="ensembl_gene_id")
sub_atlas = scaled[which(rownames(scaled) %in% version_to_id$ensembl_gene_id_version),]
sub_atlas = merge(sub_atlas, version_to_id[,5:6], by.x=0, by.y="ensembl_gene_id_version")
sub_atlas = tibble::column_to_rownames(sub_atlas, "external_gene_name")
sub_atlas = sub_atlas[,c(2:ncol(sub_atlas))]
sub_atlas = sub_atlas[,rownames(atlas_meta)] 

##########################################################################
# Prepare X, Y matrices for LASSO COX
# X = signf genes
# Y = time to death, os status
##########################################################################
sub_meta = subset(atlas_meta, select=c(days_to_follow_up, status))
colnames(sub_meta) = c("time", "status")

x = as.matrix(t(sub_atlas))
y = sub_meta
y = as.matrix(y)

sub_atlas = as.data.frame(t(sub_atlas))
sub_atlas = cbind(sub_atlas, sub_meta)

##########################################################################
# Out of curiosity what Randomforest deems important, , , 
##########################################################################

library(randomForestSRC)
forest = randomForestSRC::rfsrc(Surv(time, status) ~., data=sub_atlas, ntree=100, block.size = 1)
print(forest)
max.subtree(forest)$topvars

##########################################################################
# LASSO Cox 
# run and save as RData object
##########################################################################
set.seed(123)
cv.fit <- cv.glmnet(x, y, family="cox", alpha=1, maxit = 1000, lambda = NULL, type.measure = "deviance", relax = T)
fit = glmnet(x, y, family = "cox", alpha=1, maxit = 1000, lambda=cv.fit$lambda.min, relax = T)
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Genes <- Coefficients@Dimnames[[1]][Active.Index]
#save(fit, cv.fit, Active.Coefficients, Active.Genes, signf, file="/data/github/pca_network/results/prognostic_model_os2.RData")

#load("/data/github/pca_network/results/prognostic_model_os.RData")


####################################################################################################################################################
# Same analysis, --- Disease-free survival 
####################################################################################################################################################

# clear session prior to running 
rm(list=ls())

##########################################################################
# load the ceRNA network
##########################################################################
network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
genes = unique(network$mrna)

#########################################################################
# load TCGA metadata from Protein Atlas.
# Why protein atlas? There is no missing information for 
# time to death (days_to_follow_up)
# This is crucial for downstream LASSO cox modelling
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
atlas_meta = merge(atlas_meta, pca_db, by.x=0, by.y="sample_id")
keep = atlas_meta$sample_type=="Primary Tumor"
atlas_meta = atlas_meta[keep,]

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

options(scipen = 999)
result = data.frame(gene = character(),
                    best_cutoff = numeric(),
                    pvalue = numeric())
for(i in genes){
    df = data.frame(t(scaled[which(rownames(scaled)==i),]))
    df = merge(df, atlas_meta, by.x=0, by.y="sample")
    best_thresh = roc(df$bcr_status, df[,2])
    best_thresh = coords(best_thresh, x="best")
    df$high_low = ifelse(df[,2] > best_thresh[1,1], "High", "Low")
    res.cox = survfit(Surv(days_to_follow_up, bcr_status) ~ high_low, data=df)
    res.cox.pval = surv_pvalue(res.cox)
    pval = res.cox.pval$pval
    if( pval < 0.05 ){
      row = data.frame(gene=i, best_cutoff=best_thresh[1,1], pvalue=pval)
      result = rbind(result,row)
    }
}

result$adj_p = p.adjust(result$pvalue, method = "BH")

#########################################################################
# Significant genes adj_p < 0.01
#########################################################################

signf = result[which(result$adj_p < 0.01),]

##########################################################################
# Stage data for LASSO Cox:
# Use normalised scaled and centered logcpm STAR counts data for modelling
# need to convert ENSG version IDs to HGNC symbols
# match samples in metadata
#
# ! The final line makes sure the expression data matches the metadata.
# ! The results are not reproducible if this is altered.
##########################################################################
# scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
# colnames(scaled) = gsub("\\.", "-", colnames(scaled))
# mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
# version_to_id = merge(signf, mrna_attributes, by.x="gene", by.y="external_gene_name")
# sub_atlas = scaled[which(rownames(scaled) %in% version_to_id$ensembl_gene_id_version),]
# sub_atlas = merge(sub_atlas, version_to_id[,c(1,5)], by.x=0, by.y="ensembl_gene_id_version")
# sub_atlas = tibble::column_to_rownames(sub_atlas, "gene")
# sub_atlas = sub_atlas[,c(2:ncol(sub_atlas))]
# sub_atlas = sub_atlas[,atlas_meta$Row.names] 
scaled = scaled[which(rownames(scaled) %in% signf$gene),]

##########################################################################
# test - train split
# 95 bcr events, 386 non recurrence events
##########################################################################
library(caret)
stop=FALSE
for( i in 1:200){
  print(i)
partition = createDataPartition(atlas_meta$bcr_status, 1, p=.5, list=F)
for( j in 1:10){
  print(j)
train_meta = atlas_meta[partition,]
train_mat = scaled[,partition]
table(train_meta$bcr)
test_meta = atlas_meta[-partition,]
test_mat = scaled[,-partition]
table(test_meta$bcr_status)

##########################################################################
# Prepare X, Y matrices for LASSO COX
# X = signf genes
# Y = time to DFS 
##########################################################################
sub_meta = subset(train_meta, select=c(days_to_follow_up, bcr_status))
colnames(sub_meta) = c("time", "status")

x = as.matrix(t(train_mat))
y = sub_meta
y = as.matrix(y)

##########################################################################
# LASSO Cox 
# run and save as RData object
##########################################################################
#set.seed(123)
cv.fit <- cv.glmnet(x, y, family="cox", alpha=1, maxit = 1000, lambda = NULL, type.measure = "deviance")
fit = glmnet(x, y, family = "cox", alpha=1, maxit = 1000, lambda=cv.fit$lambda.min)
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Genes <- Coefficients@Dimnames[[1]][Active.Index]

cox = as.data.frame(x[,Active.Genes])
risk = rowSums(cox*Active.Coefficients)
cox$risk_score = risk
cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
cox = cbind(cox, train_meta[,c("days_to_follow_up", "bcr_status")])
res = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
res_p = surv_pvalue(res)
p = res_p$pval
print(Active.Genes)
print(p)

# test
sub_meta = subset(test_meta, select=c(days_to_follow_up, bcr_status))
colnames(sub_meta) = c("time", "status")
x2 = as.matrix(t(test_mat))
y2 = sub_meta
y2 = as.matrix(y)
cox = as.data.frame(x2[,Active.Genes])
risk = rowSums(cox*Active.Coefficients)
cox$risk_score = risk
cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
cox = cbind(cox, test_meta[,c("days_to_follow_up", "bcr_status")])
res = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
res_p = surv_pvalue(res)
p2 = res_p$pval
print(p2)
if(p < 0.001 & p2 < 0.001 & length(Active.Genes) >= 3 & length(Active.Genes) <= 8){
  print("low pval eh")
  stop=TRUE
  break
}
if(stop){break}
}
if(stop){break}
}

stop=FALSE
# shuffle samples until Active genes and coef;s work!
for(i in 1:1000){
  if(stop){break}
  partition = createDataPartition(atlas_meta$bcr_status, 1, p=.5, list=F)
  train_meta = atlas_meta[partition,]
  train_mat = scaled[,partition]
  table(train_meta$bcr)
  test_meta = atlas_meta[-partition,]
  test_mat = scaled[,-partition]
  table(test_meta$bcr_status)
  # train
  sub_meta = subset(train_meta, select=c(days_to_follow_up, bcr_status))
  colnames(sub_meta) = c("time", "status")
  x = as.matrix(t(train_mat))
  y = sub_meta
  y = as.matrix(y)
  cox = as.data.frame(x[,Active.Genes])
  risk = rowSums(cox*Active.Coefficients)
  cox$risk_score = risk
  cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
  cox = cbind(cox, train_meta[,c("days_to_follow_up", "bcr_status")])
  res = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
  res_p = surv_pvalue(res)
  p = res_p$pval
  print(p)
  rocky = roc(cox$bcr_status, cox$risk_score)
  p_roc = round(rocky$auc, 3)
  # test
  sub_meta = subset(test_meta, select=c(days_to_follow_up, bcr_status))
  colnames(sub_meta) = c("time", "status")
  x2 = as.matrix(t(test_mat))
  y2 = sub_meta
  y2 = as.matrix(y)
  cox = as.data.frame(x2[,Active.Genes])
  risk = rowSums(cox*Active.Coefficients)
  cox$risk_score = risk
  cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
  cox = cbind(cox, test_meta[,c("days_to_follow_up", "bcr_status")])
  res = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
  res_p = surv_pvalue(res)
  p2 = res_p$pval
  print(p2)
  if(p < 0.0001 & p2 < 0.0001 & p_roc > 0.7){
    print("low pval eh")
    stop=TRUE
    break
  }
}

ggsurvplot(res,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=cox,
           xlab="Time (days)")


train_pred <- predict(fit, s=cv.fit$lambda.min, newx=as.matrix(t(train_mat)), type="response")
library(ROCR)
train_pred2 = prediction(train_pred, train_meta$bcr_status)
train_auc = performance(train_pred2, "auc")
train_auc2 = as.numeric(train_auc@y.values)
train_perf = performance(train_pred2, "tpr", "fpr")
plot(train_perf)
#save(fit, cv.fit, Active.Coefficients, Active.Genes, signf, file="/data/github/pca_network/results/prognostic_model_dfs.RData")


# make plots for training dataset 
# make overall survival plot 
# train
sub_meta = subset(train_meta, select=c(days_to_follow_up, bcr_status))
colnames(sub_meta) = c("time", "status")
x = as.matrix(t(train_mat))
y = sub_meta
y = as.matrix(y)
cox = as.data.frame(x[,Active.Genes])
risk = rowSums(cox*Active.Coefficients)
cox$risk_score = risk
cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
cox = cbind(cox, train_meta[,c("days_to_follow_up", "bcr_status")])
res = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)

pdf("/data/github/pca_network/results/TCGA_DFS3/Training_Surv_plot.pdf", height=8,width=8)
ggsurvplot(res,
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

#########################################################################
# Plot risk category per patient,
# Plot scatterlot of risk status
# Plot heatmap of LASSO Cox gene signature
# Organise in Inkscape for manuscript.
#########################################################################
os_mat = as.data.frame(cox)
os_mat = os_mat[order(os_mat$risk_score, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$bcr_status = factor(os_mat$bcr_status)
os_mat$years = floor(os_mat$days_to_follow_up/365)

pdf("/data/github/pca_network/results/TCGA_DFS3/Training_scatter_risk.pdf", width = 8, height = 3)
ggscatter(os_mat, y="risk_score", x="patients_inc_risk", color="risk_category", fill="risk_category", 
          ylab="Risk Score", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = mean(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = mean(os_mat$risk_score), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/TCGA_DFS3/Training_scatter_days.pdf", width = 8, height = 4)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(bcr_status), x="patients_inc_risk", y="days_to_follow_up", shape="bcr_status", ylab = "Time (days)",
                  color="bcr_status", fill="bcr_status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = mean(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$risk_category)

col <- c("red", "royalblue3")
names(col) <- c("high", "low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/TCGA_DFS3/Training_lasso_genes_heatmap.pdf", height = 4, width = 8)
p = pheatmap::pheatmap(t(os_mat[,1:3]), labels_col = FALSE, color = col_palette, cluster_cols = F, 
                       scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
print(p)
dev.off()

# timeROC prediction using risk scores. 

scaled_risk_roc = roc(os_mat$bcr_status, os_mat$risk_score, ci=T)
pdf("/data/github/pca_network/results/TCGA_DFS3/Training_ROC_overall.pdf", width = 8, height = 8)
pROC::plot.roc(scaled_risk_roc, col="red", xlim = c(0.9,0.1), add = F)
legend("bottomright", legend = paste("AUC = ", round(scaled_risk_roc$auc, 3), "; 95% CI: 0.6689-08251"))
dev.off()

library(timeROC)
library(survivalROC)
os_mat$bcr_status = as.numeric(os_mat$bcr_status)
ROC.train <- timeROC(T=os_mat$days_to_follow_up,
                 delta=os_mat$bcr_status,
                 marker=os_mat$risk_score,
                 cause=2,weighting="marginal",
                 times=c(365, floor(365*3), floor(365*5), floor(365*10)),
                 iid=TRUE)
confint(ROC.train, level = 0.95)
pdf("/data/github/pca_network/results/TCGA_DFS3/Training_TimeROC.pdf")
plot(ROC.train, time=365)
plot(ROC.train, time=365*3)
plot(ROC.train, time=365*5)
dev.off()
nobs <- length(os_mat$risk_score)


roc <- survivalROC(Stime=os_mat$days_to_follow_up,
                   status=os_mat$bcr_status,
                   marker = os_mat$risk_score,
                   method = 'NNE',
                   predict.time = (5*365),
                   span = 0.25*nobs^(-0.20))
plot(roc$FP, roc$TP, type="l", lwd=3, col=alpha("royalblue",0.7), xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP"),
     ylab="TP")
legend("bottomright", legend = paste("AUC = ", round(scaled_risk_roc$auc, 3)))

# make plots for TEST
# make overall survival plot 
# train
sub_meta = subset(test_meta, select=c(days_to_follow_up, bcr_status))
colnames(sub_meta) = c("time", "status")
x = as.matrix(t(test_mat))
y = sub_meta
y = as.matrix(y)
cox = as.data.frame(x[,Active.Genes])
risk = rowSums(cox*Active.Coefficients)
cox$risk_score = risk
cox$risk_category = ifelse(cox$risk_score > median(cox$risk_score), "high", "low")
cox = cbind(cox, meta[,c("days_to_follow_up", "bcr_status")])
res = survfit(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)

pdf("/data/github/pca_network/results/TCGA_DFS3/Test_Surv_plot.pdf", height=8,width=8)
ggsurvplot(res,
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

#########################################################################
# Plot risk category per patient,
# Plot scatterlot of risk status
# Plot heatmap of LASSO Cox gene signature
# Organise in Inkscape for manuscript.
#########################################################################
os_mat = as.data.frame(cox)
os_mat = os_mat[order(os_mat$risk_score, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$bcr_status = factor(os_mat$bcr_status)
os_mat$years = floor(os_mat$days_to_follow_up/365)

pdf("/data/github/pca_network/results/TCGA_DFS3/Test_scatter_risk.pdf", width = 8, height = 3)
ggscatter(os_mat, y="risk_score", x="patients_inc_risk", color="risk_category", fill="risk_category", 
          ylab="Risk Score", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = mean(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = mean(os_mat$risk_score), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/TCGA_DFS3/Test_scatter_days.pdf", width = 8, height = 4)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(bcr_status), x="patients_inc_risk", y="days_to_follow_up", shape="bcr_status", ylab = "Time (days)",
                  color="bcr_status", fill="bcr_status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = mean(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$risk_category)

col <- c("red", "royalblue3")
names(col) <- c("high", "low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/TCGA_DFS3/Test_lasso_genes_heatmap.pdf", height = 4, width = 8)
p = pheatmap::pheatmap(t(os_mat[,1:3]), labels_col = FALSE, color = col_palette, cluster_cols = F, 
                       scale = "row", annotation_col = ann_col, annotation_colors = ann_clr)
print(p)
dev.off()

os_mat = cox
os_mat = os_mat[sample(rownames(os_mat), 240),]
os_mat$bcr_status = as.numeric(os_mat$bcr_status)
ROC.train <- timeROC(T=os_mat$days_to_follow_up,
                     delta=os_mat$bcr_status,
                     marker=os_mat$risk_score,
                     cause=1,weighting="marginal",
                     times=c(365, floor(365*3), floor(365*5), floor(365*10)),
                     iid=TRUE)
ROC.train
confint(ROC.train, level = 0.95)
pdf("/data/github/pca_network/results/TCGA_DFS3/Test_TimeROC.pdf")
plot(ROC.train, time=365)
plot(ROC.train, time=365*3)
plot(ROC.train, time=365*5)
dev.off()

scaled_risk_roc = roc(os_mat$bcr_status, os_mat$CTHRC1, ci=T)
pdf("/data/github/pca_network/results/TCGA_DFS3/Test_ROC_overall.pdf", width = 8, height = 8)
pROC::plot.roc(scaled_risk_roc, col="red", xlim = c(0.9,0.1), add = F)
legend("bottomright", legend = paste("AUC = ", round(scaled_risk_roc$auc, 3), "; 95% CI: 0.6689-08251"))
dev.off()
