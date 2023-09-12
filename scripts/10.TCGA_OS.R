#!/usr/bin/env Rscript 

library(dplyr)
library(tidyr)
library(survminer)
library(survival)
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
# Convert ENS to symbols and subset matrix by Active.Genes,
# order correctly (FPKM)
# FPKM matrix now ready for expression plots
#########################################################################
load("/data/github/pca_network/results/prognostic_model_os2.RData")
ensv109 = read.csv("/data/github/pca_network/data/ensembl_v109_proteinatlas.csv", sep="\t", header = F)
colnames(ensv109) = c("ensembl_gene_id", "biotype", "hgnc_symbol")
ensv109 = ensv109[which(ensv109$hgnc_symbol %in% Active.Genes),]
atlas_mat = atlas_mat[which(rownames(atlas_mat) %in% ensv109$ensembl_gene_id),]
atlas_mat = merge(atlas_mat, ensv109[,c(1,3)], by.x=0, by.y="ensembl_gene_id")
atlas_mat = tibble::column_to_rownames(atlas_mat, "hgnc_symbol")
atlas_mat = atlas_mat[,c(2:ncol(atlas_mat))]
atlas_mat = as.data.frame(t(atlas_mat))
atlas_mat = atlas_mat[rownames(atlas_meta),]

#########################################################################
# Use scaled centered logcpm STAR counts for CoxPH/Kaplan Meier
# coefficients came from this data - must be used on this again
#########################################################################
scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
colnames(scaled) = gsub("\\.", "-", colnames(scaled))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes = mrna_attributes[which(mrna_attributes$external_gene_name %in% Active.Genes),]
scaled = merge(scaled, mrna_attributes[,1:2], by.x=0, by.y="ensembl_gene_id_version")
scaled = tibble::column_to_rownames(scaled, "external_gene_name")
scaled = scaled[,c(2:ncol(scaled))]
scaled = scaled[,rownames(atlas_meta)]
scaled = as.data.frame(t(scaled))


#########################################################################
# Multiply gene * coefficients = risk score
# append time to event, OS status
#########################################################################
os_mat = scaled
risk = rowSums(Active.Coefficients*os_mat)
os_mat$risk_score = risk
os_mat$risk_category = ifelse(os_mat$risk_score > median(os_mat$risk_score), "high", "low")
os_mat = cbind(os_mat, subset(atlas_meta, select=c(days_to_follow_up, status)))

#########################################################################
# Kaplan-Meier OS TCGA
#########################################################################
surv_object <- Surv(os_mat$days_to_follow_up, os_mat$status)
cox = coxph(surv_object ~ risk_category, data=os_mat)
res = survfit(surv_object ~ risk_category, data=os_mat)
logrank = survdiff(surv_object ~ risk_category, data=os_mat)

pdf("/data/github/pca_network/results/TCGA_OS2/Risk_category_scaled_OS.pdf", height=8,width=8)
ggsurvplot(res,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=os_mat,
           xlab="Time (days)")
dev.off()

#########################################################################
# Scatter Histogram FPKM expr data
#########################################################################
scaled = scaled[,Active.Genes]
signf_os = signf[which(signf$gene %in% Active.Genes),]
for(i in 1:nrow(signf_os)){
    row = as.data.frame(signf_os[i,])
    gene = row$gene
    opt = row$best_cutoff

    mat = data.frame(atlas_mat[,which(colnames(atlas_mat)==gene)])
    colnames(mat) = gene
    rownames(mat) = rownames(atlas_meta)
    mat = merge(mat, atlas_meta, by.x=0, by.y="sample")
    mat$status = ifelse(mat$vital == "Dead", 1, 0)
    mat$high_low = ifelse(mat[,2] > opt, "high", "low")
    surv = survfit(Surv(days_to_follow_up, status) ~ high_low, data=mat)

    p <- mat %>% arrange(vital) %>% ggscatterhist(mat,  x=paste0(gene), y="years_to_follow_up", palette = c("royalblue3","red1"),
                                                     ylab = "Time after diagnosis (years)", xlab = "Expression level (FPKM)", fill = "vital",
                                                     color="vital", shape="vital", alpha = 0.9, ggtheme = theme_bw(), size = 2,
                                                     margin.params = list(fill="vital"), margin.plot = "density", legend = "top")
  
    p$sp <- p$sp + geom_vline(xintercept = opt, linetype = "dashed", color = "black")
    
    pdf(paste0("/data/github/pca_network/results/TCGA_OS2/",gene,"_fpkm_scatter.pdf"), height=5, width=8)
    print(p)
    dev.off()
}

#########################################################################
# Kaplan Meier OS LASSO Cox genes 
# FPKM
#########################################################################

for(i in 1:nrow(signf_os)){

  row = as.data.frame(signf_os[i,])
  gene = row$gene
  opt = row$best_cutoff
  
  mat = data.frame(atlas_mat[,which(colnames(atlas_mat)==gene)])
  colnames(mat) = gene
  rownames(mat) = rownames(atlas_meta)
  mat = merge(mat, atlas_meta, by.x=0, by.y="sample")
  mat$high_low = ifelse(mat[,2] > opt, "high", "low")
  surv = survfit(Surv(days_to_follow_up, status) ~ high_low, data=mat)

  p <- ggsurvplot(surv, pval = TRUE, conf.int = F, risk.table = F, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "none", # Specify median survival
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("red1", "royalblue3"),
                  data=mat, xlab="Time (days)")
  
  pdf(paste0("/data/github/pca_network/results/TCGA_OS2/",gene,"_fpkm_survival.pdf"), height=5, width=7)
  print(p)
  dev.off()
}

#########################################################################
# Kaplan Meier OS LASSO Cox genes 
# Scaled
#########################################################################

for(i in 1:nrow(signf_os)){
  
  row = as.data.frame(signf_os[i,])
  gene = row$gene
  opt = row$best_cutoff
  
  mat = data.frame(scaled[,which(colnames(scaled)==gene)])
  colnames(mat) = gene
  rownames(mat) = rownames(atlas_meta)
  mat = merge(mat, atlas_meta, by.x=0, by.y="sample")
  best_thresh = roc(mat$status, mat[,2])
  best_thresh = coords(best_thresh, x="best")
  mat$high_low = ifelse(mat[,2] > best_thresh[1,1], "High", "Low")
  surv = survfit(Surv(days_to_follow_up, status) ~ high_low, data=mat)
  
  p <- ggsurvplot(surv, pval = TRUE, conf.int = F, risk.table = F, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "none", # Specify median survival
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("red1", "royalblue3"),
                  data=mat, xlab="Time (days)")
  
  pdf(paste0("/data/github/pca_network/results/TCGA_OS2/",gene,"_scaled_survival.pdf"), height=5, width=7)
  print(p)
  dev.off()
}


#########################################################################
# ROC curve RISK ~ OS
#########################################################################
scaled_risk_roc = roc(os_mat$status, os_mat$risk_score)
pdf("/data/github/pca_network/results/TCGA_OS2/ROC_scaled_risk_score.pdf", width = 8, height = 8)
pROC::plot.roc(scaled_risk_roc, col="red", xlim = c(0.9,0.1), add = F)
legend("bottomright", legend = paste("AUC = ", round(scaled_risk_roc$auc, 3)))
dev.off()

#########################################################################
# ROC curve 1,3,5,10 years RISK ~ OS
#########################################################################

os_mat$risk_category = as.numeric(os_mat$risk_category)
ROC.train <- timeROC(T=os_mat$days_to_follow_up,
                     delta=os_mat$status,
                     marker=os_mat$risk_category,
                     cause=1,weighting="marginal",
                     times=c(365, floor(365*3), floor(365*5), floor(365*10)),
                     iid=TRUE)
ROC.train
confint(ROC.train, level = 0.95)
pdf("/data/github/pca_network/results/TCGA_OS2/TimeROC.pdf")
plot(ROC.train, time=365)
plot(ROC.train, time=365*3)
plot(ROC.train, time=365*5)
plot(ROC.train, time=365*10)
dev.off()




nobs <- length(os_mat$risk_score)
roc <- survivalROC(Stime=os_mat$days_to_follow_up,
                   status=os_mat$status,
                   marker = os_mat$risk_score,
                   method = 'NNE',
                   predict.time = 365,
                   span = 0.25*nobs^(-0.20))

roc2 <- survivalROC(Stime=os_mat$days_to_follow_up,
                    status=os_mat$status,
                    marker = os_mat$risk_score,
                    method = 'NNE',
                    predict.time = 365*3,
                    span = 0.25*nobs^(-0.20))

roc3 <- survivalROC(Stime=os_mat$days_to_follow_up,
                    status=os_mat$status,
                    marker = os_mat$risk_score,
                    method = 'NNE',
                    predict.time = 365*5,
                    span = 0.25*nobs^(-0.20))

roc4 <- survivalROC(Stime=os_mat$days_to_follow_up,
                    status=os_mat$status,
                    marker = os_mat$risk_score,
                    method = 'NNE',
                    predict.time = 365*10,
                    span = 0.25*nobs^(-0.20)) 

pdf("/data/github/pca_network/results/TCGA_OS2/ROC_10year_scaled_status-risk_score.pdf", width=8, height = 8)
plot(roc$FP, roc$TP, type="l", lwd=3, col=alpha("royalblue",0.7), xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP"),
     ylab="TP")
lines(roc2$FP, roc2$TP, type = "S", col=alpha("red",0.7), lwd=3)
lines(roc3$FP, roc3$TP, type = "S", col=alpha("grey9",0.7), lwd=3)
lines(roc4$FP, roc$TP, type="S", col=alpha("green",0.7), lwd=3)
abline(a = 0, b = 1, lty = 2)
dev.off()


#########################################################################
# Plot risk category per patient,
# Plot scatterlot of risk status
# Plot heatmap of LASSO Cox gene signature
# Organise in Inkscape for manuscript.
#########################################################################
os_mat = os_mat[order(os_mat$risk_score, decreasing = F),]
os_mat$patients_inc_risk = seq(1,nrow(os_mat),1)
os_mat$status = factor(os_mat$status)
os_mat$years = atlas_meta$years_to_follow_up

pdf("/data/github/pca_network/results/TCGA_OS2/Scatter_risk_category.pdf", width = 8, height = 4)
ggscatter(os_mat, y="risk_score", x="patients_inc_risk", color="risk_category", fill="risk_category", 
                  ylab="Risk Score", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw(), size = 1) + 
  geom_vline(xintercept = mean(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" ) + 
  geom_hline(yintercept = mean(os_mat$risk_score), linetype="dashed", color="grey10")
dev.off()

pdf("/data/github/pca_network/results/TCGA_OS2/Scatter_risk_category_years.pdf", width = 8, height = 4)
ggpubr::ggscatter(os_mat %>% dplyr::arrange(status), x="patients_inc_risk", y="years", shape="status", ylab = "Time (years)",
                  color="status", fill="status", palette = c("royalblue3","red"), ggtheme = theme_bw()) + 
  geom_vline(xintercept = mean(os_mat$patients_inc_risk), linetype = "dashed", color = "grey10" )
dev.off()

col_palette <- colorRampPalette(colors = c("royalblue3", "white", "red"))(100)
ann_col = data.frame(row.names = rownames(os_mat),
                     Group = os_mat$risk_category)

col <- c("red", "royalblue3")
names(col) <- c("high", "low")
ann_clr <- list(Group = col)
pdf("/data/github/pca_network/results/TCGA_OS2/lasso_genes_heatmap.pdf", height = 4, width = 8)
p = pheatmap::pheatmap(t(os_mat[,1:9]), labels_col = FALSE, color = col_palette, cluster_cols = F, 
                   scale = "column", annotation_col = ann_col, annotation_colors = ann_clr)
print(p)
dev.off()

