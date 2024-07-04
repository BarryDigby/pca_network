library(tidyr)
library(dplyr)
library(beanplot)
library(survival)

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
# Use logcpm STAR counts for CoxPH
# just use protein coding genes
#########################################################################
logcpm = read.csv("/data/github/pca_network/results/TCGA_mrna_logcpm.txt", header=T, row.names = 1, sep="\t")
colnames(logcpm) = gsub("\\.", "-", colnames(logcpm))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
mrna_attributes <- subset(mrna_attributes, mrna_attributes$external_gene_name != "")
logcpm = subset(logcpm, rownames(logcpm) %in% mrna_attributes$ensembl_gene_id_version)
mrna_attributes <- subset(mrna_attributes, mrna_attributes$ensembl_gene_id_version %in% rownames(logcpm))
logcpm = merge(logcpm, mrna_attributes[,1:2], by.x=0, by.y="ensembl_gene_id_version")
logcpm = tibble::column_to_rownames(logcpm, "external_gene_name")
logcpm = logcpm[,c(2:ncol(logcpm))]
logcpm = logcpm[,atlas_meta$Row.names]

#########################################################################
# bind meta data to expr dataset
#########################################################################
univ_mat = as.data.frame(t(scale(t(logcpm), scale = T, center = T)))
univ_mat = univ_mat[, which(colnames(univ_mat) %in% atlas_meta$sample)]
univ_mat = as.data.frame(t(univ_mat))
univ_mat = cbind(univ_mat, atlas_meta[,c("days_to_follow_up", "bcr_status")])
colnames(univ_mat) <- gsub("-", ".", colnames(univ_mat))
#########################################################################
# Extract my results first
#########################################################################

network_mat <- univ_mat[, c(genes, "days_to_follow_up", "bcr_status")]
model <- coxph(Surv(days_to_follow_up, bcr_status) ~ REG4 + SLC2A4 + JAG2 + CTHRC1,  network_mat)
# model summary
model_summary <- summary(model)
# Extract relevant info
pval <- log10(as.numeric(model_summary$logtest['pvalue']))
# stratify patients using linear predictors
risk_scores = predict(model, newdata = network_mat, type="lp")
cox = as.data.frame(network_mat[,c("REG4", "SLC2A4" ,"JAG2", "CTHRC1")])
cox$risk_score = risk_scores
cox = cbind(cox, atlas_meta[,c("days_to_follow_up", "bcr_status")])
q = quantile(risk_scores, c(.25, .5, .75), type=1)
cox$risk_category = ifelse(cox$risk_score >= median(cox$risk_score), "High Risk", "Low Risk")
## Hazard Ratio for High Risk Group
cox$risk_category <- relevel(factor(cox$risk_category), ref = "Low Risk")
fit_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
sum_fit <- summary(fit_hr)
hr <- sum_fit$coefficients[2]
hr_ci_low <- sum_fit$conf.int[3]
hr_ci_high <- sum_fit$conf.int[4]
hr_p <- sum_fit$coefficients[5]

result <- data.frame(Genes = "REG4 SLC2A4 JAG2 CTHRC1",
                     Logtest = pval,
                     HR = hr,
                     CI_low = hr_ci_low,
                     CI_high = hr_ci_high,
                     HR_Pval = hr_p)

#########################################################################
# logrank P; HR n times.. 
#########################################################################

set.seed(123)

for (i in seq(1,10000)) {
  
  # sample 4-gene set
  gene_sub_set <- sample(colnames(univ_mat)[1:12209], 4, replace = TRUE)
  
  # construct coxph model
  model_formula <- as.formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(gene_sub_set, collapse=" + ")))
  
  # run model
  model <- coxph(model_formula, univ_mat)
  
  # model summary
  model_summary <- summary(model)
  
  # Extract relevant info
  pval <- log10(as.numeric(model_summary$logtest['pvalue']))
  
  # stratify patients using linear predictors
  risk_scores = predict(model, newdata = univ_mat, type="lp")
  cox = as.data.frame(univ_mat[,gene_sub_set])
  cox$risk_score = risk_scores
  cox = cbind(cox, atlas_meta[,c("days_to_follow_up", "bcr_status")])
  q = quantile(risk_scores, c(.25, .5, .75), type=1)
  cox$risk_category = ifelse(cox$risk_score >= median(cox$risk_score), "High Risk", "Low Risk")
  
  ## Hazard Ratio for High Risk Group
  cox$risk_category <- relevel(factor(cox$risk_category), ref = "Low Risk")
  fit_hr <- coxph(Surv(days_to_follow_up, bcr_status) ~ risk_category, data=cox)
  sum_fit <- summary(fit_hr)
  
  hr <- sum_fit$coefficients[2]
  hr_ci_low <- sum_fit$conf.int[3]
  hr_ci_high <- sum_fit$conf.int[4]
  hr_p <- sum_fit$coefficients[5]
  
  out <- data.frame(Genes = paste(gene_sub_set, collapse=" "),
                    Logtest = pval,
                    HR = hr,
                    CI_low = hr_ci_low,
                    CI_high = hr_ci_high,
                    HR_Pval = hr_p)
  
  result <- rbind(result, out)
  
}

#pdf(
#  "/data/github/pca_network/results/bootstrap.pdf",
#  width = 6,
#  height = 4
#)
beanplot::beanplot(
  result$Logtest,
  range = 0,
  horizontal = TRUE,
  pch = ".",
  what = c(0, 1, 1, 0),
  xlab = expression("pvalue (log[10])"),
  lwd = 1.5,
  las = 1,
  cex.lab = 1.5,
  cex.axis = 1.0,
  col = rgb(220, 160, 20, maxColorValue = 255),
  border = 1
)

q <- quantile(result$Logtest, 0.05)
rect(
  -12,
  1 - 0.3,
  q,
  1 + 0.3,
  border = NA,
  col = rgb(0, 100, 0, 50, maxColorValue = 255)
)

points(result$Logtest[1],
       1,
       pch = 19,
       col = "red",
       cex = 2)
abline(v = log10(0.05),
       lwd = 3,
       col = "blue")
axis(
  side = 3,
  at = log10(0.05),
  labels = expression(p = log[10](0.05)),
  tick = T,
  col = "blue",
  cex.axis = 1.5,
  col.axis = "blue"
)
#dev.off()


