#!/usr/bin/env Rscript

# load human protein atlas FPKM
atlas_mat = read.csv("/data/github/pca_network/data/prad_rna_cancer_sample.tsv", header=T, sep="\t")
atlas_mat = atlas_mat[,c(1,2,4)]

library(tidyr)

atlas_mat <- atlas_mat %>%
  pivot_wider(names_from = Sample, values_from = FPKM, values_fill = 0)
atlas_mat = tibble::column_to_rownames(atlas_mat, "Gene")

# load tcga survival data (NA's supplemented by hand using human protein atlas)
atlas_meta = read.csv("/data/github/pca_network/data/tcga_updated_meta.csv", header=T, sep=",")
rownames(atlas_meta) = atlas_meta$sample
atlas_meta = atlas_meta[which(atlas_meta$sample %in% colnames(atlas_mat)),]

rem = !(atlas_meta$vital=="Not Reported")
atlas_meta = atlas_meta[rem,]
atlas_mat = atlas_mat[,rem]

# stage mrnas in network
network = read.csv("/data/github/pca_network/results/circ_mirna_mrna_network.txt", header=T, sep="\t")
genes = unique(network$mrna)

# protein atlas uses ensembl v109 in gene descriptions. this should explain discrepancies between my results and website results... 
ensv109 = read.csv("/data/github/pca_network/data/ensembl_v109_proteinatlas.csv", sep="\t", header = F)
colnames(ensv109) = c("ensembl_gene_id", "biotype", "hgnc_symbol")
ensv109 = ensv109[which(ensv109$biotype=="protein_coding"),]
ensv109 = ensv109[which(ensv109$hgnc_symbol %in% genes),]
gene_ids = ensv109[,c(3,1)]
# library(biomaRt)
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# gene_ids <- getBM(attributes=c("hgnc_symbol",
#                                "ensembl_gene_id"),
#                          filters = c("hgnc_symbol"),
#                          values = genes,
#                          mart = mart,
#                          useCache = T)

# loop over genes in network, find those significantly linked with prognosis of PCa.
library(survival)
library(survminer)
library(pROC)
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
      mat$status = ifelse(mat$vital == "Dead", 1, 0)
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

# use protein atlas pval cutoff - the cutoff results are concordant w atlas explorer
signf = result[which(result$pvalue < 0.01),]

# produce plots for each signif

for(i in 1:nrow(signf)){
  row = as.data.frame(signf[i,])
  gene = row$gene
  ens = row$ensembl
  opt = row$best_cutoff
  
  mat = data.frame(t(atlas_mat[which(rownames(atlas_mat)==ens),]))
  mat = merge(mat, atlas_meta, by.x=0, by.y="sample")
  mat$status = ifelse(mat$vital == "Dead", 1, 0)
  mat$hi_lo = ifelse(mat[,2] > opt, "high", "low")
  surv = survfit(Surv(days_to_follow_up, status) ~ hi_lo, data=mat)
  
  # scatterhist
  p <- mat %>% arrange(vital) %>% ggscatterhist(mat,  x=paste0(ens), y="years_to_follow_up", palette = c("royalblue3","red1"),
                                                   ylab = "Time after diagnosis (years)", xlab = "Expression level (FPKM)", fill = "vital",
                                                   color="vital", shape="vital", alpha = 0.9, ggtheme = theme_bw(), size = 2,
                                                   margin.params = list(fill="vital"), margin.plot = "density", legend = "top")
  
  p$sp <- p$sp + geom_vline(xintercept = opt, linetype = "dashed", color = "black")
  
  pdf(paste0("/data/github/pca_network/results/surv_plots/",gene,"_scatter.pdf"), height=3, width=6)
  print(p)
  dev.off()
  
  # survival
  p1 <- ggsurvplot(surv, pval = TRUE, conf.int = F, risk.table = F, # Add risk table
                   risk.table.col = "strata", # Change risk table color by groups
                   linetype = "strata", # Change line type by groups
                   surv.median.line = "none", # Specify median survival
                   ggtheme = theme_bw(), # Change ggplot2 theme
                   palette = c("red1", "royalblue3"),
                   data=mat, xlab="Time (days)")
  # tweak to widt = scat width afte rremoving rhs hist
  pdf(paste0("/data/github/pca_network/results/surv_plots/",gene,"_survival.pdf"), height=3, width=5.1)
  print(p1)
  dev.off()
}

prog_pca = read.csv("/data/github/pca_network/results/prognostic_prostate.tsv", header=T, sep="\t")
table(prog_pca$Gene %in% signf$gene)

# LASSO cox on signif genes 
library(glmnet)

# test using zFPKM not 'scale()' 
# library(zFPKM)
# library(SummarizedExperiment)
# granges = read.csv("/data/github/pca_network/data/ensembl_v109_proteinatlas_granges.csv", sep="\t", header=F)
# colnames(granges) = c("ensembl_id", "hgnc", "chr", "start", "end", "strand")
# granges = granges[,c(1, 3:6)]
# assay = merge(granges, atlas_mat, by.x="ensembl_id", by.y=0, all.Y=T)
# assay = tibble::column_to_rownames(assay, var="ensembl_id")
# se = SummarizedExperiment::makeSummarizedExperimentFromDataFrame(assay, )
# names(assays(se)) = "fpkm"
# zfpkm = zFPKM(se, assayName = "fpkm")

# test using counts from STAR (not FPKM protein ATLAS)
# Run TCGA_mrna.R to produce 'scaled' and mrna_attributes. 
scaled = read.csv("/data/github/pca_network/results/TCGA_mrna_scaled.txt", header=T, row.names = 1, sep="\t")
colnames(scaled) = gsub("\\.", "-", colnames(scaled))
mrna_attributes = read.csv("/data/github/pca_network/results/TCGA_mrna_attributes.txt", header=T, sep="\t")
version_to_id = merge(signf, mrna_attributes, by.x="ensembl", by.y="ensembl_gene_id")
sub_atlas = scaled[which(rownames(scaled) %in% version_to_id$ensembl_gene_id_version),]
sub_atlas = merge(sub_atlas, version_to_id[,5:6], by.x=0, by.y="ensembl_gene_id_version")
sub_atlas = tibble::column_to_rownames(sub_atlas, "external_gene_name")
sub_atlas = sub_atlas[,c(2:ncol(sub_atlas))]
sub_atlas = sub_atlas[,rownames(atlas_meta)]


# # subset mats with signf 'result' above SWITCH for zFPKM at start 
# sub_atlas = atlas_mat[which(rownames(atlas_mat) %in% signf$ensembl),]
# #sub_atlas = atlas_mat[which(rownames(atlas_mat) %in% signf$ensembl),]
# #sub_atlas = t(scale(t(sub_atlas)))
# sub_atlas = merge(sub_atlas, gene_ids, by.x=0, by.y="ensembl_gene_id")
# sub_atlas = tibble::column_to_rownames(sub_atlas, "hgnc_symbol")
# sub_atlas = sub_atlas[,c(2:ncol(sub_atlas))]

# Subset metadata for cox
sub_meta = atlas_meta[,c(8,7)]
sub_meta$vital = ifelse(sub_meta$vital=="Dead", 1,0)
sub_meta$vital = factor(sub_meta$vital)
colnames(sub_meta) = c("time", "status")

x = as.matrix(t(sub_atlas))
rn = rownames(sub_meta)
y = sub_meta
y$status = as.numeric(as.character(y$status))
y = as.matrix(y)

# set.seed(123)
# cv.fit <- cv.glmnet(x, y, family="cox", alpha=1, maxit = 1000, lambda = NULL, type.measure = "deviance")
# fit = glmnet(x, y, family = "cox", alpha=1, maxit = 1000, lambda=NULL)
# Coefficients <- coef(fit, s = cv.fit$lambda.min)
# Active.Index <- which(Coefficients != 0)
# Active.Coefficients <- Coefficients[Active.Index]
# Active.Genes <- Coefficients@Dimnames[[1]][Active.Index]

load("/data/github/pca_network/results/LASSO_cox.RData")

x = as.data.frame(x)
x = x[,colnames(x) %in% model_result$Active.Genes]
risk_scores <- rowSums(model_result$Active.Coefficients * x)
x = cbind(x, sub_meta)
x$risk_score = risk_scores
x$risk_score_cat = ifelse(x$risk_score > mean(x$risk_score), "high", "low")

surv_object <- Surv(x$time, as.numeric(x$status))

res = coxph(Surv(x$time, as.numeric(x$status)) ~ risk_score_cat, data=x)
res2 = survfit(surv_object ~ risk_score_cat, data=x)
logrank = survdiff(surv_object ~ risk_score_cat, data=x)

model_result = data.frame(Active.Coefficients = Active.Coefficients,
                          Active.Genes = Active.Genes)

pdf("/data/github/pca_network/results/surv_plots/risk_surv.pdf", height=8,width=8)
ggsurvplot(res2,
           pval = TRUE, conf.int = T,
           risk.table = T, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           #ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=x,xlab="Time (days)")
dev.off()


x$risk_score_cat_num = ifelse(x$risk_score_cat == "high", 1, 0)
library(precrec)
risk_roc= precrec::evalmod(scores = x$risk_score, labels = x$status)
plot(risk_roc)
risk_roc = roc(x$status, x$risk_score)
pROC::plot.roc(risk_roc, col="red", xlim = c(0.9,0.1), add = F)
legend("bottomright", legend = paste("AUC = ", round(risk_roc$auc, 3)))

plot_mt = x
plot_mt = plot_mt[order(plot_mt$risk_score, decreasing = F),]
plot_mt$patients_inc_risk = seq(1, nrow(plot_mt))
plot_mt$years = plot_mt$time/365
plot_mt$status = ifelse(plot_mt$status == 1, "Dead", "Alive")
ggpubr::ggscatter(plot_mt, y="risk_score", x="patients_inc_risk", color="risk_score_cat", fill="risk_score_cat", 
                  ylab="Risk Score", xlab = NULL, palette = c("red","royalblue3"), ggtheme = theme_bw())
ggpubr::ggscatter(plot_mt, x="patients_inc_risk", y="years", shape="status", color="status", fill="status", palette = c("royalblue3","red"))

ann_col = data.frame(row.names = rownames(plot_mt),
                     Group = plot_mt$risk_score_cat)

col <- c("red", "royalblue3")
names(col) <- c("high", "low")
ann_clr <- list(Group = col)
pheatmap::pheatmap(t(plot_mt[,1:9]), labels_col = FALSE, color = hcl.colors(100, "RdBu",rev=F), cluster_cols = T, 
                   annotation_col = ann_col, annotation_colors = ann_clr, scale = "column")

plot(fit$glmnet.fit, "lambda", label=T)
plot(fit)
coef(fit)

atlas_meta$status = ifelse(atlas_meta$vital=="Alive", 0, 1)
atlas_meta$risk = x$risk_score
mult.cox = coxph(Surv(days_to_follow_up, status) ~ ajcc_path_t, data=atlas_meta)

coxph.control(iter.max = 10)
mult.cox = survivalAnalysis::analyse_multivariate(data = atlas_meta, time_status = vars(days_to_follow_up, status), 
                                                  vars(age, pri_gleason, sec_gleason))
survivalAnalysis::forest_plot(mult.cox)

ggforest(mult.cox, data=atlas_meta)

















# Sanity checks 

FBXO45 = data.frame(t(atlas_mat[which(rownames(atlas_mat)=="ENSG00000174013"),]))
FBXO45 = merge(FBXO45, atlas_meta, by.x=0, by.y="sample")


FBXO45$status = ifelse(FBXO45$vital == "Dead", 1, 0)
library(survival)

library(pROC)
rocky = roc(FBXO45$status, FBXO45$ENSG00000174013)
plot(rocky)
thresh = coords(rocky, "best", best.method = "youden")
FBXO45$hi_lo = ifelse(FBXO45$ENSG00000174013 > thresh$threshold, "high", "low")
res = coxph(Surv(days_to_follow_up, status) ~ hi_lo, data=FBXO45)
res2 = survfit(Surv(days_to_follow_up, status) ~ hi_lo, data=FBXO45)

library(ggpubr)
library(dplyr)
plot = FBXO45 %>% arrange(vital) %>% ggscatterhist(FBXO45,  x="ENSG00000174013", y="years_to_follow_up", palette = c("royalblue3","red1"),
                     ylab = "Time after follow up (years)", xlab = "FBXO45 FPKM", fill = "vital",
                     color="vital", shape="vital", alpha = 0.9, ggtheme = theme_bw(), size = 2,
                     margin.params = list(fill="vital"), margin.plot = "density", legend = "top")
plot$sp <- plot$sp +
  geom_vline(xintercept = 4.68355, linetype = "dashed", color = "black") + 
  annotate("text", x=4.8, y=12, label="optimal cutoff", angle=90)

plot

library(survminer)
#x = survfit(res)
ggsurvplot(res2,
           pval = TRUE, conf.int = F,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red1", "royalblue3"),
           data=FBXO45)


