#!/usr/bin/env Rscript

# load human protein atlas FPKM
atlas_mat = read.csv("/data/github/pca_network/data/prad_rna_cancer_sample.tsv", header=T, sep="\t")
atlas_mat = atlas_mat[,c(1,2,4)]

library(tidyr)
library(dplyr)

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


# loop over genes in network, find those significantly linked with prognosis of PCa (OS)
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

# use protein atlas pval cutoff of 0.01 - the cutoff results are concordant w atlas explorer
signf = result[which(result$pvalue < 0.01),]
signf_os = signf
# produce plots for each signif

# for(i in 1:nrow(signf)){
#   row = as.data.frame(signf[i,])
#   gene = row$gene
#   ens = row$ensembl
#   opt = row$best_cutoff
# 
#   mat = data.frame(t(atlas_mat[which(rownames(atlas_mat)==ens),]))
#   mat = merge(mat, atlas_meta, by.x=0, by.y="sample")
#   mat$status = ifelse(mat$vital == "Dead", 1, 0)
#   mat$hi_lo = ifelse(mat[,2] > opt, "high", "low")
#   surv = survfit(Surv(days_to_follow_up, status) ~ hi_lo, data=mat)
# 
#   # scatterhist
#   p <- mat %>% arrange(vital) %>% ggscatterhist(mat,  x=paste0(ens), y="years_to_follow_up", palette = c("royalblue3","red1"),
#                                                    ylab = "Time after diagnosis (years)", xlab = "Expression level (FPKM)", fill = "vital",
#                                                    color="vital", shape="vital", alpha = 0.9, ggtheme = theme_bw(), size = 2,
#                                                    margin.params = list(fill="vital"), margin.plot = "density", legend = "top")
# 
#   p$sp <- p$sp + geom_vline(xintercept = opt, linetype = "dashed", color = "black")
# 
#   pdf(paste0("/data/github/pca_network/results/tcga_os_surv_plots/",gene,"_scatter.pdf"), height=3, width=6)
#   print(p)
#   dev.off()
# 
#   # survival
#   p1 <- ggsurvplot(surv, pval = TRUE, conf.int = F, risk.table = F, # Add risk table
#                    risk.table.col = "strata", # Change risk table color by groups
#                    linetype = "strata", # Change line type by groups
#                    surv.median.line = "none", # Specify median survival
#                    ggtheme = theme_bw(), # Change ggplot2 theme
#                    palette = c("red1", "royalblue3"),
#                    data=mat, xlab="Time (days)")
#   # tweak to widt = scat width afte rremoving rhs hist
#   pdf(paste0("/data/github/pca_network/results/tcga_os_surv_plots/",gene,"_survival.pdf"), height=3, width=5.1)
#   print(p1)
#   dev.off()
# }



##
##
##
## DFS in TCGA
##
##
##


pca_prad = readRDS("~/Downloads/TCGA-PRAD_eSet.RDS")
pca_meta = pca_prad@phenoData@data
# Remove last character from atlas meta -- They are the same patients. 
table(substr(atlas_meta$Row.names,1,nchar(atlas_meta$Row.names)-1) %in% pca_meta$sample_id)
# add information to atlas meta.
atlas_meta$pca_prad_key = substr(atlas_meta$Row.names,1,nchar(atlas_meta$Row.names)-1)
atlas_meta = merge(atlas_meta, pca_meta, by.x="pca_prad_key", by.y="sample_id") 

#
## time to BCR is the same as time to follow up, use time to follow up ~ BCR status.
#

# loop over genes in network, find those significantly linked with prognosis of PCa (DFS)
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
    best_thresh = roc(mat$bcr_status, mat[,2])
    best_thresh = coords(best_thresh, x="best")
    mat$high_low = ifelse(mat[,2] > best_thresh[1,1], "High", "Low")
    res.cox = survfit(Surv(days_to_follow_up, bcr_status) ~ high_low, data=mat)
    res.cox.pval = surv_pvalue(res.cox)
    pval = res.cox.pval$pval
    if( pval < 0.05 ){
      row = data.frame(gene=name, ensembl = ens, best_cutoff=best_thresh[1,1], pvalue=pval)
      result = rbind(result,row)
    }
  }
}


signf = result[which(result$pvalue < 0.01),]
signf_dfs = signf

# We are only interested in those that are signif according to LASSO cox.
# Script out of order here 9 <-> 10. 
load("/data/github/pca_network/results/LASSO_cox.RData")

signf = signf[which(signf$gene %in% Active.Genes),]


# for(i in 1:nrow(signf)){
#   row = as.data.frame(signf[i,])
#   gene = row$gene
#   ens = row$ensembl
#   opt = row$best_cutoff
#   
#   mat = data.frame(t(atlas_mat[which(rownames(atlas_mat)==ens),]))
#   mat = merge(mat, atlas_meta, by.x=0, by.y="sample")
#   mat$hi_lo = ifelse(mat[,2] > opt, "high", "low")
#   surv = survfit(Surv(days_to_follow_up, bcr_status) ~ hi_lo, data=mat)
#   
#   # scatterhist
#   p <- mat %>% arrange(bcr_status) %>% ggscatterhist(mat,  x=paste0(ens), y="years_to_follow_up", palette = c("royalblue3","red1"),
#                                                 ylab = "Time after diagnosis (years)", xlab = "Expression level (FPKM)", fill = "vital",
#                                                 color="vital", shape="vital", alpha = 0.9, ggtheme = theme_bw(), size = 2,
#                                                 margin.params = list(fill="vital"), margin.plot = "density", legend = "top")
#   
#   p$sp <- p$sp + geom_vline(xintercept = opt, linetype = "dashed", color = "black")
#   
#   pdf(paste0("/data/github/pca_network/results/tcga_dfs_surv_plots/",gene,"_scatter.pdf"), height=3, width=6)
#   print(p)
#   dev.off()
#   
#   # survival
#   p1 <- ggsurvplot(surv, pval = TRUE, conf.int = F, risk.table = F, # Add risk table
#                    risk.table.col = "strata", # Change risk table color by groups
#                    linetype = "strata", # Change line type by groups
#                    surv.median.line = "none", # Specify median survival
#                    ggtheme = theme_bw(), # Change ggplot2 theme
#                    palette = c("red1", "royalblue3"),
#                    data=mat, xlab="Time (days)")
#   # tweak to widt = scat width afte rremoving rhs hist
#   pdf(paste0("/data/github/pca_network/results/tcga_dfs_surv_plots/",gene,"_survival.pdf"), height=3, width=5.1)
#   print(p1)
#   dev.off()
# }


# save objects for COX in next script. Has all of the metadata needed from two sources of TCGA for same patients. 
save(atlas_mat, atlas_meta, signf_os, signf_dfs, network, file="/data/github/pca_network/results/TCGA_survival.RData")


# Sanity checks - fold into loop when happy with output

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



