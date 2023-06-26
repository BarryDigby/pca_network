#!/usr/bin/env Rscript 



# must create meta from series matrix :/ 
samples_named <- c("LNCAP_Control_1","LNCAP_ENZ.7d_1" ,"LNCAP_ENZ.14d_1","LNCAP_ENZ.21d_1","LNCAP_Control_2","LNCAP_ENZ.7d_2","LNCAP_ENZ.14d_2","LNCAP_ENZ.21d_2","LNCAP_Control_3","LNCAP_ENZ.7d_3","LNCAP_ENZ.14d_3","LNCAP_ENZ.21d_3")
samples_gsm <- c("GSM4258714","GSM4258715","GSM4258716","GSM4258717", "GSM4258718" ,"GSM4258719" , "GSM4258720","GSM4258721","GSM4258722", "GSM4258723" ,"GSM4258724" ,"GSM4258725")
time <- rep(c("0", "7", "14", "21"), 3)
treatment <- rep(c("control", "enzalutamide", "enzalutamide", "enzalutamide"), 3)
replicate <- c("1", "1", "1", "1", "2", "2", "2", "2", "3", "3", "3", "3")
group <-  gsub("_[^_]*$", "", samples_named)
meta <- data.frame(samples = samples_named,
                          accession = samples_gsm,
                          time = time, 
                          treatment = treatment,
                          replicate = replicate,
                   group = group)


exp <- read.table("/data/github/pca_network/data/GSE143408/GSE143408_series_matrix.txt", header=T, sep="\t", row.names = "ID_REF")
meta$group <- factor(meta$group, levels = c("LNCAP_Control", "LNCAP_ENZ.7d", "LNCAP_ENZ.14d", "LNCAP_ENZ.21d"))
table(colnames(exp) == meta$accession)
rownames(meta) = colnames(exp)
PCA <- prcomp(t(exp), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Group = meta$group)
library(ggplot2)
pdf("/data/github/pca_network/mrna/GSE143408/GSE143408_PCA.pdf", width = 5, height = 4)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Group", palette = c("chartreuse", "purple", "orange", "dodgerblue3"),
                  title = "PCA plot log-transformed\nRMA normalized expression data\n [GSE143408]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()
# looks really good. 

## LIMMA
# control for replicates.
genes <- read.table("/data/github/pca_network/data/GSE143408/genes.txt", header=F, sep="\t")
colnames(genes) = c("probe", "gene", "biotype")

group <- factor(meta$group, levels = c("LNCAP_Control", "LNCAP_ENZ.7d", "LNCAP_ENZ.14d", "LNCAP_ENZ.21d"))
replicates = factor(meta$replicate, levels = c("1", "2", "3"))

## SEQUENTIAL PAIRWISE (quoted in paper)

design = model.matrix( ~ 0 + group + replicates )
colnames(design) = c("LNCAP_Control", "LNCAP_ENZ.7d", "LNCAP_ENZ.14d", "LNCAP_ENZ.21d", "Re2", "Rep3")
contrast_matrix <- limma::makeContrasts(Day7 = LNCAP_ENZ.7d-LNCAP_Control,
                                        Day14 = LNCAP_ENZ.14d-LNCAP_ENZ.7d,
                                        Day21 = LNCAP_ENZ.21d-LNCAP_ENZ.14d,
                                        Day14_2 = LNCAP_ENZ.14d-(LNCAP_Control+LNCAP_ENZ.7d)/2,
                                        Day21_2 = LNCAP_ENZ.21d-((LNCAP_Control+LNCAP_ENZ.7d+LNCAP_ENZ.14d)/3),
                                        Day14_3 = LNCAP_ENZ.14d-LNCAP_Control,
                                        Day21_3 = LNCAP_ENZ.21d-LNCAP_Control,
                                        levels = design)

ann.res = function(df){
  x = merge(df, genes,by.x=0, by.y="probe")
  x = x[!(x$gene==""),]
  o = order(x$AveExp, decreasing = T)
  dup = duplicated(x$gene[o])
  x = x[o,][!dup,]
  return(x)
}

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(exp, design=design), contrast_matrix))
day7 = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="Day7")
day7 = ann.res(day7)

day14 = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="Day14")
day14 = ann.res(day14)

day21 = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="Day21")
#day21 = ann.res(day21)

day14_vs_below = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="Day14_2")
day14_vs_below = ann.res(day14_vs_below)

day21_vs_below = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="Day21_2")
day21_vs_below = ann.res(day21_vs_below)

day14_vs_control = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="Day14_3")
day14_vs_control = ann.res(day14_vs_control)

day21_vs_control = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="Day21_3")
day21_vs_control = ann.res(day21_vs_control)

write.table(day7, "/data/github/pca_network/mrna/GSE143408/7_vs_0.txt", row.names = T, quote=F, sep="\t")
write.table(day14, "/data/github/pca_network/mrna/GSE143408/14_vs_7.txt", row.names = T, quote=F, sep="\t")
write.table(day21, "/data/github/pca_network/mrna/GSE143408/21_vs_14.txt", row.names = T, quote=F, sep="\t")
write.table(day14_vs_below, "/data/github/pca_network/mrna/GSE143408/14_vs_7+0.txt", row.names = T, quote=F, sep="\t")
write.table(day21_vs_below, "/data/github/pca_network/mrna/GSE143408/21_vs_14+7+0.txt", row.names = T, quote=F, sep="\t")
write.table(day14_vs_control, "/data/github/pca_network/mrna/GSE143408/14_vs_0.txt", row.names = T, quote=F, sep="\t")
write.table(day21_vs_control, "/data/github/pca_network/mrna/GSE143408/21_vs_0.txt", row.names = T, quote=F, sep="\t")

summary = function(df){
  z = table(sign(df$logFC))
  up = z[[1]]
  down = z[[2]]
  print(paste0("Up reg: ", up, "     Down reg: ", down))
}

summary(day7)
summary(day14)
summary(day21)
summary(day14_vs_below)
summary(day14_vs_control)
summary(day21_vs_below)
summary(day21_vs_control)