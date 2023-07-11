#!/usr/bin/env Rscript
library(limma)
counts = read.table("/data/github/pca_network/mrna/LNCaP/mRNA_counts.txt", header=T, sep="\t", row.names = "gene_symbols")
meta = data.frame(row.names = colnames(counts),
                  group = factor(c(rep("Clone_1", 3), rep("Clone_9", 3), rep("Control", 3))),
                  replicate = factor(rep(seq(1:3), 3)))

design = model.matrix( ~ 0 + group + replicate, data = meta )                 

y <- edgeR::DGEList(counts)
keep <- edgeR::filterByExpr(y, design)
y <- y[keep, ]

# normalize and run voom transformation
y <- edgeR::calcNormFactors(y)
logcpm <- edgeR::cpm(y, normalized.lib.sizes = T, log=TRUE)
scaled <- t(scale(t(logcpm)))
v <- limma::voom(y, design, plot = F)

contrast <- limma::makeContrasts(
  Clone1_Control = groupClone_1 - groupControl,
  Clone9_Control = groupClone_9 - groupControl,
  levels = colnames(design))

fit = limma::lmFit(v, design = design)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design = design), contrast))

Clone1 = limma::topTable(fit, number=Inf, p.value = 0.05, adjust.method = "BH", coef="Clone1_Control")
Clone1 = tibble::rownames_to_column(Clone1, "Gene")
Clone9 = limma::topTable(fit, number=Inf, p.value = 0.05, adjust.method = "BH", coef="Clone9_Control")
Clone9 = tibble::rownames_to_column(Clone9, "Gene")

write.table(Clone1, "/data/github/pca_network/mrna/LNCaP/clone1_control.txt", sep="\t", row.names = F, quote = F)
write.table(Clone9, "/data/github/pca_network/mrna/LNCaP/clone9_control.txt", sep="\t", row.names = F, quote = F)

PCA <- prcomp(t(scaled), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Group = meta$group)
library(ggplot2)
pdf("/data/github/pca_network/mrna/LNCaP/LNCaP_mrna_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Group", palette = c("forestgreen", "red", "royalblue"),
                  title = "Log-transformed normalized expression data\n [LNCaP RNA-Seq]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = F, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()