#!/usr/bin/env Rscript 

# 1. PCA plot 
# 2. plot alignment rates etc. 

alignment = read.csv("/data/LNCaP_miRNA/Trim_Experiment/alignment_stats.txt", header=F, sep=" ")
colnames(alignment) = c("Sample", "Run", "Processed", "Aligned", "Failed", "Reported")

ggpubr::ggline(alignment, x="Run", y="Processed", color = "Sample", palette = c("red", "forestgreen", "royalblue", "purple", "orange", "black", "black", "black", "black", "black", "black"))

ggpubr::ggline(alignment, x="Run", y="Aligned", color = "Sample", palette = c("red", "forestgreen", "royalblue", "purple", "orange", "black", "black", "black", "black", "black", "black"))


mat = read.csv("/data/LNCaP_miRNA/Trim_Experiment/mirtop/mirna.tsv", sep="\t", header=T, row.names="miRNA")
# dump old Clone1 samples
mat = mat[,-c(2,5)]
meta = data.frame(row.names = colnames(mat),
                  group = factor(c(rep("Clone_1", 3), rep("Clone_9", 3), rep("Control", 3))),
                  replicates = factor(rep(seq(1:3), 3)),
                  batch = factor(c( rep("A", 3), rep("B", 6))) )

design = model.matrix( ~ 0 + group + replicates + batch, data = meta )                 

y <- edgeR::DGEList(mat)
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
Clone1 = tibble::rownames_to_column(Clone1, "miRNA")
Clone9 = limma::topTable(fit, number=Inf, p.value = 0.05, adjust.method = "BH", coef="Clone9_Control")
Clone9 = tibble::rownames_to_column(Clone9, "miRNA")

write.table(Clone1, "/data/github/pca_network/mirna/LNCaP/clone1_control.txt", sep="\t", row.names = F, quote = F)
write.table(Clone9, "/data/github/pca_network/mirna/LNCaP/clone9_control.txt", sep="\t", row.names = F, quote = F)

PCA <- prcomp(t(scaled), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Group = meta$group,
                     Batch = meta$batch)
library(ggplot2)
pdf("/data/github/pca_network/mirna/LNCaP/LNCaP_mirna_PCA.pdf", height=3, width=4)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2", shape = "Batch",
                  color = "Group", palette = c("forestgreen", "red", "royalblue"),
                  title = "Log-transformed normalized expression data\n [LNCaP miRNA-Seq]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = F, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()