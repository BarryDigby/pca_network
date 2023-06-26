#!usr/bin/env Rscript

library(DESeq2)

simple_files = list.files("/data/github/pca_network/data/GSE88752/")

gen_raw <- function(path, pattern){
  files = list.files(path, pattern, full.names=T, recursive=T, include.dirs=T)
  sample_names = sub("\\_.*", "", simple_files)
  mat = as.data.frame(do.call(cbind, lapply(files, function(x) data.table::fread(x, stringsAsFactors=FALSE))))
  ENSG = mat[,1]
  mat = mat[, seq(2, ncol(mat)+1, 2)]
  rownames(mat) = ENSG
  colnames(mat) = sample_names
  return(mat)
}

mat = gen_raw("/data/github/pca_network/data/GSE88752/", "\\.txt$")

# subset to LNCAP samples. construct meta data

samples = colnames(mat)
cell_line = c(rep("LAPC9", 10), rep("LNCaP", 12))
treatment = c(rep("LAPC9", 10), rep("LNCaP", 4), rep("LNCaP_Pri", 4), rep("LNCaP_Sec", 4))
replicates = c(rep(1:10), rep(1:4, 3))

meta = data.frame(samples = samples,
                  cell_line = cell_line,
                  treatment = treatment,
                  replicates = replicates)

mat <- mat[,c(meta$samples[which(meta$cell_line=="LNCaP")])]
meta <- meta[11:22,]

dds = DESeqDataSetFromMatrix(mat, meta, design = ~ 0 + treatment + replicates)
dds$treatment <- relevel(dds$treatment, ref="LNCaP")
dds <- DESeq(dds)

rld <- assay(rlog(dds))

PCA <- prcomp(t(rld))

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     treatment = meta$treatment)
library(ggplot2)
pdf("/data/github/pca_network/mrna/GSE88752/GSE88752_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "treatment", palette = c("chartreuse", "purple", "skyblue2"),
                  title = "PCA plot log-transformed\nRMA normalized expression data\n [GSE88752]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

mod_mat <- model.matrix(design(dds), colData(dds))

lncap <- colMeans(mod_mat[dds$treatment == "LNCaP", ])
primary <- colMeans(mod_mat[dds$treatment == "LNCaP_Pri", ])
secondary <- colMeans(mod_mat[dds$treatment == "LNCaP_Sec", ])

primary_v_control <- results(dds, filterFun=IHW::ihw, alpha=0.05, contrast = primary - lncap)
primary_v_control <- lfcShrink(dds, primary_v_control, contrast = primary - lncap, type="ashr")
primary_v_control <- as.data.frame(primary_v_control)
primary_v_control <- primary_v_control[primary_v_control$padj <= 0.05,]

# sanity check these vs no intercept.. 

# really suspicious LFC values, re-do DDS using intercept to conduct each contrast. instead of A V B V C, combine a and B for secondary vs. primary + control .