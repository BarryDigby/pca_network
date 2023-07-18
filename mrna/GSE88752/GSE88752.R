#!usr/bin/env Rscript

library(limma)
library(dplyr)

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
                  treatment = factor(treatment),
                  replicates = factor(replicates))

design = model.matrix( ~ 0 + treatment + replicates, data=meta)

y <- edgeR::DGEList(mat)

keep <- edgeR::filterByExpr(y, design)
y <- y[keep, ]

# normalize and run voom transformation
y <- edgeR::calcNormFactors(y)
logcpm <- edgeR::cpm(y, normalized.lib.sizes = T, log=TRUE)
scaled <- t(scale(t(logcpm)))
v <- limma::voom(y, design, plot = F)


PCA <- prcomp(t(scaled), scale. = F)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     treatment = meta$treatment)
library(ggplot2)
#pdf("/data/github/pca_network/mrna/GSE88752/GSE88752_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG[11:22,], x="PC1", y="PC2",
                  color = "treatment", palette = c("chartreuse", "purple", "skyblue2"),
                  title = "PCA plot log-transformed\nRMA normalized expression data\n [GSE88752]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#dev.off()

colnames(design)
contrast = limma::makeContrasts(
  primary_v_control = treatmentLNCaP_Pri - treatmentLNCaP,
  secondary_v_control = treatmentLNCaP_Sec - treatmentLNCaP,
  levels = colnames(design)
)

fit = limma::lmFit(v, design = design)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(v, design = design), contrast))

pri_v_ctl = limma::topTable(fit, number=Inf, p.value = 0.05, adjust.method = "BH", coef="primary_v_control")
pri_v_ctl = tibble::rownames_to_column(pri_v_ctl, "ensembl_gene_id_version")


# biomaRt not working here, or you were specifying incorrectly. 
# use gencode v21 GTF file to make map, select protein coding only. 
convert = read.table("/data/github/pca_network/mrna/GSE88752/ENSG_VERSION_TO_SYMBOL.txt", header=F, sep="\t")
convert = convert %>% distinct()
colnames(convert) = c("ensembl_gene_id_version", "biotype", "symbol")

pri_v_ctl = merge(pri_v_ctl, convert, by="ensembl_gene_id_version")
pri_v_ctl = pri_v_ctl[which(pri_v_ctl$biotype=="protein_coding"),]
write.table(pri_v_ctl, "/data/github/pca_network/mrna/GSE88752/primary_vs_control.txt", row.names = F, quote=F, sep="\t")

sec_v_ctl = limma::topTable(fit, number=Inf, p.value = 0.05, adjust.method = "BH", coef="secondary_v_control")
sec_v_ctl = tibble::rownames_to_column(sec_v_ctl, "ensembl_gene_id_version")

sec_v_ctl = merge(sec_v_ctl, convert, by="ensembl_gene_id_version")
sec_v_ctl = sec_v_ctl[which(sec_v_ctl$biotype=="protein_coding"),]
write.table(sec_v_ctl, "/data/github/pca_network/mrna/GSE88752/secondary_vs_control.txt", row.names = F, quote=F, sep="\t")
