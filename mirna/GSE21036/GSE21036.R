#!/usr/bin/env Rscript 

library(limma)

raw_data_dir <- "/data/github/pca_network/data/GSE21036/"
sdrf_location <- file.path(raw_data_dir, "E-GEOD-21036.sdrf.txt")
SDRF <- read.delim(sdrf_location)

rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- Biobase::AnnotatedDataFrame(SDRF)
meta <- SDRF@data

files <- list.files(raw_data_dir, "^G.*\\.txt$")

x <- limma::read.maimages(files, source = "agilent", path = raw_data_dir, green.only = T)
y <- limma::backgroundCorrect(x, method="normexp")
y <- limma::normalizeBetweenArrays(y, method="quantile")
Control <- y$genes$ControlType==1L
NoSymbol <- is.na(y$genes$GeneName)
yfilt <- y[!Control & !NoSymbol,]

meta$tumor_type <- meta$FactorValue..TUMOR.TYPE.
meta$tumor_type <- gsub("not specified", "Normal", meta$tumor_type)
meta$tumor_type <- gsub("Primary tumor", "Tumor", meta$tumor_type)

rownames(meta) <- meta$Hybridization.Name

index <- match(colnames(yfilt), rownames(meta))
meta_order <- meta[index,]
meta <- meta_order

table(rownames(meta) == colnames(yfilt))

# remove the 'vcap' cell line sample
vcap_meta <- meta$Characteristics.cell.line.=="Vcap"
vcap_yfilt <- colnames(yfilt)=="GSM526532"

yfilt <- yfilt[,!vcap_yfilt]
meta <- meta[!vcap_meta,]

# double check
table(rownames(meta) == colnames(yfilt))

PCA <- prcomp(t(yfilt$E), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Status = meta$tumor_type)
library(ggplot2)
pdf("/data/github/pca_network/mirna/GSE21036/GSE21036_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Status", palette = c("chartreuse", "purple", "orange"),
                  title = "Log-transformed normalized expression data\n [GSE21036]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = TRUE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()
# looks good, the vcap sample was down by metastasis samples  

# really suspicious gene annotation, force recovery:
y = yfilt$genes
y$key <- paste(y$Row, y$Col, y$GeneName, sep="_")
rownames(yfilt$E) = y$key

## LIMMA
status <- factor(meta$tumor_type, levels=c("Normal", "Tumor", "Metastasis"))

design = model.matrix( ~ 0 + status) 

contrast_matrix <- limma::makeContrasts(tumor_vs_normal = statusTumor-statusNormal, 
                                        metastasis_vs_normal = statusMetastasis-statusNormal,
                                        metastasis_vs_tumor = statusMetastasis-statusTumor,
                                        tumor_vs_all = statusTumor-(statusNormal+statusMetastasis)/2,
                                        metastasis_vs_all = statusMetastasis-(statusNormal+statusTumor)/2,
                                        normal_vs_all = statusNormal-(statusTumor+statusMetastasis)/2,
                                        levels = design)

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(yfilt, design = design), contrast_matrix))

ann_res = function(df){
  keep <- c("Row", "Col", "SystematicName", "logFC", "AveExpr", "t" ,"P.Value","adj.P.Val","B")
  o <- order(df$AveExpr, decreasing=TRUE)
  dup <- duplicated(df$GeneName[o])
  df <- df[o,][!dup,]
  df <- df[,keep]
  return(df)
}

tumor_vs_normal <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "tumor_vs_normal", adjust.method = "BH")
tumor_vs_normal = ann_res(tumor_vs_normal)

metastasis_vs_normal <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "metastasis_vs_normal", adjust.method = "BH")
metastasis_vs_normal = ann_res(metastasis_vs_normal)

metastasis_vs_tumor <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "metastasis_vs_tumor", adjust.method = "BH")
metastasis_vs_tumor = ann_res(metastasis_vs_tumor)

metastasis_vs_all <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "metastasis_vs_all", adjust.method = "BH")
metastasis_vs_all = ann_res(metastasis_vs_all)

tumor_vs_all <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "tumor_vs_all", adjust.method = "BH")
tumor_vs_all = ann_res(tumor_vs_all)

normal_vs_all <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "normal_vs_all", adjust.method = "BH")
normal_vs_all = ann_res(normal_vs_all)

# names for boxplot sanity checks which work as expected now

write.table(tumor_vs_normal[,c(3:9)], "/data/github/pca_network/mirna/GSE21036/tumor_vs_normal.txt", row.names = F, quote=F, sep="\t")
write.table(metastasis_vs_normal[,c(3:9)], "/data/github/pca_network/mirna/GSE21036/metastasis_vs_normal.txt", row.names = F, quote=F, sep="\t")
write.table(metastasis_vs_tumor[,c(3:9)], "/data/github/pca_network/mirna/GSE21036/metastasis_vs_tumor.txt", row.names = F, quote=F, sep="\t")
write.table(metastasis_vs_all[,c(3:9)], "/data/github/pca_network/mirna/GSE21036/metastasis_vs_all.txt", row.names = F, quote=F, sep="\t")
write.table(tumor_vs_all[,c(3:9)], "/data/github/pca_network/mirna/GSE21036/tumor_vs_all.txt", row.names = F, quote=F, sep="\t")
write.table(normal_vs_all[,c(3:9)], "/data/github/pca_network/mirna/GSE21036/normal_vs_all.txt", row.names = F, quote=F, sep="\t")
