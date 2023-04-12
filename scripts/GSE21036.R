#!/usr/bin/env Rscript 

# add arguments if and when directories are moved around
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Pass the path to the raw CEL files and the corresponding SDRF file", call.=FALSE)
}

raw_data_dir <- as.character(args[1])
sdrf_filename <- as.character(args[2])

#raw_data_dir <- "/data/github/pca_network/data/GSE21036/"
sdrf_location <- file.path(raw_data_dir, sdrf_filename)
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
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Status", palette = c("chartreuse", "purple", "orange"),
                  title = "PCA plot log-transformed RMA normalized expression data\n [GSE21036]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = TRUE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# looks good, the vcap sample was down by metastasis samples  

## LIMMA
status <- factor(meta$tumor_type, levels=c("Normal", "Tumor", "Metastasis"))

design = model.matrix( ~ 0 + status) 

contrast_matrix <- limma::makeContrasts(Tumor_vs_Normal = statusTumor-statusNormal, 
                                        Metastasis_vs_Normal = statusMetastasis-statusNormal,
                                        Metastasis_vs_Tumor = statusMetastasis-statusTumor,
                                        levels = design)

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(yfilt, design = design), contrast_matrix))

keep <- c("SystematicName", "logFC", "AveExpr", "t" ,"P.Value","adj.P.Val","B")

tumor_vs_normal <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "Tumor_vs_Normal", adjust.method = "BH")
o <- order(tumor_vs_normal$AveExpr, decreasing=TRUE)
dup <- duplicated(tumor_vs_normal$GeneName[o])
tumor_vs_normal <- tumor_vs_normal[o,][!dup,]
tumor_vs_normal <- tumor_vs_normal[,keep]
# sanity boxplot
# yfilt$genes[yfilt$genes$GeneName=="hsa-miR-96",] pick correspinding probe?UID from top table
#boxplot(yfilt$E[53,] ~ status) # have to dig and find hsa-mir-96 manually to get rowindex


metastasis_vs_normal <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "Metastasis_vs_Normal", adjust.method = "BH")
o <- order(metastasis_vs_normal$AveExpr, decreasing=TRUE)
dup <- duplicated(metastasis_vs_normal$GeneName[o])
metastasis_vs_normal <- metastasis_vs_normal[o,][!dup,]
metastasis_vs_normal <- metastasis_vs_normal[,keep]
# sanity boxpot doesnt look great?
#yfilt$genes[yfilt$genes$GeneName=="hsa-miR-663",]
#boxplot(yfilt$E[234,] ~ status)


metastasis_vs_tumor <- limma::topTable(fit, number=Inf, p.value=0.05, coef = "Metastasis_vs_Tumor", adjust.method = "BH")
o <- order(metastasis_vs_tumor$AveExpr, decreasing=TRUE)
dup <- duplicated(metastasis_vs_tumor$GeneName[o])
metastasis_vs_tumor <- metastasis_vs_tumor[o,][!dup,]
metastasis_vs_tumor <- metastasis_vs_tumor[,keep]
#boxplot(yfilt$E[234,]~ status)


write.table(tumor_vs_normal, "/data/github/pca_network/results/GSE21036_tumor_vs_normal.txt", row.names = T, quote=F, sep="\t")
write.table(metastasis_vs_normal, "/data/github/pca_network/results/GSE21036_metastasis_vs_normal.txt", row.names = T, quote=F, sep="\t")
write.table(metastasis_vs_tumor, "/data/github/pca_network/results/GSE21036_metastasis_vs_tumor.txt", row.names = T, quote=F, sep="\t")

