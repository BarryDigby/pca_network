#!/usr/bin/env Rscript 

raw_data_dir <- "/data/github/pca_network/data/GSE23022"
sdrf_location <- file.path(raw_data_dir, "E-GEOD-23022.sdrf.txt")
SDRF <- read.delim(sdrf_location)

rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- Biobase::AnnotatedDataFrame(SDRF)

# automatically installs 'pd.mirna.1.0' for us
raw_data <- oligo::read.celfiles(filenames=file.path(raw_data_dir, SDRF@data$Array.Data.File),
                                 verbose=F,
                                 phenoData = SDRF)

# drop this first column after populating DF
meta <- data.frame(matrix(NA, nrow = as.numeric(nrow(SDRF@data))))
meta$sample <- SDRF@data$Hybridization.Name
meta$stage <- c(SDRF@data$Characteristics.tumor.stage.[1:20], rep("Normal", 20))
meta$patient <- gsub("\\D", "", SDRF@data$Comment..Sample_description..1)
meta$status <- sapply(strsplit(SDRF@data$Comment..Sample_description..1, " "), function(x) x[1])
meta <- meta[,2:ncol(meta)]
rownames(meta) <- SDRF@data$Array.Data.File

index <- match(colnames(raw_data), rownames(meta))
meta_order <- meta[index,]
meta <- meta_order

Biobase::pData(raw_data) <- meta
eset <- oligo::rma(raw_data, normalize = TRUE)
exp <- Biobase::exprs(eset)

table(rownames(meta) == colnames(exp))

PCA <- prcomp(t(exp), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Status = meta$status,
                     Stage  = meta$stage)
library(ggplot2)
pdf("/data/github/pca_network/mirna/GSE23022/GSE23022_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Status", palette = c("royalblue", "orangered"),
                  title = "Log-transformed normalized expression data\n [GSE23022]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


# subset for Human miRNAs..?
#probes <- oligo::getProbeInfo(raw_data)
#probes <- probes[grep("^hsa-", probes$man_fsetid),]

# dont subset for this exp
#exp <- subset(exp, rownames(exp) %in% probes$man_fsetid)

## LIMMA
status <- factor(meta$status, levels=c("Tumor", "Normal"))
patient <- as.factor(meta$patient)

design = model.matrix( ~ 0 + status + patient ) 

contrast_matrix <- limma::makeContrasts(tumor_vs_normal = statusTumor-statusNormal, levels = design )

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(eset, design = design), contrast_matrix))
tt <- limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH")

# sanity check for let-7a - upregulated 0.6 FC in tumor vs normal
#boxplot(eset@assayData$exprs["hsa-let-7a_st",] ~ status)

# pretty good agreement with results in paper, roughly 5 miRs are different. 
# They used graph pad prism, we used limma which I am more comfortable citing. 
hsa_tt <- tt[grep("^hsa", row.names(tt)),]
hsa_tt <- tibble::rownames_to_column(hsa_tt, "miRNA")
write.table(hsa_tt, "/data/github/pca_network/mirna/GSE23022/tumor_vs_normal.txt", row.names = F, quote=F, sep="\t")



