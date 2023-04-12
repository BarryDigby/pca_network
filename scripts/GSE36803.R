#!/usr/bin/env Rscript 

# Benign and PCA tissue selection, and RNA/DNA extraction: We utilized tissues from clinically
# localized PCA patients who underwent radical prostatectomy as a primary therapy. The collection of
# samples from castration resistant metastatic PCA patients was previously described

# add arguments if and when directories are moved around
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Pass the path to the raw CEL files and the corresponding SDRF file", call.=FALSE)
}

raw_data_dir <- as.character(args[1])
sdrf_filename <- as.character(args[2])

raw_data_dir <- "/data/github/pca_network/data/GSE36803"
sdrf_location <- file.path(raw_data_dir, "E-GEOD-36803.sdrf.txt")
SDRF <- read.delim(sdrf_location)

# remove miR treated entries in SDRF
SDRF <- SDRF[1:42,]

rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- Biobase::AnnotatedDataFrame(SDRF)

# automatically installs 'pd.mirna.1.0' for us
raw_data <- oligo::read.celfiles(filenames=file.path(raw_data_dir, SDRF@data$Array.Data.File),
                                 verbose=F,
                                 phenoData = SDRF)

# drop this first column after populating DF
meta <- data.frame(matrix(NA, nrow = as.numeric(nrow(SDRF@data))))
meta$status <- SDRF@data$FactorValue..ORGANISM.PART.
meta$patient <- gsub("\\D", "", SDRF@data$Comment..Sample_source_name.)
meta$check <- SDRF@data$Comment..Sample_source_name.
meta <- meta[,2:3]
meta$status <- ifelse(meta$status == "Benign, prostate tissue", "Benign", "PCA")
rownames(meta) <- SDRF@data$Array.Data.File

index <- match(colnames(raw_data), rownames(meta))
meta_order <- meta[index,]
meta <- meta_order

Biobase::pData(raw_data) <- meta
eset <- oligo::rma(raw_data, normalize = TRUE)
exp <- Biobase::exprs(eset)

table(rownames(meta) == colnames(exp))

# PCA <- prcomp(t(exp), scale = FALSE)
# 
# percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
# sd_ratio <- sqrt(percentVar[2] / percentVar[1])
# 
# dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
#                      Status = meta$status)
# library(ggplot2)
# ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
#                   color = "Status", palette = c("chartreuse", "purple"),
#                   title = "PCA plot log-transformed RMA normalized expression data\n [GSE23022]",
#                   xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
#                   ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
#                   ellipse = FALSE, star.plot = FALSE,
#                   ggtheme = theme_bw()) +
#   theme(legend.position = "right") +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# subset for Human miRNAs..?
#probes <- oligo::getProbeInfo(raw_data)
#probes <- probes[grep("^hsa-", probes$man_fsetid),]

# dont subset for this exp
#exp <- subset(exp, rownames(exp) %in% probes$man_fsetid)

## LIMMA
status <- factor(meta$status, levels=c("PCA", "Benign"))
status <- stats::relevel(status, ref="Benign")
patient <- factor(meta$patient)

design = model.matrix( ~ 0 + status + patient ) 

contrast_matrix <- limma::makeContrasts(PCA_vs_Benign = statusPCA-statusBenign, levels = design )

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(eset, design = design), contrast_matrix))
tt <- limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH")

# sanity check - upregulated in tumor vs normal
#boxplot(eset@assayData$exprs["hsa-miR-153_st",] ~ status)

# pretty good agreement with results in paper, roughly 5 miRs are different. 
# They used graph pad prism, we used limma which I am more comfortable citing. 
hsa_tt <- tt[grep("^hsa", row.names(tt)),]
write.table(hsa_tt, "/data/github/pca_network/results/GSE36806_PCA_vs_Benign.txt", row.names = T, quote=F, sep="\t")


