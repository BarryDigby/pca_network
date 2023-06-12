#!/usr/bin/env Rscript 

# add arguments if and when directories are moved around
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Pass the path to the raw CEL files and the corresponding SDRF file", call.=FALSE)
}

raw_data_dir <- as.character(args[1])
sdrf_filename <- as.character(args[2])

# uses weird illumina bead chip probes, use GEO query instead to obtain log2 norm vals 
library(GEOquery)
gset <- getGEO("GSE78201", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

ex <- na.omit(ex) # eliminate rows with NAs

raw_data_dir = "/data/github/pca_network/data/GSE78201/"
sdrf_location <- file.path(raw_data_dir, "E-GEOD-78201.sdrf.txt")
SDRF <- read.delim(sdrf_location)

# we only want LNCaP 
SDRF <- SDRF[SDRF$Characteristics..cell.line.=="LNCAP",]
  
#exp <- read.table("/data/github/pca_network/data/GSE78201/GSE78201_series_matrix.txt", header=T, row.names = "ID_REF",sep="\t")
ex <- ex[,SDRF$Assay.Name]
rownames(SDRF) = SDRF$Assay.Name
SDRF$group = c(rep("CTL", 4), rep("48HR", 4), rep("RESISTANT", 4))
SDRF$replicate = rep(c("4", "3", "2", "1"), 3)

table(rownames(SDRF) == colnames(exp))

PCA <- prcomp(t(ex), scale = T)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Group = SDRF$group,
                     Rep  = SDRF$replicate)
library(ggplot2)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Group", palette = c("chartreuse", "purple", "skyblue2"),
                  title = "PCA plot log-transformed RMA normalized expression data\n [GSE78201]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# subset for Human miRNAs..?
#probes <- oligo::getProbeInfo(raw_data)
#probes <- probes[grep("^hsa-", probes$man_fsetid),]

# dont subset for this exp
#exp <- subset(exp, rownames(exp) %in% probes$man_fsetid)

## LIMMA
status <- factor(SDRF$group, levels = c("CTL", "48HR", "RESISTANT"))
replicate <- factor(SDRF$replicate, levels = as.character(rep(1:4)))

design = model.matrix( ~ 0 + status + replicate ) 

# MODERATE == 48HR after treat MDV1000
colnames(design) = c("CTL", "MODERATE", "RESISTANT", "REP2", "REP3", "REP4")

contrast_matrix <- limma::makeContrasts(MODERATE_vs_CTL = MODERATE - CTL,
                                        RESISTANT_vs_CTL = RESISTANT - CTL,
                                        MODERATE_v_all = MODERATE - (CTL + RESISTANT)/2,
                                        RESISTANT_v_all = RESISTANT - (CTL + MODERATE)/2,
                                        levels = design )

probes = read.table("/data/github/pca_network/data/illumina_symbol_probe.txt", header=T, sep="\t")

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(ex, design = design), contrast_matrix))
MOD_V_CTL = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="MODERATE_vs_CTL")
MOD_V_CTL = merge(MOD_V_CTL, probes, by.x=0, by.y="Probe_Id")

# sort by padj and remove duplicate genes. 
MOD_V_CTL = MOD_V_CTL[order(MOD_V_CTL$adj.P.Val),]
MOD_V_CTL = MOD_V_CTL[!duplicated(MOD_V_CTL$Symbol),]
