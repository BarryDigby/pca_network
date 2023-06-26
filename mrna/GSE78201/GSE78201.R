#!/usr/bin/env Rscript 

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
# changes design cols, 
SDRF$group = c(rep("CTL", 4), rep("MODERATE", 4), rep("RESISTANT", 4))
SDRF$replicate = rep(c("4", "3", "2", "1"), 3)
SDRF$plot_names = c(rep("CTL", 4), rep("48HR", 4), rep("6MO", 4))

table(rownames(SDRF) == colnames(exp))

PCA <- prcomp(t(ex), scale = T)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Group = SDRF$plot_names,
                     Rep  = SDRF$replicate)
library(ggplot2)
pdf("/data/github/pca_network/mrna/GSE78201/GSE78201_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Group", palette = c("chartreuse", "purple", "skyblue2"),
                  title = "PCA plot log-transformed\nRMA normalized expression data\n [GSE78201]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


## LIMMA
status <- factor(SDRF$group)
replicate <- factor(SDRF$replicate)

design = model.matrix( ~ 0 + status + replicate ) 

# MODERATE == 48HR after treat MDV1000
colnames(design) = c("CTL", "MODERATE", "RESISTANT", "REP2", "REP3", "REP4")

contrast_matrix <- limma::makeContrasts(MODERATE_vs_CTL = MODERATE - CTL,
                                        RESISTANT_vs_CTL = RESISTANT - CTL,
                                        MODERATE_v_all = MODERATE - (CTL + RESISTANT)/2,
                                        RESISTANT_v_all = RESISTANT - (CTL + MODERATE)/2,
                                        levels = design )

probes = read.table("/data/github/pca_network/data/illumina_symbol_probe.txt", header=T, sep="\t")
ann <- function(df){
  df = merge(df, probes, by.x=0, by.y="Probe_Id")
  df = df[order(df$adj.P.Val),]
  df = df[!duplicated(df$Symbol),]
  names = colnames(df[2:length(colnames(df))])
  names = c(names[length(names)], names[-length(names)])
  df = subset(df, select=c(paste0(names)))
  return(df)
}

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(ex, design = design), contrast_matrix))
MOD_V_CTL = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="MODERATE_vs_CTL")
MOD_V_CTL = ann(MOD_V_CTL)
write.table(MOD_V_CTL, "/data/github/pca_network/mrna/GSE78201/48h_vs_control.txt", sep="\t", quote = F, row.names = FALSE)

RES_V_CTL = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="RESISTANT_vs_CTL")
RES_V_CTL = ann(RES_V_CTL)
write.table(RES_V_CTL, "/data/github/pca_network/mrna/GSE78201/6mo_vs_control.txt", sep="\t", quote = F, row.names = FALSE)

MOD_V_ALL = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="MODERATE_v_all")
MOD_V_ALL = ann(MOD_V_ALL)
write.table(MOD_V_ALL, "/data/github/pca_network/mrna/GSE78201/48h_vs_all.txt", sep="\t", quote = F, row.names = FALSE)

RES_V_ALL = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "BH", coef="RESISTANT_v_all")
RES_V_ALL = ann(RES_V_ALL)
write.table(RES_V_ALL, "/data/github/pca_network/mrna/GSE78201/6mo_vs_all.txt", sep="\t", quote = F, row.names = FALSE)


summary = function(df){
  z = table(sign(df$logFC))
  up = z[[1]]
  down = z[[2]]
  print(paste0("Up reg: ", up, "     Down reg: ", down))
}


summary(MOD_V_CTL)
summary(RES_V_CTL)
summary(MOD_V_ALL)
summary(RES_V_ALL)