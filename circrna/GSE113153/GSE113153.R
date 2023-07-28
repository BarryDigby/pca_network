library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE113153", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21825", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

meta <- data.frame(Gleason = c(rep("High", 5), rep("Low", 5)),
                   Replicate = rep(c("1", "2", "3", "4", "5"), 2))
rownames(meta) = colnames(exprs(gset))
exp = exprs(gset)

PCA <- prcomp(t(exp), scale = F)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Gleason = meta$Gleason,
                     Replicate = meta$Replicate)
library(ggplot2)
pdf("/data/github/pca_network/circrna/GSE113153/GSE113153_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Gleason", palette = c("orangered", "royalblue"),
                  title = "Log-transformed normalized expression data\n[GSE113153]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = F, star.plot = F,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


probes = read.table("/data/github/pca_network/data/arraystar_updated_circbaseID.csv", sep="\t", header=T)

ann_res = function(res){
  probe = subset(probes, select=c(probeID, circRNA, circbaseID, GeneSymbol))
  x = merge(res, probe, by.x=0, by.y="probeID")
  return(x)
}

design <- model.matrix( ~ 0 + meta$Gleason  )
colnames(design) <- c( "High", "Low")

contrast = makeContrasts(high_vs_low = High-Low,
                         levels = design)

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(exp, design = design), contrast))

high_vs_low = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH", coef = "high_vs_low")
high_vs_low = ann_res(high_vs_low)
write.table(high_vs_low, "/data/github/pca_network/circrna/GSE113153/high_vs_low.txt", row.names = F, quote=F, sep="\t")
