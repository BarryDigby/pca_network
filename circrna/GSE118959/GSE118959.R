library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE118959", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21825", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

meta <- data.frame(Treatment = c(rep("Control", 3), rep("Clone1", 3), rep("Clone9", 3)),
                   Replicate = rep(c("1", "2", "3"), 3))
rownames(meta) = colnames(exprs(gset))
exp = exprs(gset)


PCA <- prcomp(t(exp), scale = T)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Treatment = meta$Treatment,
                     Replicate = meta$Replicate)
library(ggplot2)
pdf("/data/github/pca_network/circrna/GSE118959/GSE118959_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Treatment", palette = c("chartreuse4", "orangered", "navyblue"),
                  title = "Log-transformed normalized expression data\n[GSE118959]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
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

design <- model.matrix( ~ 0 + meta$Treatment )
colnames(design) <- c( "Clone1", "Clone9", "Control")
# design <- cbind(design, Control=meta$Treatment=="Control")
# design <- cbind(design, Clone1=meta$Treatment=="Clone1")
# design <- cbind(design, Clone9=meta$Treatment=="Clone9")

contrast = makeContrasts(Clone1_vs_Control = Clone1-Control,
                         Clone9_vs_Control = Clone9-Control,
                         Clone1_vs_Clone9  = Clone1-Clone9,
                         Clone1_vs_all     = Clone1-((Clone9+Control)/2), # correct ordering for expected FC :) 
                         Clone9_vs_all     = Clone9-((Clone1+Control)/2),
                         Control_vs_all    = Control-((Clone1+Clone9)/2),
                         levels = design)

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(exp, design = design), contrast))
clone1_vs_control = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH", coef = "Clone1_vs_Control")
clone1_vs_control = ann_res(clone1_vs_control)

clone9_vs_control = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH", coef = "Clone9_vs_Control")
clone9_vs_control = ann_res(clone9_vs_control)

clone1_vs_clone9 = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH", coef = "Clone1_vs_Clone9")
clone1_vs_clone9 = ann_res(clone1_vs_clone9)

clone1_vs_all = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH", coef = "Clone1_vs_all")
clone1_vs_all = ann_res(clone1_vs_all)

clone9_vs_all = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH", coef = "Clone9_vs_all")
clone9_vs_all = ann_res(clone9_vs_all)

control_vs_all = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH", coef = "Control_vs_all")
control_vs_all = ann_res(control_vs_all)

write.table(clone1_vs_control, "/data/github/pca_network/circrna/GSE118959/clone1_vs_control.txt", row.names = F, quote=F, sep="\t")
write.table(clone9_vs_control, "/data/github/pca_network/circrna/GSE118959/clone9_vs_control.txt", row.names = F, quote=F, sep="\t")
write.table(clone1_vs_clone9, "/data/github/pca_network/circrna/GSE118959/clone1_vs_clone9.txt", row.names = F, quote=F, sep="\t")
write.table(clone1_vs_all, "/data/github/pca_network/circrna/GSE118959/clone1_vs_all.txt", row.names = F, quote=F, sep="\t")
write.table(clone9_vs_all, "/data/github/pca_network/circrna/GSE118959/clone9_vs_all.txt", row.names = F, quote=F, sep="\t")
write.table(control_vs_all, "/data/github/pca_network/circrna/GSE118959/control_vs_all.txt", row.names = F, quote=F, sep="\t")
