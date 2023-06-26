#!/usr/bin/env Rscript 

library(GEOquery)
library(limma)

# manually make meta. 
sample_id = as.character(paste0("GSM", rep(3616205:3616228)))
status <-  c("ADT_ENZ", "ADT", "CTL", "ENZ", "ADT_XRT", "ADT_ENZ_XRT", "ENZ_XRT", "XRT")
replicate <- as.character(rep(1:3, each=8))
files = c("GSM3616205_LNCaP-N1-ADT_ENZA_MTH004B01B01_01.CEL.gz",
          "GSM3616206_LNCaP-N1-ADT_MTH003C05C02_01.CEL.gz",
          "GSM3616207_LNCaP-N1-CTR_MTH001C04C01_01.CEL.gz",
          "GSM3616208_LNCaP-N1-ENZA_MTH002A02A02_01.CEL.gz",
          "GSM3616209_LNCaP-N1-RT_ADT_MTH007G02G02_01.CEL.gz",
          "GSM3616210_LNCaP-N1-RT_ENZA_ADT_MTH008A06A03_01.CEL.gz",
          "GSM3616211_LNCaP-N1-RT_ENZA_MTH006C03C03_01.CEL.gz",
          "GSM3616212_LNCaP-N1-RT_MTH005F06F03_01.CEL.gz",
          "GSM3616213_LNCaP-N2-ADT_ENZA_MTH020E03E03_01.CEL.gz",
          "GSM3616214_LNCaP-N2-ADT_MTH019F05F02_01.CEL.gz",
          "GSM3616215_LNCaP-N2-CTR_MTH017D06D03_01.CEL.gz",
          "GSM3616216_LNCaP-N2-ENZA_MTH018E01E01_01.CEL.gz",
          "GSM3616217_LNCaP-N2-RT_ADT_MTH023A04A01_01.CEL.gz",
          "GSM3616218_LNCaP-N2-RT_ENZA_ADT_MTH024A03A03_01.CEL.gz",
          "GSM3616219_LNCaP-N2-RT_ENZA_MTH022C02C02_01.CEL.gz",
          "GSM3616220_LNCaP-N2-RT_MTH021D05D02_01.CEL.gz",
          "GSM3616221_LNCaP-N3-ADT_ENZA_MTH036C06C03_01.CEL.gz",
          "GSM3616222_LNCaP-N3-ADT_MTH035B03B03_01.CEL.gz",
          "GSM3616223_LNCaP-N3-CTR_MTH033D04D01_01.CEL.gz",
          "GSM3616224_LNCaP-N3-ENZA_MTH034H04H01_01.CEL.gz",
          "GSM3616225_LNCaP-N3-RT_ADT_MTH039A01A01_01.CEL.gz",
          "GSM3616226_LNCaP-N3-RT_ENZA_ADT_MTH040H05H02_01.CEL.gz",
          "GSM3616227_LNCaP-N3-RT_ENZA_MTH038G03G03_01.CEL.gz",
          "GSM3616228_LNCaP-N3-RT_MTH037B02B02_01.CEL.gz")

meta = data.frame(files = files,
                  sample = sample_id,
                  status = status,
                  replicate = replicate)

rownames(meta) <- meta$files

# load CEL files
raw_data_dir = "/data/github/pca_network/data/GSE126881/"
# automatically installs 'pd.mirna.1.0' for us
raw_data <- oligo::read.celfiles(filenames=file.path(raw_data_dir, meta$files),
                                 verbose=F)

Biobase::pData(raw_data) <- meta
eset <- oligo::rma(raw_data, normalize = TRUE)
exp <- Biobase::exprs(eset)

if(table(colnames(exp)==meta$files)){
  colnames(exp) = meta$sample
}

PCA <- prcomp(t(exp), scale = T)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     status = meta$status,
                     rep = meta$replicate)
library(ggplot2)
pdf("/data/github/pca_network/mrna/GSE126881/GSE126881_PCA.pdf", height = 4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "status",
                  shape = "rep",
                  title = "PCA plot log-transformed\nRMA normalized expression data\n [GSE126881]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = F, star.plot = F,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


# exp groups in paper:
# (ENZA + XRT vs. XRT, ADT+ XRT vs. XRT, ENZA+ ADT + XRT vs. XRT) i
# and
# ENZA vs. control [CTR], ADT vs. CTR, ENZA + ADT vs. CTR

## LIMMA
design = model.matrix( ~ 0 + factor(meta$status) + factor(meta$replicate) ) 
colnames(design) <- c("ADT", "ADT_ENZ", "ADT_ENZ_XRT", "ADT_XRT", "CTL", "ENZ", "ENZ_XRT", "XRT", "Rep2", "Rep3")

contrast_matrix <- limma::makeContrasts(ENZ_XRT_XRT = ENZ_XRT - XRT, 
                                        ADT_XRT_XRT = ADT_XRT - XRT,
                                        ADT_ENZ_XRT_XRT = ADT_ENZ_XRT - XRT,
                                        ENZ_CTL = ENZ - CTL,
                                        ADT_CTL = ADT - CTL,
                                        ADT_ENZ_CTL = ADT_ENZ - CTL,
                                        levels = design)

probes <- read.table("~/Downloads/GPL24324-69524.txt", header=T, sep="\t")
probes <- probes[,c(1, ncol(probes))]
probes$NM_ID <- sub("\\//.*", "", probes$SPOT_ID.1)
probes <- probes[,c(1, ncol(probes))]

library(biomaRt)
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
info <- biomaRt::getBM(attributes=c("refseq_mrna",
                                    "hgnc_symbol"),
                       filters = c("refseq_mrna"),
                       values = probes$NM_ID,
                       mart = mart,
                       useCache=T)

x = merge(probes, info, by.y="refseq_mrna", by.x="NM_ID") 

probes = probes[,c(1,2)]
write.table("/data/github/pca_network/data/pd.clariom.s.ht_converted.txt", row.names = F, quote = F, sep="\t")

ann_res = function(df){
  x = merge(df, probes, by.x=0, by.y)
}

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(exp, design = design), contrast_matrix))

ENZ_XRT = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method = "none", coef="ENZ_XRT_XRT")

