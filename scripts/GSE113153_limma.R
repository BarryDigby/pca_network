library(GEOquery)
library(limma)
library(umap)

# first argument FC, second Pval
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Pass numerical values for Fold Change and Pvalue filtering as arguments.", call.=FALSE)
}

fc <- as.numeric(args[1])
pval <- as.numeric(args[2])

# load series and platform data from GEO

gset <- getGEO("GSE113153", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21825", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000011111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("High","Low"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

nall <- nrow(gset)
gset <- gset[complete.cases(exprs(gset)), ]

# calculate precision weights and show plot of mean-variance trend
v <- vooma(gset, design, plot=F)
# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, pval)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000)

key  <- intersect(rownames(tT)[which(tT$logFC>=fc)], rownames(tT)[which(tT$P.Value<=pval)])
key2 <- intersect(rownames(tT)[which(tT$logFC<=-fc)], rownames(tT)[which(tT$P.Value<=pval)])
high_circs <- c(key,key2)
high_tt <- subset(tT, rownames(tT) %in% high_circs)
write.table(high_tt, "/data/github/pca_network/results/GSE113153_limma.txt", sep="\t", row.names = F, quote=FALSE)
