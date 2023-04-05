#!/usr/bin/env Rscript

library(GEOquery)

# first argument FC, second Pval
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Pass numerical values for Fold Change and Pvalue filtering as arguments.", call.=FALSE)
}

fc <- as.numeric(args[1])
pval <- as.numeric(args[2])

# Arraystar FC ratio
fold_change <- function(ratio){
  FC <- ifelse(ratio < 1, -1/ratio, ratio)
  return(FC)
}

# Load data
gset <- getGEO("GSE113153", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21825", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- as.data.frame(exprs(gset))
rename <- c("H-PCa1", "H-PCa2", "H-PCa3", "H-PCa4", "H-PCa5",
            "L-PCa1", "L-PCa2", "L-PCa3", "L-PCa4", "L-PCa5")
colnames(ex) <- rename

# Fold Change ratio
for(i in 1:nrow(ex)){
  
  row <- ex[i,]
  
  high <- as.numeric(row[,1:5])
  low <- as.numeric(row[,6:10])
  
  t.test_high <- t.test(high, low, paired=T)
  
  high_p.val <- t.test_high$p.value
  
  high_ratio <- 2^((mean(high))-(mean(low)))
  
  ex$High_AvExp[i] <- mean(high)
  ex$Low_AvExp[i] <- mean(low)
  
  ex$High_FC[i] <- fold_change(high_ratio)
  ex$High_pvalue[i] <- high_p.val 
}

# FDR
ex$High_FDR <- p.adjust(ex$High_pvalue, method = "fdr")

key  <- intersect(rownames(ex)[which(ex$High_FC>=fc)], rownames(ex)[which(ex$High_pvalue<=pval)])
key2 <- intersect(rownames(ex)[which(ex$High_FC<=-fc)], rownames(ex)[which(ex$High_pvalue<=pval)])
high_circs <- c(key,key2)
high <- subset(ex, rownames(ex) %in% high_circs)
high$Probe_ID <- rownames(high)
high <- subset(high, select=c(Probe_ID, High_FC, High_pvalue, High_FDR))
write.table(high, "/data/github/pca_network/results/GSE113153.txt", quote = F, sep="\t", row.names = F)


