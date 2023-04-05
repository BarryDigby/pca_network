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
gset <- getGEO("GSE118959", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21825", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- as.data.frame(exprs(gset))
rename <- c("Control_1", "Control_2", "Control_3",
            "Clone1_1", "Clone1_2", "Clone1_3",
            "Clone9_1", "Clone9_2", "Clone9_3")
colnames(ex) <- rename

# subset for Clone9
ex <- as.data.frame(cbind(ex[,1:3], ex[,7:9]))

# Fold Change ratio
for(i in 1:nrow(ex)){
  
  row <- ex[i,]
  
  control <- as.numeric(row[,1:3])
  clone9 <- as.numeric(row[,4:6])
  
  t.test_clone9 <- t.test(clone9, control, paired=T)
  
  clone9_p.val <- t.test_clone9$p.value
  
  clone9_ratio <- 2^((mean(clone9))-(mean(control)))
  
  ex$Clone9_AvExp[i] <- mean(clone9)
  ex$Control_AvExp[i] <- mean(control)
  
  ex$Clone9_FC[i] <- fold_change(clone9_ratio)
  ex$Clone9_pvalue[i] <- clone9_p.val 
}

# FDR
ex$Clone9_FDR <- p.adjust(ex$Clone9_pvalue, method = "fdr")

key  <- intersect(rownames(ex)[which(ex$Clone9_FC>=fc)], rownames(ex)[which(ex$Clone9_pvalue<=pval)])
key2 <- intersect(rownames(ex)[which(ex$Clone9_FC<=-fc)], rownames(ex)[which(ex$Clone9_pvalue<=pval)])
clone9_circs <- c(key,key2)
clone9 <- subset(ex, rownames(ex) %in% clone9_circs)
clone9$Probe_ID <- rownames(clone9)
clone9 <- subset(clone9, select=c(Probe_ID, Clone9_FC, Clone9_pvalue, Clone9_FDR))
write.table(clone9, "/data/github/pca_network/results/GSE118959_Clone9.txt", quote = F, sep="\t", row.names = F)


