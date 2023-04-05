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

# subset for Clone1
ex <- as.data.frame(ex[,c(1:6)])

# Fold Change ratio
for(i in 1:nrow(ex)){
  
  row <- ex[i,]
  
  control <- as.numeric(row[,1:3])
  clone1 <- as.numeric(row[,4:6])
  
  t.test_clone1 <- t.test(clone1, control, paired=T)
  
  clone1_p.val <- t.test_clone1$p.value
  
  clone1_ratio <- 2^((mean(clone1))-(mean(control)))
  
  ex$Clone1_AvExp[i] <- mean(clone1)
  ex$Control_AvExp[i] <- mean(control)
  
  ex$Clone1_FC[i] <- fold_change(clone1_ratio)
  ex$Clone1_pvalue[i] <- clone1_p.val 
}

# FDR
ex$Clone1_FDR <- p.adjust(ex$Clone1_pvalue, method = "fdr")

key  <- intersect(rownames(ex)[which(ex$Clone1_FC>=fc)], rownames(ex)[which(ex$Clone1_pvalue<=pval)])
key2 <- intersect(rownames(ex)[which(ex$Clone1_FC<=-fc)], rownames(ex)[which(ex$Clone1_pvalue<=pval)])
clone1_circs <- c(key,key2)
clone1 <- subset(ex, rownames(ex) %in% clone1_circs)
clone1$Probe_ID <- rownames(clone1)
clone1 <- subset(clone1, select=c(Probe_ID, Clone1_FC, Clone1_pvalue, Clone1_FDR))
write.table(clone1, "/data/github/pca_network/results/GSE118959_Clone1.txt", quote = F, sep="\t", row.names = F)


