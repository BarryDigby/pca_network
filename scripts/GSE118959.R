#!/usr/bin/env Rscript

library(GEOquery)

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
rename <-  c("Control_1", "Control_2", "Control_3",
             "Clone1_1", "Clone1_2", "Clone1_3",
             "Clone9_1", "Clone9_2", "Clone9_3")
colnames(ex) <- rename

# Fold Change ratio
for(i in 1:nrow(ex)){
  
  row <- ex[i,]
  
  control <- as.numeric(row[,1:3])
  clone1 <- as.numeric(row[,4:6])
  clone9 <- as.numeric(row[,7:9])
  
  t.test_clone1 <- t.test(clone1, control, paired=T)
  t.test_clone9 <- t.test(clone9, control, paired=T)
  
  clone1_p.val <- t.test_clone1$p.value
  clone9_p.val <- t.test_clone9$p.value
  
  clone1_ratio <- 2^((mean(clone1))-(mean(control)))
  clone9_ratio <- 2^((mean(clone9))-(mean(control)))
  
  ex$Control_AvExp[i] <- mean(control)
  ex$Clone1_AvExp[i] <- mean(clone1)
  ex$Clone9_AvExp[i] <- mean(clone9)
  
  ex$Clone1_FC[i] <- fold_change(clone1_ratio)
  ex$Clone1_pvalue[i] <- clone1_p.val 
  
  ex$Clone9_FC[i] <- fold_change(clone9_ratio)
  ex$Clone9_pvalue[i] <- clone9_p.val
}

# Calculate FDR
ex$Clone1_FDR <- p.adjust(ex$Clone1_pvalue, method = "fdr")
ex$Clone9_FDR <- p.adjust(ex$Clone9_pvalue, method = "fdr")

# Order results
probes <- read.csv("/data/github/pca_network/data/arraystar_probes.csv", header=T, sep="\t")

ex <- merge(ex, probes, by.x="row.names", by.y="probeID")

colnames(ex)[1] = "Probe_ID"

order <- c("Probe_ID", "circRNA", "Alias", "GeneSymbol",
           "Control_AvExp", "Clone1_AvExp", "Clone1_FC", "Clone1_pvalue", "Clone1_FDR",
           "Control_AvExp", "Clone9_AvExp", "Clone9_FC", "Clone9_pvalue", "Clone9_FDR",
           "Control_1", "Control_2", "Control_3",
           "Clone1_1", "Clone1_2", "Clone1_3",
           "Clone9_1", "Clone9_2", "Clone9_3")

ex <- ex[,order]

ex[ex==""] <- NA

write.table(ex, "/data/github/pca_network/results/GSE118959.txt", quote = F, sep="\t", row.names = F)


