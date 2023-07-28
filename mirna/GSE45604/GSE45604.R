#!/usr/bin/env Rscript 

library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE45604", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL14613", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

meta = data.frame(row.names = colnames(ex),
                  status = factor(c(rep("PCa", 28), rep("Normal", 10), rep("PCa", 22))))

PCA <- prcomp(t(ex), scale = T)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Status = meta$status)

library(ggplot2)
pdf("/data/github/pca_network/mirna/GSE45604/GSE45604_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Status", palette = c("chartreuse4", "orangered"),
                  title = "Log-transformed normalized expression data\n[GSE45604]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

design <- model.matrix( ~ 0 + meta$status )
colnames(design) = c("control", "pca")

contrast = makeContrasts(t_v_n = pca-control,
                         levels = design)

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(ex, design = design), contrast))
t_v_n = limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH", coef = "t_v_n")                  


t_v_n = tibble::rownames_to_column(t_v_n, "mirna")
t_v_n = t_v_n[which(grepl("hsa-", t_v_n$mirna)),]

# annotate correctly.. 

alias = read.csv("/data/github/pca_network/data/mirbase_aliases_tabbed.txt", sep="\t", header=F)

t_v_n$mirna = gsub("_st", "", t_v_n$mirna)
t_v_n$mirna = gsub("_x", "", t_v_n$mirna)
t_v_n$mirna = gsub("hp_", "", t_v_n$mirna)
t_v_n$mirna = gsub("-star", "\\*", t_v_n$mirna)


## Use updated -3p/-5p notation instead of *
update_id = function(df) {
  out_vec = c()
  for (index in 1:length(df$mirna)) {
    ## stage miRNA from tt
    mir = df$mirna[index]
    ## if its end in * add \\* for str_detect REGEX
    if (grepl("\\*", mir)) {
      mir <- gsub("\\*", "\\\\*", mir)
    }
    # grab the corresponding row in Alias
    row = alias %>% filter_all(any_vars(str_detect(., paste0("\\b", gsub("\\*", "\\\\\\*", mir), "\\b"))))
    if (nrow(row) == 0) {
      # no match with alias, keep original label as it is valid.
      out_vec = c(out_vec, mir)
    } else{
      # Vectorise these.. so you can do n + 1
      if (nrow(row) == 1) {
        # match does not like escaping chars
        row_vec = unname(unlist(row[1, ]))
        idx = match(gsub("\\\\", "", mir), row_vec)
        idx = idx + 1
        out_st = row_vec[idx]
        # if it now contains a *, move to next 'latest version'. this is rare. 
        if(grepl("\\*", out_st)){out_st = row_vec[idx+1]}
        # sometimes the matching string is correct, but also matches alias in V3 (no adjacent text)
        if(out_st == "" || is.na(out_st)){out_st=gsub("\\\\", "", mir)}
        out_vec = c(out_vec, out_st)
      } else if (nrow(row)==2){
        row_vec =  c(unname(unlist(row[1, ])), unname(unlist(row[2, ])))
        idx = match(gsub("\\\\", "", mir), row_vec)
        idx = idx + 1
        out_st = row_vec[idx]
        # if it now contains a *, move to next 'latest version'. this is rare. 
        if(grepl("\\*", out_st)){out_st = row_vec[idx+1]}
        # sometimes the matching string is correct, but also matches alias in V3 (no adjacent text)
        if(out_st == "" || is.na(out_st)){out_st=gsub("\\\\", "", mir)}
        out_vec = c(out_vec, out_st)
      } else {
        # edeg cases that draw three rows
        row_vec =  c(unname(unlist(row[1, ])), unname(unlist(row[2, ])), unname(unlist(row[3, ])))
        idx = match(gsub("\\\\", "", mir), row_vec)
        idx = idx + 1
        out_st = row_vec[idx]
        # if it now contains a *, move to next 'latest version'. this is rare. 
        if(grepl("\\*", out_st)){out_st = row_vec[idx+1]}
        # sometimes the matching string is correct, but also matches alias in V3 (no adjacent text)
        if(out_st == "" || is.na(out_st)){out_st=gsub("\\\\", "", mir)}
        out_vec = c(out_vec, out_st)
      }
    }
  }
  return(out_vec)
}

library(dplyr)
library(stringr)
t_v_n$new_id = update_id(t_v_n)

t_v_n$mirna = t_v_n$new_id
colnames(t_v_n)[1] = "miRNA"
t_v_n = t_v_n[,c(1:7)]
write.table(t_v_n, "/data/github/pca_network/mirna/GSE45604/tumor_vs_normal.txt", quote=F, sep="\t", row.names = F)
