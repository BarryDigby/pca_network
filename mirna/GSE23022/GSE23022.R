#!/usr/bin/env Rscript 

raw_data_dir <- "/data/github/pca_network/data/GSE23022"
sdrf_location <- file.path(raw_data_dir, "E-GEOD-23022.sdrf.txt")
SDRF <- read.delim(sdrf_location)

rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- Biobase::AnnotatedDataFrame(SDRF)

# automatically installs 'pd.mirna.1.0' for us
raw_data <- oligo::read.celfiles(filenames=file.path(raw_data_dir, SDRF@data$Array.Data.File),
                                 verbose=F,
                                 phenoData = SDRF)

# drop this first column after populating DF
meta <- data.frame(matrix(NA, nrow = as.numeric(nrow(SDRF@data))))
meta$sample <- SDRF@data$Hybridization.Name
meta$stage <- c(SDRF@data$Characteristics.tumor.stage.[1:20], rep("Normal", 20))
meta$patient <- gsub("\\D", "", SDRF@data$Comment..Sample_description..1)
meta$status <- sapply(strsplit(SDRF@data$Comment..Sample_description..1, " "), function(x) x[1])
meta <- meta[,2:ncol(meta)]
rownames(meta) <- SDRF@data$Array.Data.File

index <- match(colnames(raw_data), rownames(meta))
meta_order <- meta[index,]
meta <- meta_order

Biobase::pData(raw_data) <- meta
eset <- oligo::rma(raw_data, normalize = TRUE)
exp <- Biobase::exprs(eset)

table(rownames(meta) == colnames(exp))

PCA <- prcomp(t(exp), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Status = meta$status,
                     Stage  = meta$stage)
library(ggplot2)
pdf("/data/github/pca_network/mirna/GSE23022/GSE23022_PCA.pdf", height=4, width=5)
ggpubr::ggscatter(dataGG, x="PC1", y="PC2",
                  color = "Status", palette = c("royalblue", "orangered"),
                  title = "Log-transformed normalized expression data\n [GSE23022]",
                  xlab = paste0("PC1, VarExp: ", percentVar[1], "%"),
                  ylab = paste0("PC2, VarExp: ", percentVar[2], "%"),
                  ellipse = FALSE, star.plot = FALSE,
                  ggtheme = theme_bw()) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


# subset for Human miRNAs..?
#probes <- oligo::getProbeInfo(raw_data)
#probes <- probes[grep("^hsa-", probes$man_fsetid),]

# dont subset for this exp
#exp <- subset(exp, rownames(exp) %in% probes$man_fsetid)

## LIMMA
status <- factor(meta$status, levels=c("Tumor", "Normal"))
patient <- as.factor(meta$patient)

design = model.matrix( ~ 0 + status + patient ) 

contrast_matrix <- limma::makeContrasts(tumor_vs_normal = statusTumor-statusNormal, levels = design )

fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(eset, design = design), contrast_matrix))
tt <- limma::topTable(fit, number=Inf, p.value=0.05, adjust.method ="BH")

# sanity check for let-7a - upregulated 0.6 FC in tumor vs normal
#boxplot(eset@assayData$exprs["hsa-let-7a_st",] ~ status)

# pretty good agreement with results in paper, roughly 5 miRs are different. 
# They used graph pad prism, we used limma which I am more comfortable citing. 
hsa_tt <- tt[grep("^hsa", row.names(tt)),]
hsa_tt <- tibble::rownames_to_column(hsa_tt, "miRNA")
hsa_tt$miRNA <- gsub("_st", "", hsa_tt$miRNA)

## Use updated -3p/-5p notation instead of *
alias = read.csv("/data/github/pca_network/data/mirbase_aliases_tabbed.txt", sep="\t", header=F)

## Use updated -3p/-5p notation instead of *
update_id = function(df) {
  out_vec = c()
  for (index in 1:length(df$miRNA)) {
    ## stage miRNA from tt
    mir = df$miRNA[index]
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


ann_res = function(df){
  df$miRNA <- gsub("-star", "\\*", df$miRNA)
  df$new_id <- update_id(df)
  return(df)
}

hsa_tt = ann_res(hsa_tt)

write.table(hsa_tt, "/data/github/pca_network/mirna/GSE23022/tumor_vs_normal.txt", row.names = F, quote=F, sep="\t")

# export count matrix for heatmaps
# revert probes back to original naming convention: (undo prior code.. )
hsa_tt$miRNA = paste(hsa_tt$miRNA, "_st", sep="")
hsa_tt$miRNA = gsub("\\*", "-star", hsa_tt$miRNA)
# correctly gets original 25
t_v_n = exp[which(rownames(exp) %in% hsa_tt$miRNA),]
tumor_v_normal = subset(hsa_tt, select=c(miRNA, new_id))
t_v_n = merge(t_v_n, tumor_v_normal, by.x=0, by.y="miRNA")
t_v_n = t_v_n[,2:ncol(t_v_n)]
t_v_n = t_v_n[,c(ncol(t_v_n),1:(ncol(t_v_n)-1))]
colnames(t_v_n)[1] = "SystematicName"
write.table(t_v_n, "/data/github/pca_network/mirna/GSE23022/heatmap_counts.txt", sep="\t", quote=F, row.names = F)
dataGG = subset(dataGG, select=(Status))
write.table(dataGG, "/data/github/pca_network/mirna/GSE23022/heatmap_meta.txt", sep="\t", quote=F, row.names = F)

