# How do the external PCa signatures perform on the unseen, external validation dataset?
library(tidyr)
library(dplyr)
library(survival)
library(RegParallel)
library(survminer)
library(timeROC)
library(ggplot2)
library(glmnet)
library(mfp)
library(CoxBoost)
library(boot)
library(rms)
library(snowfall)
library(survAUC)
library(meta)
library(metafor)
library(survcomp)
library(randomForestSRC)
library(mRMRe)

outdir <- "/data/github/pca_network/results/2025-07-01/"
#outdir <- paste(outdir, Sys.Date(), sep="")
if(!dir.exists(outdir)){
  dir.create(outdir)
}


signatures <- readxl::read_xlsx(("/data/github/pca_network/data/pcadb_signatures/PCaDB_Prognostic_Signatures_Gene_List.xlsx"))

## unseen, external dataset validation
## external val
library(GEOquery)
library(hgu133plus2.db)

GSE46602 <- getGEO("GSE46602")
meta <- GSE46602$GSE46602_series_matrix.txt.gz@phenoData@data
meta <- meta %>%
  dplyr::rename(
    `days_to_follow_up` = `bcr_free_time:ch1`, 
    `bcr_status` = `bcr:ch1`, 
    `gleason_score` = `gleason_grade:ch1`, 
    `preop_psa` = `preop_psa:ch1`) %>% 
  dplyr::select(days_to_follow_up, bcr_status, gleason_score, preop_psa)%>%
  mutate(across(everything(), ~ if_else(. == "NA", NA_character_, .)))%>%
  mutate(bcr_status = ifelse(bcr_status == "YES", 1, 0))%>%
  mutate(days_to_follow_up = 30.44 * as.numeric(days_to_follow_up))%>%
  mutate(across(everything(), as.numeric))

meta <- na.omit(meta)

ext_mat <- GSE46602$GSE46602_series_matrix.txt.gz@assayData$exprs
ext_mat <- as.data.frame(ext_mat[, rownames(meta)])

x <- mapIds(hgu133plus2.db, keys = rownames(ext_mat), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
x <- data.frame(probe = names(x), symbol = x, stringsAsFactors = FALSE)
#x <- x %>% filter(symbol %in% genes)

ext_mat <- ext_mat[x$probe, ]
ext_mat$Amean <- rowMeans(ext_mat)
ext_mat <- ext_mat %>%
  tibble::rownames_to_column("probe") %>%
  left_join(x, by = "probe")

ext_mat <- ext_mat %>%
  arrange(desc(Amean)) %>%
  distinct(symbol, .keep_all = TRUE)

ext_mat <- na.omit(ext_mat)
rownames(ext_mat) <- NULL

ext_mat <- ext_mat %>%
  dplyr::select(-Amean, -probe) %>%
  tibble::column_to_rownames("symbol")


ext_mat <- data.frame(t(ext_mat))
ext_mat <- data.frame(scale(ext_mat, scale = T, center = T))

ext_mat <- cbind(ext_mat, meta)

ext_sig_results <- data.frame()
set.seed(123)

for(author in unique(signatures$Signature)){
  print(author)
  if(author == "Wu"){
    next
  }
  sig_genes <- signatures$`HGNC Symbol`[which(signatures$Signature %in% author)]
  missing_genes <- setdiff(sig_genes, colnames(ext_mat))
  sig_genes <- intersect(sig_genes, colnames(ext_mat))

  # External, unseen dataset results
  fx <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(sig_genes, collapse = " + ")))
  ext_fit <- coxph(fx, data = ext_mat)

  # fit is coxph on all 9 datasets above
  ext_mat$lp <- predict(ext_fit, ext_mat, "lp")
  ext_mat$prob <- predict(ext_fit, ext_mat, "risk")


  fit_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, ext_mat)
  out_val <- c(round(as.numeric(coef(fit_ext)), 4), round(confint(fit_ext)[,1], 4), round(confint(fit_ext)[,2], 4), summary(fit_ext)$coef[,5])

  z_val <- stats::qnorm(1 - (1-0.95)/2)
  ext_expected_events <- mean(ext_mat$prob) * length(ext_mat$prob)
  
  KM_ext <- summary(survival::survfit(Surv(ext_mat$days_to_follow_up, ext_mat$bcr_status) ~ 1), times = max(ext_mat$days_to_follow_up))
  ext_log_OE_ratio <- log(KM_ext$n.event) - log(ext_expected_events)
  ext_log_OE_ratio_se <- sqrt(KM_ext$surv / KM_ext$n.event)
  ext_OE_ratio <- exp(ext_log_OE_ratio)
  ext_OE_ratio_se <- sqrt(KM_ext$n.event * (1 - KM_ext$n.event / ext_expected_events) / ext_expected_events)
  ext_OE_ratio_lower <- exp(ext_log_OE_ratio - (z_val * ext_log_OE_ratio_se))
  ext_OE_ratio_upper <- exp(ext_log_OE_ratio + (z_val * ext_log_OE_ratio_se))
  
  out_val <- c(out_val, round(ext_OE_ratio, 4), round(ext_OE_ratio_lower, 4), round(ext_OE_ratio_upper, 4))
  
  form <- formula(paste("Surv(days_to_follow_up, bcr_status) ~ ", paste(sig_genes, collapse = "+"), " + offset(lp)"))
  fit_ext <- coxph(formula = form, data = ext_mat)

  out_val <- c(out_val, round(2*(diff(fit_ext$loglik)), 2), round(1-pchisq(2*(diff(fit_ext$loglik)), 8), 3))
  
  ext_c <- survcomp::concordance.index(x = ext_mat$lp, surv.time = ext_mat$days_to_follow_up, surv.event = ext_mat$bcr_status)
  out_val <- c(out_val, round(ext_c$c.index, 4), round(ext_c$lower, 4), round(ext_c$upper, 4))
  
  ext_mat$group <- as.numeric(ext_mat$lp > median(ext_mat$lp))
  fit_ext <- survfit(Surv(days_to_follow_up, bcr_status) ~ group, data = ext_mat)
  ext_mat$group <- as.factor(ext_mat$group)
  ext_mat$group <- relevel(ext_mat$group, ref = "0")
  fit_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ group, data = ext_mat)
  
  out_val <- c(out_val, round(as.numeric(coef(fit_ext)), 4), round(as.numeric(summary(fit_ext)$coef[,3]), 4),  paste(round(confint(fit_ext)[,1], 2), round(confint(fit_ext)[,2], 2), sep = "-"), summary(fit_ext)$coef[,5])

  ROC.1 <- timeROC(T=ext_mat$days_to_follow_up,
                   delta=ext_mat$bcr_status,
                   marker=ext_mat$lp,
                   cause=1,weighting="marginal",
                   times=c(365),
                   iid=TRUE)
  
  ROC.2  <- timeROC(T=ext_mat$days_to_follow_up,
                    delta=ext_mat$bcr_status,
                    marker=ext_mat$lp,
                    cause=1,weighting="marginal",
                    times=c(floor(365*3)),
                    iid=TRUE)
  
  ROC.3 <- timeROC(T=ext_mat$days_to_follow_up,
                   delta=ext_mat$bcr_status,
                   marker=ext_mat$lp,
                   cause=1,weighting="marginal",
                   times=c(floor(365*5)),
                   iid=TRUE)

  out_val <- c(out_val,
               round(as.numeric(ROC.1$AUC[2]), 2), paste(confint(ROC.1)$CI_AUC[,1], confint(ROC.1)$CI_AUC[,2], sep = "-"),
               round(as.numeric(ROC.2$AUC[2]), 2), paste(confint(ROC.2)$CI_AUC[,1], confint(ROC.2)$CI_AUC[,2], sep = "-"),
               round(as.numeric(ROC.3$AUC[2]), 2), paste(confint(ROC.3)$CI_AUC[,1], confint(ROC.3)$CI_AUC[,2], sep = "-"))
  
  out_val <- c(out_val, length(sig_genes))
  out_val <- c(out_val, paste(names(coef(ext_fit)), collapse = ", "))
  out_val <- c(out_val, paste(missing_genes, collapse = ", "))
  
  ext_sig_results <- rbind(ext_sig_results, t(as.data.frame(out_val)))
  rownames(ext_sig_results)[nrow(ext_sig_results)] <- paste(author)

}

colnames(ext_sig_results) <- c(
  "Slope Coef",
  "Slope LCI",
  "Slope UCI",
  "Slope Pval",
  "Slope ITL Coef",
  "Slope ITL LCI",
  "Slope ITL UCI",
  "ChiSq",
  "ChiSq Pval",
  "Harrells C",
  "Harrells C LCI",
  "Harrells C UCI",
  "Hazard Ratio",
  "HR SE",
  "HR CI",
  "HR Pval",
  "AUC 1yr",
  "CI 1yr",
  "AUC 3yr",
  "CI 3yr",
  "AUC 5yr",
  "CI 5yr",
  "n genes",
  "genes",
  "missing genes"
)

ext_sig_results <- ext_sig_results%>%
  tibble::rownames_to_column("Dataset")%>%
  dplyr::select(Dataset, everything())

write.table(ext_sig_results, file = paste(outdir, "/PCaDB_ext_GSE46602.txt", sep = ""), sep = "\t", quote = F, row.names = F)

# forest plot stuff
two_ext <- read.table("/data/github/pca_network/results/2025-07-01/coxboost_external.txt", sep="\t", header = T, check.names = FALSE)
two_ext$Dataset <- "Digby"
ext_sig_results <- na.omit(ext_sig_results)
ext_sig_results <- ext_sig_results %>% dplyr::select(c(1:23,25))
ext_sig_results <- rbind(ext_sig_results, two_ext)
ext_sig_results <- ext_sig_results %>% arrange(Dataset)


forest_plot <- function(model, effect, LCI, UCI, rlab){
  print(paste0("Filtering for ", model))
  rownames(ext_sig_results) <- NULL
  df <- ext_sig_results %>% tibble::column_to_rownames(., "Dataset")
  df <- df %>% dplyr::select(one_of(effect, LCI, UCI)) %>% mutate_all(as.numeric)
  met <- metagen(TE = df[, effect], lower = df[, LCI], upper = df[, UCI], studlab = rownames(df))
  ref_line <- ifelse(rlab == "Calibration Slope", 1, ifelse(rlab == "Harrell's C", NA, 0))
  pad = 0.1
  if (rlab == "Calibration Slope") {
    xlim_vals <- c(floor(min(met$lower, na.rm = TRUE)),
                   ceiling(max(met$upper, na.rm = TRUE)))
  } else if (rlab == "Calibration ITL") {
    xlim_vals <- c(-0.2,
                   max(met$upper, na.rm = TRUE) + pad)
  } else if (rlab == "Harrell's C") {
    xlim_vals <- c(max(0, min(met$lower, na.rm = TRUE) - pad), 1)
  }
  return(meta::forest(met, common = FALSE,rightcols=c("effect", "ci"), rightlabs=c(rlab, "95% CI"), leftcols=c("studlab"), leftlabs=c("Author"),addrows.below.overall=3, ref=ref_line, xlim=xlim_vals, text.predict=FALSE, overall.hetstat=TRUE, text.random=""))
}

pdf(file = paste(outdir, "/PCaDB_ext_slope.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Slope Coef", "Slope LCI", "Slope UCI", "Calibration Slope")
dev.off()

pdf(file = paste(outdir, "/PCaDB_ext_itl.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Slope ITL Coef", "Slope ITL LCI", "Slope ITL UCI", "Calibration ITL")
dev.off()

pdf(file = paste(outdir, "/PCaDB_ext_harrells.pdf", sep=""), width = 8, height = 6)
forest_plot("_val", "Harrells C", "Harrells C LCI", "Harrells C UCI", "Harrell's C")
dev.off()

source_file <- "/data/github/pca_network/scripts/10-PCa-sigs.R"
file.copy(source_file, outdir, overwrite = T)

# overlaps we have with others.
gene_list <- lapply(ext_sig_results$genes, function(g) trimws(unlist(strsplit(g, ","))))

# Name each entry in the list by the corresponding Dataset
names(gene_list) <- ext_sig_results$Dataset

for(study in names(gene_list)){
  
  print(gene_list$Digby[gene_list$Digby %in% gene_list[[study]]])
  
}
