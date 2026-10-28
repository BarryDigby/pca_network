# ==============================================================================
# Benchmarking Clinical-Only vs. Molecular Signature vs. Combined Models
# Responding to Reviewer 2, Major Comment 6:
# "Comparisons between clinical-only models (e.g., PSA and Gleason score etc.) 
#  and combined clinical-plus-molecular models should be performed..."
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(timeROC)
  library(survcomp)
  library(dplyr)
  library(ggplot2)
})

outdir <- "/data/github/pca_network/results/reviewer2_clinical_comparison/"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 1. Load Derivation and External Validation Datasets with Precomputed LP
# ------------------------------------------------------------------------------
cat("Loading derivation and external validation datasets...\n")

dev_mat <- readRDS("/data/github/pca_network/results/dev_mat_with_lp.RDS")
ext_mat <- readRDS("/data/github/pca_network/results/ext_mat_with_lp.RDS")

cat("Derivation set: ", nrow(dev_mat), " patients (", sum(dev_mat$bcr_status == 1), " BCR events)\n", sep="")
cat("External validation (GSE46602): ", nrow(ext_mat), " patients (", sum(ext_mat$bcr_status == 1), " BCR events)\n\n", sep="")

# Helper function to compute time-dependent AUCs
calc_aucs <- function(time, event, marker) {
  r1 <- timeROC(T = time, delta = event, marker = marker, cause = 1, weighting = "marginal", times = 365, iid = TRUE)
  r3 <- timeROC(T = time, delta = event, marker = marker, cause = 1, weighting = "marginal", times = floor(365 * 3), iid = TRUE)
  r5 <- timeROC(T = time, delta = event, marker = marker, cause = 1, weighting = "marginal", times = floor(365 * 5), iid = TRUE)
  c(
    auc_1 = round(as.numeric(r1$AUC[2]), 3),
    auc_3 = round(as.numeric(r3$AUC[2]), 3),
    auc_5 = round(as.numeric(r5$AUC[2]), 3)
  )
}

# ------------------------------------------------------------------------------
# 2. Benchmarking in the Derivation Dataset (n = 1,406)
# ------------------------------------------------------------------------------
cat("====================================================================\n")
cat("PART A: Derivation Dataset Benchmarking (n = 1,406)\n")
cat("====================================================================\n")

# Model 1: Clinical only (Gleason + Preop PSA)
fit_clin_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ gleason_score + preop_psa, data = dev_mat)

# Model 2: Molecular 23-gene signature only (lp)
fit_gene_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, data = dev_mat)

# Model 3: Combined model (Clinical + 23-gene signature)
fit_comb_dev <- coxph(Surv(days_to_follow_up, bcr_status) ~ gleason_score + preop_psa + lp, data = dev_mat)

# Linear predictors
dev_mat$lp_clin <- predict(fit_clin_dev, dev_mat, type = "lp")
dev_mat$lp_gene <- dev_mat$lp
dev_mat$lp_comb <- predict(fit_comb_dev, dev_mat, type = "lp")

# Harrell's C-index
c_clin_dev <- survcomp::concordance.index(x = dev_mat$lp_clin, surv.time = dev_mat$days_to_follow_up, surv.event = dev_mat$bcr_status, method = "noether")
c_gene_dev <- survcomp::concordance.index(x = dev_mat$lp_gene, surv.time = dev_mat$days_to_follow_up, surv.event = dev_mat$bcr_status, method = "noether")
c_comb_dev <- survcomp::concordance.index(x = dev_mat$lp_comb, surv.time = dev_mat$days_to_follow_up, surv.event = dev_mat$bcr_status, method = "noether")

# Statistical comparison of C-indices
comp_comb_vs_clin_dev <- survcomp::cindex.comp(c_comb_dev, c_clin_dev)
comp_gene_vs_clin_dev <- survcomp::cindex.comp(c_gene_dev, c_clin_dev)

# Likelihood Ratio Test: Does adding lp to clinical model improve fit?
lrt_dev <- anova(fit_clin_dev, fit_comb_dev)
lrt_dev_chisq <- lrt_dev[["Chisq"]][2]
lrt_dev_p <- lrt_dev[["P(>|Chi|)"]][2]

# AUCs
auc_clin_dev <- calc_aucs(dev_mat$days_to_follow_up, dev_mat$bcr_status, dev_mat$lp_clin)
auc_gene_dev <- calc_aucs(dev_mat$days_to_follow_up, dev_mat$bcr_status, dev_mat$lp_gene)
auc_comb_dev <- calc_aucs(dev_mat$days_to_follow_up, dev_mat$bcr_status, dev_mat$lp_comb)

dev_table <- data.frame(
  Model = c("Clinical Only (Gleason + PSA)", "23-Gene Signature Only (lp)", "Combined (Clinical + 23-Gene lp)"),
  C_Index = c(
    sprintf("%.3f (%.3f-%.3f)", c_clin_dev$c.index, c_clin_dev$lower, c_clin_dev$upper),
    sprintf("%.3f (%.3f-%.3f)", c_gene_dev$c.index, c_gene_dev$lower, c_gene_dev$upper),
    sprintf("%.3f (%.3f-%.3f)", c_comb_dev$c.index, c_comb_dev$lower, c_comb_dev$upper)
  ),
  AUC_1yr = c(auc_clin_dev["auc_1"], auc_gene_dev["auc_1"], auc_comb_dev["auc_1"]),
  AUC_3yr = c(auc_clin_dev["auc_3"], auc_gene_dev["auc_3"], auc_comb_dev["auc_3"]),
  AUC_5yr = c(auc_clin_dev["auc_5"], auc_gene_dev["auc_5"], auc_comb_dev["auc_5"]),
  AIC = c(round(AIC(fit_clin_dev), 1), round(AIC(fit_gene_dev), 1), round(AIC(fit_comb_dev), 1)),
  BIC = c(round(BIC(fit_clin_dev), 1), round(BIC(fit_gene_dev), 1), round(BIC(fit_comb_dev), 1))
)

print(dev_table)

cat("\nLikelihood Ratio Test (Combined vs. Clinical Only):\n")
cat("Chi-square = ", round(lrt_dev_chisq, 2), ", df = 1, P = ", format.pval(lrt_dev_p, digits = 4), "\n", sep="")

cat("C-index Improvement (Combined vs. Clinical Only):\n")
cat("Delta C = +", round(c_comb_dev$c.index - c_clin_dev$c.index, 3), ", P = ", format.pval(comp_comb_vs_clin_dev$p.value, digits = 4), "\n\n", sep="")

cat("Multivariable Model Coefficients (Combined Model in Derivation Cohort):\n")
comb_coef_summary <- summary(fit_comb_dev)$coefficients
comb_coef_df <- data.frame(
  Feature = rownames(comb_coef_summary),
  coef = round(comb_coef_summary[, "coef"], 4),
  HR = round(exp(comb_coef_summary[, "coef"]), 3),
  LCI = round(exp(confint(fit_comb_dev)[, 1]), 3),
  UCI = round(exp(confint(fit_comb_dev)[, 2]), 3),
  P_val = format.pval(comb_coef_summary[, "Pr(>|z|)"], digits = 4)
)
print(comb_coef_df)


# ------------------------------------------------------------------------------
# 3. Benchmarking in the External Validation Dataset (GSE46602, n = 36)
# ------------------------------------------------------------------------------
cat("\n====================================================================\n")
cat("PART B: External Validation Benchmarking (GSE46602, n = 36)\n")
cat("====================================================================\n")

# Model 1: Clinical only in external validation
fit_clin_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ gleason_score + preop_psa, data = ext_mat)

# Model 2: Locked 23-gene signature in external validation
fit_gene_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ lp, data = ext_mat)

# Model 3: Combined model in external validation
fit_comb_ext <- coxph(Surv(days_to_follow_up, bcr_status) ~ gleason_score + preop_psa + lp, data = ext_mat)

ext_mat$lp_clin <- predict(fit_clin_ext, ext_mat, type = "lp")
ext_mat$lp_gene <- ext_mat$lp
ext_mat$lp_comb <- predict(fit_comb_ext, ext_mat, type = "lp")

# Harrell's C-index
c_clin_ext <- survcomp::concordance.index(x = ext_mat$lp_clin, surv.time = ext_mat$days_to_follow_up, surv.event = ext_mat$bcr_status, method = "noether")
c_gene_ext <- survcomp::concordance.index(x = ext_mat$lp_gene, surv.time = ext_mat$days_to_follow_up, surv.event = ext_mat$bcr_status, method = "noether")
c_comb_ext <- survcomp::concordance.index(x = ext_mat$lp_comb, surv.time = ext_mat$days_to_follow_up, surv.event = ext_mat$bcr_status, method = "noether")

comp_comb_vs_clin_ext <- survcomp::cindex.comp(c_comb_ext, c_clin_ext)
comp_gene_vs_clin_ext <- survcomp::cindex.comp(c_gene_ext, c_clin_ext)

# Likelihood Ratio Test in External Validation
lrt_ext <- anova(fit_clin_ext, fit_comb_ext)
lrt_ext_chisq <- lrt_ext[["Chisq"]][2]
lrt_ext_p <- lrt_ext[["P(>|Chi|)"]][2]

# AUCs
auc_clin_ext <- calc_aucs(ext_mat$days_to_follow_up, ext_mat$bcr_status, ext_mat$lp_clin)
auc_gene_ext <- calc_aucs(ext_mat$days_to_follow_up, ext_mat$bcr_status, ext_mat$lp_gene)
auc_comb_ext <- calc_aucs(ext_mat$days_to_follow_up, ext_mat$bcr_status, ext_mat$lp_comb)

ext_table <- data.frame(
  Model = c("Clinical Only (Gleason + PSA)", "23-Gene Signature Only (lp)", "Combined (Clinical + 23-Gene lp)"),
  C_Index = c(
    sprintf("%.3f (%.3f-%.3f)", c_clin_ext$c.index, c_clin_ext$lower, c_clin_ext$upper),
    sprintf("%.3f (%.3f-%.3f)", c_gene_ext$c.index, c_gene_ext$lower, c_gene_ext$upper),
    sprintf("%.3f (%.3f-%.3f)", c_comb_ext$c.index, c_comb_ext$lower, c_comb_ext$upper)
  ),
  AUC_1yr = c(auc_clin_ext["auc_1"], auc_gene_ext["auc_1"], auc_comb_ext["auc_1"]),
  AUC_3yr = c(auc_clin_ext["auc_3"], auc_gene_ext["auc_3"], auc_comb_ext["auc_3"]),
  AUC_5yr = c(auc_clin_ext["auc_5"], auc_gene_ext["auc_5"], auc_comb_ext["auc_5"]),
  AIC = c(round(AIC(fit_clin_ext), 1), round(AIC(fit_gene_ext), 1), round(AIC(fit_comb_ext), 1)),
  BIC = c(round(BIC(fit_clin_ext), 1), round(BIC(fit_gene_ext), 1), round(BIC(fit_comb_ext), 1))
)

print(ext_table)

cat("\nLikelihood Ratio Test in External Validation (Combined vs. Clinical Only):\n")
cat("Chi-square = ", round(lrt_ext_chisq, 2), ", df = 1, P = ", format.pval(lrt_ext_p, digits = 4), "\n", sep="")

cat("C-index Improvement in External Validation (Combined vs. Clinical Only):\n")
cat("Delta C = +", round(c_comb_ext$c.index - c_clin_ext$c.index, 3), ", P = ", format.pval(comp_comb_vs_clin_ext$p.value, digits = 4), "\n\n", sep="")

cat("Multivariable Model Coefficients (Combined Model in External Cohort):\n")
comb_ext_summary <- summary(fit_comb_ext)$coefficients
comb_ext_df <- data.frame(
  Feature = rownames(comb_ext_summary),
  coef = round(comb_ext_summary[, "coef"], 4),
  HR = round(exp(comb_ext_summary[, "coef"]), 3),
  LCI = round(exp(confint(fit_comb_ext)[, 1]), 3),
  UCI = round(exp(confint(fit_comb_ext)[, 2]), 3),
  P_val = format.pval(comb_ext_summary[, "Pr(>|z|)"], digits = 4)
)
print(comb_ext_df)


# ------------------------------------------------------------------------------
# 4. Internal-External Cross-Validation by Study: Clinical vs. Combined
# ------------------------------------------------------------------------------
cat("\n====================================================================\n")
cat("PART C: Internal-External Cross-Validation (Leave-One-Study-Out)\n")
cat("====================================================================\n")

studies <- unique(dev_mat$study)
iecv_comp_list <- list()

for (s in studies) {
  val_fold <- dev_mat %>% filter(study == s)
  dev_fold <- dev_mat %>% filter(study != s)
  
  # Fit models on derivation fold
  m_clin <- coxph(Surv(days_to_follow_up, bcr_status) ~ gleason_score + preop_psa, data = dev_fold)
  m_comb <- coxph(Surv(days_to_follow_up, bcr_status) ~ gleason_score + preop_psa + lp, data = dev_fold)
  
  # Predict on held-out validation study
  lp_clin_val <- predict(m_clin, val_fold, type = "lp")
  lp_comb_val <- predict(m_comb, val_fold, type = "lp")
  
  c_clin_val <- survcomp::concordance.index(x = lp_clin_val, surv.time = val_fold$days_to_follow_up, surv.event = val_fold$bcr_status)$c.index
  c_comb_val <- survcomp::concordance.index(x = lp_comb_val, surv.time = val_fold$days_to_follow_up, surv.event = val_fold$bcr_status)$c.index
  
  iecv_comp_list[[s]] <- data.frame(
    Study = s,
    Events_Total = sprintf("%d/%d", sum(val_fold$bcr_status == 1), nrow(val_fold)),
    C_Clinical = round(c_clin_val, 3),
    C_Combined = round(c_comb_val, 3),
    Delta_C = round(c_comb_val - c_clin_val, 3)
  )
}

iecv_comp_df <- do.call(rbind, iecv_comp_list)
rownames(iecv_comp_df) <- NULL
print(iecv_comp_df)

cat("\nAverage IECV C-index across 9 held-out cohorts:\n")
cat("Clinical Only:  ", round(mean(iecv_comp_df$C_Clinical), 3), "\n")
cat("Combined Model: ", round(mean(iecv_comp_df$C_Combined), 3), "\n")
cat("Mean Improvement (Delta C): +", round(mean(iecv_comp_df$Delta_C), 3), "\n")


# ------------------------------------------------------------------------------
# 5. Export Tables & Plot
# ------------------------------------------------------------------------------
write.table(dev_table, file = paste0(outdir, "derivation_benchmarking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ext_table, file = paste0(outdir, "external_benchmarking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(iecv_comp_df, file = paste0(outdir, "iecv_clinical_vs_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comb_coef_df, file = paste0(outdir, "derivation_multivariable_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comb_ext_df, file = paste0(outdir, "external_multivariable_coefficients.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Generate a comparison barplot of C-indices
plot_df <- data.frame(
  Cohort = factor(c("Derivation Cohort (n=1,406)", "Derivation Cohort (n=1,406)", "Derivation Cohort (n=1,406)",
                    "External Validation (n=36)", "External Validation (n=36)", "External Validation (n=36)"),
                  levels = c("Derivation Cohort (n=1,406)", "External Validation (n=36)")),
  Model = factor(rep(c("Clinical Only (Gleason + PSA)", "23-Gene Signature Only", "Combined (Clinical + Signature)"), 2),
                 levels = c("Clinical Only (Gleason + PSA)", "23-Gene Signature Only", "Combined (Clinical + Signature)")),
  C_Index = c(c_clin_dev$c.index, c_gene_dev$c.index, c_comb_dev$c.index,
              c_clin_ext$c.index, c_gene_ext$c.index, c_comb_ext$c.index),
  LCI = c(c_clin_dev$lower, c_gene_dev$lower, c_comb_dev$lower,
          c_clin_ext$lower, c_gene_ext$lower, c_comb_ext$lower),
  UCI = c(c_clin_dev$upper, c_gene_dev$upper, c_comb_dev$upper,
          c_clin_ext$upper, c_gene_ext$upper, c_comb_ext$upper)
)

p <- ggplot(plot_df, aes(x = Cohort, y = C_Index, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), position = position_dodge(0.8), width = 0.25) +
  scale_fill_manual(values = c("#6baed6", "#fc9272", "#74c476")) +
  ylim(0, 1.0) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  labs(title = "Incremental Value: Clinical vs. Molecular vs. Combined Models",
       y = "Harrell's C-index (95% CI)", x = "") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 1))

ggsave(paste0(outdir, "cindex_comparison.pdf"), plot = p, width = 9.5, height = 6)
ggsave(paste0(outdir, "cindex_comparison.png"), plot = p, width = 9.5, height = 6, dpi = 300)

cat("\nExecution complete! Plots and tables saved in:\n", outdir, "\n")
