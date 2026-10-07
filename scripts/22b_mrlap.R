## 22b_mrlap.R  (post-AHG editor comment 1)
## MRlap (Mounier & Kutalik 2023): corrects IVW for sample overlap, weak-instrument
## bias and winner's curse simultaneously, using cross-trait LD score regression.
## Instruments fixed to the manuscript's 50 P<1e-5 SNPs (do_pruning = FALSE),
## MR_threshold = 1e-5 so the winner's-curse correction uses the actual selection threshold.
suppressPackageStartupMessages({ library(MRlap); library(data.table) })
iv  <- readRDS("data/perio_instruments_finngen2025_relaxed.rds")
exp <- readRDS("ext_data/mrlap_exposure_perio.rds")
hm3 <- fread("ext_data/eur_w_ld_chr/w_hm3.snplist")

ad <- fread("ext_data/bellenguez/35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
            select = c("hm_rsid","hm_chrom","hm_pos","hm_effect_allele","hm_other_allele","hm_beta","standard_error","n_cas","n_con"))
setnames(ad, c("rsid","chr","pos","alt","ref","beta","se","n_cas","n_con"))
ad <- ad[!is.na(beta) & !is.na(se) & rsid %in% c(hm3$SNP, iv$SNP)][!duplicated(rsid)]
ad[, N := n_cas + n_con][, c("n_cas","n_con") := NULL]   # per-SNP N (max 487,511 incl. proxy cases)
ad <- ad[!is.na(N)]
cat("Outcome rows:", nrow(ad), " instruments present:", sum(iv$SNP %in% ad$rsid), "\n")

run <- function(N_exp, label) {
  e <- copy(exp)[, N := N_exp]
  r <- MRlap(exposure = as.data.frame(e), exposure_name = "periodontitis",
             outcome = as.data.frame(ad), outcome_name = "AD",
             ld = "ext_data/eur_w_ld_chr", hm3 = "ext_data/eur_w_ld_chr/w_hm3.snplist",
             do_pruning = FALSE, user_SNPsToKeep = iv$SNP, MR_threshold = 1e-5)
  saveRDS(r, sprintf("results/R_AHG/MRlap_%s.rds", label))
  m <- r$MRcorrection; l <- r$LDSC
  data.table(label, N_exp, nIV = m$m_IVs,
             observed_b = m$observed_effect, observed_se = m$observed_effect_se, observed_p = m$observed_effect_p,
             corrected_b = m$corrected_effect, corrected_se = m$corrected_effect_se, corrected_p = m$corrected_effect_p,
             diff_p = m$p_difference,
             h2_exp = l$h2_exp, h2_out = l$h2_out, rg = l$rg, cross_trait_intercept = l$int_crosstrait,
             cross_trait_intercept_se = l$int_crosstrait_se, rg_se = NA)
}
res <- rbind(run(250000, "Nexp250k"), run(473681, "Nexp474k"))
res[, `:=`(OR_obs = exp(observed_b), OR_cor = exp(corrected_b),
           OR_cor_LCI = exp(corrected_b - 1.96 * corrected_se), OR_cor_UCI = exp(corrected_b + 1.96 * corrected_se),
           OR_cor_doubling = exp(corrected_b * log(2)))]
fwrite(res, "results/R_AHG/MRlap_summary.csv")
print(t(res))
