## 10_primary_total_effect_relaxed.R
## PRIMARY analysis (per strategic decision: total effect is the headline
## finding, mediation reframed as exploratory/secondary). Uses the 50-SNP
## relaxed-threshold (P<1e-5) instrument from Salminen 2025 FinnGen
## periodontitis-diagnosis GWAS (09_new_instruments_finngen_relaxed.R),
## which gives the best available combination of exposure specificity and
## aggregate power (R2=0.002357, power=50.8% at OR=1.15 vs AD).
suppressPackageStartupMessages({
  library(TwoSampleMR); library(dplyr); library(data.table)
})

AD <- "ebi-a-GCST90027158"       # Bellenguez 2022 (primary)
AD_SENS <- "ieu-b-2"             # Kunkle 2019 IGAP (sensitivity)

perio_iv <- readRDS("data/perio_instruments_finngen2025_relaxed.rds")
perio_iv$ncontrol.exposure <- perio_iv$samplesize.exposure - perio_iv$ncase.exposure
perio_iv$units.exposure <- "log odds"
cat("Instruments: n =", nrow(perio_iv), "\n")

cat("\nExtracting AD (Bellenguez 2022) outcome data...\n")
out_ad <- extract_outcome_data(snps = perio_iv$SNP, outcomes = AD)
dat <- harmonise_data(perio_iv, out_ad, action = 2)
dat$units.outcome <- "log odds"
cat("SNPs retained after harmonisation:", nrow(dat), "\n")

cat("\n=== PRIMARY: Total effect, Periodontitis (Salminen 2025, P<1e-5, n=50) -> AD (Bellenguez) ===\n")
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
print(res)
fwrite(res, "results/PRIMARY_total_effect_perio_to_AD_relaxed.csv")

cat("\n=== Heterogeneity (Cochran's Q) ===\n")
het <- mr_heterogeneity(dat)
print(het)
fwrite(het, "results/PRIMARY_heterogeneity_relaxed.csv")

cat("\n=== Pleiotropy (MR-Egger intercept) ===\n")
ple <- mr_pleiotropy_test(dat)
print(ple)
fwrite(ple, "results/PRIMARY_pleiotropy_relaxed.csv")

cat("\n=== Steiger directionality test ===\n")
dir_test <- tryCatch(directionality_test(dat), error = function(e) { cat(conditionMessage(e),"\n"); NULL })
if (!is.null(dir_test)) { print(dir_test); fwrite(dir_test, "results/PRIMARY_steiger_relaxed.csv") }

cat("\n=== Leave-one-out ===\n")
loo <- mr_leaveoneout(dat)
print(loo[order(loo$p), c("SNP","b","se","p")])
fwrite(loo, "results/PRIMARY_loo_relaxed.csv")

cat("\n=== Sensitivity: AD outcome = Kunkle 2019 IGAP ===\n")
out_adk <- tryCatch(extract_outcome_data(snps = perio_iv$SNP, outcomes = AD_SENS), error=function(e) NULL)
if (!is.null(out_adk) && nrow(out_adk) > 0) {
  dat_k <- harmonise_data(perio_iv, out_adk, action = 2)
  res_k <- mr(dat_k, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  print(res_k)
  fwrite(res_k, "results/PRIMARY_total_effect_perio_to_AD_relaxed_Kunkle.csv")
}

saveRDS(list(dat=dat, res=res, het=het, ple=ple, dir_test=dir_test, loo=loo),
        "data/PRIMARY_total_perio_ad_relaxed.rds")
cat("\n[OK] Done.\n")
