## 15_mrpresso_primary.R
## MR-PRESSO on the exact harmonized dataset underlying the PRIMARY total-effect
## analysis (10_primary_total_effect_relaxed.R) -- 50-SNP relaxed-threshold
## (P<1e-5) Salminen 2025 FinnGen periodontitis instrument -> Bellenguez 2022 AD.
## Run because MR-Egger intercept was nominally significant (P=0.021),
## motivating an outlier-corrected pleiotropy check.
suppressPackageStartupMessages({
  library(MRPRESSO); library(data.table)
})

obj <- readRDS("data/PRIMARY_total_perio_ad_relaxed.rds")
dat <- obj$dat
dat <- dat[dat$mr_keep, ]  # exclude the 1 SNP with an unresolvable allele mismatch (matches mr() behavior)
cat("Harmonised SNPs available (mr_keep only):", nrow(dat), "\n")

set.seed(12345)
presso <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome", SdExposure = "se.exposure",
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  data = dat, NbDistribution = 5000, SignifThreshold = 0.05
)

cat("\n=== MR-PRESSO Main MR results ===\n")
print(presso$`Main MR results`)

cat("\n=== MR-PRESSO Global Test ===\n")
print(presso$`MR-PRESSO results`$`Global Test`)

outlier_test <- presso$`MR-PRESSO results`$`Outlier Test`
if (!is.null(outlier_test)) {
  cat("\n=== Outlier Test (per-SNP) ===\n")
  outlier_test$SNP <- dat$SNP
  print(outlier_test[order(outlier_test$Pvalue), ])
  fwrite(outlier_test, "results/PRIMARY_mrpresso_outliertest_relaxed.csv")
}

main_res <- as.data.frame(presso$`Main MR results`)
fwrite(main_res, "results/PRIMARY_mrpresso_mainresults_relaxed.csv")

global_p <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
cat("\nGlobal Test RSSobs:", presso$`MR-PRESSO results`$`Global Test`$RSSobs,
    "| P-value:", global_p, "\n")

distortion <- presso$`MR-PRESSO results`$`Distortion Test`
if (!is.null(distortion) && !all(is.na(distortion))) {
  cat("\n=== Distortion Test ===\n")
  print(distortion)
}

saveRDS(presso, "data/PRIMARY_mrpresso_relaxed.rds")
cat("\n[OK] Done.\n")
