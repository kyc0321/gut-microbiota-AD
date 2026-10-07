## 12_reverse_mr_relaxed.R -- Reverse MR: AD -> Periodontitis (Salminen 2025)
## using local summary stats (no API needed; outcome already downloaded).
suppressPackageStartupMessages({ library(TwoSampleMR); library(data.table); library(dplyr) })

ad_iv <- fread("data/ad_instruments.csv")
raw <- fread("../FinGen Download/salminen_2025_sumstats/periodontitis_diagnosis_results.txt.gz")

out <- raw[rsids %in% ad_iv$SNP]
out <- out[!duplicated(rsids)]
cat("AD SNPs found in periodontitis outcome:", nrow(out), "of", nrow(ad_iv), "\n")

outcome_dat <- data.frame(
  SNP = out$rsids, chr = out$chrom, pos = out$pos,
  effect_allele.outcome = out$alt, other_allele.outcome = out$ref,
  eaf.outcome = out$af_alt, beta.outcome = out$beta, se.outcome = out$sebeta,
  pval.outcome = out$pval, outcome = "Periodontitis (Salminen 2025)", id.outcome = "Salminen2025_perio_diagnosis"
)

exposure_dat <- ad_iv
class(exposure_dat) <- "data.frame"

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
cat("SNPs retained after harmonisation:", nrow(dat), "\n")
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
print(res)
fwrite(res, "results/PRIMARY_reverse_AD_to_perio_relaxed.csv")
