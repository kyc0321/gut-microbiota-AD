## 09_new_instruments_finngen_relaxed.R
## Two-tier instrument strategy for the new Salminen 2025 FinnGen
## periodontitis-diagnosis exposure:
##   - PRIMARY: P<1e-5 (relaxed, same threshold convention as original
##     analysis), clumped -- restores enough instrument SNPs for adequate
##     aggregate R2/power while keeping the clinically valid exposure.
##   - SENSITIVITY: P<5e-8 genome-wide-significant subset (3 SNPs, see
##     07_new_instruments_finngen.R) -- cleanest instruments, addresses
##     weak-instrument concern directly.
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(TwoSampleMR); library(ieugwasr)
})

RAW <- "../FinGen Download/salminen_2025_sumstats/periodontitis_diagnosis_results.txt.gz"
N_CASES <- 38157
N_TOTAL_APPROX <- 473681

cat("Reading full GWAS and filtering to P<1e-5...\n")
dt <- fread(RAW, sep = "\t", header = TRUE)
relaxed <- dt[pval < 1e-5]
cat("Rows at P<1e-5:", nrow(relaxed), "\n")

relaxed[, rsid := rsids]
clump_input <- relaxed[, .(rsid, pval)]

cat("\nRunning LD clumping (r2=0.001, kb=10000)...\n")
clumped <- ld_clump(
  clump_input,
  clump_kb = 10000, clump_r2 = 0.001, clump_p = 1e-5,
  pop = "EUR"
)
cat("Independent instruments after clumping:", nrow(clumped), "\n")

iv <- relaxed[rsids %in% clumped$rsid]
iv <- iv[!duplicated(rsids)]

perio_iv_relaxed <- data.frame(
  SNP = iv$rsids,
  chr.exposure = iv$chrom,
  pos.exposure = iv$pos,
  effect_allele.exposure = iv$alt,
  other_allele.exposure = iv$ref,
  eaf.exposure = iv$af_alt,
  beta.exposure = iv$beta,
  se.exposure = iv$sebeta,
  pval.exposure = iv$pval,
  samplesize.exposure = N_TOTAL_APPROX,
  ncase.exposure = N_CASES,
  exposure = "Periodontitis (Salminen 2025 FinnGen diagnosis-based, P<1e-5)",
  id.exposure = "Salminen2025_perio_diagnosis_relaxed"
)
perio_iv_relaxed$F_stat <- (perio_iv_relaxed$beta.exposure / perio_iv_relaxed$se.exposure)^2

cat("\n=== Final relaxed-threshold instrument set ===\n")
print(perio_iv_relaxed[, c("SNP","chr.exposure","beta.exposure","se.exposure","pval.exposure","F_stat")])
cat(sprintf("\nn=%d  Mean F: %.1f   Min F: %.1f   Max F: %.1f\n",
            nrow(perio_iv_relaxed), mean(perio_iv_relaxed$F_stat),
            min(perio_iv_relaxed$F_stat), max(perio_iv_relaxed$F_stat)))

# Aggregate R2 (sum of per-SNP R2, standard approach)
R2_per_snp <- with(perio_iv_relaxed, (2 * beta.exposure^2 * eaf.exposure * (1-eaf.exposure)) /
  (2 * beta.exposure^2 * eaf.exposure * (1-eaf.exposure) + se.exposure^2 * 2 * samplesize.exposure * eaf.exposure * (1-eaf.exposure)))
R2_agg <- sum(R2_per_snp)
cat(sprintf("Aggregate R2 (sum across SNPs): %.6f\n", R2_agg))

fwrite(perio_iv_relaxed, "data/perio_instruments_finngen2025_relaxed.csv")
saveRDS(perio_iv_relaxed, "data/perio_instruments_finngen2025_relaxed.rds")
cat("\n[OK] Saved data/perio_instruments_finngen2025_relaxed.csv\n")
