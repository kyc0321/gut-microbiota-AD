## 24_mrcml_steiger_tost.R
##  (a) MR-cML-BIC and MR-cML-DP (Xue et al. 2021): robust to correlated and uncorrelated pleiotropy
##  (b) Steiger-filtered IVW
##  (c) Equivalence testing (two one-sided tests) per doubling of odds of periodontitis,
##      margins OR 0.90-1.10 (and 0.95-1.05)
suppressPackageStartupMessages({ library(TwoSampleMR); library(MRcML); library(data.table) })
P <- readRDS("data/PRIMARY_total_perio_ad_relaxed.rds"); dat <- subset(P$dat, mr_keep)
n_out <- 487511; n_exp <- 250000
set.seed(20261007)
cml <- mr_cML(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, n = min(n_exp, n_out), random_start = 10)
dp  <- mr_cML_DP(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, n = min(n_exp, n_out), random_start = 10, num_pert = 200)
out <- data.table(method = c("MR-cML-BIC", "MR-cML-DP"),
                  b = c(cml$BIC_theta, dp$BIC_DP_theta), se = c(cml$BIC_se, dp$BIC_DP_se),
                  p = c(cml$BIC_p, dp$BIC_DP_p), n_invalid = c(length(cml$BIC_invalid), NA))
cat("cML invalid IVs (BIC):", paste(dat$SNP[cml$BIC_invalid], collapse = ", "), "\n")

## Steiger filtering (binary traits: use liability-scale r via get_r_from_lor)
dat$ncase.outcome <- 39106 + 46828; dat$ncontrol.outcome <- 401577
dat$prevalence.exposure <- 0.10; dat$prevalence.outcome <- 0.05
dat$r.exposure <- get_r_from_lor(dat$beta.exposure, dat$eaf.exposure, dat$ncase.exposure, dat$ncontrol.exposure, 0.10)
dat$r.outcome  <- get_r_from_lor(dat$beta.outcome, dat$eaf.outcome, dat$ncase.outcome, dat$ncontrol.outcome, 0.05)
st <- steiger_filtering(dat)
cat("Steiger: SNPs with correct direction:", sum(st$steiger_dir), "/", nrow(st), "\n")
r_st <- mr(subset(st, steiger_dir), method_list = "mr_ivw")
out <- rbind(out, data.table(method = sprintf("IVW, Steiger-filtered (%d SNPs)", r_st$nsnp), b = r_st$b, se = r_st$se, p = r_st$pval, n_invalid = NA))

## IVW reference
ivw <- mr(dat, method_list = "mr_ivw")
out <- rbind(data.table(method = "IVW (reference)", b = ivw$b, se = ivw$se, p = ivw$pval, n_invalid = NA), out)
out[, `:=`(OR = exp(b), LCI = exp(b - 1.96*se), UCI = exp(b + 1.96*se),
           OR_doubling = exp(b*log(2)), LCI_d = exp((b - 1.96*se)*log(2)), UCI_d = exp((b + 1.96*se)*log(2)))]

## TOST on the per-doubling scale
tost_p <- function(b, se, m) { bd <- b*log(2); sd <- se*log(2)
  p_lower <- pnorm((bd - log(1/m)) / sd, lower.tail = FALSE)   # H0: OR <= 1/m
  p_upper <- pnorm((bd - log(m)) / sd)                        # H0: OR >= m
  max(p_lower, p_upper) }
out[, TOST_p_1.10 := mapply(tost_p, b, se, 1.10)]
out[, TOST_p_1.05 := mapply(tost_p, b, se, 1.05)]
fwrite(out, "results/R_AHG/mrcml_steiger_tost.csv")
print(out[, .(method, OR = sprintf("%.3f (%.3f-%.3f)", OR, LCI, UCI), OR_doubling = sprintf("%.3f (%.3f-%.3f)", OR_doubling, LCI_d, UCI_d),
              p = signif(p, 2), TOST_1.10 = signif(TOST_p_1.10, 2), TOST_1.05 = signif(TOST_p_1.05, 2))])
