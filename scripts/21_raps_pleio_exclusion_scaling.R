## 21_raps_pleio_exclusion_scaling.R  (post-AHG editor comments 1-3)
##  (a) MR-RAPS (robust adjusted profile score; consistent under many weak instruments)
##  (b) IVW / WM / Egger re-run after excluding instruments with GW-significant
##      associations with AD risk factors (from 20_pleiotropy_lookup.R), and
##      after excluding ALL instruments with any other GW-significant association
##  (c) Effect scaling: per 1-unit log-odds increase in liability, and per
##      doubling of the odds of periodontitis (beta x ln2; Burgess & Labrecque 2018)
suppressPackageStartupMessages({ library(TwoSampleMR); library(mr.raps); library(data.table); library(dplyr) })
source("scripts/_helpers.R")
P  <- readRDS("data/PRIMARY_total_perio_ad_relaxed.rds")
dat <- subset(P$dat, mr_keep)
summ <- fread("results/R_AHG/instrument_pleiotropy_summary.csv")
excl_ad  <- summ[n_ad_relevant > 0, SNP]
excl_any <- summ[n_assoc > 0, SNP]

fit <- function(d, label) {
  r <- mr(d, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  rp <- mr.raps.overdispersed.robust(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome,
                loss.function = "huber")
  rbind(data.table(set = label, method = r$method, nsnp = r$nsnp, b = r$b, se = r$se, p = r$pval),
        data.table(set = label, method = "MR-RAPS (overdispersed, Huber)", nsnp = nrow(d),
                   b = rp$beta.hat, se = rp$beta.se, p = 2 * pnorm(-abs(rp$beta.hat / rp$beta.se))))
}
res <- rbind(
  fit(dat, "All instruments"),
  fit(subset(dat, !SNP %in% excl_ad),  sprintf("Excluding %d SNPs associated with AD risk factors", length(excl_ad))),
  fit(subset(dat, !SNP %in% excl_any), sprintf("Excluding %d SNPs with any other GW-significant association", length(excl_any))))
res[, `:=`(OR_per_logodds = exp(b), LCI = exp(b - 1.96 * se), UCI = exp(b + 1.96 * se),
           OR_per_doubling = exp(b * log(2)), LCI_d = exp((b - 1.96 * se) * log(2)), UCI_d = exp((b + 1.96 * se) * log(2)))]
fwrite(res, "results/R_AHG/raps_pleio_exclusion_scaling.csv")
print(res[, .(set = substr(set, 1, 40), method, nsnp, OR = sprintf("%.3f (%.3f-%.3f)", OR_per_logodds, LCI, UCI),
              OR_doubling = sprintf("%.3f (%.3f-%.3f)", OR_per_doubling, LCI_d, UCI_d), p = signif(p, 2))])

## Egger intercept in reduced sets
for (s in list(excl_ad, excl_any)) print(mr_pleiotropy_test(subset(dat, !SNP %in% s))[, c("egger_intercept","se","pval")])
