## 29_conmix_lasso.R — contamination mixture (Burgess et al. 2020) and MR-Lasso (Rees et al. 2019)
suppressPackageStartupMessages({ library(MendelianRandomization); library(data.table) })
P <- readRDS("data/PRIMARY_total_perio_ad_relaxed.rds"); dat <- subset(P$dat, mr_keep)
inp <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome, snps = dat$SNP)
cm <- mr_conmix(inp); la <- mr_lasso(inp)
out <- data.table(method = c("MR-ConMix", "MR-Lasso"), nsnp = c(cm@SNPs, la@SNPs),
                  b = c(cm@Estimate, la@Estimate), LCI_b = c(cm@CILower[1], la@CILower), UCI_b = c(cm@CIUpper[length(cm@CIUpper)], la@CIUpper),
                  p = c(cm@Pvalue, la@Pvalue), valid_snps = c(NA, la@Valid))
out[, `:=`(OR = exp(b), LCI = exp(LCI_b), UCI = exp(UCI_b), OR_doubling = exp(b*log(2)))]
cat("ConMix CI is", ifelse(length(cm@CILower) > 1, "multimodal (union of intervals)", "a single interval"), "\n")
fwrite(out, "results/R_AHG/conmix_lasso.csv")
print(out[, .(method, nsnp, valid_snps, OR = sprintf("%.3f (%.3f-%.3f)", OR, LCI, UCI), p = signif(p, 2))])
