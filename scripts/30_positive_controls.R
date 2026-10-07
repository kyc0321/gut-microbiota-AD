## 30_positive_controls.R — positive-control MR (Bashir et al. 2025 checklist item):
## established causes of periodontitis -> periodontitis, using the same pipeline, to check that the
## periodontitis GWAS capture a phenotype that responds to known risk factors.
##   Exposures: smoking initiation (Liu 2019, GSCAN excl. 23andMe), type 2 diabetes (Xue 2018)
##              instruments P<5e-8, clumped r2<0.001 / 10 Mb (from 23a_mvmr_data_local.R)
##   Outcomes:  periodontitis, FinnGen (Salminen 2025) and GLIDE (Shungin 2019, EUR excl. HCHS/SOL)
suppressPackageStartupMessages({ library(TwoSampleMR); library(mr.raps); library(data.table) })
RFD <- "ext_data/mvmr_rf"
expo <- list(
  smoking = list(f = "SmkInit_Liu2019.txt.gz", c = c(SNP = "RSID", ea = "ALT", oa = "REF", eaf = "AF", b = "BETA", se = "SE", p = "PVALUE"),
                 clumped = "ext_data/mvmr_clump_smoking.clumped", label = "Smoking initiation (Liu 2019)"),
  T2D     = list(f = "T2D.h.tsv.gz", c = c(SNP = "hm_rsid", ea = "hm_effect_allele", oa = "hm_other_allele", eaf = "hm_effect_allele_frequency",
                 b = "hm_beta", se = "standard_error", p = "p_value"),
                 clumped = "ext_data/mvmr_clump_T2D.clumped", label = "Type 2 diabetes (Xue 2018)"))

## outcomes
fg <- fread("../FinGen Download/salminen_2025_sumstats/periodontitis_diagnosis_results.txt.gz",
            select = c("rsids","ref","alt","beta","sebeta","pval","af_alt"))
fg[, SNP := sub(",.*", "", rsids)]; fg <- fg[!duplicated(SNP)]
GLM <- "results/R_AHG/glide_mapped_for_controls.csv"
if (!file.exists(GLM)) {   # same mapping rules as 28_glide_replication.R
  g <- fread("ext_data/glide/EUR_perio_excl_HCHSSOL.txt")
  setnames(g, c("MarkerName","Allele1","Allele2","Effect","StdErr","P-value"), c("id","a1","a2","beta","se","p"))
  g <- g[N >= 0.8 * max(N)][, `:=`(a1 = toupper(a1), a2 = toupper(a2))]
  g <- g[!((a1 %in% c("A","T") & a2 %in% c("A","T")) | (a1 %in% c("C","G") & a2 %in% c("C","G")))]
  bim <- fread("ext_data/ld_ref/EUR.bim", col.names = c("chr","rsid","cm","pos","b1","b2"))[, id := paste(chr, pos, sep = ":")]
  g <- merge(g, bim[, .(id, rsid, b1, b2)], by = "id")
  g <- g[(a1 == b1 & a2 == b2) | (a1 == b2 & a2 == b1)][!duplicated(rsid)]
  fwrite(g[, .(rsid, a1, a2, beta, se, p, N)], GLM)
}
gl <- fread(GLM)

res <- list()
for (k in names(expo)) {
  e <- expo[[k]]
  ivs <- fread(e$clumped)$SNP
  d <- fread(file.path(RFD, e$f), select = unname(e$c)); setnames(d, unname(e$c), names(e$c))
  d <- d[SNP %in% ivs][!duplicated(SNP)]
  exp_dat <- data.frame(SNP = d$SNP, beta.exposure = as.numeric(d$b), se.exposure = as.numeric(d$se),
                        effect_allele.exposure = toupper(d$ea), other_allele.exposure = toupper(d$oa),
                        eaf.exposure = as.numeric(d$eaf), pval.exposure = as.numeric(d$p),
                        exposure = e$label, id.exposure = k, stringsAsFactors = FALSE)
  outs <- list(
    FinnGen = { o <- fg[SNP %in% exp_dat$SNP]
      data.frame(SNP = o$SNP, beta.outcome = o$beta, se.outcome = o$sebeta, effect_allele.outcome = o$alt,
                 other_allele.outcome = o$ref, eaf.outcome = o$af_alt, pval.outcome = o$pval,
                 outcome = "Periodontitis (FinnGen)", id.outcome = "FinnGen", stringsAsFactors = FALSE) },
    GLIDE = { o <- gl[rsid %in% exp_dat$SNP]
      data.frame(SNP = o$rsid, beta.outcome = o$beta, se.outcome = o$se, effect_allele.outcome = o$a1,
                 other_allele.outcome = o$a2, eaf.outcome = NA_real_, pval.outcome = o$p,
                 outcome = "Periodontitis (GLIDE)", id.outcome = "GLIDE", stringsAsFactors = FALSE) })
  for (oname in names(outs)) {
    # GLIDE has no allele frequencies and palindromes were removed upstream -> action = 1; FinnGen -> action = 2
    h <- harmonise_data(exp_dat, outs[[oname]], action = ifelse(oname == "GLIDE", 1, 2))
    h <- subset(h, mr_keep)
    r <- mr(h, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
    rp <- mr.raps.overdispersed.robust(h$beta.exposure, h$beta.outcome, h$se.exposure, h$se.outcome, loss.function = "huber")
    het <- mr_heterogeneity(h, method_list = "mr_ivw"); ple <- mr_pleiotropy_test(h)
    tab <- rbind(as.data.table(r)[, .(method, nsnp, b, se, p = pval)],
                 data.table(method = "MR-RAPS", nsnp = nrow(h), b = rp$beta.hat, se = rp$beta.se, p = 2 * pnorm(-abs(rp$beta.hat / rp$beta.se))))
    tab[, `:=`(exposure = e$label, outcome = oname, Q_p = het$Q_pval, egger_int_p = ple$pval)]
    res[[length(res) + 1]] <- tab
  }
}
res <- rbindlist(res)
## Scale: smoking initiation beta is on the log-odds scale of the GWAS (per unit log odds of smoking initiation);
## T2D beta is log OR -> OR of periodontitis per unit increase in log odds of liability to the exposure.
res[, `:=`(OR = exp(b), LCI = exp(b - 1.96 * se), UCI = exp(b + 1.96 * se))]
setcolorder(res, c("exposure", "outcome"))
fwrite(res, "results/R_AHG/positive_controls.csv")
print(res[, .(exposure = substr(exposure, 1, 18), outcome, method, nsnp, OR = sprintf("%.2f (%.2f-%.2f)", OR, LCI, UCI),
              p = signif(p, 2), Q_p = signif(Q_p, 2), egger_int_p = signif(egger_int_p, 2))])
