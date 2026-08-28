## 02_sensitivity_diagnostics.R -- pleiotropy, heterogeneity, outlier, and
## direction diagnostics for every nominally significant (IVW P<0.05) taxon.
##
## Computes, per taxon: MR-Egger intercept (pleiotropy), Cochran's Q
## (heterogeneity), MR-PRESSO global test (outlier-driven pleiotropy), and
## Steiger directionality test. Table 2 in the manuscript reports this full
## battery for the lead (smallest-P) signal only; this script computes it for
## all nominal taxa so the full evidence base is reproducible, not just the
## single tabulated row.

source("00_setup.R")

ivw_tbl <- fread("results/ivw_by_taxon.csv")
nominal <- ivw_tbl[beta_p < 0.05][order(beta_p)]
cat("Running full sensitivity battery for", nrow(nominal), "nominally significant taxa.\n")

run_sensitivity <- function(dat) {
  list(
    pleio = TwoSampleMR::mr_pleiotropy_test(dat),
    het   = TwoSampleMR::mr_heterogeneity(dat),
    dir   = TwoSampleMR::directionality_test(dat),
    presso = tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = as.data.frame(dat),
        NbDistribution = 1000, SignifThreshold = 0.05),
      error = function(e) list(error = conditionMessage(e)))
  )
}

sens_results <- list()
for (i in seq_len(nrow(nominal))) {
  tid <- nominal$taxon_id[i]; tname <- nominal$taxon[i]
  cached <- readRDS(file.path("data/taxa", paste0(tid, ".rds")))
  if (is.null(cached$dat)) next
  cat(sprintf("[%2d/%d] %s\n", i, nrow(nominal), tname))
  sens_results[[tid]] <- c(list(taxon = tname, taxon_id = tid), run_sensitivity(cached$dat))
}
saveRDS(sens_results, "results/sensitivity_full.rds")

summary_tbl <- bind_rows(lapply(sens_results, function(s) {
  data.frame(
    taxon = s$taxon, taxon_id = s$taxon_id,
    egger_intercept   = if (nrow(s$pleio)) s$pleio$egger_intercept[1] else NA,
    egger_intercept_p = if (nrow(s$pleio)) s$pleio$pval[1] else NA,
    Q_ivw   = if (nrow(s$het)) s$het$Q[s$het$method == "Inverse variance weighted"][1] else NA,
    Q_ivw_p = if (nrow(s$het)) s$het$Q_pval[s$het$method == "Inverse variance weighted"][1] else NA,
    steiger_correct_dir = if (nrow(s$dir)) s$dir$correct_causal_direction[1] else NA,
    steiger_p           = if (nrow(s$dir)) s$dir$steiger_pval[1] else NA,
    presso_global_p = if (!is.null(s$presso$`MR-PRESSO results`))
      s$presso$`MR-PRESSO results`$`Global Test`$Pvalue else NA
  )
}))
fwrite(summary_tbl, "results/sensitivity_summary.csv")
cat("[OK] 02_sensitivity_diagnostics.R complete.\n")
