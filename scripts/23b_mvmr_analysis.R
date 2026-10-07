## 23b_mvmr_analysis.R — MVMR: direct effect of periodontitis on AD adjusting for each AD risk
## factor separately and for all jointly. Estimators: MV-IVW, MV-Egger, MV-median (MendelianRandomization);
## conditional F-statistics (Sanderson et al. 2019/2021; MVMR package, gencov = 0 approximation).
suppressPackageStartupMessages({ library(MendelianRandomization); library(MVMR); library(data.table) })
M <- fread("results/R_AHG/mvmr_harmonised_matrix.csv")
RF <- intersect(c("smoking","T2D","BMI","LDL","EA","T1D"), sub("^b_", "", grep("^b_", names(M), value = TRUE)))
fitset <- function(rfs, label) {
  ivcol <- c("perio_iv", paste0(rfs, "_iv"))
  d <- M[Reduce(`|`, lapply(ivcol, function(v) get(v)))]
  bx <- c("b_perio", paste0("b_", rfs)); sx <- c("se_perio", paste0("se_", rfs))
  d <- d[complete.cases(d[, c(bx, sx, "b_AD", "se_AD"), with = FALSE])]
  inp <- mr_mvinput(bx = as.matrix(d[, ..bx]), bxse = as.matrix(d[, ..sx]), by = d$b_AD, byse = d$se_AD,
                    exposure = c("periodontitis", rfs), outcome = "AD")
  ivw <- mr_mvivw(inp); egg <- mr_mvegger(inp); med <- mr_mvmedian(inp, iterations = 1000)
  F <- tryCatch({
    fm <- format_mvmr(BXGs = as.matrix(d[, ..bx]), BYG = d$b_AD, seBXGs = as.matrix(d[, ..sx]), seBYG = d$se_AD, RSID = d$rsid)
    as.numeric(strength_mvmr(fm, gencov = 0)) }, error = function(e) rep(NA, length(bx)))
  data.table(model = label, nsnp = nrow(d),
             method = c("MV-IVW","MV-Egger","MV-median"),
             b = c(ivw@Estimate[1], egg@Estimate[1], med@Estimate[1]),
             se = c(ivw@StdError[1], egg@StdError.Est[1], med@StdError[1]),
             p = c(ivw@Pvalue[1], egg@Pvalue.Est[1], med@Pvalue[1]),
             egger_intercept_p = c(NA, egg@Pvalue.Int, NA),
             condF_perio = F[1], condF_others = paste(round(F[-1], 1), collapse = "/"),
             rf_estimates = paste(sprintf("%s:%.3f(P=%.2g)", rfs, ivw@Estimate[-1], ivw@Pvalue[-1]), collapse = "; "))
}
res <- rbindlist(c(lapply(RF, function(k) fitset(k, paste("periodontitis +", k))),
                   list(fitset(RF, sprintf("periodontitis + all %d", length(RF))))))
res[, `:=`(OR = exp(b), LCI = exp(b - 1.96*se), UCI = exp(b + 1.96*se), OR_doubling = exp(b*log(2)))]
fwrite(res, "results/R_AHG/mvmr_results.csv")
print(res[, .(model, method, nsnp, OR = sprintf("%.3f (%.3f-%.3f)", OR, LCI, UCI), p = signif(p, 2), condF_perio = round(condF_perio, 1), condF_others)])
