## 23c_mvmr_qhet.R — weak-instrument-robust MVMR (Q-statistic minimisation; Sanderson et al. 2021,
## MVMR::qhet_mvmr) for the two-exposure models, because conditional F for periodontitis was < 10.
## Exposure GWAS come from different samples, so the phenotypic correlation matrix is set to identity.
suppressPackageStartupMessages({ library(MVMR); library(data.table) })
M <- fread("results/R_AHG/mvmr_harmonised_matrix.csv")
RF <- intersect(c("smoking","T2D","BMI","LDL","EA","T1D"), sub("^b_", "", grep("^b_", names(M), value = TRUE)))
set.seed(20261007)
res <- rbindlist(lapply(RF, function(k) {
  d <- M[perio_iv | get(paste0(k, "_iv"))]
  bx <- c("b_perio", paste0("b_", k)); sx <- c("se_perio", paste0("se_", k))
  d <- d[complete.cases(d[, c(bx, sx, "b_AD", "se_AD"), with = FALSE])]
  fm <- format_mvmr(BXGs = as.matrix(d[, ..bx]), BYG = d$b_AD, seBXGs = as.matrix(d[, ..sx]), seBYG = d$se_AD, RSID = d$rsid)
  q <- suppressWarnings(qhet_mvmr(fm, pcor = diag(2), CI = TRUE, iterations = 500, ncores = 4))
  ci <- as.numeric(regmatches(q[1, 2], regexec("^(-?[0-9.]+)-(-?[0-9.]+)$", q[1, 2]))[[1]][2:3])   # "lo-hi" string
  data.table(model = paste("periodontitis +", k), nsnp = nrow(d), b = q[1, 1], LCI_b = ci[1], UCI_b = ci[2])
}))
res[, `:=`(OR = exp(b), LCI = exp(LCI_b), UCI = exp(UCI_b))]
fwrite(res, "results/R_AHG/mvmr_qhet.csv")
print(res[, .(model, nsnp, OR = sprintf("%.3f (%.3f-%.3f)", OR, LCI, UCI))])
