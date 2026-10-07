## 18_power_corrected_neff.R -- Corrected analytic power (Brion et al. 2013).
## Earlier versions took OpenGWAS metadata for ebi-a-GCST90027158 at face value
## (ncase 39,106 / ncontrol 46,828). The GWAS Catalog record shows 46,828 are
## PROXY cases; controls are 401,577. Conservative effective N here uses the
## clinically diagnosed cases only versus all controls.
power_mr <- function(b, R2_X, N_Y, alpha = 0.05) {
  ncp <- abs(b) * sqrt(R2_X * N_Y); z_a <- qnorm(1 - alpha/2)
  pnorm(ncp - z_a) + pnorm(-ncp - z_a)
}
mde_mr <- function(R2_X, N_Y, power = 0.8, alpha = 0.05)
  (qnorm(1 - alpha/2) + qnorm(power)) / sqrt(R2_X * N_Y)

N_ad  <- 4 * 39106 * 401577 / (39106 + 401577)
N_mb  <- 14306
R2 <- c(prior = 0.001538, primary = 0.002357, gws = 0.000207)
cat(sprintf("AD effective N (clinical cases vs controls) = %.0f\n\n", N_ad))
for (nm in names(R2)) {
  cat(sprintf("%-8s R2=%.6f  MDE OR=%.3f  power OR1.15=%.1f%%  OR1.25=%.1f%%  Step1 power(alpha=0.092)=%.1f%%  Step1 MDE=%.3f\n",
      nm, R2[nm], exp(mde_mr(R2[nm], N_ad)),
      100 * power_mr(log(1.15), R2[nm], N_ad), 100 * power_mr(log(1.25), R2[nm], N_ad),
      100 * power_mr(0.092, R2[nm], N_mb), mde_mr(R2[nm], N_mb)))
}
