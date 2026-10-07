## 26_lcv.R — Latent causal variable model (O'Connor & Price 2018): genetic causality
## proportion (GCP) of periodontitis on AD. HapMap3 SNPs, MHC excluded, LD scores from
## LDSC eur_w_ld_chr; SNPs ordered by position for block jackknife.
suppressPackageStartupMessages(library(data.table))
hm3 <- fread("ext_data/eur_w_ld_chr/w_hm3.snplist")
ld  <- rbindlist(lapply(1:22, function(c) fread(sprintf("ext_data/eur_w_ld_chr/%d.l2.ldscore.gz", c))))
ex  <- readRDS("ext_data/mrlap_exposure_perio.rds")
ad  <- fread("ext_data/bellenguez/35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
             select = c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error"))
setnames(ad, c("rsid","ea","oa","b","se"))
m <- merge(hm3, ld[, .(SNP, CHR, BP, L2)], by = "SNP")
m <- merge(m, ex[, .(SNP = rsid, alt, ref, b1 = beta, s1 = se)], by = "SNP")
m <- merge(m, ad[!is.na(b), .(SNP = rsid, ea, oa, b2 = b, s2 = se)], by = "SNP")
m[, z1 := fifelse(alt == A1 & ref == A2, b1/s1, fifelse(alt == A2 & ref == A1, -b1/s1, NA_real_))]
m[, z2 := fifelse(ea == A1 & oa == A2, b2/s2, fifelse(ea == A2 & oa == A1, -b2/s2, NA_real_))]
m <- m[!is.na(z1) & !is.na(z2) & !(CHR == 6 & BP > 25e6 & BP < 34e6)]
m <- m[!(CHR == 19 & BP > 44.4e6 & BP < 46.5e6)]   # APOE region: extreme AD chi2 dominates moments
setorder(m, CHR, BP)
cat("SNPs:", nrow(m), "\n")
owd <- setwd("ext_data/LCV"); source("RunLCV.R"); 
set.seed(1)
L <- RunLCV(m$L2, m$z1, m$z2)
setwd(owd)
out <- data.table(n_snps = nrow(m), rho = L$rho.est, rho_se = L$rho.err, gcp = L$gcp.pm, gcp_se = L$gcp.pse,
                  p_fullycausal = L$pval.fullycausal[1], p_gcp0 = L$pval.gcpzero.2tailed,
                  h2_zscore_perio = L$h2.zscore[1], h2_zscore_AD = L$h2.zscore[2])
fwrite(out, "results/R_AHG/lcv_results.csv"); print(t(out))
