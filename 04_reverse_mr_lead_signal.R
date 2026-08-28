## 04_reverse_mr_lead_signal.R -- reverse MR (AD as exposure, taxon abundance
## as outcome) for the single taxon with the smallest forward-direction P-value,
## using genome-wide-significant (P<5e-8) AD instruments.

source("00_setup.R")

fdr_tbl <- fread("results/ivw_by_taxon_with_fdr.csv")[order(beta_p)]
lead_cluster <- fdr_tbl$cluster_key[1]
lead_taxa <- fdr_tbl[cluster_key == lead_cluster]
cat("Lead signal (smallest P):\n"); print(lead_taxa[, .(taxon, taxon_id, beta, beta_p, q)])

ad_iv <- TwoSampleMR::extract_instruments(AD_OUTCOME_ID, p1 = AD_REVERSE_INSTR_P,
                                           clump = TRUE, r2 = 0.001, kb = 10000)
cat("AD instruments (P<5e-8):", nrow(ad_iv), "\n")
fwrite(ad_iv, "data/ad_instruments.csv")

out_rev <- TwoSampleMR::extract_outcome_data(snps = ad_iv$SNP,
                                              outcomes = lead_taxa$taxon_id, proxies = TRUE)
dat_rev <- TwoSampleMR::harmonise_data(ad_iv, out_rev, action = 2)
res_rev <- TwoSampleMR::mr(dat_rev, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
res_rev$taxon <- lead_taxa$taxon[match(res_rev$id.outcome, lead_taxa$taxon_id)]

cat("\n=== Reverse MR: AD -> lead signal taxon ===\n")
print(res_rev[, c("taxon", "method", "nsnp", "b", "se", "pval")])
fwrite(res_rev, "results/reverse_mr_lead_signal.csv")
cat("[OK] 04_reverse_mr_lead_signal.R complete.\n")
