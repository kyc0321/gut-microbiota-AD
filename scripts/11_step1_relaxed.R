## 11_step1_relaxed.R — Step 1 (Periodontitis[Salminen2025,relaxed] -> 211 taxa)
## EXPLORATORY/SECONDARY per strategic decision: total effect (script 10) is
## the primary finding; this mediation pathway is reported as exploratory,
## explicitly caveated by the microbiome GWAS power ceiling (N=14,306,
## Step1 power ~8.3% -- see power_recalc2.R), and multiplicity-corrected
## using the Meff (effective independent tests) from 06_taxa_independence_meff.R
## rather than naive Bonferroni/211.
suppressPackageStartupMessages({
  library(TwoSampleMR); library(ieugwasr); library(dplyr); library(data.table)
})
source("scripts/_helpers.R")

idx <- readRDS("data/datasets_index.rds"); gut <- idx$gut
perio_iv <- readRDS("data/perio_instruments_finngen2025_relaxed.rds")
perio_iv$ncontrol.exposure <- perio_iv$samplesize.exposure - perio_iv$ncase.exposure
perio_iv$units.exposure <- "log odds"

CHUNK <- 20
chunks <- split(gut$id, ceiling(seq_along(gut$id) / CHUNK))
cat("Step1 (relaxed instrument, n=", nrow(perio_iv), " SNPs) in", length(chunks), "chunks\n")

dir.create("data/step1_chunks_relaxed", showWarnings = FALSE)
for (i in seq_along(chunks)) {
  fp <- sprintf("data/step1_chunks_relaxed/chunk_%02d.rds", i)
  if (file.exists(fp)) { cat(sprintf("[%2d/%d] cached\n", i, length(chunks))); next }
  cat(sprintf("[%2d/%d] %d outcomes ... ", i, length(chunks), length(chunks[[i]])))
  t0 <- Sys.time()
  res <- tryCatch(
    TwoSampleMR::extract_outcome_data(snps = perio_iv$SNP,
                                      outcomes = chunks[[i]], proxies = TRUE),
    error = function(e) { cat("ERR:", conditionMessage(e)); NULL })
  if (!is.null(res)) cat(sprintf("rows=%d  %.0fs\n", nrow(res), as.numeric(Sys.time()-t0,units="secs")))
  saveRDS(res, fp)
  Sys.sleep(1)
}

all_chunks <- lapply(list.files("data/step1_chunks_relaxed", full.names = TRUE), readRDS)
out_all <- bind_rows(Filter(Negate(is.null), all_chunks))
cat("\nCombined Step1 outcome rows:", nrow(out_all),
    "  unique taxa:", length(unique(out_all$id.outcome)), "\n")
saveRDS(out_all, "data/step1_outcomes_relaxed.rds")

dat1 <- TwoSampleMR::harmonise_data(perio_iv, out_all, action = 2)
res1 <- TwoSampleMR::mr(dat1, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
res1$taxon <- gut$trait[match(res1$id.outcome, gut$id)]
fwrite(res1, "results/EXPLORATORY_step1_perio_to_taxa_relaxed.csv")

alpha_tbl <- res1 %>% filter(method == "Inverse variance weighted") %>%
  transmute(taxon, taxon_id = id.outcome, alpha = b, alpha_se = se,
            alpha_p = pval, alpha_nsnp = nsnp)

## Apply M_eff-based correction instead of naive Bonferroni/211. M_eff comes from
## 14_taxa_independence_meff_relaxed.R (Li & Ji 32.3 -> alpha 1.5e-3), which needs
## this script's Step 1 outputs; run 11 -> 14 -> 11 so the flag uses that file.
## (An earlier version read the M_eff file from the superseded instrument,
## alpha 3.1e-3; no taxon is significant under either threshold, min P = 0.019.)
meff_file <- "results/taxa_meff_summary_relaxed.csv"
alpha_meff <- if (file.exists(meff_file)) {
  meff_tbl <- fread(meff_file)
  meff_tbl$bonferroni_alpha[meff_tbl$method == "Li & Ji (2005)"]
} else NA_real_
alpha_naive <- 0.05/211
alpha_tbl$sig_naive_211 <- alpha_tbl$alpha_p < alpha_naive
alpha_tbl$sig_meff_LiJi <- alpha_tbl$alpha_p < alpha_meff
fwrite(alpha_tbl, "results/EXPLORATORY_step1_alpha_IVW_relaxed.csv")

cat("\nStep1 IVW: taxa tested =", nrow(alpha_tbl),
    " | nominal-sig (p<0.05) =", sum(alpha_tbl$alpha_p < 0.05, na.rm = TRUE),
    " | sig at naive Bonferroni (a=", signif(alpha_naive,3), ") =", sum(alpha_tbl$sig_naive_211, na.rm=TRUE),
    " | sig at Meff-corrected (a=", signif(alpha_meff,3), ") =", sum(alpha_tbl$sig_meff_LiJi, na.rm=TRUE), "\n")

cat("\n=== Top 10 taxa by p-value ===\n")
print(head(alpha_tbl[order(alpha_tbl$alpha_p), c("taxon","alpha","alpha_p","alpha_nsnp")], 10))

cat("\n[OK] 11_step1_relaxed.R complete.\n")
