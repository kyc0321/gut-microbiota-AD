## 01_instruments_and_mr.R -- instrument extraction + primary MR for all 211 taxa
##
## For each MiBioGen taxon: extract instruments at P<1e-5 (r2<0.001, 10,000kb
## clumping, European 1000 Genomes reference), harmonise against the AD outcome
## GWAS (palindromic SNPs with intermediate allele frequency excluded), and
## estimate the causal effect via IVW, MR-Egger, and weighted median.
## Incremental/resumable: caches each taxon's result to data/taxa/<id>.rds.

source("00_setup.R")

taxa_index <- ieugwasr::gwasinfo() %>%
  filter(grepl("MiBioGen|Kurilshikov", paste(consortium, author), ignore.case = TRUE)) %>%
  select(id, trait)
cat("MiBioGen taxa found:", nrow(taxa_index), "\n")
fwrite(taxa_index, "data/taxa_index.csv")

dir.create("data/taxa", showWarnings = FALSE)

for (i in seq_len(nrow(taxa_index))) {
  tid <- taxa_index$id[i]; tname <- taxa_index$trait[i]
  fp <- file.path("data/taxa", paste0(tid, ".rds"))
  if (file.exists(fp)) next
  cat(sprintf("[%3d/%d] %s ... ", i, nrow(taxa_index), tid))

  iv <- tryCatch(
    TwoSampleMR::extract_instruments(tid, p1 = 1e-5, clump = TRUE, r2 = 0.001, kb = 10000),
    error = function(e) NULL)
  if (is.null(iv) || nrow(iv) < 3) {
    saveRDS(list(taxon = tname, taxon_id = tid, status = "no_iv"), fp)
    cat("skip (IVs<3)\n"); next
  }
  iv$F_stat <- (iv$beta.exposure / iv$se.exposure)^2

  out <- tryCatch(
    TwoSampleMR::extract_outcome_data(snps = iv$SNP, outcomes = AD_OUTCOME_ID, proxies = TRUE),
    error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0) {
    saveRDS(list(taxon = tname, taxon_id = tid, status = "no_outcome", iv = iv), fp)
    cat("no outcome SNPs\n"); next
  }
  dat <- tryCatch(TwoSampleMR::harmonise_data(iv, out, action = 2), error = function(e) NULL)
  if (is.null(dat) || sum(dat$mr_keep) < 3) {
    saveRDS(list(taxon = tname, taxon_id = tid, status = "harmonise_lt3", iv = iv, out = out), fp)
    cat("harmonise<3\n"); next
  }

  r <- TwoSampleMR::mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  r$taxon <- tname; r$taxon_id <- tid
  saveRDS(list(taxon = tname, taxon_id = tid, status = "ok", iv = iv, dat = dat, mr = r), fp)
  cat(sprintf("nIV=%d  IVW b=%.4f p=%.4g\n",
              sum(dat$mr_keep),
              r$b[r$method == "Inverse variance weighted"],
              r$pval[r$method == "Inverse variance weighted"]))
}

# --- Aggregate ---
files <- list.files("data/taxa", full.names = TRUE)
all_taxa <- lapply(files, readRDS)
ok <- Filter(function(x) x$status == "ok", all_taxa)
cat("\nTestable taxa (>=3 instruments, harmonised):", length(ok), "of", nrow(taxa_index), "\n")

res_all <- bind_rows(lapply(ok, `[[`, "mr"))
fwrite(res_all, "results/mr_all_methods.csv")

ivw_tbl <- res_all %>%
  filter(method == "Inverse variance weighted") %>%
  transmute(taxon, taxon_id, beta = b, beta_se = se, beta_p = pval, beta_nsnp = nsnp)
fwrite(ivw_tbl, "results/ivw_by_taxon.csv")

cat("Nominal (P<0.05) IVW associations:", sum(ivw_tbl$beta_p < 0.05, na.rm = TRUE), "\n")
cat("[OK] 01_instruments_and_mr.R complete.\n")
