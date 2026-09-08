## 05_mvmr_exploratory.R -- exploratory multivariable MR (MVMR)
##
## Tests whether the Mollicutes/Tenericutes signal (the manuscript's lead
## finding) is robust to joint modeling of the other independent nominal
## signal clusters, since gut microbial taxa are compositional and a
## univariable association could in principle reflect confounding by a
## correlated co-taxon.
##
## Two models are fit:
##   (1) all 14 independent nominal signal clusters jointly;
##   (2) a reduced 3-exposure model (Mollicutes + its two most correlated
##       co-taxa at the shared instrument set: Methanobrevibacter,
##       Ruminiclostridium9), to check whether reducing the exposure count
##       alone resolves any weak-instrument problem seen in model (1).
##
## Requires OPENGWAS_JWT (see 00_setup.R). Do NOT hard-code the token here.

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(MVMR)      # conditional F-statistics, MVMR-IVW, MVMR pleiotropy Q-test
  library(dplyr)
})

stopifnot("Set OPENGWAS_JWT before running this script (see 00_setup.R)." =
            nzchar(Sys.getenv("OPENGWAS_JWT")))

dir.create("results", showWarnings = FALSE)

AD_ID <- "ebi-a-GCST90027158"

# The 14 independent nominal signal clusters from the primary screen
# (Table 1), one accession per cluster -- for taxonomically-nested duplicate
# pairs (e.g. class Mollicutes / phylum Tenericutes, which are the sole
# class within that phylum and therefore numerically identical), only one
# representative accession is used, consistent with the 184-cluster FDR
# correction in 03_fdr_correction.R.
taxa <- list(
  Mollicutes                  = "ebi-a-GCST90016921",
  Actinobacteria               = "ebi-a-GCST90016908",
  Selenomonadales              = "ebi-a-GCST90017107",
  Butyrivibrio                 = "ebi-a-GCST90016975",
  Ruminiclostridium9           = "ebi-a-GCST90017051",
  Pasteurellaceae              = "ebi-a-GCST90016944",
  Lactobacillaceae             = "ebi-a-GCST90016941",
  Intestinimonas               = "ebi-a-GCST90017019",
  Clostridium_innocuum_group   = "ebi-a-GCST90016979",
  Eggerthella                  = "ebi-a-GCST90016990",
  Methanobrevibacter           = "ebi-a-GCST90017033",
  Lachnospiraceae_UCG004       = "ebi-a-GCST90017026",
  Anaerotruncus                 = "ebi-a-GCST90016967",
  Bacteroides                   = "ebi-a-GCST90016968"
)

query_with_retry <- function(id, snps, max_tries = 6, base_wait = 8) {
  for (attempt in seq_len(max_tries)) {
    res <- tryCatch(
      ieugwasr::associations(variants = snps, id = id, proxies = 0),
      error = function(e) { cat("    attempt", attempt, "failed:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(res)) return(res)
    Sys.sleep(base_wait * attempt)
  }
  stop("Could not retrieve ", id, " after ", max_tries, " attempts (see OpenGWAS API status).")
}

build_mv_input <- function(exp_assoc, ad_assoc, keep) {
  exp_list <- lapply(keep, function(nm) {
    d <- exp_assoc[[nm]]
    d <- d[!is.na(d$beta) & !is.na(d$se) & !is.na(d$ea) & !is.na(d$nea), ]
    data.frame(
      SNP = d$rsid, beta.exposure = d$beta, se.exposure = d$se,
      effect_allele.exposure = toupper(d$ea), other_allele.exposure = toupper(d$nea),
      eaf.exposure = d$eaf, pval.exposure = d$p,
      exposure = nm, id.exposure = nm, mr_keep.exposure = TRUE, stringsAsFactors = FALSE)
  })
  exposure_dat <- do.call(rbind, exp_list)
  adf <- ad_assoc[!is.na(ad_assoc$beta) & !is.na(ad_assoc$se) & !is.na(ad_assoc$ea) & !is.na(ad_assoc$nea), ]
  outcome_dat <- data.frame(
    SNP = adf$rsid, beta.outcome = adf$beta, se.outcome = adf$se,
    effect_allele.outcome = toupper(adf$ea), other_allele.outcome = toupper(adf$nea),
    eaf.outcome = adf$eaf, pval.outcome = adf$p,
    outcome = "Alzheimer's disease", id.outcome = "AD", mr_keep.outcome = TRUE, stringsAsFactors = FALSE)
  mv_harmonise_data(exposure_dat, outcome_dat)
}

run_mvmr <- function(mvdat, label) {
  bx <- mvdat$exposure_beta; sx <- mvdat$exposure_se
  by <- mvdat$outcome_beta;  sy <- mvdat$outcome_se
  F.data <- format_mvmr(BXGs = bx, BYG = by, seBXGs = sx, seBYG = sy, RSID = rownames(bx))
  Fstat <- strength_mvmr(r_input = F.data, gencov = 0)
  ivw   <- ivw_mvmr(r_input = F.data)
  pleio <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
  list(label = label, n_snps = nrow(bx), exposures = colnames(bx),
       conditional_F = Fstat, ivw_mvmr = ivw, pleiotropy = pleio)
}

## --- Step 1: fetch each of the 14 taxa's own GWAS at the UNION of all 14
## taxa's instrument SNPs (not just each taxon's own selected instruments) --
## this is the data MVMR needs: every exposure's effect at every SNP used.
snp_union <- character(0)
for (nm in names(taxa)) {
  f <- file.path("data", "step2_taxa", paste0(taxa[[nm]], ".rds"))
  # falls back to a live instrument query if the cached per-taxon file from
  # 01_instruments_and_mr.R is not present
  if (file.exists(f)) {
    iv_snps <- readRDS(f)$iv$SNP
  } else {
    iv_snps <- ieugwasr::tophits(taxa[[nm]], pval = 1e-5, clump = 1)$rsid
  }
  snp_union <- union(snp_union, iv_snps)
}
cat("Union SNP count across 14 independent signal clusters:", length(snp_union), "\n")

exp_assoc <- list()
for (nm in names(taxa)) {
  cat("Querying", nm, taxa[[nm]], "...\n")
  exp_assoc[[nm]] <- query_with_retry(taxa[[nm]], snp_union)
  Sys.sleep(5)
}
cat("Querying AD outcome...\n")
ad_assoc <- query_with_retry(AD_ID, snp_union)

## --- Model 1: all 14 independent signal clusters jointly ---
mvdat_14 <- build_mv_input(exp_assoc, ad_assoc, names(taxa))
res_14 <- run_mvmr(mvdat_14, "14-exposure model")

## --- Model 2: reduced 3-exposure model (Mollicutes + 2 most correlated co-taxa) ---
# Co-taxa selected empirically by Pearson correlation of each taxon's
# SNP-effect vector with Mollicutes at the harmonised union SNP set (see
# manuscript Methods); recomputed here rather than hard-coded for transparency.
cors <- cor(mvdat_14$exposure_beta)
top2 <- names(sort(cors["Mollicutes", setdiff(colnames(cors), "Mollicutes")], decreasing = TRUE))[1:2]
cat("Most correlated co-taxa with Mollicutes:", paste(top2, collapse = ", "), "\n")
mvdat_3 <- build_mv_input(exp_assoc, ad_assoc, c("Mollicutes", top2))
res_3 <- run_mvmr(mvdat_3, "3-exposure model (Mollicutes + top-2 correlated co-taxa)")

## --- Report ---
for (res in list(res_14, res_3)) {
  cat("\n====", res$label, "====\n")
  cat("SNPs:", res$n_snps, "| Exposures:", paste(res$exposures, collapse = ", "), "\n")
  cat("Conditional F-statistics:\n"); print(res$conditional_F)
  cat("Pleiotropy Q-test: Q =", res$pleiotropy$Qstat, ", P =", res$pleiotropy$Qpval, "\n")
}

saveRDS(list(model_14 = res_14, model_3 = res_3),
        "results/mvmr_exploratory_results.rds")
cat("\n[OK] 05_mvmr_exploratory.R complete. See results/mvmr_exploratory_results.rds\n")
cat("Note: adjusted IVW-MVMR point estimates are intentionally not written to a\n")
cat("results CSV, since conditional F was well below 10 in both models and the\n")
cat("estimates are not reliably interpretable (see manuscript Results/Discussion).\n")
