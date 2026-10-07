## 27_cause.R — CAUSE (Morrison et al. 2020): compares a causal model with a model of
## correlated horizontal pleiotropy through a shared factor. Uses genome-wide HapMap3 SNPs.
## Nuisance parameters from all merged SNPs; model fitted on SNPs LD-pruned locally with
## PLINK 1.9 (1000 Genomes EUR; r2 < 0.01, 10 Mb, periodontitis P < 1e-3).
suppressPackageStartupMessages({ library(cause); library(data.table) })
## Compatibility patch: loo >= 2.8 returns loo_compare() as a data.frame whose first column is the
## model name, so cause 1.2.0's positional indexing (comp[2, 1]) fails. Index by column name instead.
local({
  f <- cause:::in_sample_elpd_loo
  b <- deparse(body(f))
  b <- gsub('rownames(comp)[1] == "model2"', 'comp$model[1] == "model2"', b, fixed = TRUE)
  b <- gsub("comp[2, 1]", "comp$elpd_diff[2]", b, fixed = TRUE)
  b <- gsub("comp[2, 2]", "comp$se_diff[2]", b, fixed = TRUE)
  body(f) <- parse(text = paste(b, collapse = "\n"))[[1]]
  assignInNamespace("in_sample_elpd_loo", f, ns = "cause")
})
ex <- readRDS("ext_data/mrlap_exposure_perio.rds")
ad <- fread("ext_data/bellenguez/35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
            select = c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error"))
setnames(ad, c("rsid","ea","oa","b","se")); ad <- ad[!is.na(b) & rsid %in% ex$rsid][!duplicated(rsid)]
X <- gwas_merge(as.data.frame(ex), as.data.frame(ad),
                snp_name_cols = c("rsid","rsid"), beta_hat_cols = c("beta","b"), se_cols = c("se","se"),
                A1_cols = c("alt","ea"), A2_cols = c("ref","oa"))
cat("Merged SNPs:", nrow(X), "\n")
set.seed(100)
params <- est_cause_params(X, sample(X$snp, min(1e6, nrow(X))))
## local LD pruning
X$p1 <- 2 * pnorm(-abs(X$beta_hat_1 / X$seb1))
fwrite(data.table(SNP = X$snp, P = X$p1), "ext_data/cause_clump_input.txt", sep = "\t")
system2("ext_data/plink/plink", c("--bfile", "ext_data/ld_ref/EUR", "--clump", "ext_data/cause_clump_input.txt",
        "--clump-p1", "1e-3", "--clump-p2", "1", "--clump-r2", "0.01", "--clump-kb", "10000",
        "--out", "ext_data/cause_clump"), stdout = FALSE)
top <- fread("ext_data/cause_clump.clumped")$SNP
cat("Pruned variants:", length(top), "\n")
res <- cause(X = X, variants = top, param_ests = params)
s <- summary(res, ci_size = 0.95)
print(s)
saveRDS(list(res = res, params = params, n_merged = nrow(X), n_pruned = length(top)), "results/R_AHG/cause_results.rds")
tab <- as.data.table(s$tab); fwrite(tab, "results/R_AHG/cause_posteriors.csv")
fwrite(as.data.table(res$elpd), "results/R_AHG/cause_elpd.csv")
print(res$elpd)
