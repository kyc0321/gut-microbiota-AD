## 23a_mvmr_data_local.R — assemble harmonised SNP x exposure matrix for MVMR from local
## full summary statistics (GWAS Catalog FTP; OpenGWAS API unavailable on 2026-10-07).
## Exposures: periodontitis (FinnGen, 50 P<1e-5 instruments) + AD risk factors (P<5e-8):
##   smoking initiation  Liu 2019, GSCAN excl. 23andMe (GCST007474; effect allele = ALT)
##   type 2 diabetes     Xue 2018 (GCST006867, harmonised)
##   body mass index     Yengo 2018, GIANT + UKB (GCST006900)
##   LDL cholesterol     Graham 2021, GLGC European (GCST90239658, harmonised)
##   education (years)   UK Biobank (GCST90029013, harmonised)
##   type 1 diabetes     Chiou 2021 (GCST90014023, harmonised)
## Outcome: Bellenguez 2022 (local harmonised file). Clumping: PLINK 1.9, 1000G EUR; per exposure r2<0.001, 10 Mb;
## union r2<0.01, 1 Mb.
suppressPackageStartupMessages(library(data.table))
RFD <- "ext_data/mvmr_rf"
spec <- list(
  smoking = list(f = "SmkInit_Liu2019.txt.gz", c = c(rsid = "RSID", ea = "ALT", oa = "REF", beta = "BETA", se = "SE", p = "PVALUE")),
  T2D     = list(f = "T2D.h.tsv.gz", c = c(rsid = "hm_rsid", ea = "hm_effect_allele", oa = "hm_other_allele", beta = "hm_beta", se = "standard_error", p = "p_value")),
  BMI     = list(f = "BMI_Yengo2018.txt.gz", c = c(rsid = "SNP", ea = "Tested_Allele", oa = "Other_Allele", beta = "BETA", se = "SE", p = "P")),
  LDL     = list(f = "LDL.h.tsv.gz", c = c(rsid = "rsid", ea = "effect_allele", oa = "other_allele", beta = "beta", se = "standard_error", p = "p_value")),
  EA      = list(f = "EA.h.tsv.gz", c = c(rsid = "hm_rsid", ea = "hm_effect_allele", oa = "hm_other_allele", beta = "hm_beta", se = "standard_error", p = "p_value")),
  T1D     = list(f = "T1D.h.tsv.gz", c = c(rsid = "hm_rsid", ea = "hm_effect_allele", oa = "hm_other_allele", beta = "hm_beta", se = "standard_error", p = "p_value"))
)
read_rf <- function(k) {
  s <- spec[[k]]
  d <- fread(file.path(RFD, s$f), select = unname(s$c), showProgress = FALSE)
  setnames(d, unname(s$c), names(s$c))
  d[, `:=`(ea = toupper(ea), oa = toupper(oa), beta = as.numeric(beta), se = as.numeric(se), p = as.numeric(p))]
  d[grepl("^rs", rsid) & !is.na(beta) & !is.na(se) & se > 0][!duplicated(rsid)]
}
clump <- function(dt, tag, p1, r2 = "0.001", kb = "10000") {
  inp <- sprintf("ext_data/mvmr_clump_%s.txt", tag)
  fwrite(dt[, .(SNP = rsid, P = p)], inp, sep = "\t")
  system2("ext_data/plink/plink", c("--bfile", "ext_data/ld_ref/EUR", "--clump", inp, "--clump-p1", p1,
          "--clump-p2", "1", "--clump-r2", r2, "--clump-kb", kb, "--out", sub(".txt$", "", inp)), stdout = FALSE)
  out <- sub(".txt$", ".clumped", inp)
  if (!file.exists(out)) return(character(0))
  fread(out)$SNP
}

RF <- names(spec)
skip <- commandArgs(trailingOnly = TRUE)   # e.g. "T1D" while its download is incomplete
avail <- setdiff(RF[file.exists(file.path(RFD, sapply(spec, `[[`, "f")))], skip)
rf <- setNames(lapply(avail, read_rf), avail)
iv <- readRDS("data/perio_instruments_finngen2025_relaxed.rds")

## 1. RF instruments, clumped separately
ivs <- lapply(avail, function(k) { s <- clump(rf[[k]][p < 5e-8], k, "5e-8"); cat(k, length(s), "instruments\n"); s })
names(ivs) <- avail

## 2. Union, then joint clumping on the minimum P across exposures
cand <- unique(c(iv$SNP, unlist(ivs)))
pmin_dt <- rbindlist(c(list(data.table(rsid = iv$SNP, p = iv$pval.exposure)),
                       lapply(avail, function(k) rf[[k]][rsid %in% cand, .(rsid, p)])))[, .(p = min(p)), by = rsid]
## joint clumping at r2 < 0.01 within 1 Mb: the stricter single-exposure setting (r2 < 0.001, 10 Mb)
## removed 43 of 50 periodontitis instruments in favour of far stronger risk-factor variants
snps <- clump(pmin_dt, "union", "1", r2 = "0.01", kb = "1000")
cat("Union after joint clumping:", length(snps), "\n")

## 3. Outcome and periodontitis at union SNPs
ad <- fread("ext_data/bellenguez/35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
            select = c("hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","hm_effect_allele_frequency"))
setnames(ad, c("rsid","ea","oa","b_AD","se_AD","eaf"))
M <- ad[rsid %in% snps & !is.na(b_AD)][!duplicated(rsid)]
pe <- fread("../FinGen Download/salminen_2025_sumstats/periodontitis_diagnosis_results.txt.gz",
            select = c("rsids","ref","alt","beta","sebeta","pval"))
pe[, rsid := sub(",.*", "", rsids)]
pe <- pe[rsid %in% snps][!duplicated(rsid)]
align <- function(b, ea, oa, rea, roa) fifelse(ea == rea & oa == roa, b, fifelse(ea == roa & oa == rea, -b, NA_real_))
M <- merge(M, pe[, .(rsid, pea = alt, poa = ref, b_perio = beta, se_perio = sebeta)], by = "rsid")
M[, b_perio := align(b_perio, pea, poa, ea, oa)][, c("pea","poa") := NULL]
for (k in avail) {
  x <- rf[[k]][rsid %in% snps, .(rsid, kea = ea, koa = oa, kb = beta, kse = se)]
  M <- merge(M, x, by = "rsid", all.x = TRUE)
  M[, (paste0("b_", k)) := align(kb, kea, koa, ea, oa)][, (paste0("se_", k)) := kse]
  M[, c("kea","koa","kb","kse") := NULL]
}
pal <- (M$ea %in% c("A","T") & M$oa %in% c("A","T")) | (M$ea %in% c("C","G") & M$oa %in% c("C","G"))
M <- M[!(pal & eaf > 0.42 & eaf < 0.58)]
M[, perio_iv := rsid %in% iv$SNP]
for (k in avail) M[, (paste0(k, "_iv")) := rsid %in% ivs[[k]]]
fwrite(M, "results/R_AHG/mvmr_harmonised_matrix.csv")
cat("Final matrix:", nrow(M), "SNPs; periodontitis IVs retained:", sum(M$perio_iv), "\n")
print(M[, lapply(.SD, function(z) sum(!is.na(z))), .SDcols = patterns("^b_")])
